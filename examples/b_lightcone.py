from multiprocessing.pool import Pool
import pandas as pd
from pandas.core.frame import DataFrame
import numpy as np
import time as ti
import random
from scipy.integrate import quad
import scipy.optimize as opt
import tables
import h5py
import math


def z2t(z):
    func1 = lambda x: 1 / ((1 + x) * (O_M * ((1 + x) ** 3) + O_K * ((1 + x) ** 2) + O_L) ** 0.5)
    v1 = quad(func1, 0, z)
    t_L = t_H * v1[0]  # look back time in Gyr/h
    return t_L


def z2r(z):
    func2 = lambda y: 1 / (O_M * ((1 + y) ** 3) + O_K * ((1 + y) ** 2) + O_L) ** 0.5
    v2 = quad(func2, 0, z)
    r_c = c * v2[0] / 1e5  # comoving distance in Mpc/h
    return r_c


def fss(begin, end):  # 找盒子位于哪几个snapshot之间,返回两个参数
    ben = 0
    enn = 0
    for l1 in range(ssn_end + 1, ssn):
        if begin == 0:
            ben = ssn - 1
            break
        elif rcp[l1] <= begin:
            ben = l1  # begin snapshot number
            break
    for l2 in range(ssn - 1, ssn_end, -1):
        if rcp[l2] >= end:
            enn = l2
            break  # 总是先找到begin，才能用
    return ben, enn


def check_z(z):  # 检查一下是不是多加了1？应该没问题
    for ch in range(ssn - 1, ssn_end, -1):
        if z < zcp[ch]:
            return ch


def car_to_RA(x, y, z):
    dec = math.pi / 2 - math.atan2(math.sqrt(x ** 2 + y ** 2), z)
    RA = math.atan2(y, x)
    if RA < 0:
        RA += math.pi * 2
    return RA * 180 / np.pi, dec * 180 / np.pi


def judge_angle_pi(v0, v1):  # judge two vectors' angle between 0 and pi
    ang = math.atan2(np.linalg.norm(np.cross(v0, v1)), np.dot(np.transpose(v0), v1))
    return ang  # return in rad. not in degree


def in_cone_snap(posc, snapnum_now):  # 判断000点会有一点问题
    # if snapnum_now == ssn - 1:
    #     return False
    rc = np.linalg.norm(posc)

    if cone_geo == 0:
        if judge_angle_pi(los, posc) < hfa and rc < r_max and rcp[snapnum_now - 1] >= rc > rcp[snapnum_now]:
            return True
        else:
            return False
    elif cone_geo == 1:
        ra, dec = car_to_RA(posc[0], posc[1], posc[2])
        if ra_min < ra < ra_max and dec_min < dec < dec_max and rc < r_max and rcp[snapnum_now - 1] >= rc > rcp[
            snapnum_now]:
            return True
        else:
            return False
    elif cone_geo == 2:
        if rc < r_max and rcp[snapnum_now - 1] >= rc > rcp[snapnum_now]:
            return True
        else:
            return False
    else:
        return True


def magcalc(Mag, distance):
    mag_now = Mag + 25 + 5 * np.log10(distance)
    return mag_now


def interpo(var0, var1, delta):
    """
    :param var0: the start point,
    :param var1: the end point,
    :param delta: time between 0 and 1,
    :return: linear interpolated value
    """
    var_now = var0 + (var1 - var0) * delta
    return var_now


def compare(var0, var1):
    if var0 <= var1:
        return 1
    else:
        return 0


def cone(paras):
    zi = paras[0]
    yi = paras[1]
    xi = paras[2]
    a = paras[3]
    b = paras[4]
    id2 = paras[5]
    snap = ga[id2].snapnum
    st = snap[0]

    ret = []
    # 加入random tiling 20220816
    phase = []  # 3 digits, present mirror x, y or z axis
    for loop99 in range(3):
        # if random.random() < 0.5:
        #     phase.append([0, 1])  # do not mirror
        # else:
        phase.append([0, 1])  # do mirror

    incone = []
    for loop13 in range(a, b - 1, -1):  # 红移由低到高
        snap_index = np.argwhere(snap == loop13).squeeze()
        if snap_index.size:
            # p0 = np.zeros(3)
            # exec('p0 = ga[id2].%s[snap_index]' % pos_name[0])
            p0 = ga[id2].ComovingAveragePosition[snap_index]
            p0[0] = (phase[0][0] * boxsize + phase[0][1] * p0[0]) % boxsize + xi * boxsize - obs[0]
            p0[1] = (phase[1][0] * boxsize + phase[1][1] * p0[1]) % boxsize + yi * boxsize - obs[1]
            p0[2] = (phase[2][0] * boxsize + phase[2][1] * p0[2]) % boxsize + zi * boxsize - obs[2]
            if in_cone_snap(p0, loop13):  # 在红移空间的中点？
                incone.append(loop13)
            else:
                incone.append(0)
        else:
            incone.append(0)  # 保证len(incone) == a-b+1
    incone = np.array(incone)

    if len(incone) >= 2 and np.linalg.norm(incone):
        for loop1 in range(len(incone) - 1):
            sit = np.argwhere(snap == incone[loop1]).squeeze()
            if sit.size and sit:  ###
                ssnow = int(sit + st)
                if incone[loop1 + 1] == 0:  # incone[loop1] != 0 and incone[loop1+1] == 0:  # 可以确定在这个时间范围内被看到
                    func0 = lambda tcos: \
                        (phase[0][0] * boxsize + phase[0][1] * (
                                    ga[id2].coeffX[ssnow - st - 1][0] * tcos ** 3 + ga[id2].coeffX[ssnow - st - 1][
                                1] * tcos ** 2 + ga[id2].coeffX[ssnow - st - 1][2] * tcos +
                                    ga[id2].coeffX[ssnow - st - 1][3]) + xi * boxsize - obs[0]) ** 2 + \
                        (phase[1][0] * boxsize + phase[1][1] * (
                                    ga[id2].coeffY[ssnow - st - 1][0] * tcos ** 3 + ga[id2].coeffY[ssnow - st - 1][
                                1] * tcos ** 2 + ga[id2].coeffY[ssnow - st - 1][2] * tcos +
                                    ga[id2].coeffY[ssnow - st - 1][3]) + yi * boxsize - obs[1]) ** 2 + \
                        (phase[2][0] * boxsize + phase[2][1] * (
                                    ga[id2].coeffZ[ssnow - st - 1][0] * tcos ** 3 + ga[id2].coeffZ[ssnow - st - 1][
                                1] * tcos ** 2 + ga[id2].coeffZ[ssnow - st - 1][2] * tcos +
                                    ga[id2].coeffZ[ssnow - st - 1][3]) + zi * boxsize - obs[2]) ** 2 - \
                        (coefs_t2r[incone[loop1]][0] * tcos ** 3 + coefs_t2r[incone[loop1]][1] * tcos ** 2 +
                         coefs_t2r[incone[loop1]][2] * tcos + coefs_t2r[incone[loop1]][3]) ** 2
                    t_root = np.array(opt.fsolve(func0, x0=ga[id2].time[ssnow - st - 1]))

                    if t_root.size and t_cos - tcp[incone[loop1] - 1] <= t_root[0] <= t_cos - tcp[
                        incone[loop1]]:  # 找到了符合要求的解
                        aa = zcp[incone[loop1] - 1]  # start point, high redshift
                        bb = zcp[incone[loop1]]  # end point, low redshift
                        tcp0 = z2t(aa)
                        tcp1 = z2t(bb)
                        delta_t = ((t_cos - t_root) - tcp0) / (tcp1 - tcp0)  # both negative values
                        if delta_t < 0 or delta_t > 1:
                            print('warning: delta_t = ', delta_t)
                            print(id2, t_root, tcp0, tcp1)
                        findout = np.zeros(3)
                        findout[0] = phase[0][0] * boxsize + phase[0][1] * (
                                    ga[id2].coeffX[ssnow - st - 1][0] * t_root ** 3 + ga[id2].coeffX[ssnow - st - 1][
                                1] * t_root ** 2 + ga[id2].coeffX[ssnow - st - 1][2] * t_root +
                                    ga[id2].coeffX[ssnow - st - 1][3]) + xi * boxsize - obs[0]
                        findout[1] = phase[1][0] * boxsize + phase[1][1] * (
                                    ga[id2].coeffY[ssnow - st - 1][0] * t_root ** 3 + ga[id2].coeffY[ssnow - st - 1][
                                1] * t_root ** 2 + ga[id2].coeffY[ssnow - st - 1][2] * t_root +
                                    ga[id2].coeffY[ssnow - st - 1][3]) + yi * boxsize - obs[1]
                        findout[2] = phase[2][0] * boxsize + phase[2][1] * (
                                    ga[id2].coeffZ[ssnow - st - 1][0] * t_root ** 3 + ga[id2].coeffZ[ssnow - st - 1][
                                1] * t_root ** 2 + ga[id2].coeffZ[ssnow - st - 1][2] * t_root +
                                    ga[id2].coeffZ[ssnow - st - 1][3]) + zi * boxsize - obs[2]
                        rc_true = np.linalg.norm(findout)
                        if in_cone_snap(findout, incone[loop1]):
                            ra_true, dec_true = car_to_RA(findout[0], findout[1], findout[2])
                            zfind = coefs_r2z[incone[loop1]][0] * rc_true ** 3 + coefs_r2z[incone[loop1]][
                                1] * rc_true ** 2 + coefs_r2z[incone[loop1]][2] * rc_true + coefs_r2z[incone[loop1]][3]
                            # return properties
                            trackID = ga[id2].TrackId[ssnow - st]
                            hosthaloID = ga[id2].HostHaloId[ssnow - st]
                            rank = ga[id2].Rank[ssnow - st]
                            posx = findout[0]
                            posy = findout[1]
                            posz = findout[2]
                            velx = interpo(ga[id2].PhysicalAverageVelocity[ssnow - st - 1][0],
                                           ga[id2].PhysicalAverageVelocity[ssnow - st][0], delta_t)[0]  ##
                            vely = interpo(ga[id2].PhysicalAverageVelocity[ssnow - st - 1][1],
                                           ga[id2].PhysicalAverageVelocity[ssnow - st][1], delta_t)[0]
                            velz = interpo(ga[id2].PhysicalAverageVelocity[ssnow - st - 1][2],
                                           ga[id2].PhysicalAverageVelocity[ssnow - st][2], delta_t)[0]
                            d_comoving = rc_true
                            v_lfs = (posx * velx + posy * vely + posz * velz) / d_comoving
                            ra = ra_true
                            dec = dec_true
                            vmax = ga[id2].VmaxPhysical[ssnow - st]
                            PeakMass = ga[id2].LastMaxMass[ssnow - st]
                            PeakVmax = ga[id2].LastMaxVmaxPhysical[ssnow - st]
                            shMbound = interpo(ga[id2].Mbound[ssnow - st - 1], ga[id2].Mbound[ssnow - st], delta_t)[0]
                            shBoundM200Crit = interpo(ga[id2].BoundM200Crit[ssnow - st - 1],
                                                      ga[id2].BoundM200Crit[ssnow - st], delta_t)[0]
                            redshift_true = zfind
                            redshift_obs = zfind + v_lfs * (1 + zfind) * 1e3 / c
                            snapNum = ssnow
                            iso = np.argwhere(snap == last_iso_list[int(trackID % 5e5)]).squeeze()  # 对于低红移光锥，可能有点问题
                            if iso.size:
                                shMbound_at_ac = ga[id2].Mbound[iso]
                            else:
                                shMbound_at_ac = ga[id2].Mbound[0]

                            for loop1 in range(len(output_list)):
                                exec('ret.append(%s)' % output_list[loop1])
                            return ret

    # elif len(incone) == 1 and incone[0]:  # 只有一个点，还在光锥内，不能丢掉，直接待在原地吧
    #     sit = np.argwhere(snap == incone[0]).squeeze()
    #     ssnow = int(sit + st)  # 在class中真实的位置
    #     p0 = np.zeros(3)
    #     for loop4 in range(3):
    #         exec('p0[loop4] = ga[id2].%s[ssnow - st][loop4]' % pos_name[0])  ##
    #     p0[0] = phase[0][0] * boxsize + phase[0][1] * p0[0] + xi * boxsize - obs[0]
    #     p0[1] = phase[1][0] * boxsize + phase[1][1] * p0[1] + yi * boxsize - obs[1]
    #     p0[2] = phase[2][0] * boxsize + phase[2][1] * p0[2] + zi * boxsize - obs[2]
    #     d_now = np.linalg.norm(p0)
    #     ra_true = 180 * math.atan(p0[1] / p0[0]) / np.pi  # 角度制
    #     dec_true = 180 * math.atan(p0[2] / np.sqrt(p0[1] ** 2 + p0[0] ** 2)) / np.pi
    #     halo_mass = ga[id2].mvir[ssnow - st]
    #     stellar_mass = ga[id2].stellarMass[ssnow - st]
    #     if ga[id2].type[ssnow - st] == 0:
    #         is_central = 1
    #     else:
    #         is_central = 0
    #     zfind = coefs_r2z[incone[0]][0] * d_now ** 3 + coefs_r2z[incone[0]][1] * d_now ** 2 + \
    #             coefs_r2z[incone[0]][2] * d_now + coefs_r2z[incone[0]][3]  # 告别 r2z时代！ 20220823
    #     D_L = d_now * (1 + zfind) / h
    #     for loop1 in range(len(output_list)):
    #         exec('ret.append(%s)' % output_list[loop1])
    #     return ret


class galaxy():
    def __init__(self, trackid):
        trackid = int(trackid % 5e5)
        start_ssn = int(born[k][trackid])
        for p in fields:  # p 为名字
            exec('%s = []' % p)
        live = []
        livetime = []
        for ii in range(k, ssn):  # ssn 从小到大读、存  0对应刚出生的ssn  # 改为k
            if ii == 69 or ii == 76:  # ?
                continue
            if 0 < deadlist[trackid] < ii:  # 捕捉死亡信息
                break
            for p in range(len(fields)):  # p为数字
                exec("%s.append(data[ii][trackid]['%s'])" % (fields[p], fields[p]))
            live.append(ii)
            livetime.append(t_cos - tcp[ii])  # 便于插值，换算为以z=3000时为时间零点计算
        for p in fields:
            exec('self.%s = np.array(%s)' % (p, p))  # 2022.7.1新增np.array
        self.snapnum = np.array(live)
        self.time = np.array(livetime)
        self.coeffX = fit_coeff(self.time, self.ComovingAveragePosition[:, 0], self.PhysicalAverageVelocity[:, 0])
        self.coeffY = fit_coeff(self.time, self.ComovingAveragePosition[:, 1], self.PhysicalAverageVelocity[:, 1])
        self.coeffZ = fit_coeff(self.time, self.ComovingAveragePosition[:, 2], self.PhysicalAverageVelocity[:, 2])


def fit_coeff(times, pos, vel):  # unit of times: Gyr/h, unit of pos: Mpc/h, unit of vel: km/s
    paras = []
    for loop1 in range(len(times) - 1):
        x0 = pos[loop1]
        x1 = pos[loop1 + 1]
        t0 = times[loop1]
        t1 = times[loop1 + 1]
        delta_x = x1 - x0
        if delta_x > boxsize / 2:
            x1 -= boxsize
        elif delta_x < - boxsize / 2:
            x1 += boxsize
        paras.append(solve_coefficients(t0, t1, x0, x1, vel[loop1], vel[loop1 + 1]))
    return np.array(paras)


def read_data(id3):  # 并行读
    return galaxy(id3)


def solve_coefficients(t1, t2, x1, x2, x1_dot, x2_dot):
    left = np.array([[t1 ** 3, t1 ** 2, t1, 1],
                     [3 * t1 ** 2, 2 * t1, 1, 0],
                     [t2 ** 3, t2 ** 2, t2, 1],
                     [3 * t2 ** 2, 2 * t2, 1, 0]])
    right = np.array([x1, x1_dot / 978.5, x2, x2_dot / 978.5])  # 978.5km/s = 1 Mpc/Gyr
    coeff = np.linalg.solve(left, right)
    return coeff


if __name__ == '__main__':
    start = ti.time()
    timenow = ti.gmtime()
    if timenow.tm_hour < 16:
        st_time = int(timenow.tm_mon * 1e4 + timenow.tm_mday * 1e2 + timenow.tm_hour + 8)
    else:
        st_time = int(timenow.tm_mon * 1e4 + (timenow.tm_mday + 1) * 1e2 + timenow.tm_hour - 16)
    print('Program start at : ', st_time)
    t_H = 9785641806 * 1e-9  # Hubble time in Gyr/h 要根据模拟输出调整
    c = 299792458
    pc = 3.08568025e16

    # parallel computing by galaxy bins
    iii = 2661  ###
    ssn_start = 118  ###

    # paths
    program_path = '/home/zltan/9t/'
    loadpath = '/home/cossim/Jiutian/M1000/hbt-tab/'
    redshift_path = program_path + 'redshift.csv'
    boxes_path = program_path + 'boxcounts_2_2_pacth.npy'
    save_path = '/home/cossim/Jiutian/M1000/lightcones/fullsky_z2/patch/9tmock_z2_%d.hdf5' % iii

    # 模拟盒子信息，或许可以从文件读出来
    cos_case = 3  # 0 for Millennium-1, 1 for scPlanck
    O_M = 0.3111
    O_L = 1 - O_M
    O_K = 0
    boxsize = 1000  # boxsize in Mpc/h
    ssn = 128  # snapshot number(可以通过LEN(AT)读出来)

    # 距离，坐标统一单位Mpc/h，时间统一单位1e9*year
    pos_interpo_method = 0  # 0 for cubic interpolation, 1 for cubic spline interpolation
    props_interpo_method = 0  # 0 for linear interpolation, 1 for cubic spline interpolation
    process = 72  ### 并行核数

    # 依据需要的光锥条件调整
    t_cos = z2t(3000)  # 宇宙年龄
    z_MAX = 10  # 链表能够索引的最大红移，需要比截断红移略大  ###
    z_max = 2  # 输入，截断红移 ###

    # cone_geo = 0  # 0 for lfs+hfa judge, 1 for R.A.+dec judge, 2 for full sky
    cone_geo = 2
    hfa = -1
    los = -1
    ra_min, ra_max, dec_min, dec_max = -1, -1, -1, -1
    obs = [500, 500, 500]

    r_max = z2r(z_max)
    N = int(r_max / boxsize) + 2

    boxes = np.load(boxes_path)
    boxes_len = len(boxes)
    print('Z_max =', z_max, 'N =', N, 'Box count =', boxes_len)
    '''最终输出文件的列名'''
    output_list = ['trackID', 'hosthaloID', 'rank', 'posx', 'posy', 'posz', 'velx', 'vely', 'velz', 'v_lfs', 'shMbound',
                   'd_comoving', 'ra', 'dec', 'vmax', 'PeakMass', 'PeakVmax', 'shBoundM200Crit', 'redshift_true',
                   'redshift_obs', 'shMbound_at_ac', 'snapNum']
    column = ['TrackId', 'Nbound', 'Mbound', 'HostHaloId', 'Rank', 'Depth', 'LastMaxMass', 'SnapshotIndexOfLastMaxMass',
              'SnapshotIndexOfLastIsolation', 'SnapshotIndexOfBirth', 'SnapshotIndexOfDeath', 'SnapshotIndexOfSink',
              'RmaxComoving', 'VmaxPhysical', 'LastMaxVmaxPhysical', 'SnapshotIndexOfLastMaxVmax', 'R2SigmaComoving',
              'RHalfComoving', 'BoundR200CritComoving', 'BoundM200Crit', 'SpecificSelfPotentialEnergy',
              'SpecificSelfKineticEnergy', 'SpecificAngularMomentum', 'InertialEigenVector',
              'InertialEigenVectorWeighted', 'InertialTensor', 'InertialTensorWeighted', 'ComovingAveragePosition',
              'PhysicalAverageVelocity', 'ComovingMostBoundPosition', 'PhysicalMostBoundVelocity',
              'MostBoundParticleId', 'SinkTrackId']
    fields = ['TrackId', 'Mbound', 'HostHaloId', 'Rank', 'LastMaxMass', 'SnapshotIndexOfBirth', 'SnapshotIndexOfDeath',
              'RmaxComoving', 'VmaxPhysical', 'LastMaxVmaxPhysical', 'RHalfComoving', 'BoundR200CritComoving',
              'BoundM200Crit', 'ComovingAveragePosition', 'PhysicalAverageVelocity', 'SnapshotIndexOfLastIsolation']

    # peak-mass -> LastMaxMass?
    # peak-vmax -> LastMaxVmaxPhysical
    # mass -> BoundM200Crit
    # not_used_name = ['Nbound', 'Mbound', 'Depth', 'SnapshotIndexOfLastMaxMass', 'SnapshotIndexOfLastIsolation', 'SnapshotIndexOfSink',
    #                  'R2SigmaComoving', 'SpecificSelfPotentialEnergy', 'SpecificSelfKineticEnergy',
    #                  'SpecificAngularMomentum', 'InertialEigenVector', 'InertialEigenVectorWeighted', 'InertialTensor',
    #                  'InertialTensorWeighted', 'ComovingMostBoundPosition', 'PhysicalMostBoundVelocity',
    #                  'MostBoundParticleId', 'SinkTrackId', 'SnapshotIndexOfLastMaxVmax', ]
    # Mag_name = ['SDSSu_Dust', 'SDSSg_Dust', 'SDSSr_Dust', 'SDSSi_Dust', 'SDSSz_Dust', 'Y_Dust', 'GNUV_Dust']
    # Mag_bulge_name = ['SDSSu_BulgeDust', 'SDSSg_BulgeDust', 'SDSSr_BulgeDust', 'SDSSi_BulgeDust', 'SDSSz_BulgeDust', 'Y_BulgeDust', 'GNUV_BulgeDust']

    redshifts = pd.read_table(redshift_path, header=0)
    zcp = redshifts['z']  # 从大到小
    zcp[69] = zcp[68]
    zcp[76] = zcp[75]

    ssn_end = max(min(boxes[:, 4]), ssn_start)  ### 截止读取、计算的ssn
    ssn_st = min(max(boxes[:, 3]), ssn)  ### 20230202 新增
    # 空列表准备读数据
    rcp = np.zeros(ssn)
    tcp = np.zeros(ssn)
    P = [[[] for row1 in range(3)] for row in range(ssn)]
    trackIDs = [[] for row2 in range(ssn)]
    born = [[] for row5 in range(ssn)]
    pos_name = ['ComovingAveragePosition']
    vel_name = ['PhysicalAverageVelocity']

    data = [[] for row in range(ssn)]
    # 读、存数据都是倒序的，变量最后存的是ssn最大的内容
    lines = len(column)
    ### 如果读数据读爆内存，可以在每一个data中删除一些不用的数据
    coefs_t2r = [[] for row2 in range(ssn)]
    coefs_r2z = [[] for row2 in range(ssn)]
    for i in range(ssn_end, ssn):
        rcp[i] = z2r(zcp[i])
        tcp[i] = z2t(zcp[i])
        if i == 69 or i == 76:  #
            continue

        rr = np.array([z2r(i2) for i2 in np.arange(zcp[i], zcp[i - 1], 1e-4)])  # z是在增加的
        tt = np.array([z2t(i2) for i2 in np.arange(zcp[i], zcp[i - 1], 1e-4)])
        zz = np.array([i2 for i2 in np.arange(zcp[i], zcp[i - 1], 1e-4)])
        tt = t_cos - tt  # 转换为宇宙年龄
        coefs_t2r[i] = np.polyfit(tt, rr, 3)  # 宇宙年龄-共动距离的关系拟合，第63个存的snapshot63到62的拟合公式
        coefs_r2z[i] = np.polyfit(rr, zz, 3)

        data[i] = h5py.File(loadpath + "./%03d/SubSnap%03d_%d.hdf5" % (i, i, iii), "r")['Subhalos'][:][fields]
        if i == ssn - 1:  #
            livesnaps = np.array(data[i]['SnapshotIndexOfDeath'] - data[i]['SnapshotIndexOfBirth']).squeeze()
            livelist = (livesnaps < 0) ^ (livesnaps > 2)  # short-aged subhalos are filtered
            print('Not short aged subhalos conuts: ', np.sum(livelist))
            deadlist = data[i]['SnapshotIndexOfDeath']
            last_iso_list = data[i]['SnapshotIndexOfLastIsolation']

        trackIDs[i] = np.array(data[i]['TrackId'])
        born[i] = np.array(data[i]['SnapshotIndexOfBirth'])
        print('len(data[%d]): %d' % (i, len(trackIDs[i])))
    print(zcp, rcp, tcp)

    ga = []
    for k in range(ssn_end, ssn_st):  ### ssn
        if k == 76 or k == 69:
            continue
        need_read = []
        Ng = len(ga)  # N galaxies
        if k == ssn_end:
            for j in range(len(trackIDs[k])):
                if int(born[k][j]) <= k and int(born[k][j]) != -1 and livelist[j] and (
                        deadlist[j] >= k or deadlist[j] == -1):
                    need_read.append(trackIDs[k][j])
        else:
            for j in range(len(trackIDs[k])):
                if int(born[k][j]) == k and livelist[j] and (deadlist[j] >= k or deadlist[j] == -1):
                    need_read.append(trackIDs[k][j])
        if len(need_read):
            print('ssn=', k, 'need read=', len(need_read))
            pool = Pool(processes=process)
            gal = pool.map(read_data, need_read)
            pool.close()
            pool.join()
            if Ng and len(gal):
                ga = np.concatenate((ga, gal), axis=0)
            elif len(gal):
                ga = gal
    ga_len = len(ga)
    print('读取结束：', ga_len, ', 共计%d个盒子' % boxes_len, '耗时：%.1f' % (ti.time() - start))
    if not ga_len:
        print('No objects to make lightcone!')
        exit(0)
    # 开始 build lightcone
    # box_now = 0
    fake_id = [[] for i in range(ga_len)]
    for i in range(ga_len):
        fake_id[i].append(i)
    time_now = ti.time()
    para = [[] for i in range(boxes_len)]
    for jj in range(boxes_len):
        box = np.tile(boxes[jj], (ga_len, 1))
        new = np.concatenate((box, fake_id), axis=1)
        para[jj] = new
    print('para[jj] = new, used time: %.1f s' % (ti.time() - time_now))

    para = np.array(para)
    para = para.reshape(boxes_len * ga_len, 6)  # 生成一批6参数数组， 然后并行处理
    print('para = para.reshape(boxes_len * ga_len, 6), used time: %.1f s' % (ti.time() - time_now))

    pool = Pool(processes=process)
    res = pool.map(cone, para)  # 只传递id就可以了，但是要传长度最小的，否则调用会出错  ### np.array
    pool.close()
    pool.join()
    print('res = pool.map(cone, para), used time: %.1f s' % (ti.time() - time_now))

    gal_tot = [i2 for i2 in res if i2 is not None]  # res = list(filter(None, res))

    print('gal_tot = [i2 for i2 in res if i2 is not None], used time: %.1f s' % (ti.time() - time_now), len(gal_tot))

    gal_tot = DataFrame(gal_tot, columns=output_list)  # , dtype=object)
    gal_tot = pd.DataFrame(gal_tot).to_records(index=None)
    print('gal_tot = pd.DataFrame(gal_tot).to_records(index=None), used time: %.1f s' % (ti.time() - time_now))

    with h5py.File(save_path, "w") as output_file:
        output_file.create_dataset('Subhalos', len(gal_tot),
                                   dtype=[('trackID', '<i8'), ('hosthaloID', '<i8'), ('rank', '<i8'), ('posx', '<f4'),
                                          ('posy', '<f4'), ('posz', '<f4'), ('velx', '<f4'), ('vely', '<f4'),
                                          ('velz', '<f4'), ('v_lfs', '<f4'), ('shMbound', '<f4'),
                                          ('d_comoving', '<f4'), ('ra', '<f4'), ('dec', '<f4'), ('vmax', '<f4'),
                                          ('PeakMass', '<f4'), ('PeakVmax', '<f4'),
                                          ('shBoundM200Crit', '<f4'), ('redshift_true', '<f4'),
                                          ('redshift_obs', '<f4'), ('shMbound_at_ac', '<f4'),
                                          ('snapNum', '<i4')])
        output_file['Subhalos'][:] = gal_tot
    output_file.close()
    print('Len: ', len(gal_tot), 'Total used: %.1f s' % (ti.time() - start))
    print('output_docu=', save_path)
