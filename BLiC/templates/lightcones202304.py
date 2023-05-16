# 20230428: 加入不插值版本光锥（直接拼接）
# interpo_method = 0 !!!

import sys
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
from pathlib import Path
import glob
import pickle
import os


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

    # if cone_geo == 0:
    #     if judge_angle_pi(los, posc) < hfa and rc < r_max and rcp[snapnum_now - 1] >= rc > rcp[snapnum_now]:
    #         return True
    #     else:
    #         return False

    if full_sky:
        if rc < r_max and rcp[snapnum_now - 1] >= rc > rcp[snapnum_now]:
            return True
        else:
            return False
    else:
        ra, dec = car_to_RA(posc[0], posc[1], posc[2])
        if ra_min < ra < ra_max and dec_min < dec < dec_max and rc < r_max and rcp[snapnum_now - 1] >= rc > rcp[
            snapnum_now]:
            return True
        else:
            return False


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
    zi = int(paras[0])
    yi = int(paras[1])
    xi = int(paras[2])
    a = int(paras[3])
    b = int(paras[4])
    rdseed = [paras[5], paras[6], paras[7]]
    id2 = int(paras[8])

    snap = ga[id2].snapnum
    st = snap[0]

    ret = []
    # 加入random tiling 20220816
    phase = []  # 3 digits, present mirror x, y or z axis
    for loop99 in range(3):
        if rdseed[loop99] < 0.5:
            phase.append([0, 1])  # do not mirror
        else:
            phase.append([1, -1])  # do mirror

    incone = []
    for loop13 in range(a, b - 1, -1):  # 红移由低到高
        snap_index = np.argwhere(snap == loop13).squeeze()
        if snap_index.size:
            # p0 = np.zeros(3)
            if len(cb.pos_name) == 1:
                exec("p00 = ga[id2].%s[snap_index] * cb.pos_unit" % cb.pos_name[0])
            else:
                for ii in range(3):
                    p00 = np.zeros(3)
                    exec("p00[ii] = ga[id2].%s[snap_index] * cb.pos_unit" % cb.pos_name[ii])
            p0 = locals()['p00']
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
                    if interpo_method[-1]:
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
                    else:
                        if len(cb.pos_name) == 1:
                            exec("p00 = ga[id2].%s[ssnow - st] * cb.pos_unit" % cb.pos_name[0])
                        else:
                            for ii in range(3):
                                exec("p00[ii] = ga[id2].%s[ssnow - st] * cb.pos_unit" % cb.pos_name[ii])
                        p0 = locals()['p00']
                        p0[0] = (phase[0][0] * boxsize + phase[0][1] * p0[0]) % boxsize + xi * boxsize - obs[0]
                        p0[1] = (phase[1][0] * boxsize + phase[1][1] * p0[1]) % boxsize + yi * boxsize - obs[1]
                        p0[2] = (phase[2][0] * boxsize + phase[2][1] * p0[2]) % boxsize + zi * boxsize - obs[2]
                        func0 = lambda tcos: \
                            (p0[0]) ** 2 + \
                            (p0[1]) ** 2 + \
                            (p0[2]) ** 2 - \
                            (coefs_t2r[incone[loop1]][0] * tcos ** 3 + coefs_t2r[incone[loop1]][1] * tcos ** 2 +
                             coefs_t2r[incone[loop1]][2] * tcos + coefs_t2r[incone[loop1]][3]) ** 2
                        t_root = np.array(opt.fsolve(func0, x0=ga[id2].time[ssnow - st - 1]))

                    if t_root.size and t_cos - tcp[incone[loop1] - 1] <= t_root[0] <= t_cos - tcp[
                        incone[loop1]]:  # 找到了符合要求的解
                        aa = zcp[incone[loop1] - 1]  # start point, high redshift
                        bb = zcp[incone[loop1]]  # end point, low redshift
                        tcp0 = z2t(aa)
                        tcp1 = z2t(bb)
                        # 20230221
                        if interpo_method[-1]:
                            delta_t = np.array(
                                ((t_cos - t_root[0]) - tcp0) / (tcp1 - tcp0)).squeeze()  # both negative values
                        else:
                            delta_t = 1

                        if delta_t < 0 or delta_t > 1:
                            print('warning: delta_t = ', delta_t)
                            print(id2, t_root, tcp0, tcp1)
                        findout = np.zeros(3)
                        if interpo_method[-1]:
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
                        else:
                            findout = p0
                        rc_true = np.linalg.norm(findout)
                        if in_cone_snap(findout, incone[loop1]):
                            # return properties
                            RA, Dec = car_to_RA(findout[0], findout[1], findout[2])
                            z_cos = coefs_r2z[incone[loop1]][0] * rc_true ** 3 + coefs_r2z[incone[loop1]][
                                1] * rc_true ** 2 + coefs_r2z[incone[loop1]][2] * rc_true + coefs_r2z[incone[loop1]][3]
                            snapnum = ssnow
                            d_comoving = rc_true

                            for i, iprop in enumerate(cb.props):
                                if iprop == cb.pos_name or iprop in cb.pos_name or iprop == cb.vel_name or iprop in cb.vel_name:
                                    continue
                                elif cb.interpo[i]:
                                    # only need to change this if more ways of interpolation can be chosen
                                    # why a '[0]' at end ?? need to check
                                    exec("%s = interpo(ga[id2].%s[ssnow - st - 1], ga[id2].%s[ssnow - st], delta_t)" % (
                                        iprop, iprop, iprop))
                                else:
                                    exec("%s = ga[id2].%s[ssnow - st]" % (iprop, iprop))

                            VELOCITOES = np.zeros(3)
                            if len(cb.pos_name) == 1:
                                exec("%s = np.array([findout[0], findout[1], findout[2]])" % cb.pos_name[0])
                                for loop2 in range(3):
                                    exec(
                                        "VELOCITOES[loop2] = interpo(ga[id2].%s[ssnow - st - 1][loop2], ga[id2].%s[ssnow - st][loop2], delta_t)" % (
                                            cb.vel_name[0], cb.vel_name[0]))  # 为什么[ssnow-st] 后面有[0]？ 已经改为loop2
                                exec("%s = VELOCITOES" % cb.vel_name[0])
                            else:
                                for ii in range(3):
                                    exec("%s = findout[ii]" % cb.pos_name[ii])
                                    exec(
                                        "VELOCITOES[ii] = interpo(ga[id2].%s[ssnow - st - 1], ga[id2].%s[ssnow - st], delta_t)" % (
                                            cb.vel_name[ii], cb.vel_name[ii]))
                                    exec("%s = VELOCITOES[ii]" % cb.vel_name[ii])

                            v_los = (findout[0] * VELOCITOES[0] + findout[1] * VELOCITOES[1] + findout[2] * VELOCITOES[
                                2]) / d_comoving
                            z_obs = z_cos + v_los * (1 + z_cos) * 1e3 / cb.c

                            for loop3 in output_list:
                                exec('ret.append(%s)' % loop3)
                            return ret


class galaxy():
    def __init__(self, trackid):
        trackid = int(trackid)
        start_ssn = k
        for p in fields:  # p 为名字
            exec('%s = []' % p)
        live = []
        livetime = []
        for ii in range(k, ssn):  # ssn 从小到大读、存  0对应刚出生的ssn  # 改为k
            # if ii == 69 or ii == 76:  # ?
            #     continue
            for p in range(len(fields)):  # p为数字
                exec("%s.append(data[ii][trackid]['%s'])" % (fields[p], fields[p]))
            live.append(ii)
            livetime.append(t_cos - tcp[ii])  # 便于插值，换算为以z=3000时为时间零点计算
        for p in fields:
            exec('self.%s = np.array(%s)' % (p, p))  # 2022.7.1新增np.array
        self.snapnum = np.array(live)
        self.time = np.array(livetime)
        # fitting traces
        if len(cb.pos_name) == 1 and len(cb.vel_name) == 1:
            exec("self.coeffX = fit_coeff(self.time, self.%s[:, 0], self.%s[:, 0])" % (cb.pos_name[0], cb.vel_name[0]))
            exec("self.coeffY = fit_coeff(self.time, self.%s[:, 1], self.%s[:, 1])" % (cb.pos_name[0], cb.vel_name[0]))
            exec("self.coeffZ = fit_coeff(self.time, self.%s[:, 2], self.%s[:, 2])" % (cb.pos_name[0], cb.vel_name[0]))
        elif len(cb.pos_name) == 3 and len(cb.vel_name) == 3:
            exec("self.coeffX = fit_coeff(self.time, self.%s, self.%s)" % (cb.pos_name[0], cb.vel_name[0]))
            exec("self.coeffY = fit_coeff(self.time, self.%s, self.%s)" % (cb.pos_name[1], cb.vel_name[1]))
            exec("self.coeffZ = fit_coeff(self.time, self.%s, self.%s)" % (cb.pos_name[2], cb.vel_name[2]))
        elif len(cb.vel_name) == 0:
            if len(cb.pos_name) == 1:
                exec("self.coeffX = linear_fit(self.time, self.%s[:, 0])" % cb.pos_name[0])
                exec("self.coeffY = linear_fit(self.time, self.%s[:, 1])" % cb.pos_name[0])
                exec("self.coeffZ = linear_fit(self.time, self.%s[:, 2])" % cb.pos_name[0])
            elif len(cb.pos_name) == 3:
                exec("self.coeffX = linear_fit(self.time, self.%s)" % cb.pos_name[0])
                exec("self.coeffY = linear_fit(self.time, self.%s)" % cb.pos_name[1])
                exec("self.coeffZ = linear_fit(self.time, self.%s)" % cb.pos_name[2])
        else:
            print('Please check pos_name and vel_name!')
            sys.exit(0)


def linear_fit(times, pos):
    paras = []
    for loop1 in range(len(times) - 1):
        x0 = pos[loop1] * cb.pos_unit
        x1 = pos[loop1 + 1] * cb.pos_unit
        t0 = times[loop1]
        t1 = times[loop1 + 1]
        delta_x = x1 - x0
        if delta_x > boxsize / 2:
            x1 -= boxsize
        elif delta_x < - boxsize / 2:
            x1 += boxsize
        paras.append([0, 0, (x1 - x0) / (t1 - t0), x0 - (x1 - x0) * t0 / (t1 - t0)])
    return np.array(paras)


def fit_coeff(times, pos, vel):  # unit of times: Gyr/h, unit of pos: Mpc/h, unit of vel: km/s
    paras = []
    for loop1 in range(len(times) - 1):
        x0 = pos[loop1] * cb.pos_unit
        x1 = pos[loop1 + 1] * cb.pos_unit
        t0 = times[loop1]
        t1 = times[loop1 + 1]
        delta_x = x1 - x0
        if delta_x > boxsize / 2:
            x1 -= boxsize
        elif delta_x < - boxsize / 2:
            x1 += boxsize
        paras.append(solve_coefficients(t0, t1, x0, x1, vel[loop1] * cb.vel_unit, vel[loop1 + 1] * cb.vel_unit))
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
    c = 299792458
    pc = 3.08568025e16
    t_H = 9785641806 * 1e-9  # Hubble time in Gyr/h
    # 距离，坐标统一单位Mpc/h，时间统一单位1e9*year

    user_path = os.path.abspath(os.path.dirname(__file__)) + '/../user_generated/'
    with open(user_path + '/cb.pkl', 'rb') as f:
        cb = pickle.loads(f.read())
    f.close()

    start = ti.time()
    print('Program start at : ', ti.ctime())

    loadpath = cb.data_path
    save_path = cb.output_path
    interpo_method = cb.interpo  ###
    full_sky = cb.full_sky
    ra_min, ra_max, dec_min, dec_max = cb.ra_min, cb.ra_max, cb.dec_min, cb.dec_max
    obs = cb.observer_pos
    r_max = cb.r_max
    z_max = cb.z_max
    N = int(cb.N)
    O_M = cb.O_M
    O_L = cb.O_L
    O_K = cb.O_K
    boxsize = cb.boxsize  # boxsize in Mpc/h
    ssn = int(cb.snap_tot)
    cores = int(cb.cores)
    boxes = cb.boxes
    fields = cb.props
    t_cos = z2t(3000)  # age of universe

    boxes_len = len(boxes)
    print('Z_max =', z_max, 'N =', N, 'Box count =', boxes_len)

    # define output path
    if not cb.output_path:
        path = os.path.abspath(os.path.dirname(__file__))
        mkpath = path + '/../BLiC_out'
        if not Path(mkpath).is_dir():
            Path(mkpath).mkdir()
        cb.output_path = path + '/../BLiC_out/'
    if cb.max_branch:
        # parallel computing by file branches
        ibranch = 2661  ###
        if Path(cb.output_path).is_dir():
            cb.output_path += 'Mockcatalog_%d.hdf5' % ibranch
        else:
            cb.output_path += '.%d.hdf5' % ibranch
        for i in range(0, ssn):
            isnap = i
            exec("iname = %s" % cb.file_name_with_branchs)
            if Path(loadpath + iname).is_file():
                ssn_start = i
                break
            elif i == ssn - 1:
                print('No input files are found, please check the path and names !!!')
    else:
        if Path(cb.output_path).is_dir():
            cb.output_path += 'mockcatalog.hdf5'
        else:
            cb.output_path += '.hdf5'

        ## 20230212
        for i in range(0, ssn):
            isnap = i
            exec("iname = %s" % cb.file_name)
            if Path(loadpath + iname).is_file():
                ssn_start = i
                break
            elif i == ssn - 1:
                print('No input files are found, please check the path and names !!!')

    '''最终输出文件的列名'''
    fields = cb.props
    output_list = cb.output_list
    zcp = cb.zcp

    ###
    # zcp[69] = zcp[68]
    # zcp[76] = zcp[75]
    ###

    ssn_end = int(max(min(boxes[:, 4]), ssn_start))  # snapshot cut
    ssn_st = int(min(max(boxes[:, 3]), ssn))  ### 20230202 新增
    # 空列表准备读数据
    rcp = np.zeros(ssn)
    tcp = np.zeros(ssn)
    pos_name = cb.pos_name
    vel_name = cb.vel_name

    data = [[] for row in range(ssn)]
    # 读、存数据都是倒序的，变量最后存的是ssn最大的内容
    # 如果读数据读爆内存，可以在每一个data中删除一些不用的数据
    coefs_t2r = [[] for row2 in range(ssn)]
    coefs_r2z = [[] for row2 in range(ssn)]
    colname = []
    for i in range(ssn_end, ssn):
        rcp[i] = z2r(zcp[i])
        tcp[i] = z2t(zcp[i])
        ###
        # if i == 69 or i == 76:
        #     continue
        ###
        rr = np.array([z2r(i2) for i2 in np.arange(zcp[i], zcp[i - 1], 1e-4)])  # z是在增加的
        tt = np.array([z2t(i2) for i2 in np.arange(zcp[i], zcp[i - 1], 1e-4)])
        zz = np.array([i2 for i2 in np.arange(zcp[i], zcp[i - 1], 1e-4)])
        tt = t_cos - tt  # 转换为宇宙年龄
        coefs_t2r[i] = np.polyfit(tt, rr, 3)  # 宇宙年龄-共动距离的关系拟合，第63个存的snapshot63到62的拟合公式
        coefs_r2z[i] = np.polyfit(rr, zz, 3)

        isnap = i

        if cb.max_branch:
            exec("iname = %s" % cb.file_name_with_branchs)
        else:
            exec("iname = %s" % cb.file_name)

        if '.hdf5' in os.path.splitext(cb.file_name)[-1] and len(cb.hdf_database):  # for .hdf5 with dataset
            data[i] = h5py.File(cb.data_path + '/' + iname, "r")[cb.hdf_database][:]
        else:
            data[i] = pd.read_table(cb.data_path + '/' + iname, header=0).to_records(index=False)  # for other tables
            colname = list(data[i].dtype.names)
            # if not colname[0] == str(colname[0]):
            #     data[i] = pd.read_table(cb.data_path + iname, header=None).to_records(index=False)
            #     colname = np.linspace(0, len(data[i]) - 1, len(data[i]), dtype='int')  # nums as col names

        # trackIDs[i] = np.array(data[i]['TrackId'])
        # born[i] = np.array(data[i]['SnapshotIndexOfBirth'])
        print('len(data[%d]): %d' % (i, len(data[i])))
    print(zcp, rcp, tcp)

    ga = []
    for k in range(ssn_end, ssn):  ### ssn_st
        # if k == 76 or k == 69:
        #     continue
        if k == ssn_end:
            need_read = np.linspace(0, len(data[k]) - 1, len(data[k]), dtype='int')
        else:
            need_read = np.linspace(len(data[k - 1]), len(data[k]) - 1, len(data[k]) - len(data[k - 1]), dtype='int')
        Ng = len(ga)  # N galaxies

        if len(need_read):
            print('ssn=', k, 'need read=', len(need_read))
            pool = Pool(processes=cb.cores)
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
    para = para.reshape(boxes_len * ga_len, 6 + 3)  # 生成一批6 + 3(random tiling)参数数组， 然后并行处理
    print('para = para.reshape(boxes_len * ga_len, 6), used time: %.1f s' % (ti.time() - time_now))

    pool = Pool(processes=cores)
    res = pool.map(cone, para)  # 只传递id就可以了，但是要传长度最小的，否则调用会出错  ### np.array
    pool.close()
    pool.join()
    print('res = pool.map(cone, para), used time: %.1f s' % (ti.time() - time_now))

    gal_tot = [i2 for i2 in res if i2 is not None]  # res = list(filter(None, res))

    print('gal_tot = [i2 for i2 in res if i2 is not None], used time: %.1f s' % (ti.time() - time_now), len(gal_tot))

    gal_tot = DataFrame(gal_tot, columns=output_list)  # , dtype=object)
    gal_tot = pd.DataFrame(gal_tot).to_records(index=None)
    print('gal_tot = pd.DataFrame(gal_tot).to_records(index=None), used time: %.1f s' % (ti.time() - time_now))

    op = cb.output_path.replace('.hdf5', '.npy')
    np.save(op, gal_tot)
    # with h5py.File(cb.output_path, "w") as output_file:
    #     if cb.dtypes:
    #         output_file.create_dataset('Subhalos', len(gal_tot),
    #                                    dtype=cb.dtypes)  ### 有问题
    #     else:
    #         output_file.create_dataset('Subhalos', gal_tot.shape)  ###len(gal_tot)
    #     output_file['Subhalos'][:] = gal_tot
    # output_file.close()
    print('Len: ', len(gal_tot), 'Total used: %.1f s' % (ti.time() - start))
    print('output_docu=', save_path)
