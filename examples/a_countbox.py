import pandas as pd
import numpy as np
import time as ti
from scipy.integrate import quad
import math
import sys


def z2r(z):
    func2 = lambda y: 1 / (O_M * ((1 + y) ** 3) + O_K * ((1 + y) ** 2) + O_L) ** 0.5
    v2 = quad(func2, 0, z)
    r_c = c * v2[0] / 1e5  # comoving distance in Mpc/h
    return r_c


def fss(begin, end):  # 找盒子位于哪几个snapshot之间,返回两个参数
    ben = 0
    enn = 0
    for l1 in range(snap_end + 1, snap_tot):
        if begin == 0:
            ben = snap_tot - 1
            break
        elif rcp[l1] <= begin:
            ben = l1  # begin snapshot number
            break
    for l2 in range(snap_tot - 1, snap_end, -1):
        if rcp[l2] >= end:
            enn = l2
            break  # 总是先找到begin，才能用
    return ben, enn


def judge_angle_2pi(v0, v1):  # judge two vectors' angle between 0 and 2*pi
    ang = math.atan2(np.linalg.norm(np.cross(v0, v1)), np.dot(np.transpose(v0), v1))
    # if these two lines used, output will between 0 and 2*pi
    if np.cross(v0, v1)[2] < 0:
        ang = 2 * np.pi - ang
    return ang  # return in rad. not in degree


def judge_angle_pi(v0, v1):  # judge two vectors' angle between 0 and pi
    ang = math.atan2(np.linalg.norm(np.cross(v0, v1)), np.dot(np.transpose(v0), v1))
    return ang  # return in rad. not in degree


def judge_box(angles):  # in 2 * pi
    wrong_points = 0
    for ii in range(len(angles)):
        if 0 <= angles[ii] < np.pi / 2 or 3 * np.pi / 2 <= angles[ii] <= 2 * np.pi:
            wrong_points += 0
        else:
            wrong_points += 1
    if 0 <= min(angles) < np.pi / 2 and 3 * np.pi / 2 <= max(angles) <= 2 * np.pi:
        wrong_points += 0
    else:
        wrong_points += 1

    if wrong_points:
        return False
    else:
        return True


def car_to_RA(x, y, z):
    dec = math.pi / 2 - math.atan2(math.sqrt(x ** 2 + y ** 2), z)
    RA = math.atan2(y, x)
    if RA < 0:
        RA += math.pi * 2
    return RA * 180 / np.pi, dec * 180 / np.pi


def count_box(cone_geo, ra_min=-1, ra_max=-1, dec_min=-1, dec_max=-1, hfa=-1, los=-1):
    boxes = []
    if cone_geo == 0 and hfa < np.pi:
        for i in range(len(check_points)):  # check every box
            check_now = boxpoints + check_points[i] - obs_unit
            angles_2pi = np.zeros(8)
            angles_pi = np.zeros(8)
            for j in range(len(check_now)):
                angles_2pi[j] = judge_angle_2pi(los, check_now[j])
                angles_pi[j] = judge_angle_pi(los, check_now[j])
            norm = np.linalg.norm(check_now, axis=1)
            if min(angles_pi) < hfa:  # 某一个角在物理锥内
                if min(norm) < r_max / boxsize:
                    ob = min(norm) * boxsize  # 盒子最近处
                    oe = max(norm) * boxsize  # 最远
                    a, b = fss(ob, oe)
                    boxes.append([check_points[i][2], check_points[i][1], check_points[i][0], a, b])  # Z-Y-X
            elif judge_box(angles_2pi):
                if min(norm) < r_max / boxsize:
                    ob = min(norm) * boxsize  # 盒子最近处
                    oe = max(norm) * boxsize  # 最远
                    a, b = fss(ob, oe)
                    boxes.append([check_points[i][2], check_points[i][1], check_points[i][0], a, b])
            elif np.linalg.norm(check_points[i]) == 0:  # 000处盒子一定会经过
                ob = 0  # 盒子最近处
                oe = max(norm) * boxsize  # 最远
                a, b = fss(ob, oe)
                boxes.append([check_points[i][2], check_points[i][1], check_points[i][0], a, b])

    elif cone_geo == 1 and 0 <= ra_min < ra_max < 360 and -90 <= dec_min < dec_max <= 90:
        for i in range(len(check_points)):  # check every box
            check_now = boxpoints + check_points[i] - obs_unit
            if np.linalg.norm(check_points[i]) == 0:  # 000处盒子一定会经过
                norm = np.linalg.norm(check_now, axis=1)
                ob = 0  # 盒子最近处
                oe = max(norm) * boxsize  # 最远
                a, b = fss(ob, oe)
                boxes.append([check_points[i][2], check_points[i][1], check_points[i][0], a, b])
                continue

            ra_now = np.zeros(8)
            dec_now = np.zeros(8)
            for j in range(8):
                ra_now[j], dec_now[j] = car_to_RA(check_now[j][0], check_now[j][1], check_now[j][2])

                if ra_min < ra_now[j] < ra_max and dec_min < dec_now[j] < dec_max:
                    # 如果某一个顶点在天区内
                    norm = np.linalg.norm(check_now, axis=1)
                    if min(norm) < r_max / boxsize:
                        ob = min(norm) * boxsize  # 盒子最近处
                        oe = max(norm) * boxsize  # 最远
                        a, b = fss(ob, oe)
                        boxes.append([check_points[i][2], check_points[i][1], check_points[i][0], a, b])
                        break

            if max(ra_now) > ra_max and min(ra_now) < ra_min and max(dec_now) > dec_max and min(dec_now) < dec_min:
                # 如果盒子横跨整个天区范围
                norm = np.linalg.norm(check_now, axis=1)
                if min(norm) < r_max / boxsize:
                    ob = min(norm) * boxsize  # 盒子最近处
                    oe = max(norm) * boxsize  # 最远
                    a, b = fss(ob, oe)
                    boxes.append([check_points[i][2], check_points[i][1], check_points[i][0], a, b])
                    continue

    elif cone_geo == 2:  # full-sky map
        for i in range(len(check_points)):  # check every box
            check_now = boxpoints + check_points[i] - obs_unit

            # 对于在坐标轴上的盒子，到盒子的最小距离应该是到某一平面的距离而不是到端点！
            # 共6个面
            if check_points[i][0] != 0 and check_points[i][1] == 0 and check_points[i][2] == 0:
                check_now = np.vstack((check_now, np.array([check_points[i][0] - obs_unit[0], 0, 0])))
                check_now = np.vstack((check_now, np.array([check_points[i][0] + 1 - obs_unit[0], 0, 0])))

            elif check_points[i][0] == 0 and check_points[i][1] != 0 and check_points[i][2] == 0:
                check_now = np.vstack((check_now, np.array([0, check_points[i][1] - obs_unit[1], 0])))
                check_now = np.vstack((check_now, np.array([0, check_points[i][1] + 1 - obs_unit[1], 0])))
            elif check_points[i][0] == 0 and check_points[i][1] == 0 and check_points[i][2] != 0:
                check_now = np.vstack((check_now, np.array([0, 0, check_points[i][2] - obs_unit[2]])))
                check_now = np.vstack((check_now, np.array([0, 0, check_points[i][2] + 1 - obs_unit[2]])))

            # 对于在坐标轴平面上的盒子，到盒子的最小距离应该是到某一棱的距离而不是到端点！
            # 共12条棱
            for jj in range(3):
                index_jj = np.delete(np.array([0, 1, 2]), jj)
                if check_points[i][jj] == 0 and check_points[i][index_jj[0]] * check_points[i][index_jj[1]] > 0:
                    need_append = np.zeros(3)
                    need_append[jj] = 0
                    need_append[index_jj[0]] = check_points[i][index_jj[0]] - obs_unit[index_jj[0]]
                    need_append[index_jj[1]] = check_points[i][index_jj[1]] - obs_unit[index_jj[1]]
                    check_now = np.vstack((check_now, need_append))
                    need_append = np.zeros(3)
                    need_append[jj] = 0
                    need_append[index_jj[0]] = check_points[i][index_jj[0]] - obs_unit[index_jj[0]] + 1
                    need_append[index_jj[1]] = check_points[i][index_jj[1]] - obs_unit[index_jj[1]] + 1
                    check_now = np.vstack((check_now, need_append))
                elif check_points[i][jj] == 0 and check_points[i][index_jj[0]] * check_points[i][index_jj[1]] < 0:
                    need_append = np.zeros(3)
                    need_append[jj] = 0
                    need_append[index_jj[0]] = check_points[i][index_jj[0]] - obs_unit[index_jj[0]] + 1
                    need_append[index_jj[1]] = check_points[i][index_jj[1]] - obs_unit[index_jj[1]]
                    check_now = np.vstack((check_now, need_append))
                    need_append = np.zeros(3)
                    need_append[jj] = 0
                    need_append[index_jj[0]] = check_points[i][index_jj[0]] - obs_unit[index_jj[0]]
                    need_append[index_jj[1]] = check_points[i][index_jj[1]] - obs_unit[index_jj[1]] + 1
                    check_now = np.vstack((check_now, need_append))

            norm = np.linalg.norm(check_now, axis=1)
            if np.linalg.norm(check_points[i]) == 0:  # 000处盒子一定会经过
                ob = 0  # 盒子最近处
                oe = max(norm) * boxsize  # 最远
                a, b = fss(ob, oe)
                # print([check_points[i][2], check_points[i][1], check_points[i][0], a, b])
                boxes.append([check_points[i][2], check_points[i][1], check_points[i][0], a, b])
            elif min(norm) < r_max / boxsize:
                # 只要有一个最近都角在球内
                ob = min(norm) * boxsize  # 盒子最近处
                oe = max(norm) * boxsize  # 最远
                a, b = fss(ob, oe)
                boxes.append([check_points[i][2], check_points[i][1], check_points[i][0], a, b])

    else:
        print("Wrong value(s): Please check cone_geo, hfa, R.A. or dec.")
        sys.exit()

    return np.array(boxes)


if __name__ == '__main__':
    '''
    check which box will sit in light-cone and skip the outside boxes
    '''

    start = ti.time()
    c = 299792458
    pc = 3.08568025e16
    t_H = 9785641806 * 1e-9  # Hubble time in Gyr/h

    # 模拟盒子信息，或许可以从文件读出来
    cos_paras = np.array([0.3111, 0.6889, 0])
    boxsize = 1000  # Mpc/h
    snap_tot = 128
    snap_end = 30

    O_M = cos_paras[0]  ###
    O_L = cos_paras[1]
    O_K = cos_paras[2]

    z_max = 2  # 输入，截断红移 ###
    obs = np.array([500, 500, 500]) % boxsize  # 不应该超出第一个盒子的范围

    Mi_redshifts = pd.read_table('/home/zltan/9t/redshift.csv', header=0)
    zcp = Mi_redshifts['z']

    r_max = z2r(z_max)
    if np.linalg.norm(obs):
        obs_unit = obs / boxsize  # np.linalg.norm(obs)
    else:
        obs = np.array([1e-5, 1e-5, 1e-5])
        obs_unit = obs
    N = int(r_max / boxsize) + 2  # +1+1是因为观测者可以在第一个盒子内任意移动
    rcp = np.zeros(snap_tot)
    for i in range(len(zcp)):
        rcp[i] = z2r(zcp[i])

    check_points = []
    tiles = N
    for i in np.arange(-tiles, tiles + 1):
        for j in np.arange(-tiles, tiles + 1):
            for k in np.arange(-tiles, tiles + 1):
                check_points.append([i, j, k])
    check_points = np.array(check_points)

    boxpoints = []
    for i in np.arange(2):
        for j in np.arange(2):
            for k in np.arange(2):
                boxpoints.append([i, j, k])

    # cone_geo = 0  # 0 for lfs+hfa judge, 1 for R.A.+dec judge, 2 for full sky
    # counted = count_box(0, los=np.array([1, 1, 1]), hfa=np.pi / 3)

    # x axis direction means RA=0 and dec=0
    # RA ranges from 0 to 360 (deg), dec. ranges from -90 to 90 (deg)
    # counted = count_box(1, ra_min=50, ra_max=55, dec_min=50, dec_max=55)
    cone_geo = 2
    counted = count_box(cone_geo)
    ncount = len(counted)
    np.save('boxcounts_%d_%d_pacth.npy' % (cone_geo, int(z_max)), counted)  ###
    print('Input parameters are:\n', 'Omega_m = %f\n' % O_M, 'Omega_lambda = %f\n' % O_L, 'Omega_k = %f\n' % O_K,
          'Boxsize = %f\n' % boxsize)
    print('Redshift cut = %f\n' % z_max, 'cone_geo = %d\n' % cone_geo, 'Comoving distance cut = %.2f Mpc/h\n' % r_max,
          'box count = %d\n' % ncount)
    #     print('min snapshot num: ', min(boxcount[:, 4]))
    print('time used:', ti.time() - start)
