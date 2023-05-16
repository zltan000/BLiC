# this file calculates the boxes and snapshots needed for the whole light-cone
# check which box will sit in light-cone and skip the outside boxes
import pandas as pd
import numpy as np
import time as ti
from scipy.integrate import quad
import math
import sys
import os
import h5py
from random import random
import pickle


class CalcBoxes:
    def __init__(self):
        init_names = ['data_path', 'file_name', 'pos_name', 'vel_name', 'props', 'redshift_table_path', 'ran_tiling',
                      'redshift_col_name', 'O_M', 'O_L', 'O_K', 'boxsize', 'z_max', 'observer_pos', 'full_sky',
                      'ra_min', 'ra_max', 'dec_min', 'dec_max', 'cores', 'max_branch', 'file_name_with_branchs',
                      'max_running', 'c', 'pc', 't_H', 'output_path']
        for loop0 in range(len(init_names)):
            exec("self.%s = %s" % (init_names[loop0], init_names[loop0]))

        data = pd.read_table(self.redshift_table_path, header=0, sep=',')  # may need to edit
        # get z check points
        self.zcp = []
        if self.redshift_col_name == str(self.redshift_col_name):
            self.zcp = np.array(data[self.redshift_col_name])
        elif isinstance(self.redshift_col_name, int):
            self.zcp = np.array(data.iloc[:, int(self.redshift_col_name)])
        if not len(self.zcp):
            print('ERROR: no redshift list read! check the redshift_table_path and the redshift_col_name value')
        self.snap_tot = len(self.zcp)

        # get column names
        self.output_list = ['RA', 'Dec', 'd_comoving', 'z_cos', 'z_obs', 'v_los', 'snapnum']
        self.hdf_database = 'Subhalos'
        self.dtypes = []
        if self.props[0] == 'all':  # all properties are needed
            isnap = self.snap_tot - 1
            if self.max_branch:
                ibranch = self.max_branch - 1
                exec("iname = %s" % self.file_name_with_branchs)
            else:
                exec("iname = %s" % self.file_name)

            if os.path.splitext(self.file_name)[-1] == '.hdf5' and len(self.hdf_database):  # for .hdf5 with dataset
                data = h5py.File(self.data_path + iname, "r")[self.hdf_database][:]
                self.props = list(data.dtype.names)
                self.dtypes = data.dtype
                self.output_list += self.props
            else:
                data = pd.read_table(self.data_path + iname, header=0)  # for other tables
                colname = list(data.columns)
                if colname[0] == str(colname[0]):
                    self.props = colname
                else:
                    data = pd.read_table(self.data_path + iname, header=None)
                    colname = np.linspace(0, len(data.iloc[0]) - 1, len(data.iloc[0]), dtype='int')  # nums as col names
                    self.props = list(colname)
                self.output_list += self.props

        elif self.props[0] == str(self.props[0]):
            self.output_list += self.props
        else:
            for i in range(len(self.props)):
                self.output_list.append(str(self.props[i]))
            # self.output_list += ['RA', 'Dec', 'd_comoving', 'z_cos', 'z_obs', 'v_los']

        self.obs = np.array(self.observer_pos) % self.boxsize  # 不应该超出第一个盒子的范围
        for loop1 in range(3):
            if self.obs[loop1] == 0:
                self.obs[loop1] == 1e-6
        self.obs_unit = self.obs / self.boxsize
        self.r_max = self.z2r(self.z_max)
        self.N = self.r_max // self.boxsize + 2  # +1+1是因为观测者可以在第一个盒子内任意移动
        self.rcp = np.zeros(self.snap_tot)
        for loop2 in range(len(self.zcp)):
            self.rcp[loop2] = self.z2r(self.zcp[loop2])
        self.snap_end = 0
        self.interpo = np.ones(len(self.props), dtype='bool')
        self.boxes = self.count_box()
        self.pos_unit = 1
        self.vel_unit = 1

    def count_box(self):
        check_points = []
        for loop3 in np.arange(-self.N, self.N + 1):
            for loop4 in np.arange(-self.N, self.N + 1):
                for loop5 in np.arange(-self.N, self.N + 1):
                    check_points.append([loop3, loop4, loop5])
        check_points = np.array(check_points)

        boxpoints = []
        for loop3 in np.arange(2):
            for loop4 in np.arange(2):
                for loop5 in np.arange(2):
                    boxpoints.append([loop3, loop4, loop5])
        boxes = []

        # calculate full-sky map first, it's precise
        for i in range(len(check_points)):  # check every box
            check_now = boxpoints + check_points[i] - self.obs_unit
            # 对于在坐标轴上的盒子，到盒子的最小距离应该是到某一平面的距离而不是到端点！
            # 共6个面
            if check_points[i][0] != 0 and check_points[i][1] == 0 and check_points[i][2] == 0:
                check_now = np.vstack((check_now, np.array([check_points[i][0] - self.obs_unit[0], 0, 0])))
                check_now = np.vstack((check_now, np.array([check_points[i][0] + 1 - self.obs_unit[0], 0, 0])))
            elif check_points[i][0] == 0 and check_points[i][1] != 0 and check_points[i][2] == 0:
                check_now = np.vstack((check_now, np.array([0, check_points[i][1] - self.obs_unit[1], 0])))
                check_now = np.vstack((check_now, np.array([0, check_points[i][1] + 1 - self.obs_unit[1], 0])))
            elif check_points[i][0] == 0 and check_points[i][1] == 0 and check_points[i][2] != 0:
                check_now = np.vstack((check_now, np.array([0, 0, check_points[i][2] - self.obs_unit[2]])))
                check_now = np.vstack((check_now, np.array([0, 0, check_points[i][2] + 1 - self.obs_unit[2]])))
            # 对于在坐标轴平面上的盒子，到盒子的最小距离应该是到某一棱的距离而不是到端点！
            # 共12条棱
            for jj in range(3):
                index_jj = np.delete(np.array([0, 1, 2]), jj)
                if check_points[i][jj] == 0 and check_points[i][index_jj[0]] * check_points[i][index_jj[1]] > 0:
                    need_append = np.zeros(3)
                    need_append[jj] = 0
                    need_append[index_jj[0]] = check_points[i][index_jj[0]] - self.obs_unit[index_jj[0]]
                    need_append[index_jj[1]] = check_points[i][index_jj[1]] - self.obs_unit[index_jj[1]]
                    check_now = np.vstack((check_now, need_append))
                    need_append = np.zeros(3)
                    need_append[jj] = 0
                    need_append[index_jj[0]] = check_points[i][index_jj[0]] - self.obs_unit[index_jj[0]] + 1
                    need_append[index_jj[1]] = check_points[i][index_jj[1]] - self.obs_unit[index_jj[1]] + 1
                    check_now = np.vstack((check_now, need_append))
                elif check_points[i][jj] == 0 and check_points[i][index_jj[0]] * check_points[i][index_jj[1]] < 0:
                    need_append = np.zeros(3)
                    need_append[jj] = 0
                    need_append[index_jj[0]] = check_points[i][index_jj[0]] - self.obs_unit[index_jj[0]] + 1
                    need_append[index_jj[1]] = check_points[i][index_jj[1]] - self.obs_unit[index_jj[1]]
                    check_now = np.vstack((check_now, need_append))
                    need_append = np.zeros(3)
                    need_append[jj] = 0
                    need_append[index_jj[0]] = check_points[i][index_jj[0]] - self.obs_unit[index_jj[0]]
                    need_append[index_jj[1]] = check_points[i][index_jj[1]] - self.obs_unit[index_jj[1]] + 1
                    check_now = np.vstack((check_now, need_append))

            norm = np.linalg.norm(check_now, axis=1)
            if np.linalg.norm(check_points[i]) == 0:  # 000处盒子一定会经过
                ob = 0  # 盒子最近处
                oe = max(norm) * self.boxsize  # 最远
                a, b = self.fss(ob, oe)
                # print([check_points[i][2], check_points[i][1], check_points[i][0], a, b])
                boxes.append([check_points[i][2], check_points[i][1], check_points[i][0], a, b])
            elif min(norm) < self.r_max / self.boxsize:
                # 只要有一个最近都角在球内
                ob = min(norm) * self.boxsize  # 盒子最近处
                oe = max(norm) * self.boxsize  # 最远
                a, b = self.fss(ob, oe)
                boxes.append([check_points[i][2], check_points[i][1], check_points[i][0], a, b])

        # if self.cone_geo == 0 and self.hfa < np.pi:
        #     need_del = []
        #     for loop6 in range(len(check_points)):
        #         check_now = boxpoints + check_points[loop6] - self.obs_unit
        #         angles_2pi = np.zeros(8)
        #         angles_pi = np.zeros(8)
        #         for loop7 in range(len(check_now)):
        #             angles_2pi[loop7] = judge_angle_2pi(self.los, check_now[loop7])
        #             angles_pi[loop7] = judge_angle_pi(self.los, check_now[loop7])
        #         norm = np.linalg.norm(check_now, axis=1)
        #         if min(angles_pi) > self.hfa or min(norm) > self.r_max / self.boxsize:  # 某一个角在物理锥内
        #             need_del.append(loop6)
        #         elif not judge_box(angles_2pi):
        #             need_del.append(loop6)
        #     boxes = np.delete(np.array(boxes), need_del, axis=1)

        if (not self.full_sky) and 0 <= self.ra_min < self.ra_max < 360 and -90 <= self.dec_min < self.dec_max <= 90:
            need_del = []
            for i in range(len(boxes)):  # check every box
                check_now = boxpoints + np.array([boxes[i][2], boxes[i][1], boxes[i][0]]) - self.obs_unit
                ra_now = np.zeros(8)
                dec_now = np.zeros(8)
                judge_del = 0
                for j in range(8):
                    ra_now[j], dec_now[j] = car_to_RA(check_now[j][0], check_now[j][1], check_now[j][2])
                    if self.ra_min < ra_now[j] < self.ra_max and self.dec_min < dec_now[j] < self.dec_max:
                        # 如果某一个顶点在天区内
                        norm = np.linalg.norm(check_now, axis=1)
                        if min(norm) < self.r_max / self.boxsize:
                            judge_del += 1
                            break

                ra_now_max, ra_now_min = max(ra_now), min(ra_now)
                dec_now_max, dec_now_min = max(dec_now), min(dec_now)
                if 270 < ra_now_max < 360 and 0 < ra_now_min < 90:
                    for j in range(8):
                        if 270 < ra_now[j] < 360:
                            ra_now[j] -= 360
                ra_now_max, ra_now_min = max(ra_now), min(ra_now)
                if ra_now_max > self.ra_max and ra_now_min < self.ra_min and (dec_now_min < self.dec_max < dec_now_max
                                                                              or dec_now_min < self.dec_min < dec_now_max or self.dec_min < dec_now_min < dec_now_max < self.dec_max):
                    judge_del += 1
                elif dec_now_max > self.dec_max and dec_now_min < self.dec_min and (
                        ra_now_min < self.ra_max < ra_now_max
                        or ra_now_min < self.ra_min < ra_now_max or self.ra_min < ra_now_min < ra_now_max < self.ra_max):
                    judge_del += 1

                if not judge_del:
                    need_del.append(i)
            boxes = np.delete(np.array(boxes), need_del, axis=0)  #, axis=-1 20220224

        # elif self.cone_geo != 2:
        #     print("Wrong value(s): Please check cone_geo, hfa, R.A. or dec.")
        #     sys.exit()
        tiling_arr = [np.zeros(3) for i in range(len(boxes))]
        if self.ran_tiling:  # random tiling when tiling boxes
            for i in range(len(boxes)):
                for j in range(3):
                    tiling_arr[i][j] = random()
        boxes = np.concatenate((boxes, tiling_arr), axis=1)

        return np.array(boxes)

    def fss(self, begin, end):  # 找盒子位于哪几个snapshot之间,返回两个参数
        ben = 0
        enn = 0
        for l1 in range(self.snap_end + 1, self.snap_tot):
            if begin == 0:
                ben = self.snap_tot - 1
                break
            elif self.rcp[l1] <= begin:
                ben = l1  # begin snapshot number
                break
        for l2 in range(self.snap_tot - 1, self.snap_end, -1):
            if self.rcp[l2] >= end:
                enn = l2
                break  # 总是先找到begin，才能用
        return ben, enn


    def z2r(self, z):
        func2 = lambda y: 1 / (self.O_M * ((1 + y) ** 3) + self.O_K * ((1 + y) ** 2) + self.O_L) ** 0.5
        v2 = quad(func2, 0, z)
        r_c = c * v2[0] / 1e5  # comoving distance in Mpc/h
        return r_c



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

    # if wrong_points:
    #     return False
    # else:
    #     return True

    # too complicated, all return True
    return True


def car_to_RA(x, y, z):
    dec = math.pi / 2 - math.atan2(math.sqrt(x ** 2 + y ** 2), z)
    RA = math.atan2(y, x)
    if RA < 0:
        RA += math.pi * 2
    return RA * 180 / np.pi, dec * 180 / np.pi


# beginning of .ini
# All parameters (file paths) are input in this file by users, please fill in carefully
# They are directly input to python files

# (absolute) snapshot folder path
data_path = "/home/cossim/IllustrisTNG/Illustris-3/subcat/"
# the quoted name formats of snapshots with variable named 'isnap', which means number i file of snapshots
# for example: "'snapshot_%d.hdf5' % isnap"
# the input file is supposed to have column names!
file_name = "'SubSnap_%03d.hdf5' % isnap"

# the column name(s) of positions
# if only one column is given, default pos_name[0], pos_name[1], pos_name[2] refer to x, y, z position
# should be in forms like: pos_name = ['a', 'b', 'c']; pos_name = ['pos']
pos_name = ['ComovingAveragePosition']
# the column name(s) of velocities, formats as same as pos_name, the length should be equal to pos_name.
# could be empty if no vel. information
vel_name = ['PhysicalAverageVelocity']
# properties needed in lightcones, including pos. and vel. These columns will be read from input files and otherwise will not
# the number of properties will also influence the allocated memory when programs are running
# if all properties are used, input: props = ['all']
props = ['TrackId', 'Nbound', 'Mbound', 'HostHaloId', 'Rank', 'BoundR200CritComoving','BoundM200Crit', 'ComovingAveragePosition', 'PhysicalAverageVelocity']

# (absolute) redshift table file path
redshift_table_path = "/home/zltan/Illustris-3/i3redshift.csv"
# the column name of redshift in table, like 'z'; it can also be a number
# NOTICE: the length of redshift column will be treated as number of snapshots in simulation
# the redshift would be decreasing while the snapshot num is increasing in this table
redshift_col_name = "z"

# cosmic parameters
O_M = 0.2726  # omega matter
O_L = 0.7274  # omega lambda
O_K = 0  # omega curvature

# the length of side of simulation box, in (Mpc/h) unit
boxsize = 75
# the redshift cut in lightcones
z_max = 1.2
# the position of the observer in original box, should be in 3-dimensions
observer_pos = [ 5, 5, 5]
# random tiling switch, True or False
ran_tiling = False
# 'True' or 'False' to make a full-sky map, if False, ra_min, ra_max, dec_min, dec_max are needed
full_sky = False
# unit in degrees, not needed to change if full_sky = True
ra_min = 42  # min 0
ra_max = 45  # max 360
dec_min = 42  # min -90
dec_max = 45  # max 90

# parallel options: the program will be paralleled using multiprocessing package
# the max CPU cores will be used (only in one node for a server), multi-nodes' tasks are not permitted
cores = 72

# below are OPTIONAL parameters

# output file path of lightcone mock catalogues, all outputs will be in hdf5 formats, you can give a dir end with '/', or dir + filename end without '.hdf5', default path is at '/blic_out'
output_path = ""
# for huge simulations, one snapshot may be stored in different files, in this case, you are supposed to
# input the max branch number and the new file forms with variables named 'isnap' and 'ibranch'
# in ordinary case, max_branch should be the number of files of the last snapshot
# NOTICE: objects at same 'ibranch' but in lower redshift files should all be included in higher redshift files
max_branch = 0
file_name_with_branchs = ""
# this program can be paralleled in different snapshot branches (in different nodes for a server), input the max number
# of tasks you may hope to run at same time
max_running = 1

# end of .ini

c = 299792458
pc = 3.08568025e16
t_H = 9785641806 * 1e-9  # Hubble time in Gyr/h

if __name__ == '__main__':
    start = ti.time()
    cb = CalcBoxes()
    
    out_path = os.path.abspath(os.path.dirname(__file__))
    ti.sleep(5)
    save_class = open(out_path + '/cb.pkl', 'wb')
    strs = pickle.dumps(cb)
    save_class.write(strs)
    save_class.close()

    counted = cb.boxes
    ncount = len(counted)
    np.save(out_path + '/boxcounts.npy', counted)  ###
    print('Input parameters are:\n', 'Omega_m = %f\n' % cb.O_M, 'Omega_lambda = %f\n' % cb.O_L,
          'Omega_k = %f\n' % cb.O_K,
          'Boxsize = %f\n' % cb.boxsize)
    print('Redshift cut = %f\n' % cb.z_max, 'full_sky = %d\n' % cb.full_sky,
          'Comoving distance cut = %.2f Mpc/h\n' % cb.r_max,
          'box count = %d\n' % ncount)
    #     print('min snapshot num: ', min(boxcount[:, 4]))
    print('time used:', ti.time() - start)
