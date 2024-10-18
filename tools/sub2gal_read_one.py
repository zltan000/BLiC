import h5py
import numpy as np
from multiprocessing.pool import Pool
import time
# import psutil
import os
import math
import pandas as pd
from scipy.integrate import quad

print('Process start:', time.ctime())


# snapnum of the subhalo is the endpoint of the interpolation
# delta_t of the subhalo is begin at high redshift and end at lower redshift (delta_t=1)


def z2r(z):
    func2 = lambda y: 1 / (O_M * ((1 + y) ** 3) + O_K * ((1 + y) ** 2) + O_L) ** 0.5
    v2 = quad(func2, 0, z)
    r_c = c * v2[0] / 1e5  # comoving distance in Mpc/h
    return r_c


def assign_mag(isub):
    jumped = 5  # first N values are jumped to be not interpolated with mags(cause they are not mags)
    trackid = int(sub_tot['trackID'][isub])
    index_in_file = trackid % filewidth
    dt = sub_tot['delta_t'][isub]
    isnap = int(sub_tot['snapNum'][isub])
    z_true = sub_tot['redshift_true'][isub]

    if index_in_file >= lengal[isnap]:
        print('warning: mag1 index out of range!!', trackid, isnap, index_in_file)
    mag1 = np.array(list(galtot[isnap][index_in_file]))  # lower redshift
    mag1 = np.array([99. if math.isnan(x) or x == -1 else x for x in mag1[int(jumped):]])
    # print(mag1, snap_dismod[isnap])
    mag1 -= snap_dismod[isnap]
    if isnap == 70 or isnap == 77:
        if index_in_file >= lengal[isnap - 2]:
            print('warning: mag0 index out of range!!', trackid, isnap - 2, index_in_file)
            return mag1
        mag0 = np.array(list(galtot[isnap - 2][index_in_file]))  # higher redshift
        interpolated = mag0
        mag0 = np.array([99. if math.isnan(x) or x == -1 else x for x in mag0[int(jumped):]])
        mag0 -= snap_dismod[isnap - 2]
        # should not be [3:], cause 5 keys are from gal table!!!, now 'CentralMvir', 'CentralGal' are wrong valued!!
    else:
        if index_in_file >= lengal[isnap - 1]:
            print('warning: mag0 index out of range!!', trackid, isnap - 1, index_in_file)
            return mag1
        mag0 = np.array(list(galtot[isnap - 1][index_in_file]))  # higher redshift
        interpolated = mag0
        mag0 = np.array([99. if math.isnan(x) or x == -1 else x for x in mag0[int(jumped):]])
        mag0 -= snap_dismod[isnap - 1]

    interpolated[int(jumped):] = mag0 + dt * (mag1 - mag0)
    interpolated[int(jumped):] += 25 + 5 * np.log10(z2r(z_true) * (1 + z_true) / h)
    if interpolated[0] != trackid:
        print('warning: trackID not equal to galID!!', trackid, isnap, index_in_file)
    return interpolated


process = 72
pid = os.getpid()
# psutil_process = psutil.Process(pid)

zlist = np.array(pd.read_csv('/home/zltan/9t/galaxy/redshift.csv')['z'])
# setting paras
O_M = 0.3111
O_L = 1 - O_M
O_K = 0
t_H = 9785641806 * 1e-9  # Hubble time in Gyr/h 要根据模拟输出调整
c = 299792458
pc = 3.08568025e16
h = 0.6766
snap_dismod = []
for i in zlist[:-1]:
    snap_dismod.append(25 + 5 * np.log10(z2r(i) * (1 + i) / h))
snap_dismod.append(0)
snap_dismod = np.array(snap_dismod)
print('snap dis modules:', snap_dismod)

# snapnum loop
ibranch = 2688  ### 

needed_sub_fields = ['trackID', 'posx', 'posy', 'posz', 'velx', 'vely', 'velz', 'v_lfs', 'd_comoving', 'RA_deg',
                     'Dec_deg', 'redshift_true', 'redshift_obs', 'snapNum', 'delta_t']

sub_dtype = [('trackID', '<i4'), ('posx', '<f4'), ('posy', '<f4'), ('posz', '<f4'), ('velx', '<f4'), ('vely', '<f4'),
             ('velz', '<f4'), ('v_lfs', '<f4'), ('d_comoving', '<f4'), ('RA_deg', '<f4'), ('Dec_deg', '<f4'),
             ('redshift_true', '<f4'), ('redshift_obs', '<f4'), ('snapNum', '<i4'), ('delta_t', '<f4')]
gal_dtype = [('GalID', '<i4'), ('Type', '<i4'), ('StellarMass', '<f4'), ('CentralMvir', '<f4'), ('CentralGal', '<i4')]
csst_names = ['csst_nuv', 'csst_u', 'csst_g', 'csst_r', 'csst_i', 'csst_z', 'csst_y']  # 0~6
sdss_names = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0', 'wise_w1', 'wise_w2']  # 7~13
desi_names = ['DECaLS_g', 'DECaLS_r', 'DECaLS_z', 'BASS_g', 'BASS_r', 'MzLS_z_corrected']  # 14~19
output_list = []
for j in ['_nodust', '_dust']:
    for i in csst_names + sdss_names + desi_names:
        output_list.append(i + j)
data_types = ['<f4' for _ in range(len(output_list))]
mag_dtype = [(output_list[i], data_types[i]) for i in range(len(data_types))]

finished = np.load('gal_finished.npy')
for iibranch in range(ibranch, (ibranch + 1)):
    if iibranch in finished:
        print("%d: Succusess generated !" % iibranch, time.ctime())
        continue
    # load subhalos
    print('ibranch start:', iibranch)
    subdir = '/home/cossim/Jiutian/M1000/lightcones/fullsky_z2.5/9tmock_z2.5_%d.npy' % iibranch
    sub_tot = np.load(subdir)[needed_sub_fields]
    snap_st = min(sub_tot['snapNum'])
    snap_ed = max(sub_tot['snapNum'])
    ssn = 128
    filewidth = int(5e5)

    # load all galaxies in this branch
    galtot = [[] for i in range(ssn)]
    lengal = np.zeros(ssn)
    if snap_st == 70 or snap_st == 77:
        snap_st -= 2
    else:
        snap_st -= 1
    gal_fields = ['GalID', 'Type', 'StellarMass', 'CentralMvir', 'CentralGal']
    gal_fields = gal_fields + output_list
    galdir = '/home/cossim/Jiutian/M1000/lightcones/gal-tab/%03d/SubSnap%03d_%d.hdf5'
    for i in range(snap_st, snap_ed + 1):
        if i == 69 or i == 76:
            continue
        else:
            try:
                galtot[i] = h5py.File(galdir % (i, i, iibranch))['Galaxies'][:][gal_fields]
                lengal[i] = len(galtot[i])
            except:
                continue
    print('len gal:', lengal)

    # print("Memory usage:", round(psutil_process.memory_info().rss / 1024**3, 2), "Gbytes")
    # get interpolated mags
    pool = Pool(processes=process)
    res = pool.map(assign_mag, np.arange(len(sub_tot)))
    # res = pd.DataFrame(res)
    res = np.array(res)
    pool.close()
    pool.join()
    print('Pool done:', time.ctime())
    # print("Memory usage:", round(psutil_process.memory_info().rss / 1024**3, 2), "Gbytes")
    print(len(sub_tot), len(res))

    # compile data into hdf5
    sub_tot = np.array(pd.DataFrame(sub_tot))
    print(sub_tot.shape, res.shape)

    # pd_dtype = {i: int if ii == '<i4' else float for i, ii in enumerate(sub_dtype + gal_dtype + mag_dtype)}
    result = pd.DataFrame(np.concatenate((sub_tot, res), axis=1))
    # print(result.shape)
    emptyarr = np.empty(len(result), dtype=sub_dtype + gal_dtype + mag_dtype)
    i1 = 0
    for i2, i3 in (sub_dtype + gal_dtype + mag_dtype):
        emptyarr[i2] = np.array(result.iloc[:, i1])
        i1 += 1
    # print("Memory usage:", round(psutil_process.memory_info().rss / 1024**3, 2), "Gbytes")
    output_file_name = '/home/cossim/Jiutian/M1000/lightcones/fullsky_z2.5/galaxies/9tmockGal_z2.5_%d.hdf5' % iibranch
    # np.save(output_file_name, result)
    with h5py.File(output_file_name, "w") as output_file:
        output_file.create_dataset('Galaxies', data=emptyarr)
        # output_file['Galaxies'][:] = result
        print(time.ctime())  # timenow
    # del res, result, emptyarr, sub_tot, galtot
    # print("Memory usage:", round(psutil_process.memory_info().rss / 1024 ** 3, 2), "Gbytes")
    print("%d: Succusess generated !" % iibranch, time.ctime(), 'path: ', output_file_name)

