# from .templates.CalcBoxes import CalcBoxes
import os
import subprocess
from pathlib import Path
import pickle
import numpy as np


def init(path_ini=None):
    # initialize of BLiC
    # a new file (BLic.init) will be generated to get all parameters needed
    path = os.path.abspath(os.path.dirname(__file__))
    if path_ini:
        mkpath = path_ini + '/BLiC_out'
        if not Path(mkpath).is_dir():
            Path(mkpath).mkdir()
        ini_file = path_ini + '/BLiC_out/blic.ini'
    else:
        mkpath = path + '/../BLiC_out'
        if not Path(mkpath).is_dir():
            Path(mkpath).mkdir()
        ini_file = path + '/../BLiC_out/blic.ini'

    cmd = 'touch ' + ini_file
    p = subprocess.Popen(cmd, shell=True)
    p.wait()

    # read init.txt and generate new file
    f = open(path + '/templates/init_info.txt', 'r+')
    text = f.readlines()
    f.close()
    f = open(ini_file, 'w+')
    f.writelines(text)
    f.close()
    print('Initial file generated at: "%s"' % ini_file)


def findreplace(text, stra, strb):
    linepos = 0
    j = 0
    while j < len(text):
        if text[j][0:len(stra)] == stra:
            linepos = j
            break
        j += 1
    if linepos == 0 and text[linepos][0:len(stra)] != stra:
        print('Warning: replace defeat !', stra)
    else:
        text[linepos] = strb
    return text


def count_boxes(path_ini=None):
    path = os.path.abspath(os.path.dirname(__file__))
    if not path_ini:
        path_ini = path + '/../BLiC_out/blic.ini'
    f = open(path + '/templates/boxes.py', 'r+')
    flist = f.readlines()
    f.close()

    text_st = 0
    text_ed = 0
    j = 0
    while j < len(flist):
        if flist[j][0:-1] == "# beginning of .ini":
            text_st = j
        if flist[j][0:-1] == "# end of .ini":
            text_ed = j
            break
        j += 1
    flist_ed = flist[text_ed:]

    f = open(path_ini, 'r+')
    attached = f.readlines()
    f.close()

    flist[text_st + 1:] = attached
    flist[-1:] = flist_ed
    user_path = os.path.abspath(os.path.dirname(__file__)) + '/user_generated/'
    f = open(user_path + '/a_boxes.py', 'w+')
    f.writelines(flist)
    f.close()

    cmd = 'python3 ' + user_path + '/a_boxes.py'
    p = subprocess.Popen(cmd, shell=True)
    p.wait()  # run box calulator

    # generate .sub file
    cb = set_cb_class('r')
    if not cb.output_path:
        cb.output_path = path_ini.replace('blic.ini', '')
    set_cb_class('w', cb)

    f = open(path + '/templates/lightcones.sub', 'r+')
    text = f.readlines()
    f.close()

    path_ini = path_ini.replace('blic.ini', 'lightcones.sub')
    text = findreplace(text, "#PBS -l nodes=1:ppn=1\n", '#PBS -l nodes=1:ppn=%d\n' % cb.cores)
    opath = path_ini.replace('lightcones.sub', 'output_${PBS_JOBNAME}')
    text = findreplace(text, "#PBS -o \n", '#PBS -o %s\n' % opath)
    lightpath = os.path.abspath(os.path.dirname(__file__)) + '/templates/lightcones.py'
    text = findreplace(text, "python -u lightcones.py\n", 'python -u %s\n' % lightpath)
    f = open(path_ini, 'w+')
    f.writelines(text)
    f.close()
    print('Initial submit file generated at: "%s"' % path_ini)


def run(path_ini=None):  # ='./blic.ini'
    path = os.path.abspath(os.path.dirname(__file__))
    if not path_ini:
        path_ini = path + '/../BLiC_out/blic.ini'
    cb = set_cb_class('r')
    if cb.max_branch:
        cb.max_branch = 0
    else:
        subfile = path_ini.replace('blic.ini', '') + 'lightcones.sub'
        cmd = 'qsub %s' % subfile
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print('Submitted!')


# def run_example(path):
#
#
# # generate .ini
#
# # run .ini
#
#
def set_interpo(props, cases):
    
    if len(props) != len(cases):
        print("input information do not consistent! Please check inputs")
        return 0
    cb = set_cb_class('r')

    for i, iprop in enumerate(props):
        prop_index = cb.props.index(iprop)
        cb.interpo[prop_index] = cases[i]

    set_cb_class('w', cb)


def set_hdf_database(stra):
    cb = set_cb_class('r')
    cb.hdf_database = str(stra)
    set_cb_class('w', cb)


def set_unit(num):
    cb = set_cb_class('r')

    num = np.array(num)
    if num.size > 1:
        cb.pos_unit = num[0]
        cb.vel_unit = num[1]
    else:
        cb.pos_unit = num

    set_cb_class('w', cb)


def set_cb_class(mode, classcb=None):
    from .CalcBoxes import CalcBoxes
    user_path = os.path.abspath(os.path.dirname(__file__)) + '/user_generated/'
    if mode == 'r':
        with open(user_path + '/cb.pkl', 'rb') as f:
            cb = pickle.loads(f.read())
        f.close()
        return cb
    
    elif mode == 'w':
        save_class = open(user_path + '/cb.pkl', 'wb')
        strs = pickle.dumps(classcb)
        save_class.write(strs)
        save_class.close()
        

# some tools may useful

# tool 1: convert R.A. & dec to a percent of the full sky
def fullsky(RAmin, RAmax, decmin, decmax):
    """
    :param RAmin: the min value of R.A., between 0 and RAmax,
    :param RAmax: the max value of R.A., between RAmin and 360,
    :param RAmin: the min value of declination, between -90 and decmax,
    :param RAmin: the max value of declination, between decmin and 90,
    :return: the percentile of input field of view, compared with the full sky
    """
    if 0 <= RAmin < RAmax < 360 and -90 <= decmin < decmax <= 90:
        para = (-np.cos(np.deg2rad(90-decmin)) + np.cos(np.deg2rad(90-decmax))) * (np.deg2rad(RAmax) - np.deg2rad(RAmin)) / (4 * np.pi) * 100
        return para
    else:
        print('Please check inputs of the R.A. and the declination, in units of degrees.')
        return 0

