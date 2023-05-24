import subprocess
import time
import numpy as np


gal_bin = 2689
start_bin = 0
end_bin = 2689
filepath = '/home/zltan/9t/'
subfile = 'b_lightcone.sub'
start_ssn = np.load('0_start_ssn.npy')

pool_num = 0
pool_num_max = 18
check_str = 'output_docu'
subname = 'b_lightcone.py'  ###
lcname = '9t-fullsky'  ###
output_duc = './out_fs_z2/4.1_output_' + lcname + '-'  ###
i = start_bin

# while i < end_bin:  # for i in [2195, 2270, 2681, 2682, 2683, 2684, 2685, 2686, 2687, 2688]:
while i < end_bin:
    running = 0
    nofind = 0
    for j in range(0, pool_num_max*2):
        try:  # 有可能还在排队
            if i-1-j >= 0:
                f = open(filepath + output_duc + '%d' % (i-1-j), 'r+')
                flist = f.readlines()
                f.close()
                if flist[-1][0:len(check_str)] != check_str:
                    running += 1
        except:
            print('no find:', i-1-j)
            nofind = 1
            time.sleep(10)
    if running < pool_num_max and nofind == 0:
        time.sleep(3)
        f = open(filepath + subname, 'r+')
        flist = f.readlines()
        f.close()
        j = 0
        binpos = j
        while j < len(flist):
            if flist[j][4:-1] == "# parallel computing by galaxy bins":
                binpos = j
            j += 1
        if binpos == 0:
            print('Warning: galbin pos do not true!!')
        flist[binpos + 1] = '    iii = %d  ###\n' % i  # 行数 -1
        flist[binpos + 2] = '    ssn_start = %d  ###\n' % start_ssn[i]
        f = open(filepath + subname, 'w+')
        f.writelines(flist)
        f.close()

        f = open(filepath + subfile, 'r+')
        flist =f.readlines()
        f.close()
        flist[1] = '#PBS -N ' + lcname + '-%d\n' % i
        flist[24] = 'python -u  %s\n' % subname
        f = open(filepath + subfile, 'w+')
        f.writelines(flist)
        f.close()

        cmd = 'qsub %s' % subfile
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        time.sleep(5)
        cmd = 'qstat -u zltan'
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print('bin_now:', i, 'running:', running)
        i += 1
        time.sleep(5)

# 最后检查一遍有没有galbin失败报错了
for i in range(start_bin, end_bin):
    f = open(filepath + output_duc + '%d' % i, 'r+')
    flist = f.readlines()
    f.close()
    if flist[-1][0:len(check_str)] != check_str:
        print('Running or not exit! : %d' % i)