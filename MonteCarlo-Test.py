import numpy as np
import os
import random
from scipy.io import savemat


def run():
    filename = "dm_real.txt"
    # this executable file is compiled obove folder.
    command = "max_inf6.exe -t Test_dm_real.txt 100 100 GC_spread_dm_real.txt < C:\\Users\\508\source\\repos\\max_inf6\\data_example\\" + filename + " > out.txt"
    os.system(command)


def read_influence_spread():
    x = 1
    filename = "GC_spread_dm_real.txt"
    f = open(filename, 'r')
    lines = f.readlines()
    i = lines[-1]
    i = i.split()
    x = eval(i[-1])
    return x


def set_seeds(seeds):
    filename = "Test_dm_real.txt"
    f = open(filename, 'w+')
    f.write(str(len(seeds)) + '\n')
    for i in seeds:
        f.write(str(i) + '\n')


def calculate_ratio(dp, dn, a, seeds_list):
    newdp = dp
    newdn = dn
    rate = newdn / (newdn + newdp) 
    for seeds in seeds_list:
        if seeds == []:
            for i in range(a):
                if random.random() > rate:
                    newdp = newdp + 1
                else:
                    newdn = newdn + 1
        else:
            set_seeds(seeds)
            run()
            spread = read_influence_spread()
            rate = newdn / (newdn + newdp)
            for i in range(a):
                if random.random() > rate:
                    newdp = newdp + 1
                else:
                    newdn = newdn + 1
            newdn = newdn + spread

    return newdn / newdp


b = np.load("dm_real.npy", allow_pickle=True)
dp = 400
dn = 100
a = 10
data = 'dm_real'
print(np.shape(b))
print(b[0][0][0])
ratio_list_round =[]
for i in range(3):
    ratio_list_k = []
    for j in range(3):
        ratio_list_alg = []
        for k in range(6):
            ratio = calculate_ratio(dp, dn, a, b[i][j][k])
            ratio_list_alg.append(round(ratio, 4))
        ratio_list_k.append(ratio_list_alg)
    ratio_list_round.append(ratio_list_k)
savemat('MonteCarlo_' + data + '.mat', {'MonteCarlo_' + data: ratio_list_round})


