import numpy as np
import os
import random
from scipy.io import savemat


def run():
    filename = "dm_wc.txt"
    # this executable file is compiled by one of the above folder.
    # -t parameter is used to estimate the influence spread of a seed set.
    command = "max_inf6.exe -t Test.txt 100 100 GC_spread.txt < " + filename + " > out.txt"
    os.system(command)


def read_influence_spread():
    x = 1
    filename = "GC_spread.txt"
    f = open(filename, 'r')
    lines = f.readlines()
    i = lines[-1]
    i = i.split()
    x = eval(i[-1])
    return x


def set_seeds(seeds):
    filename = "Test.txt"
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

# seedSet.npy has a list in python, contains the seed set.
b = np.load("seedSet.npy", allow_pickle=True)
dp = 400
dn = 100
a = 10
data = 'dm_real'
print(np.shape(b))
print(b[0][0][0])

ratio = calculate_ratio(dp, dn, a, b[i][j][k])
# here we save the result as ".mat" because we use matlab to draw the Figure.
savemat('MonteCarlo_' + data + '.mat', {'MonteCarlo_' + data: ratio})


