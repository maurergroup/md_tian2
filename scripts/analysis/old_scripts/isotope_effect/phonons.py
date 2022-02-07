#!/usr/bin/env python

# intention: analyze trajectories due to energy transfer from slab (phonons) to projectile (energy gain of projectile)

import os
from os.path import join
import numpy as np

path = os.getcwd() #get directory
cur_dir = os.getcwd()
path2 = os.path.join(cur_dir, 'traj')
file_names = os.listdir(path) #get filenames in directory
sorted_file_names = file_names.sort() # get the right ordering, which results unfortunately not in the right order

ekin_init  = []
ekin_final = []
ekin_diff  = []

outfile = open('phonons.dat', 'w')

for name in file_names:

    if name.startswith('mxt_'):

        infile = open(name, 'r') # read out of file

        for line in infile:

            if 'ekin_p_i' in line:

                ekin_init.append(np.array([float(x) for x in line.split()[2:]]))

            if 'ekin_p_f' in line:

                ekin_final.append(np.array([float(x) for x in line.split()[2:]]))

        infile.close()


ekin_diff = np.array(ekin_init) - np.array(ekin_final)

n = 0
for ind, diff in enumerate(ekin_diff):

    if diff < 0:

        outfile.write('{} {}\n'.format(float(diff), file_names[ind]))
        n += 1

outfile.close()
print('Phonon energy transfer happens in', n, 'trajectories')
