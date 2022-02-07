#!/usr/bin/env python

# intention: converts the .dat files from md_tian to readable .POSCAR files!

import os
from os.path import join

path = os.getcwd() #get directory
file_names = os.listdir(path) #get filenames in directory

expr = ':1:' # search for expression to replace

for i, name in enumerate(file_names): # loop over filenames

    if name.startswith('poscar') and name.endswith('.dat'):

        infile = open(name, 'r') # read out of file
        out_string = name[:-4] + '.POSCAR' # cut off .dat and add .POSCAR
        outfile = open(out_string, 'w') # write to file

        for ind, line in enumerate(infile):

            if expr in line:

                change = line[4:] # if expression to replace is found replace
                outfile.write('{}\n'.format(change.strip())) # write to file

            else:
                outfile.write('{}'.format(line)) # anything else will be the same

        infile.close() # close readfile
        os.remove(name)
        outfile.close() # close writefile
