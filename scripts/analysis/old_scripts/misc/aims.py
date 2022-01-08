#!/usr/bin/env python

# intention: in case of H atom adds the spin to the aims data file

import os
from os.path import join

path = os.getcwd() #get directory
file_names = os.listdir(path) #get filenames in directory

expr = 'H' # search for expression

for i, name in enumerate(file_names): # loop over filenames
    with open(name, "r+") as f:
        lines = f.read() # read everything in the file
        for line in lines:
            f.write(line)
            if expr in line:
                f.write("        initial_moment 1")
    f.close()




'''

for i, name in enumerate(file_names): # loop over filenames

    if name.endswith('.in'):

        infile = open(name, 'r') # read out of file
        out_string = name[:-3] + '.' # cut off .dat and add .POSCAR
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

'''
