#!/usr/bin/env python

# intention: creates given number of input files for md_tian

import sys
import numpy as np

infile_string=sys.argv[1] # for example md_tian.inp
file_number=sys.argv[2] # number of input files

expr = 'start'

#change = []
#end = []

n = 0

for i in range(1, int(file_number)+1):
    infile = open(infile_string, 'r')
    outfile_string = infile_string + '.' + str(i) # creates for example md_tian.inp.1
    outfile = open(outfile_string, 'w')

    result = ((n * 10000) + 1)

    for line in infile:
        if expr in line:
            change = line.split()
            change[-1] = result
            outfile.write('{} {}\n'.format(change[0],change[1]))
        else:
            outfile.write('{}'.format(line))

    outfile.close()
    infile.close()
    n += 1
