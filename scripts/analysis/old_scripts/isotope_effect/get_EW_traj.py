#!/usr/bin/env python

# intention: used to get all the trajectories with extrapolation warnings in the output.dat file

infile_name = 'output.dat'
outfile_name = 'recalc_EW_traj.dat'

traj_counter = 1
warning_counter = 0
total_counter = 0

infile = open(infile_name, 'r')
outfile = open(outfile_name, 'w')

for line in infile:

    if 'Eref' in line:

        if (warning_counter > 0):

            outfile.write('Trajectory Number: {:5d} Number of Extrapolation Warnings: {:5d}\n'.format(traj_counter, warning_counter))
            total_counter += 1

        warning_counter = 0
        traj_counter += 1

    if 'WARNING' in line:

        warning_counter += 1


outfile.write('\nTotal number of trajectories with Extrapolation Warnings: {}'.format(total_counter))

infile.close()
outfile.close()


#2do:
#
# - loop over trajectory folders to be independent of number and don't need a bash script anymore!
