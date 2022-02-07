#!/usr/bin/env python

# intention: plots temperature and energies (total, potential, kinetic) of system, projectile as well as slab for each trajectory

import os
import matplotlib.pyplot as plt
import numpy as np

path = os.getcwd() # get directory
file_names = os.listdir(path) # get file names
sorted_file_names = file_names.sort() # file names are not sorted by default, which is very weird and confusing
plt_form = ".pdf"

for i, name in enumerate(file_names): # loop over filenames

    if name.endswith('.dat'): # we only want files which end on .dat (output by md_tian2)

        infile = open(name, 'r') # read out of file
        
        time = [] # create strings
        temp = []
        ekin_p = []
        ekin_l = []
#        epot = []
        e_total = []
        ekin = []
#        eges = []
        etot = []

        for ind,line in enumerate(infile): # loop over lines in file

            time.append(np.asarray(line.split())[0]) # extract data
            temp.append(np.asarray(line.split())[1])
            ekin_p.append(np.asarray(line.split())[5])
            ekin_l.append(np.asarray(line.split())[6])
#            epot.append(np.asarray(line.split())[11])
            e_total.append(np.asarray(line.split())[13])
            etot.append(abs(float(e_total[ind]) - float(e_total[0])) * 1000)
            ekin.append(float(ekin_p[ind]) + float(ekin_l[ind]))
#            eges.append(float(ekin_p[ind]) + float(ekin_l[ind]) + float(epot[ind]))

        change_temp = name[0:-4] + "_T+E" + plt_form # create plot filename out of data filename
        fig, ax1 = plt.subplots() # multiple plots
        ax1.plot(time, temp, 'r-')
        ax1.set_xlabel(r'$t$ / fs') # use LaTeX syntax
        ax1.set_ylabel(r'$T$ / K', color='r')
        ax1.tick_params('y', colors='r') # set axis numbers red
        ax2 = ax1.twinx() # multiple axis
#        ax2.plot(time, ekin_p,  'b--', label="Proj") # b-- for dashed lines
#        ax2.plot(time, ekin_l, label="Latt", color='b')
        ax2.plot(time, ekin,  'b-')
        ax2.set_ylabel(r'$E_{\mathrm{kin}}$ / eV', color='b')
        ax2.tick_params('y', colors='b')
#        ax2.ylim(1,3.3)
        #plt.legend() # without no legend will be printed
        fig.tight_layout() # not sure what this will do
        plt.xlim(0,float(time[-1])) # set range of x axis
        plt.savefig(change_temp) # save plot
        plt.close() # close plot


        change_e = name[0:-4] + "_E" + plt_form # repeat for next data
#        plt.plot(time, eges, 'r--', label=r"$E_{\mathrm{ges}}$")
        plt.plot(time, etot, 'r-')
        plt.xlabel(r'$t$ / fs')
        plt.ylabel(r'$\Delta E_{\mathrm{tot}}$ / meV')
        plt.xlim(0,float(time[-1]))
        
#        plt.legend()
        plt.savefig(change_e)
        plt.close()

        infile.close() # close file if not needed anymore
