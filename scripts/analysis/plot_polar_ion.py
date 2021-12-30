#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import griddata


parametername = "analysis/plot_parameter.txt" # conditions to write in plot
filename_ion  = "analysis/2d-ang-dist_ion_imaging_norm.txt"

parameterfile = open(parametername,"r")

for line in parameterfile:
    if not "Einc" in line: # skip first line
        einc = float(line.split()[0])
        vinc = float(line.split()[1])
        ainc = float(line.split()[2])
        temp = float(line.split()[3])
        bins = int(line.split()[5])

size_val = bins # should be equal to BINS in script
cmap = "nipy_spectral"
mintheta = 0
maxtheta = 90
angles = [0.0,30.0,60.0,90.0] # own ticks
angles.append(ainc)
labels = angles
levels = 20

x, y, z = np.loadtxt(filename_ion,unpack=True) # x=E_f / E_i, y=theta, z=flux

mintheta = -90
maxtheta = 90
angles = [-90.0,-60.0,-30.0,0.0,30.0,60.0,90.0] # own ticks
angles.append(ainc)
labels = angles
size_val_2d = 2*size_val # double the fun, double the gun (bin)

xi, yi = np.linspace(0, 1.1, size_val_2d), np.linspace(-90, 90, size_val_2d) 
xi, yi = np.meshgrid(xi, yi)
print("Plot polar_2d_ion.png")

zi = griddata((x, y), z, (xi, yi),method='linear',rescale=False)

fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')
ax.set_thetamin(mintheta)
ax.set_thetamax(maxtheta)
ax.set_thetagrids(angles, labels=labels)
ax.get_xticklabels()[-1].set_color("red")
ax.set_theta_direction(-1)
ax.set_theta_zero_location("N")
ax.set_xlabel(r'$E_{\mathrm{s}}\, / \, E_{\mathrm{i}}$')
ax.set_ylabel(r'$\vartheta\, / \, \circ$',rotation=0)
ax.yaxis.set_label_coords(0.3, 0.9)
ax.grid(linewidth=0.5,linestyle="--")
ax.set_yticks(np.linspace(0, 1, 3, endpoint=True))
ax.set_ylim(0,1.1)
ax.grid(linewidth=0.5,linestyle="--")
minor_locator_y = AutoMinorLocator(2)
ax.yaxis.set_minor_locator(minor_locator_y)
minor_locator_x = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minor_locator_x)
ax.text(0.05, 1.00, r'$T$ = {} K'.format(temp),verticalalignment='bottom', horizontalalignment='right',transform=ax.transAxes,color='k', fontsize=10,bbox={'facecolor': 'grey', 'alpha': 0.8, 'pad': 10})
ax.text(0.06, 0.75, r'$E_{{\mathrm{{i}}}}$ = {} eV'.format(einc),verticalalignment='bottom', horizontalalignment='right',transform=ax.transAxes,color='k', fontsize=10,bbox={'facecolor': 'grey', 'alpha': 0.8, 'pad': 10})
plt.grid(which='minor',linestyle="--")
ctf = ax.contourf(yi/360*2*np.pi, xi, zi, 20, cmap=cmap, levels=levels) # theta = yi/360*2*np.pi, r = xi, values = zi
plt.colorbar(ctf, pad=0.10)

plt.savefig("analysis/polar_2d_ion.png", dpi=600)
