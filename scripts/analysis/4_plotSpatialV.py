#!/usr/bin/env python

# intetnion: plot spatial velocity distribution

# use like: python <scriptname>

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm

xmax = 8.5476
ymax = 7.4024
zmax = 13.0055 # this only for completeness


dat = np.genfromtxt('analysis/spatial_v.txt',skip_header=1)
X_dat = dat[:,0]
Y_dat = dat[:,1]
Z_dat = dat[:,2]*100



X, Y, Z, = np.array([]), np.array([]), np.array([])
for i in range(len(X_dat)):
    X = np.append(X, X_dat[i])
    Y = np.append(Y, Y_dat[i])
    Z = np.append(Z, Z_dat[i])

# create x-y points to be used in heatmap
xi = np.linspace(X.min(), X.max(), 1000)
yi = np.linspace(Y.min(), Y.max(), 1000)


# Interpolate for plotting
zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='cubic')

# I control the range of my colorbar by removing data 
# outside of my range of interest
#zmin = min(Z_dat)
zmin = 0
zmax = max(Z_dat)

zi[(zi<zmin) | (zi>zmax)] = None

'''
# Create the contour plot
CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow, vmax=zmax, vmin=zmin)



plt.plot([0, xmax], [0, 0], color='k', linestyle='-', linewidth=2)
plt.plot([0, xmax], [ymax, ymax], color='k', linestyle='-', linewidth=2)
plt.plot([0, 0], [0, ymax], color='k', linestyle='-', linewidth=2)
plt.plot([xmax, xmax], [0, ymax], color='k', linestyle='-', linewidth=2)

plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('spatial_v.pdf', format='pdf')
'''
# Create the contour plot
plt.figure()

origin = 'lower'

fig1, ax2 = plt.subplots(constrained_layout=True)
CS1 = ax2.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow, origin=origin, vmax=zmax, vmin=zmin)

#CS2 = ax2.contour(CS1, levels=CS.levels[::2], colors='r', origin=origin)

cbar = fig1.colorbar(CS1)
cbar.ax.set_ylabel('km / s')
# Add the contour line levels to the colorbar
#cbar.add_lines(CS2)


plt.plot([0, xmax], [0, 0], color='k', linestyle='-', linewidth=2)
plt.plot([0, xmax], [ymax, ymax], color='k', linestyle='-', linewidth=2)
plt.plot([0, 0], [0, ymax], color='k', linestyle='-', linewidth=2)
plt.plot([xmax, xmax], [0, ymax], color='k', linestyle='-', linewidth=2)

plt.xlabel('x')
plt.ylabel('y')
plt.savefig('analysis/spatial_v.pdf', format='pdf')

plt.figure()
fig1, ax2 = plt.subplots(constrained_layout=True)
plt.hist2d(X, Y, bins=200, cmap=plt.cm.rainbow,norm=LogNorm())
plt.colorbar(label='counts / bin')

plt.xlabel('x')
plt.ylabel('y')
plt.savefig('analysis/spatial_v_hist.pdf', format='pdf')
