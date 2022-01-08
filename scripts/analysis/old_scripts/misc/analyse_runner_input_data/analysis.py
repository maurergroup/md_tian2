import numpy as np
import matplotlib.pyplot as plt
import statistics as stats
from math import sqrt
import sys

def readcolumn(filename, column, mode='float'):
    data=np.loadtxt(filename, usecols=(column-1,), dtype=mode)
    #data=data[:,column-1]
    return data

def get_forces(filename, element):
    #format: atom x y z element ?? ?? fx fy fz
    element_array=readcolumn(filename, 5, mode='S2')
    #print(element_array)

    force_array=np.array([readcolumn(filename, 8), readcolumn(filename, 9), readcolumn(filename, 10)])
    force_array=force_array*ha2ev/b2a
    #print(force_array)

    #probably a way to do this with array masking, more efficiently
    del_indx=[]
    i=0
    for elm in element_array:
        if (elm.astype(str) != element):
            del_indx.append(i)
        i=i+1
    #print(del_indx)
    force_array=np.array([np.delete(force_array[0], del_indx), np.delete(force_array[1], del_indx), np.delete(force_array[2], del_indx)])

    force_tot=np.array([force_array[0]**2 + force_array[1]**2 + force_array[2]**2])
    force_tot=np.sqrt(force_tot)
    force_tot=np.concatenate(force_tot)

    #print("COMPONENTS")
    #print(force_array)
    #print("FORCE")
    #print(force_tot)

    return force_array[0], force_array[1], force_array[2], force_tot

def get_energy(filename, elem_array, atomic_energy_array):
    infile=open(filename, mode='r')
    natoms=0
    energy=0.0
    e_array=[]
    n=0
    dE=0.0
    for line in infile:
        spline=line.split()
        if "begin" in spline[0]:
            natoms=0
            energy=0.0
            dE=0.0
            #n+=1
            #print(n)
        if "atom" in spline[0]:
            natoms+=1
            dE+=-1.0*atomic_energy_array[elem_array.index(spline[4])]
        if "energy" in spline[0]:
            energy=float(spline[1])
        if "end" in spline[0]:
            energy=((energy*ha2ev)+dE)/natoms
            e_array.append(energy)
    infile.close()
    return e_array

def get_distances(filename, chosen_elem, same_element=True):
    #returns first neighbor distances, ie, minimum distances
    #can be with any atom, or only those of the same element

    infile=open(filename, mode="r")
    dist_array=[]
    n=0
    natoms=0

    for line in infile:
        spline=line.split()
        if "begin" in spline[0]: #reset stuff
            pos_array=[]
            elem_array=[]
            latt_counter=0
            latt=np.array([0.0, 0.0, 0.0])
            n+=1
            natoms=0

        if "lattice" in spline[0]:
            #spline=line.split()
            l=b2a*float(spline[latt_counter+2-1])
            latt[latt_counter]=l
            latt_counter += 1

        if "atom" in spline[0]:
            #spline=line.split()
            #print(spline)
            elm=spline[5-1]
            x=b2a*float(spline[2-1])
            y=b2a*float(spline[3-1])
            z=b2a*float(spline[4-1])
            pos_array.append(np.array([x, y, z]))
            elem_array.append(elm)
            natoms+=1
            
        if "end" in spline[0]: #process
            if natoms>1:
                for (a1, e1) in zip(pos_array, elem_array):
                    min_dist=10000.0
                    if e1==chosen_elem:
                        for (a2, e2) in zip(pos_array, elem_array):
                            #dist=np.linalg.norm(a2-a1)
                            #mic
                            d=(a2-a1)
                            d=d-np.rint(d/latt)*latt
                            dist=np.linalg.norm(d)
                            #print("A1: "+str(a1)+" "+e1+" A2: "+str(a2)+" "+e2+" Dist: "+ \
                            #str(dist)+" "+str(min_dist)+" Elem Lock "+str(elem_lock(e1, e2, same_element)))
                            if (dist<min_dist) and (dist>0.0001) and (elem_lock(e1, e2, same_element)) and (e1 == chosen_elem): 
                            #prevents comparing distance with oneself
                                min_dist=dist
                        dist_array.append(min_dist)
                        if min_dist>100.0:
                            print(n)

    infile.close()
    #print(len(dist_array))
    #print(dist_array)

    return dist_array

def elem_lock(e1, e2, same_element):
    res=True
    if same_element:
        if (e1 == e2):
            res=True
        else:
            res=False
    return res

def do_stats(data, count, mean_color='black', median_color='blue', stdev_color='red', plot=True):
    mean  =stats.mean(data)
    median=stats.median(data)
    stdev =stats.pstdev(data)
    print("Max: "+str(max(data))+" Min: "+str(min(data))+" Mean: "+str(mean)+" Median: "+str(median)+" Std. Dev.:"+str(stdev))

    if plot:
        plt.plot((mean, mean),     (0.0, max(count)), mean_color)
        plt.plot((median, median), (0.0, max(count)), median_color)
        plt.plot((mean+stdev, mean+stdev),     (0.0, max(count)), stdev_color)
        plt.plot((mean-stdev, mean-stdev),     (0.0, max(count)), stdev_color)



#files
atom_file_dir="atom_info.out"
infile_dir=sys.argv[1]

#unit conversion
ha2ev=27.2114
b2a=0.52917

do_xyz=False
do_energy=True
do_force=True
do_distance=True
elem_array=["C", "H"]                                                                            	## VanSibner: Change element here
atomic_energy_array=[0.0, 0.0] #isolated atom energy for each, in eV       				## VanSibner: Change atomic energy here (if you had them just put them in outherwise fill in zeros; energy convert is crucial to get rid of extremely large numbers and spare 
ha_energy=[0.0, 0.0]											## 	      caused float errors)

#export to xyz
if do_xyz:
    print("EXPORTING TO XYZ")
    outfile_dir="structure.lammpstrj"
    #infile_dir="ruNNer.input.small"
    outfile=open(outfile_dir, mode="w")
    infile=open(infile_dir, mode="r")

    nsteps=0
    natoms=0
    for line in infile:
        spline=line.split()
        if "begin" in spline[0]:
            nsteps+=1
            natoms=0
            atom_array=[]
            latt_counter=0
            latt_array=[]
            outfile.write("ITEM: TIMESTEP\n")
            outfile.write("{}\n".format(nsteps))
        if "lattice" in spline[0]:
            latt_array.append(spline[latt_counter+1])
            latt_counter+=1
        if "atom" in spline[0]:
            elm=spline[5-1]
            x=str(float(spline[2-1])*b2a)
            y=str(float(spline[3-1])*b2a)
            z=str(float(spline[4-1])*b2a)
            atom_array.append([elm, x, y, z])
            natoms+=1
        if "end" in spline[0]:
            outfile.write("ITEM: NUMBER OF ATOMS\n")
            outfile.write("{}\n".format(natoms))
            outfile.write("ITEM: BOX BOUNDS pp pp pp\n")
            outfile.write("0.0 {}\n".format(latt_array[0]))
            outfile.write("0.0 {}\n".format(latt_array[1]))
            outfile.write("0.0 {}\n".format(latt_array[2]))
            outfile.write("ITEM: ATOMS id element x y z\n")
            for i,a in enumerate(atom_array):
                outfile.write("{} {} {} {} {}\n".format(i, a[0], a[1], a[2], a[3]))

    infile.close()
    outfile.close()
            
if do_energy:
    #get data
    #energy=readcolumn(energy_file_dir, 2)
    #energy=energy*ha2ev
    energy=get_energy(infile_dir, elem_array, atomic_energy_array)

    #plot data
    ##ENERGY
    print("ENERGY")
    plt.figure()

    n, bins, patches = plt.hist(energy, bins=25, normed=False, facecolor='green', alpha=0.75)           ## VanSibner: bins are the number of boxes in the histogram, should be not too large as well as not too small
    plt.suptitle("Energy Histogram")
    plt.axis([min(energy), max(energy), 0, max(n)])
    plt.xlabel('Energy (eV/atom)')
    plt.ylabel('Count')
    do_stats(energy, n)
    plt.savefig("01_energy_histogram.png")

if do_force:
    ##FORCES
    print("FORCES")
    for elem in elem_array:
        f_x, f_y, f_z, f_tot=get_forces(atom_file_dir, element=elem)
        if len(f_x)>0:
            plt.figure()
            n, bins, patches = plt.hist([f_x, f_y, f_z], bins=31, normed=False, stacked=False, label=["fx", "fy", "fz"])#,facecolor='blue', alpha=0.75)
            plt.axis([-6, 6, 0, max((max(n[0]), max(n[1]), max(n[2])))])    ## VanSibner: change range if necessary
            plt.xlabel('{} Force Component (eV/A)'.format(elem))
            plt.ylabel('Count')
            plt.suptitle("{} Force Component Histogram".format(elem))
            plt.legend()
            print("{} Individual Force Components".format(elem))
            print("F_x")
            do_stats(f_x, n[0], plot=False)
            print("F_y")
            do_stats(f_y, n[1], plot=False)
            print("F_z")
            do_stats(f_z, n[2], plot=False)
            plt.savefig("02_{}_force_comp_histogram.png".format(elem))

            plt.figure()
            n, bins, patches = plt.hist(f_tot, bins=100, normed=False, stacked=False)#,facecolor='blue', alpha=0.75)
            plt.axis([min(f_tot), max(f_tot), 0, max(n)])
            plt.xlabel('{} Total Force (eV/A)'.format(elem))
            plt.ylabel('Count')
            print("{} Total Force".format(elem))
            do_stats(f_tot, n)
            plt.suptitle("{} Total Force Histogram".format(elem))
            plt.savefig("03_{}_total_force_histogram.png".format(elem))

if do_distance:
    ##DISTANCES
    print("DISTANCES")
    for elem in elem_array:
        dist=get_distances(infile_dir, chosen_elem=elem, same_element=False)
        if len(dist)>0:
            plt.figure()
            n, bins, patches = plt.hist(dist, bins=500, normed=False, stacked=False,facecolor='blue', alpha=0.75)
            plt.axis([min(dist), max(dist), 0, max(n)])
            plt.xlabel('{} Nearest Neighbor Distance (A)'.format(elem))
            plt.ylabel('Count')
            print("{} Nearest Neighbor Distances".format(elem))
            do_stats(dist, n)
            plt.suptitle("{} Nearest Neighbor Distance Histogram".format(elem))
            plt.savefig("05_{}_nn_dist_histogram.png".format(elem))

#plt.show()
