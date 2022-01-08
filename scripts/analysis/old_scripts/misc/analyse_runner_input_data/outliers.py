import math, sys
import numpy as np
def get_distances(filename):
    #returns first neighbor distances, ie, minimum distances
    #can be with any atom, or only those of the same element

    infile=open(filename, mode="r")
    dist_array=[]
    natoms=0

    for line in infile:
        spline=line.split()
        if "begin" in spline[0]: #reset stuff
            pos_array=[]
            elem_array=[]
            latt_counter=0
            latt=np.array([0.0, 0.0, 0.0])
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
            min_dist=10000.0
            if natoms>1:
                for (a1, e1) in zip(pos_array, elem_array):
                    for (a2, e2) in zip(pos_array, elem_array):
                        #mic
                        d=(a2-a1)
                        d=d-np.rint(d/latt)*latt
                        dist=np.linalg.norm(d)
                        if (dist<min_dist) and (dist>0.0001):
                        #prevents comparing distance with oneself
                            min_dist=dist
            dist_array.append(min_dist)

    infile.close()
    #print(len(dist_array))
    #print(dist_array)
    return dist_array

#format: atom x y z element ?? ?? fx fy fz
#unit conversion
ha2ev=27.2114
b2a=0.52917

nstructs=50324
outfile1_dir="inliers.data"
outfile2_dir="outliers.data"
infile_dir=sys.argv[1]
outfile1=open(outfile1_dir, mode="w")
outfile2=open(outfile2_dir, mode="w")
#infile=open(infile_dir, mode="r")

#process command line input
nopts=len(sys.argv)
print(nopts)
mode=[]
thr=[]
n=2
while n<nopts:
    mode.append(sys.argv[n])
    thr.append(float(sys.argv[n+1]))
    n+=2
print(mode, thr)

#find outliers
out_index=[]
out=False #avoids multiple counting of same structure
n=0

if "force" in mode:
    infile=open(infile_dir, mode="r")
    for line in infile:
        spline=line.split()
        if "begin" in spline[0]:
            out=False
            n=n+1

        if "atom" in spline[0]:
            #spline=line.split()
            #print(spline)
            elm=spline[5-1]
            #x=str(b2a*float(spline[2-1]))
            #y=str(b2a*float(spline[3-1]))
            #z=str(b2a*float(spline[4-1]))
            #outfile.write(elm+" "+x+" "+y+" "+z+"\n")
            fx=float(spline[8-1])*(ha2ev/b2a)       
            fy=float(spline[9-1])*(ha2ev/b2a)
            fz=float(spline[10-1])*(ha2ev/b2a)
            ft=math.sqrt(fx**2+fy**2+fz**2)

            if ft>thr[mode.index("force")]:
                print("Structure: "+str(n)+" Atom: "+elm+" Force too big: "+str(fx)+" "+str(fy)+" "+str(fz)+" "+str(ft)) 
                if not out:
                    out=True
                    if not (n in out_index):
                        out_index.append(n)
    infile.close()

if "distance" in mode:
    dist_array=get_distances(infile_dir)
    for n,d in enumerate(dist_array):
        if d<thr[mode.index("distance")]:
            print("Structure: "+str(n+1)+" Distance too small: "+str(d))
            if not ((n+1) in out_index):
                out_index.append(n+1)

if "energy" in mode:
    n=0
    natoms=0
    infile=open(infile_dir, mode="r")
    for line in infile:
        spline=line.split()
        if "begin" in spline[0]:
            n=n+1
            natoms=0
        if "atom" in spline[0]:
            natoms+=1
        if "energy" in spline[0]:
            e=float(spline[1])*ha2ev
        if "end" in spline[0]:
            e=e/float(natoms)
            if e>thr[mode.index("energy")]:
                print("Structure: "+str(n)+" Energy too big: "+str(e)) 
                if not (n in out_index):
                    out_index.append(n)
    infile.close()

                    
print("Outliers "+str(len(out_index))+ " out of "+str(nstructs)+" Percentage: "+str(100.0*len(out_index)/float(nstructs)))

#output
#infile.seek(0) 
infile=open(infile_dir, mode="r")
n=0
for line in infile:
    if "begin" in line:
        n=n+1

    if (n in out_index):
        outfile2.write(line)
    else:
        outfile1.write(line)

infile.close()
outfile1.close()
outfile2.close()
