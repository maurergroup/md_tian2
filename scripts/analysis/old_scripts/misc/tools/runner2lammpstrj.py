import sys
import triclinic_functions as trif
#Unit Conversion
bohr2angst=0.529177
angst2bohr=1.0/bohr2angst
thr=0.000000000001 #threshold for detecting values close to zero for non diagonal elements

infile_name=sys.argv[1]
outfile_name=infile_name.split(".")[0]+".lammpstrj"
infile=open(infile_name, mode="r")
outfile=open(outfile_name, mode="w")

#process file
nstruct=0
for line in infile:
    spline=line.lstrip().split()
    if spline[0]=="begin":
        atom_array=[]
        bbox=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        latt_count=0
        natoms=0
        nstruct+=1
    elif spline[0]=="lattice":
        for i in range(0, 3):
            bbox[latt_count][i]=float(spline[i+1])*bohr2angst
        latt_count+=1

    elif spline[0]=="atom":
        natoms+=1
        atom_array.append([bohr2angst*float(spline[1]), bohr2angst*float(spline[2]), bohr2angst*float(spline[3]), spline[4]])
    elif spline[0]=="end":
        outfile.write("ITEM: TIMESTEP\n")
        outfile.write("{}\n".format(nstruct))
        outfile.write("ITEM: NUMBER OF ATOMS\n")
        outfile.write("{}\n".format(natoms))
        
        #Check for non orthogonal box:
        if trif.check_triclinic(bbox):
            outfile.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
            xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound, xy, xz, yz=trif.cellvectors_to_lammpstrjformat(bbox)
            outfile.write("{} {} {}\n".format(xlo_bound, xhi_bound, xy))
            outfile.write("{} {} {}\n".format(ylo_bound, yhi_bound, xz))
            outfile.write("{} {} {}\n".format(zlo_bound, zhi_bound, yz))
        else: #orthogonal, not much to do
            outfile.write("ITEM: BOX BOUNDS pp pp pp\n")
            for i in range(0, 3):
                outfile.write("0.0 {}\n".format(bbox[i][i]))

        outfile.write("ITEM: ATOMS id element x y z\n")
        for i in range(0, natoms):
            outfile.write("{0} {1[3]} {1[0]} {1[1]} {1[2]}\n".format(i+1, atom_array[i]))
infile.close()
outfile.close()
