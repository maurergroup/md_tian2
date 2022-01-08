#!/usr/bin/env/ python

# intention: get grid structures

# use like: python <scriptname> <input.POSCAR> <files_prefix>

import sys
import numpy as np

# give ranges (start stop step) for C and H:
start_C = -2.0
stop_C  = 2.5
step_C  = 0.1

start_H = 0.8
stop_H  = 6.0
step_H  = 0.1


##################### NO CHANGES BEYOND THIS LINE ##################

inpfile_name=sys.argv[1]
name_input=sys.argv[2]


slab  = 0
thres = 0.000000001
name  = str(name_input)


class atom():
    '''
    this part is copied as is from RuNNerUC.py
    '''
    def __init__(self, **kwargs):
        self.default_f=False
        try: 
            self.element=kwargs['element']
            #print(self.element)
        except:
            print('Atom was initialized without element type')
            sys.exit()
        try:
            self.xyz=[kwargs['x'],kwargs['y'],kwargs['z']]
            #print(self.xyz)
        except:
            print('Atom was initialized without xyz')
            sys.exit()
        fx=kwargs.pop('fx',99.999999)
        fy=kwargs.pop('fy',99.999999)
        fz=kwargs.pop('fz',99.999999)
        if abs(fx-99.999999)<thres or abs(fy-99.999999)<thres or abs(fz-99.999999)<thres:
            self.default_f=True
        self.force=[fx,fy,fz]
        atom_charge=kwargs.pop('atom_charge',0.)
        atom_energy=kwargs.pop('atom_energy',0.)
        self.atom_charge=atom_charge
        self.atom_energy=atom_energy       


def read_poscar(inpfile_name):
    # open file
    inpfile=open(inpfile_name, 'r')
    #comment line
    comments=str(inpfile.readline())
    #scaling factor
    sf=float(inpfile.readline().strip())
    #lattices
    lattice=[None,]
    l=[]
    for i in range(0, 3):
        line=inpfile.readline().strip()
        spline=line.split()
        l.append([sf*float(spline[j]) for j in range(0, 3)]) #notice that we need to multiply by the scaling factor
    lattice=np.array(l)
    # elements
    elementline = str(inpfile.readline().strip())
    # number of each element
    countline = str(inpfile.readline().strip())
    count = countline.split()
    count=[int(c) for c in count if c[0]!="!"]
    # Cartesian, Direct etc.
    typeline = str(inpfile.readline().strip())
    # position of atoms
    atoms = []
    natoms = sum(count)
    for i in range(0, natoms):
        spline=inpfile.readline().strip().split()
        pos=[sp for sp in spline[0:3]]
        pos=[sf*float(p) for p in pos]
        thisatom=atom(element="XX", x=pos[0], y=pos[1], z=pos[2])
        atoms.append(thisatom)
    # close file
    inpfile.close()
    #return values
    return(comments, sf, lattice, elementline, countline, typeline, atoms)


def write_poscar(name, slab, comments, sf, lattice, elementline, countline, typeline, atoms):
    outfile_name = name+"_"+str(slab)
    outfile = open(outfile_name, 'w')
    # comment
    outfile.write(comments)
    # scaling factor
    outfile.write(str(sf)+"\n")
    # lattice
    for l in lattice:
        outfile.write("{:12.6f} {:12.6f} {:12.6f}\n".format(l[0]*1, l[1]*1, l[2]*1))
    # elements
    outfile.write(elementline+"\n")
    # number of each element
    outfile.write(countline+"\n")
    # Cartesian, Direct etc.
    outfile.write(typeline+"\n")
    # positions of atoms
    for at in atoms:
        pos = at.xyz
        outfile.write("{:12.6f} {:12.6f} {:12.6f}\n".format(pos[0], pos[1], pos[2]))

    outfile.close()


#name = "slab"
#slab=1
#write_poscar(name, slab, comments, sf, lattice, elementline, countline, typeline, atoms)
#exit()    

# read poscar file
comments, sf, lattice, elementline, countline, typeline, atoms = read_poscar(inpfile_name)

#for h_z in np.linspace(0.8, 6.0, num=53, endpoint=True):
final_stop_C = stop_C + step_C # start stop step, end has to be stop+step!!
final_stop_H = stop_H + step_H # start stop step, end has to be stop+step!!

if step_C == 0 and step_H != 0:
    for h_z in np.arange(start_H, final_stop_H, step_H):
        slab += 1
        for ind,at in enumerate(atoms):
            pos = at.xyz
            if ind == 0:
                pos[2] = float(h_z)
        write_poscar(name, slab, comments, sf, lattice, elementline, countline, typeline, atoms)

elif step_C != 0 and step_H == 0:
    for c_z in np.arange(start_C, final_stop_C, step_C):
        slab += 1
        for ind,at in enumerate(atoms):
            pos = at.xyz
            if ind == 2: # I always place the H above the second atom!
                pos[2] = float(c_z)
        write_poscar(name, slab, comments, sf, lattice, elementline, countline, typeline, atoms)

elif step_C != 0 and step_H != 0:
    for c_z in np.arange(start_C, final_stop_C, step_C):
        for h_z in np.arange(start_H, final_stop_H, step_H):
            slab += 1
            for ind,at in enumerate(atoms):
                pos = at.xyz
                if ind == 0:
                    pos[2] = float(h_z)
                if ind == 2:
                    pos[2] = float(c_z)
            write_poscar(name, slab, comments, sf, lattice, elementline, countline, typeline, atoms)
else:
    print("Error: No proper step given for neither H nor C!")
    sys.exit()
