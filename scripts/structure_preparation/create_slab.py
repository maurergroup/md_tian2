#!/usr/bin/env python3

# Intention: script to create fcc, bcc and diamond slabs with given miller indices for surface type

#use script like: 
#./scriptname <element symbol> <lattice type> <miller index x-axis> <miller index y-axis> <miller index z-axis> <surface layers> <lattice constant> <repeat in x> <repeat in y> <vacuum>

#example: ./create_slab.py Ag fcc 1 1 1 6 4.17 2 2 6.5 (Create a POSCAR file of an 2x2 Ag fcc(111) surface with 6 layers, a lattice constant of 4.17 Ang and 6.5 Ang vacuum)

# 2do: include general procedure for hcp metals?, convert orthorhombic to orthogonal lattice, 

import sys
from ase.build import bulk
from ase.build import surface
from ase.io import write
#from ase.visualize import view

miller_index = [] # x y z miller indices

symbol = sys.argv[1] # element symbol
lattice_type = sys.argv[2] # fcc, bcc, diamond
miller_index.append(sys.argv[3])
miller_index.append(sys.argv[4])
miller_index.append(sys.argv[5])
layer = sys.argv[6]
lattice_constant = sys.argv[7] # lattice constant in Ang
repeat_x = sys.argv[8]
repeat_y = sys.argv[9]
vacuum = sys.argv[10]

print("Create a POSCAR file of an {}x{} {} {}({}{}{}) slab with {} layer, a lattice constant of {} Ang and {} Ang vacuum!".format(repeat_x, repeat_y, symbol, lattice_type, miller_index[0], miller_index[1], miller_index[2], layer, lattice_constant, vacuum))

bulk = bulk(symbol, lattice_type, a=lattice_constant, cubic=True)
#view(bulk)

slab = surface(bulk, (int(miller_index[0]),int(miller_index[1]),int(miller_index[2])), int(layer), vacuum=int(vacuum))
#view(slab)
file_name = symbol+"_"+lattice_type+miller_index[0]+miller_index[1]+miller_index[2]+".POSCAR"

print(slab.cell)
print(len(slab))
for atom in slab:
    print(atom.symbol, atom.position)

slab = slab.repeat([int(repeat_x),int(repeat_y),1])
#view(slab)
write(file_name, slab, format='vasp')
