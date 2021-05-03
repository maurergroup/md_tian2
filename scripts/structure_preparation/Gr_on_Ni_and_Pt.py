#!/usr/bin/env python3

# intention: create input structures for graphene (Gr) on Ni and Pt

# use: python3 <scriptname>

#2do: add needed comments, mention paper with algo for lattice mismatch, include MLP w/ Mail as reference, mention/show how to come up with lattice vectors for super structure that respects both lattice constants

from ase.lattice.hexagonal import Graphene #for creating the Gr slab
from ase.build import fcc111 #for creating Ni fcc(111) surface
from ase.io import write #for writing to POSCAR
from ase.build import add_adsorbate  #for putting adsorbates Gr and H to slab!
from ase.build import sort #
from ase.build import surface
from ase import Atoms, Atom # add H
from ase.visualize import view # view the structures
import numpy as np
from math import sqrt, radians, sin, cos
from numpy.linalg import det

# class and functions

class create_layer():
    def __init__(self, structure, rotation, l, lprime, c1wave, c2wave):
        self.structure=Atoms(pbc=[1, 1, 1])
        self.original_cell=structure.get_cell()
        self.basis=get_basis(structure)
        self.rotation=rotation
        self.l=l
        self.lprime=lprime
        self.c1wave=c1wave
        self.c2wave=c2wave

    def make_cell(self, mark_extremes=False):
        l=self.l
        lprime=self.lprime
        lpluslprime=l+lprime
        origin=np.array([0, 0])
        for l1 in range(min(l[0], lprime[0], lpluslprime[0], origin[0]), max(l[0], lprime[0], lpluslprime[0], origin[0])+1):
            for l2 in range(min(l[1], lprime[1], lpluslprime[1], origin[1]), max(l[1], lprime[1], lpluslprime[1], origin[1])+1):
                if (l1==l[0] and l2==l[1]) or (l1==lprime[0] and l2==lprime[1]) or (l1==lpluslprime[0] and l2==lpluslprime[1]): #avoid the endpoints
                    continue
                if is_inside_parallelogram(l, lprime, np.array([l1, l2])):
                    lattice_pos=np.dot(self.original_cell.transpose(), np.array([l1, l2, 0]))
                    for atom in self.basis:
                        pos=lattice_pos+atom[1] #for eacha tom in the basis
                        self.structure.append(Atom(atom[0], pos))

        #self.structure.rotate(self.rotation, v="-z", rotate_cell=True)
        self.structure.set_cell([self.c1wave, self.c2wave, self.original_cell[2]]) #keep the z vector
        self.structure.rotate(self.rotation, v="-z", rotate_cell=True)

        #self.structure.wrap(pbc=[1, 1, 0]) #don"t wrap in the z direction

        #mark_extremes=True
        if mark_extremes:
            mark_atom="Li"
            lattice_pos=np.dot(self.original_cell.transpose(), np.array([  0,    0, 0]))
            self.structure.append(Atom(mark_atom, rotate(self.rotation, lattice_pos)))
            lattice_pos=np.dot(self.original_cell.transpose(), np.array([l[0], l[1], 0]))
            self.structure.append(Atom(mark_atom, rotate(self.rotation, lattice_pos)))
            lattice_pos=np.dot(self.original_cell.transpose(), np.array([lprime[0], lprime[1], 0]))
            self.structure.append(Atom(mark_atom, rotate(self.rotation, lattice_pos)))
            lattice_pos=np.dot(self.original_cell.transpose(), np.array([lpluslprime[0], lpluslprime[1], 0]))
            self.structure.append(Atom(mark_atom, rotate(self.rotation, lattice_pos)))
        
        #self.structure.rotate(self.rotation, v="-z", rotate_cell=True)
        #self.structure.set_cell([self.c1wave, self.c2wave, self.original_cell[2]]) #keep the z vector
        #self.structure.wrap(pbc=[1, 1, 0]) #don"t wrap in the z direction

    def expand(self, newc1wave, newc2wave):
        self.structure.set_cell([newc1wave, 
                                 newc2wave, 
                                 self.structure.get_cell()[2]],
                                scale_atoms=True)

def is_inside_parallelogram(A, B, C, thr=0.0001):
    n= det(np.array([A, C]))/det(np.array([A, B]))
    m=-det(np.array([B, C]))/det(np.array([A, B]))
    #print(n, m)
    #m= det(B, C)/det(A,B)
    if m<0.0 or m>(1.0-thr) or n<0.0 or n>(1.0-thr):
        return False
    else:
        return True

def platinum_cell(a=3.9242, miller=[1, 1, 1], vac=0.0, nlayers=1, rep=[1, 1, 1]):
    v1=np.array([0.50,            0.50,   0.00])
    v2=np.array([0.00,            0.50,   0.50])
    v3=np.array([0.50,            0.00,   0.50])

    v1=v1*a
    v2=v2*a
    v3=v3*a

    Pt = Atoms('Pt',
                pbc=(1,1,1),
                scaled_positions=[(0.000,0.000,0.000)],
                cell=[v1,v2,v3]
               )
    Pt=surface(Pt, indices=miller, layers=nlayers, vacuum=vac)
    Pt=Pt.repeat(rep)
    return Pt

def get_basis(atoms_object):
    basis_atoms=[]
    for atom in atoms_object:
        basis_atoms.append([atom.symbol, atom.position])
    return basis_atoms

def displace(atoms_object, displ=np.array([0.0, 0.0, 0.0])):
    ret_obj=atoms_object.copy()
    for atom in ret_obj:
        atom.position+=displ
    return ret_obj

def rotate(angle, vector):
    rad_angle=radians(angle)
    Mrot=np.array([[ cos(rad_angle), sin(rad_angle), 0.0], 
                   [-sin(rad_angle), cos(rad_angle), 0.0], 
                   [            0.0,            0.0, 0.0]])
    rot_vector=np.dot(Mrot, vector.transpose()).transpose()
    return rot_vector

# set variables
max_layer = 7
Ni_vacuum = 10 # 9+?
Pt_vacuum = 12 # 10+?
Pt_rotation = 0.0

# get lattice parameters of Gr on Pt super cell
aprime_coeff=np.array([1,  1])
bprime_coeff=np.array([1, -2])
platinum_cell_parameter=platinum_cell(nlayers=1)
a,b,c=platinum_cell_parameter.get_cell()
aprime=aprime_coeff[0]*a+aprime_coeff[1]*b
bprime=bprime_coeff[0]*a+bprime_coeff[1]*b

# parameters for Gr
CCdist=1.42 #carbon carbon distance in Ang
graphene_a=sqrt(3)*CCdist #definition of Gr lattice constant (see paper, geometrical conditions)
graphene_c=graphene_a # geometrical condition for Gr

# arrays for different metal layer
Ni_layer=[]
Pt_layer=[]
Pt_super_layer=[]
layer_list = list(range(max_layer))

# create Gr for Ni and Gr for Pt super structure
Gr = Graphene(symbol='C',latticeconstant=(graphene_a, graphene_c)) # building block has 2 atoms in x

# Pt:
Gr_super = Gr.repeat([2,2,1]) # 2x2 super structure for Pt (8 atoms total)

# Ni:
a,b,c=Gr.get_cell()
Gr = Gr.repeat([1,2,1])
Gr.set_cell([a, 2*b+a, c], scale_atoms=False) # orthorhombic to orthogonal conversion (only works for 60 deg)
Gr.wrap()
Gr.rotate(60, 'z') # without rotation the hcp sites are occupied; now the fcc sites are occupied by Gr
Gr.wrap()
Gr = Gr.repeat([2,1,1]) # 8 atoms total

# generate first type Ni slabs (Gr on fcc Ni):
for layer,slab in enumerate(layer_list): # loop over Ni layer
    slab = fcc111('Ni', (2,2,layer+1), orthogonal=True, vacuum=Ni_vacuum)
    Ni_cell = slab.get_cell()
    Gr_cell  = Gr.get_cell()
    Gr.set_cell([Ni_cell[0], Ni_cell[1], Gr_cell[2]], scale_atoms=True) # adjust metal lattice vectors to Gr (bulk less flexible)
    add_adsorbate(slab, Gr, 2.1, 'ontop') # ontop, bridge, fcc, hcp possible
    add_adsorbate(slab, 'H', 6.1, 'ontop')
    slab.center()
    Ni_tags=slab.numbers
    slab=sort(slab, tags=Ni_tags) # required when the element sorting should be according to the mass (PES), otherwise alphabetical
    print(slab.cell)
    print(len(slab))
    for i,atom in enumerate(slab):
        print(atom.symbol, atom.position)
    print(slab.numbers)
    Ni_layer.append(slab)
    write("Ni_"+str(layer+1)+"l"+".POSCAR", slab, format='vasp') # vasp, vasp-out, vasp-xdatcar, vasp-xml possible

#view(Ni_layer)


# generate Pt superstructure slabs:
for layer,slab in enumerate(layer_list): # loop over Pt layer
    Pt_layer=create_layer(platinum_cell(nlayers=layer+1, vac=Pt_vacuum), Pt_rotation, aprime_coeff, bprime_coeff, aprime, bprime)
    Pt_layer.make_cell()
    Pt_layer.structure.wrap()
    slab=Pt_layer.structure
    slab.rotate(180, 'y', center=(0,0,0), rotate_cell=True) # needed so that the surface layer stays the same (add layer top-down); without bottom-up
    Pt_cell = slab.get_cell()
    Gr_super_rot = Gr_super.copy() # without, every time it will be rotated again!
    #Gr_super_temp.rotate(66.6, 'z', center=(0,0,0), rotate_cell=True)
    #Gr_super_temp.wrap()
    Gr_super_cell = Gr_super_rot.get_cell()
    Gr_super_rot.set_cell([Pt_cell[0], Pt_cell[1], Gr_super_cell[2]], scale_atoms=True) # adjust super structure lattice vectors to Gr
    Gr_super_rot.wrap()
    #Gr_super_rot.rotate(13.5, 'z', rotate_cell=True)
    #Gr_super_rot.wrap()
    #view(Gr_super_temp)
    #exit()
    
    add_adsorbate(slab, Gr_super_rot, 3.3, (0,0)) # for Pt maybe extend() is better for combining when to rotate and translate the Gr sheet
    slab.center()
    a,b,c = slab.get_cell() 
    slab = slab.repeat([1,2,1]) # crucial to make cell and pbc right
    slab.set_cell([a, 2*b+a, c], scale_atoms=False) # convert orthorhombic to orthogonal cell, result in doubling the atoms in y direction (therefore repeat is crucial)
    slab.wrap()
    add_adsorbate(slab, 'H', 7.3, (0,0))
    Pt_tags=slab.numbers
    slab=sort(slab, tags=Pt_tags)
    print(slab.cell)
    print(len(slab))
    for i,atom in enumerate(slab):
        print(atom.symbol, atom.position)
    print(slab.numbers)
    Pt_super_layer.append(slab)
    write("Pt_"+str(layer+1)+"l"+".POSCAR", slab, format='vasp')

#view(Pt_super_layer)
