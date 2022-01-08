#!/bin/env python

import os

path = os.getcwd()

file_names = os.listdir(path)

for name in file_names:
    if name.endswith('.lammpstrj'):
        infile_name = name # get full name and use first part with _# as file name
        outfile_name_start=infile_name.split(".")[0]


#infile_name  = "CH.lammpstrj" # loop over files in folder -> if entry.endswith(.lammpstrj) => get file name!
outfile_files = [] # array, assign names when number is known
timesteps = []
number_of_atoms = []
lattice_vectors = []
atoms = []



with open(infile_name, "r") as input_file:

    for line in input_file:

        if line.startswith("ITEM:"):

            if "TIMESTEP" in line:
                timesteps.append(next(input_file, '').strip())

                #print("Timestep:", timesteps)

            elif "NUMBER OF ATOMS" in line:
                number_of_atoms.append(next(input_file, '').strip())

                #print("Number of atoms:", number_of_atoms)

            elif "BOX BOUNDS pp pp pp" in line:
                lattice_vectors.append(next(input_file, '').strip())
                lattice_vectors.append(next(input_file, '').strip())
                lattice_vectors.append(next(input_file, '').strip())

                #print("Lattice vectors:", lattice_vectors)

            elif "ATOMS id element x y z" in line:

                for atom_id in range(int(number_of_atoms[-1])):
                    atoms.append(next(input_file, '').strip())

                #print("Atoms:", atoms)
                

print("timesteps:", len(timesteps))
print("number of atoms entries:", len(number_of_atoms))
print("lattice entries:", len(lattice_vectors))
print("atoms entries:", len(atoms))
