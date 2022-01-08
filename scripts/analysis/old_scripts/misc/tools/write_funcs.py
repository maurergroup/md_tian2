from math import sqrt, cos, radians

bohr2ang=0.529177
ang2bohr=1.0/0.529177

#LAMMPS manual page 6.13 describes how to go from 3 cell vectors to LAMMPS triclinic format, and from 3 lengths + 3 angles to LAMMPS triclinic format

def write_runner(filename, atom_object, comment="#", write_mode="w"):
    dump=open(filename, mode=write_mode)
    dump.write("begin\n")
    dump.write("comment Structure generated automatically by script \n")
    dump.write("comment"+" "+comment+"\n")

    pbc=atom_object.get_pbc()
    cell=atom_object.get_cell()
    for i in range(3):
            if (pbc[i] == True):
                    #cellstr=str(cell[i][0])+" "+str(cell[i][1])+" "+str(cell[i][2])
                    cellstr="{0[0]:12.9f} {0[1]:12.9f} {0[2]:12.9f}".format(ang2bohr*cell[i])
                    dump.write("lattice "+cellstr+"\n")

    for atom in atom_object:
            atomstr="{0[0]:12.9f} {0[1]:12.9f} {0[2]:12.9f} {1:2}".format(ang2bohr*atom.position, atom.symbol)
            dump.write("atom    "+atomstr+" 0.0 0.0 0.0 0.0 0.0\n")

    dump.write("energy 0.00\n")
    dump.write("charge 0.00\n")
    dump.write("end\n")
    dump.close()

def write_lammpstrj(filename, atom_object, write_mode="w", timestep=0):
    outfile=open(filename, mode=write_mode)

    triclinic=check_triclinic(atom_object)

    outfile.write("ITEM: TIMESTEP\n")
    outfile.write("{}\n".format(timestep))
    outfile.write("ITEM: NUMBER OF ATOMS\n")
    outfile.write("{}\n".format(len(atom_object)))
    if triclinic:
        xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz=get_triclinic(atom_object) #[xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz]
        bounds=[xlo+min([0.0, xy, xz, xy+xz]),
                xhi+max([0.0, xy, xz, xy+xz]),
                ylo+min([0.0, yz]),
                yhi+max([0.0, yz]),
                zlo,
                zhi]
        tilts=[xy, xz, yz]
        outfile.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
        for i in range(0, 3):
            outfile.write("{} {} {}\n".format(bounds[i*2], bounds[i*2+1], tilts[i]))
    else: #orthogonal normal box
        outfile.write("ITEM: BOX BOUNDS pp pp pp\n")
        for i in range(0, 3):
            outfile.write("0.0 {}\n".format(atom_object.get_cell()[i][i]))
    outfile.write("ITEM: ATOMS id element x y z\n")
    for i, atom in enumerate(atom_object):
        outfile.write("{0} {1} {2[0]} {2[1]} {2[2]}\n".format(i+1, atom.symbol, atom.position))
    
    outfile.close()


def write_lammps_data(filename, atom_object, comment="#", masses=[["H", 1.0], ], round_pos=False):
    outfile=open(filename, mode='w')

    triclinic=check_triclinic(atom_object)

    outfile.write("{}\n".format(comment))
    outfile.write("\n")

    outfile.write("{} atoms\n".format(len(atom_object)))
    outfile.write("\n")

    outfile.write("{} atom types\n".format(len(masses)))
    outfile.write("\n")

    if triclinic:
        params=get_triclinic(atom_object) #[xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz]
        bound_strings=[["xlo", "xhi"], ["ylo", "yhi"], ["zlo", "zhi"]]
        for i in range(0, 3):
            outfile.write("0.0 {} {} {}\n".format(params[1+2*i], bound_strings[i][0], bound_strings[i][1]))
        outfile.write("{} {} {} xy xz yz\n".format(params[6], params[7], params[8]))
        outfile.write("\n")

    else: #orthogonal normal box
        bound_strings=[["xlo", "xhi"], ["ylo", "yhi"], ["zlo", "zhi"]]
        for i in range(0, 3):
            outfile.write("0.0 {} {} {}\n".format(atom_object.get_cell()[i][i], bound_strings[i][0], bound_strings[i][1]))
        outfile.write("\n")

    outfile.write("Masses\n")
    outfile.write("\n")
    for i, mass_pair in enumerate(masses):
        outfile.write("{} {}\n".format(i+1, mass_pair[1]))
    outfile.write("\n")

    outfile.write("Atoms\n")
    outfile.write("\n")
    for i, atom in enumerate(atom_object):
        for j, elt in enumerate(masses):
            if elt[0]==atom.symbol:
                type_num=j+1
                break
        pos=atom.position
        for n, x in enumerate(pos): #avoid very small almost zero coordinates
            if round_pos:
                pos[n]=round(pos[n], 5)
            if abs(x)<0.000001:
                pos[n]=0.0
        outfile.write("{0} {1} {2[0]} {2[1]} {2[2]}\n".format(i+1, type_num, pos))
    
    outfile.close()

def check_triclinic(atom_object, tol=0.01):
    #check if we need to output the structure as a triclinic object in LAMMPS
    cell=atom_object.get_cell_lengths_and_angles() #[len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)]
    if abs(cell[3]-90.0)>tol or abs(cell[4]-90.0)>tol or abs(cell[5]-90.0)>tol:
        return True
    else:
        return False

def get_triclinic(atom_object):
    #return the necessary values for a triclinic output
    a, b, c, alpha, beta, gamma=atom_object.get_cell_lengths_and_angles() #[len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)]
    #print([a, b, c, alpha, beta, gamma])
    alpha=radians(alpha)
    beta=radians(beta)
    gamma=radians(gamma)
    xlo=ylo=zlo=0.0
    xhi=a
    xy=b*cos(gamma)
    xz=c*cos(beta)
    yhi=sqrt(b*b-xy*xy)
    yz=(b*c*cos(alpha)-xy*xz)/(yhi)
    zhi=sqrt(c*c-xz*xz-yz*yz)
    if abs(xy)<0.001:
        xy=0.0
    if abs(xz)<0.001:
        xz=0.0
    if abs(yz)<0.001:
        yz=0.0
    return [xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz]
