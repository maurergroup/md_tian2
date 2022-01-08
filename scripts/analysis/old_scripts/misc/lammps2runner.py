import sys
import os
from math import sqrt

# intention: Convert LAMMPS trajectory to RuNNer input.data format
# use: python <scriptname> <LAMMPS input>

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

def lammpsformat_to_cellvectors(xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound, xy, xz, yz, tol=0.000001):
    amod, bmod, cmod, cosalpha, cosbeta, cosgamma=lammpsformat_to_lengthsandangles(xlo_bound, xhi_bound, 
                                                                                   ylo_bound, yhi_bound, 
                                                                                   zlo_bound, zhi_bound, 
                                                                                   xy, xz, yz, tol=0.000001)
    #print(amod, bmod, cmod, cosalpha, cosbeta, cosgamma)
    return lengthsandangles_to_cellvectors(amod, bmod, cmod, cosalpha, cosbeta, cosgamma)

def get_displ_vector(xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound, xy, xz, yz, tol=0.000001):
    xlo = xlo_bound - min(0.0,xy,xz,xy+xz)
    #xhi = xhi_bound - max(0.0,xy,xz,xy+xz)
    ylo =ylo_bound - min(0.0,yz)
    #yhi =yhi_bound - max(0.0,yz)
    zlo = zlo_bound
    #zhi = zhi_bound
    displ=[-xlo, -ylo, -zlo]
    for i in range(0, 3):
        if abs(displ[i])<tol:
            displ[i]=0.0
    return displ

def lammpsformat_to_lengthsandangles(xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound, xy, xz, yz, tol=0.000001):

    xlo = xlo_bound - min(0.0,xy,xz,xy+xz)
    xhi = xhi_bound - max(0.0,xy,xz,xy+xz)
    ylo =ylo_bound - min(0.0,yz)
    yhi =yhi_bound - max(0.0,yz)
    zlo = zlo_bound
    zhi = zhi_bound

    if abs(xlo)<tol:
        xlo=0.0
    if abs(ylo)<tol:
        ylo=0.0
    if abs(zlo)<tol:
        zlo=0.0

    if abs(xy)<tol:
        xy=0.0
    if abs(xz)<tol:
        xz=0.0
    if abs(yz)<tol:
        yz=0.0

    lx=xhi-xlo
    ly=yhi-ylo
    lz=zhi-zlo

    amod=lx
    bmod=sqrt(ly*ly+xy*xy)
    cmod=sqrt(lz*lz+xz*xz+yz*yz)

    cosalpha=(xy*xz+ly*yz)/(bmod*cmod)
    cosbeta =xz/cmod
    cosgamma=xy/bmod

    return amod, bmod, cmod, cosalpha, cosbeta, cosgamma

def lengthsandangles_to_cellvectors(amod, bmod, cmod, cosalpha, cosbeta, cosgamma, tol=0.000001):
    a=[0.0, 0.0, 0.0]
    b=[0.0, 0.0, 0.0]
    c=[0.0, 0.0, 0.0]

    #A
    ax=amod
    a[0]=ax

    #B
    bx=cosgamma*amod*bmod/ax
    if abs(bx)<tol:
        bx=0.0
    #print(bx)
    by=sqrt(bmod*bmod-bx*bx)
    b[0]=bx
    b[1]=by

    #C
    cx=cosbeta*amod*cmod/ax
    if abs(cx)<tol:
        cx=0.0
    cy=(cosalpha*bmod*cmod-bx*cx)/by
    if abs(cy)<tol:
        cy=0.0
    cz=sqrt(cmod*cmod-cx*cx-cy*cy)
    c[0]=cx
    c[1]=cy
    c[2]=cz

    return a, b, c

#Unit Conversion
bohr2angst=0.529177
angst2bohr=1.0/bohr2angst
#angst2bohr=1.0

infile_name=sys.argv[1]
outfile_name=infile_name.split(".")[0]+".runner"
infile=open(infile_name, mode="r")
outfile=open(outfile_name, mode="w")

#process file
checked_positions=False
for line in infile:
    spline=line.split()
    if spline[1]=="TIMESTEP":
        ts=next(infile).split()[0]
        outfile.write("begin\n")
        outfile.write("comment from file {} timestep {}\n".format(infile_name, ts))
        #outfile.write("comment files in {} user {} at host {}\n".format(os.getcwd(), os.environ['USER'], os.environ['HOSTNAME']))
        outfile.write("comment files in {} user {}\n".format(os.getcwd(), os.environ['USER']))

    elif spline[1]=="NUMBER":
        tot_natoms=int(next(infile))
    elif spline[1]=="BOX": #WARNING: This only works for orthogonal simulation boxes, which start at 0.0,0.0,0.0 
        box=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        #xhi, xlo, yhi, ylo, zhi, zlo, xy, xz, yz=0.0
        ilobound=[0.0, 0.0, 0.0] #not directly xlo, xhi, etc; but "bounds" according to the LAMMPS manual
        ihibound=[0.0, 0.0, 0.0]
        tilt=[0.0, 0.0, 0.0]
        triclinic=False
        if ("xy" or "xz" or "yz") in spline:
            triclinic=True
        for i in range(0, 3):
            spline=next(infile).split()
            ilobound[i]=float(spline[0])#*angst2bohr
            ihibound[i]=float(spline[1])#*angst2bohr
            if triclinic:
                tilt[i]=float(spline[2])
        box[0], box[1], box[2]=lammpsformat_to_cellvectors(ilobound[0], ihibound[0], 
                                                           ilobound[1], ihibound[1], 
                                                           ilobound[2], ihibound[2], 
                                                           tilt[0], tilt[1], tilt[2])
        #get displacement if ilo's are not zero
        displ=get_displ_vector(ilobound[0], ihibound[0],
                               ilobound[1], ihibound[1],
                               ilobound[2], ihibound[2],
                               tilt[0], tilt[1], tilt[2])
        displ=[di*angst2bohr for di in displ]

        for i in range(0, 3):
            box[i]=[bi*angst2bohr for bi in box[i]]
            outfile.write("lattice {0[0]:12.9f} {0[1]:12.9f} {0[2]:12.9f}\n".format(box[i]))

    elif spline[1]=="ATOMS":

        if checked_positions==False: #check in which position the element, x y and z tags are
            pos_x=spline.index("x")-2
            pos_y=spline.index("y")-2
            pos_z=spline.index("z")-2
            pos_el=spline.index("element")-2

        #loop thru all the atoms
        for i in range(1, tot_natoms+1):
            spline=next(infile).split()
            data=[angst2bohr*float(spline[pos_x])+displ[0], 
                  angst2bohr*float(spline[pos_y])+displ[1], 
                  angst2bohr*float(spline[pos_z])+displ[2], 
                  spline[pos_el]]
            #atom x y z element q e_atom fx fy fz
            outfile.write("atom {0[0]:12.9f} {0[1]:12.9f} {0[2]:12.9f} {0[3]:2} 0.0 0.0 0.0 0.0 0.0\n".format(data))

        #close the current runner block
        outfile.write("energy 0.0\n")
        outfile.write("charge 0.0\n")
        outfile.write("end\n")

infile.close()
outfile.close()
