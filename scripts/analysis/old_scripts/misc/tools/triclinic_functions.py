from math import sqrt, cos, acos, radians, degrees

#A bunch of triclinic cell utility functions for converting to and from the dreaded LAMMPS format
#Ideally should not use numpy or ase, for portability
#More info, check section 8.3.2 of the How To in the LAMMPS manual: https://lammps.sandia.gov/doc/Howto_triclinic.html

#Assumes cell vectors are passed in as:
#[ax, ay, az]
#[bx, by, bz]
#[cx, cy, cz]

#LAMMPS requires the vectors to be:
#[ ax, 0.0, 0.0], with ax>0.0
#[ bx,  by, 0.0], with by>0.0
#[ cx,  cy,  cz], with cz>0.0
#This forms a right handed basis
#This requires some rotations and cell vector swaps sometimes (which also leads to atomic position being changed)

#Sometimes notations use the Transpose of this matrix, so be careful!

#An additional complication, is that LAMMPS uses a different notation for its data files, and for its trj files 9where the described cell is a "boudning box"

def modulus(vector):
    length=sqrt(sum([vi**2 for vi in vector]))
    return length

def angle(v, u):
    dot_product=dot(v, u)
    len_v, len_u=[modulus(v), modulus(u)]
    angle_vu=acos(dot_product/(len_v*len_u))
    return angle_vu

def dot(v, u):
    dot_product=sum([vi*ui for vi, ui in zip(v, u)])
    return dot_product

def cross(v, u):
    x, y, z=0, 1, 2 #components, for ease of indexing
    cross_vector=[v[y]*u[z]-v[z]*u[y],
                -(v[x]*u[z]-v[z]*u[x]),
                  v[x]*u[y]-v[y]*u[x]
                 ]

    test=[dot(v, cross_vector), dot(u, cross_vector)]
    if (test[0]>0.00001) or (test[1]>0.00001):
        print("ERROR: Cross product not orthogonal to original vectors, formula wrong")
        exit(1)
    return cross_vector

def is_lammps_format(cell_vectors):
    
    cv=cell_vectors
    a, b, c=0, 1, 2
    x, y, z=0, 1, 2

    #positive components, and zero components
    if ((cv[a][x]>0.0 and cv[b][y]>0.0 and cv[c][z]>0.0) and
        (abs(cv[a][y])<0.00001 and abs(cv[a][z])<0.00001 and abs(cv[b][z])<0.00001)):
        return True
    else:
        return False

def is_right_handed_basis(cell_vectors):
    #right handed if signed_V>0
    sv=signed_volume(cell_vectors)
    if sv>0.0:
        return True
    else:
        return False
    
def signed_volume(cell_vectors):
    #We need to calculate (v1Xv2).v3=dot(cross(a, b), c)
    a, b, c=cell_vectors
    return dot(cross(a, b), c)

def volume(cell_vectors):
    return abs(signed_volume(cell_vectors))

def get_cell_lengths_and_angles(cell_vectors):
    l_a, l_b, l_c=[modulus(vi) for vi in cell_vectors]

    cv=cell_vectors
    angles=[angle(cv[1], cv[2]), angle(cv[0], cv[2]), angle(cv[0], cv[1])]
    angle_bc, angle_ac, angle_ab=[degrees(ai) for ai in angles]

    return l_a, l_b, l_c, angle_bc, angle_ac, angle_ab

def check_triclinic(cell_vectors, tol=0.01):
    #check if we need to output the structure as a triclinic object in LAMMPS
    cell=get_cell_lengths_and_angles(cell_vectors) #[len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)]
    if abs(cell[3]-90.0)>tol or abs(cell[4]-90.0)>tol or abs(cell[5]-90.0)>tol:
        return True
    else:
        return False

def get_triclinic(cell_vectors):
    #return the necessary values for a triclinic output
    a, b, c, alpha, beta, gamma=get_cell_lengths_and_angles(cell_vectors) #[len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)]
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

def lammpstrjformat_to_cellvectors(xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound, xy, xz, yz, tol=0.000001):
    amod, bmod, cmod, cosalpha, cosbeta, cosgamma=lammpstrjformat_to_lengthsandangles(xlo_bound, xhi_bound, 
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

def lammpstrjformat_to_lengthsandangles(xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound, xy, xz, yz, tol=0.000001):

    xlo = xlo_bound - min(0.0,xy,xz,xy+xz)
    xhi = xhi_bound - max(0.0,xy,xz,xy+xz)
    ylo = ylo_bound - min(0.0,yz)
    yhi = yhi_bound - max(0.0,yz)
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

def lammpsdataformat_to_cellvectors():
    pass

def cellvectors_to_lammpsdataformat():
    pass

def cellvectors_to_lammpstrjformat(cell_vectors):
    #Check format is correct, otherwise we have to do some annoying transformations
    if not(is_right_handed_basis(cell_vectors)):
        print("ERROR: Cell vectors don't form a right handed basis, not prepared")
        exit(1)
    if not(is_lammps_format(cell_vectors)):
        print("ERROR: Cell vectors are not in lammps format, not prepared")
        exit(1)

    #From now on, format should be correct, just need to calculate some values

    xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz=get_triclinic(cell_vectors)

    xlo_bound = xlo + min(0.0,xy,xz,xy+xz)
    xhi_bound = xhi + max(0.0,xy,xz,xy+xz)
    ylo_bound = ylo + min(0.0,yz)
    yhi_bound = yhi + max(0.0,yz)
    zlo_bound = zlo
    zhi_bound = zhi

    return xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound, xy, xz, yz
