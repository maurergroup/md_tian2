#!/usr/bin/python/

import os, sys, math, copy, numpy, time

### edit here ###

#METAL_TYPE = "Pd"
#SHOT_THRU_LIMIT = -6.8256    # yes, this should be a negative number

#METAL_TYPE = "Au"
#SHOT_THRU_LIMIT = -7.27636607765

METAL_TYPE = "C"
SHOT_THRU_LIMIT = 0.0

SPECULAR_RADIUS = 1.5
### do not edit anything below this line ###

VERSION_ID = 2.07

class Point3D:
	def __init__(self,x,y,z):
		self.x = x
		self.y = y
		self.z = z

	def add(self,other):
		self.x += other.x
		self.y += other.y
		self.z += other.z
	
	def minus(self,other):
                self.x -= other.x
                self.y -= other.y
                self.z -= other.z
	
	def multiply(self, factor):
		self.x *= factor
		self.y *= factor
		self.z *= factor

	def distance(self, other):
		return math.sqrt( (self.x-other.x)**2 + \
			(self.y-other.y)**2 +(self.z-other.z)**2 )

	def normalize(self):
                l = self.length()
                self.x = self.x / l
                self.y = self.y / l
                self.z = self.z / l
	
	def project_onto_plane(self, plane):
    		d = self.dot_product(plane) / plane.length()
		plane.normalize()
    		plane.multiply(d)
    		self.minus(plane)


	def rotate_about_z(self, angle):
		rangle = math.radians(float(angle))
		self.x, self.y = math.cos(rangle)*self.x - math.sin(rangle)*self.y, math.sin(rangle)*self.x + math.cos(rangle)*self.y

	def rotate_about_y(self, angle):
                rangle = math.radians(float(angle))
                self.x, self.z = math.cos(rangle)*self.x + math.sin(rangle)*self.z, -math.sin(rangle)*self.x + math.cos(rangle)*self.z

	def dot_product(self, other):
                return self.x*other.x + self.y*other.y + self.z*other.z

	def cross_product(self, other):
		return [self.y*other.z - self.z*other.y, self.z*other.x - self.x*other.z, self.x*other.y - self.y*other.x]

	def cross_product2(self, other):
		cp = self.cross_product(other)
		return Point3D(cp[0], cp[1], cp[2])

	def length(self):
		return math.sqrt( self.x**2 + self.y**2 + self.z**2)
	
	def __str__(self):
		return str(self.x) + '  ' + str(self.y) + '  ' + str(self.z)

class Traj:


	def __init__(self, ekin_p_i, ekin_l_i, epot_i, etotal_i,\
			r_p_i, polar_i, azi_i, ekin_p_f, ekin_l_f, epot_f,   \
                        etotal_f, r_p_f, polar_f, azi_f, time, turn_pnts,    \
			cl_appr, r_p_min):
		self.ekin_p_i = ekin_p_i
		self.ekin_l_i = ekin_l_i
		self.epot_i   = epot_i
		self.etotal_i = etotal_i
		self.r_p_i    = r_p_i
		self.polar_i  = polar_i
		self.azi_i    = azi_i
		self.ekin_p_f = ekin_p_f
		self.ekin_l_f = ekin_l_f
		self.epot_f   = epot_f
		self.etotal_f = etotal_f
		self.r_p_f    = r_p_f
		self.polar_f  = polar_f
		self.azi_f    = azi_f
		self.time     = time
		self.turn_pnts = turn_pnts
		self.cl_appr  = cl_appr
		self.r_p_min  = r_p_min
		self.eloss    = ekin_p_i - ekin_p_f
		self.has_scattered = r_p_f.z > r_p_i.z
		self.has_transmitted = r_p_f.z < SHOT_THRU_LIMIT
		self.has_adsorbed = not (self.has_scattered or self.has_scattered)
		self.delta_azi = min(360-abs(azi_f-azi_i), abs(azi_f-azi_i))
		self.in_spec = math.sqrt( (polar_f-polar_i)**2 + (azi_f-azi_i)**2 ) < SPECULAR_RADIUS
		self.in_plane = self.delta_azi < SPECULAR_RADIUS

	def spec_scattering_vector(self):
		inc = copy.deepcopy(self.v_p_i)
                inc.z *= -1.
                inc.normalize()
		return inc

	def angle_with_spec_scattering_vector(self):
		inc = self.spec_scattering_vector()
		out = copy.deepcopy(self.v_p_f)
		out.normalize()
		dp = inc.dot_product(out)
		angle_in_rad = math.acos(dp)
		return math.degrees( angle_in_rad )

	def angle_with_vector(self, vec):
                out = copy.deepcopy(self.v_p_f)
                out.normalize()
		vec = Point3D(vec[0], vec[1], vec[2])
		vec.normalize()
                dp = out.dot_product(vec)
                angle_in_rad = math.acos(dp)
                return math.degrees( angle_in_rad )


def angle_between(fv, sv):
        fv1 = Point3D(fv[0], fv[1], fv[2])
	sv1 = Point3D(sv[0], sv[1], sv[2])
        fv1.normalize()
	sv1.normalize()
        dp = fv1.dot_product(sv1)
        angle_in_rad = math.acos(dp)
        return math.degrees( angle_in_rad )

def convert_index_to_float(line, idx):
	l = line.strip(' \n\t\r').split()
	return float(l[idx])

def convert_line_to_point3d(line):
	l = line.strip(' \n\t\r').split()
	return Point3D(float(l[0]), float(l[1]), float(l[2]))

def obtain_rot_mat(vec, angle):
	# vec describes line thru origin around which the rotation is performed
	# rotation angle needs to be given in degrees
        # normalize rotation vector
        vec.normalize()

        rangle = math.radians(float(angle))

        rot = [[None,None,None], [None,None,None], [None,None,None]]
        rot[0][0] = vec.x*vec.x*(1.-math.cos(rangle)) +        math.cos(rangle)
        rot[1][0] = vec.x*vec.y*(1.-math.cos(rangle)) - vec.z*math.sin(rangle)
        rot[2][0] = vec.x*vec.z*(1.-math.cos(rangle)) + vec.y*math.sin(rangle)

        rot[0][1] = vec.y*vec.x*(1.-math.cos(rangle)) + vec.z*math.sin(rangle)
        rot[1][1] = vec.y*vec.y*(1.-math.cos(rangle)) +        math.cos(rangle)
        rot[2][1] = vec.y*vec.z*(1.-math.cos(rangle)) - vec.x*math.sin(rangle)

        rot[0][2] = vec.z*vec.x*(1.-math.cos(rangle)) - vec.y*math.sin(rangle)
        rot[1][2] = vec.z*vec.y*(1.-math.cos(rangle)) + vec.x*math.sin(rangle)
        rot[2][2] = vec.z*vec.z*(1.-math.cos(rangle)) +        math.cos(rangle)
        return rot

def matmul(mat, vec):
	v1 = mat[0][0]*vec[0] + mat[0][1]*vec[1] + mat[0][2]*vec[2]
	v2 = mat[1][0]*vec[0] + mat[1][1]*vec[1] + mat[1][2]*vec[2]
	v3 = mat[2][0]*vec[0] + mat[2][1]*vec[1] + mat[2][2]*vec[2]
	return [v1, v2, v3]

def initialize():
	ntrajs = sum(1 for line in open("MXT2Summary.txt", "r")) -1 	# first line is commment
	print "Reading %d trajectories" % ntrajs
	traj_list = []					# init list
	mxt_file = open("MXT2Summary.txt", "r")
	counter = 0
	scattered = 0
	absorbed = 0
	transmitted = 0
	for line in mxt_file:
		if line.startswith("#"):				# is comment line
			continue
		if (counter % (ntrajs/10) == 0):
                        print (100*counter+1)/ntrajs, "%"
		sl = line.strip("\n\t\r ").split()
		ekin_p_i = float(sl[0])						#	e_kin_p_i = 3.33000
		ekin_l_i = float(sl[1])						#	e_kin_l_i = 4.07264
		epot_i   = float(sl[2])						#	epot_i    = 30.09554
		etotal_i = float(sl[3])						#       e_total_i = 37.49818
		r_p_i    = Point3D(float(sl[4]), float(sl[5]), float(sl[6]))	#	r_i       = 14.88455  -2.57611  6.00000
		polar_i  = float(sl[7])                                         #       polar_i   = 50.00000
		azi_i    = float(sl[8])                                         #       azi_i     = 0.000000

		ekin_p_f = float(sl[9])						#	ekin_p_f  = 0.05978
		ekin_l_f = float(sl[10])					#	ekin_l_f  = 5.26957
		epot_f   = float(sl[11])					#	epot_f    = 28.73945
		etotal_f = float(sl[12])					#	etotal_f  = 34.06880
		r_p_f    = Point3D(float(sl[13]), float(sl[14]), float(sl[15]))	#	r_f       = 13.70619  1.33464  -1.07431
		polar_f  = float(sl[16])                                        #       polar_f   = 27.23457
                azi_f    = float(sl[17])                                        # 	azi_f     = 4.23478	
		
		time     = float(sl[18])					#	time      = 978.70000
		turn_pnts = int(sl[19])						#       turn_pnts = 14
		cl_appr  = float(sl[20])					#	cl_appr  =      0.9846501
		r_p_min  = Point3D(float(sl[21]), float(sl[22]), float(sl[23])) #	r_min_p  = 33.4630699  31.9529501  0.9322836
	
		this_traj = Traj(ekin_p_i, ekin_l_i, epot_i, etotal_i, r_p_i, polar_i, azi_i, ekin_p_f, \
					ekin_l_f, epot_f, etotal_f, r_p_f, polar_f, azi_f, time, turn_pnts, \
					cl_appr, r_p_min)

		traj_list.append(this_traj)
		counter += 1

		if this_traj.ekin_p_f > 1.4*this_traj.ekin_p_i:
			print "Warning: a projectile gained more than 40% of its initial kinetic energy", this_traj.ekin_p_f


	mxt_file.close()

	if not os.path.exists("analysis"):
    		os.makedirs("analysis")

	for traj in traj_list:
		if traj.has_scattered:
			scattered += 1
		elif traj.has_transmitted:
			transmitted += 1
		else:
			absorbed += 1

	return traj_list, scattered, absorbed, transmitted



def analyze(trajs):

        ### ANGULAR DISTRIBUTION ###
	print "Calculating angular energy loss"
	# get trajectories that are within specular radius in azimuth direction
	energy_collect  = []
	angle_collect   = []
	
	for traj in trajs:
		einc = traj.ekin_p_i
		inc_angle = traj.polar_i
		if traj.has_scattered and traj.in_plane:
			energy_collect.append(traj.ekin_p_f)
			angle_collect.append(traj.polar_f)

	ang_dist_file = open("analysis/2D_matrix.txt", "w")

	angle_step = 2.5
	energy_step = 0.04

	ang_dist_file.write("# Columns represent scattering angles. First column is 0 degrees and it continues in steps of %f degrees\n" % angle_step)
	ang_dist_file.write("# Rows represent final kinetic energy. First row is 0 eV and in continues in steps of %f eV.\n" % energy_step)
	ang_dist_file.write("# An entry means that this probability lies between this_value and this_value+step.\n")
	ang_dist_file.write("# Incidence energy: %f eV    incidence angle: %f degrees    detector radius %f degrees\n" % (einc, inc_angle, SPECULAR_RADIUS))
	ang_dist_file.write("#\n")
	
	

	if len(energy_collect) != 0:
		angle_eloss_hist, xedges, yedges = numpy.histogram2d(energy_collect, angle_collect,  bins=(numpy.arange(0, 2.51, energy_step), numpy.arange(0, 90, angle_step)) , normed=True)

		# OUTPUT

		ang_dist_file.write("0.000000  ")

		for i in range(len(yedges)-1):
			if (0.5*(yedges[i]+yedges[i+1])) < 10:
				ang_dist_file.write("%f  " % (0.5*(yedges[i]+yedges[i+1])))
			else:
				ang_dist_file.write("%6.5f  " % (0.5*(yedges[i]+yedges[i+1])))
		ang_dist_file.write("\n")

		for i in range(len(angle_eloss_hist)):
			ang_dist_file.write("%f  " % (0.5*(xedges[i]+xedges[i+1])/einc))
			for j in range(len(angle_eloss_hist[i])):
				ang_dist_file.write("%f  " % angle_eloss_hist[i][j])
			ang_dist_file.write("\n")
	else:
		ang_dist_file.write("%f %f %d\n"   % (0.1, 0.5, 0))
                ang_dist_file.write("%f %f %d\n\n" % (0.1, 1.0, 1))
		ang_dist_file.write("%f %f %d\n"   % (0.2, 1.0, 1))
		ang_dist_file.write("%f %f %d\n"   % (0.2, 0.5, 2))
		
	ang_dist_file.close()



	
	
	

# Average energy of H-atoms in bulk 0.03969474585 eV.
# Average energy of H-atoms reflected 1.93905890836 eV.
# Peak energy of H-atoms reflected 2.770000 eV
#


### READ IN TRAJS ###
traj_collection, SCATTERED, ABSORBED, TRANSMITTED = initialize()

### CALCULATE USEFUL CONSTANTS ###
NTRAJS = len(traj_collection)
FRAC_SCATTERED = float(SCATTERED)/NTRAJS
FRAC_ABSORBED = float(ABSORBED)/NTRAJS
FRAC_TRANSMITTED = float(TRANSMITTED)/NTRAJS

### OUTPUT ###
analyze(traj_collection)

