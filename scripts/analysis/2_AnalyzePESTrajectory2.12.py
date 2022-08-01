#!/usr/bin/env python3

# Original version 2.09 was done by Marvin Kammler
# Version 2.10 was modified by Sebastian Wille
# Version 2.11 was modified by Sebastian Wille
# Version 2.12 was modified by Nils Hertl

# Modifications v2.10:
# traj_id added in Traj class
# logfile will be written
# set readfile name
# fast and slow component analysis added (H@Gr)

# Modifications v2.11:
# velocities added in Traj class
# velocity analysis files will be written
# logfile entries added for velocities

# Modifications v2.12:
# angle correction added in polar plots

# intention: analyze compressed traj file to generate data files needed for plotting

# use like: python3 <scriptname> or ./<scriptname> # check first command line

import os, sys, math, copy, numpy, time

### edit here ###

#METAL_TYPE = "Pd"
#SHOT_THRU_LIMIT = -6.8256    # yes, this should be a negative number

#METAL_TYPE = "Au"
#SHOT_THRU_LIMIT = -7.27636607765

#METAL_TYPE = "C"
#SHOT_THRU_LIMIT = 0.0 

#SPECULAR_RADIUS = 1.5 # should match experimental settings (RAT detector radius)
#ION_IMAGING_AZI = 5   # should match experimental settings (ion imaging detector settings)
#AZIMUTHAL_ANGLE = 42  # should match experimental settings (ion imaging surface shift according to LEED)

#ANGLE_MAX = 90  # maximum angle in degrees
#ANGLE_MIN = -90 # minimum angle in degrees

#BINS = int((ANGLE_MAX-ANGLE_MIN)/5.0) # Denominator defining bin width, default is 2.5

inpname     = "MXT2Summary.txt" # output file
logfilename = "AnalyzePESTrajectory.log" # logfile
mdtinpname  = "md_tian.inp" # md_tian.inp file
settingname = "plot_settings.dat" # settings file

# H@Gr related
defnrgname = "deformation_energy_trajids.txt"


### CHANGES BELOW THIS LINE DEVELOPERS ONLY###

VERSION_ID = 2.12

### READ SETTINGS ###

settingfile = open(settingname, 'r')

for line in settingfile:
    if "METAL_TYPE" in line:
        METAL_TYPE      = str(line.split()[-1])
    if "SHOT_THRU_LIMIT" in line:
        SHOT_THRU_LIMIT = float(line.split()[-1])
    if "SPECULAR_RADIUS" in line:
        SPECULAR_RADIUS = float(line.split()[-1])
    if "ION_IMAGING_AZI" in line:
        ION_IMAGING_AZI = float(line.split()[-1])
    if "AZIMUTHAL_ANGLE" in line:
        AZIMUTHAL_ANGLE = float(line.split()[-1])
    if "ANGLE_MAX" in line:
        ANGLE_MAX       = int(line.split()[-1])
    if "ANGLE_MIN" in line:
        ANGLE_MIN       = int(line.split()[-1])

BINS = int((ANGLE_MAX-ANGLE_MIN)/5.0) # Denominator defining bin width, default is 2.5

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
			r_p_i, v_p_i, polar_i, azi_i, ekin_p_f, ekin_l_f, epot_f,   \
                        etotal_f, r_p_f, v_p_f, polar_f, azi_f, time, turn_pnts,    \
			cl_appr, cl_appr_t, r_p_min, traj_id):
                self.ekin_p_i   = ekin_p_i
                self.ekin_l_i   = ekin_l_i
                self.epot_i     = epot_i
                self.etotal_i   = etotal_i
                self.r_p_i      = r_p_i
                self.v_p_i      = v_p_i
                self.polar_i    = polar_i
                self.azi_i      = azi_i
                self.ekin_p_f   = ekin_p_f
                self.ekin_l_f   = ekin_l_f
                self.epot_f     = epot_f
                self.etotal_f   = etotal_f
                self.r_p_f      = r_p_f
                self.v_p_f      = v_p_f
                self.polar_f    = polar_f
                self.azi_f      = azi_f
                self.time       = time
                self.turn_pnts  = turn_pnts
                self.cl_appr    = cl_appr
                self.cl_appr_t  = cl_appr_t
                self.r_p_min    = r_p_min
                self.traj_id    = traj_id
                self.eloss      = ekin_p_i - ekin_p_f
                self.efrac      = ekin_p_f / ekin_p_i
                self.vloss      = length(v_p_i) - length(v_p_f)
                self.has_scattered = r_p_f.z > r_p_i.z
                self.has_transmitted = r_p_f.z < SHOT_THRU_LIMIT
                self.has_adsorbed = not (self.has_scattered or self.has_transmitted)
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

def length(self):
	return math.sqrt( (self.x)**2 + (self.y)**2 +(self.z)**2 )

def initialize(inpname,logfile):
        ntrajs = sum(1 for line in open(inpname, "r")) -1 	# first line is commment
        print("Reading {} trajectories".format(ntrajs))
        logfile.write("Reading {} trajectories\n".format(ntrajs))
        traj_list = []					# init list
        inp_file = open(inpname, "r")
        counter = 0
        scattered = 0
        #sc_count = 0
        absorbed = 0
        transmitted = 0
        for line in inp_file:
                if line.startswith("#"):				# is comment line
                	continue
                if (counter % (ntrajs/10) == 0):
                        print("{}%".format(100*counter/ntrajs+1))
                        logfile.write("{}%\n".format(100*counter/ntrajs+1))
                sl = line.strip("\n\t\r ").split()
                traj_id  = str(sl[0])                                           #       traj_id   = 00000001
                ekin_p_i = float(sl[1])						#	e_kin_p_i = 3.33000
                ekin_l_i = float(sl[2])						#	e_kin_l_i = 4.07264
                epot_i   = float(sl[3])						#	epot_i    = 30.09554
                etotal_i = float(sl[4])						#       e_total_i = 37.49818
                r_p_i    = Point3D(float(sl[5]), float(sl[6]), float(sl[7]))	#	r_i       = 14.88455  -2.57611  6.00000
                v_p_i    = Point3D(float(sl[8]), float(sl[9]), float(sl[10]))	#	r_i       = 14.88455  -2.57611  6.00000
                polar_i  = float(sl[11])                                         #       polar_i   = 50.00000
                azi_i    = float(sl[12])                                         #       azi_i     = 0.000000

                ekin_p_f = float(sl[13])					#	ekin_p_f  = 0.05978
                ekin_l_f = float(sl[14])					#	ekin_l_f  = 5.26957
                epot_f   = float(sl[15])					#	epot_f    = 28.73945
                etotal_f = float(sl[16])					#	etotal_f  = 34.06880
                r_p_f    = Point3D(float(sl[17]), float(sl[18]), float(sl[19]))	#	r_f       = 13.70619  1.33464  -1.07431
                v_p_f    = Point3D(float(sl[20]), float(sl[21]), float(sl[22]))	#	r_f       = 13.70619  1.33464  -1.07431
                polar_f  = float(sl[23])                                        #       polar_f   = 27.23457
                azi_f    = float(sl[24])                                        # 	azi_f     = 4.23478	
                
                time     = float(sl[25])					#	time      = 978.70000
                turn_pnts = int(sl[26])						#       turn_pnts = 14
                cl_appr  = float(sl[27])					#	cl_appr  =      0.9846501
                cl_appr_t  = float(sl[28])					#	cl_appr_t = 128
                r_p_min  = Point3D(float(sl[29]), float(sl[30]), float(sl[31])) #	r_min_p  = 33.4630699  31.9529501  0.9322836
        
                this_traj = Traj(ekin_p_i, ekin_l_i, epot_i, etotal_i, r_p_i, v_p_i, polar_i, azi_i, ekin_p_f, \
                			ekin_l_f, epot_f, etotal_f, r_p_f, v_p_f, polar_f, azi_f, time, turn_pnts, \
                			cl_appr, cl_appr_t, r_p_min, traj_id)

                traj_list.append(this_traj)
                #if counter == 0:
                #        einc = ekin_p_i
                #        print("Einc_set = {}".format(einc))
                counter += 1

                if this_traj.ekin_p_f > 1.4*this_traj.ekin_p_i:
                	print("Warning in traj {}: a projectile with final kinetic energy of {} gained more \
                            than 40% of its initial kinetic energy!".format(this_traj.traj_id,this_traj.ekin_p_f))
                	logfile.write("Warning in traj {}: a projectile with final kinetic energy of {} gained \
                            more than 40% of its initial kinetic energy!\n".format(this_traj.traj_id,this_traj.ekin_p_f))

                #if this_traj.has_transmitted is True:
                #        print("Particle was transmitted in traj {}".format(this_traj.traj_id))
                #        logfile.write("Particle was transmitted in traj {}".format(this_traj.traj_id))





        inp_file.close()

        if not os.path.exists("analysis"):
        	os.makedirs("analysis")

        for traj in traj_list:
                if traj.has_scattered:
                        scattered += 1
                        #if traj.cl_appr < 1.4:
                                #if traj.turn_pnts == 1:
                                        #def_nrg_file.write(traj.traj_id)
                                        #def_nrg_file.write(' ')
                                #print("Analyze slow component in traj {}".format(traj.traj_id))
                                #logfile.write("Analyze slow component in traj {}\n".format(traj.traj_id))
                                #sc_count += 1
                elif traj.has_transmitted:
                        transmitted += 1
                        print("Particle was transmitted in traj {}".format(traj.traj_id))
                        logfile.write("Particle was transmitted in traj {}\n".format(traj.traj_id))
                else:
                	absorbed += 1

        #print("total traj: {}".format(ntrajs))
        #logfile.write("total traj: {}\n".format(ntrajs))
        #print("# traj with scattering after barrier: {} out of {} scattered traj ({}%)".format(sc_count,scattered,float(sc_count)*100/float(scattered)))
        #logfile.write("# traj with scattering after barrier: {} out of {} scattered traj ({}%)\n".format(sc_count,scattered,float(sc_count)*100/float(scattered)))
        print("trajs with scattering: {} out of total {} traj ({:4.2f}%)".format(scattered,ntrajs,float(scattered)*float(100)/float(ntrajs)))
        print("trajs with adsorption: {} out of total {} traj ({:4.2f}%)".format(absorbed,ntrajs,float(absorbed)*float(100)/float(ntrajs)))
        print("trajs transmitted: {} out of total {} traj ({:4.2f}%)".format(transmitted,ntrajs,float(transmitted)*float(100)/float(ntrajs)))
        logfile.write("trajs with scattering: {} out of total {} traj ({:4.2f}%)\n".format(scattered,ntrajs,float(scattered)*float(100)/float(ntrajs)))
        logfile.write("trajs with adsorption: {} out of total {} traj ({:4.2f}%)\n".format(absorbed,ntrajs,float(absorbed)*float(100)/float(ntrajs)))
        logfile.write("trajs transmitted: {} out of total {} traj ({:4.2f}%)\n".format(transmitted,ntrajs,float(transmitted)*float(100)/float(ntrajs)))

        return traj_list, scattered, absorbed, transmitted

def numbins(inp):
	if isinstance(inp, list):
		return int(3*(len(inp)**(1./3)))
	elif isinstance(inp, int):
		return int(3*((inp)**(1./3)))
	elif isinstance(inp, float):
		return int(3*((inp)**(1./3)))
	else:
		sys.exit("Unknown type from which to compute number of bins in histogram")

def analyze(trajs,logfile):

        ### BOUNCES ###
        print("Calculating bounces.")
        logfile.write("Calculating bounces.\n")

        # LOOP 
        all_bounces    = [traj.turn_pnts for traj in trajs]
        scat_bounces   = [traj.turn_pnts for traj in trajs if traj.has_scattered]
        abso_bounces   = [traj.turn_pnts for traj in trajs if traj.has_adsorbed]
        transm_bounces = [traj.turn_pnts for traj in trajs if traj.has_transmitted]
        
        # ANALYSIS
        all_bounce_hist,    all_bounce_edges    = numpy.histogram(all_bounces,    bins=max(all_bounces), range=(0,max(all_bounces)), density=True)
        scat_bounce_hist,   scat_bounce_edges   = numpy.histogram(scat_bounces,   bins=max(all_bounces), range=(0,max(all_bounces)), density=True)
        abso_bounce_hist,   abso_bounce_edges   = numpy.histogram(abso_bounces,   bins=max(all_bounces), range=(0,max(all_bounces)), density=True)
        transm_bounce_hist, transm_bounce_edges = numpy.histogram(transm_bounces, bins=max(all_bounces), range=(0,max(all_bounces)), density=True)

        # OUTPUT
        bounce_file = open("analysis/bounces.txt", "w")
        bounce_file.write("# bounces  all  scattered  absorbed  transmitted\n")
        for i in range(len(all_bounce_hist)):
        	bounce_file.write("%d %f %f %f %f\n" % ( 0.5*(all_bounce_edges[i]+all_bounce_edges[i+1]), all_bounce_hist[i], FRAC_SCATTERED*scat_bounce_hist[i], FRAC_ABSORBED*abso_bounce_hist[i], FRAC_TRANSMITTED*transm_bounce_hist[i]))
        bounce_file.close()


        ### TOTAL ENERGY LOSS ###
        print("Calculating total energy loss.")
        logfile.write("Calculating total energy loss.\n")
        # LOOP 
        all_eloss = [traj.eloss for traj in trajs if traj.has_scattered] 
        #all_efrac = [traj.efrac for traj in trajs if traj.has_scattered] 
        one_b     = [traj.eloss for traj in trajs if traj.has_scattered and traj.turn_pnts == 1]
        two_b     = [traj.eloss for traj in trajs if traj.has_scattered and traj.turn_pnts == 3]
        mul_b     = [traj.eloss for traj in trajs if traj.has_scattered and traj.turn_pnts >= 5]
        
        absorbed_eloss = [traj.eloss for traj in trajs if traj.has_adsorbed]

        # ANALYSIS
        all_eloss_hist, all_eloss_edges = numpy.histogram(all_eloss, bins=numbins(SCATTERED), range=(min(all_eloss), max(all_eloss)), density=True)
        #all_efrac_hist, all_efrac_edges = numpy.histogram(all_efrac, bins=numbins(SCATTERED), range=(min(all_efrac), max(all_efrac)), density=True)
        one_b_hist,     one_b_edges     = numpy.histogram(one_b,     bins=numbins(SCATTERED), range=(min(all_eloss), max(all_eloss)), density=True)
        two_b_hist,     two_b_edges     = numpy.histogram(two_b,     bins=numbins(SCATTERED), range=(min(all_eloss), max(all_eloss)), density=True)
        mul_b_hist,     mul_b_edges     = numpy.histogram(mul_b,     bins=numbins(SCATTERED), range=(min(all_eloss), max(all_eloss)), density=True)
        frac_one_b = float(len(one_b))/SCATTERED
        frac_two_b = float(len(two_b))/SCATTERED
        frac_mul_b = float(len(mul_b))/SCATTERED

        # OUTPUT 
        eloss_file = open("analysis/eloss.txt", "w")
        eloss_file.write("# eloss/eV  all  single bounce  double bounce  multi bounce\n")
        for i in range(len(all_eloss_hist)):
        	eloss_file.write("%f %f %f %f %f\n" % (0.5*(all_eloss_edges[i]+all_eloss_edges[i+1]), all_eloss_hist[i], frac_one_b*one_b_hist[i], frac_two_b*two_b_hist[i], frac_mul_b*mul_b_hist[i]))
        eloss_file.close()


        ### SPECULAR ENERGY LOSS ###
        print("Calculating specular energy loss.")
        logfile.write("Calculating specular energy loss.\n")

        # INIT
        spec_all_eloss = []
        #spec_all_efrac = []
        spec_one_b     = []
        spec_two_b     = []
        spec_mul_b     = []

        # LOOP
        for traj in trajs:	# List comprehension simply need too much time. This is ugly, but fast.
                if traj.in_spec and traj.has_scattered:
                        spec_all_eloss.append(traj.eloss)
                        #spec_all_efrac.append(traj.efrac)
                        if traj.turn_pnts == 1:
                        	spec_one_b.append(traj.eloss)
                        elif traj.turn_pnts == 3:
                        	spec_two_b.append(traj.eloss)
                        else:
                        	spec_mul_b.append(traj.eloss)

        spec_eloss_file = open("analysis/spec_eloss.txt", "w")
        spec_eloss_file.write("# eloss/eV  all  single bounce  double bounce  multi bounce\n")
        if len(spec_all_eloss) > 0:	
                # ANALYSIS 
                spec_all_eloss_hist, spec_all_eloss_edges = numpy.histogram(spec_all_eloss, bins=numbins(spec_all_eloss), range=(min(spec_all_eloss), max(spec_all_eloss)), density=True)
                #spec_all_efrac_hist, spec_all_efrac_edges = numpy.histogram(spec_all_efrac, bins=numbins(spec_all_efrac), range=(min(spec_all_efrac), max(spec_all_efrac)), density=True)
                spec_one_b_hist,     spec_one_b_edges     = numpy.histogram(spec_one_b,     bins=numbins(spec_all_eloss), range=(min(spec_all_eloss), max(spec_all_eloss)), density=True)
                spec_two_b_hist,     spec_two_b_edges     = numpy.histogram(spec_two_b,     bins=numbins(spec_all_eloss), range=(min(spec_all_eloss), max(spec_all_eloss)), density=True)
                spec_mul_b_hist,     spec_mul_b_edges     = numpy.histogram(spec_mul_b,     bins=numbins(spec_all_eloss), range=(min(spec_all_eloss), max(spec_all_eloss)), density=True)
                spec_frac_one_b = float(len(spec_one_b))/len(spec_all_eloss)
                spec_frac_two_b = float(len(spec_two_b))/len(spec_all_eloss)
                spec_frac_mul_b = float(len(spec_mul_b))/len(spec_all_eloss)

                # OUTPUT 
                spec_eloss_file = open("analysis/spec_eloss.txt", "w")
                spec_eloss_file.write("# eloss/eV  all  single bounce  double bounce  multi bounce\n")
                for i in range(len(spec_all_eloss_hist)):
	                spec_eloss_file.write("%f %f %f %f %f\n" % (0.5*(spec_all_eloss_edges[i]+spec_all_eloss_edges[i+1]), spec_all_eloss_hist[i], spec_frac_one_b*spec_one_b_hist[i], spec_frac_two_b*spec_two_b_hist[i], spec_frac_mul_b*spec_mul_b_hist[i]))
        else:
        	spec_eloss_file.write("%f %f %f %f %f\n" % ( 0.0, 0.0, 0.0, 0.0, 0.0))
        spec_eloss_file.close()


        ### IN PLANE ENERGY LOSS ###
        print("Calculating in-plane energy loss.")
        logfile.write("Calculating in-plane energy loss.\n")

        pp_file = open("analysis/plot_parameter.txt", "w") # here the parameters needed for plotting are stored
        mdtinpfile = open(mdtinpname, "r")
        for line in mdtinpfile:
            if not line.startswith("!"): # skip comment lines
                if "Tsurf" in line:
                    temp = float(line.split()[-1]) # "Tsurf 300"
        pp_file.write("Einc Vinc Ainc Temp Detector_radius Bins\n{} {} {} {:f} {} {}".format(trajs[0].ekin_p_i,length(trajs[0].v_p_i),trajs[0].polar_i,temp,SPECULAR_RADIUS,BINS))
        pp_file.close()


        # INIT
        in_plane_all_eloss = []
        in_plane_one_b     = []
        in_plane_two_b     = []
        in_plane_mul_b     = []

        # LOOP
        for traj in trajs:	# List comprehension simply need too much time. This is ugly, but fast.
                if traj.in_plane and traj.has_scattered:

                        in_plane_all_eloss.append(traj.eloss)

                        if traj.turn_pnts == 1:
                        	in_plane_one_b.append(traj.eloss)
                        elif traj.turn_pnts == 3:
                        	in_plane_two_b.append(traj.eloss)
                        else:
                        	in_plane_mul_b.append(traj.eloss)


        in_plane_eloss_file = open("analysis/in_plane_eloss.txt", "w")
        in_plane_eloss_file.write("# eloss/eV  all  single bounce  double bounce  multi bounce\n")
        if len(in_plane_all_eloss) > 0:	
                # ANALYSIS 
                in_plane_all_eloss_hist, in_plane_all_eloss_edges = numpy.histogram(in_plane_all_eloss, bins=numbins(in_plane_all_eloss), range=(min(in_plane_all_eloss), max(in_plane_all_eloss)), density=True)
                in_plane_one_b_hist,     in_plane_one_b_edges     = numpy.histogram(in_plane_one_b,     bins=numbins(in_plane_all_eloss), range=(min(in_plane_all_eloss), max(in_plane_all_eloss)), density=True)
                in_plane_two_b_hist,     in_plane_two_b_edges     = numpy.histogram(in_plane_two_b,     bins=numbins(in_plane_all_eloss), range=(min(in_plane_all_eloss), max(in_plane_all_eloss)), density=True)
                in_plane_mul_b_hist,     in_plane_mul_b_edges     = numpy.histogram(in_plane_mul_b,     bins=numbins(in_plane_all_eloss), range=(min(in_plane_all_eloss), max(in_plane_all_eloss)), density=True)
                in_plane_frac_one_b = float(len(in_plane_one_b))/len(in_plane_all_eloss)
                in_plane_frac_two_b = float(len(in_plane_two_b))/len(in_plane_all_eloss)
                in_plane_frac_mul_b = float(len(in_plane_mul_b))/len(in_plane_all_eloss)

                # OUTPUT 
                in_plane_eloss_file = open("analysis/in_plane_eloss.txt", "w")
                in_plane_eloss_file.write("# eloss/eV  all  single bounce  double bounce  multi bounce\n")
                for i in range(len(in_plane_all_eloss_hist)):
                        in_plane_eloss_file.write("%f %f %f %f %f\n" % (0.5*(in_plane_all_eloss_edges[i]+in_plane_all_eloss_edges[i+1]), in_plane_all_eloss_hist[i], in_plane_frac_one_b*in_plane_one_b_hist[i], in_plane_frac_two_b*in_plane_two_b_hist[i], in_plane_frac_mul_b*in_plane_mul_b_hist[i]))
        else:
        	in_plane_eloss_file.write("%f %f %f %f %f\n" % ( 0.0, 0.0, 0.0, 0.0, 0.0))
        in_plane_eloss_file.close()



        ### Z-POSITION ###
        print("Calculating final z positions.")
        logfile.write("Calculating final z positions.\n")
        # LOOP
        final_z = [traj.r_p_f.z for traj in trajs if traj.has_adsorbed]

        if (len(final_z) > 0):
        	# ANALYSIS
        	final_z_hist, final_z_edges = numpy.histogram(final_z, bins=numbins(final_z), range=(min(final_z), max(final_z)), density=True)

        	# OUTPUT
        	final_z_file = open("analysis/final_z.txt", "w")
        	final_z_file.write("# z/A  probability density\n")
        	for i in range(len(final_z_hist)):
        		final_z_file.write("%f %f\n" % (0.5*(final_z_edges[i]+final_z_edges[i+1]), final_z_hist[i]))
        	final_z_file.close()
        	

        ### BOUNCES VS ELOSS ###
        print("Calculating bounces/energy loss correlation.")
        logfile.write("Calculating bounces/energy loss correlation.\n")
        # ANALYSIS
        bounce_vs_eloss_hist, xedges, yedges = numpy.histogram2d(scat_bounces, all_eloss, bins=(max(scat_bounces),numbins(all_eloss)), range=[[0,max(scat_bounces)],[min(all_eloss),max(all_eloss)]], normed=False)

        # OUTPUT
        out = open("analysis/bounces_vs_eloss.txt", "w")
        out.write("# bounces  eloss/eV  counts\n")
        for i in range(len(bounce_vs_eloss_hist)):
        	for j in range(len(bounce_vs_eloss_hist[i])):
        		out.write("%d %f %d\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), bounce_vs_eloss_hist[i][j]))
        	out.write("\n")
        out.close()


        ### ANGULAR DISTRIBUTION ###
        print("Calculating angular energy loss.")
        logfile.write("Calculating angular energy loss.\n")
        # get trajectories that are within specular radius in azimuth direction
        energy_collect  = []
        efrac_collect   = []
        angle_collect   = []
        ps_dist_collect = []
        polar_scatt_azi_int_energy = []
        polar_scatt_azi_int_efrac  = []
        polar_scatt_azi_int_angle  = []
        #total_counts = 0
        for traj in trajs:
        	if traj.has_scattered:
        		polar_scatt_azi_int_energy.append(traj.eloss)
        		polar_scatt_azi_int_efrac.append(traj.efrac)
        		polar_scatt_azi_int_angle.append(traj.polar_f)

        		if traj.in_plane:
        			energy_collect.append(traj.eloss)
        			efrac_collect.append(traj.efrac)
        			angle_collect.append(traj.polar_f)
        			ps_dist_collect.append(traj.cl_appr)

        ang_dist_file          = open("analysis/ang_res_eloss.txt", "w")
        ang_dist_file_norm     = open("analysis/ang_res_eloss_norm.txt", "w")
        ang_dist_mat_file      = open("analysis/ang_res_eloss_matrix.txt", "w")
        ang_dist_mat_file_norm = open("analysis/ang_res_eloss_matrix_norm.txt", "w")
        occurence_file         = open("analysis/ang_res_occurrence.txt", "w")
        if len(energy_collect) != 0:
                angle_eloss_hist, xedges, yedges = numpy.histogram2d(energy_collect, angle_collect, bins=(numbins(energy_collect)), density=False)
                angle_efrac_hist, xefedges, yefedges = numpy.histogram2d(efrac_collect, angle_collect, bins=BINS, range=[[0, 1.1],[0, 90]], density=False) # bins=90
                #occurence_hist, occ_edges = numpy.histogram(angle_collect, bins=numpy.arange(91), normed=False)
                occurence_hist, occ_edges = numpy.histogram(angle_collect, bins=numpy.arange(91), density=False)

                # OUTPUT
                for i in range(len(angle_eloss_hist)):
                    for j in range(len(angle_eloss_hist[i])):
                        ang_dist_file.write("%f %f %d\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), int(angle_eloss_hist[i][j])/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))))
                	#ang_dist_file.write("%f %f %d\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), angle_eloss_hist[i][j]))
                        ang_dist_file.write("\n")

                ang_dist_mat_file.write("# x-range describing energy loss in eV (left to right) from %f to %f in steps of %f\n" % (0.5*(xedges[0]+xedges[1]), 0.5*(xedges[-2]+xedges[-1]), abs(xedges[0]-xedges[1])))
                ang_dist_mat_file.write("# y-range describing scattering angle in degrees (top to bottom) from %f to %f in steps of %f\n" % (0.5*(yedges[0]+yedges[1]), 0.5*(yedges[-2]+yedges[-1]), abs(yedges[0]-yedges[1])))
                ang_dist_mat_file.write("# specular scattering angle is %f degrees and detector radius is %f degrees\n" % (trajs[0].polar_i, SPECULAR_RADIUS))
                ang_dist_mat_file.write("{:8.4f}".format(0)); [ang_dist_mat_file.write("{:8.4f}".format((0.5*(xedges[i]+xedges[i+1])))) for i in range(len(xedges)-1)]; ang_dist_mat_file.write("\n")
                ang_dist_mat_file.write("{:8.4f}".format(0)); [ang_dist_mat_file.write("{:8.4f}".format((trajs[0].ekin_p_i - (0.5*(xedges[i]+xedges[i+1]))))) for i in range(len(xedges)-1)]; ang_dist_mat_file.write("\n")
                for i in range(len(angle_eloss_hist[0])):
                        ang_dist_mat_file.write("{:8.4f}".format(0.5*(yedges[i]+yedges[i+1])))
                        for j in range(len(angle_eloss_hist)):
                            ang_dist_mat_file.write("{:4d}".format(int(angle_eloss_hist[j][i])))
                            #total_counts += int(angle_eloss_hist[j][i])
                        ang_dist_mat_file.write("\n")
                
                # OUTPUT without and with norm
                for i in range(len(angle_efrac_hist)):
                    for j in range(len(angle_efrac_hist[i])):
                        abs_val = int(angle_efrac_hist[i][j]/abs(numpy.sin(0.5*(yefedges[j]+yefedges[j+1])/360*2*numpy.pi))) # correct for angle in experiment
                        #abs_val = int(angle_eloss_hist[i][j]/numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))
                        rel_val = float(abs_val)/float(angle_efrac_hist.sum()) # get flux (integrated over all angles)
                        ang_dist_file_norm.write("%f %f %f\n" % (0.5*(xefedges[i]+xefedges[i+1]), 0.5*(yefedges[j]+yefedges[j+1]), rel_val))

                #print("Total counts in-plane {}".format(total_counts))
                #logfile.write("Total counts in-plane {}\n".format(total_counts))

                ang_dist_mat_file_norm.write("# x-range describing energy in eV (left to right) from %f to %f in steps of %f\n" % (0.5*(xedges[0]+xedges[1]), 0.5*(xedges[-2]+xedges[-1]), abs(xedges[0]-xedges[1])))
                ang_dist_mat_file_norm.write("# y-range describing scattering angle in degrees (top to bottom) from %f to %f in steps of %f\n" % (0.5*(yedges[0]+yedges[1]), 0.5*(yedges[-2]+yedges[-1]), abs(yedges[0]-yedges[1])))
                ang_dist_mat_file_norm.write("# specular scattering angle is {} degrees, incidence kinetic energy is {} eV and temperature is {} K\n".format(trajs[0].polar_i, trajs[0].ekin_p_i, temp))
                ang_dist_mat_file_norm.write("# detector radius is %f degrees and total number of counts is %d \n" % (SPECULAR_RADIUS, angle_eloss_hist.sum()))
                #ang_dist_mat_file_norm.write("{:8.4f}".format(0)); [ang_dist_mat_file_norm.write("{:8.4f}".format((0.5*(xedges[i]+xedges[i+1])))) for i in range(len(xedges)-1)]; ang_dist_mat_file_norm.write("\n")
                ang_dist_mat_file_norm.write("{:8.4f}".format(0)); [ang_dist_mat_file_norm.write("{:8.4f}".format((trajs[0].ekin_p_i - (0.5*(xedges[i]+xedges[i+1]))))) for i in range(len(xedges)-1)]; ang_dist_mat_file_norm.write("\n")
                for i in range(len(angle_eloss_hist[0])):
                         ang_dist_mat_file_norm.write("{:8.4f}".format(0.5*(yedges[i]+yedges[i+1])))
                         for j in range(len(angle_eloss_hist)):
                             ang_dist_mat_file_norm.write("{:8.4f}".format(float(angle_eloss_hist[j][i]/angle_eloss_hist.sum())))
                         ang_dist_mat_file_norm.write("\n")
                ang_dist_mat_file.write("Total count: {}".format(angle_eloss_hist.sum()))
                		

                occurence_file.write("# number of trajs in plane in unit polar angle, unit azimuthal angle\n")
                for i in range(len(occurence_hist)):
                	this_angle = 0.5*(occ_edges[i]+occ_edges[i+1])
                	occurence_file.write("%f %f\n" % (this_angle, 1.0*occurence_hist[i]/sum(occurence_hist)))
        else:
                ang_dist_file.write("%f %f %d\n"   % (0.1, 0.5, 0))
                ang_dist_file.write("%f %f %d\n\n" % (0.1, 1.0, 1))
                ang_dist_file.write("%f %f %d\n"   % (0.2, 1.0, 1))
                ang_dist_file.write("%f %f %d\n"   % (0.2, 0.5, 2))
                
                #ang_dist_file_norm.write("%f %f %d\n"   % (0.1, 0.5, 0)) # don't know if this is needed and what the numbers should be
                #ang_dist_file_norm.write("%f %f %d\n\n" % (0.1, 1.0, 1))
                #ang_dist_file_morm.write("%f %f %d\n"   % (0.2, 1.0, 1))
                #ang_dist_file_norm.write("%f %f %d\n"   % (0.2, 0.5, 2))

                occurence_file.write("%f %f\n" % (0.1, 0.1))

        ang_dist_mat_file.write("\n")
        occurence_file.close()
        ang_dist_file.close()
        ang_dist_mat_file.close()
        ang_dist_mat_file_norm.close()


        # Graphene bounce events #
        #print("Calculating graphene bounce events.")
        #logfile.write("Calculating graphene bounce events.\n")
        #fast_c = open("analysis/component_fast.txt", "w")
        #slow_c_sb = open("analysis/component_slow_single.txt", "w")
        #slow_c_mb = open("analysis/component_slow_multi.txt", "w")
        #for traj in trajs:
        	#if traj.has_scattered and traj.in_plane:
        		#if traj.turn_pnts == 1:
        			#if traj.cl_appr > 1.4:
        				#fast_c.write("%f %f\n" % (traj.ekin_p_f/traj.ekin_p_i, traj.polar_f))
        			#if traj.cl_appr < 1.4:
        				#slow_c_sb.write("%f %f\n" % (traj.ekin_p_f/traj.ekin_p_i, traj.polar_f))
        		#elif traj.turn_pnts > 1 and traj.cl_appr < 1.4:
        				#slow_c_mb.write("%f %f\n" % (traj.ekin_p_f/traj.ekin_p_i, traj.polar_f))
        #slow_c_mb.close(); slow_c_sb.close(); fast_c.close()
        			
        			
        			



        # INTERGRATED OVER ALL AZIMUTH ANGLES #
        polar_scatt_azi_file = open("analysis/polar_scatt_azi_int.txt", "w")
        #total_counts_all = 0
        if len(polar_scatt_azi_int_energy) != 0:
        	polar_scatt_azi_int_hist, xedges, yedges = numpy.histogram2d(polar_scatt_azi_int_energy, polar_scatt_azi_int_angle, bins=(numbins(polar_scatt_azi_int_energy)), density=False)
        	
        	for i in range(len(polar_scatt_azi_int_hist)):
        		for j in range(len(polar_scatt_azi_int_hist[i])):
                                polar_scatt_azi_file.write("%f %f %d\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), int(polar_scatt_azi_int_hist[i][j])/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))))
                                #total_counts_all += int(polar_scatt_azi_int_hist[i][j])

        		polar_scatt_azi_file.write("\n")
        else:
        	polar_scatt_azi_file.write("%f %f %d\n"   % (0.1, 0.5, 0))
        	polar_scatt_azi_file.write("%f %f %d\n\n" % (0.1, 1.0, 1))
        	polar_scatt_azi_file.write("%f %f %d\n"   % (0.2, 1.0, 1))
        	polar_scatt_azi_file.write("%f %f %d\n"   % (0.2, 0.5, 2))

        polar_scatt_azi_file.close()
        #print("Total counts all angles: {}".format(int(total_counts_all)))
        #logfile.write("Total counts all angles: {}\n".format(int(total_counts_all)))

        polar_scatt_azi_file_norm = open("analysis/polar_scatt_azi_int_norm.txt", "w")
        if len(polar_scatt_azi_int_efrac) != 0:
                #polar_scatt_azi_int_hist, xedges, yedges = numpy.histogram2d(polar_scatt_azi_int_energy, polar_scatt_azi_int_angle, bins=(numbins(polar_scatt_azi_int_energy)), normed=False)
                polar_scatt_azi_int_hist_f, xefedges, yefedges = numpy.histogram2d(polar_scatt_azi_int_efrac, polar_scatt_azi_int_angle, bins=BINS, range=[[0, 1.1],[0, 90]], density=False)

                for i in range(len(polar_scatt_azi_int_hist_f)):
                        for j in range(len(polar_scatt_azi_int_hist_f[i])):
                                abs_val = int(polar_scatt_azi_int_hist_f[i][j]/abs(numpy.sin(0.5*(yefedges[j]+yefedges[j+1])/360*2*numpy.pi)))
                                rel_val = float(abs_val)/float(polar_scatt_azi_int_hist_f.sum())
                                polar_scatt_azi_file_norm.write("%f %f %f\n" % (0.5*(xefedges[i]+xefedges[i+1]), 0.5*(yefedges[j]+yefedges[j+1]), rel_val))

                        #polar_scatt_azi_file_norm.write("\n")
        else:
                polar_scatt_azi_file_norm.write("%f %f %d\n"   % (0.1, 0.5, 0))
                polar_scatt_azi_file_norm.write("%f %f %d\n\n" % (0.1, 1.0, 1))
                polar_scatt_azi_file_norm.write("%f %f %d\n"   % (0.2, 1.0, 1))
                polar_scatt_azi_file_norm.write("%f %f %d\n"   % (0.2, 0.5, 2))



        ### LOSS TO EHP AND PHONONS ###
        print("Calculating loss to ehps and phonons.")
        logfile.write("Calculating loss to ehps and phonons.\n")
        # LOOP
        loss_to_ehps = [ traj.etotal_i - traj.etotal_f for traj in trajs if traj.has_scattered ]
        loss_to_ehps_spec = [ traj.etotal_i - traj.etotal_f for traj in trajs if traj.has_scattered and traj.in_spec ]
        
        # ANALYSIS
        loss_to_ehps_hist, loss_to_ehps_edges = numpy.histogram(loss_to_ehps, bins=numbins(loss_to_ehps), range=(min(loss_to_ehps), max(loss_to_ehps)), density=True)
        if len(loss_to_ehps_spec) > 0:
        	loss_to_ehps_spec_hist, loss_to_ehps_spec_edges = numpy.histogram(loss_to_ehps_spec, bins=numbins(loss_to_ehps_spec), range=(min(loss_to_ehps_spec), max(loss_to_ehps_spec)), density=True)

        # OUTPUT
        loss_to_ehps_file = open("analysis/eloss_to_ehps.txt", "w")
        for i in range(len(loss_to_ehps_hist)):
        	loss_to_ehps_file.write("%f %f\n" % (0.5*(loss_to_ehps_edges[i]+loss_to_ehps_edges[i+1]), loss_to_ehps_hist[i]))
        loss_to_ehps_file.close()

        
        loss_to_ehps_spec_file = open("analysis/eloss_to_ehps_spec.txt", "w")
        if len(loss_to_ehps_spec) > 0:
                for i in range(len(loss_to_ehps_spec_hist)):
                        loss_to_ehps_spec_file.write("%f %f\n" % (0.5*(loss_to_ehps_spec_edges[i]+loss_to_ehps_spec_edges[i+1]), loss_to_ehps_spec_hist[i]))
        else:
        	loss_to_ehps_spec_file.write("%f %f\n" % (0.0, 0.0))
        loss_to_ehps_spec_file.close()


        ### SPHERICAL SYMMETRY ###
        # LOOP
        print("Calculating spherical symmetry.")
        logfile.write("Calculating spherical symmetry.\n")
        abs_azi = []
        rel_azi = []
        yvals = []
        this_azi = None
        for traj in trajs:
        	if traj.has_scattered:
        		delta_azi = traj.azi_f-traj.azi_i
        		if -180 <= delta_azi <= 180:
        			rel_azi.append(delta_azi)
        		elif delta_azi < -180:
        			rel_azi.append(delta_azi+360)
        		elif delta_azi > 180:
        			rel_azi.append(delta_azi-360)
        		else:
                                print("Weird angle in spherical symmetry.")
                                logfile.write("Weird angle in spherical symmetry.\n")
                                sys.exit()

        		abs_azi.append(traj.azi_f)
        		yvals.append(traj.polar_f)

        # ANALYSIS
        rel_spherical_hist, rel_xedges, yedges = numpy.histogram2d(rel_azi, yvals, bins=numbins(rel_azi), normed=False)
        abs_spherical_hist, abs_xedges, yedges = numpy.histogram2d(abs_azi, yvals, bins=numbins(abs_azi), normed=False)

        # OUTPUT
        spherical_file = open("analysis/rel_spherical_symmetry.txt", "w")
        for i in range(len(rel_spherical_hist)):
        	for j in range(len(rel_spherical_hist[i])):
        		spherical_file.write("%f %f %d\n" % (0.5*(rel_xedges[i]+rel_xedges[i+1]), -0.5*(yedges[j]+yedges[j+1]), int(rel_spherical_hist[i][j])/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))))
        	spherical_file.write("\n")
        spherical_file.close()

        spherical_file = open("analysis/abs_spherical_symmetry.txt", "w")
        for i in range(len(abs_spherical_hist)):
        	for j in range(len(abs_spherical_hist[i])):
        		spherical_file.write("%f %f %d\n" % (0.5*(abs_xedges[i]+abs_xedges[i+1]), -0.5*(yedges[j]+yedges[j+1]), int(abs_spherical_hist[i][j])/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))))
        	spherical_file.write("\n")
        spherical_file.close()

        ### 1D ANGULAR DISTRIBUTION
        print("Calculate 1D angular distribution")
        angplane_collect = []
        for traj in trajs:
            if traj.azi_f < 0:
                inv_azi = AZIMUTHAL_ANGLE - 180
                if traj.has_scattered and abs(inv_azi - traj.azi_f) < SPECULAR_RADIUS:
                    angplane_collect.append(-traj.polar_f)
                elif traj.has_scattered and abs(AZIMUTHAL_ANGLE - traj.azi_f) < SPECULAR_RADIUS:
                    angplane_collect.append(traj.polar_f)
            if traj.azi_f > 0:
                inv_azi = AZIMUTHAL_ANGLE + 180
                if traj.has_scattered and abs(inv_azi - traj.azi_f) < SPECULAR_RADIUS:
                    angplane_collect.append(-traj.polar_f)
                elif traj.has_scattered and abs(AZIMUTHAL_ANGLE - traj.azi_f) < SPECULAR_RADIUS:
                    angplane_collect.append(traj.polar_f)

        ang_dist_file = open("analysis/ang_dist.txt", "w")
        if len(angplane_collect) > 0:
            angle_hist, angplane_collect_edges = numpy.histogram(angplane_collect, bins=BINS, range=(-90, 90), density=False)
            for i in range(len(angle_hist)):
                angle_hist[i] = angle_hist[i]/abs(numpy.sin(0.5*(angplane_collect_edges[i]+angplane_collect_edges[i+1])/360*2*numpy.pi)) # needs to be divided by sin(x) to correct geometry of experiment.

            #OUTPUT
            for i in range(len(angle_hist)):
                ang_dist_file.write("%f %f %f\n" % (0.5*(angplane_collect_edges[i]+angplane_collect_edges[i+1]), angle_hist[i], angle_hist[i]/float(max(angle_hist))))
            ang_dist_file.close()


        ### 2D ANGULAR DISTRIBUTION
        print("Calculate 2D angular distribution")
        angle_all_collect = [] # all polar angle from -90 to 90 in deg
        #angle_inp_collect = [] # in-plane polar angle from -90 to 90 in deg
        efrac_all_collect = [] # all E_s / E_i in eV
        #efrac_inp_collect = [] # in-plane E_s / E_i in eV

        for traj in trajs:
            if traj.azi_f < 0:
                inv_azi = AZIMUTHAL_ANGLE - 180
                if traj.has_scattered and abs(inv_azi - traj.azi_f) < SPECULAR_RADIUS:
                    angle_all_collect.append(-traj.polar_f)
                    efrac_all_collect.append(traj.efrac)
                elif traj.has_scattered and abs(AZIMUTHAL_ANGLE - traj.azi_f) < SPECULAR_RADIUS:
                    angle_all_collect.append(traj.polar_f)
                    efrac_all_collect.append(traj.efrac)
            if traj.azi_f > 0:
                inv_azi = AZIMUTHAL_ANGLE + 180
                if traj.has_scattered and abs(inv_azi - traj.azi_f) < SPECULAR_RADIUS:
                    angle_all_collect.append(-traj.polar_f)
                    efrac_all_collect.append(traj.efrac)
                elif traj.has_scattered and abs(AZIMUTHAL_ANGLE - traj.azi_f) < SPECULAR_RADIUS:
                    angle_all_collect.append(traj.polar_f)
                    efrac_all_collect.append(traj.efrac)

        ang_dist_nrg_ang      = open("analysis/2d-ang-dist.txt", "w")
        ang_dist_nrg_ang_norm = open("analysis/2d-ang-dist_norm.txt", "w")

        if len(efrac_all_collect) != 0 and len(angle_all_collect) != 0:
            #angle_efrac_hist, xedges, yedges = numpy.histogram2d(efrac_all_collect, angle_all_collect,  bins=(numbins(efrac_all_collect)), density=False)
            #angle_efrac_hist, xedges, yedges = numpy.histogram2d(efrac_all_collect, angle_all_collect, bins=BINS, range=[[0, 1.1],[-90, 90]], density=False) # fixed bin size of 2 deg; try 36,72,180
            angle_efrac_hist, xedges, yedges = numpy.histogram2d(efrac_all_collect, angle_all_collect, bins=BINS, range=[[0, 1.1],[ANGLE_MIN, ANGLE_MAX]], density=False) # fixed bin size of 2 deg; try 36,72,180

            for i in range(len(angle_efrac_hist)):
                for j in range(len(angle_efrac_hist[i])):
                    abs_val = int(angle_efrac_hist[i][j]/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))) # needs to be divided by sin(x) to correct geometry of experiment.
                    rel_val = float(abs_val)/float(angle_efrac_hist.sum()) # current bin devided by sum of all bins to get "flux"
                    ang_dist_nrg_ang.write("%f %f %d\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), abs_val)) # write data with bins
                    ang_dist_nrg_ang_norm.write("%f %f %f\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), rel_val)) # write data with "flux""

        ang_dist_nrg_ang.close()
        ang_dist_nrg_ang_norm.close()


        ### Projectile-Surface distance ###
        print("Calculating projectile-surface distance.")
        logfile.write("Calculating projectile-surface distance.\n")
        # LOOP
        ps_dist = [traj.cl_appr for traj in trajs if traj.has_scattered]

        # ANALYSIS
        ps_hist, xedges = numpy.histogram(ps_dist, bins=numbins(ps_dist), range=(min(ps_dist), max(ps_dist)), density=True)
        
        # OUTPUT
        ps_file = open("analysis/ps_dist.txt", "w")
        for i in range(len(ps_hist)):
                ps_file.write("%f %f\n" % (0.5*(xedges[i]+xedges[i+1]), ps_hist[i]))
        ps_file.close()



        ### Eloss vs Projectile-Surface distance ###
        print("Calculating energy loss projectile-surface distance relationship.")
        logfile.write("Calculating energy loss projectile-surface distance relationship.\n")
        # ANALYSIS
        eloss_psd_hist, xedges, yedges = numpy.histogram2d(all_eloss, ps_dist, bins=numbins(all_eloss), normed=False)
        
        # OUTPUT
        eloss_psd_file = open("analysis/eloss_psd.txt", "w")
        for i in range(len(eloss_psd_hist)):
                for j in range(len(eloss_psd_hist[i])):
                        eloss_psd_file.write("%f %f %d\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), eloss_psd_hist[i][j]))
                eloss_psd_file.write("\n")
        eloss_psd_file.close()



        ### Scattering polar angle vs Projectile-Surface distance in-plane ###
        print("Calculating scattering angle projectile-surface distance relationship.")
        logfile.write("Calculating scattering angle projectile-surface distance relationship.\n")
        # ANALYSIS
        polar_psd_file = open("analysis/polar_psd.txt", "w")
        if len(ps_dist_collect) > 0:
        	polar_psd_hist, xedges, yedges = numpy.histogram2d(ps_dist_collect, angle_collect, bins=numbins(ps_dist_collect), normed=False)

        	# OUTPUT
        	for i in range(len(polar_psd_hist)):
                        for j in range(len(polar_psd_hist[i])):
                                polar_psd_file.write("%f %f %d\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), int(polar_psd_hist[i][j])/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))))
                        polar_psd_file.write("\n")
        else:
                polar_psd_file.write("%f %f %d\n"   % (0.1, 0.5, 0))
                polar_psd_file.write("%f %f %d\n\n" % (0.1, 1.0, 1))
                polar_psd_file.write("%f %f %d\n"   % (0.2, 1.0, 1))
                polar_psd_file.write("%f %f %d\n"   % (0.2, 0.5, 2))
        polar_psd_file.close()


        ### Eloss vs Projectile-Surface distance in plane ###
        print("Calculating in-plane energy loss projectile-surface distance relationship.")
        logfile.write("Calculating in-plane energy loss projectile-surface distance relationship.\n")
        # ANALYSIS
        eloss_psd_in_plane_file = open("analysis/eloss_psd_in_plane.txt", "w")
        if len(ps_dist_collect) > 0:
        	eloss_psd_in_plane_hist, xedges, yedges = numpy.histogram2d(energy_collect, ps_dist_collect, bins=numbins(ps_dist_collect), normed=False)

        	# OUTPUT
        	for i in range(len(eloss_psd_in_plane_hist)):
                        for j in range(len(eloss_psd_in_plane_hist[i])):
                                eloss_psd_in_plane_file.write("%f %f %d\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), int(eloss_psd_in_plane_hist[i][j])/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))))
                        eloss_psd_in_plane_file.write("\n")
        else:
                eloss_psd_in_plane_file.write("%f %f %d\n"   % (0.1, 0.5, 0))
                eloss_psd_in_plane_file.write("%f %f %d\n\n" % (0.1, 1.0, 1))
                eloss_psd_in_plane_file.write("%f %f %d\n"   % (0.2, 1.0, 1))
                eloss_psd_in_plane_file.write("%f %f %d\n"   % (0.2, 0.5, 2))
        eloss_psd_in_plane_file.close()
        

        ### TOTAL VELOCITY LOSS ###
        print("Calculating total velocity loss.")
        logfile.write("Calculating total velocity loss.\n")
        # LOOP 
        all_vloss  = [traj.vloss for traj in trajs if traj.has_scattered] 
        one_vb     = [traj.vloss for traj in trajs if traj.has_scattered and traj.turn_pnts == 1]
        two_vb     = [traj.vloss for traj in trajs if traj.has_scattered and traj.turn_pnts == 3]
        mul_vb     = [traj.vloss for traj in trajs if traj.has_scattered and traj.turn_pnts >= 5]

        all_vf     = [length(traj.v_p_f) for traj in trajs if traj.has_scattered]
        one_vfb    = [length(traj.v_p_f) for traj in trajs if traj.has_scattered and traj.turn_pnts == 1]
        two_vfb    = [length(traj.v_p_f) for traj in trajs if traj.has_scattered and traj.turn_pnts == 3]
        mul_vfb    = [length(traj.v_p_f) for traj in trajs if traj.has_scattered and traj.turn_pnts >= 5]
        
        absorbed_vloss = [traj.vloss for traj in trajs if traj.has_adsorbed]

        # ANALYSIS
        all_vloss_hist, all_vloss_edges  = numpy.histogram(all_vloss,   bins=numbins(SCATTERED), range=(min(all_vloss), max(all_vloss)), density=True)
        one_vb_hist,     one_vb_edges     = numpy.histogram(one_vb,     bins=numbins(SCATTERED), range=(min(all_vloss), max(all_vloss)), density=True)
        two_vb_hist,     two_vb_edges     = numpy.histogram(two_vb,     bins=numbins(SCATTERED), range=(min(all_vloss), max(all_vloss)), density=True)
        mul_vb_hist,     mul_vb_edges     = numpy.histogram(mul_vb,     bins=numbins(SCATTERED), range=(min(all_vloss), max(all_vloss)), density=True)

        all_vf_hist,      all_vf_edges      = numpy.histogram(all_vf,     bins=numbins(SCATTERED), range=(min(all_vf), max(all_vf)), density=True)
        one_vfb_hist,     one_vfb_edges     = numpy.histogram(one_vfb,     bins=numbins(SCATTERED), range=(min(all_vf), max(all_vf)), density=True)
        two_vfb_hist,     two_vfb_edges     = numpy.histogram(two_vfb,     bins=numbins(SCATTERED), range=(min(all_vf), max(all_vf)), density=True)
        mul_vfb_hist,     mul_vfb_edges     = numpy.histogram(mul_vfb,     bins=numbins(SCATTERED), range=(min(all_vf), max(all_vf)), density=True)

        frac_one_vb = float(len(one_vb))/SCATTERED
        frac_two_vb = float(len(two_vb))/SCATTERED
        frac_mul_vb = float(len(mul_vb))/SCATTERED

        frac_one_vfb = float(len(one_vfb))/SCATTERED
        frac_two_vfb = float(len(two_vfb))/SCATTERED
        frac_mul_vfb = float(len(mul_vfb))/SCATTERED

        # OUTPUT 
        vloss_file = open("analysis/vloss.txt", "w")
        vloss_file.write("# vloss/Ang*fs^-1  all  single bounce  double bounce  multi bounce\n")
        for i in range(len(all_vloss_hist)):
        	vloss_file.write("%f %f %f %f %f\n" % (0.5*(all_vloss_edges[i]+all_vloss_edges[i+1]), all_vloss_hist[i], frac_one_vb*one_vb_hist[i], frac_two_vb*two_vb_hist[i], frac_mul_vb*mul_vb_hist[i]))
        vloss_file.close()

        # write final velocities
        v_final_file = open("analysis/all_final_v.txt", "w")
        v_final_file.write("# final v/Ang*fs^-1  all  single bounce  double bounce  multi bounce\n")
        for i in range(len(all_vf_hist)):
                v_final_file.write("%f %f %f %f %f\n" % (0.5*(all_vf_edges[i]+all_vf_edges[i+1]), all_vf_hist[i], frac_one_vfb*one_vfb_hist[i], frac_two_vfb*two_vfb_hist[i], frac_mul_vfb*mul_vfb_hist[i]))
        v_final_file.close()


        spatial_file = open("analysis/spatial_v.txt", "w")
        spatial_file.write("# x-pos y-pos v_scat\n")
        for traj in trajs:
            if traj.has_scattered:
                spatial_file.write("{:10.4f} {:10.4f} {:10.5f}\n".format(traj.r_p_f.x,traj.r_p_f.y,length(traj.v_p_f)))
        spatial_file.close()



        ### ANGULAR DISTRIBUTION ###
        print("Calculating angular velocity loss.")
        logfile.write("Calculating angular velocity loss.\n")
        # get trajectories that are within specular radius in azimuth direction
        velocity_collect  = []
        angle_collect   = []
        ps_dist_collect = []
        polar_scatt_azi_int_velocity = []
        polar_scatt_azi_int_angle  = []
        for traj in trajs:
        	if traj.has_scattered:
        		polar_scatt_azi_int_velocity.append(traj.vloss)
        		polar_scatt_azi_int_angle.append(traj.polar_f)

        		if traj.in_plane:
        			velocity_collect.append(traj.vloss)
        			angle_collect.append(traj.polar_f)
        			ps_dist_collect.append(traj.cl_appr)

        ang_dist_file_v = open("analysis/ang_res_vloss.txt", "w")
        ang_dist_mat_file_v = open("analysis/ang_res_vloss_matrix.txt", "w")
        occurence_file_v    = open("analysis/ang_res_occurrence_v.txt", "w")
        if len(velocity_collect) != 0:
                angle_vloss_hist, xedges, yedges = numpy.histogram2d(velocity_collect, angle_collect,  bins=(numbins(velocity_collect)), normed=False)
                #occurence_hist_v, occ_edges = numpy.histogram(angle_collect, bins=numpy.arange(91), normed=False)
                occurence_hist_v, occ_edges = numpy.histogram(angle_collect, bins=numpy.arange(91), density=False)
                	
                # OUTPUT
                for i in range(len(angle_vloss_hist)):
                	for j in range(len(angle_vloss_hist[i])):
                		ang_dist_file_v.write("%f %f %d\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), int(angle_vloss_hist[i][j])/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))))
                	ang_dist_file_v.write("\n")

                ang_dist_mat_file_v.write("# x-range describing velocity loss in Ang/fs (left to right) from %f to %f in steps of %f\n" % (0.5*(xedges[0]+xedges[1]), 0.5*(xedges[-2]+xedges[-1]), abs(xedges[0]-xedges[1])))
                ang_dist_mat_file_v.write("# y-range describing scattering angle in degrees (top to bottom) from %f to %f in steps of %f\n" % (0.5*(yedges[0]+yedges[1]), 0.5*(yedges[-2]+yedges[-1]), abs(yedges[0]-yedges[1])))
                ang_dist_mat_file_v.write("# specular scattering angle is %f degrees and detector radius is %f degrees\n" % (trajs[0].polar_i, SPECULAR_RADIUS))
                ang_dist_mat_file_v.write("0.0 "); [ang_dist_mat_file_v.write("%f " % (0.5*(xedges[i]+xedges[i+1]))) for i in range(len(xedges)-1)]; ang_dist_mat_file_v.write("\n")
                for i in range(len(angle_vloss_hist[0])):
                        ang_dist_mat_file_v.write("%f " % (0.5*(yedges[i]+yedges[i+1])))
                        for j in range(len(angle_vloss_hist)):
                                ang_dist_mat_file_v.write("%d " % angle_vloss_hist[j][i])
                        ang_dist_mat_file_v.write("\n")

                occurence_file_v.write("# number of trajs in plane in unit polar angle, unit azimuthal angle\n")
                for i in range(len(occurence_hist_v)):
                        this_angle = 0.5*(occ_edges[i]+occ_edges[i+1])
                        occurence_file_v.write("%f %f\n" % (this_angle, 1.0*occurence_hist_v[i]/sum(occurence_hist_v)))
        else:
                ang_dist_file_v.write("%f %f %d\n"   % (0.1, 0.5, 0))
                ang_dist_file_v.write("%f %f %d\n\n" % (0.1, 1.0, 1))
                ang_dist_file_v.write("%f %f %d\n"   % (0.2, 1.0, 1))
                ang_dist_file_v.write("%f %f %d\n"   % (0.2, 0.5, 2))

                occurence_file_v.write("%f %f\n" % (0.1, 0.1))
        	
        occurence_file_v.close()
        ang_dist_file_v.close()
        ang_dist_mat_file_v.close()

        # INTERGRATED OVER ALL AZIMUTH ANGLES #
        polar_scatt_azi_file_v = open("analysis/polar_scatt_azi_int_v.txt", "w")
        if len(polar_scatt_azi_int_velocity) != 0:
                polar_scatt_azi_int_hist_v, xedges, yedges = numpy.histogram2d(polar_scatt_azi_int_velocity, polar_scatt_azi_int_angle, bins=(numbins(polar_scatt_azi_int_velocity)), normed=False)
        	
                for i in range(len(polar_scatt_azi_int_hist_v)):
                        for j in range(len(polar_scatt_azi_int_hist_v[i])):
                                polar_scatt_azi_file_v.write("%f %f %d\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), int(polar_scatt_azi_int_hist_v[i][j])/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))))
                        polar_scatt_azi_file_v.write("\n")
        else:
                polar_scatt_azi_file_v.write("%f %f %d\n"   % (0.1, 0.5, 0))
                polar_scatt_azi_file_v.write("%f %f %d\n\n" % (0.1, 1.0, 1))
                polar_scatt_azi_file_v.write("%f %f %d\n"   % (0.2, 1.0, 1))
                polar_scatt_azi_file_v.write("%f %f %d\n"   % (0.2, 0.5, 2))

        polar_scatt_azi_file_v.close()

        ### SPHERICAL SYMMETRY ###
        # LOOP
        print("Calculating spherical symmetry.")
        logfile.write("Calculating spherical symmetry.\n")
        abs_azi = []
        rel_azi = []
        yvals = []
        this_azi = None
        for traj in trajs:
                if traj.has_scattered:
                        delta_azi = traj.azi_f-traj.azi_i
                        if -180 <= delta_azi <= 180:
                                rel_azi.append(delta_azi)
                        elif delta_azi < -180:
                                rel_azi.append(delta_azi+360)
                        elif delta_azi > 180:
                                rel_azi.append(delta_azi-360)
                        else:
                                sys.exit("Weird angle in spherical symmetry.")

                        abs_azi.append(traj.azi_f)
                        yvals.append(traj.polar_f)

        # ANALYSIS
        rel_spherical_hist_v, rel_xedges, yedges = numpy.histogram2d(rel_azi, yvals, bins=numbins(rel_azi), normed=False)
        abs_spherical_hist_v, abs_xedges, yedges = numpy.histogram2d(abs_azi, yvals, bins=numbins(abs_azi), normed=False)

        # OUTPUT
        spherical_file_rel_v = open("analysis/rel_spherical_symmetry_v.txt", "w")
        for i in range(len(rel_spherical_hist_v)):
                for j in range(len(rel_spherical_hist_v[i])):
                        spherical_file_rel_v.write("%f %f %d\n" % (0.5*(rel_xedges[i]+rel_xedges[i+1]), -0.5*(yedges[j]+yedges[j+1]), int(rel_spherical_hist_v[i][j])/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))))
                spherical_file_rel_v.write("\n")
        spherical_file_rel_v.close()

        spherical_file_abs_v = open("analysis/abs_spherical_symmetry_v.txt", "w")
        for i in range(len(abs_spherical_hist_v)):
                for j in range(len(abs_spherical_hist_v[i])):
                        spherical_file_abs_v.write("%f %f %d\n" % (0.5*(abs_xedges[i]+abs_xedges[i+1]), -0.5*(yedges[j]+yedges[j+1]), int(abs_spherical_hist_v[i][j])/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))))
                spherical_file_abs_v.write("\n")
        spherical_file_abs_v.close()









        ### SUMMARY ###
        # ANALYSIS
        energy_won = sum( [1 for e_final in all_eloss if e_final < 0] )

        # OUTPUT 
        out = open("analysis/Summary.txt", "w")
        out.write("Created by version %4.2f\n" % VERSION_ID)
        out.write("Scattered:   %d (%f%%)\n" % (SCATTERED,   100.*FRAC_SCATTERED))
        out.write("Absorbed:    %d (%f%%)\n" % (ABSORBED,    100.*FRAC_ABSORBED))
        out.write("Transmitted: %d (%f%%)\n\n" % (TRANSMITTED, 100.*FRAC_TRANSMITTED))

        out.write("%d (%f%%) of the scattered projectiles won kinetic energy.\n" % (energy_won, 100.*energy_won/SCATTERED))
        out.write("%f%% of scattered trajectories were within +-%f degrees in plane.\n" % (100.*len(in_plane_all_eloss)/SCATTERED, SPECULAR_RADIUS))
        out.write("%f%% of scattered trajectories were within +-%f degrees to specular scattering angle.\n\n" % (100.*len(spec_all_eloss)/SCATTERED, SPECULAR_RADIUS))

        out.write("Average energy loss of H-atoms in specular scattering angle %f eV.\n"   % numpy.mean(spec_all_eloss))
        try:
                out.write("Peak energy loss of H-atoms in specular scattering angle    %f eV.\n\n" % spec_all_eloss_edges[numpy.argmax(spec_all_eloss_hist)])
        except(UnboundLocalError):
                out.write("Peak energy loss of H-atoms in specular scattering angle    %s eV.\n\n" % "No atoms in specular scattering angle")
        out.write("Average energy loss of H-atoms reflected %f eV.\n"     % numpy.mean(all_eloss))
        out.write("Peak energy loss of H-atoms reflected    %f eV.\n\n"   % all_eloss_edges[numpy.argmax(all_eloss_hist)])
        out.write("Average energy loss of H-atoms in bulk   %f eV.\n" % numpy.mean(absorbed_eloss))
        out.write("Average energy loss of reflected H-atoms to ehps %f eV.\n" % numpy.mean(loss_to_ehps))
        out.write("Average energy loss of specularly reflected H-atoms to ehps %f eV.\n\n" % numpy.mean(loss_to_ehps_spec))
        
        out.write("%refl  %in bulk  %shot_thru  %E_won  %in_spec  avg_E_in_spec  peak_E_in_spec  avg_E_in_bulk  avg_E_refl  peak_E_refl  avg_ehp_loss  avg_ehp_loss_spec Trajs\n")
        try:
                out.write("%f %f %f %f %f %f %f %f %f %f %f %f %d\n" % (100.*FRAC_SCATTERED, 100.*FRAC_ABSORBED, 100.*FRAC_TRANSMITTED, 100.*energy_won/SCATTERED, 100.*len(spec_all_eloss)/SCATTERED, trajs[0].ekin_p_i-numpy.mean(spec_all_eloss), trajs[0].ekin_p_i-spec_all_eloss_edges[numpy.argmax(spec_all_eloss_hist)], trajs[0].ekin_p_i-numpy.mean(absorbed_eloss), trajs[0].ekin_p_i-numpy.mean(all_eloss), trajs[0].ekin_p_i-all_eloss_edges[numpy.argmax(all_eloss_hist)], numpy.mean(loss_to_ehps), numpy.mean(loss_to_ehps_spec), len(trajs)))
        except(UnboundLocalError):
	        out.write("%f %f %f %f %f %f %f %f %f %f %f %d\n" % (100.*FRAC_SCATTERED, 100.*FRAC_ABSORBED, 100.*FRAC_TRANSMITTED, 100.*energy_won/SCATTERED, 100.*len(spec_all_eloss)/SCATTERED, trajs[0].ekin_p_i-numpy.mean(spec_all_eloss), trajs[0].ekin_p_i-numpy.mean(absorbed_eloss), trajs[0].ekin_p_i-numpy.mean(all_eloss), trajs[0].ekin_p_i-all_eloss_edges[numpy.argmax(all_eloss_hist)], numpy.mean(loss_to_ehps), numpy.mean(loss_to_ehps_spec), len(trajs)))
        out.close()

	
# Average energy of H-atoms in bulk 0.03969474585 eV.
# Average energy of H-atoms reflected 1.93905890836 eV.
# Peak energy of H-atoms reflected 2.770000 eV
#


def graphene_bounce_events(trajs,logfile):

    print("Calculating graphene bounce events.")
    logfile.write("Calculating graphene bounce events.\n")

    fast_c    = open("analysis/component_fast.txt", "w")
    slow_c_sb = open("analysis/component_slow_single.txt", "w")
    slow_c_mb = open("analysis/component_slow_multi.txt", "w")

    for traj in trajs:
        if traj.has_scattered and traj.in_plane:
            if traj.turn_pnts == 1:
                if traj.cl_appr > 1.4:
                    fast_c.write("%f %f\n" % (traj.ekin_p_f/traj.ekin_p_i, traj.polar_f))

                if traj.cl_appr < 1.4:
                    slow_c_sb.write("%f %f\n" % (traj.ekin_p_f/traj.ekin_p_i, traj.polar_f))

            elif traj.turn_pnts > 1 and traj.cl_appr < 1.4:
                slow_c_mb.write("%f %f\n" % (traj.ekin_p_f/traj.ekin_p_i, traj.polar_f))

    slow_c_mb.close()
    slow_c_sb.close()
    fast_c.close()

def analyze_angles(trajs,logfile):
        for traj in trajs:	# List comprehension simply need too much time. This is ugly, but fast.
                if traj.in_plane and traj.has_scattered:

                        if traj.cl_appr < 1.4: # our structural parameter for the barrie # our structural parameter for the barrierr
                            outfile_string = "slow_component.log"
                        else:
                            outfile_string = "fast_component.log"

                        outfile = open(outfile_string,'a+')
                        if 14 <= traj.polar_f <= 16:
                            if 1.44 <= traj.ekin_p_f:
                                outfile.write("15+-1, 1.44 <= E_s: trajid {} and closest approach {}\n".format(traj.traj_id,traj.cl_appr))
                            if 0.960 <= traj.ekin_p_f < 1.44:
                                outfile.write("15+-1, 0.96 <= E_s < 1.44: trajid {} and closest approach {}\n".format(traj.traj_id,traj.cl_appr))
                            if traj.ekin_p_f < 0.960:
                                outfile.write("15+-1, E_s < 0.96: trajid {} and closest approach {}\n".format(traj.traj_id,traj.cl_appr))
                        if 29 <= traj.polar_f <= 31:
                            if 1.44 <= traj.ekin_p_f:
                                outfile.write("30+-1, 1.44 <= E_s: trajid {} and closest approach {}\n".format(traj.traj_id,traj.cl_appr))
                            if 0.960 <= traj.ekin_p_f < 1.44:
                                outfile.write("30+-1, 0.96 <= E_s < 1.44: trajid {} and closest approach {}\n".format(traj.traj_id,traj.cl_appr))
                            if traj.ekin_p_f < 0.960:
                                outfile.write("30+-1, E_s < 0.96: trajid {} and closest approach {}\n".format(traj.traj_id,traj.cl_appr))
                        if 44 <= traj.polar_f <= 46:
                            if 1.44 <= traj.ekin_p_f:
                                outfile.write("45+-1, 1.44 <= E_s: trajid {} and closest approach {}\n".format(traj.traj_id,traj.cl_appr))
                            if 0.960 <= traj.ekin_p_f < 1.44:
                                outfile.write("45+-1, 0.96 <= E_s < 1.44: trajid {} and closest approach {}\n".format(traj.traj_id,traj.cl_appr))
                            if traj.ekin_p_f < 0.960:
                                outfile.write("45+-1, E_s < 0.96: trajid {} and closest approach {}\n".format(traj.traj_id,traj.cl_appr))
                        if 59 <= traj.polar_f <= 61:
                            if 1.44 <= traj.ekin_p_f:
                                outfile.write("60+-1, 1.44 <= E_s: trajid {} and closest approach {}\n".format(traj.traj_id,traj.cl_appr))
                            if 0.960 <= traj.ekin_p_f < 1.44:
                                outfile.write("60+-1, 0.96 <= E_s < 1.44: trajid {} and closest approach {}\n".format(traj.traj_id,traj.cl_appr))
                            if traj.ekin_p_f < 0.960:
                                outfile.write("60+-1, E_s < 0.96: trajid {} and closest approach {}\n".format(traj.traj_id,traj.cl_appr))

                        outfile.close()


def get_traj(trajs,logfile):

    def_nrg_file  = open(defnrgname, "w")

    ntrajs        = len(trajs)
    scattered_ctr = 0
    absorbed_ctr  = 0
    slow_comp_ctr = 0
    fast_comp_ctr = 0
    back_scat_ctr = 0
    traj_after    = 0
    traj_before   = 0
    cl_appr_dist  = 1.4

    for traj in trajs:
        if traj.has_scattered:
            scattered_ctr += 1

            if traj.cl_appr <= cl_appr_dist: # to get slow component
                traj_after += 1
                if traj.turn_pnts == 1:
                    slow_comp_ctr += 1

                    def_nrg_file.write(traj.traj_id)
                    def_nrg_file.write(' ')

                    print("Analyze slow component in traj {}".format(traj.traj_id))
                    logfile.write("Analyze slow component in traj {}\n".format(traj.traj_id))

            else: # traj.cl_appr > 1.4
                traj_before += 1
                if traj.turn_pnts == 1:
                    fast_comp_ctr += 1

                    print("Analyze fast component in traj {}".format(traj.traj_id))
                    logfile.write("Analyze fast component in traj {}\n".format(traj.traj_id))

            if traj.azi_f < -100: # for backscattering
                back_scat_ctr += 1

                print("potential backscattering in traj {}".format(traj.traj_id))
                logfile.write("potential backscattering in traj {}\n".format(traj.traj_id))

        else:
            absorbed_ctr += 1
            print("particle adsorbed in traj {}".format(traj.traj_id))
            logfile.write("particle adsorbed in traj {}\n".format(traj.traj_id))

    print("trajs with scattering after barrier (single bounce): {} out of {} scattered traj ({}%)".format(slow_comp_ctr,scattered_ctr,float(slow_comp_ctr)*100/float(scattered_ctr)))
    logfile.write("trajs with scattering after barrier (single bounce): {} out of {} scattered traj ({}%)\n".format(slow_comp_ctr,scattered_ctr,float(slow_comp_ctr)*100/float(scattered_ctr)))

    print("trajs with adsorption: {} out of total {} traj ({}%)".format(absorbed_ctr,ntrajs,float(absorbed_ctr)*float(100)/float(ntrajs)))
    logfile.write("trajs with adsorption: {} out of total {} traj ({}%)\n".format(absorbed_ctr,ntrajs,float(absorbed_ctr)*float(100)/float(ntrajs)))

    print("trajs with scattering before barrier: {} out of {} scattered traj ({}%)".format(traj_before,scattered_ctr,float(traj_before)*float(100)/float(scattered_ctr)))
    logfile.write("trajs with scattering before barrier: {} out of {} scattered traj ({}%)\n".format(traj_before,scattered_ctr,float(traj_before)*float(100)/float(scattered_ctr)))

    print("trajs with scattering after barrier: {} out of {} scattered traj ({}%)".format(traj_after,scattered_ctr,float(traj_after)*float(100)/float(scattered_ctr)))
    logfile.write("trajs with scattering after barrier: {} out of {} scattered traj ({}%)\n".format(traj_after,scattered_ctr,float(traj_after)*float(100)/float(scattered_ctr)))



def rat_analysis(trajs,logfile):

        # 2D ANGULAR DISTRIBUTION RAT
        print("Calculate 2D angular distribution for RAT experiment")
        angle_rat_collect = []
        efrac_rat_collect = []
        AZIMUTHAL_ANGLE_RAT_D1 = 13.5
        AZIMUTHAL_ANGLE_RAT_D2 = -13.5

       # DOMAIN 1 with azi_i = +13.5 degree
        for traj in trajs:
            if traj.azi_i == AZIMUTHAL_ANGLE_RAT_D1:
                if traj.azi_f < 0:
                    inv_azi = AZIMUTHAL_ANGLE_RAT_D1 - 180
                    if traj.has_scattered and abs(inv_azi - traj.azi_f) < SPECULAR_RADIUS:
                        angle_rat_collect.append(-traj.polar_f)
                        efrac_rat_collect.append(traj.efrac)
                    elif traj.has_scattered and abs(AZIMUTHAL_ANGLE_RAT_D1 - traj.azi_f) < SPECULAR_RADIUS:
                        angle_rat_collect.append(traj.polar_f)
                        efrac_rat_collect.append(traj.efrac)
                if traj.azi_f > 0:
                    inv_azi = AZIMUTHAL_ANGLE_RAT_D1 + 180
                    if traj.has_scattered and abs(inv_azi - traj.azi_f) < SPECULAR_RADIUS:
                        angle_rat_collect.append(-traj.polar_f)
                        efrac_rat_collect.append(traj.efrac)
                    elif traj.has_scattered and abs(AZIMUTHAL_ANGLE_RAT_D1 - traj.azi_f) < SPECULAR_RADIUS:
                        angle_rat_collect.append(traj.polar_f)
                        efrac_rat_collect.append(traj.efrac)

       # DOMAIN 2 with azi_i = -13.5 degree
        for traj in trajs:
            if traj.azi_i == AZIMUTHAL_ANGLE_RAT_D2:
                if traj.azi_f < 0:
                    inv_azi = AZIMUTHAL_ANGLE_RAT_D2 - 180
                    if traj.has_scattered and abs(inv_azi - traj.azi_f) < SPECULAR_RADIUS:
                        angle_rat_collect.append(-traj.polar_f)
                        efrac_rat_collect.append(traj.efrac)
                    elif traj.has_scattered and abs(AZIMUTHAL_ANGLE_RAT_D2 - traj.azi_f) < SPECULAR_RADIUS:
                        angle_rat_collect.append(traj.polar_f)
                        efrac_rat_collect.append(traj.efrac)
                if traj.azi_f > 0:
                    inv_azi = AZIMUTHAL_ANGLE_RAT_D2 + 180
                    if traj.has_scattered and abs(inv_azi - traj.azi_f) < SPECULAR_RADIUS:
                        angle_rat_collect.append(-traj.polar_f)
                        efrac_rat_collect.append(traj.efrac)
                    elif traj.has_scattered and abs(AZIMUTHAL_ANGLE_RAT_D2 - traj.azi_f) < SPECULAR_RADIUS:
                        angle_rat_collect.append(traj.polar_f)
                        efrac_rat_collect.append(traj.efrac)

        ang_dist_nrg_rat_ang          = open("analysis/2d-ang-dist_rat.txt", "w")
        ang_dist_nrg_rat_ang_norm_sum = open("analysis/2d-ang-dist_rat_norm_sum.txt", "w")
        ang_dist_nrg_rat_ang_norm_max = open("analysis/2d-ang-dist_rat_norm_max.txt", "w")
        
        ang_dist_nrg_rat_ang_norm_bin_sum = open("analysis/2d-ang-dist_rat_norm_bin_sum.txt", "w")
        ang_dist_nrg_rat_ang_norm_bin_max = open("analysis/2d-ang-dist_rat_norm_bin_max.txt", "w")

        bin_arr = []

        if len(efrac_rat_collect) != 0 and len(angle_rat_collect) != 0:
            #angle_efrac_hist, xedges, yedges = numpy.histogram2d(efrac_all_collect, angle_all_collect,  bins=(numbins(efrac_all_collect)), density=False)
            #angle_efrac_hist, xedges, yedges = numpy.histogram2d(efrac_all_collect, angle_all_collect, bins=BINS, range=[[0, 1.1],[-90, 90]], density=False) # fixed bin size of 2 deg; try 36,72,180
            angle_efrac_hist, xedges, yedges = numpy.histogram2d(efrac_rat_collect, angle_rat_collect, bins=BINS, range=[[0, 1.1],[ANGLE_MIN, ANGLE_MAX]], density=False) # fixed bin size of 2 deg; try 36,72,180

            for i in range(len(angle_efrac_hist)):
                for j in range(len(angle_efrac_hist[i])):
                    val_abs = int(angle_efrac_hist[i][j]/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))) # needs to be divided by sin(x) to correct geometry of experiment.
                    bin_arr.append(float(val_abs))
                    val_sum = float(val_abs)/float(angle_efrac_hist.sum()) # current bin divided by sum of all values to get "flux"
                    val_max = float(val_abs)/float(angle_efrac_hist.max()) # current bin divided by max value to get "flux"
                    ang_dist_nrg_rat_ang.write("%f %f %d\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), val_abs)) # write data with bins
                    ang_dist_nrg_rat_ang_norm_sum.write("%f %f %f\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), val_sum)) # write data with area integrtated "flux""
                    ang_dist_nrg_rat_ang_norm_max.write("%f %f %f\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), val_max)) # write data with normalized "flux""

            for i in range(len(angle_efrac_hist)):
                for j in range(len(angle_efrac_hist[i])):
                    val_abs = int(angle_efrac_hist[i][j]/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi)))
                    val_sum = float(val_abs)/float(sum(bin_arr))
                    val_max = float(val_abs)/float(max(bin_arr))
                    ang_dist_nrg_rat_ang_norm_bin_sum.write("%f %f %f\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), val_sum))
                    ang_dist_nrg_rat_ang_norm_bin_max.write("%f %f %f\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), val_max))

        ang_dist_nrg_rat_ang.close()
        ang_dist_nrg_rat_ang_norm_sum.close()
        ang_dist_nrg_rat_ang_norm_max.close()

        ang_dist_nrg_rat_ang_norm_bin_sum.close()
        ang_dist_nrg_rat_ang_norm_bin_max.close()


def ion_imaging_analysis(trajs,logfile):

        print("Calculate 2D angular distribution for ion imaging experiment")
        logfile.write("Calculate 2D angular distribution for ion imaging experiment\n")
        angle_ion_collect = [] # all polar angle from -90 to 90 in deg
        efrac_ion_collect = [] # all E_s / E_i in eV
 
        for traj in trajs:
            if traj.azi_f < 0:
                inv_azi = AZIMUTHAL_ANGLE - 180
                if traj.has_scattered and abs(inv_azi - traj.azi_f) < ION_IMAGING_AZI:
                    angle_ion_collect.append(-traj.polar_f)
                    efrac_ion_collect.append(traj.efrac)
                elif traj.has_scattered and abs(AZIMUTHAL_ANGLE - traj.azi_f) < ION_IMAGING_AZI:
                    angle_ion_collect.append(traj.polar_f)
                    efrac_ion_collect.append(traj.efrac)
            if traj.azi_f > 0:
                inv_azi = AZIMUTHAL_ANGLE + 180
                if traj.has_scattered and abs(inv_azi - traj.azi_f) < ION_IMAGING_AZI:
                    angle_ion_collect.append(-traj.polar_f)
                    efrac_ion_collect.append(traj.efrac)
                elif traj.has_scattered and abs(AZIMUTHAL_ANGLE - traj.azi_f) < ION_IMAGING_AZI:
                    angle_ion_collect.append(traj.polar_f)
                    efrac_ion_collect.append(traj.efrac)
 
        ang_dist_nrg_ang_ion          = open("analysis/2d-ang-dist_ion_imaging.txt", "w")
        ang_dist_nrg_ang_ion_norm_sum = open("analysis/2d-ang-dist_ion_imaging_norm_sum.txt", "w")
        ang_dist_nrg_ang_ion_norm_max = open("analysis/2d-ang-dist_ion_imaging_norm_max.txt", "w")
        
        ang_dist_nrg_ang_ion_norm_bin_sum = open("analysis/2d-ang-dist_ion_imaging_norm_bin_sum.txt", "w")
        ang_dist_nrg_ang_ion_norm_bin_max = open("analysis/2d-ang-dist_ion_imaging_norm_bin_max.txt", "w")

        bin_arr = []

        if len(efrac_ion_collect) != 0 and len(angle_ion_collect) != 0:
            angle_efrac_hist, xedges, yedges = numpy.histogram2d(efrac_ion_collect, angle_ion_collect, bins=BINS, range=[[0, 1.1],[ANGLE_MIN, ANGLE_MAX]], density=False) # fixed bin size of 2 deg; try 36,72,180

            for i in range(len(angle_efrac_hist)):
                for j in range(len(angle_efrac_hist[i])):
                    val_abs = int(angle_efrac_hist[i][j]/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))) # needs to be divided by sin(x) to correct geometry of experiment.
                    bin_arr.append(float(val_abs))
                    val_sum = float(val_abs)/float(angle_efrac_hist.sum()) # current bin divided by sum of all values to get "flux"; area normed
                    val_max = float(val_abs)/float(angle_efrac_hist.max()) # current bin divided by maximum value of all values to get "flux", max value normed
                    ang_dist_nrg_ang_ion.write("%f %f %d\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), val_abs)) # write data with bins
                    ang_dist_nrg_ang_ion_norm_sum.write("%f %f %f\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), val_sum)) # write data with area normed "flux"
                    ang_dist_nrg_ang_ion_norm_max.write("%f %f %f\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), val_max)) # write data with maximum value normed "flux"

            for i in range(len(angle_efrac_hist)):
                for j in range(len(angle_efrac_hist[i])):
                    val_abs = int(angle_efrac_hist[i][j]/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi)))
                    val_sum = float(val_abs)/float(sum(bin_arr)) # current bin divided by sum of all bins
                    val_max = float(val_abs)/float(max(bin_arr)) # current bin divided by maximum value of all bins
                    ang_dist_nrg_ang_ion_norm_bin_sum.write("%f %f %f\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), val_sum))
                    ang_dist_nrg_ang_ion_norm_bin_max.write("%f %f %f\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), val_max))

        ang_dist_nrg_ang_ion.close()
        ang_dist_nrg_ang_ion_norm_sum.close()
        ang_dist_nrg_ang_ion_norm_max.close()

        print(angle_efrac_hist.max())
        print(angle_efrac_hist.sum())

        ang_dist_nrg_ang_ion_norm_bin_sum.close()
        ang_dist_nrg_ang_ion_norm_bin_max.close()


def cmd_analysis(trajs,logfile):
        print("Calculate 2D angular distribution for cMD simulations\n")
        logfile.write("Calculate 2D angular distribution for cMD simulations\n")
        angle_cmd_collect = [] # all polar angle from -90 to 90 in deg
        efrac_cmd_collect = [] # all E_s / E_i in eV
 
        cmd_azi = 0
        cmd_azi_det = 90
 
        for traj in trajs:
            if traj.azi_f < 0:
                inv_azi = cmd_azi - 180
                if traj.has_scattered and abs(inv_azi - traj.azi_f) < cmd_azi_det:
                    angle_cmd_collect.append(-traj.polar_f)
                    efrac_cmd_collect.append(traj.efrac)
                elif traj.has_scattered and abs(cmd_azi - traj.azi_f) < cmd_azi_det:
                    angle_cmd_collect.append(traj.polar_f)
                    efrac_cmd_collect.append(traj.efrac)
            if traj.azi_f > 0:
                inv_azi = cmd_azi + 180
                if traj.has_scattered and abs(inv_azi - traj.azi_f) < cmd_azi_det:
                    angle_cmd_collect.append(-traj.polar_f)
                    efrac_cmd_collect.append(traj.efrac)
                elif traj.has_scattered and abs(cmd_azi - traj.azi_f) < cmd_azi_det:
                    angle_cmd_collect.append(traj.polar_f)
                    efrac_cmd_collect.append(traj.efrac)

        ang_dist_nrg_ang_cmd      = open("analysis/2d-ang-dist_cmd.txt", "w")
        ang_dist_nrg_ang_cmd_norm = open("analysis/2d-ang-dist_cmd_norm.txt", "w")

        if len(efrac_cmd_collect) != 0 and len(angle_cmd_collect) != 0:
            angle_efrac_hist, xedges, yedges = numpy.histogram2d(efrac_cmd_collect, angle_cmd_collect, bins=BINS, range=[[0, 1.1],[ANGLE_MIN, ANGLE_MAX]], density=False)

            for i in range(len(angle_efrac_hist)):
                for j in range(len(angle_efrac_hist[i])):
                    abs_val = int(angle_efrac_hist[i][j]/abs(numpy.sin(0.5*(yedges[j]+yedges[j+1])/360*2*numpy.pi))) # needs to be divided by sin(x) to correct geometry of experiment.
                    rel_val = float(abs_val)/float(angle_efrac_hist.sum()) # current bin devided by sum of all bins to get "flux"
                    ang_dist_nrg_ang_cmd.write("%f %f %d\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), abs_val)) # write data with bins
                    ang_dist_nrg_ang_cmd_norm.write("%f %f %f\n" % (0.5*(xedges[i]+xedges[i+1]), 0.5*(yedges[j]+yedges[j+1]), rel_val)) # write data with "flux"

        ang_dist_nrg_ang_cmd.close()
        ang_dist_nrg_ang_cmd_norm.close()


def get_movies(trajs,logfile):

    for traj in trajs:
        if traj.has_scattered:
            #if 1.80 <= traj.
            if 0.09 <= traj.eloss <= 0.11:
                print(traj.traj_id)



###### SCRIPT ######

# open logfile
logfile = open(logfilename, 'w')

print("Screen output will be automatically written to {}!".format(logfilename))

print("Created by version %4.2f" % VERSION_ID)
logfile.write("Created by version %4.2f\n" % VERSION_ID)


### READ IN TRAJS ###
traj_collection, SCATTERED, ABSORBED, TRANSMITTED = initialize(inpname,logfile)

### CALCULATE USEFUL CONSTANTS ###
NTRAJS = len(traj_collection)
FRAC_SCATTERED = float(SCATTERED)/NTRAJS
FRAC_ABSORBED = float(ABSORBED)/NTRAJS
FRAC_TRANSMITTED = float(TRANSMITTED)/NTRAJS

### RAT ###
#rat_analysis(traj_collection,logfile)

### ION IMAGING ###
ion_imaging_analysis(traj_collection,logfile)

### constrained MD ###
#cmd_analysis(traj_collection,logfile)

### OUTPUT ###
analyze(traj_collection,logfile)



### H@Gr related functions
#graphene_bounce_events(traj_collection,logfile)
#analyze_angles(traj_collection,logfile)
#get_traj(traj_collection,logfile) # get number of trajs for backscattering
#get_movies(traj_collection,logfile) # get traj id for movies


logfile.close()
