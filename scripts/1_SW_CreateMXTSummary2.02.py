#!/usr/bin/env python

# original version 2.00 was done by Marvin Kammler
# version 2.01 is modified by Sebastian Wille
# version 2.02 is modified by Sebastian Wille

# Modifications v2.01:
# trajid added in Traj class
# logfile will be written
# set readfile name
# depending on final polar angle additional data files are written

# Modifications v2.02:
# initial and final velocities will also be given
# added time at closest approach

# intention: analyze all traj/mxt_fin files and create the MXt2Summary file

# use like: python <scriptname>


# 2DO:
# merge this script with the analysis script and in the end use matplotlib to generate all plots using a single script



import os, sys, glob 

# set names for output and log file
outname     = "MXT2Summary.txt"
logfilename = "CreateMXTSummary.log"

# add range of scattered angle to look at

############# NO CHANGES BEYOND THIS LINE ######################################

VERSION_ID = 2.02

class Traj:
	def __init__(self, fname, ekin_p_i, ekin_l_i, epot_i, etotal_i,\
			r_p_i, v_p_i, polar_i, azi_i, ekin_p_f, ekin_l_f, epot_f,   \
                        etotal_f, r_p_f, v_p_f, polar_f, azi_f, time, turn_pnts, \
			cl_appr, cl_appr_t, r_p_min, traj_id):
		self.fname     = fname
		self.ekin_p_i  = ekin_p_i
		self.ekin_l_i  = ekin_l_i
		self.epot_i    = epot_i
		self.etotal_i  = etotal_i
		self.r_p_i     = r_p_i
		self.v_p_i     = v_p_i
		self.polar_i   = polar_i
		self.azi_i     = azi_i
		self.ekin_p_f  = ekin_p_f
		self.ekin_l_f  = ekin_l_f
		self.epot_f    = epot_f
		self.etotal_f  = etotal_f
		self.r_p_f     = r_p_f
		self.v_p_f     = v_p_f
		self.polar_f   = polar_f
		self.azi_f     = azi_f
		self.time      = time
		self.turn_pnts = turn_pnts
		self.cl_appr   = cl_appr
		self.cl_appr_t = cl_appr_t
		self.r_p_min   = r_p_min
                self.traj_id   = traj_id

	def containsNaN(self):
		return "NaN" in self.ekin_p_i or "NaN" in self.ekin_l_i or "NaN" in self.epot_i or "NaN" in self.etotal_i or \
			"NaN" in self.ekin_p_f or "NaN" in self.ekin_l_f or "NaN" in self.epot_f or "NaN" in self.etotal_f

	def __str__(self):
		return "   fname " + str(self.fname) +\
			"\nekin_p_i " + str(self.ekin_p_i) +\
                        "\nekin_l_i " + str(self.ekin_l_i) +\
                        "\nepot_i " + str(self.epot_i) +\
                        "\netotal_i " + str(self.etotal_i) +\
                        "\nr_p_i " + str(self.r_p_i) +\
                        "\nv_p_i " + str(self.v_p_i) +\
                        "\npolar_i " + str(self.polar_i) +\
			"\nazi_i " + str(self.azi_i) +\
                        "\nekin_p_f " + str(self.ekin_p_f) +\
                        "\nekin_l_f " + str(self.ekin_l_f) +\
                        "\nepot_f " + str(self.epot_f) +\
                        "\netotal_f " + str(self.etotal_f) +\
                        "\nr_p_f " + str(self.r_p_f) +\
                        "\nv_p_f " + str(self.v_p_f) +\
                        "\npolar_f " + str(self.polar_f) +\
                        "\nazi_f " + str(self.azi_f) +\
                        "\ntime " + str(self.time) +\
                        "\nturn_pnts " + str(self.turn_pnts) +\
			"\ncl_appr " + str(self.cl_appr) +\
			"\ncl_appr_t " + str(self.cl_appr_t) +\
			"\nr_p_min " + str(self.r_p_min) +\
                        "\ntraj_id " + str(self.traj_id)



#def read_traj_id():
#        os.chdir("traj/")
#        traj_dir = os.getcwd()
#        folder_list = sorted(glob.glob('mxt_*'))
#        num_folders = len(folder_list)
#        print "Reading traj ids..."
#        traj_id = []
#        for folder in folder_list:
#                traj_id.append(folder[7:15])
#
#        os.chdir("../")
#        return traj_id

def file_len(fname):
        i = -1 # to handle empty files
        with open(fname) as f:
                for i, l in enumerate(f):
                        pass
        return i + 1



def read_in_mxt_fins(logfile):
	os.chdir("traj/")
	traj_dir = os.getcwd()
	folder_list = sorted(glob.glob('mxt_*'))
	#folder_list = [folder for folder in os.listdir(traj_dir)]
	num_folders = len(folder_list)
#	traj_list = num_folders*[None]
        traj_list = []
	print("Reading {} trajs...".format(num_folders))
        logfile.write("Reading {} trajs...\n".format(num_folders))
	counter = 0
        traj_id = []
	for folder in folder_list:
		if (counter % (num_folders/10) == 0):
			print 100*counter/num_folders+1, "%"
                        logfile.write("{}%\n".format(100*counter/num_folders+1))
		infile = open(folder, 'r')						        	####    Reference Trajectory Output   ####

                traj_id = folder[7:15]

                checkline = file_len(folder) # check if number of entries match finished traj file

                if checkline == 23:
		        for line in infile:
			        try:
			        	ekin_p_i = line.strip(' \n\t\r').split()[-1]; line=infile.next()	#	ekin_p_i =      1.9200000
				        ekin_l_i = line.strip(' \n\t\r').split()[-1]; line=infile.next()	#	ekin_l_i =      0.9553543
				        epot_i   = line.strip(' \n\t\r').split()[-1]; line=infile.next()	#	epot_i   =   -392.2094094
				        etotal_i = line.strip(' \n\t\r').split()[-1]; line=infile.next()	#       etotal_i =   -369.3305490
			        	r_p_i    = line.strip(' \n\t\r').split()[-3:] ; line=infile.next()	#	r_i      =     30.1309358    32.4013593     3.5000000
			        	v_p_i    = line.strip(' \n\t\r').split()[-3:] ; line=infile.next()	#	v_i      =     0.0000000     0.0000000     -0.1934888
		        		polar_i  = line.strip(' \n\t\r').split()[-1]; line=infile.next()	#	polar_i  =     50.0000000
	        			azi_i    = line.strip(' \n\t\r').split()[-1]; line=infile.next()	#	azi_i    =      0.0000000
		        		line=infile.next()							#
			        	ekin_p_f = line.strip(' \n\t\r').split()[-1]; line=infile.next()	#       ekin_p_f =      1.5942738
		        		ekin_l_f = line.strip(' \n\t\r').split()[-1]; line=infile.next()	#	ekin_l_f =      0.9610213
       	                	        epot_f   = line.strip(' \n\t\r').split()[-1]; line=infile.next()	#	epot_f   =   -391.8832698
       	        	                etotal_f = line.strip(' \n\t\r').split()[-1]; line=infile.next()	#	etotal_f =   -369.3305177
       	        	                r_p_f    = line.strip(' \n\t\r').split()[-3:]; line=infile.next()	#	r_f      =     33.1131979    32.3776858     1.3585263
       	        	                v_p_f    = line.strip(' \n\t\r').split()[-3:]; line=infile.next()	#	r_f      =      0.0026410     0.0117354     0.1513796
       	 		                polar_f  = line.strip(' \n\t\r').split()[-1]; line=infile.next()	#	polar_f  =    122.6747121
       		       	                azi_f    = line.strip(' \n\t\r').split()[-1]; line=infile.next()	#	azi_f    =     -5.5773256
				        line=infile.next()
				        time     = line.strip(' \n\t\r').split()[-1]; line=infile.next()	#	time     =     20.0000000
				        turn_pnts= line.strip(' \n\t\r').split()[-1]; line=infile.next()	#	turn_pnts =           1
				        cl_appr  = line.strip(' \n\t\r').split()[-1]; line=infile.next()        #       cl_appr  =      0.9694690
				        cl_appr_t  = line.strip(' \n\t\r').split()[-1]; line=infile.next()        #       cl_appr_time  =          278
				        r_p_min  = line.strip(' \n\t\r').split()[-3:]				#	r_min_p  =     51.3246421    33.1600141     0.9996723

				        traj = Traj(folder, ekin_p_i, ekin_l_i, epot_i, etotal_i,  \
                      			    r_p_i, v_p_i, polar_i, azi_i, ekin_p_f, ekin_l_f, epot_f,\
                       			    etotal_f, r_p_f, v_p_f, polar_f, azi_f, time, turn_pnts,\
					    cl_appr, cl_appr_t, r_p_min, traj_id)
			        except:
				        print "Error in file ", folder, " in line ", line
                                        logfile.write("Error in file {} in line {}\n".format(folder,line))
				        sys.exit()	

		        #infile.close()
		        #traj_list[counter] = traj
                        traj_list.append(traj)
                        counter += 1

                else:
                        print "Skipping unfinished traj", folder
                        logfile.write("Skipping unfinished traj {}\n".format(folder))
                        #break
                        continue

                infile.close()

        
	os.chdir("../")
	return traj_list


def write_summary(logfile, outfile_name, traj_list):
	outfile = open(outfile_name, "w")
	outfile.write("# traj_id E_kin_p   E_kin_l        E_pot      E_total r_p(    x,        y,         z) v_p(    x,        y,         z)      polar       azi   E_kin_p  E_kin_l       E_pot       E_total     r_p(x,       y,         z) v_p(x,       y,         z)        polar       azi     simtime turn_pnts   cl_appr   cl_appr_t   r_p_min\n")
	for traj in traj_list:
		try:
			outfile.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" \
				% ( traj.traj_id, traj.ekin_p_i, traj.ekin_l_i, traj.epot_i, traj.etotal_i, traj.r_p_i[0], traj.r_p_i[1], traj.r_p_i[2], \
                                traj.v_p_i[0], traj.v_p_i[1], traj.v_p_i[2], \
				traj.polar_i, traj.azi_i, traj.ekin_p_f, traj.ekin_l_f, traj.epot_f, traj.etotal_f, traj.r_p_f[0], \
				traj.r_p_f[1], traj.r_p_f[2], traj.v_p_f[0], traj.v_p_f[1], traj.v_p_f[2], traj.polar_f, traj.azi_f, \
                                traj.time, traj.turn_pnts, traj.cl_appr, traj.cl_appr_t, \
				traj.r_p_min[0], traj.r_p_min[1], traj.r_p_min[2]))
		except (IndexError):
			outfile.close()
                        logfile.write("ERROR: There's something wrong with file " + str(traj))
                        logfile.close()
			sys.exit("ERROR: There's something wrong with file " + str(traj))
			

	outfile.close()

# the following is not needed anymore since empty or unfinished files are skipped; might be useful for other applications
def remove_unfinished_traj(logfile):
    os.chdir("traj/")
    traj_dir = os.getcwd()
    folder_list = sorted(glob.glob('mxt_*'))
    num_folders = len(folder_list)
    print("Found {} trajs".format(num_folders))
    logfile.write("Found {} trajs\n".format(num_folders))
    count_files = 0
    for folder in folder_list:
        inpfile = open(folder, 'r')
        count_line = 0
        for line in inpfile:
            count_line += 1
        if count_line != 20:
            count_files += 1
            print("WARNING: the unfinished trajectory {} was removed".format(folder))
            logfile.write("WARNING: the unfinished trajectory {} was removed\n".format(folder))
            os.remove(folder)

        inpfile.close()

    print("In total, {} trajectories were removed".format(count_files))
    logfile.write("In total, {} trajectories were removed\n".format(count_files))

    os.chdir("../")

def write_traj_to_file(traj,this_outfile):
        outfile = open(this_outfile, "w")
        outfile.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" \
                       % ( traj.traj_id, traj.ekin_p_i, traj.ekin_l_i, traj.epot_i, traj.etotal_i, traj.r_p_i[0], \
                       traj.r_p_i[1], traj.r_p_i[2], traj.v_p_i[0], traj.v_p_i[1], traj.v_p_i[2], traj.polar_i, \
                       traj.azi_i, traj.ekin_p_f, traj.ekin_l_f, traj.epot_f, traj.etotal_f, traj.r_p_f[0], 
                       traj.r_p_f[1], traj.r_p_f[2], traj.v_p_f[0], traj.v_p_f[1], traj.v_p_f[2], traj.polar_f, \
                       traj.azi_f, traj.time, traj.turn_pnts, traj.cl_appr, traj.cl_appr_t, traj.r_p_min[0], \
                       traj.r_p_min[1], traj.r_p_min[2]))
        outfile.close()


# 2DO: include write function to spare unnecessary code!!
def write_statistics(traj_list):
        outfile_15 = "MXT2Summary_15.txt"
        outfile_30 = "MXT2Summary_30.txt"
        outfile_45 = "MXT2Summary_45.txt"
        outfile_60 = "MXT2Summary_60.txt"
        for traj in traj_list:
                if 7.5 <= float(traj.polar_f) <= 22.5:
                        write_traj_to_file(traj,outfile_15)
                        #outfile_15.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" \
                        #        % ( traj.traj_id, traj.ekin_p_i, traj.ekin_l_i, traj.epot_i, traj.etotal_i, traj.r_p_i[0], traj.r_p_i[1], traj.r_p_i[2], \
                        #        traj.polar_i, traj.azi_i, traj.ekin_p_f, traj.ekin_l_f, traj.epot_f, traj.etotal_f, traj.r_p_f[0], \
                        #        traj.r_p_f[1], traj.r_p_f[2], traj.polar_f, traj.azi_f, traj.time, traj.turn_pnts, traj.cl_appr, \
                        #        traj.r_p_min[0], traj.r_p_min[1], traj.r_p_min[2]))
                if 22.5 <= float(traj.polar_f) <= 37.5:
                        write_traj_to_file(traj,outfile_30)
                        #outfile_30.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" \
                        #        % ( traj.traj_id, traj.ekin_p_i, traj.ekin_l_i, traj.epot_i, traj.etotal_i, traj.r_p_i[0], traj.r_p_i[1], traj.r_p_i[2], \
                        #        traj.polar_i, traj.azi_i, traj.ekin_p_f, traj.ekin_l_f, traj.epot_f, traj.etotal_f, traj.r_p_f[0], \
                        #        traj.r_p_f[1], traj.r_p_f[2], traj.polar_f, traj.azi_f, traj.time, traj.turn_pnts, traj.cl_appr, \
                        #        traj.r_p_min[0], traj.r_p_min[1], traj.r_p_min[2]))
                if 37.5 <= float(traj.polar_f) <= 52.5:
                        write_traj_to_file(traj,outfile_45)
                        #outfile_45.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" \
                        #        % ( traj.traj_id, traj.ekin_p_i, traj.ekin_l_i, traj.epot_i, traj.etotal_i, traj.r_p_i[0], traj.r_p_i[1], traj.r_p_i[2], \
                        #        traj.polar_i, traj.azi_i, traj.ekin_p_f, traj.ekin_l_f, traj.epot_f, traj.etotal_f, traj.r_p_f[0], \
                        #        traj.r_p_f[1], traj.r_p_f[2], traj.polar_f, traj.azi_f, traj.time, traj.turn_pnts, traj.cl_appr, \
                        #        traj.r_p_min[0], traj.r_p_min[1], traj.r_p_min[2]))
                if 52.5 <= float(traj.polar_f) <= 67.5:
                        write_traj_to_file(traj,outfile_60)
                        #outfile_60.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" \
                        #        % ( traj.traj_id, traj.ekin_p_i, traj.ekin_l_i, traj.epot_i, traj.etotal_i, traj.r_p_i[0], traj.r_p_i[1], traj.r_p_i[2], \
                        #        traj.polar_i, traj.azi_i, traj.ekin_p_f, traj.ekin_l_f, traj.epot_f, traj.etotal_f, traj.r_p_f[0], \
                        #        traj.r_p_f[1], traj.r_p_f[2], traj.polar_f, traj.azi_f, traj.time, traj.turn_pnts, traj.cl_appr, \
                        #        traj.r_p_min[0], traj.r_p_min[1], traj.r_p_min[2]))


        #outfile_15.close()
        #outfile_30.close()
        #outfile_45.close()
        #outfile_60.close()







###### SCRIPT ######

# open logfile
logfile = open(logfilename, 'w')

print("Screen output will be automatically written to {}!".format(logfilename))

print("Created by version %4.2f" % VERSION_ID)
logfile.write("Created by version %4.2f\n" % VERSION_ID)



# WARNING: the following is outdated since now empty or unfinished trajectory files are just skipped which is wanted!!
# check if unfinished trajectories are there and remove them
#remove_unfinished_traj(logfile)

traj_list = read_in_mxt_fins(logfile)
write_summary(logfile, outname, traj_list)

# in case of further analysis
write_statistics(traj_list)

logfile.close()
