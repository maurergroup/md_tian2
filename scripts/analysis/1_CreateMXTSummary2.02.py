#!/usr/bin/python3

# original version 2.00 and previous versions were done by Marvin Kammler
# version 2.01 was modified by Sebastian Wille
# version 2.02 was modified by Sebastian Wille

# Modifications v2.01:
# trajid added in Traj class
# logfile will be written
# set readfile name
# depending on final polar angle additional data files are written

# Modifications v2.02:
# initial and final velocities will also be given
# added time at closest approach
# added function to select traj according to settings from ion imaging experiment

# intention: analyze all traj/mxt_fin files and create the MXt2Summary file

# use like: python3 <scriptname> or ./<scriptname> # check first command line


# 2DO:
# merge this script with the analysis script and in the end use matplotlib to generate all plots using a single script



import os, sys, glob 

# set names for output and log file
outname     = "MXT2Summary"
logfilename = "CreateMXTSummary.log"
settingname = "plot_settings.dat"

METAL_TYPE = "C"
SHOT_THRU_LIMIT = 0.0

SPECULAR_RADIUS = 1.5 # should match experimental settings (RAT detector radius)
ION_IMAGING_AZI = 5.0   # should match experimental settings (ion imaging detector settings)
AZIMUTHAL_ANGLE = 42.0  # should match experimental settings (ion imaging surface shift according to LEED)

ANGLE_MAX = 90  # maximum angle in degrees
ANGLE_MIN = -90 # minimum angle in degrees

# add range of scattered angle to look at

############# NO CHANGES BELOW THIS LINE ######################################

### WRITE SETTINGS TO FILE ###
settingfile = open(settingname, 'w')

settingfile.write("METAL_TYPE {}\n".format(METAL_TYPE))
settingfile.write("SHOT_THRU_LIMIT {}\n".format(SHOT_THRU_LIMIT))
settingfile.write("SPECULAR_RADIUS {}\n".format(SPECULAR_RADIUS))
settingfile.write("ION_IMAGING_AZI {}\n".format(ION_IMAGING_AZI))
settingfile.write("AZIMUTHAL_ANGLE {}\n".format(AZIMUTHAL_ANGLE))
settingfile.write("ANGLE_MAX {}\n".format(ANGLE_MAX))
settingfile.write("ANGLE_MIN {}\n".format(ANGLE_MIN))

settingfile.close()


VERSION_ID = 2.02

class Traj:
        def __init__(self, fname, ekin_p_i, ekin_l_i, epot_i, etotal_i,\
                r_p_i, v_p_i, polar_i, azi_i, ekin_p_f, ekin_l_f, epot_f,   \
                        etotal_f, r_p_f, v_p_f, polar_f, azi_f, time, turn_pnts, \
                        cl_appr, cl_appr_t, r_p_min, traj_id):
                self.fname     = fname
                self.traj_id   = traj_id
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

        def containsNaN(self):
	        return "NaN" in self.ekin_p_i or "NaN" in self.ekin_l_i or "NaN" in self.epot_i or "NaN" in self.etotal_i or \
		        "NaN" in self.ekin_p_f or "NaN" in self.ekin_l_f or "NaN" in self.epot_f or "NaN" in self.etotal_f

        def __str__(self):
                return "   fname " + str(self.fname) +\
                        "\ntraj_id " + str(self.traj_id) +\
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
			"\nr_p_min " + str(self.r_p_min)



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
                        print("{}%".format(100*counter/num_folders+1))
                        logfile.write("{}%\n".format(100*counter/num_folders+1))

                infile = open(folder, 'r')
                traj_id = folder[7:15]
                checkline = file_len(folder) # check if number of entries match finished traj file

                if checkline == 23: # only finished traj files have all needed entries, this will skip unfinished trajs
                        for line in infile:
                                sline = line.split()
                                #try:
                                if "ekin_p_i" in sline:
                                    ekin_p_i = sline[-1]	#	ekin_p_i =      1.9200000
                                elif "ekin_l_i" in sline:
                                    ekin_l_i = sline[-1]	#	ekin_l_i =      0.9553543
                                elif "epot_i" in line:
                                    epot_i  = sline[-1]         #       epot_i   =   -392.2094094
                                elif "etotal_i" in line:
                                    etotal_i = sline[-1]        #       etotal_i =   -369.3305490
                                if line.startswith("r_i"):
                                    r_p_i = sline[-3:]          #       r_i      =     30.1309358    32.4013593     3.5000000
                                if "v_i" in line:
                                    v_p_i = sline[-3:]          #       v_i      =     0.0000000     0.0000000     -0.1934888 
                                if "polar_i" in line:
                                    polar_i = sline[-1]         #       polar_i  =     50.0000000
                                if "azi_i" in line:
                                    azi_i = sline[-1]           #       azi_i    =      0.0000000
                                if "ekin_p_f" in line:
                                    ekin_p_f = sline[-1]        #       ekin_p_f =      1.5942738
                                if "ekin_l_f" in line:
                                    ekin_l_f = sline[-1]        #       ekin_l_f =      0.9610213
                                if "epot_f" in line:
                                    epot_f = sline[-1]          #       epot_f   =   -391.8832698
                                if "etotal_f" in line:
                                    etotal_f = sline[-1]        #       etotal_f =   -369.3305177
                                if line.startswith("r_f"):
                                    r_p_f = sline[-3:]          #       r_f      =     33.1131979    32.3776858     1.3585263
                                if "v_f" in line:
                                    v_p_f = sline[-3:]          #       v_f      =      0.0026410     0.0117354     0.1513796
                                if "polar_f" in line:
                                    polar_f = sline[-1]         #       polar_f  =    122.6747121
                                if "azi_f" in line:
                                    azi_f = sline[-1]           #       azi_f    =     -5.5773256
                                if "time" in line:
                                    time = sline[-1]            #       time     =     20.0000000
                                if "turn_pnts" in line:
                                    turn_pnts = sline[-1]       #       turn_pnts =           1
                                if "cl_appr " in line:
                                    cl_appr = sline[-1]         #       cl_appr  =      0.9694690
                                if "cl_appr_time" in line:
                                    cl_appr_t = sline[-1]       #       cl_appr_time  =          278
                                if "r_min_p" in line:
                                    r_p_min = sline[-3:]        #       r_min_p  =     51.3246421    33.1600141     0.9996723

                        traj = Traj(folder, ekin_p_i, ekin_l_i, epot_i, etotal_i,  \
                                    r_p_i, v_p_i, polar_i, azi_i, ekin_p_f, ekin_l_f, epot_f,\
                                    etotal_f, r_p_f, v_p_f, polar_f, azi_f, time, turn_pnts,\
                                    cl_appr, cl_appr_t, r_p_min, traj_id)
                                #except:
                                        #print("Error in file {} in line {}".format(folder,line))
                                        #logfile.write("Error in file {} in line {}\n".format(folder,line))
                                        #sys.exit()	

		        #infile.close()
		        #traj_list[counter] = traj
                        traj_list.append(traj)
                        counter += 1

                else:
                        print("Skipping unfinished traj {}".format(folder))
                        logfile.write("Skipping unfinished traj {}\n".format(folder))
                        #break
                        continue

                infile.close()

        
        os.chdir("../")
        return traj_list


def write_summary(logfile, outfile_name_tmp, traj_list):
        outfile_name = outfile_name_tmp + ".txt" 
        outfile = open(outfile_name, "w")
        outfile.write("# traj_id E_kin_p   E_kin_l        E_pot      E_total r_p(    x,        y,         z) v_p(    x,        y,         z)      polar       azi   E_kin_p  E_kin_l       E_pot       E_total     r_p(    x,       y,         z) v_p(    x,       y,         z)        polar       azi     simtime turn_pnts   cl_appr   cl_appr_t   r_p_min\n")
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
                        print("ERROR: There's something wrong with file {}\n".format(str(traj)))
                        logfile.write("ERROR: There's something wrong with file {}\n".format(str(traj)))
                        sys.exit()

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

def write_traj_to_file(this_traj,this_outfile):
        this_outfile.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" \
                       % ( this_traj.traj_id, this_traj.ekin_p_i, this_traj.ekin_l_i, this_traj.epot_i, 
                           this_traj.etotal_i, this_traj.r_p_i[0], this_traj.r_p_i[1], this_traj.r_p_i[2], \
                           this_traj.v_p_i[0], this_traj.v_p_i[1], this_traj.v_p_i[2], this_traj.polar_i, \
                           this_traj.azi_i, this_traj.ekin_p_f, this_traj.ekin_l_f, this_traj.epot_f, \
                           this_traj.etotal_f, this_traj.r_p_f[0], this_traj.r_p_f[1], this_traj.r_p_f[2], \
                           this_traj.v_p_f[0], this_traj.v_p_f[1], this_traj.v_p_f[2], this_traj.polar_f, \
                           this_traj.azi_f, this_traj.time, this_traj.turn_pnts, this_traj.cl_appr, this_traj.cl_appr_t, \
                           this_traj.r_p_min[0], this_traj.r_p_min[1], this_traj.r_p_min[2]))


def write_angle_files(traj_list):
        # make the following more general
        outfilestr_15 = "MXT2Summary_15.txt"
        outfilestr_30 = "MXT2Summary_30.txt"
        outfilestr_45 = "MXT2Summary_45.txt"
        outfilestr_60 = "MXT2Summary_60.txt"

        outfile_15 = open(outfilestr_15, "w")
        outfile_30 = open(outfilestr_30, "w")
        outfile_45 = open(outfilestr_45, "w")
        outfile_60 = open(outfilestr_60, "w")

        for traj in traj_list:
                if 7.5 <= float(traj.polar_f) <= 22.5:
                        write_traj_to_file(traj,outfile_15)
                if 22.5 <= float(traj.polar_f) <= 37.5:
                        write_traj_to_file(traj,outfile_30)
                if 37.5 <= float(traj.polar_f) <= 52.5:
                        write_traj_to_file(traj,outfile_45)
                if 52.5 <= float(traj.polar_f) <= 67.5:
                        write_traj_to_file(traj,outfile_60)


        outfile_15.close()
        outfile_30.close()
        outfile_45.close()
        outfile_60.close()


def traj_in_ion_imaging(traj_list):

        foldername = "ion_imaging"

        if not os.path.exists(foldername):
            os.makedirs(foldername)

        ion_imaging_filename  =  foldername + "/" + outname + ".txt"
        ion_imaging_filenamet =  foldername + "/" + outname + "_test.txt"
        
        ion_imaging_file     = open(ion_imaging_filename, 'w')
        ion_imaging_filet      = open(ion_imaging_filenamet, 'w')
        
        ion_imaging_file.write("# traj_id E_kin_p   E_kin_l        E_pot      E_total r_p(    x,        y,         z) v_p(    x,        y,         z)      polar       azi   E_kin_p  E_kin_l       E_pot       E_total     r_p(    x,       y,         z) v_p(    x,       y,         z)        polar       azi     simtime turn_pnts   cl_appr   cl_appr_t   r_p_min\n")
        ion_imaging_filet.write("# traj_id E_kin_p   E_kin_l        E_pot      E_total r_p(    x,        y,         z) v_p(    x,        y,         z)      polar       azi   E_kin_p  E_kin_l       E_pot       E_total     r_p(    x,       y,         z) v_p(    x,       y,         z)        polar       azi     simtime turn_pnts   cl_appr   cl_appr_t   r_p_min\n")

        for traj in traj_list:
            if float(traj.polar_f) <= 31:
                if float(traj.azi_f) < 0:
                    inv_azi = AZIMUTHAL_ANGLE - 180
                    if inv_azi - ION_IMAGING_AZI <= abs(inv_azi - float(traj.azi_f)) <= inv_azi + ION_IMAGING_AZI:
                        write_traj_to_file(traj,ion_imaging_filet)
                    elif AZIMUTHAL_ANGLE - ION_IMAGING_AZI <= abs(AZIMUTHAL_ANGLE - float(traj.azi_f)) <= AZIMUTHAL_ANGLE + ION_IMAGING_AZI:
                        write_traj_to_file(traj,ion_imaging_filet)
                if float(traj.azi_f) > 0:
                    inv_azi = AZIMUTHAL_ANGLE + 180
                    if inv_azi - ION_IMAGING_AZI <= abs(inv_azi - float(traj.azi_f)) <= inv_azi + ION_IMAGING_AZI:
                        write_traj_to_file(traj,ion_imaging_filet)
                    elif AZIMUTHAL_ANGLE - ION_IMAGING_AZI <= abs(AZIMUTHAL_ANGLE - float(traj.azi_f)) <= AZIMUTHAL_ANGLE + ION_IMAGING_AZI:
                        write_traj_to_file(traj,ion_imaging_filet)

        
        for traj in traj_list:
            if float(traj.polar_f) <= 31:
                if float(traj.azi_f) < 0:
                    inv_azi = AZIMUTHAL_ANGLE - 180
                    if abs(inv_azi - float(traj.azi_f)) <= ION_IMAGING_AZI:
                        write_traj_to_file(traj,ion_imaging_file)
                    elif abs(AZIMUTHAL_ANGLE - float(traj.azi_f)) <= ION_IMAGING_AZI:
                        write_traj_to_file(traj,ion_imaging_file)
                if float(traj.azi_f) > 0:
                    inv_azi = AZIMUTHAL_ANGLE + 180
                    if abs(inv_azi - float(traj.azi_f)) <= ION_IMAGING_AZI:
                        write_traj_to_file(traj,ion_imaging_file)
                    elif abs(AZIMUTHAL_ANGLE - float(traj.azi_f)) <= ION_IMAGING_AZI:
                        write_traj_to_file(traj,ion_imaging_file)



        ion_imaging_file.close()
        ion_imaging_filet.close()





###### SCRIPT ######

# open logfile
logfile = open(logfilename, 'w')

print("Screen output will be automatically written to {}!".format(logfilename))

print("Created by version %4.2f" % VERSION_ID)
logfile.write("Created by version %4.2f\n" % VERSION_ID)



# WARNING: the following function is outdated since now empty or unfinished trajectory files are just skipped which is wanted!!
# check if unfinished trajectories are there and remove them
#remove_unfinished_traj(logfile)

traj_list = read_in_mxt_fins(logfile)
write_summary(logfile, outname, traj_list)

# Ion Imaging Experiment
traj_in_ion_imaging(traj_list)


# H@Gr related functions
#write_angle_files(traj_list)
#traj_in_es(traj_list)


logfile.close()
