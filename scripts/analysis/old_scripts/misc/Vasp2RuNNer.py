#!/usr/bin/python
# -*- coding: iso-8859-15 -*-
'''
Created on Feb 27, 2017

@author: mkammle
'''

from copy import copy, deepcopy
from sys import exit
from numpy import reshape, asarray
from sys import argv
from os import getcwd

RUNTYPE = None
 
    
class RuNNerInput:

    def __init__(self):
        self.latticeVectors = 3*[3*[None]]
        self.atomicCharge = 0.0
        self.atomicEnergy = 0.0
        self.atomicCoords = []
        self.atomicForces = []
        self.atomicSymbol = []
        self.energy = 0.0
        self.charge = 0.0

    def write(self, f):
        f.write("begin\n")
        for i in range(len(self.latticeVectors)):
            f.write("lattice %15.8f %15.8f %15.8f\n" % (self.latticeVectors[i][0],
                                                     self.latticeVectors[i][1],
                                                     self.latticeVectors[i][2]))
        for i in range(len(self.atomicCoords)):
            f.write("atom %15.8f %15.8f %15.8f %s %15.8f %15.8f %15.8f %15.8f %15.8f\n" % 
                        (self.atomicCoords[i][0], self.atomicCoords[i][1], self.atomicCoords[i][2],
                         self.atomicSymbol[i], self.atomicCharge, self.atomicEnergy,
                         self.atomicForces[i][0], self.atomicForces[i][1], self.atomicForces[i][2]))
        f.write("energy %15.8f\n" % self.energy)
        f.write("charge %15.8f\n" % self.charge)
        f.write("end\n")
        
        
class VaspOutput:
    
    massOf = {"Ag" : 107.868,
              "Al" : 26.981,
              "Au" : 196.966,
              "C"  : 12.011,
              "Cu" : 63.546,
              "H"  : 1.000,
              "Ir" : 192.220,
              "Ni" : 58.690,
              "Pd" : 106.420,
              "Pt" : 195.080,
              "Rh" : 102.906}
    
    ang2bohr = 1.889725989
    eV2Hartree = 0.0367493
    amu2kg = 1.660539040e-27
    angstrom2meter = 1e-10
    second2fs = 1e15
    joule2hartree = 2.293710449e17
    bohr2meter = 5.29177249e-11
    VASP2RuNNerForce = amu2kg * angstrom2meter * second2fs**2 * joule2hartree * bohr2meter
    
    def __init__(self):
        self.latticeVectors = 3*[3*[None]]
        self.prevAtomicCoords = []
        self.atomicCoords = []
        self.atomicCoordsWoPBCs = None
        self.nextAtomicCoords = []
        self.atomicSymbol = []
        self.atomicCount = []
        self.atomicForces = []
        self.timestep = 0.0
        self.numSteps = 0
        self.energy = []
        
        
    def __generateSymbols(self):
        
        # delete temporary save
        symbolline = copy(self.atomicSymbol)
        self.atomicSymbol = []
        
        for i in range(len(symbolline)):
            for _ in range(self.atomicCount[i]):
                self.atomicSymbol.append(symbolline[i])
                
    def __forcesFromCentralDifferences(self):
        
        if self.numSteps < 3:
            exit("Dynamic run must contain at least three entries. You have %d" % self.numSteps)
            
        for step in range(1, len(self.atomicCoords)-1):
            for atom in range(sum(self.atomicCount)):
                try:
                    mass = VaspOutput.massOf[self.atomicSymbol[atom]]
                except(KeyError):
                    exit("Mass of element %s is unknown. Please add to dictionary." % self.atomicSymbol[atom])
                
                self.atomicForces.append([ mass *
                                          (    self.atomicCoordsWoPBCs[step+1][atom][dim]
                                           - 2*self.atomicCoordsWoPBCs[step]  [atom][dim]
                                           +   self.atomicCoordsWoPBCs[step-1][atom][dim]) / self.timestep**2 
                                          for dim in range(3)])
        # this is now in amu*angstrom/fs**2, but we need it in hartree/bohr
        #
        #  amu * A     kg * m * fs²    J        J     Ha * m     Ha
        # --------- * ------------- = ---  ;;  --- * -------- = ----
        #    fs²       amu * A * s²    m        m     J * a0     a0
        self.atomicForces = [[VaspOutput.VASP2RuNNerForce * dim for dim in force] for force in self.atomicForces]
        
#        tmp = reshape(self.atomicForces, ((len(self.atomicCoords)-2)*25*3))
#        print max(abs(tmp))
        
        # delete first and last set of atomic coordinates since we don't have forces for them
        self.atomicCoords = self.atomicCoords[1:-1]
        self.numSteps -= 2


    def __removePBCs(self):
        self.atomicCoords = reshape(self.atomicCoords, (self.numSteps, sum(self.atomicCount), 3))
        self.atomicCoordsWoPBCs = deepcopy(self.atomicCoords)
        
        for ts in range(1, self.numSteps):
            for atom in range(sum(self.atomicCount)):
                for dim in range(3):
                    diff = self.atomicCoordsWoPBCs[ts][atom][dim] - self.atomicCoordsWoPBCs[ts-1][atom][dim]
                    if abs(diff) < 0.5:
                        continue
                    self.atomicCoordsWoPBCs[ts][atom][dim] -= round(diff)
        
        
        self.atomicCoords       = reshape(self.atomicCoords,       (self.numSteps*sum(self.atomicCount), 3))
        self.atomicCoordsWoPBCs = reshape(self.atomicCoordsWoPBCs, (self.numSteps*sum(self.atomicCount), 3))
        
    
       
                
                
    def process(self):
        # create array with atomic symbols
        self.__generateSymbols()
        
        # remove periodic boundaries
        self.__removePBCs()
        
        # convert to cartesian
        self.atomicCoords       = [asarray(self.latticeVectors).dot(asarray(coords)) for coords in self.atomicCoords]
        self.atomicCoordsWoPBCs = [asarray(self.latticeVectors).dot(asarray(coords)) for coords in self.atomicCoordsWoPBCs]
        
        # reshape atomic coordinates
        self.atomicCoords       = reshape(self.atomicCoords,       (self.numSteps, sum(self.atomicCount), 3))
        self.atomicCoordsWoPBCs = reshape(self.atomicCoordsWoPBCs, (self.numSteps, sum(self.atomicCount), 3))
        
        if (RUNTYPE == "dynamic"):
            self.__forcesFromCentralDifferences()
            self.atomicForces = reshape(self.atomicForces, (self.numSteps, sum(self.atomicCount), 3))
        elif (RUNTYPE == "static"):
            self.atomicForces = [[VaspOutput.eV2Hartree/VaspOutput.ang2bohr * dim for dim in force] for force in self.atomicForces]
        else:
            exit("Unknown runtype %s" % RUNTYPE)
            

        
    def read(self):
        
        lc = 1
        f = None
        
        ### XDATCAR ###
        # file might be called any of these names
        try:
            f = open("XDATCAR", "r")
        except(IOError):
            try:
                f = open("comb_XDATCAR", "r")
            except(IOError):
                exit("No XDATCAR file found in %s." % getcwd())
                      
            
        # read file XDATCAR
        for line in f:
            if lc < 3:
                pass
            elif 3 <= lc < 6:
                sl = line.strip(" \n\r\t").split()
                self.latticeVectors[lc-3] = [float(entry) for entry in sl]
            elif lc == 6:
                self.atomicSymbol = line.strip(" \n\r\t").split()  # will be adjusted later
            elif lc == 7:
                sl = line.strip(" \n\r\t").split()
                self.atomicCount = [int(entry) for entry in sl]
            elif lc > 7:
                if line.strip(" \n\r\t").startswith("Direct configuration="):
                    continue
                sl = line.strip(" \n\r\t").split()
                self.atomicCoords.append([float(entry) for entry in sl])
            lc += 1
        f.close()
        f = None
        
        
        ### OSZICAR ###
        # OSZICAR might be called any of these names
        try:
            f = open("OSZICAR", "r")
        except(IOError):
            try:
                f = open("comb_OSZICAR", "r")
            except(IOError):
                exit("No OSZICAR file found in %s." % getcwd())

            
            
        # read file OSZICAR
        for line in f:
            if "E0=" in line:
                sl = line.strip(" \n\r\t").split()
                e0TagPosition = sl.index("E0=")
                self.energy.append(float(sl[e0TagPosition+1]))
        f.close()
        f = None     
        
        
        # check if OSZICAR matches XDATCAR
        self.numSteps = len(self.atomicCoords)/sum(self.atomicCount)
        if self.numSteps != len(self.energy):
            exit("OSZICAR and XDATCAR do not match. I found %d energies for %d timesteps" % (len(self.energy), self.numSteps))
            
            
        if (RUNTYPE == "static"):
            readingForces = False
            # if static, read OUTCAR
            try:
                f = open("OUTCAR", "r")
            except(IOError):
                exit("No OSZICAR file found in %s. This file is required for a static calculation." % getcwd())
            
            # read file OUTCAR
            for line in f:
                if "POSITION" in line and "TOTAL-FORCE" in line:
                    readingForces = True
                elif line.strip(" \n\r\t").startswith("----------------------"):
                    continue
                elif line.strip(" \n\r\t").startswith("total drift"):
                    break
                elif readingForces:
                    sl = line.strip(" \n\r\t").split()
                    self.atomicForces.append([float(entry) for entry in sl[3:]])
            f.close()
            
        
        # if dynamic, try to extract the timestep
        if (RUNTYPE == "dynamic"):
            try:
                f = open("OUTCAR", "r")
            except(IOError):
                self.timestep = 0.1 # fs
                print "No OUTCAR found, assuming timestep is 0.1 fs"
            else:
                for line in f:
                    if "POTIM" in line:
                        sl = line.strip(" \n\r\t").split()
                        self.timestep = float(sl[2])
                        print "Found timestep of %f fs." % self.timestep
                        break
                f.close()
                
        

                    
                    

if __name__ == '__main__':
    
    if len(argv) < 2:
        exit("Need to specify -static or -dynamic as argument.")
    elif (argv[1] == "-static" or
          argv[1] == "-dynamic"):
        RUNTYPE = argv[1]
    else:
        exit("Unknown flag.")
    RUNTYPE = RUNTYPE[1:]   # delete the -
        
    vo = VaspOutput()
    vo.read()
    vo.process()
    
    outfile = open("ruNNer.input", "w")
    for i in range(len(vo.atomicCoords)):
        ri = RuNNerInput()
        ri.latticeVectors = [[VaspOutput.ang2bohr * dim for dim in latticeVector] for latticeVector in vo.latticeVectors]
        ri.atomicCoords = [[VaspOutput.ang2bohr * dim for dim in config] for config in vo.atomicCoords[i]]
        ri.atomicForces = vo.atomicForces[i]
        ri.atomicSymbol = vo.atomicSymbol
        ri.energy = VaspOutput.eV2Hartree * vo.energy[i]
        ri.write(outfile)
    outfile.close()
        

          
