!######################################################################
! This routine is part of
! RuNNer - RuNNer Neural Network Energy Representation
! (c) 2008-2019 Prof. Dr. Joerg Behler 
! Georg-August-Universitaet Goettingen, Germany
!
! This program is free software: you can redistribute it and/or modify it 
! under the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or 
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
! or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
! for more details.
!
! You should have received a copy of the GNU General Public License along 
! with this program. If not, see http://www.gnu.org/licenses. 
!######################################################################
      module globaloptions 

      implicit none

!! cutoff_type = functional form of cutoff function (cutoff_type=1: f_c= 0.5(cos(pi*R_ij/R_c)+1)  , cutoff_type=1: f_c=(tanh(1.d0-rij/funccutoff_local(i2,iindex)))**3 )
      integer cutoff_type
      real*8    cutoff_alpha ! modified by kenko
!! nran = integer switch for selecting the type of random number generator
      integer nran

!! nblock = number of structure in memory at the same time (modes 1 + 2) or number of atoms in memory at the same time (mode 3)
      integer nblock

!! max_num_atoms: number of atoms in the largest structure of all structures in input.data
      integer max_num_atoms 
      integer max_num_pairs

!! elementindex: yields number of element in system when nuclear charge of atom is provided. Example: for pure ZnO we have elementindex(6)=1 and elementindex(30)=2, all others are 0. 
      integer elementindex(102)

!! pairindex
      integer pairindex(102,102)
!! maxnum_weightsshort:
      integer maxnum_weights_short_atomic
      integer maxnum_weights_elec
      integer maxnum_weights_short_pair

!! maxnum_layers_short_atomic = number of short range hidden layers + 1 (atomic NN)
!! determined in main/initialization/getdimensions
      integer maxnum_layers_short_atomic

!! maxnum_layersewald = number of electrostatic NN hidden layers + 1 (atomic and pair NN)
!! determined in main/initialization/getdimensions
      integer maxnum_layers_elec

!! maxnum_layerspair = number of short range hidden layers + 1 (pair NN)
!! determined in main/initialization/getdimensions
      integer maxnum_layers_short_pair

      integer maxnum_funcvalues_short_atomic
      integer maxnum_funcvalues_elec
      integer maxnum_funcvalues_short_pair

      integer paramode
      integer listdim
      integer ewaldkmax

!! determined in getdimensions
      integer npairs
!! enforced value for max_num_neighbors
      integer max_num_neighbors_atomic_input

!! determined in getdimensions
      integer nelem ! number of elements in input.nn

!! determined in main/initialization/structurecount
      integer totnum_structures ! total number of structures in input.data file

      integer, dimension(:)  , allocatable :: nucelem
      integer, dimension(:,:), allocatable :: elempair

      real*8 ewaldalpha
      real*8 ewaldcutoff
      real*8 rscreen_cut
      real*8 rscreen_onset
      real*8, dimension(:) , allocatable :: fixedcharge
      real*8 rmin
      real*8, dimension(:)     , allocatable :: atomrefenergies
      real*8, dimension(:)     , allocatable :: dmin_element
      real*8 kalman_epsilon 
!! Noise matrix  and modified by kenko      
      real*8 kalman_q
      real*8 kalman_qtau
      real*8 kalman_qmin
      real*8 kalman_q0

      logical lscalesym
      logical lcentersym
      logical lremoveatomenergies
      logical lnormnodes
      logical lscreen
      logical lreadunformatted
      logical lwriteunformatted
      logical lfinetime
      logical lfinetimeepoch
      logical lsilent
      logical ldebug
      logical lompmkl
      logical lcheckf
      logical luseatomcharges
      logical luseatomenergies
      logical luseforces
      logical lenforcemaxnumneighborsatomic
      logical lmd
      logical lusenoisematrix
!! calculate correlation of symmetry functions
      logical lpearson_correlation

!! calculate the sensitivity
      logical lsens

      logical ldostress
      logical ldohessian

      logical ldetect_saturation
      real*8 saturation_threshold

!! give statistics of atomic environments
      logical lenvironmentanalysis

      character*2, dimension(:), allocatable :: element
      character*20 pstring

      end module globaloptions 

