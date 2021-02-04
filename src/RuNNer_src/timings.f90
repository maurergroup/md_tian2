!######################################################################
! This routine is part of
! RuNNer - RuNNer Neural Network Energy Representation
! (c) 2008-2020 Prof. Dr. Joerg Behler 
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
      module timings 

      implicit none

!! for main.f90
      integer runtimeday
      real*8 runtimestart
      real*8 runtimeend

      integer dayinitnn                                
      real*8 timeinitnnstart 
      real*8 timeinitnnend   
      real*8 timeinitnn

      integer dayreadinput                                
      real*8 timereadinputstart 
      real*8 timereadinputend   
      real*8 timereadinput

      integer daymode1                                
      real*8 timemode1start 
      real*8 timemode1end   
      real*8 timemode1

      integer daymode2                                
      real*8 timemode2start 
      real*8 timemode2end   
      real*8 timemode2

      integer daymode3                                
      real*8 timemode3start 
      real*8 timemode3end   
      real*8 timemode3

      integer dayfinalize                                
      real*8 timefinalizestart 
      real*8 timefinalizeend   
      real*8 timefinalize

!! times still to be sorted

      integer dayoutput                               
      real*8 timeoutputstart                                            
      real*8 timeoutputend                                              
      real*8 timeoutput                                                 

!! timings for detailed_timing_epoch for mode 1

!! timings for detailed_timing_epoch for mode 2

      integer dayepochini                               
      real*8 timeepochinistart              
      real*8 timeepochiniend                
      real*8 timeepochini                          

      integer dayepoch                               
      real*8 timeepochstart                                            
      real*8 timeepochend                                              
      real*8 timeepoch                                                 

      integer daygeterror  
      real*8 timegeterrorstart
      real*8 timegeterrorend    
      real*8 timegeterror     
   
      integer daycalconefunction  
      real*8 timecalconefunctionstart
      real*8 timecalconefunctionend    
      real*8 timecalconefunction     
   
      integer daygetshortforces  
      real*8 timegetshortforcesstart
      real*8 timegetshortforcesend    
      real*8 timegetshortforces     
   
      integer daymix  
      real*8 timemixstart  
      real*8 timemixend      
      real*8 timemix

      integer dayio                                   
      real*8 timeiostart                                                
      real*8 timeioend                                                  
      real*8 timeio                                                     

      integer dayefitting
      real*8 timeefittingstart 
      real*8 timeefittingend  
      real*8 timeefitting    
 
      integer dayqfitting
      real*8 timeqfittingstart 
      real*8 timeqfittingend  
      real*8 timeqfitting    
 
      integer dayefittingrepeat
      real*8 timeefittingrepeatstart 
      real*8 timeefittingrepeatend  
      real*8 timeefittingrepeat    
 
      integer dayffitting
      real*8 timeffittingstart  
      real*8 timeffittingend   
      real*8 timeffitting     

      integer dayshortfit
      real*8 timeshortfitstart 
      real*8 timeshortfitend  
      real*8 timeshortfit    
 
      integer dayelecfit
      real*8 timeelecfitstart 
      real*8 timeelecfitend  
      real*8 timeelecfit    
 
      integer daydeshortdw                                
      real*8 timedeshortdwstart 
      real*8 timedeshortdwend   
      real*8 timedeshortdw

      integer daydeshortdwrepeat                                
      real*8 timedeshortdwrepeatstart 
      real*8 timedeshortdwrepeatend   
      real*8 timedeshortdwrepeat

      integer daydfshortdw 
      real*8 timedfshortdwstart
      real*8 timedfshortdwend 
      real*8 timedfshortdw  

      integer dayeupdate 
      real*8 timeeupdatestart
      real*8 timeeupdateend   
      real*8 timeeupdate 
    
      integer dayeupdaterepeat 
      real*8 timeeupdaterepeatstart
      real*8 timeeupdaterepeatend   
      real*8 timeeupdaterepeat 
    
      integer dayfupdate  
      real*8 timefupdatestart 
      real*8 timefupdateend     
      real*8 timefupdate           

      integer dayferror 
      real*8 timeferrorstart 
      real*8 timeferrorend  
      real*8 timeferror   

      integer dayeerror  
      real*8 timeeerrorstart 
      real*8 timeeerrorend    
      real*8 timeeerror     
   
      integer dayeerrorrepeat  
      real*8 timeeerrorrepeatstart 
      real*8 timeeerrorrepeatend    
      real*8 timeeerrorrepeat     
   
      integer dayfsym  
      real*8 timefsymstart  
      real*8 timefsymend      
      real*8 timefsym

      integer dayqerror  
      real*8 timeqerrorstart 
      real*8 timeqerrorend    
      real*8 timeqerror     
   
      integer daydqdw  
      real*8 timedqdwstart 
      real*8 timedqdwend    
      real*8 timedqdw     
   
      integer dayqupdate  
      real*8 timequpdatestart 
      real*8 timequpdateend    
      real*8 timequpdate     
   
!! timings for detailed_timing_epoch for mode 3

      integer dayshort                                
      real*8 timeshortstart                                             
      real*8 timeshortend                                               
      real*8 timeshort                                                  

      integer dayallocshort                             
      real*8 timeallocshortstart                       
      real*8 timeallocshortend                        
      real*8 timeallocshort                                 

      integer daysymshort                             
      real*8 timesymshortstart                                          
      real*8 timesymshortend                                            
      real*8 timesymshort                                               

      integer dayextrapolationshort                             
      real*8 timeextrapolationshortstart                   
      real*8 timeextrapolationshortend             
      real*8 timeextrapolationshort                                   

      integer dayextrapolationewald                             
      real*8 timeextrapolationewaldstart                   
      real*8 timeextrapolationewaldend             
      real*8 timeextrapolationewald                                   

      integer dayscalesymshort                             
      real*8 timescalesymshortstart                             
      real*8 timescalesymshortend                              
      real*8 timescalesymshort                                     

      integer dayscalesymewald                             
      real*8 timescalesymewaldstart                             
      real*8 timescalesymewaldend                              
      real*8 timescalesymewald                                     

      integer dayscaledsfuncshort                             
      real*8 timescaledsfuncshortstart                             
      real*8 timescaledsfuncshortend                              
      real*8 timescaledsfuncshort                                     

      integer dayeshort                               
      real*8 timeeshortstart                                            
      real*8 timeeshortend                                              
      real*8 timeeshort                                                 

      integer dayfshort                               
      real*8 timefshortstart                                            
      real*8 timefshortend                                              
      real*8 timefshort                                                 

      integer daysshort                               
      real*8 timesshortstart                                            
      real*8 timesshortend                                              
      real*8 timesshort                                                 

      integer daycharge                               
      real*8 timechargestart                                            
      real*8 timechargeend                                              
      real*8 timecharge                                                 

      integer daycomm1                                
      real*8 timecomm1start                                             
      real*8 timecomm1end                                               
      real*8 timecomm1                                                  

      integer dayelec                                 
      real*8 timeelecstart                                              
      real*8 timeelecend                                                
      real*8 timeelec                                                   

      integer daysymelec1                             
      real*8 timesymelec1start                                          
      real*8 timesymelec1end                                            
      real*8 timesymelec1                                               

      integer daysymelec2                             
      real*8 timesymelec2start                                          
      real*8 timesymelec2end                                            
      real*8 timesymelec2                                               

      integer dayeelec                             
      real*8 timeeelecstart                                    
      real*8 timeeelecend                                
      real*8 timeeelec                                   

! variables for subroutine date_and_time
      character*8 fulldate
      character*10 fulltime
      character*5 zone
      integer*4 timevalues(8)

      ! Emir (testing)
      integer daytest
      real*8 teststart
      real*8 testend

      end module timings 

