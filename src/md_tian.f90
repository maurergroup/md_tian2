program md_tian
    ! Purpose:
    !       Do molecular dynamics, Langevin dynamics, Ring Polymer dynamics
    !
    ! Date          	Author          	    History of Revison
    ! ====          	======          	    ==================
    ! 30.03.2017    	Marvin Kammler		    new data structure
    !                   Sascha Kandratsenka     and propagation mathods
    !
    ! 18.02.2014    	Svenja M. Janke		    Original
    !			        Sascha Kandratsenka
    !			        Dan J. Auerbach
    !


use md_init
use atom_class
use pes_lj_module

implicit none

type(atoms) :: teil, slab


call simbox_init(teil, slab)
call compute_lj(teil, slab, energy_and_force)



end program md_tian
