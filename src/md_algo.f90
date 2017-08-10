module md_algo

    use universe_mod

    implicit none



contains

    subroutine propagate_1(atoms)

        type(universe), intent(inout) :: atoms


    end subroutine propagate_1


    subroutine propagate_2(atoms)
        type(universe), intent(inout) :: atoms


    end subroutine propagate_2

    subroutine ring_polymer_step(atoms)
        type(universe), intent(inout) :: atoms

    end subroutine ring_polymer_step

end module md_algo
