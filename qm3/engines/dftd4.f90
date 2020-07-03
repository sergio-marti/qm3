!
! gfortran -shared|-dynamiclib -o libdftd4.so -Ilib dftd4.f90 libdftd4.a
!
module qm3
    use mctc_environment
    use class_molecule
    use class_set
    use class_param
    use class_results
    use dispersion_calculator
    use dftd4
    implicit none
    public
    real*8, parameter  :: AA__Bohr = 1.0d0 / 0.52917726d0
    type(molecule) :: mol
    type(dispersion_model), allocatable :: dispm
    type(dftd_parameter) :: dparam
    type(dftd_options),  parameter :: opt = dftd_options ( &
        &  lmbd = p_mbd_approx_atm, refq = p_refq_goedecker, &
        &  wf = 6.0d0, g_a = 3.0d0, g_c = 2.0d0, &
        &  lmolpol=.false., lenergy=.true., lgradient=.true., lhessian=.false., &
        &  print_level = 0)
end module qm3


subroutine qm3_dftd4_init( nat, siz, dat )
	use qm3
	implicit none
	integer, intent( in ) :: nat, siz
	real*8, dimension(0:siz-1), intent( inout ) :: dat
    integer :: i

    call mol%allocate( nat, .false. )
    mol%chrg = dat(0)
    do i = 1, nat
        mol%at(i) = dint( dat(i) )
    end do
    dparam = dftd_parameter( s6 = dat(nat+1), s8 = dat(nat+2), a1 = dat(nat+3), a2 = dat(nat+4) )
end subroutine qm3_dftd4_init


subroutine qm3_dftd4_calc( nat, siz, dat )
	use qm3
	implicit none
	integer, intent( in ) :: nat, siz
	real*8, dimension(0:siz-1), intent( inout ) :: dat
    type(mctc_logger) :: env
    type(dftd_results) :: dresul
	integer :: i, j, n

	do i = 1, nat
		j = 3 * ( i - 1 )
		mol%xyz(1,i) = dat(j  ) * AA__Bohr
		mol%xyz(2,i) = dat(j+1) * AA__Bohr
		mol%xyz(3,i) = dat(j+2) * AA__Bohr
	end do

    call d4_calculation( 999, env, opt, mol, dparam, dresul )

	dat(0) = dresul%energy
	do i = 1, nat
		j = 3 * ( i - 1 )
		dat(j+1) = dresul%gradient(1,i)
		dat(j+2) = dresul%gradient(2,i)
		dat(j+3) = dresul%gradient(3,i)
	end do

    call dresul%deallocate

end subroutine qm3_dftd4_calc
