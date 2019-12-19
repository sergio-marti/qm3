!
! gfortran -shared|-dynamiclib -o libdftd3.so -Ilib dftd3.f90 libdftd3.a
!
module qm3
    use dftd3_api
    implicit none
    public
    real*8, parameter  :: AA__Bohr = 1.0d0 / 0.52917726d0
    type( dftd3_calc ) :: calc
end module qm3


! s6,rs6,s18,rs18,alp
subroutine qm3_dftd3_init( prm )
	use qm3
	implicit none
    real*8, dimension(0:5), intent( in ) ::prm 
    type( dftd3_input ) :: input
    real*8, dimension(1:5) :: xprm
    integer :: xver

    call dftd3_init( calc, input )
!    call dftd3_set_functional( calc, func = trim( fun ), version = ver, tz = ( qtz == 1 ) )
    xver = dint( prm(0) )
    xprm(1:5) = prm(1:5)
    call dftd3_set_params( calc, xprm, xver )
end subroutine qm3_dftd3_init


subroutine qm3_dftd3_calc( nat, siz, dat )
	use qm3
	implicit none
	integer, intent( in ) :: nat, siz
	real*8, dimension(0:siz-1), intent( inout ) :: dat

	integer :: i, j, n
! - allocatable?
    integer :: anum(nat)
	real*8  :: ener
	real*8  :: coor(3,nat), grad(3,nat)

	do i = 1, nat
        anum(i) = dint( dat(i-1) )
    end do
	do i = 1, nat
		j = nat + 3 * ( i - 1 )
		coor(1,i) = dat(j  ) * AA__Bohr
		coor(2,i) = dat(j+1) * AA__Bohr
		coor(3,i) = dat(j+2) * AA__Bohr
	end do

    call dftd3_dispersion( calc, coor, anum, ener, grad )

	dat(0) = ener 
	do i = 1, nat
		j = 3 * ( i - 1 )
		dat(j+1) = grad(1,i)
		dat(j+2) = grad(2,i)
		dat(j+3) = grad(3,i)
	end do

end subroutine qm3_dftd3_calc
