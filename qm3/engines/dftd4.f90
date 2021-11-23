!
! gfortran -shared|-dynamiclib -o libdftd4.so -Ilib dftd4.f90 libdftd4.a
!
module qm3
    use dftd4, only: structure_type, new, d4_model, new_d4_model, rational_damping_param, get_dispersion, realspace_cutoff
    implicit none
    public
    real*8, parameter  :: AA__Bohr = 1.0d0 / 0.52917726d0
    type( structure_type ) :: mol
    type( d4_model ) :: d4
    type( rational_damping_param ) :: param
end module qm3


subroutine qm3_dftd4_init( nat, siz, dat )
	use qm3
	implicit none
	integer, intent( in ) :: nat, siz
	real*8, dimension(0:siz-1), intent( inout ) :: dat
    real*8, dimension(:,:), allocatable :: xyz
    integer, dimension(:), allocatable :: atn
    integer :: i

    allocate( atn(1:nat), xyz(1:3,1:nat) )
    do i = 1, nat
        atn(i) = dint( dat(i) )
    end do
    xyz = 0.0d0
    call new( mol, atn, xyz )
    mol%charge = dat(0)
    call new_d4_model( d4, mol, ga = 3.0d0, gc = 2.0d0, wf = 6.0d0 )
    param = rational_damping_param( s6 = dat(nat+1), s8 = dat(nat+2), s9 = 1.0d0, a1 = dat(nat+3), a2 = dat(nat+4), alp = 16.0d0 )
    deallocate( atn, xyz )
end subroutine qm3_dftd4_init


subroutine qm3_dftd4_calc( nat, siz, dat )
	use qm3
	implicit none
	integer, intent( in ) :: nat, siz
	real*8, dimension(0:siz-1), intent( inout ) :: dat
    real*8, dimension(:,:), allocatable :: grad, sigm
    real*8 :: ener
	integer :: i, j, n

	do i = 1, nat
		j = 3 * ( i - 1 )
		mol%xyz(1,i) = dat(j  ) * AA__Bohr
		mol%xyz(2,i) = dat(j+1) * AA__Bohr
		mol%xyz(3,i) = dat(j+2) * AA__Bohr
	end do

    allocate( grad(1:3,1:nat), sigm(1:3,1:3) )
    call get_dispersion( mol, d4, param, realspace_cutoff(), ener, grad, sigm )

	dat(0) = ener
	do i = 1, nat
		j = 3 * ( i - 1 )
		dat(j+1) = grad(1,i)
		dat(j+2) = grad(2,i)
		dat(j+3) = grad(3,i)
	end do

    deallocate( grad, sigm )

end subroutine qm3_dftd4_calc
