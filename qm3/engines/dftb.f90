!
! gfortran -shared|-dynamiclib -o dftb.so -I{DFTBPLUS_git}/_build/api/mm dftb.f90 libdftb+.a libxmlf90.a libdftd3.a -framework Accelerate -fopenmp
!
module qm3
	use dftbplus
	use dftbp_constants, only : AA__Bohr
	implicit none
	public
	type( TDftbPlus ) :: dftbp
	integer :: devnull
end module qm3



subroutine qm3_dftb_init
	use qm3
	implicit none
	type( TDftbPlusInput ) :: input

	open( newunit = devnull, file = "dftb_in.log", action = "write" )
!	open( newunit = devnull, file = "/dev/null", action = "write" )
	call TDftbPlus_init( dftbp, outputUnit = devnull )
	call dftbp%getInputFromFile( "dftb_in.hsd", input )
	call dftbp%setupCalculator( input )

end subroutine qm3_dftb_init



subroutine qm3_dftb_calc( nQM, nMM, siz, dat )
	use qm3
	implicit none
	integer, intent( in ) :: nQM, nMM, siz
	real*8, dimension(0:siz-1), intent( inout ) :: dat

	integer :: i, j
! - allocatable?
	real*8 :: ener, qm_chg(nQM), mm_chg(nMM)
	real*8 :: qm_crd(3,nQM), mm_crd(3,nMM), qm_grd(3,nQM), mm_grd(3,nMM)

	do i = 1, nQM
		j = 3 * ( i - 1 )
		qm_crd(1,i) = dat(j)   * AA__Bohr
		qm_crd(2,i) = dat(j+1) * AA__Bohr
		qm_crd(3,i) = dat(j+2) * AA__Bohr
	end do
	do i = 1, nMM
		j = 3 * ( nQM + i - 1 )
		mm_crd(1,i) = dat(j)   * AA__Bohr
		mm_crd(2,i) = dat(j+1) * AA__Bohr
		mm_crd(3,i) = dat(j+2) * AA__Bohr
	end do
	do i = 1, nMM
		j = 3 * ( nQM + nMM ) + i - 1
		mm_chg(i) = dat(j)
	end do

	call dftbp%setGeometry( qm_crd )
	call dftbp%setExternalCharges( mm_crd, mm_chg )
	call dftbp%getEnergy( ener )
	call dftbp%getGrossCharges( qm_chg )
	call dftbp%getGradients( qm_grd )
	call dftbp%getExtChargeGradients( mm_grd )

	dat(0) = ener 
	do i = 1, nQM
		dat(i) = qm_chg( i )
	end do
	do i = 1, nQM
		j = 1 + nQM + ( i - 1 ) * 3
		dat(j)   = qm_grd(1,i)
		dat(j+1) = qm_grd(2,i)
		dat(j+2) = qm_grd(3,i)
	end do
	do i = 1, nMM 
		j = 1 + 4 * nQM + ( i - 1 ) * 3
		dat(j)   = mm_grd(1,i)
		dat(j+1) = mm_grd(2,i)
		dat(j+2) = mm_grd(3,i)
	end do

end subroutine qm3_dftb_calc



subroutine qm3_dftb_stop
	use qm3
	implicit none

	call TDftbPlus_destruct( dftbp )
	close( devnull )

end subroutine qm3_dftb_stop
