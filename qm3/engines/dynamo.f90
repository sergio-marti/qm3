!
! compile to a library (.so) by adding "-shared" on linux, or "-dynamiclib" on macOS
!
subroutine qm3_initialize
	use dynamo
	implicit none
	logical, dimension(:), allocatable :: flg

	open( unit=output, file="dynamo.log", status="replace", access="stream", form="formatted" )
	call dynamo_header

	call mm_file_process( "borra", "opls" )
	call mm_system_construct( "borra", "seq" )
	call coordinates_read( "crd" )

	open( unit=999, file="mm_data", action="write", form="unformatted" )
	write(999) atmmas
	write(999) atmchg
	write(999) atmeps * 0.5d0 ! >> [geom] sqrt( kJ/mol )
	write(999) ( atmsig ** 2 ) * 0.5612310241546865d0 ! >> [arith] rmin/2 (A)
	close(999)

	allocate( flg(1:natoms) )
	flg = atom_selection( subsystem = (/ "ACS" /) )
	call mopac_setup( method = "AM1", charge = 0, selection = flg )
	call energy_initialize
	call energy_non_bonding_options( &
		list_cutoff   = 15.5_dp, &
		outer_cutoff  = 13.5_dp, &
		inner_cutoff  = 11.5_dp, &
		minimum_image = .true. )
	deallocate( flg )
end subroutine qm3_initialize



subroutine qm3_update_coor( coor )
	use dynamo
	implicit none
	real*8, dimension(1:3,1:natoms), intent( in ) :: coor
	atmcrd = coor
end subroutine qm3_update_coor

subroutine qm3_update_chrg( chrg )
	use dynamo
	implicit none
	real*8, dimension(1:natoms), intent( in ) :: chrg
	atmchg = chrg
end subroutine qm3_update_chrg

subroutine qm3_get_func( coor, func )
	use dynamo
	implicit none
	real*8, dimension(1:3,1:natoms), intent( in ) :: coor
	real*8, intent( out ) :: func
	call qm3_update_coor( coor )
	call energy
	func = etotal
end subroutine qm3_get_func

subroutine qm3_get_grad( coor, func, grad )
	use dynamo
	implicit none
	real*8, dimension(1:3,1:natoms), intent( in ) :: coor
	real*8, intent( out ) :: func
	real*8, dimension(1:3,1:natoms), intent( out ) :: grad
	call qm3_update_coor( coor )
	call gradient
	func = etotal
	grad = atmder
end subroutine qm3_get_grad
