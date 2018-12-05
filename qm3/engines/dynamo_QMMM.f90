subroutine driver( acs, flg )
    use dynamo
    implicit none
    logical, dimension(1:natoms), intent( in ) :: acs, flg
    character( len = 256 ) :: str
    integer :: i, j, n
    real*8, dimension(:,:), allocatable :: crd, der
    integer, dimension(:), allocatable :: idx

    n = count( acs )
    allocate( crd(1:3,1:n), der(1:3,1:n), idx(1:n) )
    j = 1
    do i = 1, natoms
        if( acs(i) ) then
            idx(j) = i
            j = j + 1
        end if
    end do

    call getarg( 1, str )
    write(*,"(a,a,a)") "[", trim( str ), "]"
    open( file = trim( str ), unit = 998, action = "read", form = "formatted" )
    read( 998, "(a)" ) str
    write(*,"(a,a,a)") "[", trim( str ), "]"
    do while( trim( str ) /= "exit" )
        if( trim( str ) == "coordinates" ) then
            open( file = "dynamo.crd", unit = 999, action = "read", form = "unformatted" )
            read( 999 ) crd(1:3,1:n)
            close( 999 )
            do i = 1, n
                atmcrd(1:3,idx(i)) = crd(1:3,i)
            end do
        end if
        if( trim( str ) == "energy" ) then
            call atoms_fix( .not. flg )
            call optimize_lbfgs( nupd = 9, print_frequency = 100, step_number = 1000, &
                                 step_size = 0.01d0, gradient_tolerance = 1.0d0 )
            call coordinates_write( "last.crd" )
            call atoms_fix( .not. acs )
            call energy
            write( *, "(6f20.10)" ) eqm, qmlj
            open( file = "dynamo.dat", unit = 999, action = "write", form = "unformatted" )
            write( 999 ) eqm + qmlj
            close( 999 )
        end if
        if( trim( str ) == "gradient" ) then
            call atoms_fix( .not. flg )
            call optimize_lbfgs( nupd = 9, print_frequency = 100, step_number = 1000, &
                                 step_size = 0.01d0, gradient_tolerance = 1.0d0 )
            call coordinates_write( "last.crd" )
            call atoms_fix( .not. acs )
            call gradient
            write( *, "(6f20.10)" ) eqm, qmlj
            do i = 1, n
                der(1:3,i) = atmder(1:3,idx(i))
            end do
            open( file = "dynamo.dat", unit = 999, action = "write", form = "unformatted" )
            write( 999 ) eqm + qmlj
            write( 999 ) der(1:3,1:n)
            close( 999 )
        end if
        read( 998, "(a)" ) str
        write(*,"(a,a,a)") "[", trim( str ), "]"
    end do
    close( 998 )
    deallocate( crd, der )
end subroutine driver


program slave
    use dynamo
    implicit none

    logical, dimension(:), allocatable :: flg, acs
    integer :: i

    call dynamo_header

    call mm_file_process( "borra", "opls" )
    call mm_system_construct( "borra", "seq" )
    call coordinates_read( "seed" )

    allocate( flg(1:natoms), acs(1:natoms) )
    acs = atom_selection ( subsystem = (/ "X" /) )
    call my_sele( flg )

    call mopac_setup( "AM1", charge = -2, selection = acs )

    call energy_initialize
    call energy_non_bonding_options( &
        list_cutoff   = 18.0_dp, &
        outer_cutoff  = 16.0_dp, &
        inner_cutoff  = 14.5_dp, &
        minimum_image = .true. )

    call driver( acs, flg )

    deallocate( acs, flg )
end program
include "dentro.f90"
