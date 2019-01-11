subroutine driver
    use dynamo
    implicit none
    character( len = 256 ) :: str
    integer :: i, j

    call getarg( 1, str )
    write(*,"(a,a,a)") "[", trim( str ), "]"
    open( file = trim( str ), unit = 998, action = "read", form = "formatted" )
    read( 998, "(a)" ) str
    write(*,"(a,a,a)") "[", trim( str ), "]"
    do while( trim( str ) /= "exit" )
        if( trim( str ) == "charges" ) then
            open( file = "dynamo.chrg", unit = 999, action = "read", form = "unformatted" )
            read( 999 ) atmchg
            close( 999 )
        end if
        if( trim( str ) == "coordinates" ) then
            open( file = "dynamo.crd", unit = 999, action = "read", form = "unformatted" )
            read( 999 ) atmcrd(1:3,1:natoms)
            close( 999 )
        end if
        if( trim( str ) == "energy" ) then
            call energy
            open( file = "dynamo.dat", unit = 999, action = "write", form = "unformatted" )
            write( 999 ) etotal
            close( 999 )
        end if
        if( trim( str ) == "gradient" ) then
            if( .not. allocated( atmder ) ) allocate( atmder(1:3,1:natoms ) )
            call gradient
            open( file = "dynamo.dat", unit = 999, action = "write", form = "unformatted" )
            write( 999 ) etotal
            write( 999 ) atmder
            close( 999 )
        end if
        read( 998, "(a)" ) str
        write(*,"(a,a,a)") "[", trim( str ), "]"
    end do
    close( 998 )
end subroutine driver


program slave
    use dynamo
    implicit none

    logical, dimension(:), allocatable :: flg
    integer :: i

    call dynamo_header

    call mm_file_process( "borra", "opls" )
    call mm_system_construct( "borra", "seq" )
    call coordinates_read( "crd" )

    allocate( flg(1:natoms) )
    flg = .false.
    flg = atom_selection( subsystem = (/ "W" /), residue_number = (/ 1 /) )
    call atoms_fix( flg )
    do i = 1, natoms
        if( flg(i) ) then
            atmchg(i) = 0.0_dp
            atmchg14(i) = 0.0_dp
        end if
    end do

    call energy_initialize
    call energy_non_bonding_options( &
        list_cutoff   = 18.0_dp, &
        outer_cutoff  = 16.0_dp, &
        inner_cutoff  = 14.0_dp, &
        minimum_image = .false. )

    call driver

    deallocate( flg )
end program
