program de_mierda
    use dynamo
    implicit none
    integer :: i, t
    logical, dimension(:), allocatable :: acs

    call dynamo_header

    call mm_file_process( "borra", "opls" )
    call mm_system_construct( "borra", "seq" )

    allocate( acs(1:natoms) )

    acs = atom_selection( subsystem = (/ "A" /), residue_number = (/ 144 /), &
        atom_name = (/ "C", "O" /) )
    acs = acs .or. atom_selection( subsystem = (/ "A" /), residue_number = (/ 145 /), &
        atom_name = (/ "N  ", "H43", "CA ", "CB ", "SG ", "C34", "C24", "H63", "N23", &
        "H62", "C21", "O22", "C26", "H64", "H65", "C27", "H66", "C29", "H67", "H68", &
        "C30", "H69", "H70", "N31", "H71", "C32", "O33", "O35", "H72", "C36", "H73", &
        "H74", "O37", "H75", "H47", "H48", "H46", "C  ", "O  " /) )
    acs = acs .or. atom_selection( subsystem = (/ "A" /), residue_number = (/ 146 /), &
        atom_name = (/ "N  ", "H  ", "CA ", "HA3", "HA2" /) )
    acs = acs .or. atom_selection( subsystem = (/ "A" /), residue_number = (/ 41 /), &
        atom_name = (/ "CB ", "HB3", "HB2", "CG ", "ND1", "CE1", "HE1", "NE2", "HE2", &
        "CD2", "HD2" /) )

    open( unit = 999, file = "data", action = "write", form = "unformatted" )
    write( 999 ) atmchg
    ! - ATMEPS(I)   = 2.0_DP * SQRT ( TYPES(MM)%EPSILON )
    write( 999 ) ( atmeps * 0.5d0 ) ** 2
! ------------------------------------------------------------
!    ! - ATMSIG(I)   = 0.5_DP * TYPES(MM)%SIGMA
!    write( 999 ) ( atmsig * 2.0d0 ) * 0.5612310241546865d0
! ------------------------------------------------------------
    ! - ATMSIG(I)   = SQRT ( TYPES(MM)%SIGMA )
    write( 999 ) ( atmsig ** 2 ) * 0.5612310241546865d0
! ------------------------------------------------------------
    close( 999 )
    ! - exlcusions -
    write( *, "(a)" ) "        self.exc = []"
    write( *, "(a)" ) "        # ------------------------------"
    do i = 1, nbonds
        t = 0
        if( acs( bonds(i)%i ) ) t = t + 1
        if( acs( bonds(i)%j ) ) t = t + 1
        if( t == 0 .or. t == 2 ) cycle
        write( *, "(a,f8.1,a,f8.3,a,i6,a,i6,a)" ) &
                "        self.exc.append( qm3.engines.mmres.distance(", &
                bonds(i)%fc, ", ", bonds(i)%eq, &
                ", [", bonds(i)%i - 1, ", ", bonds(i)%j - 1, " ] ) )"
        write( *, "(a)" ) "        self.exc[-1].ffac = 0.0"
        if( acs( bonds(i)%i ) ) then
            write( *, "(a)" ) "        self.exc[-1].gfac = [ 1.0, 0.0 ]"
        else
            write( *, "(a)" ) "        self.exc[-1].gfac = [ 0.0, 1.0 ]"
        end if
    end do
    write( *, "(a)" ) "        # ------------------------------"
    do i = 1, nangles
        t = 0
        if( acs( angles(i)%i ) ) t = t + 1
        if( acs( angles(i)%j ) ) t = t + 1
        if( acs( angles(i)%k ) ) t = t + 1
        if( t == 0 .or. t == 2 .or. t == 3 ) cycle
        write( *, "(a,f8.1,a,f8.3,a,i6,a,i6,a,i6,a)" ) &
                "        self.exc.append( qm3.engines.mmres.angle(", &
                angles(i)%fc, ", ", angles(i)%eq, &
                ", [", angles(i)%i - 1, &
                ", ", angles(i)%j - 1, &
                ", ", angles(i)%k - 1, " ] ) )"
        write( *, "(a)" ) "        self.exc[-1].ffac = 0.0"
        if( acs( angles(i)%i ) ) then
            write( *, "(a)" ) "        self.exc[-1].gfac = [ 1.0, 0.0, 0.0 ]"
        else if( acs( angles(i)%j ) ) then
            write( *, "(a)" ) "        self.exc[-1].gfac = [ 0.0, 1.0, 0.0 ]"
        else
            write( *, "(a)" ) "        self.exc[-1].gfac = [ 0.0, 0.0, 1.0 ]"
        end if
    end do
    write( *, "(a)" ) "        # ------------------------------"
    do i = 1, ndihedrals
        t = 0
        if( acs( dihedrals(i)%i ) ) t = t + 1
        if( acs( dihedrals(i)%j ) ) t = t + 1
        if( acs( dihedrals(i)%k ) ) t = t + 1
        if( acs( dihedrals(i)%l ) ) t = t + 1
        if( t == 0 .or. t == 3 .or. t == 4 ) cycle
        write( *, "(a,i3,a,f8.4,a,f6.1,a,a,i6,a,i6,a,i6,a,i6,a)" ) &
                "        self.exc.append( qm3.engines.mmres.dihedral( { ", &
                dihedrals(i)%period, ": [ ", dihedrals(i)%fc, ", ", 180.d0 * dihedrals(i)%phase," ] }", &
                ", [", dihedrals(i)%i - 1, &
                ", ", dihedrals(i)%j - 1, &
                ", ", dihedrals(i)%k - 1, &
                ", ", dihedrals(i)%l - 1, " ] ) )"
        write( *, "(a)" ) "        self.exc[-1].ffac = 0.0"
        write( *, "(a)", advance = "no" ) "        self.exc[-1].gfac = [ "
        if( acs( dihedrals(i)%i ) ) then
            write( *, "(a)", advance = "no" ) "1.0, "
        else
            write( *, "(a)", advance = "no" ) "0.0, "
        end if
        if( acs( dihedrals(i)%j ) ) then
            write( *, "(a)", advance = "no" ) "1.0, "
        else
            write( *, "(a)", advance = "no" ) "0.0, "
        end if
        if( acs( dihedrals(i)%k ) ) then
            write( *, "(a)", advance = "no" ) "1.0, "
        else
            write( *, "(a)", advance = "no" ) "0.0, "
        end if
        if( acs( dihedrals(i)%l ) ) then
            write( *, "(a)" ) "1.0 ]"
        else
            write( *, "(a)" ) "0.0 ]"
        end if
    end do

end
