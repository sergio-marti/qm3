subroutine qm3_lio_init
    use garcha_mod, only: ntatom, natom, nsol, Iz, mulliken, writeforces, writexyz
    use fileio_data, only: verbose

    implicit none
    integer :: ierr
    character(len=20) :: inpfile

    call lio_defaults()
    writexyz = .false.
    inpfile = "lio.inp"
    call read_options(inpfile, ierr)
    inpfile = "lio.xyz"
    verbose = 0
    mulliken = .false.
    writeforces = .false.
    call read_coords( inpfile )
    call init_lio_common( natom, Iz, nsol, 0 )
end subroutine qm3_lio_init


subroutine qm3_lio_calc( siz, dat )
    use garcha_mod, only: ntatom, natom, nsol, Smat, RealRho, sqsm, &
                           Eorbs, Eorbs_b, d, r, Iz, Em, Rm, pc, rqm, Pmat_vec
    use basis_data, only: M, MM, Nuc
    use faint_cpu, only: int1
    use SCF_aux, only: fix_densmat
    implicit none
    integer, intent( in ) :: siz
    real*8, dimension(0:siz-1), intent( inout ) :: dat

    real*8, dimension(:,:), allocatable :: qm_grd, mm_grd
    real*8, dimension(:), allocatable :: qm_chg, Fock_1e, Hmat
    real*8 :: ene, tmp
    integer :: i, j, k

    if( .not. allocated( Smat ) )    allocate( Smat(M,M) )
    if( .not. allocated( RealRho ) ) allocate( RealRho(M,M) )
    if( .not. allocated( sqsm ) )    allocate( sqsm(M,M) )
    if( .not. allocated( Eorbs ) )   allocate( Eorbs(M) )
    if( .not. allocated( Eorbs_b ) ) allocate( Eorbs_b(M) )

    allocate( qm_chg(1:natom), qm_grd(1:3,1:natom) )
    if( nsol > 0 ) allocate( mm_grd(1:3,1:nsol) )

    k = 0
    do i = 1, natom
        do j = 1, 3
            r(i,j)   = dat(k) / 0.529177d0
            rqm(i,j) = dat(k) / 0.529177d0
            k = k + 1
        end do
    end do
    do i = 1, nsol
        do j = 1, 3
            r(natom+i,j) = dat(k) / 0.529177d0
            k = k + 1
        end do
    end do

    call recenter_coords( rqm, r, natom, nsol )

    call SCF( ene )

    allocate( Fock_1e(1:MM), Hmat(1:MM) )
    call int1( tmp, Fock_1e, Hmat, Smat, d, r, Iz, natom, ntatom )
    call spunpack( 'L', M, Fock_1e, Smat )
    call spunpack( 'L', M, Pmat_vec(1), RealRho )
    deallocate( Fock_1e, Hmat )
    call fix_densmat( RealRho )
    do i = 1, natom
         qm_chg(i) = real( Iz(i), kind=8 )
    enddo
    call mulliken_calc( natom, M, RealRho, Smat, Nuc, qm_chg )

    call dft_get_qm_forces( qm_grd )
    if( nsol > 0 ) call dft_get_mm_forces( mm_grd, qm_grd )

    dat(0) = ene
    do i = 1, natom
        dat(i) = qm_chg( i )
    end do
    do i = 1, natom
        j = 1 + natom + ( i - 1 ) * 3
        dat(j)   = qm_grd(1,i)
        dat(j+1) = qm_grd(2,i)
        dat(j+2) = qm_grd(3,i)
    end do
    do i = 1, nsol
        j = 1 + 4 * natom + ( i - 1 ) * 3
        dat(j)   = mm_grd(1,i)
        dat(j+1) = mm_grd(2,i)
        dat(j+2) = mm_grd(3,i)
    end do

    deallocate( qm_chg, qm_grd )
    if( nsol > 0 ) deallocate( mm_grd )
end subroutine qm3_lio_calc


subroutine qm3_lio_clean
    call lio_finalize()
end subroutine qm3_lio_clean
