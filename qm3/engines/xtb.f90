!
! gfortran -shared|-dynamiclib -o libxtb.so -Ixtb-mod xtb.f90 libxtb.a -framework Accelerate
!
module qm3
    use xtb_type_environment
    use xtb_type_molecule
    use xtb_type_wavefunction
    use xtb_type_basisset
    use xtb_type_param
    use xtb_type_data
    use xtb_type_pcem
    use xtb_setparam
    use xtb_basis
    use xtb_scf
    use xtb_scc_core
    use xtb_paramset
    use xtb_xtb_data
    use xtb_xtb_gfn2

    use xtb_eeq
    use xtb_chargemodel
    use xtb_disp_ncoord, only: ncoord_erf
    use xtb_type_solvation, only : TSolvation

    implicit none
    public
    integer, parameter  :: maxiter = 1000
    integer, parameter  :: prlevel = 0
    logical, parameter  :: lgrad = .true.
    logical             :: restart = .false.
    real*8, parameter   :: et = 300.0d0
    real*8, parameter   :: acc = 1.0d0
    type(TEnvironment)  :: env
    type(TMolecule)     :: mol
    type(scc_results)   :: res
    type(TBasisset)     :: basis
    type(TWavefunction) :: wfn
    type(tb_pcem)       :: pcem
    type(TxTBData)      :: xtbData
    class(TSolvation), allocatable :: solvation
    real*8              :: ene, gap
    type(TxTBParameter) :: globpar
    logical             :: okpar, okbas, okrun
    real*8, parameter   :: AA__Bohr = 1.0d0 / 0.52917726d0
end module qm3


subroutine qm3_xtb_calc( nQM, nMM, siz, dat )
    use qm3
    implicit none
    integer, intent( in ) :: nQM, nMM, siz
    real*8, dimension(0:siz-1), intent( inout ) :: dat
   
    integer, dimension(nQM)  :: atn
    real*8, dimension(3,nQM) :: xyz, grd
    integer                  :: i, i3, j
    real*8, allocatable      :: cn(:)
    type(chrg_parameter)     :: chrgeq

    ! 3 + nQM [QM_chg] + 3 * nQM [QM_crd/grd] + nQM [QM_mul] + nMM [MM_chg] + 3 * nMM [MM_crd/grd]
    do i = 1, nQM
        atn(i) = dint( dat(2+i) )
        j = 3 + nQM + 3 * ( i - 1 )
        xyz(1,i) = dat(j)   * AA__Bohr
        xyz(2,i) = dat(j+1) * AA__Bohr
        xyz(3,i) = dat(j+2) * AA__Bohr
    end do
    if( nMM > 0 ) then
        if( .not. restart ) call pcem%allocate( nMM )
        pcem%gam = 999.0d0
        do i = 1, nMM
            pcem%q(i) = dat(2+5*nQM+i)
            j = 3 + 5 * nQM + nMM + 3 * ( i - 1 )
            pcem%xyz(1,i) = dat(j)   * AA__Bohr 
            pcem%xyz(2,i) = dat(j+1) * AA__Bohr
            pcem%xyz(3,i) = dat(j+2) * AA__Bohr
        end do
    end if

    if( .not. restart ) then
        set%gfn_method = 2
        call init( env )
        call init( mol, atn, xyz )
        wfn%nel = dint( sum( mol%z ) -  dat(1) )
        wfn%nopen = dint( dat(2) )
        call use_parameterset( "param_gfn2-xtb.txt", globpar, xtbData, okpar )
        call newBasisset( xtbData, mol%n, mol%at, basis, okbas )
        call wfn%allocate( mol%n, basis%nshell, basis%nao )
        mol%chrg = dat(1)
        mol%uhf  = dint( dat(2) )
        !wfn%q = mol%chrg / real( mol%n, kind=8 )

        ! EEQ guess
        allocate( cn(mol%n), source = 0.0d0 )
        call new_charge_model_2019( chrgeq, mol%n, mol%at )
        call ncoord_erf( mol%n, mol%at, mol%xyz, cn )
        call eeq_chrgeq( mol, env, chrgeq, cn, wfn%q )
        deallocate( cn )

        call iniqshell( xtbData, mol%n, mol%at, mol%z, basis%nshell, wfn%q, wfn%qsh, set%gfn_method )
    else
        mol%xyz = xyz
    end if

    grd = 0.0d0
    if( nMM > 0 ) pcem%grd = 0.0d0
    call scf( env, mol, wfn, basis, pcem, xtbData, solvation, &
        gap, et, maxiter, prlevel, restart, lgrad, acc, ene, grd, res )
    call env%check( okrun )

    dat(0) = ene
    do i = 1, nQM
        j = 3 + nQM + 3 * ( i - 1 )
        dat(j)   = grd(1,i)
        dat(j+1) = grd(2,i)
        dat(j+2) = grd(3,i)
        dat(3+4*nQM+i-1) = wfn%q(i)
    end do
    if( nMM > 0 ) then
        do i = 1, nMM
            j = 3 + 5 * nQM + nMM + 3 * ( i - 1 )
            dat(j)   = pcem%grd(1,i)
            dat(j+1) = pcem%grd(2,i)
            dat(j+2) = pcem%grd(3,i)
        end do
    end if
    if( .not. restart ) restart = .true.

end subroutine qm3_xtb_calc


subroutine qm3_xtb_clean
    use qm3
    implicit none
    call mol%deallocate
    call wfn%deallocate
    call basis%deallocate
    call pcem%deallocate
end subroutine qm3_xtb_clean
