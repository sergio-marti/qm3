module qm3
   use qm2_pm6_hof_module, only : strlen
   use file_io_dat, only : MAX_FN_LEN
   implicit none
   public
   real*8 x(3000), f(3000), escf
   character(len=8) atnam(1000)
   real*8 born_radii(1000), one_born_radii(1000)
   real*8 intdiel, extdiel, Arad
   integer natom, ier, atnum(1000)
   character(len=80) arg
   integer ntpr
   character(len=MAX_FN_LEN) mdin, mdout 
   real*8 excharge(40000)
   integer chgatnum(10000)
   character(len=8) chgnam(10000)
   integer ncharge
   integer :: igb, maxcyc
   real*8  :: grms_tol
   real*8  :: total_energy
   logical :: master=.true.
end module qm3


subroutine qm3_sqm_init
   use qmmm_module, only : qmmm_nml, qmmm_struct, qmmm_mpi, qm2_struct
   use sqm_qmmm_read_and_alloc, only : read_qmmm_nm_and_alloc
   use constants, only : KCAL_TO_EV, EV_TO_KCAL
   use file_io_dat, only : MAX_FN_LEN
   use UtilitiesModule, only : print
   use qm3
   implicit none
   qmmm_struct%qm_mm_first_call = .true.
   qmmm_struct%fock_first_call = .true.
   qmmm_struct%fock2_2atm_first_call = .true.
   ! qmmm_struct%qm2_deriv_qm_analyt_first_call = .true.
   qmmm_struct%qm2_allocate_e_repul_first_call = .true.
   ! qmmm_struct%qm2_rotate_qmqm_first_call = .true.
   qmmm_struct%qm2_calc_rij_eqns_first_call = .true.
   qmmm_struct%qm2_scf_first_call = .true.
   qmmm_struct%zero_link_charges_first_call = .true.
   qmmm_struct%adj_mm_link_pair_crd_first_call = .true.
   mdin   = 'mdin'
   mdout  = 'mdout'
   igb = 0
   call amopen( 5, mdin, 'O', 'F', 'R' )
   call x_getsqmx(natom,x,atnam,atnum,ncharge,excharge,chgnam,chgatnum)
   call read_qmmm_nm_and_alloc(natom,igb,atnum,maxcyc,grms_tol,ntpr, &
                               ncharge,excharge,chgatnum )
   call qm_assign_atom_types
   qmmm_mpi%commqmmm_master = master
   qmmm_mpi%numthreads = 1
   qmmm_mpi%mytaskid = 0
   qmmm_mpi%natom_start = 1
   qmmm_mpi%natom_end = natom
   qmmm_mpi%nquant_nlink_start = 1
   qmmm_mpi%nquant_nlink_end = qmmm_struct%nquant_nlink
   call allocate_qmgb(qmmm_struct%nquant_nlink)
   allocate( qmmm_struct%dxyzqm(3, qmmm_struct%nquant_nlink) )
   allocate ( qm2_struct%scf_mchg(qmmm_struct%nquant_nlink) )
   close( 5 )
end subroutine qm3_sqm_init



subroutine qm3_sqm_calc( siz, dat )
   use qmmm_module, only : qm_gb, qmmm_nml, qmmm_struct, qmmm_mpi, qm2_struct
   use qm2_dftb_module, only : ks_struct
   use qm2_pm6_hof_module, only : cct, nsp2, print, strlen
   use constants, only : EV_TO_KCAL, KCAL_TO_EV, zero, one, alpb_alpha
   use UtilitiesModule, only : print
   use qm3
   implicit none
   integer, intent( in ) :: siz
   real*8, dimension(0:siz-1), intent( inout ) :: dat
   integer :: i, j, i3, m
   real*8 :: alpb_beta
   real*8 , dimension(:), allocatable :: vectmp1, vectmp2, vectmp3, vectmp4

   do i = 1, 3 * ( natom - ncharge )
      x(i) = dat(i-1)
   end do
   do i = 1, ncharge
      j = 3 * ( natom - ncharge + i - 1 )
      qmmm_struct%qm_xcrd(1,i) = dat(j)
      qmmm_struct%qm_xcrd(2,i) = dat(j+1)
      qmmm_struct%qm_xcrd(3,i) = dat(j+2)
   end do

   call amopen( 6, mdout,'R', 'F', 'W' )
   if (qmmm_struct%qm_mm_first_call) then
     allocate ( qmmm_struct%qm_coords(3,qmmm_struct%nquant_nlink) )
     if (qmmm_nml%qmgb == 2) then
       if (qm_gb%alpb_on) then
         alpb_beta=alpb_alpha*(intdiel/extdiel)
         qm_gb%intdieli = one/(intdiel*(one + alpb_beta))
         qm_gb%extdieli = one/(extdiel*(one + alpb_beta))
         qm_gb%one_Arad_beta = alpb_beta/Arad
       else
         qm_gb%intdieli = 1.0d0/intdiel
         qm_gb%extdieli = 1.0d0/extdiel
       end if
       qm_gb%mmcut2 = 999.d0
     end if
   end if
   i3 = 0
   do i=1,qmmm_struct%nquant_nlink
      qmmm_struct%qm_coords(1,i) = x(i3+1)
      qmmm_struct%qm_coords(2,i) = x(i3+2)
      qmmm_struct%qm_coords(3,i) = x(i3+3)
      i3 = i3 + 3
   end do
   if(qmmm_struct%qm_mm_first_call) call qm2_load_params_and_allocate(.false.)
   call qm2_calc_rij_and_eqns(qmmm_struct%qm_coords, qmmm_struct%nquant_nlink, &
          qmmm_struct%qm_xcrd, natom, qmmm_struct%qm_mm_pairs)
   call qm2_energy(escf, qm2_struct%scf_mchg, natom, born_radii, one_born_radii)
   call qm2_calc_dipole(x)
   qmmm_struct%qm_mm_first_call = .false.
   total_energy = qmmm_struct%elec_eng +  qmmm_struct%enuclr_qmqm + qmmm_struct%enuclr_qmmm

   qmmm_struct%dxyzqm=zero
   if (qmmm_nml%qmtheory%DFTB) then
      call qm2_dftb_get_qm_forces(qmmm_struct%dxyzqm)
      qmmm_struct%dxyzcl = zero
      allocate( vectmp1(qmmm_struct%qm_mm_pairs), vectmp2(qmmm_struct%qm_mm_pairs), &
                vectmp3(qmmm_struct%qm_mm_pairs), vectmp4(qmmm_struct%qm_mm_pairs) )
      call qm2_dftb_get_qmmm_forces(qmmm_struct%dxyzcl,qmmm_struct%dxyzqm, vectmp1,vectmp2,vectmp3,vectmp4)
      deallocate( vectmp1, vectmp2, vectmp3, vectmp4 )
   else
      call qm2_get_qm_forces(qmmm_struct%dxyzqm)
      qmmm_struct%dxyzcl = zero
      call qm2_get_qmmm_forces(qmmm_struct%dxyzqm,qmmm_struct%qm_xcrd,qmmm_struct%dxyzcl,qm2_struct%scf_mchg)
   end if
   close( 6 )

   dat(0) = escf
   do i = 1, natom - ncharge
      dat(i) = qm2_struct%scf_mchg( i )
   end do
   do i = 1, natom - ncharge
      j = 1 + ( natom - ncharge ) + ( i - 1 ) * 3
      dat(j)   = qmmm_struct%dxyzqm(1,i)
      dat(j+1) = qmmm_struct%dxyzqm(2,i)
      dat(j+2) = qmmm_struct%dxyzqm(3,i)
   end do
   do i = 1, ncharge
      j = 1 + 4 * ( natom - ncharge ) + ( i - 1 ) * 3
      dat(j)   = qmmm_struct%dxyzcl(1,i)
      dat(j+1) = qmmm_struct%dxyzcl(2,i)
      dat(j+2) = qmmm_struct%dxyzcl(3,i)
   end do
end subroutine qm3_sqm_calc



subroutine qm3_sqm_clean
   use qmmm_module, only : qmmm_nml, qmmm_struct, qmmm_vsolv, qm2_params, deallocate_qmmm
   implicit none
   call deallocate_qmmm(qmmm_nml, qmmm_struct, qmmm_vsolv, qm2_params)
   call mexit(6,0)
end subroutine qm3_sqm_clean



subroutine qm3_sqm_writedens
   use qmmm_module, only : qm2_struct
   implicit none
   open( unit = 666, file = "sqm_dens", action = "write", form = "unformatted" )
   write( 666 ) qm2_struct%den_matrix(1:qm2_struct%matsize)
   close( 666 )
end subroutine qm3_sqm_writedens



subroutine x_getsqmx(natom,x,atnam,atnum,ncharge,excharge,chgnam,chgatnum)
   implicit none
   real*8 x(*)
   integer i,i3,lun
   integer natom,atnum(*)
   character(len=8) atnam(*)
   character(len=80) line
   real*8 excharge(*)
   integer chgatnum(*)
   character(len=8) chgnam(*)
   integer ncharge
   integer ia, ic, ihead, iend
   logical mdin_external_charge
   lun = 5 
   mdin_external_charge = .false.
   ncharge = 0
   ihead=0
   iend=0
   do i=1,11999
      read(lun,'(a)',end=10) line
      if (line(1:1) == "#") then
         if (line(1:80) == "#EXCHARGES") then
            mdin_external_charge = .true.
            ihead = ihead + 1
         else if (line(1:80) == "#END") then
            iend = iend + 1
         else
            write(0,*) 'Unrecognized header name'
            write(0,*) line(1:80)
            call mexit(6,1)
         end if
      end if
   end do
   10 if (iend < ihead) then
      write(0,*) 'Missing "#END" termination sign, exit program'
      call mexit(6,1)
   end if
   rewind(lun)
   do i=1,20
      read(5,'(a)') line
      if( line(1:2) == " /" ) go to 11
   end do
   write(0,*) 'Error in finding end of qmmm namelist'
   call mexit(6,1)
   11 i3=0
   ia=0
   do i=1,999
      read(lun,'(a)',end=12) line
      if (line(1:80) /= "") then
         if (line(1:80) /= "#EXCHARGES") then
            ia = ia + 1
            read(line,*,err=15) atnum(ia),atnam(ia),x(i3+1),x(i3+2),x(i3+3)
            atnum(ia) = abs( atnum(ia) )
            i3 = i3 + 3
         else
            go to 12
         end if
      end if
   end do
   12 natom = ia
   if (mdin_external_charge) then
      i3=0
      ic=0
      do i=1,9999
         read(lun,'(a)',end=14) line
         if (line(1:80) /= "") then
            if (line(1:80) /= "#END") then
               ic = ic + 1
               read(line,*,err=16) chgatnum(ic), chgnam(ic),excharge(i3+1:i3+4) 
               i3 = i3 + 4
            else
               go to 13
            end if
         end if
      end do
   13 ncharge = ic
   end if
   return
   14 write(0,*) 'The termination sign "#END" is missing'
   call mexit(6,1)
   15 write(0,*) 'Error in reading QM atoms'
   call mexit(6,1)
   16 write(0,*) 'Error in reading external charges'
   call mexit(6,1)
end subroutine x_getsqmx 
