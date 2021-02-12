!-- MOPAC taken (and cleaned) from fDynamo-2.2
! =======================================================================================
!program proof
!    use mopac
!    implicit none
!    real*8 :: nuc, scf
!    real*8, dimension(:), allocatable :: chg
!    natm = 6
!    naqm = 3
!    allocate( atmchg(1:natm), atmcrd(1:3,1:natm), atmder(1:3,1:natm), chg(1:naqm) )
!    atmchg = (/ 8.0d0, 1.0d0, 1.0d0, -0.834d0, 0.417d0, 0.417d0 /)
!    atmcrd(1:3,1) = (/  0.12109d0,  0.06944d0, -0.22458d0 /)
!    atmcrd(1:3,2) = (/ -0.52694d0,  0.16499d0, -0.94583d0 /)
!    atmcrd(1:3,3) = (/  0.70159d0, -0.63565d0, -0.54677d0 /)
!    atmcrd(1:3,4) = (/ -0.45114d0,  1.12675d0,  2.21102d0 /)
!    atmcrd(1:3,5) = (/ -0.29157d0,  0.59483d0,  1.39876d0 /)
!    atmcrd(1:3,6) = (/  0.05804d0,  1.92714d0,  2.01036d0 /)
!
!    call mopac_setup( "AM1" )
!    call mopac_scf_options()
!
!    atmder = 0.0d0
!    call mopac_integrals( nuc, atmder )
!    call mopac_scf( scf, 200, .true. )
!    write( *, "(f20.10)" ) ev_to_kj * ( scf + nuc ) + atheat
!    call mopac_gradients( atmder )
!    write( *, "(3f20.10)" ) atmder
!    write( *, "(f20.10)" ) dsqrt( sum( atmder ** 2 ) / ( 3.0d0 * real( natm, 8 ) ) )
!    call density_write( "dens.bin" )
!    call mopac_charges( chg )
!    write( *, "(3f20.10)" ) chg
!end
! ---------------------------------------------------------------------------------------
!     -272.0275070986
!       16.9469154876       18.2125899843       60.0178061824
!      -29.6607486482       -0.0216901028      -42.3696920893
!       13.6228622679      -21.9914204108      -19.8080426929
!       12.7198650069      -15.1441325253      -36.0956811568
!      -12.3670427743        7.2988550752       26.0667237078
!       -1.2618513400       11.6457979794       12.1888860490
!       24.5210375472
!       -0.4299752861        0.2144379951        0.2155372910
! =======================================================================================
MODULE MOPAC

IMPLICIT NONE
PRIVATE
PUBLIC :: NATM, NAQM, ATMCHG, ATMCRD, ATMDER, ATHEAT, EV_TO_KJ, &
          MOPAC_INTEGRALS, MOPAC_SETUP, MOPAC_SCF_OPTIONS, MOPAC_SCF, &
          MOPAC_GRADIENTS, MOPAC_CHARGES, MOPAC_DATA_INITIALIZE, &
          DENSITY_WRITE, DENSITY_READ, CUT_OFF, CUT_ON
SAVE

REAL*8, PARAMETER :: AU_TO_EV = 27.21D0
REAL*8, PARAMETER :: KCAL_TO_KJ = 4.184D0
REAL*8, PARAMETER :: EV_TO_KJ = KCAL_TO_KJ * 23.061D0
REAL*8, PARAMETER :: BOHRS_TO_ANGSTROMS = 0.529177249D0
REAL*8, PARAMETER :: ANGSTROMS_TO_BOHRS = 1.0D0 / BOHRS_TO_ANGSTROMS
REAL*8, PARAMETER :: PI32 = 5.56832799683170D0, RLN10 = 2.30258D0, SQRT3 = 1.73205080756888D0

INTEGER, PARAMETER :: NELEMENTS = 100
INTEGER, DIMENSION(1:NELEMENTS) :: NATORB = 0
LOGICAL, DIMENSION(1:NELEMENTS) :: SEPAR = .FALSE.
REAL*8, DIMENSION(1:NELEMENTS) :: AD    = 0.0D0, ALP   = 0.0D0, AM    = 0.0D0, AQ   = 0.0D0, &
    BETAD = 0.0D0, BETAP = 0.0D0, BETAS = 0.0D0, CORE = 0.0D0, &
    DD    = 0.0D0, EHEAT = 0.0D0, EISOL = 0.0D0, GDD  = 0.0D0, &
    GPD   = 0.0D0, GSD   = 0.0D0, GP2   = 0.0D0, GPP  = 0.0D0, &
    GSP   = 0.0D0, GSS   = 0.0D0, HSP   = 0.0D0, QQ   = 0.0D0, &
    USS   = 0.0D0, UPP   = 0.0D0, UDD   = 0.0D0, ZD   = 0.0D0, &
    ZP    = 0.0D0, ZS    = 0.0D0
REAL*8, DIMENSION(1:NELEMENTS,1:4) :: FN1 = 0.0D0, FN2 = 0.0D0, FN3 = 0.0D0
REAL*8, DIMENSION(1:NELEMENTS,1:2) :: PDDGC = 0.0D0, PDDGE = 0.0D0

INTEGER :: NATM = 0, NAQM = 0
REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ATMCRD, ATMDER
REAL*8, ALLOCATABLE, DIMENSION(:) :: ATMCHG
INTEGER, ALLOCATABLE, DIMENSION(:) :: BFIRST, BLAST
CHARACTER ( LEN = 4 ) :: HAMILTONIAN
INTEGER :: MULTIP = 1, NALPHA = 0, NBETA = 0, NELEC = 0, N2ELEC = 0, TOTCHG = 0
REAL*8 :: ATHEAT  = 0.0D0
INTEGER :: NBASIS = 0, NBASTR = 0
INTEGER, ALLOCATABLE, DIMENSION(:) :: BFINDEX
REAL*8, ALLOCATABLE, DIMENSION(:) :: BETA, HCORE, OVERLAP, SETEI, USPD
REAL*8, ALLOCATABLE, DIMENSION(:) :: DENMAT, DENMATA, DENMATB
INTEGER, DIMENSION(1:16,1:4,1:4,1:2) :: JIND2
INTEGER, DIMENSION(1:4,1:4,1:16,1:2) :: JIND3
INTEGER, DIMENSION(1:16,1:4,1:4,1:2) :: KIND2
INTEGER :: NSHELL = 0
INTEGER, ALLOCATABLE, DIMENSION(:) :: KSTART, KATOM, KTYPE, KNG, KLOC, KMIN, KMAX
REAL*8, ALLOCATABLE, DIMENSION(:) :: EX, CS, CP, CD, CF, CG
REAL*8, DIMENSION(1:21), PARAMETER :: H =  (/ 0.0D0,     -0.707106781186548D0, &
                                                          0.707106781186548D0, -1.22474487139159D0,  &
                                                          0.0D0,                1.22474487139159D0,  &
                                                         -1.65068012388578D0,  -0.524647623275290D0, &
                                                          0.524647623275290D0,  1.65068012388578D0,  &
                                                         -2.02018287045609D0,  -0.958572464613819D0, &
                                                          0.0D0,                0.958572464613819D0, &
                                                          2.02018287045609D0,  -2.350604973674D0,    &
                                                         -1.335849074014D0,    -0.436077411928D0,    &
                                                          0.436077411928D0,     1.335849074014D0,    &
                                                          2.350604973674D0 /)
REAL*8, DIMENSION(1:21), PARAMETER :: W =  (/ 1.77245385090552D0,   0.8862269254528D0,    &
                                                          0.8862269254528D0,  0.2954089751509D0,  &
                                                          1.181635900604D0,   0.2954089751509D0,  &
                                                          8.131283544725D-02, 8.049140900055D-01, &
                                                          8.049140900055D-01, 8.131283544725D-02, &
                                                          1.995324205905D-02, 3.936193231522D-01, &
                                                          9.453087204829D-01, 3.936193231522D-01, &
                                                          1.995324205905D-02, 4.530009905509D-03, &
                                                          1.570673203229D-01, 7.246295952244D-01, &
                                                          7.246295952244D-01, 1.570673203229D-01, &
                                                          4.530009905509D-03 /)

INTEGER :: NQMINT
REAL*8 :: CUT_OFF = 990000.0D0, CUT_ON = 980000.0D0
REAL*8 :: CUTOFFB = 0.0D0, CUTONB = 0.0D0, GAMMAB = 1.0D0, R2OFFB = 0.0D0, R2ONB = 0.0D0, &
          T0FAC   = 0.0D0, T2FAC  = 0.0D0, T4FAC  = 0.0D0, T6FAC  = 0.0D0
REAL*8, ALLOCATABLE, DIMENSION(:)   :: DXE1BA, DYE1BA, DZE1BA
REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DXE1BB, DYE1BB, DZE1BB
INTEGER            :: IEXTRP    = 15,         NDIIS = 19,     NDIISP = 20 ! = NDIIS + 1
LOGICAL            :: QFORCEUHF = .FALSE., QSCFBOMB = .TRUE.
REAL*8 :: ACURCY    = 1.0D-8,  DAMPF = 0.0D0, DAMP0 = 0.0D0, SHIFTO = 0.0D0, SHIFTV = 0.0D0
INTEGER                                         :: MATNUM
INTEGER,            ALLOCATABLE, DIMENSION(:)   :: MATIND
REAL*8, ALLOCATABLE, DIMENSION(:,:) :: BCOEFF

CONTAINS

   FUNCTION NORMALIZE( VECTOR )
   REAL*8, DIMENSION(:), INTENT(IN) :: VECTOR
   REAL*8, DIMENSION(1:SIZE(VECTOR)) :: NORMALIZE
   INTEGER :: N
   REAL*8 :: FACT

   N = SIZE( VECTOR )
   FACT = SQRT( DOT_PRODUCT( VECTOR(1:N), VECTOR(1:N) ) )
   IF ( FACT > 0.0D0 ) NORMALIZE = VECTOR / FACT
   END FUNCTION NORMALIZE


   FUNCTION RANDOM_VECTOR( N )
   INTEGER, INTENT(IN) :: N
   REAL*8, DIMENSION(1:N) :: RANDOM_VECTOR
   INTEGER :: I

   CALL RANDOM_SEED
   DO I = 1, N
       CALL RANDOM_NUMBER( RANDOM_VECTOR(I) )
   END DO
   END FUNCTION RANDOM_VECTOR


   SUBROUTINE SYMMETRIC_UPPER ( MATRIX, EIGENVALUES, EIGENVECTORS )
   REAL*8, DIMENSION(:), INTENT(INOUT) :: MATRIX
   REAL*8, DIMENSION(:), INTENT(OUT)   :: EIGENVALUES
   REAL*8, DIMENSION(:,:), INTENT(OUT), OPTIONAL :: EIGENVECTORS
   INTEGER :: IFAIL, N
   LOGICAL :: QOK
   REAL*8, ALLOCATABLE, DIMENSION(:)   :: WORK
   REAL*8,              DIMENSION(1:1) :: EVTEMP

   N = SIZE ( EIGENVALUES )
   QOK = ( ( N * ( N + 1 ) ) / 2 ) == SIZE ( MATRIX )
   IF ( PRESENT ( EIGENVECTORS ) ) THEN
      QOK = QOK .AND. ( SIZE ( EIGENVECTORS, 1 ) == N ) &
                .AND. ( SIZE ( EIGENVECTORS, 2 ) == N )
   END IF
   IF ( .NOT. QOK ) THEN
       WRITE( *, * ) ">> Array dimension error"
       STOP
   END IF
   ALLOCATE ( WORK(1:3*N) )
   IFAIL = 0
   IF ( PRESENT ( EIGENVECTORS ) ) THEN
      CALL DSPEV ( "V", "U", N, MATRIX, EIGENVALUES, EIGENVECTORS, N, WORK, IFAIL )
   ELSE
      CALL DSPEV ( "N", "U", N, MATRIX, EIGENVALUES, EVTEMP,       N, WORK, IFAIL )
   END IF
   IF ( IFAIL /= 0 ) THEN
       WRITE( *, * ) ">> Diagonalization error"
       STOP
   END IF
   DEALLOCATE ( WORK )
   END SUBROUTINE SYMMETRIC_UPPER
   

   SUBROUTINE LINEAR_EQUATIONS_SOLVE ( N, A, IA, B, IERR )
   INTEGER, INTENT(IN)  :: IA, N
   INTEGER, INTENT(OUT) :: IERR
   REAL*8, DIMENSION(1:IA,1:N), INTENT(INOUT) :: A
   REAL*8, DIMENSION(1:N),      INTENT(INOUT) :: B
   INTEGER :: INFO
   INTEGER, DIMENSION(1:N) :: IPIV

   INFO = 0
   CALL DGESV ( N, 1, A(1:IA,1:N), IA, IPIV(1:N), B(1:N), N, INFO )
   IERR = INFO
   END SUBROUTINE LINEAR_EQUATIONS_SOLVE


   SUBROUTINE GAUSSIAN_INITIALIZE
   NSHELL   = 0 
   IF( ALLOCATED ( KSTART  ) ) DEALLOCATE( KSTART  )
   IF( ALLOCATED ( KATOM   ) ) DEALLOCATE( KATOM   )
   IF( ALLOCATED ( KTYPE   ) ) DEALLOCATE( KTYPE   )
   IF( ALLOCATED ( KNG     ) ) DEALLOCATE( KNG     )
   IF( ALLOCATED ( KLOC    ) ) DEALLOCATE( KLOC    )
   IF( ALLOCATED ( KMIN    ) ) DEALLOCATE( KMIN    )
   IF( ALLOCATED ( KMAX    ) ) DEALLOCATE( KMAX    )
   IF( ALLOCATED ( EX      ) ) DEALLOCATE( EX      )
   IF( ALLOCATED ( CS      ) ) DEALLOCATE( CS      )
   IF( ALLOCATED ( CP      ) ) DEALLOCATE( CP      )
   IF( ALLOCATED ( CD      ) ) DEALLOCATE( CD      )
   IF( ALLOCATED ( CF      ) ) DEALLOCATE( CF      )
   IF( ALLOCATED ( CG      ) ) DEALLOCATE( CG      )
   END SUBROUTINE GAUSSIAN_INITIALIZE


   SUBROUTINE GAUSSIAN_OVERLAP( OVERLAP )
   REAL*8, DIMENSION(:), INTENT(OUT) :: OVERLAP
   INTEGER :: NI, NJ
   REAL*8 :: XINT, YINT, ZINT, T, X0, Y0, Z0, XI, YI, ZI, XJ, YJ, ZJ
   CALL GET_OVERLAP( OVERLAP, 1, NSHELL )
   CONTAINS
      SUBROUTINE GET_OVERLAP( OVERLAP, NSTART, NSTOP )
      INTEGER, INTENT(IN) :: NSTART, NSTOP
      REAL*8, DIMENSION(:), INTENT(OUT) :: OVERLAP
      INTEGER :: I, IG, II, IJ, IN, ITOL, I1, I2, J, JG, JGMAX, JJ, JN, &
                 J1, J2, LI, LIT, LJ, LJT, LOCI, LOCJ, MAX, MAXI,       &
                 MAXJ, MINI, MINJ, N, NN, NX, NY, NZ
      LOGICAL :: DOUBLE, IANDJ
      REAL*8 :: AA, AA1, AI, AJ, ARRI, AX, AXI, AY, AYI, AZ, AZI,   &
                CDI, CDJ, CPI, CPJ, CSI, CSJ, DUM, DUM1, DUM2, FAC, &
                RR, TOL, T1, T2, YZ
      INTEGER, DIMENSION(1:10) :: IX, IY, IZ, JX, JY, JZ
      INTEGER, DIMENSION(1:36) :: IJX, IJY, IJZ
      REAL*8, DIMENSION(1:3)  :: RI, RJ
      REAL*8, DIMENSION(1:36) :: DIJ, S
      REAL*8, DIMENSION(1:27) :: XIN, YIN, ZIN
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0/
      DATA IX / 1, 4, 1, 1, 7, 1, 1, 4, 4, 1/
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1/
      DATA IY / 1, 1, 4, 1, 1, 7, 1, 4, 1, 4/
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1/
      DATA IZ / 1, 1, 1, 4, 1, 1, 7, 1, 4, 4/

      ITOL = 15
      TOL  = RLN10 * REAL( ITOL, 8 )
      OVERLAP = 0.0D0
      DO 9000 II = NSTART,NSTOP
      RI = ATMCRD(1:3,KATOM(II)) * ANGSTROMS_TO_BOHRS
      XI = RI(1)
      YI = RI(2)
      ZI = RI(3)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
      DO 8000 JJ = 1,II
      RJ = ATMCRD(1:3,KATOM(JJ)) * ANGSTROMS_TO_BOHRS
      XJ = RJ(1)
      YJ = RJ(2)
      ZJ = RJ(3)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      IANDJ=II.EQ.JJ
      IJ=0
      MAX=MAXJ
      DO I=MINI,MAXI
         NX=IX(I)
         NY=IY(I)
         NZ=IZ(I)
         IF(IANDJ) MAX=I
         DO J=MINJ,MAX
            IJ=IJ+1
            IJX(IJ)=NX+JX(J)
            IJY(IJ)=NY+JY(J)
            IJZ(IJ)=NZ+JZ(J)
         END DO
      END DO
      DO 60 I=1,IJ
      S(I)=0.0D0
      60 CONTINUE
      JGMAX=J2
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      IF(IANDJ) JGMAX=IG
      DO 6000 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AA1=1.0D0/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC=DEXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
      DOUBLE=IANDJ.AND.IG.NE.JG
      MAX=MAXJ
      NN=0
      DO I=MINI,MAXI
         GO TO (70,80,110,110,90,110,110,100,110,110),I
         70 DUM1=CSI*FAC
         GO TO 110
         80 DUM1=CPI*FAC
         GO TO 110
         90 DUM1=CDI*FAC
         GO TO 110
         100 DUM1=DUM1*SQRT3
         110 IF(IANDJ) MAX=I
         DO J=MINJ,MAX
            GO TO (120,130,160,160,140,160,160,150,160,160),J
            120 DUM2=DUM1*CSJ
            IF(.NOT.DOUBLE) GO TO 160
            IF(I.GT.1) GO TO 125
            DUM2=DUM2+DUM2
            GO TO 160
            125 DUM2=DUM2+CSI*CPJ*FAC
            GO TO 160
            130 DUM2=DUM1*CPJ
            IF(DOUBLE) DUM2=DUM2+DUM2
            GO TO 160
            140 DUM2=DUM1*CDJ
            IF(DOUBLE) DUM2=DUM2+DUM2
            GO TO 160
            150 DUM2=DUM2*SQRT3
            160 NN=NN+1
            DIJ(NN)=DUM2
         END DO
      END DO
      T= SQRT(AA1)
      T1=-2.0D0*AJ*AJ*T
      T2=-0.5D0*T
      X0=AX
      Y0=AY
      Z0=AZ
      IN=-3
      DO I=1,LIT
         IN=IN+3
         NI=I
         DO J=1,LJT
            JN=IN+J
            NJ=J
            CALL STVXYZ
            XIN(JN)=XINT*T
            YIN(JN)=YINT*T
            ZIN(JN)=ZINT*T
            NJ=J+2
            CALL STVXYZ
            XIN(JN+9)=XINT*T1
            YIN(JN+9)=YINT*T1
            ZIN(JN+9)=ZINT*T1
            NJ=J-2
            IF(NJ.GT.0) GO TO 320
            XINT=0.0D0
            YINT=0.0D0
            ZINT=0.0D0
            GO TO 330
            320 CALL STVXYZ
            330 N=(J-1)*(J-2)
            DUM= FLOAT(N)*T2
            XIN(JN+18)=XINT*DUM
            YIN(JN+18)=YINT*DUM
            ZIN(JN+18)=ZINT*DUM
         END DO
      END DO
      DO 350 I=1,IJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      YZ=YIN(NY)*ZIN(NZ)
      DUM  = YZ*XIN(NX)
      DUM1 = (XIN(NX+9)+XIN(NX+18))*YZ+(YIN(NY+9)+YIN(NY+18))*XIN(NX)*ZIN(NZ)+(ZIN(NZ+9)+ZIN(NZ+18))*XIN(NX)*YIN(NY)
      S(I)=S(I)+DIJ(I)*DUM
       350 CONTINUE
      6000 CONTINUE
      7000 CONTINUE
      MAX=MAXJ
      NN=0
      DO I=MINI,MAXI
         LI=LOCI+I
         IN=BFINDEX(LI)
         IF(IANDJ) MAX=I
         DO J=MINJ,MAX
            LJ=LOCJ+J
            JN=LJ+IN
            NN=NN+1
            OVERLAP(JN) = S(NN)
         END DO
      END DO
         8000 CONTINUE
      9000 CONTINUE
      END SUBROUTINE GET_OVERLAP


      SUBROUTINE STVXYZ
      INTEGER :: I, IMAX, IMIN, NPTS
      REAL*8 :: AX, AY, AZ, BX, BY, BZ, DUM, PTX, PTY, PTZ, PX, PY, PZ
      INTEGER, DIMENSION(1:6) :: MIN, MAX
      DATA MIN /1, 2, 4, 7, 11, 16/
      DATA MAX /1, 3, 6, 10, 15, 21/

      XINT=0.0D0
      YINT=0.0D0
      ZINT=0.0D0
      NPTS=(NI+NJ-2)/2+1
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO I=IMIN,IMAX
         DUM=W(I)
         PX=DUM
         PY=DUM
         PZ=DUM
         DUM=H(I)*T
         PTX=DUM+X0
         PTY=DUM+Y0
         PTZ=DUM+Z0
         AX=PTX-XI
         AY=PTY-YI
         AZ=PTZ-ZI
         BX=PTX-XJ
         BY=PTY-YJ
         BZ=PTZ-ZJ
         GO TO (5,4,3,2,1),NI
         1 PX=PX*AX
         PY=PY*AY
         PZ=PZ*AZ
         2 PX=PX*AX
         PY=PY*AY
         PZ=PZ*AZ
         3 PX=PX*AX
         PY=PY*AY
         PZ=PZ*AZ
         4 PX=PX*AX
         PY=PY*AY
         PZ=PZ*AZ
         5 GO TO (12,11,10,9,8,7,6),NJ
         6 PX=PX*BX
         PY=PY*BY
         PZ=PZ*BZ
         7 PX=PX*BX
         PY=PY*BY
         PZ=PZ*BZ
         8 PX=PX*BX
         PY=PY*BY
         PZ=PZ*BZ
         9 PX=PX*BX
         PY=PY*BY
         PZ=PZ*BZ
         10 PX=PX*BX
         PY=PY*BY
         PZ=PZ*BZ
         11 PX=PX*BX
         PY=PY*BY
         PZ=PZ*BZ
         12 CONTINUE
         XINT=XINT+PX
         YINT=YINT+PY
         ZINT=ZINT+PZ
      END DO
      END SUBROUTINE STVXYZ
   END SUBROUTINE GAUSSIAN_OVERLAP


   SUBROUTINE GAUSSIAN_OVERLAP_DERIVATIVES( SX, SY, SZ )
   REAL*8, DIMENSION(:), INTENT(OUT) ::   SX,   SY,   SZ
   INTEGER :: NI, NJ
   REAL*8 :: T, XI, XINT, XJ, X0, YI, YINT, YJ, Y0, ZI, ZINT, ZJ, Z0
   CALL GET_OVERLAP_DERIVATIVES( SX, SY, SZ, 1, NSHELL )
   CONTAINS
      SUBROUTINE GET_OVERLAP_DERIVATIVES( SX, SY, SZ, NSTART, NSTOP )
      INTEGER, INTENT(IN) :: NSTART, NSTOP
      REAL*8, DIMENSION(:), INTENT(OUT) :: SX, SY, SZ
      INTEGER :: I, IG, II, IJ, IN, ITOL, I1, I2, J, JG, JJ, JN, J1, J2, LIT, LOCI, LJT, LOCJ, &
                 MAXI, MAXJ, MINI, MINJ, N, NN, NX, NY, NZ, N0
      LOGICAL :: IANDJ
      REAL*8 :: AA, AA1, AI, AJ, ARRI, AX, AXI, AY, AYI, AZ, AZI, CDI, CDJ, CFI, CPI, CPJ, CSI, CSJ, &
                DUM, DUM1, DUM2, FAC, RR, TOL
      INTEGER, DIMENSION(1:78) :: IJG, IJX, IJY, IJZ
      LOGICAL, DIMENSION(1:20) :: ISKIP
      REAL*8, DIMENSION(1:3)  :: RI, RJ
      REAL*8, DIMENSION(1:78) :: DIJ, S
      REAL*8, DIMENSION(1:36) :: XIN, YIN, ZIN
      INTEGER, DIMENSION(1:40) :: INDX
      INTEGER, DIMENSION(1:20) :: IX, IY, IZ
      INTEGER, DIMENSION(1:10) :: JX, JY, JZ
      DATA INDX / 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0, &
                  -0, 1, 2, 3,-0,-0,-0,-0,-0,-0,4, 5, 6, 7, 8, 9,10,11,12,13/
      DATA IX   / 1, 4, 1, 1, 7, 1, 1, 4, 4, 1,10, 1, 1, 7, 7, 4, 1, 4, 1, 4/
      DATA IY   / 1, 1, 4, 1, 1, 7, 1, 4, 1, 4,1,10, 1, 4, 1, 7, 7, 1, 4, 4/
      DATA IZ   / 1, 1, 1, 4, 1, 1, 7, 1, 4, 4,1, 1,10, 1, 4, 1, 4, 7, 7, 4/
      DATA JX   / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0/
      DATA JY   / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1/
      DATA JZ   / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1/

      SX = 0.0D0 ; SY = 0.0D0 ; SZ = 0.0D0
      ITOL = 15
      TOL  = RLN10 * REAL( ITOL, 8 )
      DO 9000 II=NSTART,NSTOP
      RI = ATMCRD(1:3,KATOM(II)) * ANGSTROMS_TO_BOHRS
      XI = RI(1)
      YI = RI(2)
      ZI = RI(3)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)+1
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
      ISKIP(1:20)=.TRUE.
      DO 20 I=MINI,MAXI
      GOTO (11,13,20,20,15,20,20,20,20,20),I
   11 ISKIP(2:4)=.FALSE.
      GOTO 20
   13 ISKIP(5:10)=.FALSE.
      ISKIP(1)=.FALSE.
      GOTO 20
   15 ISKIP(2:4)=.FALSE.
      ISKIP(11:20)=.FALSE.
   20 CONTINUE
      DO 8000 JJ=1,II
      RJ = ATMCRD(1:3,KATOM(JJ)) * ANGSTROMS_TO_BOHRS
      XJ = RJ(1)
      YJ = RJ(2)
      ZJ = RJ(3)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      IANDJ=II.EQ.JJ
      N0=0
      IF(LIT.EQ.4) N0=20
      IJ=0
      DO 60 I=1,20
      IF(ISKIP(I)) GOTO 60
      IN=INDX(I+N0)
      NX=IX(I)
      NY=IY(I)
      NZ=IZ(I)
      DO 50 J=MINJ,MAXJ
      IJ=IJ+1
      IJX(IJ)=NX+JX(J)
      IJY(IJ)=NY+JY(J)
      IJZ(IJ)=NZ+JZ(J)
      IJG(IJ)=IN+13*(J-MINJ)
   50 CONTINUE
   60 CONTINUE
      DO I=1,IJ
         N=IJG(I)
         S(N)=0.0D0
      END DO
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      DUM=AI+AI
      CSI=CP(IG)
      CPI=CS(IG)*DUM
      IF(LIT.EQ.4) CPI=CD(IG)
      CDI=CP(IG)*DUM
      CFI=CD(IG)*DUM
      DO 6000 JG=J1,J2
      AJ=EX(JG)
      AA =AI+AJ
      AA1=1.0D0/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GOTO 6000
      FAC=DEXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
      NN=0
      DO 180 I=1,20
      IF(ISKIP(I)) GOTO 180
      GOTO ( 70, 80,110,110, 90,110,110,110,110,110, 100,110,110,110,110,110,110,110,110,110),I
   70 DUM1=CSI*FAC
      GOTO 110
   80 DUM1=CPI*FAC
      GOTO 110
   90 DUM1=CDI*FAC
      GOTO 110
  100 DUM1=CFI*FAC
  110 CONTINUE
      DO J=MINJ,MAXJ
         GOTO (120,130,160,160,140,160,160,150,160,160),J
  120    DUM2=DUM1*CSJ
         GOTO 160
  130    DUM2=DUM1*CPJ
         GOTO 160
  140    DUM2=DUM1*CDJ
         GOTO 160
  150    DUM2=DUM2*SQRT3
  160    NN=NN+1
         DIJ(NN)=DUM2
      END DO
  180 CONTINUE
      T=SQRT(AA1)
      X0=AX
      Y0=AY
      Z0=AZ
      IN=-3
      DO I=1,LIT
         IN=IN+3
         NI=I
         DO J=1,LJT
            JN=IN+J
            NJ=J
            CALL DERXYZ
            XIN(JN)=XINT*T
            YIN(JN)=YINT*T
            ZIN(JN)=ZINT*T
         END DO
      END DO
      DO 350 I=1,IJ
      N=IJG(I)
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      S(N)=S(N)+DIJ(I)*XIN(NX)*YIN(NY)*ZIN(NZ)
  350 CONTINUE
 6000 CONTINUE
 7000 CONTINUE
      NN=0
      N=1
      DO 7301 J=MINJ,MAXJ
      IF(MINI.GT.1) GOTO 7100
      NN=NN+1
      XIN(NN)= S(N+ 1)
      YIN(NN)= S(N+ 2)
      ZIN(NN)= S(N+ 3)
      IF(MAXI.EQ.1) GOTO 7300
 7100 IF(MINI.GT.2) GOTO 7200
      NN=NN+1
      XIN(NN)=(S(N+ 4)-S(N   ))
      YIN(NN)= S(N+ 7)
      ZIN(NN)= S(N+ 8)
      NN=NN+1
      XIN(NN)= S(N+ 7)
      YIN(NN)=(S(N+ 5)-S(N   ))
      ZIN(NN)= S(N+ 9)
      NN=NN+1
      XIN(NN)= S(N+ 8)
      YIN(NN)= S(N+ 9)
      ZIN(NN)=(S(N+ 6)-S(N   ))
      IF(MAXI.EQ.4) GOTO 7300
 7200 CONTINUE
      NN=NN+1
      XIN(NN)=(S(N+ 3)-S(N   )-S(N   ))
      YIN(NN)= S(N+ 6)
      ZIN(NN)= S(N+ 7)
      NN=NN+1
      XIN(NN)= S(N+ 8)
      YIN(NN)=(S(N+ 4)-S(N+ 1)-S(N+ 1))
      ZIN(NN)= S(N+ 9)
      NN=NN+1
      XIN(NN)= S(N+10)
      YIN(NN)= S(N+11)
      ZIN(NN)=(S(N+ 5)-S(N+ 2)-S(N+ 2))
      NN=NN+1
      DUM=SQRT3
      XIN(NN)=DUM*(S(N+ 6)-S(N+ 1))
      YIN(NN)=DUM*(S(N+ 8)-S(N   ))
      ZIN(NN)=DUM* S(N+12)
      NN=NN+1
      XIN(NN)=DUM*(S(N+ 7)-S(N+ 2))
      YIN(NN)=DUM* S(N+12)
      ZIN(NN)=DUM*(S(N+10)-S(N   ))
      NN=NN+1
      XIN(NN)=DUM* S(N+12)
      YIN(NN)=DUM*(S(N+ 9)-S(N+ 2))
      ZIN(NN)=DUM*(S(N+11)-S(N+ 1))
 7300 N=N+13
 7301 CONTINUE
      N=0
      DO 7500 J=MINJ,MAXJ
      JN=LOCJ+J
      DO 7400 I=MINI,MAXI
      N=N+1
      IN=LOCI+I
      IF(JN.GT.IN) GOTO 7400
      NN=BFINDEX(IN)+JN
      SX(NN)=XIN(N)
      SY(NN)=YIN(N)
      SZ(NN)=ZIN(N)
 7400 CONTINUE
 7500 CONTINUE
 8000 CONTINUE
 9000 CONTINUE
      END SUBROUTINE GET_OVERLAP_DERIVATIVES

      SUBROUTINE DERXYZ
      INTEGER ::  I, IMIN, IMAX, NPTS
      REAL*8 :: AX, AY, AZ, BX, BY, BZ, DUM, PTX, PTY, PTZ, PX, PY, PZ
      INTEGER, DIMENSION(1:6) :: MAX, MIN
      DATA MAX / 1, 3, 6, 10, 15, 21/
      DATA MIN / 1, 2, 4,  7, 11, 16/

      XINT=0.0D0
      YINT=0.0D0
      ZINT=0.0D0
      NPTS=(NI+NJ-2)/2+1
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO 13 I=IMIN,IMAX
      DUM=W(I)
      PX=DUM
      PY=DUM
      PZ=DUM
      DUM=H(I)*T
      PTX=DUM+X0
      PTY=DUM+Y0
      PTZ=DUM+Z0
      AX=PTX-XI
      AY=PTY-YI
      AZ=PTZ-ZI
      BX=PTX-XJ
      BY=PTY-YJ
      BZ=PTZ-ZJ
      GOTO (5,4,3,2,1),NI
    1 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    2 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    3 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    4 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    5 GOTO (12,11,10,9,8,7,6),NJ
    6 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
    7 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
    8 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
    9 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   10 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   11 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   12 CONTINUE
      XINT=XINT+PX
      YINT=YINT+PY
      ZINT=ZINT+PZ
   13 CONTINUE
      END SUBROUTINE DERXYZ
   END SUBROUTINE GAUSSIAN_OVERLAP_DERIVATIVES

   SUBROUTINE GAUSSIAN_SETUP( NHEAVY, NLIGHT )
   INTEGER, INTENT(IN) :: NHEAVY, NLIGHT
   CALL GAUSSIAN_INITIALIZE
   CALL GET_BASIS_SET   
   CALL NORMALIZE_BASIS_FUNCTIONS
   CONTAINS
      SUBROUTINE GET_BASIS_SET
      INTEGER :: I, IGAUSS, IQM, LOC, NGAUSS, NI
      REAL*8 :: EE, FACP, FACS, SCALE, ZETAP, ZETAS
      REAL*8, DIMENSION(1:6) :: CCP, CCS, EXP, EXS
      NSHELL = ( 2 * NHEAVY ) + NLIGHT
      NGAUSS = 6 * NSHELL
      ALLOCATE ( KSTART(1:NSHELL), KATOM(1:NSHELL), KTYPE(1:NSHELL), KNG(1:NSHELL), &
                 KLOC(1:NSHELL),   KMIN(1:NSHELL),  KMAX(1:NSHELL) )
      ALLOCATE ( EX(1:NGAUSS), CS(1:NGAUSS), CP(1:NGAUSS), CD(1:NGAUSS), CF(1:NGAUSS), CG(1:NGAUSS) )
      CD     = 0.0D0
      CF     = 0.0D0
      CG     = 0.0D0
      CP     = 0.0D0
      CS     = 0.0D0
      EX     = 0.0D0
      LOC    = 0
      NGAUSS = 0
      NSHELL = 0
      DO IQM = 1,NAQM
         NI    = ANINT( ATMCHG(IQM) )
         IF( NATORB(NI) > 0 ) THEN
            ZETAS = ZS(NI)
            ZETAP = ZP(NI)
            SELECT CASE ( NI )
               ! . 1s orbitals.
               CASE ( 1:2 )
                  ! . 1s expansion.
                  EXS(1) =  2.310303149D+01
                  CCS(1) =  9.163596280D-03
                  EXS(2) =  4.235915534D+00
                  CCS(2) =  4.936149294D-02
                  EXS(3) =  1.185056519D+00
                  CCS(3) =  1.685383049D-01
                  EXS(4) =  4.070988982D-01
                  CCS(4) =  3.705627997D-01
                  EXS(5) =  1.580884151D-01
                  CCS(5) =  4.164915298D-01
                  EXS(6) =  6.510953954D-02
                  CCS(6) =  1.303340841D-01
               ! . 2s and 2p orbitals.
               CASE ( 3:10 )
                  ! . 2s expansion.
                  EXS(1) =  2.768496241D+01
                  CCS(1) = -4.151277819D-03
                  EXS(2) =  5.077140627D+00
                  CCS(2) = -2.067024148D-02
                  EXS(3) =  1.426786050D+00
                  CCS(3) = -5.150303337D-02
                  EXS(4) =  2.040335729D-01
                  CCS(4) =  3.346271174D-01
                  EXS(5) =  9.260298399D-02
                  CCS(5) =  5.621061301D-01
                  EXS(6) =  4.416183978D-02
                  CCS(6) =  1.712994697D-01
                  ! . 2p expansion.
                  EXP(1) =  5.868285913D+00
                  CCP(1) =  7.924233646D-03
                  EXP(2) =  1.530329631D+00
                  CCP(2) =  5.144104825D-02
                  EXP(3) =  5.475665231D-01
                  CCP(3) =  1.898400060D-01
                  EXP(4) =  2.288932733D-01
                  CCP(4) =  4.049863191D-01
                  EXP(5) =  1.046655969D-01
                  CCP(5) =  4.012362861D-01
                  EXP(6) =  4.948220127D-02
                  CCP(6) =  1.051855189D-01
               ! . 3s and 3p orbitals.
               CASE ( 11:18 )
                  ! . 3s expansion.
                  EXS(1) =  3.273031938D+00
                  CCS(1) = -6.775596947D-03
                  EXS(2) =  9.200611311D-01
                  CCS(2) = -5.639325779D-02
                  EXS(3) =  3.593349765D-01
                  CCS(3) = -1.587856086D-01
                  EXS(4) =  8.636686991D-02
                  CCS(4) =  5.534527651D-01
                  EXS(5) =  4.797373812D-02
                  CCS(5) =  5.015351020D-01
                  EXS(6) =  2.724741144D-02
                  CCS(6) =  7.223633674D-02
                  ! . 3p expansion.
                  EXP(1) =  5.077973607D+00
                  CCP(1) = -3.329929840D-03
                  EXP(2) =  1.340786940D+00
                  CCP(2) = -1.419488340D-02
                  EXP(3) =  2.248434849D-01
                  CCP(3) =  1.639395770D-01
                  EXP(4) =  1.131741848D-01
                  CCP(4) =  4.485358256D-01
                  EXP(5) =  6.076408893D-02
                  CCP(5) =  3.908813050D-01
                  EXP(6) =  3.315424265D-02
                  CCP(6) =  7.411456232D-02
               ! . 4s and 4p orbitals.
               CASE ( 19:36 )
                  ! . 4s expansion.
                  EXS(1) =  3.232838646D+00
                  CCS(1) =  1.374817488D-03
                  EXS(2) =  3.605788802D-01
                  CCS(2) = -8.666390043D-02
                  EXS(3) =  1.717905487D-01
                  CCS(3) = -3.130627309D-01
                  EXS(4) =  5.277666487D-02
                  CCS(4) =  7.812787397D-01
                  EXS(5) =  3.163400284D-02
                  CCS(5) =  4.389247988D-01
                  EXS(6) =  1.874093091D-02
                  CCS(6) =  2.487178756D-02
                  ! . 4p expansion.
                  EXP(1) =  2.389722618D+00
                  CCP(1) = -1.665913575D-03
                  EXP(2) =  7.960947826D-01
                  CCP(2) = -1.657464971D-02
                  EXP(3) =  3.415541380D-01
                  CCP(3) = -5.958513378D-02
                  EXP(4) =  8.847434525D-02
                  CCP(4) =  4.053115554D-01
                  EXP(5) =  4.958248334D-02
                  CCP(5) =  5.433958189D-01
                  EXP(6) =  2.816929784D-02
                  CCP(6) =  1.204970491D-01
               ! . 5s and 5p orbitals.
               CASE ( 37:54 )
                  ! . 5s expansion.
                  EXS(1) =  1.410128298D+00
                  CCS(1) =  2.695439582D-03
                  EXS(2) =  5.077878915D-01
                  CCS(2) =  1.850157487D-02
                  EXS(3) =  1.847926858D-01
                  CCS(3) = -9.588628125D-02
                  EXS(4) =  1.061070594D-01
                  CCS(4) = -5.200673560D-01
                  EXS(5) =  3.669584901D-02
                  CCS(5) =  1.087619490D+00
                  EXS(6) =  2.213558430D-02
                  CCS(6) =  3.103964343D-01
                  ! . 5p expansion.
                  EXP(1) =  3.778623374D+00
                  CCP(1) =  1.163246387D-04
                  EXP(2) =  3.499121109D-01
                  CCP(2) = -2.920771322D-02
                  EXP(3) =  1.683175469D-01
                  CCP(3) = -1.381051233D-01
                  EXP(4) =  5.404070736D-02
                  CCP(4) =  5.706134877D-01
                  EXP(5) =  3.328911801D-02
                  CCP(5) =  4.768808140D-01
                  EXP(6) =  2.063815019D-02
                  CCP(6) =  6.021665516D-02
               ! . Higher shells.
               CASE DEFAULT
                  WRITE( *, * ) ">> Gaussian expansions unavailable for atomic number > 54"
                  STOP
               END SELECT
               NSHELL = NSHELL + 1
               KATOM  (NSHELL) = IQM
               KLOC   (NSHELL) = LOC + 1
               KMAX   (NSHELL) = 1
               KMIN   (NSHELL) = 1
               KNG    (NSHELL) = 6
               KSTART (NSHELL) = NGAUSS + 1
               KTYPE  (NSHELL) = 1
               NGAUSS = NGAUSS + 6
               IGAUSS = KSTART(NSHELL) - 1
               SCALE  = ZETAS * ZETAS
               CS(IGAUSS+1:IGAUSS+6) = CCS(1:6)
               EX(IGAUSS+1:IGAUSS+6) = SCALE * EXS(1:6)
               LOC = LOC + 1
               DO I = 1,6
                  EE   = EX(IGAUSS+I) + EX(IGAUSS+I)
                  FACS = PI32 / ( EE * SQRT ( EE ) )
                  CS(IGAUSS+I) = CS(IGAUSS+I) / SQRT ( FACS )
               END DO
               IF( NATORB(NI) > 1 ) THEN
                  NSHELL = NSHELL + 1
                  KATOM  (NSHELL) = IQM
                  KLOC   (NSHELL) = LOC + 1
                  KMAX   (NSHELL) = 4
                  KMIN   (NSHELL) = 2
                  KNG    (NSHELL) = 6
                  KSTART (NSHELL) = NGAUSS + 1
                  KTYPE  (NSHELL) = 2
                  NGAUSS = NGAUSS + 6
                  IGAUSS = KSTART(NSHELL) - 1
                  SCALE  = ZETAP * ZETAP
                  CP(IGAUSS+1:IGAUSS+6) = CCP(1:6)
                  EX(IGAUSS+1:IGAUSS+6) = SCALE * EXP(1:6)
                  LOC = LOC + 3
                  DO I = 1,6
                     EE   = EX(IGAUSS+I) + EX(IGAUSS+I)
                     FACS = PI32 / ( EE * SQRT ( EE ) )
                     FACP = 0.5D0 * FACS / EE
                     CP(IGAUSS+I) = CP(IGAUSS+I) / SQRT ( FACP )
                  END DO
               END IF
         END IF
      END DO
      END SUBROUTINE GET_BASIS_SET

      SUBROUTINE NORMALIZE_BASIS_FUNCTIONS
      INTEGER :: I, IG, II, I1, I2, JG, LIT, LOCI, MAXI, MINI
      REAL*8 :: DNRM, DNRMI, DUM, DUMD, DUMF, DUMG, DUMP, DUMS, EE, &
                FNRM, FNRMI, GNRM, GNRMI, PNRM, PNRMI, SNRM, SNRMI

      DO II = 1,NSHELL
         I1 = KSTART(II)
         I2 = I1 + KNG(II) - 1
         LIT  = KTYPE(II)
         MINI = KMIN(II)
         MAXI = KMAX(II)
         LOCI = KLOC(II) - MINI
         SNRM = 0.0D0
         PNRM = 0.0D0
         DNRM = 0.0D0
         FNRM = 0.0D0
         GNRM = 0.0D0
         DO I = MINI,MAXI
            SELECT CASE ( I )
            CASE (  1 ) ; SNRM = 1.0D0
            CASE (  2 ) ; PNRM = 1.0D0
            CASE (  5 ) ; DNRM = 1.0D0
            CASE ( 11 ) ; FNRM = 1.0D0
            CASE ( 21 ) ; GNRM = 1.0D0
            END SELECT
         END DO
         SNRMI = 0.0D0
         PNRMI = 0.0D0
         DNRMI = 0.0D0
         FNRMI = 0.0D0
         GNRMI = 0.0D0
         DO IG = I1,I2
            DO JG = I1,IG
               EE   = EX(IG) + EX(JG)
               DUM  = EE * SQRT (EE)
               DUMS =             CS(IG) * CS(JG) / DUM
               DUMP = 0.5D0    * CP(IG) * CP(JG) / ( EE      * DUM )
               DUMD = 0.75D0   * CD(IG) * CD(JG) / ( EE ** 2 * DUM )
               DUMF = 1.875D0  * CF(IG) * CF(JG) / ( EE ** 3 * DUM )
               DUMG = 6.5625D0 * CG(IG) * CG(JG) / ( EE ** 4 * DUM )
               IF ( JG /= IG ) THEN
                  DUMS = DUMS + DUMS
                  DUMP = DUMP + DUMP
                  DUMD = DUMD + DUMD
                  DUMF = DUMF + DUMF
                  DUMG = DUMG + DUMG
               END IF
               SNRMI = SNRMI + DUMS
               PNRMI = PNRMI + DUMP
               DNRMI = DNRMI + DUMD
               FNRMI = FNRMI + DUMF
               GNRMI = GNRMI + DUMG
            END DO
         END DO
         IF ( SNRMI > 1.0D-10 ) SNRM = 1.0D0 / SQRT ( SNRMI * PI32 )
         IF ( PNRMI > 1.0D-10 ) PNRM = 1.0D0 / SQRT ( PNRMI * PI32 )
         IF ( DNRMI > 1.0D-10 ) DNRM = 1.0D0 / SQRT ( DNRMI * PI32 )
         IF ( FNRMI > 1.0D-10 ) FNRM = 1.0D0 / SQRT ( FNRMI * PI32 )
         IF ( GNRMI > 1.0D-10 ) GNRM = 1.0D0 / SQRT ( GNRMI * PI32 )
         DO IG = I1,I2
            CS(IG) = CS(IG) * SNRM
            CP(IG) = CP(IG) * PNRM
            CD(IG) = CD(IG) * DNRM
            CF(IG) = CF(IG) * FNRM
            CG(IG) = CG(IG) * GNRM
         END DO
      END DO
      END SUBROUTINE NORMALIZE_BASIS_FUNCTIONS
   END SUBROUTINE GAUSSIAN_SETUP


   SUBROUTINE MOPAC_DATA_INITIALIZE
   HAMILTONIAN = "AM1 "
   MULTIP      = 1
   NBASIS      = 0
   NBASTR      = 0
   NALPHA      = 0
   NBETA       = 0
   NELEC       = 0
   N2ELEC      = 0
   TOTCHG      = 0
   ATHEAT      = 0.0D0
   IF( ALLOCATED( BFINDEX ) ) DEALLOCATE( BFINDEX )
   IF( ALLOCATED( BETA    ) ) DEALLOCATE( BETA    )
   IF( ALLOCATED( DENMAT  ) ) DEALLOCATE( DENMAT  )
   IF( ALLOCATED( DENMATA ) ) DEALLOCATE( DENMATA )
   IF( ALLOCATED( DENMATB ) ) DEALLOCATE( DENMATB )
   IF( ALLOCATED( HCORE   ) ) DEALLOCATE( HCORE   )
   IF( ALLOCATED( OVERLAP ) ) DEALLOCATE( OVERLAP )
   IF( ALLOCATED( SETEI   ) ) DEALLOCATE( SETEI   )
   IF( ALLOCATED( USPD    ) ) DEALLOCATE( USPD    )
   IF( ALLOCATED( BFIRST  ) ) DEALLOCATE( BFIRST  )
   IF( ALLOCATED( BLAST   ) ) DEALLOCATE( BLAST   )
   END SUBROUTINE MOPAC_DATA_INITIALIZE


   SUBROUTINE DENSITY_CALCULATE ( DENSITY, MORBS, NOCC, OCCNUM )
   INTEGER, INTENT(IN) :: NOCC
   REAL*8, INTENT(IN) :: OCCNUM
   REAL*8, DIMENSION(1:NBASTR), INTENT(OUT) :: DENSITY
   REAL*8, DIMENSION(1:NBASIS,1:NOCC), INTENT(IN) :: MORBS
   INTEGER :: I, IBASIS, JBASIS

   I = 0
   DO IBASIS = 1, NBASIS
      DO JBASIS = 1, IBASIS
         I = I + 1
         DENSITY(I) = SUM( MORBS(IBASIS,1:NOCC) * MORBS(JBASIS,1:NOCC) )
      END DO
   END DO
   IF ( OCCNUM /= 1.0D0 ) DENSITY = OCCNUM * DENSITY
   END SUBROUTINE DENSITY_CALCULATE


   SUBROUTINE DENSITY_GUESS
   INTEGER :: NELEBAS
   
   IF( ALLOCATED ( DENMAT  ) ) DEALLOCATE( DENMAT  )
   IF( ALLOCATED ( DENMATA ) ) DEALLOCATE( DENMATA )
   IF( ALLOCATED ( DENMATB ) ) DEALLOCATE( DENMATB )
   ALLOCATE( DENMAT(1:NBASTR) )
   IF ( NALPHA /= NBETA ) ALLOCATE( DENMATA(1:NBASTR), DENMATB(1:NBASTR) )
   NELEBAS = NALPHA + NBETA + TOTCHG
   IF( NALPHA == NBETA ) THEN
      CALL GET_DENSITY( DENMAT, NALPHA + NBETA, NELEBAS, 1.0D0 )
   ELSE
      CALL GET_DENSITY( DENMATA, NALPHA, NELEBAS, 0.5D0 )
      CALL GET_DENSITY( DENMATB, NBETA,  NELEBAS, 0.5D0 )
      DENMAT = DENMATA + DENMATB
   END IF
   CONTAINS
      SUBROUTINE GET_DENSITY( DMATRIX, NELETOT, NELEBAS, OCCFAC )
      INTEGER, INTENT(IN) :: NELEBAS, NELETOT
      REAL*8, INTENT(IN) :: OCCFAC
      REAL*8, DIMENSION(:), INTENT(OUT) :: DMATRIX
      REAL*8, PARAMETER :: SMALL = 0.01D0
      INTEGER :: IBASIS, II, IQM, NI, NORBS
      REAL*8 :: FRACTQ, ORBCHG

      DMATRIX(:) = SMALL * RANDOM_VECTOR( NBASTR )
      FRACTQ = ( REAL ( NELETOT, 8 ) - OCCFAC * NELEBAS ) / REAL( NBASIS, 8 )
      DO IQM = 1, NAQM
         NI    = ANINT( ATMCHG(IQM) )
         NORBS = NATORB(NI)
         IF( NORBS > 0 ) THEN
            ORBCHG = ( OCCFAC * CORE(NI) ) / REAL ( NORBS, 8 ) - FRACTQ
            DO IBASIS = BFIRST(IQM), BLAST(IQM)
               II = ( IBASIS * ( IBASIS + 1 ) ) / 2
               DMATRIX(II) = ORBCHG
            END DO
         END IF
      END DO
      END SUBROUTINE GET_DENSITY
   END SUBROUTINE DENSITY_GUESS


   SUBROUTINE MOPAC_SETUP( METHOD, CHARGE, MULTIPLICITY )
   CHARACTER ( LEN = * ), INTENT(IN) :: METHOD
   INTEGER, INTENT(IN), OPTIONAL :: CHARGE, MULTIPLICITY
   INTEGER :: I, IFIRST, ILAST, IQM, J, K, KK, L, LL, &
                 NHEAVY, NI, NLIGHT, NUMORB, N0HEAVY
   REAL*8 :: EAT, ELECS
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: JIND4, KIND4

   CALL MOPAC_DATA_INITIALIZE
   HAMILTONIAN = METHOD
   CALL MOPAC_PARAMETERS_INITIALIZE( HAMILTONIAN )
   IF( PRESENT ( CHARGE       ) ) TOTCHG = CHARGE
   IF( PRESENT ( MULTIPLICITY ) ) MULTIP = MULTIPLICITY
   EAT       = 0.0D0
   ELECS     = 0.0D0
   N0HEAVY   = 0
   NHEAVY    = 0
   NLIGHT    = 0
   DO IQM = 1, NAQM
      NI     = ANINT( ATMCHG(IQM) )
      ATHEAT = ATHEAT + EHEAT(NI)
      EAT    = EAT    + EISOL(NI)
      ELECS  = ELECS  +  CORE(NI)
      NUMORB = NATORB(NI)
      IF( NUMORB == 0 ) THEN
          N0HEAVY = N0HEAVY + 1
      ELSE IF( NUMORB == 1 ) THEN
         NLIGHT = NLIGHT + 1
      ELSE
         NHEAVY = NHEAVY + 1
      END IF
   END DO
   NELEC = - TOTCHG + NINT( ELECS )
   NALPHA = ( NELEC + MULTIP - 1 ) / 2
   NBETA  = ( NELEC - MULTIP + 1 ) / 2
   IF( ( NALPHA + NBETA ) /= NELEC ) THEN
       WRITE( *, * ) ">> Impossible multiplicity for the system"
       STOP
   END IF
   ATHEAT = ( KCAL_TO_KJ * ATHEAT ) - ( EV_TO_KJ * EAT )
   NBASIS = ( 4 * NHEAVY ) + NLIGHT
   NBASTR = ( NBASIS * ( NBASIS + 1 ) ) / 2
   N2ELEC = (50*NHEAVY*(NHEAVY-1)) + (10*NHEAVY*NLIGHT) + (NLIGHT*(NLIGHT-1))/2
   ALLOCATE( BETA(1:NBASIS), BFIRST(1:NAQM), BLAST(1:NAQM), USPD(1:NBASIS) )
   IFIRST = 1
   ILAST  = 0
   DO IQM = 1, NAQM
      NI          = ANINT( ATMCHG(IQM) )
      NUMORB      = NATORB(NI)
      ILAST       = IFIRST + NUMORB - 1
      BFIRST(IQM) = IFIRST
      BLAST(IQM)  = ILAST
      IF( NUMORB > 0 ) THEN
         BETA(IFIRST) = BETAS(NI)
         USPD(IFIRST) = USS(NI)
         IF( IFIRST < ILAST ) THEN
            BETA((IFIRST+1):ILAST) = BETAP(NI)
            USPD((IFIRST+1):ILAST) = UPP(NI)
          END IF
          IFIRST = ILAST + 1
      END IF
   END DO
   CALL GAUSSIAN_SETUP( NHEAVY, NLIGHT )
   ALLOCATE( BFINDEX(1:NBASIS) )
   DO I = 1,NBASIS
      BFINDEX(I) = ( I * ( I - 1 ) ) / 2
   END DO
   ALLOCATE ( JIND4(1:4,1:4,1:4,1:4,1:2), KIND4(1:4,1:4,1:4,1:4,1:2) )
   KK = 0
   LL = 0
   DO I = 1,4
      DO J = 1,I
         LL = LL + 1
         JIND4(1,1,I,J,2) = LL
         JIND4(1,1,J,I,2) = LL
         JIND4(I,J,1,1,2) = LL
         JIND4(J,I,1,1,2) = LL
         KIND4(I,1,1,J,2) = LL
         KIND4(J,1,1,I,2) = LL
         KIND4(I,1,J,1,2) = LL
         KIND4(J,1,I,1,2) = LL
         DO K = 1,4
            DO L = 1,K
               KK = KK + 1
               KIND4(K,I,J,L,1) = KK
               KIND4(K,J,I,L,1) = KK
               KIND4(L,I,J,K,1) = KK
               KIND4(L,J,I,K,1) = KK
               JIND4(K,L,I,J,1) = KK
               JIND4(L,K,I,J,1) = KK
               JIND4(K,L,J,I,1) = KK
               JIND4(L,K,J,I,1) = KK
            END DO
         END DO
      END DO
   END DO
   JIND2 = RESHAPE( JIND4, SHAPE( JIND2 ) )
   JIND3 = RESHAPE( JIND4, SHAPE( JIND3 ) )
   KIND2 = RESHAPE( KIND4, SHAPE( KIND2 ) )
   DEALLOCATE ( JIND4, KIND4 )
   CALL DENSITY_GUESS
   CALL MOPAC_SCF_INITIALIZE
   END SUBROUTINE MOPAC_SETUP


   SUBROUTINE MOPAC_INTEGRALS( ENUCLEAR, GRADIENT )
   REAL*8, INTENT(OUT) :: ENUCLEAR
   REAL*8, DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: GRADIENT
   INTEGER            :: I, IBASIS, IFIRST, INDTEI, ILAST, IQM, J, JBASIS, &
                         JFIRST, JLAST, JQM, NI, NJ, NQMINTH
   LOGICAL            :: QAM1PM3, QGRADIENT
   REAL*8 :: BETAIJ, ENUC, ENUCLR
   REAL*8, DIMENSION(1:3) :: RI, RJ
   REAL*8, DIMENSION(1:10) :: E1B, E2A
   REAL*8, ALLOCATABLE, DIMENSION(:) :: OVERLAP

   QAM1PM3 = ( HAMILTONIAN == "AM1 " ) .OR. ( HAMILTONIAN == "PDDG " ) .OR. ( HAMILTONIAN == "PM3 " ) .OR. ( HAMILTONIAN == "RM1 " )
   QGRADIENT = PRESENT ( GRADIENT )

   IF( ALLOCATED ( HCORE   ) ) DEALLOCATE( HCORE   )
   IF( ALLOCATED ( SETEI   ) ) DEALLOCATE( SETEI   )
   ALLOCATE ( HCORE(1:NBASTR), OVERLAP(1:NBASTR), SETEI(1:N2ELEC) )
   HCORE = 0.0D0
   I = 0
   DO IBASIS = 1,NBASIS
      I = I + IBASIS
      HCORE(I) = USPD(IBASIS)
   END DO
   CALL GAUSSIAN_OVERLAP( OVERLAP )
   I = 0
   DO IBASIS = 1,NBASIS
      DO JBASIS = 1,(IBASIS-1)
         BETAIJ = 0.5D0 * ( BETA(IBASIS) + BETA(JBASIS) )
         I = I + 1
         OVERLAP(I) = BETAIJ * OVERLAP(I)
      END DO
      I = I + 1
      OVERLAP(I) = 0.0D0
   END DO
   HCORE(1:NBASTR) = HCORE(1:NBASTR) + OVERLAP(1:NBASTR)
   DEALLOCATE( OVERLAP )
   INDTEI  = 1
   NQMINTH = 0
   ENUCLR  = 0.0D0
   DO IQM = 1,NAQM
      NI     = ANINT( ATMCHG(IQM) )
      RI     = ATMCRD(1:3,IQM)
      IFIRST = BFIRST(IQM)
      ILAST  = BLAST(IQM)
      DO JQM = 1,(IQM-1)
         NJ     = ANINT( ATMCHG(JQM) )
         RJ     = ATMCRD(1:3,JQM)
         JFIRST = BFIRST(JQM)
         JLAST  = BLAST(JQM)
         CALL ROTATEI( NI, NJ, RI, RJ, SETEI(INDTEI:), INDTEI, E1B, E2A, ENUC, QAM1PM3 )
         ENUCLR = ENUCLR + ENUC
         IF( NATORB(NI) > 0 ) THEN
            J = 0
            DO IBASIS = IFIRST,ILAST
               I = (IBASIS*(IBASIS-1))/2 + IFIRST-1
               DO JBASIS = IFIRST,IBASIS
                  I = I + 1
                  J = J + 1
                  HCORE(I) = HCORE(I) + E1B(J)
               END DO
            END DO
         END IF
         IF( NATORB(NJ) > 0 ) THEN
            J = 0
            DO IBASIS = JFIRST,JLAST
               I = (IBASIS*(IBASIS-1))/2 + JFIRST-1
               DO JBASIS = JFIRST,IBASIS
                  I = I + 1
                  J = J + 1
                  HCORE(I) = HCORE(I) + E2A(J)
               END DO
            END DO
         END IF
      END DO
      IF ( NATORB(NI) > 1 ) NQMINTH = NQMINTH + 1
   END DO
   NQMINT = NAQM * ( NATM - NAQM )
   NQMINTH = NQMINTH * ( NATM - NAQM )
   IF( NAQM < NATM ) THEN
      IF ( QGRADIENT ) THEN
         ALLOCATE( DXE1BA(1:NQMINT), DYE1BA(1:NQMINT), DZE1BA(1:NQMINT), &
                   DXE1BB(2:10,1:NQMINTH), DYE1BB(2:10,1:NQMINTH), DZE1BB(2:10,1:NQMINTH) )
      END IF
      CALL SEINTC
   END IF
   ENUCLEAR = ENUCLR
   CONTAINS
      SUBROUTINE SEINTC
      INTEGER            :: I, II, IINT, J, JJ, N1, N2, MATOM, NINTQ, NQM, NTERM, QATOM
      REAL*8 :: ALPHA, CGQM, D1, D2, D4, F1, F2, F3, GAMMA, RHO0, RHO1, RHO2, R2OFF, R2ON, SWITCH, XQ, YQ, ZQ
      REAL*8 :: SS,   SPZ,   PZPZ,   PPPP,  &
                            DSS,  DSPZ,  DPZPZ,  DPPPP, &
                            SFCT1M,   SFCT1P,  SFCT20,  SFCT2M,  SFCT2P,  SFCT4P, &
                            DSFCT1M, DSFCT1P, DSFCT20, DSFCT2M, DSFCT2P, DSFCT4P
      REAL*8 :: CGL, DXENUC, DYENUC, DZENUC, DXPPPP, DYPPPP, DZPPPP, &
                            DXPZPZ, DYPZPZ, DZPZPZ, DXSPZ, DYSPZ, DZSPZ,         &
                            DXSS, DYSS, DZSS, ENUC, RQM, RQM2, RQMB, RQMI,       &
                            TEMP1, TEMP2, TEMP3, TEMP4, TEMP5,                   &
                            XN, YN, ZN, XN2, YN2, ZN2, XQM, YQM, ZQM
      REAL*8, DIMENSION(1:3)  :: GVEC, RQ
      REAL*8, DIMENSION(1:10) :: E1B
      REAL*8, PARAMETER :: ALPMM = 5.0D0, RHO0MM = 0.0D0

      R2OFF = CUT_OFF**2
      R2ON  = CUT_ON**2
      GAMMA = ( R2OFF - R2ON )**3
      CUTOFFB = ANGSTROMS_TO_BOHRS * CUT_OFF
      CUTONB  = ANGSTROMS_TO_BOHRS * CUT_ON
      R2OFFB  = CUTOFFB**2
      R2ONB   = CUTONB**2
      GAMMAB  = ( R2OFFB - R2ONB )**3
      T6FAC =   2.0D0 / GAMMAB
      T4FAC = - 3.0D0 * ( R2OFFB + R2ONB ) / GAMMAB
      T2FAC =   6.0D0 * R2OFFB * R2ONB / GAMMAB
      T0FAC =   R2OFFB * R2OFFB * ( R2OFFB - 3.0D0 * R2ONB ) / GAMMAB

      IF ( QGRADIENT ) THEN
         DXE1BA = 0.0D0 ; DXE1BB = 0.0D0
         DYE1BA = 0.0D0 ; DYE1BB = 0.0D0
         DZE1BA = 0.0D0 ; DZE1BB = 0.0D0
      END IF
      IINT  = 0
      NINTQ = 0
      DO QATOM = 1, NAQM
         NQM   = ANINT( ATMCHG(QATOM) )
         RQ    = ATMCRD(1:3,QATOM)
         ALPHA = ALP(NQM)
         CGQM  = CORE(NQM)
         N1    = BFIRST(QATOM)
         N2    = BLAST(QATOM)
         XQ    = RQ(1)
         YQ    = RQ(2)
         ZQ    = RQ(3)
         RHO0  = ( 0.5D0 / AM(NQM) + RHO0MM )**2
         IF( QAM1PM3 ) NTERM = COUNT ( ABS ( FN1(NQM,1:4) ) > 0.0D0 )
         DO MATOM = NAQM+1, NATM
            CGL = ATMCHG(MATOM)
            XQM = ( XQ - ATMCRD(1,MATOM) )
            YQM = ( YQ - ATMCRD(2,MATOM) )
            ZQM = ( ZQ - ATMCRD(3,MATOM) )
            RQM2 = XQM * XQM + YQM * YQM + ZQM * ZQM
            IINT = IINT + 1
            IF ( NATORB(NQM) > 1 ) NINTQ = NINTQ + 1
            IF ( RQM2 > R2OFF ) CYCLE
            RQM  = SQRT ( RQM2 )
            RQMI = 1.0D0 / RQM
            XN   = RQMI * XQM
            YN   = RQMI * YQM
            ZN   = RQMI * ZQM
            IF ( RQM2 > R2ON ) THEN
               SWITCH = ( R2OFF - RQM2 )**2 * ( R2OFF + 2.0D0 * RQM2 - 3.0D0 * R2ON ) / GAMMA
            ELSE
               SWITCH = 1.0D0
            END IF
            IF ( NATORB(NQM) == 1 ) THEN
               RQMB = ANGSTROMS_TO_BOHRS * RQM
               CALL SEINTC_INTEGRALB ( RQMB, RHO0, SWITCH, SS, DSS )
               E1B(1)    = AU_TO_EV * SS
               E1B(2:10) = 0.0D0
               IF ( QGRADIENT ) THEN
                  DSS = ANGSTROMS_TO_BOHRS * DSS
                  DXE1BA(IINT) = XN * DSS
                  DYE1BA(IINT) = YN * DSS
                  DZE1BA(IINT) = ZN * DSS
               END IF
            ELSE
               RHO1 = ( 0.5D0 / AD(NQM) + RHO0MM )**2
               RHO2 = ( 0.5D0 / AQ(NQM) + RHO0MM )**2
               D1   = DD(NQM)
               D2   = 2.0D0 * QQ(NQM)
               D4   = D2 * D2
               XN2  = XN * XN
               YN2  = YN * YN
               ZN2  = ZN * ZN
               RQMB = ANGSTROMS_TO_BOHRS * RQM
               CALL SEINTC_INTEGRALB ( RQMB, RHO0, SWITCH, SS, DSS )
               CALL SEINTC_INTEGRALAB ( RQMB,  D1, RHO1, SWITCH, SFCT1P, DSFCT1P )
               CALL SEINTC_INTEGRALAB ( RQMB, -D1, RHO1, SWITCH, SFCT1M, DSFCT1M )
               SPZ  = 0.5D0 * (  SFCT1P -  SFCT1M )
               DSPZ = 0.5D0 * ( DSFCT1P - DSFCT1M )
               CALL SEINTC_INTEGRALB  ( RQMB,      RHO2, SWITCH, SFCT20, DSFCT20 )
               CALL SEINTC_INTEGRALAB ( RQMB,  D2, RHO2, SWITCH, SFCT2P, DSFCT2P )
               CALL SEINTC_INTEGRALAB ( RQMB, -D2, RHO2, SWITCH, SFCT2M, DSFCT2M )
               PZPZ  =  SS + 0.25D0 * (  SFCT2P +  SFCT2M ) - 0.5D0 *  SFCT20
               DPZPZ = DSS + 0.25D0 * ( DSFCT2P + DSFCT2M ) - 0.5D0 * DSFCT20
               CALL SEINTC_INTEGRALB ( RQMB, ( D4 + RHO2 ), SWITCH, SFCT4P, DSFCT4P )
               PPPP  =  SS + 0.5D0 * (  SFCT4P -  SFCT20 )
               DPPPP = DSS + 0.5D0 * ( DSFCT4P - DSFCT20 )
               E1B(1)  = SS
               E1B(2)  = XN * SPZ
               E1B(3)  = XN2 * PZPZ + ( YN2 + ZN2 ) * PPPP
               E1B(4)  = YN * SPZ
               E1B(5)  = XN * YN * ( PZPZ - PPPP )
               E1B(6)  = YN2 * PZPZ + ( XN2 + ZN2 ) * PPPP
               E1B(7)  = ZN * SPZ
               E1B(8)  = XN * ZN * ( PZPZ - PPPP )
               E1B(9)  = YN * ZN * ( PZPZ - PPPP )
               E1B(10) = ZN2 * PZPZ + ( XN2 + YN2 ) * PPPP
               E1B(1)    = AU_TO_EV * E1B(1)
               E1B(2:10) = AU_TO_EV * CGL * E1B(2:10)
               IF ( QGRADIENT ) THEN
                  DSS   = ANGSTROMS_TO_BOHRS * DSS
                  DSPZ  = ANGSTROMS_TO_BOHRS * DSPZ
                  DPZPZ = ANGSTROMS_TO_BOHRS * DPZPZ
                  DPPPP = ANGSTROMS_TO_BOHRS * DPPPP
                  DXSS  = XN * DSS
                  DYSS  = YN * DSS
                  DZSS  = ZN * DSS
                  DXSPZ = XN * DSPZ
                  DYSPZ = YN * DSPZ
                  DZSPZ = ZN * DSPZ
                  DXPZPZ = XN * DPZPZ
                  DYPZPZ = YN * DPZPZ
                  DZPZPZ = ZN * DPZPZ
                  DXPPPP = XN * DPPPP
                  DYPPPP = YN * DPPPP
                  DZPPPP = ZN * DPPPP
                  TEMP1 = PZPZ - PPPP
                  TEMP2 = YN2 + ZN2
                  TEMP3 = XN2 + ZN2
                  TEMP4 = RQMI * SPZ
                  TEMP5 = XN2 + YN2
                  DXE1BA(IINT) = DXSS
                  DYE1BA(IINT) = DYSS
                  DZE1BA(IINT) = DZSS
                  DXE1BB(2,NINTQ) = TEMP2 * TEMP4 + XN * DXSPZ
                  DYE1BB(2,NINTQ) = -XN * (TEMP4 * YN - DYSPZ)
                  DZE1BB(2,NINTQ) = -XN * (TEMP4 * ZN - DZSPZ)
                  DXE1BB(3,NINTQ) = XN2 * DXPZPZ + TEMP2 * DXPPPP + 2.0D0 * XN * RQMI * TEMP2 * TEMP1
                  DYE1BB(3,NINTQ) = XN2 * DYPZPZ + TEMP2 * DYPPPP - 2.0D0 * YN * RQMI * XN2 * TEMP1
                  DZE1BB(3,NINTQ) = XN2 * DZPZPZ + TEMP2 * DZPPPP - 2.0D0 * ZN * RQMI * XN2 * TEMP1
                  DXE1BB(4,NINTQ) = -YN * (XN * TEMP4 - DXSPZ)
                  DYE1BB(4,NINTQ) = TEMP4 * TEMP3 + YN * DYSPZ
                  DZE1BB(4,NINTQ) = -YN * (ZN * TEMP4 - DZSPZ)
                  DXE1BB(5,NINTQ) = YN * (RQMI * TEMP1 * (TEMP2 - XN2) + XN * (DXPZPZ - DXPPPP))
                  DYE1BB(5,NINTQ) = XN * (RQMI * TEMP1 * (TEMP3 - YN2) + YN * (DYPZPZ - DYPPPP))
                  DZE1BB(5,NINTQ) = XN * YN * ((DZPZPZ - DZPPPP) - 2.0D0 * ZN * RQMI * TEMP1)
                  DXE1BB(6,NINTQ) = YN2 * DXPZPZ + TEMP3 * DXPPPP - 2.0D0 * XN * YN2 * RQMI * TEMP1
                  DYE1BB(6,NINTQ) = YN2 * DYPZPZ + TEMP3 * (DYPPPP + 2.0D0 * YN * RQMI * TEMP1)
                  DZE1BB(6,NINTQ) = YN2 * DZPZPZ + TEMP3 * DZPPPP - 2.0D0 * ZN * YN2 * RQMI * TEMP1
                  DXE1BB(7,NINTQ) = -ZN * (XN * TEMP4 - DXSPZ)
                  DYE1BB(7,NINTQ) = -ZN * (YN * TEMP4 - DYSPZ)
                  DZE1BB(7,NINTQ) = TEMP5 * TEMP4 + ZN * DZSPZ
                  DXE1BB(8,NINTQ) = ZN * (RQMI * TEMP1 * (TEMP2 - XN2) + XN * (DXPZPZ - DXPPPP))
                  DYE1BB(8,NINTQ) = XN * ZN * ((DYPZPZ - DYPPPP) - 2.0D0 * YN * RQMI * TEMP1)
                  DZE1BB(8,NINTQ) = XN * (RQMI * TEMP1 * (TEMP5 - ZN2) + ZN * (DZPZPZ - DZPPPP))
                  DXE1BB(9,NINTQ) = YN * ZN * ((DXPZPZ - DXPPPP) - 2.0D0 * XN * RQMI * TEMP1)
                  DYE1BB(9,NINTQ) = ZN * (RQMI * TEMP1 * (TEMP3 - YN2) + YN * (DYPZPZ - DYPPPP))
                  DZE1BB(9,NINTQ) = YN * (RQMI * TEMP1 * (TEMP5 - ZN2) + ZN * (DZPZPZ - DZPPPP))
                  DXE1BB(10,NINTQ) = ZN2 * DXPZPZ + TEMP5 * DXPPPP - 2.0D0 * XN * ZN2 * RQMI * TEMP1
                  DYE1BB(10,NINTQ) = ZN2 * DYPZPZ + TEMP5 * DYPPPP - 2.0D0 * YN * ZN2 * RQMI * TEMP1
                  DZE1BB(10,NINTQ) = ZN2 * DZPZPZ + TEMP5 * (DZPPPP + 2.0D0 * ZN * RQMI * TEMP1)
                  DXE1BB(2:10,NINTQ) = CGL * DXE1BB(2:10,NINTQ)
                  DYE1BB(2:10,NINTQ) = CGL * DYE1BB(2:10,NINTQ)
                  DZE1BB(2:10,NINTQ) = CGL * DZE1BB(2:10,NINTQ)
               END IF
            END IF
            TEMP1 = E1B(1)
            TEMP2 = CGQM * CGL
            TEMP3 = ABS ( TEMP2 ) * EXP ( - ALPHA * RQM )
            TEMP4 = ABS ( TEMP2 ) * EXP ( - ALPMM * RQM )
            ENUC  = TEMP1 * ( TEMP2 + TEMP3 + TEMP4 )
            IF ( QGRADIENT ) THEN
               DXENUC = AU_TO_EV * ( TEMP2 + TEMP3 + TEMP4 ) * DXE1BA(IINT) - XN * TEMP1 * ( ALPHA * TEMP3 + ALPMM * TEMP4 )
               DYENUC = AU_TO_EV * ( TEMP2 + TEMP3 + TEMP4 ) * DYE1BA(IINT) - YN * TEMP1 * ( ALPHA * TEMP3 + ALPMM * TEMP4 )
               DZENUC = AU_TO_EV * ( TEMP2 + TEMP3 + TEMP4 ) * DZE1BA(IINT) - ZN * TEMP1 * ( ALPHA * TEMP3 + ALPMM * TEMP4 )
            END IF
            IF ( QAM1PM3 ) THEN
               DO I = 1,NTERM
                  F1    = FN1(NQM,I)
                  F2    = FN2(NQM,I)
                  F3    = FN3(NQM,I)
                  TEMP1 = EXP ( MAX ( -30.0D0, - F2 * ( RQM - F3 )**2 ) )
                  ENUC  = ENUC + CGQM * CGL * RQMI * F1 * TEMP1
                  IF ( QGRADIENT ) THEN
                     TEMP2 = CGQM * CGL / RQM2
                     DXENUC = DXENUC - F1 * TEMP1 * TEMP2 * XN * ( 1.0D0 + 2.0D0 * F2 * RQM * ( RQM - F3 ) )
                     DYENUC = DYENUC - F1 * TEMP1 * TEMP2 * YN * ( 1.0D0 + 2.0D0 * F2 * RQM * ( RQM - F3 ) )
                     DZENUC = DZENUC - F1 * TEMP1 * TEMP2 * ZN * ( 1.0D0 + 2.0D0 * F2 * RQM * ( RQM - F3 ) )
                  END IF
               END DO
            END IF
            E1B(1) = CGL * E1B(1)
            IF( QGRADIENT ) THEN
               DXE1BA(IINT) = CGL * DXE1BA(IINT)
               DYE1BA(IINT) = CGL * DYE1BA(IINT)
               DZE1BA(IINT) = CGL * DZE1BA(IINT)
            END IF
            ENUCLR = ENUCLR + ENUC
            JJ = 0
            DO I = N1,N2
               II = (I * (I - 1)) / 2 + N1 - 1
               DO J = N1,I
                  II = II + 1
                  JJ = JJ + 1
                  HCORE(II) = HCORE(II) - E1B(JJ)
               END DO
            END DO
            IF( QGRADIENT ) THEN
               TEMP1 = EV_TO_KJ
               GVEC(1) = TEMP1 * DXENUC
               GVEC(2) = TEMP1 * DYENUC
               GVEC(3) = TEMP1 * DZENUC
               GRADIENT(1:3,MATOM) = GRADIENT(1:3,MATOM) - GVEC
               GRADIENT(1:3,QATOM) = GRADIENT(1:3,QATOM) + GVEC
            END IF
         END DO
      END DO
      IF( QGRADIENT ) THEN
         DXE1BA = AU_TO_EV * EV_TO_KJ * DXE1BA
         DYE1BA = AU_TO_EV * EV_TO_KJ * DYE1BA
         DZE1BA = AU_TO_EV * EV_TO_KJ * DZE1BA
         DXE1BB = AU_TO_EV * EV_TO_KJ * DXE1BB
         DYE1BB = AU_TO_EV * EV_TO_KJ * DYE1BB
         DZE1BB = AU_TO_EV * EV_TO_KJ * DZE1BB
      END IF
      END SUBROUTINE SEINTC
   END SUBROUTINE MOPAC_INTEGRALS


   SUBROUTINE REPP ( NI, NJ, RIJ, RI, CORINT )
   ! -----------------------------------------------------------------------------
   !     REPP calculates the two-electron integrals and electron-nuclear
   !     attraction terms in the local frame between the atoms I and J.
   !
   !     Input:   NI, NJ        the atomic numbers.
   !              RIJ           interatomic distance in atomic units.
   !              AM, AD, AQ    arrays of two-electron one-centre integrals.
   !              CORE          array on atomic nuclear charges.
   !              DD            array of dipole charge separations.
   !              QQ            array of quadrupole charge separations.
   !
   !    Output:   RI            two-electron repulsion integrals.
   !              CORINT          array of electron-nuclear attraction integrals.
   ! -----------------------------------------------------------------------------
   INTEGER,            INTENT(IN) :: NI, NJ
   REAL*8, INTENT(IN) :: RIJ
   REAL*8, DIMENSION(1:4,1:2), INTENT(OUT) :: CORINT
   REAL*8, DIMENSION(1:22),    INTENT(OUT) :: RI
   REAL*8 :: ADD, ADE, ADQ, AED, AEE, AEQ, AQD, AQE, AQQ, DA, DB,  &
                         DXDX, DXQXZ, DZDZ, DZE, DZQXX, DZQZZ, EDZ, EE, EQXX,  &
                         EQZZ, QA, QB, QXXDZ, QXXE, QXXQXX, QXXQYY, QXXQZZ,    &
                         QXYQXY, QXZQXZ, QXZDX, QZZDZ, QZZE, QZZQXX, QZZQZZ, R
   REAL*8, PARAMETER :: EIGHT = 8.0D0, FD = 4.0D0, FOUR = 4.0D0, OD = 1.0D0, &
                                    PT25  = 0.25D0, SIXTN = 16.0D0, TD = 2.0D0

   R      = RIJ
   CORINT = 0.0D0
   RI     = 0.0D0
   DA = DD(NI)
   DB = DD(NJ)
   QA = QQ(NI)
   QB = QQ(NJ)
   AEE = 0.25D0 * ( OD/AM(NI) + OD/AM(NJ) )**2
   EE  = OD / SQRT ( R**2+AEE )
   RI(1) = EE * AU_TO_EV
   CORINT(1,1) = CORE(NJ) * RI(1)
   CORINT(1,2) = CORE(NI) * RI(1)
   IF ( (NATORB(NI) < 3) .AND. (NATORB(NJ) < 3) ) RETURN
   IF ( (NATORB(NI) < 3) ) GOTO 10
   ADE=0.25D0*(OD/AD(NI)+OD/AM(NJ))**2
   AQE=0.25D0*(OD/AQ(NI)+OD/AM(NJ))**2
   DZE=-OD/SQRT((R+DA)**2+ADE)+OD/SQRT((R-DA)**2+ADE)
   QZZE=OD/SQRT((R-TD*QA)**2+AQE)-TD/SQRT(R**2+AQE)+OD/SQRT((R+TD*QA)**2+AQE)
   QXXE=TD/SQRT(R**2+FD*QA**2+AQE)-TD/SQRT(R**2+AQE)
   DZE=DZE/2.0D0
   QXXE=QXXE/FOUR
   QZZE=QZZE/FOUR
   RI(2)=-DZE
   RI(3)=EE+QZZE
   RI(4)=EE+QXXE
   IF ( NATORB(NJ) < 3 ) GOTO 20
   10 CONTINUE
   AED=0.25D0*(OD/AM(NI)+OD/AD(NJ))**2
   AEQ=0.25D0*(OD/AM(NI)+OD/AQ(NJ))**2
   EDZ=-OD/SQRT((R-DB)**2+AED)+OD/SQRT((R+DB)**2+AED)
   EQZZ=OD/SQRT((R-TD*QB)**2+AEQ)-TD/SQRT(R**2+AEQ)+OD/SQRT((R+TD*QB)**2+AEQ)
   EQXX=TD/SQRT(R**2+FD*QB**2+AEQ)-TD/SQRT(R**2+AEQ)
   EDZ=EDZ/2.0D0
   EQXX=EQXX/FOUR
   EQZZ=EQZZ/FOUR
   RI(5)=-EDZ
   RI(11)=EE+EQZZ
   RI(12)=EE+EQXX
   IF ( NATORB(NI) < 3 ) GOTO 20
   ADD=0.25D0*(OD/AD(NI)+OD/AD(NJ))**2
   ADQ=0.25D0*(OD/AD(NI)+OD/AQ(NJ))**2
   AQD=0.25D0*(OD/AQ(NI)+OD/AD(NJ))**2
   AQQ=0.25D0*(OD/AQ(NI)+OD/AQ(NJ))**2
   DXDX=TD/SQRT(R**2+(DA-DB)**2+ADD)-TD/SQRT(R**2+(DA+DB)**2+ADD)
   DZDZ=OD/SQRT((R+DA-DB)**2+ADD)+OD/SQRT((R-DA+DB)**2+ADD)-OD/SQRT((R-DA-DB)**2+ADD)-OD/SQRT((R+DA+DB)**2+ADD)
   DZQXX=-TD/SQRT((R+DA)**2+FD*QB**2+ADQ)+TD/SQRT((R-DA)**2+FD*QB**2+ADQ)+TD/SQRT((R+DA)**2+ADQ)-TD/SQRT((R-DA)**2+ADQ)
   QXXDZ=-TD/SQRT((R-DB)**2+FD*QA**2+AQD)+TD/SQRT((R+DB)**2+FD*QA**2+AQD)+TD/SQRT((R-DB)**2+AQD)-TD/SQRT((R+DB)**2+AQD)
   DZQZZ=-OD/SQRT((R+DA-TD*QB)**2+ADQ)+OD/SQRT((R-DA-TD*QB)**2+ADQ)-OD/SQRT((R+DA+TD*QB)**2+ADQ)+OD/SQRT((R-DA+TD*QB) &
          **2+ADQ)-TD/SQRT((R-DA)**2+ADQ)+TD/SQRT((R+DA)**2+ADQ)
   QZZDZ=-OD/SQRT((R+TD*QA-DB)**2+AQD)+OD/SQRT((R+TD*QA+DB)**2+AQD)-OD/SQRT((R-TD*QA-DB)**2+AQD)+OD/SQRT((R-2.0D0 &
          *QA+DB)**2+AQD)+TD/SQRT((R-DB)**2+AQD)-TD/SQRT((R+DB)**2+AQD)
   QXXQXX=TD/SQRT(R**2+FD*(QA-QB)**2+AQQ)+TD/SQRT(R**2+FD*(QA+QB)**2+ &
          AQQ)-FD/SQRT(R**2+FD*QA**2+AQQ)-FD/SQRT(R**2+FD*QB**2+AQQ)+FD/SQRT(R**2+AQQ)
   QXXQYY=FD/SQRT(R**2+FD*QA**2+FD*QB**2+AQQ)-FD/SQRT(R**2+FD*QA**2+AQQ)-FD/SQRT(R**2+FD*QB**2+AQQ)+FD/SQRT(R**2+AQQ)
   QXXQZZ=TD/SQRT((R-TD*QB)**2+FD*QA**2+AQQ)+TD/SQRT((R+TD*QB)**2+FD* &
          QA**2+AQQ)-TD/SQRT((R-TD*QB)**2+AQQ)-TD/SQRT((R+TD*QB)**2+AQQ)-FD/ &
          SQRT(R**2+FD*QA**2+AQQ)+FD/SQRT(R**2+AQQ)
   QZZQXX=TD/SQRT((R+TD*QA)**2+FD*QB**2+AQQ)+TD/SQRT((R-TD*QA)**2+FD* &
          QB**2+AQQ)-TD/SQRT((R+TD*QA)**2+AQQ)-TD/SQRT((R-TD*QA)**2+AQQ)-FD/ &
          SQRT(R**2+FD*QB**2+AQQ)+FD/SQRT(R**2+AQQ)
   QZZQZZ=OD/SQRT((R+TD*QA-TD*QB)**2+AQQ)+OD/SQRT((R+TD*QA+TD*QB)**2+ &
          AQQ)+OD/SQRT((R-TD*QA-TD*QB)**2+AQQ)+OD/SQRT((R-TD*QA+TD*QB)**2+AQQ) &
          -TD/SQRT((R-TD*QA)**2+AQQ)-TD/SQRT((R+TD*QA)**2+AQQ)-TD/SQRT((R- &
          TD*QB)**2+AQQ)-TD/SQRT((R+TD*QB)**2+AQQ)+FD/SQRT(R**2+AQQ)
   DXQXZ=-TD/SQRT((R-QB)**2+(DA-QB)**2+ADQ)+TD/SQRT((R+QB)**2+(DA-QB) &
         **2+ADQ)+TD/SQRT((R-QB)**2+(DA+QB)**2+ADQ)-TD/SQRT((R+QB)**2+(DA+QB)**2+ADQ)
   QXZDX=-TD/SQRT((R+QA)**2+(QA-DB)**2+AQD)+TD/SQRT((R-QA)**2+(QA-DB) &
         **2+AQD)+TD/SQRT((R+QA)**2+(QA+DB)**2+AQD)-TD/SQRT((R-QA)**2+(QA+DB)**2+AQD)
   QXYQXY=FD/SQRT(R**2+TD*(QA-QB)**2+AQQ)+FD/SQRT(R**2+TD*(QA+QB)**2+AQQ)-EIGHT/SQRT(R**2+TD*(QA**2+QB**2)+AQQ)
   QXZQXZ=TD/SQRT((R+QA-QB)**2+(QA-QB)**2+AQQ)-TD/SQRT((R+QA+QB)**2+( & 
          QA-QB)**2+AQQ)-TD/SQRT((R-QA-QB)**2+(QA-QB)**2+AQQ)+TD/SQRT((R-QA+  &
          QB)**2+(QA-QB)**2+AQQ)-TD/SQRT((R+QA-QB)**2+(QA+QB)**2+AQQ)+TD/SQRT &
          ((R+QA+QB)**2+(QA+QB)**2+AQQ)+TD/SQRT((R-QA-QB)**2+(QA+QB)**2+AQQ &
          )-TD/SQRT((R-QA+QB)**2+(QA+QB)**2+AQQ)
   DXDX=DXDX/FOUR
   DZDZ=DZDZ/FOUR
   DZQXX=DZQXX/EIGHT
   QXXDZ=QXXDZ/EIGHT
   DZQZZ=DZQZZ/EIGHT
   QZZDZ=QZZDZ/EIGHT
   DXQXZ=DXQXZ/EIGHT
   QXZDX=QXZDX/EIGHT
   QXXQXX=QXXQXX/16.0D0
   QXXQYY=QXXQYY/16.0D0
   QXXQZZ=QXXQZZ/16.0D0
   QZZQXX=QZZQXX/16.0D0
   QZZQZZ=QZZQZZ/16.0D0
   QXZQXZ=QXZQXZ/16.0D0
   QXYQXY=QXYQXY/16.0D0
   RI(6)=DZDZ
   RI(7)=DXDX
   RI(8)=-EDZ-QZZDZ
   RI(9)=-EDZ-QXXDZ
   RI(10)=-QXZDX
   RI(13)=-DZE-DZQZZ
   RI(14)=-DZE-DZQXX
   RI(15)=-DXQXZ
   RI(16)=EE+EQZZ+QZZE+QZZQZZ
   RI(17)=EE+EQZZ+QXXE+QXXQZZ
   RI(18)=EE+EQXX+QZZE+QZZQXX
   RI(19)=EE+EQXX+QXXE+QXXQXX
   RI(20)=QXZQXZ
   RI(21)=EE+EQXX+QXXE+QXXQYY
   RI(22)=0.5D0*(QXXQXX-QXXQYY)
   20 CONTINUE
   RI(2:22) = RI(2:22) * AU_TO_EV
   CORINT(2,1) = CORE(NJ) * RI(2)
   CORINT(3,1) = CORE(NJ) * RI(3)
   CORINT(4,1) = CORE(NJ) * RI(4)
   CORINT(2,2) = CORE(NI) * RI(5)
   CORINT(3,2) = CORE(NI) * RI(11)
   CORINT(4,2) = CORE(NI) * RI(12)
   END SUBROUTINE REPP


   SUBROUTINE ROTATEI ( NI, NJ, XI, XJ, W, KR, E1B, E2A, ENUC, QAM1PM3 )
   ! -----------------------------------------------------------------------------
   !     ROTATE calculates the two-electron integrals and the electron-nuclear
   !     attraction integrals. The integrals are evaluated first in the local
   !     frame by REPP and then transformed to the molecular frame by ROTATE.
   !
   !     Input:     NI, NJ      atomic numbers of first and second atoms.
   !                XI, XJ      coordinates of atoms (in atomic units).
   !
   !     Output:    W           array of two-electron integrals.
   !                E1B         attraction integrals for electron of atom NI
   !                            and nucleus of NJ.
   !                E2A         attraction integrals for electron of atom NJ
   !                            and nucleus of NI.
   !                ENUC        nuclear-nuclear repulsion term.
   !
   !
   !     The integrals in the local frame are stored in RI and they have the
   !     order :
   !
   !         (ss/ss) = 1,     (so/ss) = 2,   (oo/ss) = 3,   (pp/ss) = 4,
   !         (ss/os) = 5,     (so/so) = 6,   (sp/sp) = 7,   (oo/so) = 8,
   !         (pp/so) = 9,     (po/sp) = 10,  (ss/oo) = 11,  (ss/pp) = 12,
   !         (so/oo) = 13,    (so/pp) = 14,  (sp/op) = 15,  (oo/oo) = 16,
   !         (pp/oo) = 17,    (oo/pp) = 18,  (pp/pp) = 19,  (po/po) = 20,
   !                       (pp/p*p*) = 21,  (p*p/p*p) = 22
   !
   !     where o is a p-sigma orbital and p and p* are p-pi orbitals.
   !
   !     The nuclear attraction integrals are stored in CORINT(kl/ij) in the
   !     order :
   !              (ss/) = 1,  (so/) = 2, (oo/) = 3, (pp/) = 4
   !
   !     where ij = 1 if the orbitals are on atom NI and ij = 2 if they are
   !     on atom NJ.
   ! -----------------------------------------------------------------------------
   INTEGER,            INTENT(IN)    :: NI, NJ
   INTEGER,            INTENT(INOUT) :: KR
   LOGICAL,            INTENT(IN)    :: QAM1PM3
   REAL*8, INTENT(OUT)   :: ENUC
   REAL*8, DIMENSION(1:10), INTENT(OUT)   :: E1B, E2A
   REAL*8, DIMENSION(1:3),  INTENT(IN)    :: XI, XJ
   REAL*8, DIMENSION(:),    INTENT(INOUT) :: W
   REAL*8, PARAMETER :: PDDG_EXPONENT = 10.0D0
   INTEGER            :: I, IB, IG, II, IJ, IK, J, JB, JJ, JK, K, KI, KK, L, LL, NT
   LOGICAL            :: SI, SK
   REAL*8 :: A, AX, D, GAM, RIJ, SCALE, ZAF, ZBF
   REAL*8, DIMENSION(1:4,1:2) :: CORINT
   REAL*8, DIMENSION(1:22)    :: RI
   REAL*8, DIMENSION(1:3)     :: X, Y, Z
   REAL*8 :: CSS1, CSP1, CPPS1, CPPP1, CSS2, CSP2, CPPS2, CPPP2
   REAL*8 :: X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3
   REAL*8, PARAMETER :: APROX1 = 1.0D0 - 1.0D-8

   X   = ANGSTROMS_TO_BOHRS * ( XI - XJ )
   RIJ = SQRT ( DOT_PRODUCT ( X, X ) )
   CALL REPP ( NI, NJ, RIJ, RI, CORINT )
   GAM = RI(1)
   A = 1.0D0 / RIJ
   X = A * X
   Z(3)=0.0D0
   IF ( ABS ( X(3) ) > APROX1 ) THEN
      X(3)=SIGN(1.0D0,X(3))
      Y(1)=0.0D0
      Y(2)=1.0D0
      Y(3)=0.0D0
      Z(1)=1.0D0
      Z(2)=0.0D0
   ELSE
      Z(3)=SQRT(1.0D0-X(3)**2)
      A=1.0D0/Z(3)
      Y(1)=-A*X(2)*SIGN(1.0D0,X(1))
      Y(2)=ABS(A*X(1))
      Y(3)=0.0D0
      Z(1)=-A*X(1)*X(3)
      Z(2)=-A*X(2)*X(3)
   END IF
   X1 = X(1)
   X2 = X(2)
   X3 = X(3)
   Y1 = Y(1)
   Y2 = Y(2)
   Y3 = Y(3)
   Z1 = Z(1)
   Z2 = Z(2)
   Z3 = Z(3)
   IB=MIN(NATORB(NI),4)
   JB=MIN(NATORB(NJ),4)
   KI=0
   DO I=1,IB
      SI=I.EQ.1
      II=I-1
      DO J=1,I
         JJ=J-1
         IJ=0
         IF (JJ.EQ.0) IJ=-1
         IF (SI) IJ=+1
         DO K=1,JB
            KK=K-1
            SK=KK.GT.0
            DO L=1,K
               KI=KI+1
               IF (SK) GOTO 180
               IF ( IJ < 0 ) THEN
                  W(KI)=RI(2)*X(II)
               ELSE IF ( IJ == 0 ) THEN
                  W(KI)=RI(3)*X(II)*X(JJ)+RI(4)*(Y(II)*Y(JJ)+Z(II)*Z(JJ))
               ELSE
                  W(KI)=RI(1)
               END IF
               GO TO 260
180             LL=L-1
               IF (LL.GT.0) GOTO 220
               IF ( IJ < 0 ) THEN
                  W(KI)=RI(6)*X(II)*X(KK)+RI(7)*(Y(II)*Y(KK)+Z(II)*Z(KK))
               ELSE IF ( IJ == 0 ) THEN
                  W(KI)=X(KK)*(RI(8)*X(II)*X(JJ)+RI(9)*(Y(II)*Y(JJ)+Z(II)*Z(JJ))) &
                        +RI(10)*(X(II)*(Y(JJ)*Y(KK)+Z(JJ)*Z(KK))+X(JJ)*(Y(II)*Y(KK)+Z(II)*Z(KK)))
               ELSE
                  W(KI)=RI(5)*X(KK)
               END IF
               GO TO 260
220             CONTINUE
               IF ( IJ < 0 ) THEN
                  W(KI)=X(II)*(RI(13)*X(KK)*X(LL)+RI(14)*(Y(KK)*Y(LL)+Z(KK)*Z(LL) &
                        ))+RI(15)*(Y(II)*(Y(KK)*X(LL)+Y(LL)*X(KK))+Z(II)*(Z(KK)*X(LL)+Z(LL)*X(KK)))
                ELSE IF ( IJ == 0 ) THEN
                  W(KI)=(RI(16)*X(II)*X(JJ)+RI(17)*(Y(II)*Y(JJ)+Z(II)*Z(JJ)))*X(KK)*X(LL)+   &
                        RI(18)*X(II)*X(JJ)*(Y(KK)*Y(LL)+Z(KK)*Z(LL))+                        &
                        RI(19)*(Y(II)*Y(JJ)*Y(KK)*Y(LL)+Z(II)*Z(JJ)*Z(KK)*Z(LL))+RI(20)      &
                        *(X(II)*(   X(KK)*(Y(JJ)*Y(LL)+Z(JJ)*Z(LL))+X(LL)*(Y(JJ)*Y(KK)+Z(JJ) &
                        *Z(KK))   )+X(JJ)*(X(KK)*(Y(II)*Y(LL)+Z(II)*Z(LL))+X(LL)*(Y(II)*     &
                        Y(KK)+Z(II)*Z(KK))))+RI(21)*(Y(II)*Y(JJ)*Z(KK)*Z(LL)+Z(II)*Z(JJ)     &
                        *Y(KK)*Y(LL))+RI(22)*(Y(II)*Z(JJ)+Z(II)*Y(JJ))*(Y(KK)*Z(LL)+Z(KK)*Y(LL))
               ELSE
                  W(KI)=RI(11)*X(KK)*X(LL)+RI(12)*(Y(KK)*Y(LL)+Z(KK)*Z(LL))
               END IF
260             CONTINUE
            END DO
         END DO
      END DO
   END DO
!280 CONTINUE
   CSS1  = CORINT(1,1)
   CSP1  = CORINT(2,1)
   CPPS1 = CORINT(3,1)
   CPPP1 = CORINT(4,1)
   CSS2  = CORINT(1,2)
   CSP2  = CORINT(2,2)
   CPPS2 = CORINT(3,2)
   CPPP2 = CORINT(4,2)
   E1B = 0.0D0
   E2A = 0.0D0
   IF ( NATORB(NI) > 0 ) THEN
   E1B(1)=-CSS1
   IF(NI.GT.1) THEN
      E1B(2) = -CSP1 *X1
      E1B(3) = -CPPS1*X1**2-CPPP1*(Y1**2+Z1**2)
      E1B(4) = -CSP1 *X2
      E1B(5) = -CPPS1*X1*X2-CPPP1*(Y1*Y2+Z1*Z2)
      E1B(6) = -CPPS1*X2*X2-CPPP1*(Y2*Y2+Z2*Z2)
      E1B(7) = -CSP1 *X3
      E1B(8) = -CPPS1*X1*X3-CPPP1*(Y1*Y3+Z1*Z3)
      E1B(9) = -CPPS1*X2*X3-CPPP1*(Y2*Y3+Z2*Z3)
      E1B(10)= -CPPS1*X3*X3-CPPP1*(Y3*Y3+Z3*Z3)
   END IF
   END IF
   IF ( NATORB(NJ) > 0 ) THEN
   E2A(1)=-CSS2
   IF(NJ.GT.1) THEN
      E2A(2) = -CSP2 *X1
      E2A(3) = -CPPS2*X1**2-CPPP2*(Y1**2+Z1**2)
      E2A(4) = -CSP2 *X2
      E2A(5) = -CPPS2*X1*X2-CPPP2*(Y1*Y2+Z1*Z2)
      E2A(6) = -CPPS2*X2*X2-CPPP2*(Y2*Y2+Z2*Z2)
      E2A(7) = -CSP2 *X3
      E2A(8) = -CPPS2*X1*X3-CPPP2*(Y1*Y3+Z1*Z3)
      E2A(9) = -CPPS2*X2*X3-CPPP2*(Y2*Y3+Z2*Z3)
      E2A(10)= -CPPS2*X3*X3-CPPP2*(Y3*Y3+Z3*Z3)
   END IF
   END IF
      IF(ABS(CORE(NI)).GT.20.AND.ABS(CORE(NJ)).GT.20) THEN
         ENUC=0.0D0
         RETURN
      ELSE IF (RIJ.LT.1.0D0.AND.NATORB(NI)*NATORB(NJ).EQ.0) THEN
         ENUC=0.0D0
         RETURN
      END IF
      RIJ   = RIJ / ANGSTROMS_TO_BOHRS
      SCALE = EXP(-ALP(NI)*RIJ)+EXP(-ALP(NJ)*RIJ)
      NT=NI+NJ
      IF(NT.EQ.8.OR.NT.EQ.9) THEN
         IF(NI.EQ.7.OR.NI.EQ.8) SCALE=SCALE+(RIJ-1.0D0)*EXP(-ALP(NI)*RIJ)
         IF(NJ.EQ.7.OR.NJ.EQ.8) SCALE=SCALE+(RIJ-1.0D0)*EXP(-ALP(NJ)*RIJ)
      END IF
      ENUC = CORE(NI)*CORE(NJ)*GAM
      SCALE=ABS(SCALE*ENUC)
      IF ( QAM1PM3 ) THEN
         DO IG = 1,4
            IF ( ABS ( FN1(NI,IG) ) > 0.0D0 ) THEN
               SCALE = SCALE + (CORE(NI)*CORE(NJ)/RIJ)*FN1(NI,IG)*EXP ( MAX ( -30.0D0, -FN2(NI,IG)*(RIJ-FN3(NI,IG))**2 ) )
            END IF
            IF ( ABS ( FN1(NJ,IG) ) > 0.0D0 )  THEN
               SCALE = SCALE + (CORE(NI)*CORE(NJ)/RIJ)*FN1(NJ,IG)*EXP ( MAX ( -30.0D0, -FN2(NJ,IG)*(RIJ-FN3(NJ,IG))**2 ) )
            END IF
         END DO
      END IF
      IF ( HAMILTONIAN == "PDDG" ) THEN 
          ZAF = CORE(NI)/(CORE(NI) + CORE(NJ))
          ZBF = CORE(NJ)/(CORE(NI) + CORE(NJ))
          DO IK= 1,2   
             DO JK = 1,2
                D = RIJ - PDDGE(NI,IK) - PDDGE(NJ,JK)
                AX = PDDG_EXPONENT * D * D
                SCALE = SCALE + (ZAF * PDDGC(NI,IK) + ZBF * PDDGC(NJ,JK)) * EXP(-AX)
             END DO
          END DO
      END IF
      ENUC = ENUC + SCALE
   KR = KR + KI
   END SUBROUTINE ROTATEI


   SUBROUTINE SEINTC_INTEGRALAB ( R, A, B, SWITCH, F, DF )
   REAL*8, INTENT(IN)  :: R, A, B, SWITCH
   REAL*8, INTENT(OUT) :: F, DF
   REAL*8 :: IOFF, ION, IMODOFF, IMODON, IMODR, IR
   CALL SEINTC_TNAB ( CUTOFFB, A, B, IMODOFF, IOFF )
   IF ( R <= CUTONB ) THEN
      IR = 1.0D0 / SQRT ( ( R + A )**2 + B )
      CALL SEINTC_TNAB ( CUTONB, A, B, IMODON, ION )
      F = IR - ( IMODOFF - IMODON ) - ION
   ELSE
      CALL SEINTC_TNAB ( R, A, B, IMODR, IR )
      F = IMODR - IMODOFF
   END IF
   DF = - SWITCH * ( R + A ) * IR * IR * IR
   END SUBROUTINE SEINTC_INTEGRALAB


   SUBROUTINE SEINTC_INTEGRALB ( R, B, SWITCH, F, DF )
   REAL*8, INTENT(IN)  :: R, B, SWITCH
   REAL*8, INTENT(OUT) :: F, DF
   REAL*8 :: IOFF, ION, IMODOFF, IMODON, IMODR, IR
   CALL SEINTC_TNB ( CUTOFFB, B, IMODOFF, IOFF )
   IF ( R <= CUTONB ) THEN
      IR = 1.0D0 / SQRT ( R**2 + B )
      CALL SEINTC_TNB ( CUTONB, B, IMODON, ION )
      F = IR - ( IMODOFF - IMODON ) - ION
   ELSE
      CALL SEINTC_TNB ( R, B, IMODR, IR )
      F = IMODR - IMODOFF
   END IF
   DF = - SWITCH * R * IR * IR * IR
   END SUBROUTINE SEINTC_INTEGRALB


   SUBROUTINE SEINTC_TNAB ( R, A, B, IMOD, T0 )
   REAL*8, INTENT(IN)  :: R, A, B
   REAL*8, INTENT(OUT) :: IMOD, T0
   REAL*8 :: A2, B2, FACT, LOGFAC, R2, RA, RASQB, T2, T4, T6
   A2     = A * A
   B2     = B * B
   R2     = R * R
   RA     = R + A
   RASQB  = RA * RA + B
   FACT   = SQRT ( RASQB )
   LOGFAC = LOG ( RA + FACT )
   T0 = 1.0D0 / FACT
   T2 = - ( 2.0D0 * A2 + 2.0D0 * B + 4.0D0 * A * R + R2 ) * T0 + 2.0D0 * A * LOGFAC
   T4 = ( - 22.0D0 * A2 * A2 - 14.0D0 * A2 * B + 8.0D0 * B2 - 34.0D0 * A2 * A * R + 26.0D0 * A * B * R &
                     - 6.0D0 * A2 * R2 + 4.0D0 * B * R2 + 2.0D0 * A * R2 * R - R2 * R2 ) * T0 / 3.0D0 + &
                                                        2.0D0 * A * ( 2.0D0 * A2 - 3.0D0 * B ) * LOGFAC
   T6 = ( - 274.0D0 * A2 * A2 * A2 + 333.0D0 * A2 * A2 * B + 543.0D0 * A2 * B2 - 64.0D0 * B2 * B -                 &
            394.0D0 * A2 * A2 * A * R + 1207.0D0 * A2 * A * B * R - 289.0D0 * A * B2 * R - 60.0D0 * A2 * A2 * R2 + &
            223.0D0 * A2 * B * R2 - 32.0D0 * B2 * R2 + 20.0D0 * A2 * A * R2 * R - 43.0D0 * A * B * R2 * R -        &
            10.0D0 * A2 * R2 * R2 + 8.0D0 * B * R2 * R2 + 6.0D0 * A * R2 * R2 * R - 4.0D0 * R2 * R2 * R2 ) * T0 /  &
            20.0D0 + 3.0D0 * A * ( 8.0D0 * A2 * A2 - 40.0D0 * A2 * B + 15.0D0 * B2 ) * LOGFAC / 4.0D0
   IMOD = T0FAC * T0 + T2FAC * T2 + T4FAC * T4 + T6FAC * T6
   END SUBROUTINE SEINTC_TNAB

   
   SUBROUTINE SEINTC_TNB ( R, B, IMOD, T0 )
   REAL*8, INTENT(IN)  :: R, B
   REAL*8, INTENT(OUT) :: IMOD, T0
   REAL*8 :: B2, FACT, R2, RSQB, T2, T4, T6
   B2   = B * B
   R2   = R * R
   RSQB = R * R + B
   FACT = SQRT ( RSQB )
   T0 = 1.0D0 / FACT
   T2 = - ( 2.0D0 * B + R2 ) * T0
   T4 = ( 8.0D0 * B2 + 4.0D0 * B * R2 - R2 * R2 ) * T0 / 3.0D0
   T6 = ( - 64.0D0 * B2 * B - 32.0D0 * B2 * R2 + 8.0D0 * B * R2 * R2 - 4.0D0 * R2 * R2 * R2 ) * T0 / 20.0D0
   IMOD = T0FAC * T0 + T2FAC * T2 + T4FAC * T4 + T6FAC * T6
   END SUBROUTINE SEINTC_TNB

   
   SUBROUTINE FOCK_MATRIX_RESTRICTED ( DENCLS, FCKCLS )
   REAL*8, DIMENSION(1:NBASTR), INTENT(IN)  :: DENCLS
   REAL*8, DIMENSION(1:NBASTR), INTENT(OUT) :: FCKCLS
   INTEGER            :: I, IAB, IFIRST, ILAST, INDTEI, IP, IQM, I1, I2, J, JAB, JFIRST, &
                         JJ, JLAST, JP, JQM, KA, KFI, KFJ, L, M, NI, NJ
   REAL*8 :: PTPOP
   REAL*8, ALLOCATABLE, DIMENSION(:)     :: FDIAG
   REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PDIAG

   FCKCLS = 0.0D0
   ALLOCATE ( FDIAG(10*NAQM), PDIAG(1:4,1:4,1:NAQM) )
   DO IQM = 1,NAQM
      IFIRST = BFIRST(IQM)
      ILAST  = BLAST(IQM)
      NI     = ANINT( ATMCHG(IQM ) )
      IF ( NATORB(NI) == 0 ) CYCLE
      PTPOP  = 0.0D0

      ! . Hydrogen.
      IF ( NATORB(NI) <= 1 ) THEN
         CONTINUE
      ! . SP elements.
      ELSE
         PTPOP = DENCLS((ILAST*(ILAST+1))/2) + DENCLS(((ILAST-1)*(ILAST))/2) + DENCLS(((ILAST-2)*(ILAST-1))/2)
      END IF

      ! . ss diagonal terms.
      KA = (IFIRST*(IFIRST+1))/2
      FCKCLS(KA) = FCKCLS(KA) + ( 0.5D0 * DENCLS(KA) * GSS(NI) ) + PTPOP * ( GSP(NI) - 0.5D0 * HSP(NI) )
      IF ( NI < 3 ) CYCLE

      L = KA
      DO J = (IFIRST+1),ILAST
         M = L+IFIRST
         L = L+J
         ! . pp diagonal terms.
         FCKCLS(L) = FCKCLS(L) + DENCLS(KA) * (GSP(NI)-0.5D0*HSP(NI)) + 0.5D0 * DENCLS(L) * GPP(NI) + &
                     ( PTPOP - DENCLS(L) ) * GP2(NI) - 0.25D0 * (PTPOP-DENCLS(L)) * (GPP(NI)-GP2(NI))
         ! . sp terms.
         FCKCLS(M) = FCKCLS(M) + ( 2.0D0 * DENCLS(M) * HSP(NI) ) - 0.5D0 * DENCLS(M) * ( HSP(NI)+GSP(NI) )
      END DO

      ! . pp* terms.
      DO J = (IFIRST+1),(ILAST-1)
         DO L = (J+1),ILAST
            M = (L*(L-1))/2+J
            FCKCLS(M) = FCKCLS(M) + DENCLS(M) * (GPP(NI)-GP2(NI)) - 0.25D0 * DENCLS(M) * (GPP(NI)+GP2(NI))
         END DO
      END DO
   END DO

   INDTEI   = 1
   KFI      = 1
   FDIAG    = 0.0D0
   PDIAG    = 0.0D0
   DO IQM = 1,NAQM
      IFIRST = BFIRST(IQM)
      ILAST  = BLAST(IQM)
      IAB    = ILAST - IFIRST + 1
      NI     = ANINT( ATMCHG(IQM ) )
      IF ( NATORB(NI) == 0 ) CYCLE
      IF ( IAB == 1 ) THEN
         PDIAG(1,1,IQM) = DENCLS((IFIRST*(IFIRST+1))/2)
      ELSE
         I1 = ((IFIRST-1)*(IFIRST-2))/2
         IP = 0
         DO I = IFIRST,ILAST
            IP = IP + 1
            I1 = I1 + I - 1
            I2 = I1 + IFIRST - 1
            JP = 0
            DO J = IFIRST,I
               JP = JP + 1
               I2 = I2 + 1
               PDIAG(JP,IP,IQM) = DENCLS(I2)
               PDIAG(IP,JP,IQM) = DENCLS(I2)
            END DO
         END DO
      END IF
      KFJ = 1

      DO JQM = 1,(IQM-1)
         JFIRST = BFIRST(JQM)
         JLAST  = BLAST(JQM)
         JAB    = JLAST - JFIRST + 1
         NJ     = ANINT( ATMCHG(JQM) )
         IF ( NATORB(NJ) == 0 ) CYCLE

         CALL SEFOCK ( FDIAG(KFI:), FDIAG(KFJ:),   SETEI(  INDTEI:),     INDTEI, PDIAG(:,:,IQM), PDIAG(:,:,JQM) )

         ! . Increment the counter for FDIAG.
         KFJ = KFJ + (JAB*(JAB+1))/2

      END DO

      ! . Increment the counter for FDIAG.
      KFI = KFI + (IAB*(IAB+1))/2

   END DO

   ! . Add in the diagonal parts of the Fock matrix.
   JJ = 0
   DO IQM = 1,NAQM
      IFIRST = BFIRST(IQM)
      ILAST  = BLAST(IQM)
      I1     = ((IFIRST-1)*(IFIRST-2))/2
      NI     = ANINT( ATMCHG(IQM) )
      IF ( NATORB(NI) == 0 ) CYCLE
      ! . Loop over the atom basis functions.
      DO I = IFIRST,ILAST
         I1 = I1 + I - 1
         I2 = I1 + IFIRST - 1
         DO J = IFIRST,I
            I2 = I2 + 1
            JJ = JJ + 1
            FCKCLS(I2) = FCKCLS(I2) + FDIAG(JJ)
         END DO
      END DO
   END DO

   ! . Deallocate space.
   DEALLOCATE ( FDIAG, PDIAG )
   CONTAINS
      SUBROUTINE SEFOCK ( FII, FJJ, SETEI, KR, P2II, P2JJ )
      !     SEFCC1 constructs the two-centre two-electron part of the Fock
      !     matrix given :
      !
      !     F - the partial packed Fock matrix, FII - the atom ii part of the
      !     Fock matrix, FJJ - the atom jj part of the Fock matrix, W - the
      !     two-centre two-electron integrals and P - the density matrix.
      !
      !     Atom indices are IFIRST to ILAST for atom I and JFIRST to JLAST
      !     for atom J.
      !

      ! . Scalar argument declarations.
      INTEGER, INTENT(INOUT) :: KR

      ! . Array argument declarations.
      REAL*8, DIMENSION(:),    INTENT(INOUT) :: FII,  FJJ
      REAL*8, DIMENSION(1:16), INTENT(IN)    :: P2II, P2JJ
      REAL*8, DIMENSION(:),    INTENT(IN)    :: SETEI

      ! . Local scalars.
      INTEGER            :: I, IAC, IAC2, II, IJ, IJAC, J, JAC, JAC2, JJ, K, KK, KL, NTYPE
      REAL*8 :: SUM

      ! . Local arrays.
      REAL*8, DIMENSION(1:16) :: PAIJ

      ! . Establish NTYPE. If 1 both atoms in the diatomic have 4 orbitals,
      ! . if 2 one atom has 4 orbitals and the other 1 and if 3 both atoms
      ! . have single orbitals.
      IAC = ILAST - IFIRST + 1
      JAC = JLAST - JFIRST + 1
      IF ( ( IAC + JAC ) == 8 ) THEN
         NTYPE = 1
      ELSE IF ( ( IAC + JAC ) == 5 ) THEN
         NTYPE = 2
      ELSE
         IJ = (IFIRST*(IFIRST-1))/2 + JFIRST
         FCKCLS(IJ) = FCKCLS(IJ) - ( 0.5D0 * SETEI(1) * DENCLS(IJ) )
         FII(1)   = FII(1) + ( SETEI(1) * P2JJ(1) )
         FJJ(1)   = FJJ(1) + ( SETEI(1) * P2II(1) )
         KR = KR + 1
         RETURN
      END IF

      ! . Unpack off-diagonal spin density matrix elements.
      KK = ((IFIRST-1)*(IFIRST-2))/2
      DO I = IFIRST,ILAST
         KK = KK + I - 1
         JJ = ( I - IFIRST ) * JAC
         K  = KK + JFIRST - 1
         IF ( JAC == 1 ) THEN
            PAIJ(JJ+1) = 0.5D0 * DENCLS(K+1)
         ELSE
            DO J = JFIRST,JLAST
               JJ = JJ + 1
               K  = K  + 1
               PAIJ(JJ) = 0.5D0 * DENCLS(K)
            END DO
         ENDIF
      END DO

      ! . Construct the off-diagonal exchange contributions.
      IJAC = IAC * JAC
      II   = ((IFIRST-1)*(IFIRST-2))/2
      DO I = 1,IAC
         II = II + IFIRST + I - 2
         IJ = II + JFIRST - 1
         DO J = 1,JAC
            IJ  = IJ + 1
            SUM = 0.0D0
            DO KL = 1,IJAC
               SUM = SUM + ( SETEI( KIND2(KL,I,J,NTYPE) ) * PAIJ(KL) )
            END DO
            FCKCLS(IJ) = FCKCLS(IJ) - SUM
         END DO
      END DO

      ! . Construct the on-diagonal Coulombic contributions.
      ! . Atom I affected by atom J.
      JAC2 = JAC * JAC
      JJ   = 0
      DO I = 1,IAC
         DO J = 1,I
            JJ = JJ + 1
            IF ( JAC2 == 1 ) THEN
               FII(JJ) = FII(JJ) + ( SETEI( JIND2(1,I,J,2) ) * P2JJ(1) )
            ELSE
               SUM = 0.0D0
               DO KL = 1,JAC2
                  SUM = SUM + ( SETEI( JIND2(KL,I,J,NTYPE) ) * P2JJ(KL) )
               END DO
               FII(JJ) = FII(JJ) + SUM
            END IF
         END DO
      END DO

      ! . Atom J affected by atom I.
      IAC2 = IAC * IAC
      JJ   = 0
      DO I = 1,JAC
         DO J = 1,I
            JJ = JJ + 1
            IF ( IAC2 == 1 ) THEN
               FJJ(JJ) = FJJ(JJ) + ( SETEI( JIND3(I,J,1,2) ) * P2II(1) )
            ELSE
               SUM = 0.0D0
               DO KL = 1,IAC2
                  SUM = SUM + ( SETEI( JIND3(I,J,KL,NTYPE) ) * P2II(KL) )
               END DO
               FJJ(JJ) = FJJ(JJ) + SUM
            END IF
         END DO
      END DO

      ! . Increment KR counter to line up W matrix for the next pair of atoms.
      KR = KR + ( ( IAC*(IAC+1) * JAC*(JAC+1) ) / 4 )

      END SUBROUTINE SEFOCK

   END SUBROUTINE FOCK_MATRIX_RESTRICTED


   SUBROUTINE FOCK_MATRIX_UNRESTRICTED ( DENTOT, DENALP, DENBET, FCKALP, FCKBET )
   REAL*8, DIMENSION(1:NBASTR), INTENT(IN)  :: DENALP, DENBET, DENTOT
   REAL*8, DIMENSION(1:NBASTR), INTENT(OUT) :: FCKALP, FCKBET
   INTEGER            :: I, IAB, IFIRST, ILAST, INDTEI, IP, IQM, I1, I2, J, JAB, JFIRST, &
                         JJ, JLAST, JP, JQM, KA, KFI, KFJ, L, M, NI, NJ, PIINDTEI
   REAL*8 :: PAPOP, PBPOP, PTPOP
   REAL*8, ALLOCATABLE, DIMENSION(:)     :: FDIAG
   REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PDIAG

   FCKALP = 0.0D0
   FCKBET = 0.0D0
   ALLOCATE ( FDIAG(10*NAQM), PDIAG(1:4,1:4,1:NAQM) )
   DO IQM = 1,NAQM
      IFIRST = BFIRST(IQM)
      ILAST  = BLAST(IQM)
      NI     = ANINT( ATMCHG(IQM) )
      PAPOP  = 0.0D0
      PBPOP  = 0.0D0
      PTPOP  = 0.0D0

      ! . Hydrogen.
      IF ( NATORB(NI) <= 1 ) THEN
         CONTINUE
      ! . SP elements.
      ELSE
         PAPOP = DENALP((ILAST*(ILAST+1))/2) + DENALP(((ILAST-1)*(ILAST))/2) + DENALP(((ILAST-2)*(ILAST-1))/2)
         PBPOP = DENBET((ILAST*(ILAST+1))/2) + DENBET(((ILAST-1)*(ILAST))/2) + DENBET(((ILAST-2)*(ILAST-1))/2)
         PTPOP = DENTOT((ILAST*(ILAST+1))/2) + DENTOT(((ILAST-1)*(ILAST))/2) + DENTOT(((ILAST-2)*(ILAST-1))/2)
      END IF

      ! . ss diagonal terms.
      KA = (IFIRST*(IFIRST+1))/2
      FCKALP(KA) = FCKALP(KA) + DENBET(KA) * GSS(NI) + PTPOP * GSP(NI) - PAPOP * HSP(NI)
      FCKBET(KA) = FCKBET(KA) + DENALP(KA) * GSS(NI) + PTPOP * GSP(NI) - PBPOP * HSP(NI)
      IF ( NI < 3 ) CYCLE

      L = KA
      DO J = (IFIRST+1),ILAST
         M = L+IFIRST
         L = L+J
         ! . pp diagonal terms.
         FCKALP(L) = FCKALP(L) + DENTOT(KA) * GSP(NI) - DENALP(KA) * HSP(NI) + &
                                 DENBET(L)  * GPP(NI) + ( PTPOP - DENTOT(L) ) * GP2(NI) - &
                                 0.5D0 * ( PAPOP - DENALP(L) ) * ( GPP(NI) - GP2(NI) )
         FCKBET(L) = FCKBET(L) + DENTOT(KA) * GSP(NI) - DENBET(KA) * HSP(NI) + &
                                 DENALP(L)  * GPP(NI) + ( PTPOP - DENTOT(L) ) * GP2(NI) - &
                                 0.5D0 * ( PBPOP - DENBET(L) ) * ( GPP(NI) - GP2(NI) )
     ! . sp terms.
         FCKALP(M) = FCKALP(M) + ( 2.0D0 * DENTOT(M) * HSP(NI) ) - DENALP(M) * ( HSP(NI)+GSP(NI) )
         FCKBET(M) = FCKBET(M) + ( 2.0D0 * DENTOT(M) * HSP(NI) ) - DENBET(M) * ( HSP(NI)+GSP(NI) )
      END DO

      ! . pp* terms.
      DO J = (IFIRST+1),(ILAST-1)
         DO L = (J+1),ILAST
            M = (L*(L-1))/2+J
            FCKALP(M) = FCKALP(M) + DENTOT(M) * (GPP(NI)-GP2(NI)) - 0.5D0 * DENALP(M) * (GPP(NI)+GP2(NI))
            FCKBET(M) = FCKBET(M) + DENTOT(M) * (GPP(NI)-GP2(NI)) - 0.5D0 * DENBET(M) * (GPP(NI)+GP2(NI))
         END DO
      END DO
   END DO

   INDTEI   = 1
   KFI      = 1
   PIINDTEI = 1
   FDIAG    = 0.0D0
   PDIAG    = 0.0D0
   DO IQM = 1,NAQM
      IFIRST = BFIRST(IQM)
      ILAST  = BLAST(IQM)
      IAB    = ILAST - IFIRST + 1
      NI     = ANINT( ATMCHG(IQM) ) 
      IF ( IAB == 1 ) THEN
         PDIAG(1,1,IQM) = DENTOT((IFIRST*(IFIRST+1))/2)
      ELSE
         I1 = ((IFIRST-1)*(IFIRST-2))/2
         IP = 0
         DO I = IFIRST,ILAST
            IP = IP + 1
            I1 = I1 + I - 1
            I2 = I1 + IFIRST - 1
            JP = 0
            DO J = IFIRST,I
               JP = JP + 1
               I2 = I2 + 1
               PDIAG(JP,IP,IQM) = DENTOT(I2)
               PDIAG(IP,JP,IQM) = DENTOT(I2)
            END DO
         END DO
      END IF
      KFJ = 1

      DO JQM = 1,(IQM-1)
         JFIRST = BFIRST(JQM)
         JLAST  = BLAST(JQM)
         JAB    = JLAST - JFIRST + 1
         NJ     = ANINT( ATMCHG(JQM) )

         CALL SEFOCK ( FDIAG(KFI:), FDIAG(KFJ:),   SETEI(  INDTEI:),     INDTEI, PDIAG(:,:,IQM), PDIAG(:,:,JQM) )

         ! . Increment the counter for FDIAG.
         KFJ = KFJ + (JAB*(JAB+1))/2

      END DO

      ! . Increment the counter for FDIAG.
      KFI = KFI + (IAB*(IAB+1))/2

   END DO

   ! . Add in the diagonal parts of the Fock matrix.
   JJ = 0
   DO IQM = 1,NAQM
      IFIRST = BFIRST(IQM)
      ILAST  = BLAST(IQM)
      I1     = ((IFIRST-1)*(IFIRST-2))/2
      NI     = ANINT( ATMCHG(IQM) )

      ! . Loop over the atom basis functions.
      DO I = IFIRST,ILAST
         I1 = I1 + I - 1
         I2 = I1 + IFIRST - 1
         DO J = IFIRST,I
            I2 = I2 + 1
            JJ = JJ + 1
            FCKALP(I2) = FCKALP(I2) + FDIAG(JJ)
            FCKBET(I2) = FCKBET(I2) + FDIAG(JJ)
         END DO
      END DO
   END DO

   ! . Deallocate space.
   DEALLOCATE ( FDIAG, PDIAG )
   CONTAINS
      SUBROUTINE SEFOCK ( FII, FJJ, SETEI, KR, P2II, P2JJ )
      !     SEFCC1 constructs the two-centre two-electron part of the Fock
      !     matrix given :
      !
      !     F - the partial packed Fock matrix, FII - the atom ii part of the
      !     Fock matrix, FJJ - the atom jj part of the Fock matrix, W - the
      !     two-centre two-electron integrals and P - the density matrix.
      !
      !     Atom indices are IFIRST to ILAST for atom I and JFIRST to JLAST
      !     for atom J.
      !

      ! . Scalar argument declarations.
      INTEGER, INTENT(INOUT) :: KR

      ! . Array argument declarations.
      REAL*8, DIMENSION(:),     INTENT(INOUT) :: FII,  FJJ
      REAL*8, DIMENSION(1:16),  INTENT(IN)    :: P2II, P2JJ
      REAL*8, DIMENSION(:),     INTENT(IN)    :: SETEI

      ! . Local scalars.
      INTEGER            :: I, IAC, IAC2, II, IJ, IJAC, J, JAC, JAC2, JJ, K, KK, KL, NTYPE
      REAL*8 :: SUM, SUMA, SUMB

      ! . Local arrays.
      REAL*8, DIMENSION(1:16) :: PAIJ, PBIJ

      ! . Establish NTYPE. If 1 both atoms in the diatomic have 4 orbitals,
      ! . if 2 one atom has 4 orbitals and the other 1 and if 3 both atoms
      ! . have single orbitals.
      IAC = ILAST - IFIRST + 1
      JAC = JLAST - JFIRST + 1
      IF ( ( IAC + JAC ) == 8 ) THEN
         NTYPE = 1
      ELSE IF ( ( IAC + JAC ) == 5 ) THEN
         NTYPE = 2
      ELSE
         IJ = (IFIRST*(IFIRST-1))/2 + JFIRST
         FCKALP(IJ) = FCKALP(IJ) - ( SETEI(1) * DENALP(IJ) )
         FCKBET(IJ) = FCKBET(IJ) - ( SETEI(1) * DENBET(IJ) )
         FII(1)     = FII(1)     + ( SETEI(1) * P2JJ(1)    )
         FJJ(1)     = FJJ(1)     + ( SETEI(1) * P2II(1)    )
         KR = KR + 1
         RETURN
      END IF

      ! . Unpack off-diagonal spin density matrix elements.
      KK = ((IFIRST-1)*(IFIRST-2))/2
      DO I = IFIRST,ILAST
         KK = KK + I - 1
         JJ = ( I - IFIRST ) * JAC
         K  = KK + JFIRST - 1
         IF ( JAC == 1 ) THEN
            PAIJ(JJ+1) = DENALP(K+1)
            PBIJ(JJ+1) = DENBET(K+1)
         ELSE
            DO J = JFIRST,JLAST
               JJ = JJ + 1
               K  = K  + 1
               PAIJ(JJ) = DENALP(K)
               PBIJ(JJ) = DENBET(K)
            END DO
         ENDIF
      END DO

      ! . Construct the off-diagonal exchange contributions.
      IJAC = IAC * JAC
      II   = ((IFIRST-1)*(IFIRST-2))/2
      DO I = 1,IAC
         II = II + IFIRST + I - 2
         IJ = II + JFIRST - 1
         DO J = 1,JAC
            IJ   = IJ + 1
            SUMA = 0.0D0
            SUMB = 0.0D0
            DO KL = 1,IJAC
               SUMA = SUMA + ( SETEI( KIND2(KL,I,J,NTYPE) ) * PAIJ(KL) )
               SUMB = SUMB + ( SETEI( KIND2(KL,I,J,NTYPE) ) * PBIJ(KL) )
            END DO
            FCKALP(IJ) = FCKALP(IJ) - SUMA
            FCKBET(IJ) = FCKBET(IJ) - SUMB
         END DO
      END DO

      ! . Construct the on-diagonal Coulombic contributions.
      ! . Atom I affected by atom J.
      JAC2 = JAC * JAC
      JJ   = 0
      DO I = 1,IAC
         DO J = 1,I
            JJ = JJ + 1
            IF ( JAC2 == 1 ) THEN
               FII(JJ) = FII(JJ) + ( SETEI( JIND2(1,I,J,2) ) * P2JJ(1) )
            ELSE
               SUM = 0.0D0
               DO KL = 1,JAC2
                  SUM = SUM + ( SETEI( JIND2(KL,I,J,NTYPE) ) * P2JJ(KL) )
               END DO
               FII(JJ) = FII(JJ) + SUM
            END IF
         END DO
      END DO

      ! . Atom J affected by atom I.
      IAC2 = IAC * IAC
      JJ   = 0
      DO I = 1,JAC
         DO J = 1,I
            JJ = JJ + 1
            IF ( IAC2 == 1 ) THEN
               FJJ(JJ) = FJJ(JJ) + ( SETEI( JIND3(I,J,1,2) ) * P2II(1) )
            ELSE
               SUM = 0.0D0
               DO KL = 1,IAC2
                  SUM = SUM + ( SETEI( JIND3(I,J,KL,NTYPE) ) * P2II(KL) )
               END DO
               FJJ(JJ) = FJJ(JJ) + SUM
            END IF
         END DO
      END DO

      ! . Increment KR counter to line up W matrix for the next pair of atoms.
      KR = KR + ( ( IAC*(IAC+1) * JAC*(JAC+1) ) / 4 )

      END SUBROUTINE SEFOCK
   END SUBROUTINE FOCK_MATRIX_UNRESTRICTED


   SUBROUTINE DENSITY_READ ( FILE )
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE
   INTEGER :: NMATIN
   
   OPEN ( 999, FILE = FILE, ACTION = "READ", STATUS = "OLD", FORM = "UNFORMATTED" )
   READ( 999 ) NMATIN
   IF ( NMATIN == 2 ) THEN
       READ( 999 ) DENMATA(1:NBASTR)
       READ( 999 ) DENMATB(1:NBASTR)
   ELSE
       READ( 999 ) DENMAT(1:NBASTR)
   END IF
   CLOSE( 999 )
   IF ( NMATIN == 2 ) DENMAT = DENMATA + DENMATB
   END SUBROUTINE DENSITY_READ


   SUBROUTINE DENSITY_WRITE ( FILE )
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE
   
   OPEN( 999, FILE = FILE, ACTION = "WRITE", STATUS = "UNKNOWN", FORM = "UNFORMATTED" )
   IF ( ALLOCATED ( DENMATA ) .AND. ALLOCATED ( DENMATB ) ) THEN
      WRITE( 999 ) 2
      WRITE( 999 ) DENMATA(1:NBASTR)
      WRITE( 999 ) DENMATB(1:NBASTR)
   ELSE
      WRITE( 999 ) 1
      WRITE( 999 ) DENMAT(1:NBASTR)
   END IF
   CLOSE( 999 )
   END SUBROUTINE DENSITY_WRITE


   SUBROUTINE MOPAC_SCF ( EHF, MAXIT, QPRINT )
   LOGICAL, INTENT(IN) :: QPRINT
   INTEGER, INTENT(IN) :: MAXIT
   REAL*8, INTENT(OUT) :: EHF
   REAL*8, PARAMETER :: CVGTOL = 0.005D0, DMPTOL = 1.0D0, SHFTOL = 0.2D0
   INTEGER            :: ICALL, IDAMP, ITER
   LOGICAL            :: CVGED, DMPFLG, LEVSHF, QBETA, QOPEN, QUHF, XTPFLG
   REAL*8 :: DAMP, DEHF, DEHFOLD, DIFF, DMPSAV, EHF0, OCCNUM, OSHIFT, VSHIFT
   REAL*8, ALLOCATABLE, DIMENSION(:)   :: FOCK1, FOCK2, ORBENE1, ORBENE2, TMPMAT
   REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ERRMAT1, ERRMAT2, FCKDMP1, FCKDMP2, FCKOLD1, FCKOLD2, &
                                                      MORBS1, MORBS2

   EHF = 0.0D0
   IF ( ( NBASIS <= 0 ) .OR. ( ( NALPHA <= 0 ) .AND. ( NBETA <= 0 ) ) ) RETURN
   QUHF = QFORCEUHF .OR. ( NALPHA /= NBETA )
   QBETA = QUHF .AND. ( NBETA > 0 )
   QOPEN = ( NALPHA /= NBETA )
   IF ( QUHF ) THEN
      IF ( QFORCEUHF ) THEN
         IF ( .NOT. ALLOCATED ( DENMATA ) ) THEN
            ALLOCATE ( DENMATA(1:NBASTR) ) ; DENMATA = 0.5D0 * DENMAT
         END IF
         IF ( .NOT. ALLOCATED ( DENMATB ) ) THEN
            ALLOCATE ( DENMATB(1:NBASTR) ) ; DENMATB = 0.5D0 * DENMAT
         END IF
      END IF
      IF ( NBETA == 0 ) THEN
         DENMAT  = DENMATA
         DENMATB = 0.0D0
      END IF
   ELSE
      ALLOCATE ( DENMATA(1:NBASTR) ) ; DENMATA = DENMAT ; DEALLOCATE ( DENMAT )
   END IF
   IF ( ALLOCATED ( MATIND ) ) DEALLOCATE ( MATIND )
   IF ( ALLOCATED ( BCOEFF ) ) DEALLOCATE ( BCOEFF )
   ALLOCATE ( BCOEFF(1:NDIIS,1:NDIIS), MATIND(1:NDIIS) )
   ICALL = 0
   IDAMP = 0
   DEHF  = 0.0D0
   DIFF  = 0.0D0
   EHF0  = 0.0D0
   DAMP   = DAMPF
   OSHIFT = SHIFTO
   VSHIFT = SHIFTV
   IF ( QOPEN ) THEN
       LEVSHF = ( OSHIFT > 0.0D0 ) .OR. ( VSHIFT > 0.0D0 )
   ELSE
       LEVSHF = ( VSHIFT > 0.0D0 )
   END IF
   DMPFLG = ( DAMP > DMPTOL )
   XTPFLG = ( .NOT. DMPFLG ) .AND. ( .NOT. LEVSHF )
   IF ( QUHF ) THEN
       OCCNUM = 1.0D0
   ELSE
       OCCNUM = 2.0D0
   END IF
   ALLOCATE ( FOCK1(1:NBASTR), MORBS1(1:NBASIS,1:NBASIS), ORBENE1(1:NBASIS), &
                   ERRMAT1(1:NBASTR,1:NDIIS), FCKDMP1(1:NBASTR,1:2), FCKOLD1(1:NBASTR,1:NDIIS) )
   IF ( QUHF ) THEN
      ALLOCATE ( FOCK2(1:NBASTR), MORBS2(1:NBASIS,1:NBASIS), ORBENE2(1:NBASIS), &
                   ERRMAT2(1:NBASTR,1:NDIIS), FCKDMP2(1:NBASTR,1:2), FCKOLD2(1:NBASTR,1:NDIIS) )
   END IF
   ALLOCATE ( TMPMAT(1:NBASTR) )
   IF( QPRINT ) THEN
       IF ( QUHF ) THEN
           WRITE( *, "(A)" ) "* Spin-Unrestricted SCF Calculation"
       ELSE
           WRITE( *, "(A)" ) "* Spin-Restricted SCF Calculation"
       END IF
       WRITE( *, "(A8,4A18)" ) "Cycle", "Energy", "Energy Change", "Density Change", "Damp/Shift"
       WRITE( *, "(80('-'))" )
   END IF

   ! . Perform 0SCF calculation
   IF( MAXIT <= 0 ) THEN
       IF ( QUHF ) THEN
           CALL FOCK_MATRIX_UNRESTRICTED ( DENMAT, DENMATA, DENMATB, FOCK1, FOCK2 )
       ELSE
           CALL FOCK_MATRIX_RESTRICTED ( DENMATA, FOCK1 )
       END IF
       FOCK1 = FOCK1 + HCORE
       IF ( QBETA ) FOCK2 = FOCK2 + HCORE
       EHF0 = EHF
       EHF  = 0.5D0 * ( SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATA, HCORE ) + &
                         SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATA, FOCK1 ) )
       IF ( QBETA ) THEN
           EHF = EHF + 0.5D0 * ( SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATB, HCORE ) + &
                                  SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATB, FOCK2 ) )
       END IF
       DEHFOLD = DEHF
       DEHF    = EHF - EHF0
       CVGED   = .TRUE.
       IF( QPRINT ) THEN
           IF ( DMPFLG ) THEN
               WRITE ( *, "(I8,4E18.8)" ) -1, EHF, DEHF, DIFF, DAMP
           ELSE IF ( LEVSHF ) THEN
               WRITE ( *, "(I8,4E18.8)" ) -1, EHF, DEHF, DIFF, MAX ( OSHIFT, VSHIFT )
           ELSE
               WRITE ( *, "(I8,4E18.8)" ) -1, EHF, DEHF, DIFF, 0.0D0
           END IF
       END IF
   END IF

   DO ITER = 1, MAXIT
       IF ( QUHF ) THEN
           CALL FOCK_MATRIX_UNRESTRICTED ( DENMAT, DENMATA, DENMATB, FOCK1, FOCK2 )
       ELSE
           CALL FOCK_MATRIX_RESTRICTED ( DENMATA, FOCK1 )
       END IF
       FOCK1 = FOCK1 + HCORE
       IF ( QBETA ) FOCK2 = FOCK2 + HCORE
       EHF0 = EHF
       EHF  = 0.5D0 * ( SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATA, HCORE ) + &
                         SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATA, FOCK1 ) )
       IF ( QBETA ) THEN
           EHF = EHF + 0.5D0 * ( SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATB, HCORE ) + &
                                  SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATB, FOCK2 ) )
       END IF
       DEHFOLD = DEHF
       DEHF    = EHF - EHF0
       IF( QPRINT ) THEN
           IF ( DMPFLG ) THEN
               WRITE ( *, "(I8,4E18.8)" ) ITER, EHF, DEHF, DIFF, DAMP
           ELSE IF ( LEVSHF ) THEN
               WRITE ( *, "(I8,4E18.8)" ) ITER, EHF, DEHF, DIFF, MAX ( OSHIFT, VSHIFT )
           ELSE
               WRITE ( *, "(I8,4E18.8)" ) ITER, EHF, DEHF, DIFF, 0.0D0
           END IF
       END IF
       CVGED = ( DEHF < CVGTOL ) .AND. ( DIFF < ACURCY ) .AND. ( ITER > 1 )
       IF ( CVGED ) EXIT
       IF ( DMPFLG .OR. LEVSHF ) THEN
           XTPFLG = ( ( ABS ( DEHF ) < CVGTOL ) .AND. ( DAMP < DMPTOL ) .AND. ( MAX ( OSHIFT, VSHIFT ) < SHFTOL ) ) .OR. &
                    ( ( ABS ( DEHF ) < CVGTOL ) .AND. ( DIFF < ACURCY ) )
           IF ( XTPFLG ) THEN
               DMPFLG = .FALSE.
               LEVSHF = .FALSE.
               OSHIFT = 0.0D0
              VSHIFT = 0.0D0
           END IF
       ELSE
           IF ( DEHF >= CVGTOL ) THEN
               DAMP   = DMPTOL
               DMPFLG = .TRUE.
               LEVSHF = .FALSE.
               IDAMP  = 0
               XTPFLG = .FALSE.
          END IF
       END IF
       IF ( LEVSHF .AND. ( ITER >  1 ) ) THEN
           IF ( QOPEN ) THEN
               CALL DENSITY_CALCULATE ( TMPMAT, MORBS1(1:NBASIS,NBETA+1:NALPHA), NALPHA-NBETA, 1.0D0 )
                FOCK1 = FOCK1 + OSHIFT * TMPMAT
           END IF
           CALL DENSITY_CALCULATE ( TMPMAT, MORBS1(1:NBASIS,NALPHA+1:NBASIS), NBASIS-NALPHA, 1.0D0 )
           FOCK1 = FOCK1 + VSHIFT * TMPMAT
           IF ( QBETA ) THEN
               CALL DENSITY_CALCULATE ( TMPMAT, MORBS2(1:NBASIS,NALPHA+1:NBASIS), NBASIS-NBETA, 1.0D0 )
               FOCK2 = FOCK2 + VSHIFT * TMPMAT
            END IF
       ELSE IF ( DMPFLG ) THEN
           IF ( ITER > 2 ) CALL GET_DAMP_FACTOR ( ITER, DEHF, DEHFOLD, DAMP )
           IDAMP = IDAMP + 1
           CALL DAVIDSON_DAMPING ( IDAMP, IEXTRP, DAMP, DAMP0, DMPSAV, FOCK1, FCKDMP1(:,1), FCKDMP1(:,2) )
           IF ( QBETA ) CALL DAVIDSON_DAMPING ( IDAMP, IEXTRP, DAMP, DAMP0, DMPSAV, FOCK2, FCKDMP2(:,1), FCKDMP2(:,2) )
       ELSE IF ( XTPFLG ) THEN
           IF ( QBETA ) THEN
               CALL DIIS ( ICALL, DENMATA, ERRMAT1, FOCK1, FCKOLD1, &
                                  DENMATB, ERRMAT2, FOCK2, FCKOLD2  )
           ELSE
               CALL DIIS ( ICALL, DENMATA, ERRMAT1, FOCK1, FCKOLD1 )
           END IF
       END IF
       CALL SYMMETRIC_UPPER ( FOCK1, ORBENE1, MORBS1 )
       IF ( LEVSHF .AND. ( ITER > 1 ) ) THEN
          IF ( QOPEN ) ORBENE1(NBETA+1:NALPHA) = ORBENE1(NBETA+1:NALPHA) - OSHIFT
          ORBENE1(NALPHA+1:NBASIS) = ORBENE1(NALPHA+1:NBASIS) - VSHIFT
       END IF
       TMPMAT = DENMATA
       CALL DENSITY_CALCULATE ( DENMATA, MORBS1(1:NBASIS,1:NALPHA), NALPHA, OCCNUM )
       DIFF = MAXVAL ( ABS ( DENMATA - TMPMAT ) )
       IF ( QUHF ) DENMAT = DENMATA
       IF ( QBETA ) THEN
           CALL SYMMETRIC_UPPER ( FOCK2, ORBENE2, MORBS2 )
           IF ( LEVSHF .AND. ( ITER > 1 ) ) ORBENE2(NBETA+1:NBASIS) = ORBENE2(NBETA+1:NBASIS) - VSHIFT
           TMPMAT = DENMATB
           CALL DENSITY_CALCULATE ( DENMATB, MORBS2(1:NBASIS,1:NBETA), NBETA, 1.0D0 )
           DIFF = MAX ( DIFF, MAXVAL ( ABS ( DENMATB - TMPMAT ) ) )
           DENMAT = DENMAT + DENMATB
       END IF
       IF ( CVGED ) EXIT
   END DO
   IF( QPRINT ) WRITE( *, "(80('-'))" )
   IF ( .NOT. CVGED ) THEN
       WRITE( *, *) ">> Excessive number of SCF iterations"
       IF ( QSCFBOMB ) STOP
   END IF
   DEALLOCATE ( FOCK1, ERRMAT1, FCKDMP1, FCKOLD1 )
   IF ( QBETA ) DEALLOCATE ( FOCK2, ERRMAT2, FCKDMP2, FCKOLD2 )
   IF ( QPRINT .AND. QUHF ) CALL SPIN_VALUES
   DEALLOCATE ( MORBS1, ORBENE1, TMPMAT )
   IF ( QBETA ) DEALLOCATE ( MORBS2, ORBENE2 )
   IF ( .NOT. QUHF ) THEN
       ALLOCATE ( DENMAT(1:NBASTR) ) ; DENMAT = DENMATA ; DEALLOCATE ( DENMATA )
   END IF
   END SUBROUTINE MOPAC_SCF


   SUBROUTINE MOPAC_SCF_INITIALIZE
   IEXTRP    =  15
   NDIIS     =  19
   NDIISP    = NDIIS + 1
   QFORCEUHF = .FALSE.
   ACURCY    = 1.0D-8
   DAMPF     = 0.0D0
   DAMP0     = 0.0D0
   SHIFTO    = 0.0D0
   SHIFTV    = 0.0D0
   END SUBROUTINE MOPAC_SCF_INITIALIZE


   SUBROUTINE MOPAC_SCF_OPTIONS ( ACCURACY, DAMP_FACTOR, DIIS_MATRICES, EXTRAPOLATION, FORCE_UHF, &
                                  OPENSHELL_SHIFT, VIRTUAL_SHIFT, ZERO_DAMP_FACTOR )
   INTEGER,            INTENT(IN), OPTIONAL :: DIIS_MATRICES, EXTRAPOLATION
   LOGICAL,            INTENT(IN), OPTIONAL :: FORCE_UHF
   REAL*8, INTENT(IN), OPTIONAL :: ACCURACY, DAMP_FACTOR, OPENSHELL_SHIFT, VIRTUAL_SHIFT, ZERO_DAMP_FACTOR

   IF ( PRESENT ( ACCURACY         ) ) ACURCY    = ACCURACY
   IF ( PRESENT ( DAMP_FACTOR      ) ) DAMPF     = DAMP_FACTOR
   IF ( PRESENT ( DIIS_MATRICES    ) ) NDIIS     = DIIS_MATRICES
   IF ( PRESENT ( EXTRAPOLATION    ) ) IEXTRP    = EXTRAPOLATION
   IF ( PRESENT ( FORCE_UHF        ) ) QFORCEUHF = FORCE_UHF
   IF ( PRESENT ( OPENSHELL_SHIFT  ) ) SHIFTO    = OPENSHELL_SHIFT
   IF ( PRESENT ( VIRTUAL_SHIFT    ) ) SHIFTV    = VIRTUAL_SHIFT
   IF ( PRESENT ( ZERO_DAMP_FACTOR ) ) DAMP0     = ZERO_DAMP_FACTOR   
   NDIISP = NDIIS + 1
   IF ( ( SHIFTO >  0.0D0 ) .AND. ( SHIFTV <= 0.0D0 ) ) SHIFTV = SHIFTO
   IF ( ( SHIFTO <= 0.0D0 ) .AND. ( SHIFTV >  0.0D0 ) ) SHIFTO = SHIFTV
   WRITE( *, "(A)"   ) "* Options ----------------------------"
   WRITE( *, "(A,E14.8)" ) "Accuracy Tolerance    : ", ACURCY   
   WRITE( *, "(A,L14)"   ) "Force UHF             : ", QFORCEUHF
   WRITE( *, "(A,I14)"   ) "Extrapolation Freq.   : ", IEXTRP   
   WRITE( *, "(A,E14.8)" ) "Damping Factor        : ", DAMPF    
   WRITE( *, "(A,E14.8)" ) "Zero Damping Factor   : ", DAMP0    
   WRITE( *, "(A,E14.8)" ) "Open-Shell Level Shift: ", SHIFTO   
   WRITE( *, "(A,E14.8)" ) "Virtual Level Shift   : ", SHIFTV   
   WRITE( *, "(A,I14)"   ) "No. of DIIS Matrices  : ", NDIIS    
   WRITE( *, "(38('-'))"   )
   END SUBROUTINE MOPAC_SCF_OPTIONS


   SUBROUTINE DAVIDSON_DAMPING ( ITER, IEXTRP, DAMP, DAMP0, DMPSAV, D0, D1, D2 )
   INTEGER,            INTENT(IN)    :: IEXTRP, ITER
   REAL*8, INTENT(IN)    :: DAMP, DAMP0
   REAL*8, INTENT(INOUT) :: DMPSAV
   REAL*8, DIMENSION(1:NBASTR), INTENT(INOUT) :: D0, D1, D2
   INTEGER            :: N
   REAL*8 :: CUTOFF, DUM11, DUM21, DUM22, EPS
   REAL*8, PARAMETER :: TOL = 0.01D0

   ! . On the first iteration save the first matrix only.
   IF ( ITER <= 1 ) THEN
      D1 = D0
   ELSE
      ! . Perform the damping without extrapolation for ITER > 1.
      IF ( IEXTRP <= 0 ) THEN
     D0 = ( D0 + DAMP * D1 ) / ( 1.0D0 + DAMP )
     D1 = D0
      ELSE
     ! . Perform the damping for ITER = 2.
     IF ( ITER == 2 ) THEN
            D0 = ( D0 + DAMP * D1 ) / ( 1.0D0 + DAMP )
            D2 = D1
            D1 = D0
            DMPSAV = DAMP
     ELSE
            ! . Check whether any extrapolation needs to be performed.
            N = ITER - ( ITER / IEXTRP ) * IEXTRP
            ! . Perform damping without extrapolation for ITER > 2.
            IF ( N /= 0 ) THEN
               D0 = ( D0  + DAMP * D1 ) / ( 1.0D0 + DAMP )
               D2 = D1
               D1 = D0
               DMPSAV = DAMP
            ELSE
               ! . Perform damping and extrapolation for ITER > 2.
               CUTOFF = 0.5D0 - DAMP0
               IF ( CUTOFF >= 1.0D0 ) CUTOFF = 0.99D0
               IF ( CUTOFF <= 0.0D0 ) CUTOFF = TOL
               D0    = ( D0 + DMPSAV * D1 ) / ( 1.0D0 + DMPSAV )
               DUM11 = DOT_PRODUCT ( ( D0 - D1 ), ( D0 - D1 ) )
               DUM21 = DOT_PRODUCT ( ( D1 - D2 ), ( D0 - D1 ) )
               DUM22 = DOT_PRODUCT ( ( D1 - D2 ), ( D1 - D2 ) )
               EPS   = ( DUM11 - DUM21 ) / ( DUM21 - DUM22 )
               IF ( DUM21 * DUM21 < 0.5D0 * DUM11 * DUM22 ) EPS = 0.0D0
               IF ( EPS > CUTOFF ) EPS = CUTOFF
               D0 = ( D0 - EPS * D1 ) / ( 1.0D0 - EPS )
               D2 = D1
               D1 = D0
               DMPSAV = -EPS * ( 1.0D0 + DMPSAV ) + DMPSAV
            END IF
     END IF
      END IF
   END IF
   END SUBROUTINE DAVIDSON_DAMPING


   SUBROUTINE DIIS ( NCALLS, DENMATA, ERRMATA, FCKNEWA, FCKOLDA, &
                             DENMATB, ERRMATB, FCKNEWB, FCKOLDB  )
   INTEGER, INTENT(INOUT) :: NCALLS
   REAL*8, DIMENSION(1:NBASTR),         INTENT(IN)    :: DENMATA
   REAL*8, DIMENSION(1:NBASTR),         INTENT(INOUT) :: FCKNEWA
   REAL*8, DIMENSION(1:NBASTR,1:NDIIS), INTENT(INOUT) :: ERRMATA, FCKOLDA
   REAL*8, DIMENSION(1:NBASTR),         INTENT(IN),    OPTIONAL :: DENMATB
   REAL*8, DIMENSION(1:NBASTR),         INTENT(INOUT), OPTIONAL :: FCKNEWB
   REAL*8, DIMENSION(1:NBASTR,1:NDIIS), INTENT(INOUT), OPTIONAL :: ERRMATB, FCKOLDB
   INTEGER :: I, IERR, II, IMAT, J, JJ, JMAT, NDIM, NEWEST, START, STOP
   LOGICAL :: QBAD
   REAL*8, DIMENSION(1:NDIISP,1:NDIISP) :: A
   REAL*8, DIMENSION(1:NDIISP)          :: B

   NCALLS = NCALLS + 1
   IF ( NCALLS == 1 ) THEN
      MATNUM = 0
      DO I = 1,NDIIS
         MATIND(I) = I
      END DO
   END IF
   IF ( MATNUM < NDIIS ) THEN
      MATNUM = MATNUM + 1
      NEWEST = MATIND(MATNUM)
   ELSE
      MATNUM = NDIIS
      NEWEST = MATIND(1)
      DO I = 1,(NDIIS-1)
         MATIND(I) = MATIND(I+1)
      END DO
      MATIND(NDIIS) = NEWEST
   END IF
   FCKOLDA(1:NBASTR,NEWEST) = FCKNEWA
   CALL SYMMETRIC_MATRIX_MULTIPLY_2 ( FCKNEWA, DENMATA, ERRMATA(1,NEWEST),   1.0D0, .TRUE.  )
   CALL SYMMETRIC_MATRIX_MULTIPLY_2 ( DENMATA, FCKNEWA, ERRMATA(1,NEWEST), - 1.0D0, .FALSE. )
   DO I = 1,MATNUM
      IMAT = MATIND(I)
      BCOEFF(IMAT,NEWEST) = 2.0D0 * DOT_PRODUCT ( ERRMATA(:,IMAT), ERRMATA(:,NEWEST) )
      BCOEFF(NEWEST,IMAT) = BCOEFF(IMAT,NEWEST)
   END DO
   IF ( PRESENT ( DENMATB ) ) THEN
      FCKOLDB(1:NBASTR,NEWEST) = FCKNEWB
      CALL SYMMETRIC_MATRIX_MULTIPLY_2 ( FCKNEWB, DENMATB, ERRMATB(1,NEWEST),   1.0D0, .TRUE.  )
      CALL SYMMETRIC_MATRIX_MULTIPLY_2 ( DENMATB, FCKNEWB, ERRMATB(1,NEWEST), - 1.0D0, .FALSE. )
      DO I = 1,MATNUM
         IMAT = MATIND(I)
         BCOEFF(IMAT,NEWEST) = BCOEFF(IMAT,NEWEST) + 2.0D0 * DOT_PRODUCT ( ERRMATB(:,IMAT), ERRMATB(:,NEWEST) )
         BCOEFF(NEWEST,IMAT) = BCOEFF(IMAT,NEWEST)
      END DO
   END IF
   FCKNEWA = FCKOLDA(1:NBASTR,NEWEST)
   IF ( MATNUM <= 1 ) GO TO 9999
   START = 1
   STOP  = MATNUM
   10 CONTINUE
   NDIM = STOP - START + 2
   A = 0.0D0
   B = 0.0D0
   II = 0
   DO I = START,STOP
      II   = II + 1
      IMAT = MATIND(I)
      JJ   = 0
      DO J = START,STOP
         JJ   = JJ + 1
         JMAT = MATIND(J)
         A(II,JJ) = BCOEFF(IMAT,JMAT)
      END DO
   END DO
   DO I = 1,(NDIM-1)
      A(NDIM,I) = - 1.0D0
      A(I,NDIM) = - 1.0D0
   END DO
   B(NDIM) = - 1.0D0
   IERR = 0
   CALL LINEAR_EQUATIONS_SOLVE ( NDIM, A(1:NDIISP,1:NDIM), NDIISP, B(1:NDIM), IERR )
   QBAD = ( IERR /= 0 )
   IF ( QBAD ) THEN
      START  = START + 1
      MATNUM = STOP - START + 1
      IF ( MATNUM <= 1 ) GO TO 20
      GO TO 10
   END IF
   20 CONTINUE
   IF ( START > 1 ) THEN
      II = 0
      DO I = START,STOP
         II = II + 1
         MATIND(II) = MATIND(I)
      END DO
      DO I = (MATNUM+1),NDIIS
         MATIND(I) = MATIND(I-1)+1
         IF ( MATIND(I) > NDIIS ) MATIND(I) = MATIND(I) - NDIIS
      END DO
      IF ( MATNUM <= 1 ) GOTO 9999
   END IF
   FCKNEWA = 0.0D0
   DO I = 1,MATNUM
      FCKNEWA = FCKNEWA + B(I) * FCKOLDA(:,MATIND(I))
   END DO
   IF ( PRESENT ( DENMATB ) ) THEN
      FCKNEWB = 0.0D0
      DO I = 1,MATNUM
         FCKNEWB = FCKNEWB + B(I) * FCKOLDB(:,MATIND(I))
      END DO
   END IF
   9999 CONTINUE
   END SUBROUTINE DIIS


   SUBROUTINE GET_DAMP_FACTOR ( ITER, DE, DEP, DAMP )
   INTEGER,            INTENT(IN)    :: ITER
   REAL*8, INTENT(IN)    :: DE, DEP
   REAL*8, INTENT(INOUT) :: DAMP
   REAL*8, PARAMETER :: DMPMAX = 256.0D0, FAC = 16.0D0
   REAL*8, SAVE :: DEAVG

   SELECT CASE ( ITER )
   CASE ( 1 )   ; DEAVG = 0.0D0
   CASE ( 2 )   ; DEAVG = ABS ( DE )
   CASE DEFAULT ; DEAVG = ( ABS ( DE ) + ABS ( DEP ) + ( 0.2D0 * DEAVG ) ) / 2.2D0
   END SELECT
   IF ( DE > 0.0D0 ) THEN
      IF ( DEP > 0.0D0 ) THEN
         DAMP = 4.0D0 * MAX ( DAMP, DEAVG )
         IF ( DE >= ( 4.0D0 * DEP ) ) THEN
            DAMP = FAC * MAX ( DAMP, DEAVG )
         ELSE IF ( DE <= ( 0.25D0 * DEP ) ) THEN
            DAMP = DAMP / FAC
         ELSE
            DAMP = ( DE / DEP )**2 * MAX ( DAMP, DEAVG )
         END IF
      ELSE
         DAMP = 4.0D0 * MAX ( DAMP, DEAVG )
         IF ( DE > 0.5D0 * DEAVG ) DAMP = DAMP * FAC
         IF ( ( DE - DEP ) < ( 0.2D0 * DEAVG ) ) THEN
            DAMP = DAMP / FAC
         END IF
      END IF
   ELSE
      IF ( DEP > 0.0D0 ) THEN
         DAMP = 4.0D0 * MAX ( DAMP, DEAVG )
         IF ( -DE > DEAVG ) DAMP = DAMP * FAC
         IF ( ( -DE + DEP ) < DEAVG ) THEN
            DAMP = DAMP / FAC
         END IF
      ELSE
         IF ( DE > DEP ) THEN
            IF ( DE <= ( 0.25D0 * DEP ) ) THEN
               DAMP = ( DE / DEP )**2 * MAX ( DAMP, DEAVG )
            ELSE
               DAMP = DAMP / FAC
            END IF
         ELSE
            IF ( ABS ( DE ) >= ( 2.0D0 * DEAVG ) ) THEN
               DAMP = FAC * MAX ( DAMP, DEAVG )
            ELSE
               IF ( ABS ( DE ) <= ( 0.5D0 * DEAVG ) ) DAMP = DAMP / FAC
            END IF
         END IF
      END IF
   END IF
   DAMP = MIN ( DAMP, DMPMAX )
   END SUBROUTINE GET_DAMP_FACTOR


   SUBROUTINE SPIN_VALUES
   REAL*8 :: SZ, S2
   INTEGER :: I, IFIRST, ILAST
   REAL*8, DIMENSION(1:NBASIS) :: BQTEMP

   SZ = 0.5D0 * REAL ( NALPHA - NBETA, 8 )
   S2 = SZ**2 + 0.5D0 * REAL ( NALPHA + NBETA, 8 ) - SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATA, DENMATB )
   WRITE ( *, "(A,E14.8)" ) "<Sz>:  ", SZ
   WRITE ( *, "(A,E14.8)" ) "<S^2>: ", S2
   WRITE ( *, "(A)" ) "* Spin densities:"
   DO I = 1,NBASIS
       BQTEMP(I) = DENMATA((I*(I+1))/2) - DENMATB((I*(I+1))/2)
   END DO
   DO I = 1, NAQM
       IFIRST = BFIRST(I)
       ILAST  = BLAST(I)
       WRITE( *, "(E14.8)" ) SUM ( BQTEMP(IFIRST:ILAST) )
   END DO
   END SUBROUTINE SPIN_VALUES


   FUNCTION SYMMETRIC_MATRIX_COLUMN ( COLUMN, N, MATRIX )
   INTEGER, INTENT(IN) :: COLUMN, N
   REAL*8, DIMENSION(:), INTENT(IN) :: MATRIX
   REAL*8, DIMENSION(1:N) :: SYMMETRIC_MATRIX_COLUMN
   INTEGER :: I

   SYMMETRIC_MATRIX_COLUMN(1:COLUMN) = MATRIX(BFINDEX(COLUMN)+1:BFINDEX(COLUMN)+COLUMN)
   DO I = (COLUMN+1),N
      SYMMETRIC_MATRIX_COLUMN(I) = MATRIX(BFINDEX(I)+COLUMN)
   END DO
   END FUNCTION SYMMETRIC_MATRIX_COLUMN


   FUNCTION SYMMETRIC_MATRIX_DIAGONAL ( N, MATRIX )
   INTEGER, INTENT(IN) :: N
   REAL*8, DIMENSION(:), INTENT(IN) :: MATRIX
   REAL*8, DIMENSION(1:N) :: SYMMETRIC_MATRIX_DIAGONAL
   INTEGER :: I

   DO I = 1,N
      SYMMETRIC_MATRIX_DIAGONAL(I) = MATRIX(BFINDEX(I)+I)
   END DO
   END FUNCTION SYMMETRIC_MATRIX_DIAGONAL


   FUNCTION SYMMETRIC_MATRIX_DOT_PRODUCT ( A, B )
   REAL*8 :: SYMMETRIC_MATRIX_DOT_PRODUCT
   REAL*8, DIMENSION(:), INTENT(IN) :: A, B
   SYMMETRIC_MATRIX_DOT_PRODUCT = 2.0D0 * DOT_PRODUCT ( A, B ) - &
                                           DOT_PRODUCT ( SYMMETRIC_MATRIX_DIAGONAL ( NBASIS, A ), &
                                                         SYMMETRIC_MATRIX_DIAGONAL ( NBASIS, B )  )
   END FUNCTION SYMMETRIC_MATRIX_DOT_PRODUCT


   SUBROUTINE SYMMETRIC_MATRIX_MULTIPLY_2 ( A, B, C, FACTOR, QINITIALIZE )
   LOGICAL,            INTENT(IN) :: QINITIALIZE
   REAL*8, INTENT(IN) :: FACTOR
   REAL*8, DIMENSION(1:NBASTR), INTENT(IN)    :: A, B
   REAL*8, DIMENSION(1:NBASTR), INTENT(INOUT) :: C
   INTEGER :: I, IK, K

   IF ( QINITIALIZE ) C = 0.0D0
   IF ( FACTOR == 0.0D0 ) RETURN
   IK = 0
   DO I = 1,NBASIS
      DO K = 1,I
         IK = IK + 1
         C(IK) = C(IK) + FACTOR * DOT_PRODUCT ( SYMMETRIC_MATRIX_COLUMN ( I, NBASIS, A ), &
                                            SYMMETRIC_MATRIX_COLUMN ( K, NBASIS, B )  )
      END DO
   END DO
   END SUBROUTINE SYMMETRIC_MATRIX_MULTIPLY_2


   SUBROUTINE MOPAC_GRADIENTS ( GRADIENT )
   REAL*8, DIMENSION(:,:), INTENT(INOUT) :: GRADIENT
   LOGICAL :: QAM1PM3
   QAM1PM3 = ( HAMILTONIAN == "AM1 " ) .OR. ( HAMILTONIAN == "PDDG" ) .OR. ( HAMILTONIAN == "PM3 " ) .OR. ( HAMILTONIAN == "RM1 " )
   CALL DERIVATIVES_QM
   CALL DERIVATIVES_QM_MM
   CONTAINS
      SUBROUTINE DERIVATIVES_QM
      INTEGER            :: I, IBASIS, IFIRST, IJ, ILAST, IQM, ISHELL, ISTART,    &
                            ISTOP, J, JBASIS, JFIRST, JLAST, JQM, JSHELL, JSTART, &
                            JSTOP, K, L
      LOGICAL            :: QUHF
      REAL*8                            :: BETAIJ, PFAC
      INTEGER, DIMENSION(1:2)           :: NDI
      REAL*8, DIMENSION(1:3,1:2)        :: CDI
      REAL*8, DIMENSION(1:3)            :: ENG
      REAL*8, DIMENSION(1:171)          :: PDIT, PDALP, PDBET, PDIA
      REAL*8, DIMENSION(1:3,1:NAQM)     :: GRADQM
      REAL*8, ALLOCATABLE, DIMENSION(:) ::   SX,   SY,   SZ

      GRADQM = 0.0D0
      ALLOCATE ( SX(1:NBASTR), SY(1:NBASTR), SZ(1:NBASTR) )
      CALL GAUSSIAN_OVERLAP_DERIVATIVES ( SX, SY, SZ )
      I = 0
      DO IBASIS = 1,(NBASIS)
         DO JBASIS = 1,IBASIS
            BETAIJ = 0.5D0 * ( BETA(IBASIS) + BETA(JBASIS) )
            I = I + 1
            SX(I) = BETAIJ * SX(I)
            SY(I) = BETAIJ * SY(I)
            SZ(I) = BETAIJ * SZ(I)
         END DO
      END DO
      SX = 2.0D0 * ANGSTROMS_TO_BOHRS * EV_TO_KJ * SX
      SY = 2.0D0 * ANGSTROMS_TO_BOHRS * EV_TO_KJ * SY
      SZ = 2.0D0 * ANGSTROMS_TO_BOHRS * EV_TO_KJ * SZ
      I = 0
      DO ISHELL = 1,NSHELL
         IQM    = KATOM(ISHELL)
         ISTART = KLOC(ISHELL)
         IF ( ISHELL /= NSHELL ) THEN
            ISTOP  = KLOC(ISHELL+1) - 1
         ELSE
            ISTOP = NBASIS
         END IF
         DO IBASIS = ISTART,ISTOP
            DO JSHELL = 1,ISHELL
               JQM    = KATOM(JSHELL)
               JSTART = KLOC(JSHELL)
               IF ( JSHELL /= NSHELL ) THEN
                  JSTOP  = MIN ( IBASIS, (KLOC(JSHELL+1)-1) )
               ELSE
                  JSTOP  = MIN ( IBASIS, NBASIS )
               END IF
               DO JBASIS = JSTART,JSTOP
                  I = I + 1
                  PFAC = DENMAT(I)
                  GRADQM(1,IQM) = GRADQM(1,IQM) + PFAC * SX(I)
                  GRADQM(2,IQM) = GRADQM(2,IQM) + PFAC * SY(I)
                  GRADQM(3,IQM) = GRADQM(3,IQM) + PFAC * SZ(I)
                  GRADQM(1,JQM) = GRADQM(1,JQM) - PFAC * SX(I)
                  GRADQM(2,JQM) = GRADQM(2,JQM) - PFAC * SY(I)
                  GRADQM(3,JQM) = GRADQM(3,JQM) - PFAC * SZ(I)
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE ( SX, SY, SZ )

      QUHF = ALLOCATED ( DENMATA ) .AND. ALLOCATED ( DENMATB )
      DO IQM = 1,NAQM
         IFIRST = BFIRST(IQM)
         ILAST  = BLAST(IQM)
         NDI(2) = ANINT( ATMCHG(IQM) )
         CDI(1:3,2) = ATMCRD(1:3,IQM)
         DO JQM = 1,(IQM-1)
            JFIRST = BFIRST(JQM)
            JLAST  = BLAST(JQM)
            NDI(1) = ANINT( ATMCHG(JQM) )
            CDI(1:3,1) = ATMCRD(1:3,JQM)
            IJ = 0
            DO I = JFIRST,JLAST
               K = I*(I-1)/2+JFIRST-1
               DO J = JFIRST,I
                  IJ = IJ+1
                  K  = K+1
                  PDIA(IJ) = DENMAT(K)
                  PDIT(IJ) = DENMAT(K)
               END DO
            END DO
            DO I = IFIRST,ILAST
               L = I*(I-1)/2
               K = L+JFIRST-1
               DO J = JFIRST,JLAST
                  IJ = IJ+1
                  K  = K+1
                  PDIA(IJ) = DENMAT(K)
                  PDIT(IJ) = DENMAT(K)
               END DO
               K = L+IFIRST-1
               DO L = IFIRST,I
                  IJ = IJ+1
                  K  = K+1
                  PDIA(IJ) = DENMAT(K)
                  PDIT(IJ) = DENMAT(K)
               END DO
            END DO
            IF ( QUHF ) THEN
                IJ = 0
                DO I = JFIRST,JLAST
                    K = I*(I-1)/2+JFIRST-1
                    DO J = JFIRST,I
                        IJ = IJ+1
                        K  = K+1
                        PDALP(IJ) = DENMATA(K)
                        PDBET(IJ) = DENMATB(K)
                     END DO
                END DO
                DO I = IFIRST,ILAST
                    L = I*(I-1)/2
                    K = L+JFIRST-1
                    DO J = JFIRST,JLAST
                        IJ = IJ+1
                        K  = K+1
                        PDALP(IJ) = DENMATA(K)
                        PDBET(IJ) = DENMATB(K)
                    END DO
                    K = L+IFIRST-1
                    DO L = IFIRST,I
                        IJ = IJ+1
                        K  = K+1
                        PDALP(IJ) = DENMATA(K)
                        PDBET(IJ) = DENMATB(K)
                    END DO
                 END DO
             ELSE
                PDALP(1:IJ) = 0.5D0 * PDIA(1:IJ)
                PDBET(1:IJ) = 0.5D0 * PDIA(1:IJ)
             END IF
             CALL ANALYT ( PDIT, PDIA, PDALP, PDBET, 1, CDI, NDI, JFIRST, JLAST, IFIRST, ILAST, ENG, QAM1PM3 )
             GRADQM(1:3,IQM) = GRADQM(1:3,IQM) + ENG(1:3)
             GRADQM(1:3,JQM) = GRADQM(1:3,JQM) - ENG(1:3)
         END DO
      END DO
      DO I = 1,NAQM
        GRADIENT(1:3,I) = GRADIENT(1:3,I) + GRADQM(1:3,I)
      END DO
      END SUBROUTINE DERIVATIVES_QM


      SUBROUTINE DERIVATIVES_QM_MM
      INTEGER            :: I, IATOM, II, IINT, J, JATOM, JJ, N1, NINTQ, QATOM
      REAL*8                      :: PDEN
      REAL*8, DIMENSION(1:3)      :: GVEC
      REAL*8, DIMENSION(1:10)     :: DENFAC

      DENFAC     = 2.0D0
      DENFAC(1)  = 1.0D0
      DENFAC(3)  = 1.0D0
      DENFAC(6)  = 1.0D0
      DENFAC(10) = 1.0D0
      IINT  = 0
      NINTQ = 0
      DO IATOM = 1, NAQM
          QATOM = ANINT( ATMCHG(IATOM) )
          DO JATOM = NAQM + 1, NATM
              IINT = IINT + 1
              N1   = BFIRST(IATOM)
              II   = ( N1 * ( N1 + 1 ) ) / 2
              PDEN = DENFAC(1) * DENMAT(II)
              GVEC(1) = - PDEN * DXE1BA(IINT)
              GVEC(2) = - PDEN * DYE1BA(IINT)
              GVEC(3) = - PDEN * DZE1BA(IINT)
              GRADIENT(1:3,IATOM) = GRADIENT(1:3,IATOM) + GVEC
              GRADIENT(1:3,JATOM) = GRADIENT(1:3,JATOM) - GVEC
              IF ( NATORB(QATOM) > 1 ) THEN
                 NINTQ = NINTQ + 1
                 JJ    = 1
                 DO I = (N1 + 1),BLAST(IATOM)
                     II = (I * (I - 1)) / 2 + N1 - 1
                     DO J = N1,I
                         II = II + 1
                         JJ = JJ + 1
                         PDEN = DENFAC(JJ) * DENMAT(II)
                         GVEC(1) = - PDEN * DXE1BB(JJ,NINTQ)
                         GVEC(2) = - PDEN * DYE1BB(JJ,NINTQ)
                         GVEC(3) = - PDEN * DZE1BB(JJ,NINTQ)
                         GRADIENT(1:3,IATOM) = GRADIENT(1:3,IATOM) + GVEC
                         GRADIENT(1:3,JATOM) = GRADIENT(1:3,JATOM) - GVEC
                     END DO
                 END DO
              END IF
          END DO
      END DO
      IF ( ALLOCATED ( DXE1BA ) ) DEALLOCATE ( DXE1BA )
      IF ( ALLOCATED ( DYE1BA ) ) DEALLOCATE ( DYE1BA )
      IF ( ALLOCATED ( DZE1BA ) ) DEALLOCATE ( DZE1BA )
      IF ( ALLOCATED ( DXE1BB ) ) DEALLOCATE ( DXE1BB )
      IF ( ALLOCATED ( DYE1BB ) ) DEALLOCATE ( DYE1BB )
      IF ( ALLOCATED ( DZE1BB ) ) DEALLOCATE ( DZE1BB )
      END SUBROUTINE DERIVATIVES_QM_MM

   END SUBROUTINE MOPAC_GRADIENTS


   SUBROUTINE ANALYT ( PSUMT, PSUMA, PALPHA, PBETA, NP, COORD, NAT, JFIRST, JLAST, IFIRST, ILAST, ENG, QAM1PM3 )
   INTEGER, INTENT(IN) :: IFIRST, ILAST, JFIRST, JLAST, NP
   LOGICAL, INTENT(IN) :: QAM1PM3
   INTEGER, DIMENSION(1:2),       INTENT(IN)  :: NAT
   REAL*8, DIMENSION(1:3,1:2),    INTENT(IN)  :: COORD
   REAL*8, DIMENSION(1:3),        INTENT(OUT) :: ENG
   REAL*8, DIMENSION(1:171),      INTENT(IN)  :: PSUMT
   REAL*8, DIMENSION(1:171,1:NP), INTENT(IN)  :: PALPHA, PBETA, PSUMA
   REAL*8, PARAMETER :: PDDG_EXPONENT = 10.0D0
   INTEGER :: I, IA, ID, IG, IK, ISP, IX, J, JA, JD, JK, K, KK, KL, L, LL, M, MK, ML, MN, N, NI, NJ, NK, NL
   REAL*8  :: AA, ANAM1, AX, BB, DFAC, DEL1, DM, F3, PFAC, RIJ, RR, R0, R2, TERMAA, TERMAB, TERMNC, ZAF, ZBF
   REAL*8, DIMENSION(1:4,1:2) :: CORINT
   REAL*8, DIMENSION(1:22)    :: DG, G
   REAL*8, DIMENSION(1:100)   :: DR
   REAL*8, DIMENSION(1:3)     :: EAA, EAB, ENUC

   EAA  = 0.0D0
   EAB  = 0.0D0
   ENG  = 0.0D0
   ENUC = 0.0D0

   ! . Initialize some counters.
   JD = JLAST - JFIRST + 1
   JA = 1
   ID = ILAST - IFIRST + 1 + JD
   IA = JD + 1

   ! . Get the interatomic distance.
   I  = 2
   NI = NAT(I)
   J  = 1
   NJ = NAT(J)
   R2 = SUM ( ( COORD(1:3,I) - COORD(1:3,J) )**2 )
   RIJ = SQRT ( R2 )
   R0  = ANGSTROMS_TO_BOHRS * RIJ
   RR  = R0 * R0

   ! . Calculate the two-electron integrals.
   CALL REPP ( NI, NJ, R0, G, CORINT )

   ! . Loop over the Cartesian components.
   DO IX = 1,3
      DEL1   = COORD(IX,I) - COORD(IX,J)
      TERMAA = 0.0D0
      TERMAB = 0.0D0
      ISP = 0
      CALL DELRI  ( DG, NI, NJ, R0, DEL1 )
      CALL DELMOL ( COORD, I, J, NI, NJ, IA, ID, JA, JD, IX, RIJ, DEL1, ISP, DG, DR, G )

      ! . Calculate the first derivative of the nuclear repulsion term.
      ! . Skip invalid terms.
      IF ( (RIJ.LT.1.0D0) .AND. (NATORB(NI)*NATORB(NJ).EQ.0) ) THEN
         TERMNC = 0.0D0
         GO TO 10
      END IF

      ! . Treat N-H and O-H interactions.
      IF ( (NI.EQ.1) .AND. ( (NJ.EQ.7) .OR. (NJ.EQ.8) ) ) THEN
         F3 = 1.0D0 + EXP(-ALP(1)*RIJ)+RIJ*EXP(-ALP(NJ)*RIJ)
         DFAC = DG(1)*F3-G(1)*(DEL1/RIJ)*(ALP(1)*EXP(-ALP(1)*RIJ)+(ALP(NJ)*RIJ-1.0D0)*EXP(-ALP(NJ)*RIJ))
         ELSE IF ( ((NI.EQ.7) .OR. (NI.EQ.8)) .AND. (NJ.EQ.1) ) THEN
            F3 = 1.0D0+EXP(-ALP(1)*RIJ)+RIJ*EXP(-ALP(NI)*RIJ)
            DFAC = DG(1)*F3-G(1)*(DEL1/RIJ)*(ALP(1)*EXP(-ALP(1)*RIJ)+(ALP(NI)*RIJ-1.0D0)*EXP(-ALP(NI)*RIJ))
      ! . Treat other interactions.
      ELSE
         F3 = 1.0D0+EXP(-ALP(NI)*RIJ)+EXP(-ALP(NJ)*RIJ)
         DFAC = DG(1)*F3-G(1)*(DEL1/RIJ)*(ALP(NI)*EXP(-ALP(NI)*RIJ)+ALP(NJ)*EXP(-ALP(NJ)*RIJ))
      END IF
      TERMNC = CORE(NI) * CORE(NJ) * DFAC

      ! . Calculate the extra core-core repulsion terms required for an AM1 calculation.
      IF ( QAM1PM3 ) THEN
         ANAM1 = 0.0D0
         DO IG = 1,4
            IF ( ABS ( FN1(NI,IG) ) .GT. 0.0D0 ) THEN
               ANAM1 = ANAM1+FN1(NI,IG)*(1.0D0/(RIJ*RIJ) + 2.0D0*FN2(NI,IG)*(RIJ-FN3(NI,IG))/RIJ)* &
                                               EXP(MAX(-30.0D0,-FN2(NI,IG)*(RIJ-FN3(NI,IG))**2))
            END IF
            IF ( ABS ( FN1(NJ,IG) ) .GT. 0.0D0 ) THEN
               ANAM1 = ANAM1+FN1(NJ,IG)*(1.0D0/(RIJ*RIJ) + 2.0D0*FN2(NJ,IG)*(RIJ-FN3(NJ,IG))/RIJ)* &
                                               EXP(MAX(-30.0D0,-FN2(NJ,IG)*(RIJ-FN3(NJ,IG))**2))
            END IF
         END DO
         ANAM1  = ANAM1 * CORE(NI) * CORE(NJ)
         TERMNC = TERMNC - ANAM1 * DEL1 / RIJ
      END IF

      ! . PDDG specific terms.
      IF ( HAMILTONIAN == "PDDG" ) THEN
        ANAM1 = 0.0D0
        ZAF = CORE(NI)/(CORE(NI) + CORE(NJ))
        ZBF = CORE(NJ)/(CORE(NI) + CORE(NJ))
        DO IK = 1,2
           DO JK = 1,2
              DM = RIJ - PDDGE(NI,IK) - PDDGE(NJ,JK)
              AX = PDDG_EXPONENT * DM * DM
              ANAM1 = ANAM1 + (ZAF * PDDGC(NI,IK) + ZBF * PDDGC(NJ,JK)) * &
                      2.0D0 * PDDG_EXPONENT * DM * EXP( -AX ) 
            END DO
         END DO
         TERMNC = TERMNC - ANAM1 * DEL1 / RIJ
      END IF

      10 CONTINUE

      ! . Calculate the core-electron attraction derivatives.
      ! . Atom core I affecting AOs on J.
      ISP = 0
      DO M = JA,JD
         BB = 1.0D0
         DO N = M,JD
            MN = M+N*(N-1)/2
            ISP = ISP+1
            TERMAB = TERMAB-BB*CORE(NI)*PSUMT(MN)*DR(ISP)
            BB = 2.0D0
         END DO
      END DO

      ! . Atom core J affecting AOs on I.
      K = MAX(JD-JA+1,1)
      K = (K*(K+1))/2
      ISP = -K+1
      DO M = IA,ID
         BB = 1.0D0
         DO N = M,ID
            MN = M+N*(N-1)/2
            ISP = ISP+K
            TERMAB = TERMAB-BB*CORE(NJ)*PSUMT(MN)*DR(ISP)
            BB = 2.0D0
         END DO
      END DO

      ! . Calculate the two electron terms.
      ISP = 0
      DO K = IA,ID
         AA = 1.0D0
         KK = (K*(K-1))/2
         DO L = K,ID
            LL = (L*(L-1))/2
            DO M = JA,JD
               BB = 1.0D0
               DO N = M,JD
                  ISP = ISP+1
                  KL = K+LL
                  MN = M+N*(N-1)/2

                  ! . Coulomb term.
                  PFAC   = SUM ( PSUMA(KL,1:NP) * PSUMA(MN,1:NP) )
                  TERMAA = TERMAA + AA * BB * PFAC * DR(ISP)
                  MK = M+KK
                  NK = N+KK
                  ML = M+LL
                  NL = N+LL

                  ! . Exchange term.
                  PFAC   = SUM ( PALPHA(MK,1:NP) * PALPHA(NL,1:NP) + PALPHA(NK,1:NP) * PALPHA(ML,1:NP) + &
                         PBETA (MK,1:NP) * PBETA (NL,1:NP) + PBETA (NK,1:NP) * PBETA (ML,1:NP)   )
                  TERMAA = TERMAA - 0.5D0 * AA * BB * DR(ISP) * PFAC
                  BB = 2.0D0

               END DO
            END DO
            AA = 2.0D0
         END DO
      END DO

      ! . Accumulate the results.
      EAA(IX)  = TERMAA
      EAB(IX)  = TERMAB
      ENUC(IX) = TERMNC

   END DO

   ENG = EV_TO_KJ * ( EAA + EAB + REAL ( NP, 8 ) * ENUC )

   CONTAINS

      SUBROUTINE DELMOL ( COORD, I, J, NI, NJ, IA, ID, JA, JD, IX, RIJ, DEL1, ISP, DG, DR, G )
      INTEGER,            INTENT(IN)    :: I, IA, ID, IX, J, JA, JD, NI, NJ
      INTEGER,            INTENT(INOUT) :: ISP
      REAL*8, INTENT(IN)    :: DEL1, RIJ
      REAL*8, DIMENSION(1:3,1:2), INTENT(IN)    :: COORD
      REAL*8, DIMENSION(1:22),    INTENT(IN)    :: DG, G
      REAL*8, DIMENSION(1:100),   INTENT(INOUT) :: DR
      INTEGER :: IB, JB, K, KK, L, LL, M, MM, N, NN
      REAL*8  :: TEMP1, TEMP2
      REAL*8, DIMENSION(1:3) :: TDX, TDY, TDZ, TX, TY, TZ

      IF ( ( NI.GT.1 ) .OR. ( NJ.GT.1) ) THEN
         CALL DER_ROTAT ( COORD, I, J, IX, RIJ, DEL1, 2, TDX, TDY, TDZ, TX, TY, TZ )
      ENDIF
      IB=MAX(IA,ID)
      JB=MAX(JA,JD)
      DO K=IA,IB
         KK=K-IA
         DO L=K,IB
            LL=L-IA
            DO M=JA,JB
               MM=M-JA
               DO N=M,JB
                  NN=N-JA
                  ISP=ISP+1
                  IF(NN.EQ.0)THEN
                     IF(LL.EQ.0) THEN
                        ! . (SS/SS)
                        DR(ISP)=DG(1)
                     ELSEIF(KK.EQ.0) THEN
                        ! . (SP/SS)
                        DR(ISP)=DG(2)*TX(LL)+G(2)*TDX(LL)
                     ELSE
                        ! . (PP/SS)
                        DR(ISP)=DG(3)*TX(KK)*TX(LL)+G(3)*(TDX(KK)*TX(LL)+TX(KK)*TDX(LL)) &
                               +DG(4)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))+G(4)*(TDY(KK)*TY(LL)+TY(KK)*TDY(LL) &
                               +TDZ(KK)*TZ(LL)+TZ(KK)*TDZ(LL))
                     ENDIF
                  ELSEIF(MM.EQ.0) THEN
                     IF(LL.EQ.0) THEN
                        ! . (SS/SP)
                        DR(ISP)=DG(5)*TX(NN)+G(5)*TDX(NN)
                     ELSEIF(KK.EQ.0) THEN
                        ! . (SP/SP)
                        DR(ISP)=DG(6)*TX(LL)*TX(NN)+G(6)*(TDX(LL)*TX(NN)+TX(LL)*TDX(NN)) &
                               +DG(7)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+G(7)*(TDY(LL)*TY(NN)+TY(LL)*TDY(NN) &
                               +TDZ(LL)*TZ(NN)+TZ(LL)*TDZ(NN))
                     ELSE
                        ! . (PP/SP)
                        DR(ISP)=DG(8)*TX(KK)*TX(LL)*TX(NN)+G(8)*(TDX(KK)*TX(LL)*TX(NN)+TX(KK)*TDX(LL)*TX(NN) &
                   +TX(KK)*TX(LL)*TDX(NN))+DG(9)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(NN) &
             +G(9)*((TDY(KK)*TY(LL)+TY(KK)*TDY(LL)+TDZ(KK)*TZ(LL)+TZ(KK)*TDZ(LL))*TX(NN) &
                   +(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TDX(NN))+DG(10)*(TX(KK)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN)) &
                     +TX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN)))+G(10)*(TDX(KK)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN)) &
                    +TDX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(KK)*(TDY(LL)*TY(NN)+TY(LL)*TDY(NN) &
                            +TDZ(LL)*TZ(NN)+TZ(LL)*TDZ(NN))+TX(LL)*(TDY(KK)*TY(NN)+TY(KK)*TDY(NN) &
                            +TDZ(KK)*TZ(NN)+TZ(KK)*TDZ(NN)))
                     ENDIF
                  ELSEIF(LL.EQ.0) THEN
                     ! . (SS/PP)
                     DR(ISP)=DG(11)*TX(MM)*TX(NN)+G(11)*(TDX(MM)*TX(NN)+TX(MM)*TDX(NN)) &
             +DG(12)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+G(12)*(TDY(MM)*TY(NN)+TY(MM)*TDY(NN) &
                    +TDZ(MM)*TZ(NN)+TZ(MM)*TDZ(NN))
                  ELSEIF(KK.EQ.0) THEN
                     ! . (SP/PP)
                     DR(ISP)=DG(13)*TX(LL)*TX(MM)*TX(NN)+G(13)*(TDX(LL)*TX(MM)*TX(NN)+TX(LL)*TDX(MM)*TX(NN) &
                    +TX(LL)*TX(MM)*TDX(NN))+DG(14)*TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN)) &
             +G(14)*(TDX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(LL)*(TDY(MM)*TY(NN)+TY(MM)*TDY(NN) &
                            +TDZ(MM)*TZ(NN)+TZ(MM)*TDZ(NN)))+DG(15)*(TY(LL)*(TY(MM)*TX(NN)+TY(NN)*TX(MM)) &
                     +TZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)*TX(MM)))+G(15)*(TDY(LL)*(TY(MM)*TX(NN)+TY(NN)*TX(MM)) &
                    +TDZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)*TX(MM))+TY(LL)*(TDY(MM)*TX(NN)+TY(MM)*TDX(NN) &
                            +TDY(NN)*TX(MM)+TY(NN)*TDX(MM))+TZ(LL)*(TDZ(MM)*TX(NN)+TZ(MM)*TDX(NN) &
                            +TDZ(NN)*TX(MM)+TZ(NN)*TDX(MM)))
                  ELSE
                     ! . (PP/PP)
                     DR(ISP)=DG(16)*TX(KK)*TX(LL)*TX(MM)*TX(NN)+G(16)*(TDX(KK)*TX(LL)*TX(MM)*TX(NN) &
                    +TX(KK)*TDX(LL)*TX(MM)*TX(NN)+TX(KK)*TX(LL)*TDX(MM)*TX(NN) &
                    +TX(KK)*TX(LL)*TX(MM)*TDX(NN))+DG(17)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(MM)*TX(NN) &
             +G(17)*((TDY(KK)*TY(LL)+TY(KK)*TDY(LL)+TDZ(KK)*TZ(LL)+TZ(KK)*TDZ(LL))*TX(MM)*TX(NN) &
                    +(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*(TDX(MM)*TX(NN)+TX(MM)*TDX(NN))) &
             +DG(18)*TX(KK)*TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+G(18)*((TDX(KK)*TX(LL)+TX(KK)*TDX(LL)) &
                       *(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(KK)*TX(LL)*(TDY(MM)*TY(NN)+TY(MM)*TDY(NN) &
                                   +TDZ(MM)*TZ(NN)+TZ(MM)*TDZ(NN)))
                     DR(ISP)=DR(ISP)+DG(19)*(TY(KK)*TY(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*TZ(MM)*TZ(NN)) &
             +G(19)*(TDY(KK)*TY(LL)*TY(MM)*TY(NN)+TY(KK)*TDY(LL)*TY(MM)*TY(NN)+TY(KK)*TY(LL)*TDY(MM)*TY(NN) &
                       +TY(KK)*TY(LL)*TY(MM)*TDY(NN)+TDZ(KK)*TZ(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TDZ(LL)*TZ(MM)*TZ(NN) &
                       +TZ(KK)*TZ(LL)*TDZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*TZ(MM)*TDZ(NN)) &
             +DG(20)*(TX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)*(TY(LL)*TY(MM)+TZ(LL)*TZ(MM))) &
                        +TX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM))))
                     ! . To avoid compiler difficulties the expression is divided.
                     TEMP1 = TDX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)*(TY(LL)*TY(MM)+TZ(LL)*TZ(MM)))  &
                            +TDX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM)))  &
                            +TX(KK)*(TDX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TDX(NN)*(TY(LL)*TY(MM)+TZ(LL)*TZ(MM))) &
                            +TX(LL)*(TDX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TDX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM)))
                     TEMP2 = TX(KK)*(TX(MM)*(TDY(LL)*TY(NN)+TY(LL)*TDY(NN)+TDZ(LL)*TZ(NN)+TZ(LL)*TDZ(NN)) &
                            +TX(NN)*(TDY(LL)*TY(MM)+TY(LL)*TDY(MM)+TDZ(LL)*TZ(MM)+TZ(LL)*TDZ(MM))) &
                            +TX(LL)*(TX(MM)*(TDY(KK)*TY(NN)+TY(KK)*TDY(NN)+TDZ(KK)*TZ(NN)+TZ(KK)*TDZ(NN)) &
                            +TX(NN)*(TDY(KK)*TY(MM)+TY(KK)*TDY(MM)+TDZ(KK)*TZ(MM)+TZ(KK)*TDZ(MM)))
                     DR(ISP)=DR(ISP)+G(20)*(TEMP1+TEMP2)
                     DR(ISP)=DR(ISP)+DG(21)*(TY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*TY(MM)*TY(NN)) &
             +G(21)*(TDY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TY(KK)*TDY(LL)*TZ(MM)*TZ(NN)+TY(KK)*TY(LL)*TDZ(MM)*TZ(NN) &
                       +TY(KK)*TY(LL)*TZ(MM)*TDZ(NN)+TDZ(KK)*TZ(LL)*TY(MM)*TY(NN)+TZ(KK)*TDZ(LL)*TY(MM)*TY(NN) &
                       +TZ(KK)*TZ(LL)*TDY(MM)*TY(NN)+TZ(KK)*TZ(LL)*TY(MM)*TDY(NN))
                     DR(ISP)=DR(ISP)+DG(22)*(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN)) &
             +G(22)*((TDY(KK)*TZ(LL)+TY(KK)*TDZ(LL)+TDZ(KK)*TY(LL)+TZ(KK)*TDY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN)) &
                       +(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(TDY(MM)*TZ(NN)+TY(MM)*TDZ(NN)+TDZ(MM)*TY(NN)+TZ(MM)*TDY(NN)))
                  END IF
               END DO
            END DO
         END DO
      END DO
      END SUBROUTINE DELMOL

      SUBROUTINE DELRI ( DG, NI, NJ, RR, DEL1 )
      INTEGER,            INTENT(IN) :: NI, NJ
      REAL*8, INTENT(IN) :: DEL1, RR
      REAL*8, DIMENSION(1:22), INTENT(OUT) :: DG
      REAL*8 :: ADD, ADE, ADQ, AED, AEE, AEQ, AQD, AQE, AQQ, DA, DB, &
                            DXDX, DXQXZ, DZDZ, DZE, DZQXX, DZQZZ, EDZ, EE,       &
                            EQXX, EQZZ, QA, QB, QXXDZ, QXXE, QXXQXX, QXXQYY,     &
                            QXXQZZ, QXYQXY, QXZDX, QXZQXZ, QZZDZ, QZZE, QZZQXX,  &
                            QZZQZZ, TERM

      TERM = ( ANGSTROMS_TO_BOHRS * ANGSTROMS_TO_BOHRS * AU_TO_EV * DEL1 ) / RR
      DA = DD(NI)
      DB = DD(NJ)
      QA = QQ(NI)
      QB = QQ(NJ)
      ! . Hydrogen - hydrogen.
      AEE=0.25D0*(1.0D0/AM(NI)+1.0D0/AM(NJ))**2
      EE    =-RR/(SQRT(RR**2+AEE))**3
      DG(1)=TERM*EE
      IF(NATORB(NI).LE.2.AND.NATORB(NJ).LE.2) RETURN
      IF(NATORB(NI).LE.2) GO TO 10
      ! . Heavy atom - hydrogen.
      ADE=0.25D0*(1.0D0/AD(NI)+1.0D0/AM(NJ))**2
      AQE=0.25D0*(1.0D0/AQ(NI)+1.0D0/AM(NJ))**2
      DZE   = (RR+DA)/(SQRT((RR+DA)**2+ADE))**3-(RR-DA)/(SQRT((RR-DA)**2+ADE))**3
      QZZE  =-(RR+2.0D0*QA)/(SQRT((RR+2.0D0*QA)**2+AQE))**3-(RR-2.0D0*QA)/(SQRT((RR-2.0D0*QA)**2+AQE))**3 &
             +(2.0D0*RR)/(SQRT(RR**2+AQE))**3
      QXXE  =-(2.0D0*RR)/(SQRT(RR**2+4.0D0*QA**2+AQE))**3+(2.0D0*RR)/(SQRT(RR**2+AQE))**3
      DG(2)=-(TERM*DZE)/2.0D0
      DG(3)=TERM*(EE+QZZE/4.0D0)
      DG(4)=TERM*(EE+QXXE/4.0D0)
      IF(NATORB(NJ).LE.2) RETURN
      ! . Hydrogen - Heavy atom.
   10 AED=0.25D0*(1.0D0/AM(NI)+1.0D0/AD(NJ))**2
      AEQ=0.25D0*(1.0D0/AM(NI)+1.0D0/AQ(NJ))**2
      EDZ   = (RR-DB)/(SQRT((RR-DB)**2+AED))**3-(RR+DB)/(SQRT((RR+DB)**2+AED))**3
      EQZZ  =-(RR-2.0D0*QB)/(SQRT((RR-2.0D0*QB)**2+AEQ))**3-(RR+2.0D0*QB)/(SQRT((RR+2.0D0*QB)**2+AEQ))**3 &
             +(2.0D0*RR)/(SQRT(RR**2+AEQ))**3
      EQXX  =-(2.0D0*RR)/(SQRT(RR**2+4.0D0*QB**2+AEQ))**3+(2.0D0*RR)/(SQRT(RR**2+AEQ))**3
      DG(5)=-(TERM*EDZ)/2.0D0
      DG(11)=TERM*(EE+EQZZ/4.0D0)
      DG(12)=TERM*(EE+EQXX/4.0D0)
      IF(NATORB(NI).LE.2) RETURN
      ! . Heavy atom - Heavy atom.
      ADD=0.25D0*(1.0D0/AD(NI)+1.0D0/AD(NJ))**2
      ADQ=0.25D0*(1.0D0/AD(NI)+1.0D0/AQ(NJ))**2
      AQD=0.25D0*(1.0D0/AQ(NI)+1.0D0/AD(NJ))**2
      AQQ=0.25D0*(1.0D0/AQ(NI)+1.0D0/AQ(NJ))**2
      DXDX  =-(2.0D0*RR)/(SQRT(RR**2+(DA-DB)**2+ADD))**3+(2.0D0*RR)/(SQRT(RR**2+(DA+DB)**2+ADD))**3
      DZDZ  =-(RR+DA-DB)/(SQRT((RR+DA-DB)**2+ADD))**3-(RR-DA+DB)/(SQRT((RR-DA+DB)**2+ADD))**3 &
             +(RR-DA-DB)/(SQRT((RR-DA-DB)**2+ADD))**3+(RR+DA+DB)/(SQRT((RR+DA+DB)**2+ADD))**3
      DZQXX = 2.0D0*(RR+DA)/(SQRT((RR+DA)**2+4.0D0*QB**2+ADQ))**3-2.0D0*(RR-DA)/(SQRT((RR-DA)**2+4.0D0*QB**2+ADQ))**3 &
             -2.0D0*(RR+DA)/(SQRT((RR+DA)**2+ADQ))**3+2.0D0*(RR-DA)/(SQRT((RR-DA)**2+ADQ))**3
      QXXDZ = 2.0D0*(RR-DB)/(SQRT((RR-DB)**2+4.0D0*QA**2+AQD))**3-2.0D0*(RR+DB)/(SQRT((RR+DB)**2+4.0D0*QA**2+AQD))**3 &
             -2.0D0*(RR-DB)/(SQRT((RR-DB)**2+AQD))**3+2.0D0*(RR+DB)/(SQRT((RR+DB)**2+AQD))**3
      DZQZZ = (RR+DA-2.0D0*QB)/(SQRT((RR+DA-2.0D0*QB)**2+ADQ))**3-(RR-DA-2.0D0*QB)/(SQRT((RR-DA-2.0D0*QB)**2+ADQ))**3 &
             +(RR+DA+2.0D0*QB)/(SQRT((RR+DA+2.0D0*QB)**2+ADQ))**3-(RR-DA+2.0D0*QB)/(SQRT((RR-DA+2.0D0*QB)**2+ADQ))**3 &
             +2.0D0*(RR-DA)/(SQRT((RR-DA)**2+ADQ))**3-2.0D0*(RR+DA)/(SQRT((RR+DA)**2+ADQ))**3
      QZZDZ = (RR+2.0D0*QA-DB)/(SQRT((RR+2.0D0*QA-DB)**2+AQD))**3-(RR+2.0D0*QA+DB)/(SQRT((RR+2.0D0*QA+DB)**2+AQD))**3 &
             +(RR-2.0D0*QA-DB)/(SQRT((RR-2.0D0*QA-DB)**2+AQD))**3-(RR-2.0D0*QA+DB)/(SQRT((RR-2.0D0*QA+DB)**2+AQD))**3 &
             -2.0D0*(RR-DB)/(SQRT((RR-DB)**2+AQD))**3+2.0D0*(RR+DB)/(SQRT((RR+DB)**2+AQD))**3
      QXXQXX=-(2.0D0*RR)/(SQRT(RR**2+4.0D0*(QA-QB)**2+AQQ))**3-(2.0D0*RR)/(SQRT(RR**2+4.0D0*(QA+QB)**2+AQQ))**3 &
             +(4.0D0*RR)/(SQRT(RR**2+4.0D0*QA**2+AQQ))**3+(4.0D0*RR)/(SQRT(RR**2+4.0D0*QB**2+AQQ))**3 &
             -(4.0D0*RR)/(SQRT(RR**2+AQQ))**3
      QXXQYY=-(4.0D0*RR)/(SQRT(RR**2+4.0D0*QA**2+4.0D0*QB**2+AQQ))**3+(4.0D0*RR)/(SQRT(RR**2+4.0D0*QA**2+AQQ))**3 &
             +(4.0D0*RR)/(SQRT(RR**2+4.0D0*QB**2+AQQ))**3-(4.0D0*RR)/(SQRT(RR**2+AQQ))**3
      QXXQZZ= &
           -2.0D0*(RR-2.0D0*QB)/(SQRT((RR-2.0D0*QB)**2+4.0D0*QA**2+AQQ))**3-&
            2.0D0*(RR+2.0D0*QB)/(SQRT((RR+2.0D0*QB)**2+4.0D0*QA**2+AQQ))**3 &
             +2.0D0*(RR-2.0D0*QB)/(SQRT((RR-2.0D0*QB)**2+AQQ))**3+2.0D0*(RR+2.0D0*QB)/(SQRT((RR+2.0D0*QB)**2+AQQ))**3 &
             +(4.0D0*RR)/(SQRT(RR**2+4.0D0*QA**2+AQQ))**3-(4.0D0*RR)/(SQRT(RR**2+AQQ))**3
      QZZQXX= &
           -2.0D0*(RR+2.0D0*QA)/(SQRT((RR+2.0D0*QA)**2+4.0D0*QB**2+AQQ))**3-&
            2.0D0*(RR-2.0D0*QA)/(SQRT((RR-2.0D0*QA)**2+4.0D0*QB**2+AQQ))**3 &
             +2.0D0*(RR+2.0D0*QA)/(SQRT((RR+2.0D0*QA)**2+AQQ))**3+2.0D0*(RR-2.0D0*QA)/(SQRT((RR-2.0D0*QA)**2+AQQ))**3 &
             +(4.0D0*RR)/(SQRT(RR**2+4.0D0*QB**2+AQQ))**3-(4.0D0*RR)/(SQRT(RR**2+AQQ))**3
      QZZQZZ= &
           -(RR+2.0D0*QA-2.0D0*QB)/(SQRT((RR+2.0D0*QA-2.0D0*QB)**2+AQQ))**3-&
            (RR+2.0D0*QA+2.0D0*QB)/(SQRT((RR+2.0D0*QA+2.0D0*QB)**2+AQQ))**3 &
           -(RR-2.0D0*QA-2.0D0*QB)/(SQRT((RR-2.0D0*QA-2.0D0*QB)**2+AQQ))**3-&
            (RR-2.0D0*QA+2.0D0*QB)/(SQRT((RR-2.0D0*QA+2.0D0*QB)**2+AQQ))**3 &
             +2.0D0*(RR-2.0D0*QA)/(SQRT((RR-2.0D0*QA)**2+AQQ))**3+2.0D0*(RR+2.0D0*QA)/(SQRT((RR+2.0D0*QA)**2+AQQ))**3 &
             +2.0D0*(RR-2.0D0*QB)/(SQRT((RR-2.0D0*QB)**2+AQQ))**3+2.0D0*(RR+2.0D0*QB)/(SQRT((RR+2.0D0*QB)**2+AQQ))**3 &
             -(4.0D0*RR)/(SQRT(RR**2+AQQ))**3
      DXQXZ = 2.0D0*(RR-QB)/(SQRT((RR-QB)**2+(DA-QB)**2+ADQ))**3-2.0D0*(RR+QB)/(SQRT((RR+QB)**2+(DA-QB)**2+ADQ))**3 &
             -2.0D0*(RR-QB)/(SQRT((RR-QB)**2+(DA+QB)**2+ADQ))**3+2.0D0*(RR+QB)/(SQRT((RR+QB)**2+(DA+QB)**2+ADQ))**3
      QXZDX = 2.0D0*(RR+QA)/(SQRT((RR+QA)**2+(QA-DB)**2+AQD))**3-2.0D0*(RR-QA)/(SQRT((RR-QA)**2+(QA-DB)**2+AQD))**3 &
             -2.0D0*(RR+QA)/(SQRT((RR+QA)**2+(QA+DB)**2+AQD))**3+2.0D0*(RR-QA)/(SQRT((RR-QA)**2+(QA+DB)**2+AQD))**3
      QXYQXY=-(4.0D0*RR)/(SQRT(RR**2+2.0D0*(QA-QB)**2+AQQ))**3-(4.0D0*RR)/(SQRT(RR**2+2.0D0*(QA+QB)**2+AQQ))**3 &
             +(8.0D0*RR)/(SQRT(RR**2+2.0D0*(QA**2+QB**2)+AQQ))**3
      QXZQXZ=-2.0D0*(RR+QA-QB)/(SQRT((RR+QA-QB)**2+(QA-QB)**2+AQQ))**3+2.0D0*(RR+QA+QB)/(SQRT((RR+QA+QB)**2+(QA-QB)**2+AQQ))**3 &
             +2.0D0*(RR-QA-QB)/(SQRT((RR-QA-QB)**2+(QA-QB)**2+AQQ))**3-2.0D0*(RR-QA+QB)/(SQRT((RR-QA+QB)**2+(QA-QB)**2+AQQ))**3 &
             +2.0D0*(RR+QA-QB)/(SQRT((RR+QA-QB)**2+(QA+QB)**2+AQQ))**3-2.0D0*(RR+QA+QB)/(SQRT((RR+QA+QB)**2+(QA+QB)**2+AQQ))**3 &
             -2.0D0*(RR-QA-QB)/(SQRT((RR-QA-QB)**2+(QA+QB)**2+AQQ))**3+2.0D0*(RR-QA+QB)/(SQRT((RR-QA+QB)**2+(QA+QB)**2+AQQ))**3
      DG(6)=(TERM*DZDZ)/4.0D0
      DG(7)=(TERM*DXDX)/4.0D0
      DG(8)=-TERM*(EDZ/2.0D0+QZZDZ/8.0D0)
      DG(9)=-TERM*(EDZ/2.0D0+QXXDZ/8.0D0)
      DG(10)=-(TERM*QXZDX)/8.0D0
      DG(13)=-TERM*(DZE/2.0D0+DZQZZ/8.0D0)
      DG(14)=-TERM*(DZE/2.0D0+DZQXX/8.0D0)
      DG(15)=-(TERM*DXQXZ)/8.0D0
      DG(16)=TERM*(EE+EQZZ/4.0D0+QZZE/4.0D0+QZZQZZ/16.0D0)
      DG(17)=TERM*(EE+EQZZ/4.0D0+QXXE/4.0D0+QXXQZZ/16.0D0)
      DG(18)=TERM*(EE+EQXX/4.0D0+QZZE/4.0D0+QZZQXX/16.0D0)
      DG(19)=TERM*(EE+EQXX/4.0D0+QXXE/4.0D0+QXXQXX/16.0D0)
      DG(20)=(TERM*QXZQXZ)/16.0D0
      DG(21)=TERM*(EE+EQXX/4.0D0+QXXE/4.0D0+QXXQYY/16.0D0)
      DG(22)=TERM*(QXXQXX-QXXQYY)/32.0D0
      END SUBROUTINE DELRI

      SUBROUTINE DER_ROTAT ( COORD, I, J, IX, RIJ, DEL1, IDX, TDX, TDY, TDZ, TX, TY, TZ )
      INTEGER,            INTENT(IN) :: I, IDX, IX, J
      REAL*8, INTENT(IN) :: RIJ, DEL1
      REAL*8, DIMENSION(1:3,1:2), INTENT(IN)  :: COORD
      REAL*8, DIMENSION(1:3),     INTENT(OUT) :: TDX, TDY, TDZ, TX, TY, TZ
      INTEGER            :: IJK
      REAL*8 :: RXY, RYZ, RZX, TERM, XD, YD, ZD

      XD=COORD(1,I)-COORD(1,J)
      YD=COORD(2,I)-COORD(2,J)
      ZD=COORD(3,I)-COORD(3,J)
      RXY=SQRT(XD*XD+YD*YD)
      RYZ=SQRT(YD*YD+ZD*ZD)
      RZX=SQRT(ZD*ZD+XD*XD)
      DO IJK=1,3
         TX(IJK)=0.0D0
         TY(IJK)=0.0D0
         TZ(IJK)=0.0D0
         TDX(IJK)=0.0D0
         TDY(IJK)=0.0D0
         TDZ(IJK)=0.0D0
      END DO
      IF(RXY.LT.1.0D-4) THEN
         !   MOLECULAR Z AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(3)=1.0D0
         IF(ZD.LT.0.0D0) TX(3)=-1.0D0
         TY(2)=1.0D0
         TZ(1)=TX(3)
         IF(IDX.EQ.1) RETURN
         IF(IX.EQ.1) TDX(1)=1.0D0/RIJ
         IF(IX.EQ.2) TDX(2)=1.0D0/RIJ
         IF(IX.EQ.1) TDZ(3)=-1.0D0/RIJ
         IF(IX.EQ.2) TDY(3)=-TX(3)/RIJ
      ELSEIF(RYZ.LT.1.0D-4) THEN
         !   MOLECULAR X AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(1)=1.0D0
         IF(XD.LT.0.0D0) TX(1)=-1.0D0
         TY(2)=TX(1)
         TZ(3)=1.0D0
         IF(IDX.EQ.1) RETURN
         IF(IX.EQ.2) TDX(2)=1.0D0/RIJ
         IF(IX.EQ.3) TDX(3)=1.0D0/RIJ
         IF(IX.EQ.2) TDY(1)=-1.0D0/RIJ
         IF(IX.EQ.3) TDZ(1)=-TX(1)/RIJ
      ELSEIF(RZX.LT.1.0D-4) THEN
         !   MOLECULAR Y AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(2)=1.0D0
         IF(YD.LT.0.0D0) TX(2)=-1.0D0
         TY(1)=-TX(2)
         TZ(3)=1.0D0
         IF(IDX.EQ.1) RETURN
         IF(IX.EQ.1) TDX(1)=1.0D0/RIJ
         IF(IX.EQ.3) TDX(3)=1.0D0/RIJ
         IF(IX.EQ.1) TDY(2)=1.0D0/RIJ
         IF(IX.EQ.3) TDZ(2)=-TX(2)/RIJ
      ELSE
         TX(1)=XD/RIJ
         TX(2)=YD/RIJ
         TX(3)=ZD/RIJ
         TZ(3)=RXY/RIJ
         TY(1)=-TX(2)*SIGN(+1.0D0,TX(1))/TZ(3)
         TY(2)=ABS(TX(1)/TZ(3))
         TY(3)=0.0D0
         TZ(1)=-TX(1)*TX(3)/TZ(3)
         TZ(2)=-TX(2)*TX(3)/TZ(3)
         IF(IDX.EQ.1) RETURN
         TERM=DEL1/(RIJ*RIJ)
         IF(IX.EQ.1)THEN
            TDX(1)=1.0D0/RIJ-TX(1)*TERM
            TDX(2)=-TX(2)*TERM
            TDX(3)=-TX(3)*TERM
            TDZ(3)=TX(1)/RXY-TZ(3)*TERM
         ELSEIF(IX.EQ.2) THEN
            TDX(1)=-TX(1)*TERM
            TDX(2)=1.0D0/RIJ-TX(2)*TERM
            TDX(3)=-TX(3)*TERM
            TDZ(3)=TX(2)/RXY-TZ(3)*TERM
         ELSEIF(IX.EQ.3)THEN
            TDX(1)=-TX(1)*TERM
            TDX(2)=-TX(2)*TERM
            TDX(3)=1.0D0/RIJ-TX(3)*TERM
            TDZ(3)=-TZ(3)*TERM
         ENDIF
         TDY(1)=-TDX(2)/TZ(3)+TX(2)*TDZ(3)/TZ(3)**2
         IF(TX(1).LT.0.0D0) TDY(1)=-TDY(1)
         TDY(2)=TDX(1)/TZ(3)-TX(1)*TDZ(3)/TZ(3)**2
         IF(TX(1).LT.0.0D0) TDY(2)=-TDY(2)
         TDY(3)=0.0D0
         TDZ(1)=-TX(3)*TDX(1)/TZ(3)-TX(1)*TDX(3)/TZ(3)+TX(1)*TX(3)*TDZ(3)/TZ(3)**2
         TDZ(2)=-TX(3)*TDX(2)/TZ(3)-TX(2)*TDX(3)/TZ(3)+TX(2)*TX(3)*TDZ(3)/TZ(3)**2
      END IF
      END SUBROUTINE DER_ROTAT

   END SUBROUTINE ANALYT


   SUBROUTINE MOPAC_CHARGES ( CHARGES )
   REAL*8, DIMENSION(1:NAQM), INTENT(OUT) :: CHARGES
   INTEGER :: I, IFIRST, ILAST, IQM
   REAL*8, DIMENSION(1:NBASIS) :: BQTEMP
   CHARGES = 0.0D0
   DO I = 1,NBASIS
      BQTEMP(I) = DENMAT((I*(I+1))/2)
   END DO
   DO I = 1,NAQM
       IFIRST = BFIRST(I)
       ILAST  = BLAST(I)
       IQM    = ANINT( ATMCHG(I) )
       CHARGES(I) = CHARGES(I) + ( CORE(IQM) - SUM ( BQTEMP(IFIRST:ILAST) ) )
   END DO
   END SUBROUTINE MOPAC_CHARGES


! ========================================================================================================


   SUBROUTINE MOPAC_PARAMETERS_INITIALIZE ( METHOD )
   CHARACTER ( LEN = * ), INTENT(IN) :: METHOD
   CALL PARAMETERS_COMMON

   IF ( METHOD == "AM1" ) THEN
      CALL PARAMETERS_AM1
   ELSE IF ( METHOD == "MNDO" ) THEN
      CALL PARAMETERS_MNDO
   ELSE IF ( METHOD == "PDDG" ) THEN
      CALL PARAMETERS_PDDG
   ELSE IF ( METHOD == "PM3" ) THEN
      CALL PARAMETERS_PM3
   ELSE IF ( METHOD == "RM1" ) THEN
      CALL PARAMETERS_RM1
   END IF
   END SUBROUTINE MOPAC_PARAMETERS_INITIALIZE


   SUBROUTINE PARAMETERS_COMMON
   SEPAR = .FALSE.
   AD    = 0.0D0 ; ALP   = 0.0D0 ; AM    = 0.0D0 ; AQ   = 0.0D0
   BETAD = 0.0D0 ; BETAP = 0.0D0 ; BETAS = 0.0D0 ; CORE = 0.0D0
   DD    = 0.0D0 ; EHEAT = 0.0D0 ; EISOL = 0.0D0 ; GDD  = 0.0D0
   GPD   = 0.0D0 ; GSD   = 0.0D0 ; GP2   = 0.0D0 ; GPP  = 0.0D0
   GSP   = 0.0D0 ; GSS   = 0.0D0 ; HSP   = 0.0D0 ; QQ   = 0.0D0
   USS   = 0.0D0 ; UPP   = 0.0D0 ; UDD   = 0.0D0 ; ZD   = 0.0D0
   ZP    = 0.0D0 ; ZS    = 0.0D0
   FN1   = 0.0D0 ; FN2   = 0.0D0 ; FN3   = 0.0D0
   NATORB (  1) = 1 ; CORE (  1) =  1.0D0 ; EHEAT (  1) =  52.102D0
   NATORB (  2) = 1 ; CORE (  2) =  0.0D0 ; EHEAT (  2) =   0.000D0
   NATORB (  3) = 4 ; CORE (  3) =  1.0D0 ; EHEAT (  3) =  38.410D0
   NATORB (  4) = 4 ; CORE (  4) =  2.0D0 ; EHEAT (  4) =  76.960D0
   NATORB (  5) = 4 ; CORE (  5) =  3.0D0 ; EHEAT (  5) = 135.700D0
   NATORB (  6) = 4 ; CORE (  6) =  4.0D0 ; EHEAT (  6) = 170.890D0
   NATORB (  7) = 4 ; CORE (  7) =  5.0D0 ; EHEAT (  7) = 113.000D0
   NATORB (  8) = 4 ; CORE (  8) =  6.0D0 ; EHEAT (  8) =  59.559D0
   NATORB (  9) = 4 ; CORE (  9) =  7.0D0 ; EHEAT (  9) =  18.890D0
   NATORB ( 10) = 4 ; CORE ( 10) =  0.0D0 ; EHEAT ( 10) =   0.000D0
   NATORB ( 11) = 0 ; CORE ( 11) =  1.0D0 ; EHEAT ( 11) =  25.850D0
   NATORB ( 12) = 4 ; CORE ( 12) =  2.0D0 ; EHEAT ( 12) =  35.000D0
   NATORB ( 13) = 4 ; CORE ( 13) =  3.0D0 ; EHEAT ( 13) =  79.490D0
   NATORB ( 14) = 4 ; CORE ( 14) =  4.0D0 ; EHEAT ( 14) = 108.390D0
   NATORB ( 15) = 4 ; CORE ( 15) =  5.0D0 ; EHEAT ( 15) =  75.570D0
   NATORB ( 16) = 4 ; CORE ( 16) =  6.0D0 ; EHEAT ( 16) =  66.400D0
   NATORB ( 17) = 4 ; CORE ( 17) =  7.0D0 ; EHEAT ( 17) =  28.990D0
   NATORB ( 18) = 4 ; CORE ( 18) =  0.0D0 ; EHEAT ( 18) =   0.000D0
   NATORB ( 19) = 0 ; CORE ( 19) =  1.0D0 ; EHEAT ( 19) =  21.420D0
   NATORB ( 20) = 4 ; CORE ( 20) =  2.0D0 ; EHEAT ( 20) =  42.600D0
   NATORB ( 21) = 9 ; CORE ( 21) =  3.0D0 ; EHEAT ( 21) =  90.300D0
   NATORB ( 22) = 9 ; CORE ( 22) =  4.0D0 ; EHEAT ( 22) = 112.300D0
   NATORB ( 23) = 9 ; CORE ( 23) =  5.0D0 ; EHEAT ( 23) = 122.900D0
   NATORB ( 24) = 9 ; CORE ( 24) =  6.0D0 ; EHEAT ( 24) =  95.000D0
   NATORB ( 25) = 9 ; CORE ( 25) =  7.0D0 ; EHEAT ( 25) =  67.700D0
   NATORB ( 26) = 9 ; CORE ( 26) =  8.0D0 ; EHEAT ( 26) =  99.300D0
   NATORB ( 27) = 9 ; CORE ( 27) =  9.0D0 ; EHEAT ( 27) = 102.400D0
   NATORB ( 28) = 9 ; CORE ( 28) = 10.0D0 ; EHEAT ( 28) = 102.800D0
   NATORB ( 29) = 9 ; CORE ( 29) = 11.0D0 ; EHEAT ( 29) =  80.700D0
   NATORB ( 30) = 4 ; CORE ( 30) =  2.0D0 ; EHEAT ( 30) =  31.170D0
   NATORB ( 31) = 4 ; CORE ( 31) =  3.0D0 ; EHEAT ( 31) =  65.400D0
   NATORB ( 32) = 4 ; CORE ( 32) =  4.0D0 ; EHEAT ( 32) =  89.500D0
   NATORB ( 33) = 4 ; CORE ( 33) =  5.0D0 ; EHEAT ( 33) =  72.300D0
   NATORB ( 34) = 4 ; CORE ( 34) =  6.0D0 ; EHEAT ( 34) =  54.300D0
   NATORB ( 35) = 4 ; CORE ( 35) =  7.0D0 ; EHEAT ( 35) =  26.740D0
   NATORB ( 36) = 4 ; CORE ( 36) =  0.0D0 ; EHEAT ( 36) =   0.000D0
   NATORB ( 37) = 4 ; CORE ( 37) =  1.0D0 ; EHEAT ( 37) =  19.600D0
   NATORB ( 38) = 4 ; CORE ( 38) =  2.0D0 ; EHEAT ( 38) =  39.100D0
   NATORB ( 39) = 9 ; CORE ( 39) =  3.0D0 ; EHEAT ( 39) = 101.500D0
   NATORB ( 40) = 9 ; CORE ( 40) =  4.0D0 ; EHEAT ( 40) = 145.500D0
   NATORB ( 41) = 9 ; CORE ( 41) =  5.0D0 ; EHEAT ( 41) = 172.400D0
   NATORB ( 42) = 9 ; CORE ( 42) =  6.0D0 ; EHEAT ( 42) = 157.300D0 
   NATORB ( 43) = 9 ; CORE ( 43) =  7.0D0 ; EHEAT ( 43) =   0.000D0
   NATORB ( 44) = 9 ; CORE ( 44) =  8.0D0 ; EHEAT ( 44) = 155.500D0
   NATORB ( 45) = 9 ; CORE ( 45) =  9.0D0 ; EHEAT ( 45) = 133.000D0
   NATORB ( 46) = 9 ; CORE ( 46) = 10.0D0 ; EHEAT ( 46) =  90.000D0
   NATORB ( 47) = 9 ; CORE ( 47) = 11.0D0 ; EHEAT ( 47) =  68.100D0
   NATORB ( 48) = 4 ; CORE ( 48) =  2.0D0 ; EHEAT ( 48) =  26.720D0
   NATORB ( 49) = 4 ; CORE ( 49) =  3.0D0 ; EHEAT ( 49) =  58.000D0
   NATORB ( 50) = 4 ; CORE ( 50) =  4.0D0 ; EHEAT ( 50) =  72.200D0
   NATORB ( 51) = 4 ; CORE ( 51) =  5.0D0 ; EHEAT ( 51) =  63.200D0
   NATORB ( 52) = 4 ; CORE ( 52) =  6.0D0 ; EHEAT ( 52) =  47.000D0 
   NATORB ( 53) = 4 ; CORE ( 53) =  7.0D0 ; EHEAT ( 53) =  25.517D0
   NATORB ( 54) = 4 ; CORE ( 54) =  0.0D0 ; EHEAT ( 54) =   0.000D0
   NATORB ( 55) = 2 ; CORE ( 55) =  1.0D0 ; EHEAT ( 55) =  18.700D0
   NATORB ( 56) = 2 ; CORE ( 56) =  2.0D0 ; EHEAT ( 56) =  42.500D0
   NATORB ( 57) = 8 ; CORE ( 57) =  3.0D0 ; EHEAT ( 57) =   0.000D0
   NATORB ( 58) = 8 ; CORE ( 58) =  4.0D0 ; EHEAT ( 58) = 101.300D0
   NATORB ( 59) = 8 ; CORE ( 59) =  5.0D0 ; EHEAT ( 59) =   0.000D0
   NATORB ( 60) = 8 ; CORE ( 60) =  6.0D0 ; EHEAT ( 60) =   0.000D0
   NATORB ( 61) = 8 ; CORE ( 61) =  7.0D0 ; EHEAT ( 61) =   0.000D0
   NATORB ( 62) = 8 ; CORE ( 62) =  8.0D0 ; EHEAT ( 62) =  49.400D0
   NATORB ( 63) = 8 ; CORE ( 63) =  9.0D0 ; EHEAT ( 63) =   0.000D0
   NATORB ( 64) = 8 ; CORE ( 64) = 10.0D0 ; EHEAT ( 64) =   0.000D0
   NATORB ( 65) = 8 ; CORE ( 65) = 11.0D0 ; EHEAT ( 65) =   0.000D0
   NATORB ( 66) = 8 ; CORE ( 66) = 12.0D0 ; EHEAT ( 66) =   0.000D0
   NATORB ( 67) = 8 ; CORE ( 67) = 13.0D0 ; EHEAT ( 67) =   0.000D0
   NATORB ( 68) = 8 ; CORE ( 68) = 14.0D0 ; EHEAT ( 68) =  75.800D0
   NATORB ( 69) = 8 ; CORE ( 69) = 15.0D0 ; EHEAT ( 69) =   0.000D0
   NATORB ( 70) = 8 ; CORE ( 70) = 16.0D0 ; EHEAT ( 70) =  36.350D0
   NATORB ( 71) = 9 ; CORE ( 71) =  3.0D0 ; EHEAT ( 71) =   0.000D0
   NATORB ( 72) = 9 ; CORE ( 72) =  4.0D0 ; EHEAT ( 72) = 148.000D0
   NATORB ( 73) = 9 ; CORE ( 73) =  5.0D0 ; EHEAT ( 73) = 186.900D0
   NATORB ( 74) = 9 ; CORE ( 74) =  6.0D0 ; EHEAT ( 74) = 203.100D0
   NATORB ( 75) = 9 ; CORE ( 75) =  7.0D0 ; EHEAT ( 75) = 185.000D0
   NATORB ( 76) = 9 ; CORE ( 76) =  8.0D0 ; EHEAT ( 76) = 188.000D0
   NATORB ( 77) = 9 ; CORE ( 77) =  9.0D0 ; EHEAT ( 77) = 160.000D0
   NATORB ( 78) = 9 ; CORE ( 78) = 10.0D0 ; EHEAT ( 78) = 135.200D0
   NATORB ( 79) = 9 ; CORE ( 79) = 11.0D0 ; EHEAT ( 79) =  88.000D0
   NATORB ( 80) = 4 ; CORE ( 80) =  2.0D0 ; EHEAT ( 80) =  14.690D0
   NATORB ( 81) = 4 ; CORE ( 81) =  3.0D0 ; EHEAT ( 81) =  43.550D0
   NATORB ( 82) = 4 ; CORE ( 82) =  4.0D0 ; EHEAT ( 82) =  46.620D0
   NATORB ( 83) = 4 ; CORE ( 83) =  5.0D0 ; EHEAT ( 83) =  50.100D0
   NATORB ( 84) = 4 ; CORE ( 84) =  6.0D0 ; EHEAT ( 84) =   0.000D0
   NATORB ( 85) = 4 ; CORE ( 85) =  7.0D0 ; EHEAT ( 85) =   0.000D0
   END SUBROUTINE PARAMETERS_COMMON


   SUBROUTINE PARAMETERS_AM1
   ! . Hydrogen
   SEPAR(1) = .TRUE.
   ALP(1)   =   2.8823240D0
   EISOL(1) = -11.3964270D0
   BETAS(1) =  -6.1737870D0
   ZS(1)    =   1.1880780D0
   AM(1)    =   0.4721793D0
   AD(1)    =   0.4721793D0
   AQ(1)    =   0.4721793D0
   USS(1)   = -11.3964270D0
   GSS(1)   =  12.8480000D0

   FN1(1,1) =  0.1227960D0 ; FN2(1,1) = 5.0000000D0 ; FN3(1,1) = 1.2000000D0
   FN1(1,2) =  0.0050900D0 ; FN2(1,2) = 5.0000000D0 ; FN3(1,2) = 1.8000000D0
   FN1(1,3) = -0.0183360D0 ; FN2(1,3) = 2.0000000D0 ; FN3(1,3) = 2.1000000D0

   ! . Lithium
   SEPAR(3) = .TRUE.
   ALP(3)   =   1.2501400D0
   EISOL(3) =  -5.1280000D0
   BETAS(3) =  -1.3500400D0
   BETAP(3) =  -1.3500400D0
   ZS(3)    =   0.7023800D0
   ZP(3)    =   0.7023800D0
   DD(3)    =   2.0549783D0
   QQ(3)    =   1.7437069D0
   AM(3)    =   0.2682837D0
   AD(3)    =   0.2269793D0
   AQ(3)    =   0.2614581D0
   USS(3)   =  -5.1280000D0
   UPP(3)   =  -2.7212000D0
   GSS(3)   =   7.3000000D0
   GSP(3)   =   5.4200000D0
   GPP(3)   =   5.0000000D0
   GP2(3)   =   4.5200000D0
   HSP(3)   =   0.8300000D0

   ! . Beryllium
   SEPAR(4) = .TRUE.
   ALP(4)   =   1.6694340D0
   EISOL(4) = -24.2047560D0
   BETAS(4) =  -4.0170960D0
   BETAP(4) =  -4.0170960D0
   ZS(4)    =   1.0042100D0
   ZP(4)    =   1.0042100D0
   DD(4)    =   1.4373245D0
   QQ(4)    =   1.2196103D0
   AM(4)    =   0.3307607D0
   AD(4)    =   0.3356142D0
   AQ(4)    =   0.3846373D0
   USS(4)   = -16.6023780D0
   UPP(4)   = -10.7037710D0
   GSS(4)   =   9.0000000D0
   GSP(4)   =   7.4300000D0
   GPP(4)   =   6.9700000D0
   GP2(4)   =   6.2200000D0
   HSP(4)   =   1.2800000D0

   ! . Boron
   SEPAR(5) = .TRUE.
   ALP(5)   =   2.4469090D0
   EISOL(5) = -63.7172650D0
   BETAS(5) =  -9.5991140D0
   BETAP(5) =  -6.2737570D0
   ZS(5)    =   1.6117090D0
   ZP(5)    =   1.5553850D0
   DD(5)    =   0.9107622D0
   QQ(5)    =   0.7874223D0
   AM(5)    =   0.3891951D0
   AD(5)    =   0.5045152D0
   AQ(5)    =   0.5678856D0
   USS(5)   = -34.4928700D0
   UPP(5)   = -22.6315250D0
   GSS(5)   =  10.5900000D0
   GSP(5)   =   9.5600000D0
   GPP(5)   =   8.8600000D0
   GP2(5)   =   7.8600000D0
   HSP(5)   =   1.8100000D0

   ! . Carbon
   SEPAR(6) = .TRUE.
   ALP(6)   =    2.6482740D0
   EISOL(6) = -120.8157940D0
   BETAS(6) =  -15.7157830D0
   BETAP(6) =   -7.7192830D0
   ZS(6)    =    1.8086650D0
   ZP(6)    =    1.6851160D0
   DD(6)    =    0.8236736D0
   QQ(6)    =    0.7268015D0
   AM(6)    =    0.4494671D0
   AD(6)    =    0.6082946D0
   AQ(6)    =    0.6423492D0
   USS(6)   =  -52.0286580D0
   UPP(6)   =  -39.6142390D0
   GSS(6)   =   12.2300000D0
   GSP(6)   =   11.4700000D0
   GPP(6)   =   11.0800000D0
   GP2(6)   =    9.8400000D0
   HSP(6)   =    2.4300000D0

   FN1(6,1) =  0.0113550D0 ; FN2(6,1) =  5.0000000D0 ; FN3(6,1) =  1.6000000D0
   FN1(6,2) =  0.0459240D0 ; FN2(6,2) =  5.0000000D0 ; FN3(6,2) =  1.8500000D0
   FN1(6,3) = -0.0200610D0 ; FN2(6,3) =  5.0000000D0 ; FN3(6,3) =  2.0500000D0
   FN1(6,4) = -0.0012600D0 ; FN2(6,4) =  5.0000000D0 ; FN3(6,4) =  2.6500000D0

   ! . Nitrogen
   SEPAR(7) = .TRUE.
   ALP(7)   =    2.9472860D0
   EISOL(7) = -202.4077430D0
   BETAS(7) =  -20.2991100D0
   BETAP(7) =  -18.2386660D0
   ZS(7)    =    2.3154100D0
   ZP(7)    =    2.1579400D0
   DD(7)    =    0.6433247D0
   QQ(7)    =    0.5675528D0
   AM(7)    =    0.4994487D0
   AD(7)    =    0.7820840D0
   AQ(7)    =    0.7883498D0
   USS(7)   =  -71.8600000D0
   UPP(7)   =  -57.1675810D0
   GSS(7)   =   13.5900000D0
   GSP(7)   =   12.6600000D0
   GPP(7)   =   12.9800000D0
   GP2(7)   =   11.5900000D0
   HSP(7)   =    3.1400000D0

   FN1(7,1) =  0.0252510D0 ; FN2(7,1) =  5.0000000D0 ; FN3(7,1) =  1.5000000D0
   FN1(7,2) =  0.0289530D0 ; FN2(7,2) =  5.0000000D0 ; FN3(7,2) =  2.1000000D0
   FN1(7,3) = -0.0058060D0 ; FN2(7,3) =  2.0000000D0 ; FN3(7,3) =  2.4000000D0

   ! . Oxygen
   SEPAR(8) = .TRUE.
   ALP(8)   =    4.4553710D0
   EISOL(8) = -316.0995200D0
   BETAS(8) =  -29.2727730D0
   BETAP(8) =  -29.2727730D0
   ZS(8)    =    3.1080320D0
   ZP(8)    =    2.5240390D0
   DD(8)    =    0.4988896D0
   QQ(8)    =    0.4852322D0
   AM(8)    =    0.5667034D0
   AD(8)    =    0.9961066D0
   AQ(8)    =    0.9065223D0
   USS(8)   =  -97.8300000D0
   UPP(8)   =  -78.2623800D0
   GSS(8)   =   15.4200000D0
   GSP(8)   =   14.4800000D0
   GPP(8)   =   14.5200000D0
   GP2(8)   =   12.9800000D0
   HSP(8)   =    3.9400000D0

   FN1(8,1) = 0.2809620D0 ; FN2(8,1) = 5.0000000D0 ; FN3(8,1) = 0.8479180D0
   FN1(8,2) = 0.0814300D0 ; FN2(8,2) = 7.0000000D0 ; FN3(8,2) = 1.4450710D0

   ! . Fluorine
   SEPAR(9) = .TRUE.
   ALP(9)   =    5.5178000D0
   EISOL(9) = -482.2905830D0
   BETAS(9) =  -69.5902770D0
   BETAP(9) =  -27.9223600D0
   ZS(9)    =    3.7700820D0
   ZP(9)    =    2.4946700D0
   DD(9)    =    0.4145203D0
   QQ(9)    =    0.4909446D0
   AM(9)    =    0.6218302D0
   AD(9)    =    1.2088792D0
   AQ(9)    =    0.9449355D0
   USS(9)   = -136.1055790D0
   UPP(9)   = -104.8898850D0
   GSS(9)   =   16.9200000D0
   GSP(9)   =   17.2500000D0
   GPP(9)   =   16.7100000D0
   GP2(9)   =   14.9100000D0
   HSP(9)   =    4.8300000D0

   FN1(9,1) = 0.2420790D0 ; FN2(9,1) = 4.8000000D0 ; FN3(9,1) = 0.9300000D0
   FN1(9,2) = 0.0036070D0 ; FN2(9,2) = 4.6000000D0 ; FN3(9,2) = 1.6600000D0

   ! . Sodium
   SEPAR(11) = .TRUE.
   ALP(11)   =  1.32D0
   EISOL(11) =  0.00D0
   AM(11)    =  0.50D0

   ! . Magnesium (Hutter et al., JPCB 102, 8080-8090, 1998)
   SEPAR(12) = .TRUE.
   USS  ( 12) =  -14.96959313D0
   UPP  ( 12) =  -11.56229248D0
   BETAS( 12) =   -1.25974355D0
   BETAP( 12) =   -0.77836604D0
   ZS   ( 12) =    1.22339270D0
   ZP   ( 12) =    1.02030798D0
   ALP  ( 12) =    1.67049799D0
   EISOL( 12) =  -22.43786349D0 ! From CALPAR formula. Otherwise PM3 value?
   GSS  ( 12) =    7.50132277D0
   GSP  ( 12) =    6.34591536D0
   GPP  ( 12) =    4.77534467D0
   GP2  ( 12) =    4.34017279D0
   HSP  ( 12) =    0.48930466D0
   DD   ( 12) =    1.75012115D0 ! From formula (15) page 94 Thiel + Dewar
   QQ   ( 12) =    1.64001467D0 ! From formula (16) page 95 Thiel + Dewar
   AM   ( 12) =    0.27568257D0 ! From formula (21) page 95 Thiel + Dewar (rho0 = 1/2 Gss = 1/2 AM so AM = Gss/AU_TO_EV)
   AD   ( 12) =    0.19957236D0 ! By solving numerically equation (19) page 95 Thiel + Dewar
   AQ   ( 12) =    0.26440152D0 ! By solving numerically equation (20) page 95 Thiel + Dewar

   FN1( 12,1) =  2.55017735D0 ; FN2( 12,1) = 4.29397225D0 ; FN3( 12,1) = 0.79989601D0
   FN1( 12,2) = -0.00565806D0 ; FN2( 12,2) = 2.96053910D0 ; FN3( 12,2) = 1.47499983D0
   FN1( 12,3) = -0.00610286D0 ; FN2( 12,3) = 2.61416919D0 ; FN3( 12,3) = 2.42604040D0

   ! . Aluminium
   SEPAR(13) = .TRUE.
   ALP(13)   =   1.9765860D0
   EISOL(13) = -46.4208150D0
   BETAS(13) =  -3.8668220D0
   BETAP(13) =  -2.3171460D0
   ZS(13)    =   1.5165930D0
   ZP(13)    =   1.3063470D0
   ZD(13)    =   1.0000000D0
   DD(13)    =   1.4040443D0
   QQ(13)    =   1.2809154D0
   AM(13)    =   0.2973172D0
   AD(13)    =   0.2630229D0
   AQ(13)    =   0.3427832D0
   USS(13)   = -24.3535850D0
   UPP(13)   = -18.3636450D0
   GSS(13)   =   8.0900000D0
   GSP(13)   =   6.6300000D0
   GPP(13)   =   5.9800000D0
   GP2(13)   =   5.4000000D0
   HSP(13)   =   0.7000000D0

   FN1(13,1) = 0.0900000D0 ; FN2(13,1) = 12.3924430D0 ; FN3(13,1) = 2.0503940D0

   ! . Silicon
   SEPAR(14) = .TRUE.
   ALP(14)   =   2.257816D0
   EISOL(14) = -79.0017420D0
   BETAS(14) =  -3.784852D0
   BETAP(14) =  -1.968123D0
   ZS(14)    =   1.830697D0
   ZP(14)    =   1.2849530D0
   ZD(14)    =   1.0000000D0
   DD(14)    =   1.1631107D0
   QQ(14)    =   1.3022422D0
   AM(14)    =   0.3608967D0
   AD(14)    =   0.3829813D0
   AQ(14)    =   0.3999442D0
   USS(14)   = -33.9536220D0
   UPP(14)   = -28.9347490D0
   GSS(14)   =   9.8200000D0
   GSP(14)   =   8.3600000D0
   GPP(14)   =   7.3100000D0
   GP2(14)   =   6.5400000D0
   HSP(14)   =   1.3200000D0

   ! . Phosphorus
   SEPAR(15) = .TRUE.
   ALP(15)   =    2.4553220D0
   EISOL(15) = -124.4368355D0
   BETAS(15) =   -6.3537640D0
   BETAP(15) =   -6.5907090D0
   ZS(15)    =    1.9812800D0
   ZP(15)    =    1.8751500D0
   ZD(15)    =    1.0000000D0
   DD(15)    =    1.0452022D0
   QQ(15)    =    0.8923660D0
   AM(15)    =    0.4248440D0
   AD(15)    =    0.3275319D0
   AQ(15)    =    0.4386854D0
   USS(15)   =  -42.0298630D0
   UPP(15)   =  -34.0307090D0
   GSS(15)   =   11.5600050D0
   GSP(15)   =   5.2374490D0
   GPP(15)   =    7.8775890D0
   GP2(15)   =    7.3076480D0
   HSP(15)   =    0.7792380D0

   ! . Sulphur
   SEPAR(16) = .TRUE.
   ALP(16)   =    2.4616480D0
   EISOL(16) = -191.7321930D0
   BETAS(16) =  -3.9205660D0
   BETAP(16) =  -7.9052780D0
   ZS(16)    =    2.3665150D0
   ZP(16)    =    1.6672630D0
   ZD(16)    =    1.0000000D0
   DD(16)    =    0.9004265D0
   QQ(16)    =    1.0036329D0
   AM(16)    =    0.4331617D0
   AD(16)    =    0.5907115D0
   AQ(16)    =    0.6454943D0
   USS(16)   =  -56.6940560D0
   UPP(16)   =  -48.7170490D0
   GSS(16)   =   11.7863290D0
   GSP(16)   =    8.6631270D0
   GPP(16)   =   10.0393080D0
   GP2(16)   =    7.7816880D0
   HSP(16)   =    2.5321370D0

   ! . Chlorine
   SEPAR(17) = .TRUE.
   ALP(17)   =    2.9193680D0
   EISOL(17) = -372.1984310D0
   BETAS(17) =  -24.5946700D0
   BETAP(17) =  -14.6372160D0
   ZS(17)    =    3.6313760D0
   ZP(17)    =    2.0767990D0
   ZD(17)    =    1.0000000D0
   DD(17)    =    0.5406286D0
   QQ(17)    =    0.8057208D0
   AM(17)    =    0.5523705D0
   AD(17)    =    0.7693200D0
   AQ(17)    =    0.6133369D0
   USS(17)   = -111.6139480D0
   UPP(17)   =  -76.6401070D0
   GSS(17)   =   15.0300000D0
   GSP(17)   =   13.1600000D0
   GPP(17)   =   11.3000000D0
   GP2(17)   =    9.9700000D0
   HSP(17)   =    2.4200000D0

   FN1(17,1) = 0.0942430D0 ; FN2(17,1) = 4.0000000D0 ; FN3(17,1) = 1.3000000D0
   FN1(17,2) = 0.0271680D0 ; FN2(17,2) = 4.0000000D0 ; FN3(17,2) = 2.1000000D0

   ! . Potassium
   SEPAR(19) = .TRUE.
   ALP(19)   =  1.16D0
   EISOL(19) =  0.00D0
   AM(19)    =  0.50D0

   ! . Zinc
   SEPAR(30)  = .TRUE.
   USS(30)    =    -21.0400080D0
   UPP(30)    =    -17.6555740D0
   BETAS(30)  =     -1.9974290D0
   BETAP(30)  =     -4.7581190D0
   ZS(30)     =      1.9542990D0
   ZP(30)     =      1.3723650D0
   ZD(30)     =      1.0000000D0
   ALP(30)    =      1.4845630D0
   EISOL(30)  =    -30.2800160D0
   GSS(30)    =     11.8000000D0
   GSP(30)    =     11.1820180D0
   GPP(30)    =     13.3000000D0
   GP2(30)    =     12.9305200D0
   HSP(30)    =      0.4846060D0
   DD(30)     =      1.3581113D0
   QQ(30)     =      1.5457406D0
   AM(30)     =      0.4336641D0
   AD(30)     =      0.2317423D0
   AQ(30)     =      0.2621165D0

   ! . Bromine
   SEPAR(35) = .TRUE.
   ALP(35)   =    2.5765460D0
   EISOL(35) = -352.3142087D0
   BETAS(35) =  -19.3998800D0
   BETAP(35) =   -8.9571950D0
   ZS(35)    =    3.0641330D0
   ZP(35)    =    2.0383330D0
   ZD(35)    =    1.0000000D0
   DD(35)    =    0.8458104D0
   QQ(35)    =    1.0407133D0
   AM(35)    =    0.5526071D0
   AD(35)    =    0.6024598D0
   AQ(35)    =    0.5307555D0
   USS(35)   = -104.6560630D0
   UPP(35)   =  -74.9300520D0
   GSS(35)   =   15.0364395D0
   GSP(35)   =   13.0346824D0
   GPP(35)   =   11.2763254D0
   GP2(35)   =    9.8544255D0
   HSP(35)   =    2.4558683D0

   FN1(35,1) = 0.0666850D0 ; FN2(35,1) = 4.0000000D0 ; FN3(35,1) = 1.5000000D0
   FN1(35,2) = 0.0255680D0 ; FN2(35,2) = 4.0000000D0 ; FN3(35,2) = 2.3000000D0

   ! . Iodine
   SEPAR(53) = .TRUE.
   USS(53)   = -103.5896630D0
   UPP(53)   =  -74.4299970D0
   BETAS(53) =   -8.4433270D0
   BETAP(53) =   -6.3234050D0
   ZS(53)    =    2.1028580D0
   ZP(53)    =    2.1611530D0
   ZD(53)    =    1.0000000D0
   ALP(53)   =    2.2994240D0
   EISOL(53) = -346.8642857D0
   GSS(53)   =   15.0404486D0
   GSP(53)   =   13.0565580D0
   GPP(53)   =   11.1477837D0
   GP2(53)   =    9.9140907D0
   HSP(53)   =    2.4563820D0
   DD(53)    =    1.4878778D0
   QQ(53)    =    1.1887388D0
   AM(53)    =    0.5527544D0
   AD(53)    =    0.4497523D0
   AQ(53)    =    0.4631775D0

   FN1(53,1) = 0.0043610D0 ; FN2(53,1) = 2.3000000D0 ; FN3(53,1) = 1.8000000D0
   FN1(53,2) = 0.0157060D0 ; FN2(53,2) = 3.0000000D0 ; FN3(53,2) = 2.2400000D0
   END SUBROUTINE PARAMETERS_AM1


   SUBROUTINE PARAMETERS_MNDO
   ! . Hydrogen
   SEPAR(1) = .TRUE.
   ALP(1)   =   2.5441341D0
   EISOL(1) = -11.9062760D0
   BETAS(1) =  -6.9890640D0
   ZS(1)    =   1.3319670D0
   AM(1)    =   0.4721793D0
   AD(1)    =   0.4721793D0
   AQ(1)    =   0.4721793D0
   USS(1)   = -11.9062760D0
   GSS(1)   =  12.8480000D0

   ! . Lithium
   SEPAR(3) = .TRUE.
   ALP(3)   =   1.2501400D0
   EISOL(3) =  -5.1280000D0
   BETAS(3) =  -1.3500400D0
   BETAP(3) =  -1.3500400D0
   ZS(3)    =   0.7023800D0
   ZP(3)    =   0.7023800D0
   DD(3)    =   2.0549783D0
   QQ(3)    =   1.7437069D0
   AM(3)    =   0.2682837D0
   AD(3)    =   0.2269793D0
   AQ(3)    =   0.2614581D0
   USS(3)   =  -5.1280000D0
   UPP(3)   =  -2.7212000D0
   GSS(3)   =   7.3000000D0
   GSP(3)   =   5.4200000D0
   GPP(3)   =   5.0000000D0
   GP2(3)   =   4.5200000D0
   HSP(3)   =   0.8300000D0

   ! . Beryllium
   SEPAR(4) = .TRUE.
   ALP(4)   =   1.6694340D0
   EISOL(4) = -24.2047560D0
   BETAS(4) =  -4.0170960D0
   BETAP(4) =  -4.0170960D0
   ZS(4)    =   1.0042100D0
   ZP(4)    =   1.0042100D0
   DD(4)    =   1.4373245D0
   QQ(4)    =   1.2196103D0
   AM(4)    =   0.3307607D0
   AD(4)    =   0.3356142D0
   AQ(4)    =   0.3846373D0
   USS(4)   = -16.6023780D0
   UPP(4)   = -10.7037710D0
   GSS(4)   =   9.0000000D0
   GSP(4)   =   7.4300000D0
   GPP(4)   =   6.9700000D0
   GP2(4)   =   6.2200000D0
   HSP(4)   =   1.2800000D0

   ! . Boron
   SEPAR(5) = .TRUE.
   ALP(5)   =   2.1349930D0
   EISOL(5) = -64.3159500D0
   BETAS(5) =  -8.2520540D0
   BETAP(5) =  -8.2520540D0
   ZS(5)    =   1.5068010D0
   ZP(5)    =   1.5068010D0
   DD(5)    =   0.9579073D0
   QQ(5)    =   0.8128113D0
   AM(5)    =   0.3891951D0
   AD(5)    =   0.4904730D0
   AQ(5)    =   0.5556979D0
   USS(5)   = -34.5471300D0
   UPP(5)   = -23.1216900D0
   GSS(5)   =  10.5900000D0
   GSP(5)   =   9.5600000D0
   GPP(5)   =   8.8600000D0
   GP2(5)   =   7.8600000D0
   HSP(5)   =   1.8100000D0

   ! . Carbon
   SEPAR(6) = .TRUE.
   ALP(6)   =    2.5463800D0
   EISOL(6) = -120.5006060D0
   BETAS(6) =  -18.9850440D0
   BETAP(6) =   -7.9341220D0
   ZS(6)    =    1.7875370D0
   ZP(6)    =    1.7875370D0
   DD(6)    =    0.8074662D0
   QQ(6)    =    0.6851578D0
   AM(6)    =    0.4494671D0
   AD(6)    =    0.6149474D0
   AQ(6)    =    0.6685897D0
   USS(6)   =  -52.2797450D0
   UPP(6)   =  -39.2055580D0
   GSS(6)   =   12.2300000D0
   GSP(6)   =   11.4700000D0
   GPP(6)   =   11.0800000D0
   GP2(6)   =    9.8400000D0
   HSP(6)   =    2.4300000D0

   ! . Nitrogen
   SEPAR(7) = .TRUE.
   ALP(7)   =    2.8613420D0
   EISOL(7) = -202.5812010D0
   BETAS(7) =  -20.4957580D0
   BETAP(7) =  -20.4957580D0
   ZS(7)    =    2.2556140D0
   ZP(7)    =    2.2556140D0
   DD(7)    =    0.6399037D0
   QQ(7)    =    0.5429763D0
   AM(7)    =    0.4994487D0
   AD(7)    =    0.7843643D0
   AQ(7)    =    0.8144720D0
   USS(7)   =  -71.9321220D0
   UPP(7)   =  -57.1723190D0
   GSS(7)   =   13.5900000D0
   GSP(7)   =   12.6600000D0
   GPP(7)   =   12.9800000D0
   GP2(7)   =   11.5900000D0
   HSP(7)   =    3.1400000D0

   ! . Oxygen
   SEPAR(8) = .TRUE.
   ALP(8)   =    3.1606040D0
   EISOL(8) = -317.8685060D0
   BETAS(8) =  -32.6880820D0
   BETAP(8) =  -32.6880820D0
   ZS(8)    =    2.6999050D0
   ZP(8)    =    2.6999050D0
   DD(8)    =    0.5346024D0
   QQ(8)    =    0.4536252D0
   AM(8)    =    0.5667034D0
   AD(8)    =    0.9592562D0
   AQ(8)    =    0.9495934D0
   USS(8)   =  -99.6443090D0
   UPP(8)   =  -77.7974720D0
   GSS(8)   =   15.4200000D0
   GSP(8)   =   14.4800000D0
   GPP(8)   =   14.5200000D0
   GP2(8)   =   12.9800000D0
   HSP(8)   =    3.9400000D0

   ! . Fluorine
   SEPAR(9) = .TRUE.
   ALP(9)   =    3.4196606D0
   EISOL(9) = -476.6837810D0
   BETAS(9) =  -48.2904660D0
   BETAP(9) =  -36.5085400D0
   ZS(9)    =    2.8484870D0
   ZP(9)    =    2.8484870D0
   DD(9)    =    0.5067166D0
   QQ(9)    =    0.4299633D0
   AM(9)    =    0.6218302D0
   AD(9)    =    1.0850301D0
   AQ(9)    =    1.0343643D0
   USS(9)   = -131.0715480D0
   UPP(9)   = -105.7821370D0
   GSS(9)   =   16.9200000D0
   GSP(9)   =   17.2500000D0
   GPP(9)   =   16.7100000D0
   GP2(9)   =   14.9100000D0
   HSP(9)   =    4.8300000D0

   ! . Sodium
   SEPAR(11) = .TRUE.
   ALP(11)   =  1.66D0
   EISOL(11) =  0.00D0
   AM(11)    =  0.50D0

   ! . Aluminium
   SEPAR(13) = .TRUE.
   ALP(13)   =   1.8688394D0
   EISOL(13) = -44.4840720D0
   BETAS(13) =  -2.6702840D0
   BETAP(13) =  -2.6702840D0
   ZS(13)    =   1.4441610D0
   ZP(13)    =   1.4441610D0
   ZD(13)    =   1.0000000D0
   DD(13)    =   1.3992387D0
   QQ(13)    =   1.1586797D0
   AM(13)    =   0.2973172D0
   AD(13)    =   0.2635574D0
   AQ(13)    =   0.3673560D0
   USS(13)   = -23.8070970D0
   UPP(13)   = -17.5198780D0
   GSS(13)   =   8.0900000D0
   GSP(13)   =   6.6300000D0
   GPP(13)   =   5.9800000D0
   GP2(13)   =   5.4000000D0
   HSP(13)   =   0.7000000D0

   ! . Silicon
   SEPAR(14) = .TRUE.
   ALP(14)   =   2.2053160D0
   EISOL(14) = -82.8394220D0
   BETAS(14) =  -9.0868040D0
   BETAP(14) =  -1.0758270D0
   ZS(14)    =   1.3159860D0
   ZP(14)    =   1.7099430D0
   ZD(14)    =   1.0000000D0
   DD(14)    =   1.2580349D0
   QQ(14)    =   0.9785824D0
   AM(14)    =   0.3608967D0
   AD(14)    =   0.3664244D0
   AQ(14)    =   0.4506740D0
   USS(14)   = -37.0375330D0
   UPP(14)   = -27.7696780D0
   GSS(14)   =   9.8200000D0
   GSP(14)   =   8.3600000D0
   GPP(14)   =   7.3100000D0
   GP2(14)   =   6.5400000D0
   HSP(14)   =   1.3200000D0

   ! . Phosphorus
   SEPAR(15) = .TRUE.
   ALP(15)   =    2.4152800D0
   EISOL(15) = -152.9599600D0
   BETAS(15) =   -6.7916000D0
   BETAP(15) =   -6.7916000D0
   ZS(15)    =    2.1087200D0
   ZP(15)    =    1.7858100D0
   ZD(15)    =    1.0000000D0
   DD(15)    =    1.0129699D0
   QQ(15)    =    0.9370090D0
   AM(15)    =    0.4248438D0
   AD(15)    =    0.4882420D0
   AQ(15)    =    0.4979406D0
   USS(15)   =  -56.1433600D0
   UPP(15)   =  -42.8510800D0
   GSS(15)   =   11.5600000D0
   GSP(15)   =   10.0800000D0
   GPP(15)   =    8.6400000D0
   GP2(15)   =    7.6800000D0
   HSP(15)   =    1.9200000D0

   ! . Sulphur
   SEPAR(16) = .TRUE.
   ALP(16)   =    2.4780260D0
   EISOL(16) = -226.0123900D0
   BETAS(16) =  -10.7616700D0
   BETAP(16) =  -10.1084330D0
   ZS(16)    =    2.3129620D0
   ZP(16)    =    2.0091460D0
   ZD(16)    =    1.0000000D0
   DD(16)    =    0.9189935D0
   QQ(16)    =    0.8328514D0
   AM(16)    =    0.4733554D0
   AD(16)    =    0.5544502D0
   AQ(16)    =    0.5585244D0
   USS(16)   =  -72.2422810D0
   UPP(16)   =  -56.9732070D0
   GSS(16)   =   12.8800000D0
   GSP(16)   =   11.2600000D0
   GPP(16)   =    9.9000000D0
   GP2(16)   =    8.8300000D0
   HSP(16)   =    2.2600000D0

   ! . Chlorine
   SEPAR(17) = .TRUE.
   ALP(17)   =    2.5422010D0
   EISOL(17) = -353.1176670D0
   BETAS(17) =  -14.2623200D0
   BETAP(17) =  -14.2623200D0
   ZS(17)    =    3.7846450D0
   ZP(17)    =    2.0362630D0
   ZD(17)    =    1.0000000D0
   DD(17)    =    0.4986870D0
   QQ(17)    =    0.8217603D0
   AM(17)    =    0.5523705D0
   AD(17)    =    0.8061220D0
   AQ(17)    =    0.6053435D0
   USS(17)   = -100.2271660D0
   UPP(17)   =  -77.3786670D0
   GSS(17)   =   15.0300000D0
   GSP(17)   =   13.1600000D0
   GPP(17)   =   11.3000000D0
   GP2(17)   =    9.9700000D0
   HSP(17)   =    2.4200000D0

   ! . Potassium
   SEPAR(19) = .TRUE.
   ALP(19)   =  1.16D0
   EISOL(19) =  0.00D0
   AM(19)    =  0.50D0

   ! . Chromium
   SEPAR(24) = .TRUE.
   ALP(24)   =    3.0683070D0
   EISOL(24) = -134.8187920D0
   BETAS(24) =   -0.1000000D0
   BETAP(24) =   -0.1000000D0
   BETAD(24) =   -8.7766360D0
   ZS(24)    =    1.5000000D0
   ZP(24)    =    1.5000000D0
   ZD(24)    =    2.8845490D0
   DD(24)    =    1.7320508D0
   QQ(24)    =    1.4142136D0
   AM(24)    =    0.2205072D0
   AD(24)    =    0.2711332D0
   AQ(24)    =    0.4464656D0
   USS(24)   =  -17.5170270D0
   UPP(24)   =  -12.5337290D0
   UDD(24)   =  -44.1249280D0
   GSS(24)   =    6.0000000D0
   GSP(24)   =    4.1500000D0
   GPP(24)   =    5.0000000D0
   GP2(24)   =    3.5000000D0
   HSP(24)   =    1.0000000D0
   GSD(24)   =    2.8746410D0
   GPD(24)   =    3.0000000D0
   GDD(24)   =    8.8949670D0

   ! . Zinc
   SEPAR(30) = .TRUE.
   ALP(30)   =    1.5064570D0
   EISOL(30) =  -29.8794320D0
   BETAS(30) =   -1.0000000D0
   BETAP(30) =   -2.0000000D0
   ZS(30)    =    2.0473590D0
   ZP(30)    =    1.4609460D0
   ZD(30)    =    1.0000000D0
   DD(30)    =    1.3037826D0
   QQ(30)    =    1.4520183D0
   AM(30)    =    0.4336641D0
   AD(30)    =    0.2375912D0
   AQ(30)    =    0.2738858D0
   USS(30)   =  -20.8397160D0
   UPP(30)   =  -19.6252240D0
   GSS(30)   =   11.8000000D0
   GSP(30)   =   11.1820180D0
   GPP(30)   =   13.3000000D0
   GP2(30)   =   12.9305200D0
   HSP(30)   =    0.4846060D0

   ! . Germanium
   SEPAR(32) = .TRUE.
   ALP(32)   =   1.9784980D0
   EISOL(32) = -76.2489440D0
   BETAS(32) =  -4.5164790D0
   BETAP(32) =  -1.7555170D0
   ZS(32)    =   1.2931800D0
   ZP(32)    =   2.0205640D0
   DD(32)    =   1.2556091D0
   QQ(32)    =   1.0498655D0
   AM(32)    =   0.3601617D0
   AD(32)    =   0.3643722D0
   AQ(32)    =   0.4347337D0
   USS(32)   = -33.9493670D0
   UPP(32)   = -27.4251050D0
   GSS(32)   =   9.8000000D0
   GSP(32)   =   8.3000000D0
   GPP(32)   =   7.3000000D0
   GP2(32)   =   6.5000000D0
   HSP(32)   =   1.3000000D0

   ! . Bromine
   SEPAR(35) = .TRUE.
   ALP(35)   =    2.44570510D0
   EISOL(35) = -346.68125000D0
   BETAS(35) =   -8.91710700D0
   BETAP(35) =   -9.94374000D0
   ZS(35)    =    3.85430190D0
   ZP(35)    =    2.19920910D0
   ZD(35)    =    1.00000000D0
   DD(35)    =    0.60510740D0
   QQ(35)    =    0.96458730D0
   AM(35)    =    0.55260680D0
   AD(35)    =    0.72583300D0
   AQ(35)    =    0.55745890D0
   USS(35)   =  -99.98644050D0
   UPP(35)   =  -75.67130750D0
   GSS(35)   =   15.03643948D0
   GSP(35)   =   13.03468242D0
   GPP(35)   =   11.27632539D0
   GP2(35)   =    9.85442552D0
   HSP(35)   =    2.45586832D0

   ! . Tin
   SEPAR(50) = .TRUE.
   ALP(50)   =   1.8008140D0
   EISOL(50) = -92.3241020D0
   BETAS(50) =  -3.2351470D0
   BETAP(50) =  -4.2904160D0
   ZS(50)    =   2.0803800D0
   ZP(50)    =   1.9371060D0
   DD(50)    =   1.5697766D0
   QQ(50)    =   1.3262292D0
   AM(50)    =   0.3601617D0
   AD(50)    =   0.3219998D0
   AQ(50)    =   0.3713827D0
   USS(50)   = -40.8518020D0
   UPP(50)   = -28.5602490D0
   GSS(50)   =   9.8000000D0
   GSP(50)   =   8.3000000D0
   GPP(50)   =   7.3000000D0
   GP2(50)   =   6.5000000D0
   HSP(50)   =   1.3000000D0

   ! . Iodine
   SEPAR(53) = .TRUE.
   USS(53)   = -100.00305380D0
   UPP(53)   =  -74.61146920D0
   BETAS(53) =   -7.41445100D0
   BETAP(53) =   -6.19678100D0
   ZS(53)    =    2.27296100D0
   ZP(53)    =    2.16949800D0
   ZD(53)    =    1.00000000D0
   ALP(53)   =    2.20732000D0
   EISOL(53) = -340.59836000D0
   GSS(53)   =   15.04044855D0
   GSP(53)   =   13.05655798D0
   GPP(53)   =   11.14778369D0
   GP2(53)   =    9.91409071D0
   HSP(53)   =    2.45638202D0
   DD(53)    =    1.42532330D0
   QQ(53)    =    1.18417070D0
   AM(53)    =    0.55275410D0
   AD(53)    =    0.45934510D0
   AQ(53)    =    0.45853760D0

   ! . Mercury
   SEPAR(80) = .TRUE.
   ALP(80)   =   1.3356410D0
   EISOL(80) = -28.8191480D0
   BETAS(80) =  -0.4045250D0
   BETAP(80) =  -6.2066830D0
   ZS(80)    =   2.2181840D0
   ZP(80)    =   2.0650380D0
   DD(80)    =   1.7378048D0
   QQ(80)    =   1.4608064D0
   AM(80)    =   0.3969129D0
   AD(80)    =   0.3047694D0
   AQ(80)    =   0.3483102D0
   USS(80)   = -19.8095740D0
   UPP(80)   = -13.1025300D0
   GSS(80)   =  10.8000000D0
   GSP(80)   =   9.3000000D0
   GPP(80)   =  14.3000000D0
   GP2(80)   =  13.5000000D0
   HSP(80)   =   1.3000000D0

   ! . Lead
   SEPAR(82) = .TRUE.
   ALP(82)   =    1.7283330D0
   EISOL(82) = -105.8345040D0
   BETAS(82) =   -8.0423870D0
   BETAP(82) =   -3.0000000D0
   ZS(82)    =    2.4982860D0
   ZP(82)    =    2.0820710D0
   DD(82)    =    1.5526624D0
   QQ(82)    =    1.4488558D0
   AM(82)    =    0.3601617D0
   AD(82)    =    0.3239309D0
   AQ(82)    =    0.3502057D0
   USS(82)   =  -47.3196920D0
   UPP(82)   =  -28.8475600D0
   GSS(82)   =    9.8000000D0
   GSP(82)   =    8.3000000D0
   GPP(82)   =    7.3000000D0
   GP2(82)   =    6.5000000D0
   HSP(82)   =    1.3000000D0
   END SUBROUTINE PARAMETERS_MNDO
  

   SUBROUTINE PARAMETERS_PDDG
   ! . Hydrogen
   SEPAR(1)   = .TRUE.
   USS  (  1) =  -12.8932720D0 
   BETAS(  1) =   -6.1526540D0
   ZS   (  1) =    0.9727860D0
   ALP  (  1) =    3.3816860D0
   EISOL(  1) =  -13.1205660D0
   GSS  (  1) =   14.7942080D0
   AM   (  1) =    0.5437048D0
   AD   (  1) =    0.5437048D0
   AQ   (  1) =    0.5437048D0

   FN1(  1,1) =  1.1222440D0 ; FN2(  1,1) = 4.7077900D0 ; FN3(  1,1) = 1.5470990D0
   FN1(  1,2) = -1.0697370D0 ; FN2(  1,2) = 5.8579950D0 ; FN3(  1,2) = 1.5678930D0
   PDDGC ( 1,1) =  0.0571930D0 ; PDDGE(1,1) =  0.6633950D0
   PDDGC ( 1,2) = -0.0348230D0 ; PDDGE(1,2) =  1.0819010D0

   ! . Carbon
   SEPAR(6)   = .TRUE.
   USS  (  6) =  -48.2412410D0
   UPP  (  6) =  -36.4612560D0
   BETAS(  6) =  -11.9528180D0
   BETAP(  6) =   -9.9224110D0
   ZS   (  6) =    1.5678640D0
   ZP   (  6) =    1.8466590D0
   ALP  (  6) =    2.7257720D0
   EISOL(  6) = -113.4282420D0
   GSS  (  6) =   11.2007080D0
   GSP  (  6) =   10.2650270D0
   GPP  (  6) =   10.7962920D0
   GP2  (  6) =    9.0425660D0
   HSP  (  6) =    2.2909800D0
   DD   (  6) =    0.8314130D0
   QQ   (  6) =    0.6632220D0
   AM   (  6) =    0.4116388D0
   AD   (  6) =    0.5892981D0
   AQ   (  6) =    0.7659489D0

   FN1(  6,1) = 0.0489060D0 ; FN2(  6,1) = 5.7653400D0 ; FN3(  6,1) = 1.6822320D0
   FN1(  6,2) = 0.0476970D0 ; FN2(  6,2) = 5.9737210D0 ; FN3(  6,2) = 0.8944060D0
   PDDGC ( 6,1) = -0.0007430D0 ; PDDGE(6,1) =  0.8369150D0
   PDDGC ( 6,2) =  0.0009850D0 ; PDDGE(6,2) =  1.5852360D0

   ! . Nitrogen.                      
   SEPAR(7)  = .TRUE.
   USS  (  7) =  -49.4545460D0
   UPP  (  7) =  -47.7574060D0
   BETAS(  7) =  -14.1172300D0
   BETAP(  7) =  -19.9385090D0
   ZS   (  7) =    2.0358070D0
   ZP   (  7) =    2.3243270D0
   ALP  (  7) =    2.8491240D0
   EISOL(  7) = -158.4162050D0
   GSS  (  7) =   11.9047870D0
   GSP  (  7) =    7.3485650D0
   GPP  (  7) =   11.7546720D0
   GP2  (  7) =   10.8072770D0
   HSP  (  7) =    1.1367130D0
   DD   (  7) =    0.6548550D0
   QQ   (  7) =    0.5269240D0
   AM   (  7) =    0.4375150D0
   AD   (  7) =    0.5044213D0
   AQ   (  7) =    0.7388754D0

   FN1(  7,1) =  1.5133200D0 ; FN2(  7,1) = 5.9043940D0 ; FN3(  7,1) = 1.7283760D0
   FN1(  7,2) = -1.5118920D0 ; FN2(  7,2) = 6.0300140D0 ; FN3(  7,2) = 1.7341080D0
   PDDGC( 7,1) = -0.0031600D0 ; PDDGE(7,1) =  1.0041720D0
   PDDGC( 7,2) =  0.0125010D0 ; PDDGE(7,2) =  1.5163360D0
 
 
   ! . Oxygen
   SEPAR(8)   = .TRUE.
   USS  (  8) =  -87.4125050D0
   UPP  (  8) =  -72.1830698D0
   BETAS(  8) =  -44.8745530D0
   BETAP(  8) =  -24.6019390D0
   ZS   (  8) =    3.8145650D0
   ZP   (  8) =    2.3180110D0
   ALP  (  8) =    3.2253090D0
   EISOL(  8) = -292.1887660D0
   GSS  (  8) =   15.7557600D0
   GSP  (  8) =   10.6211600D0
   GPP  (  8) =   13.6540160D0
   GP2  (  8) =   12.4060950D0
   HSP  (  8) =    0.5938830D0
   DD   (  8) =    0.4037410D0
   QQ   (  8) =    0.5283600D0
   AM   (  8) =    0.5790428D0
   AD   (  8) =    0.5340363D0
   AQ   (  8) =    0.8009086D0

   FN1(  8,1) = -1.1384550D0 ; FN2(  8,1) = 6.0000430D0 ; FN3(  8,1) = 1.6223620D0
   FN1(  8,2) =  1.1460070D0 ; FN2(  8,2) = 5.9634940D0 ; FN3(  8,2) = 1.6147880D0
   PDDGC( 8,1) = -0.0010000D0 ; PDDGE(8,1) =  1.3606850D0
   PDDGC( 8,2) = -0.0015220D0 ; PDDGE(8,2) =  1.3664070D0

   ! . Fluorine
   SEPAR(9)   = .TRUE.
   USS  (  9) = -111.4004320D0
   UPP  (  9) = -106.3952640D0
   BETAS(  9) =  -50.9373010D0
   BETAP(  9) =  -31.6369760D0
   ZS   (  9) =    5.5380330D0
   ZP   (  9) =    2.5380660D0
   ALP  (  9) =    3.2005710D0
   EISOL(  9) = -442.4571330D0
   GSS  (  9) =   10.4966670D0
   GSP  (  9) =   16.0736890D0
   GPP  (  9) =   14.8172560D0
   GP2  (  9) =   14.4183930D0
   HSP  (  9) =    0.7277630D0
   DD   (  9) =    0.2466010D0
   QQ   (  9) =    0.4825510D0
   AM   (  9) =    0.3857650D0
   AD   (  9) =    0.7878570D0
   AQ   (  9) =    0.6205000D0

   FN1(  9,1) = -0.0080790D0 ; FN2(  9,1) = 5.9389690D0 ; FN3(  9,1) = 1.8639490D0
   FN1(  9,2) = -0.0026590D0 ; FN2(  9,2) = 5.9251050D0 ; FN3(  9,2) = 2.3888640D0
   
   PDDGC(  9,1) =    -0.0128660000D0
   PDDGC(  9,2) =     0.0073150000D0
   PDDGE(  9,1) =     1.3056810000D0
   PDDGE(  9,2) =     1.8425720000D0

   ! . Silicon
   SEPAR( 14)   = .TRUE.
   USS  ( 14)   =   -26.3325220000D0
   UPP  ( 14)   =   -22.6025400000D0
   BETAS( 14)   =    -3.3764450000D0
   BETAP( 14)   =    -3.6209690000D0
   ZS   ( 14)   =     1.5863890000D0
   ZP   ( 14)   =     1.4859580000D0
   ALP  ( 14)   =     2.2151570000D0
   EISOL( 14)   =   -66.8390000000D0
   DD   ( 14)   =     1.3105150000D0
   QQ   ( 14)   =     1.1260890000D0
   AM   ( 14)   =     0.1854904888D0
   AD   ( 14)   =     0.3066060731D0
   AQ   ( 14)   =     0.5267593763D0
   FN1  ( 14,1) =    -0.0713140000D0
   FN2  ( 14,1) =     6.0000000000D0
   FN3  ( 14,1) =     0.2379950000D0
   FN1  ( 14,2) =     0.0894510000D0
   FN2  ( 14,2) =     6.0000000000D0
   FN3  ( 14,2) =     1.8977280000D0
   PDDGC( 14,1) =    -0.0919280000D0
   PDDGC( 14,2) =    -0.0407530000D0
   PDDGE( 14,1) =     1.1631900000D0
   PDDGE( 14,2) =     2.1905260000D0

   GSS  ( 14) =    5.0471960D0
   GSP  ( 14) =    5.9490570D0
   GPP  ( 14) =    6.7593670D0
   GP2  ( 14) =    5.1612970D0
   HSP  ( 14) =    0.9198320D0

   ! . Phosphorus
   SEPAR( 15)   = .TRUE.
   USS  ( 15)   =   -37.8821130000D0
   UPP  ( 15)   =   -30.3129790000D0
   BETAS( 15)   =   -12.6762970000D0
   BETAP( 15)   =    -7.0933180000D0
   ZS   ( 15)   =     2.3958820000D0
   ZP   ( 15)   =     1.7422130000D0
   ALP  ( 15)   =     2.0052940000D0
   EISOL( 15)   =  -117.2128540000D0
   DD   ( 15)   =     0.8939780000D0
   QQ   ( 15)   =     0.9604570000D0
   AM   ( 15)   =     0.2867186201D0
   AD   ( 15)   =     0.4758048477D0
   AQ   ( 15)   =     0.4135967448D0
   FN1  ( 15,1) =    -0.3980550000D0
   FN2  ( 15,1) =     1.9972720000D0
   FN3  ( 15,1) =     0.9500730000D0
   FN1  ( 15,2) =    -0.0796530000D0
   FN2  ( 15,2) =     1.9983600000D0
   FN3  ( 15,2) =     2.3369590000D0
   PDDGC( 15,1) =     0.4627410000D0
   PDDGC( 15,2) =    -0.0204440000D0
   PDDGE( 15,1) =     0.7142960000D0
   PDDGE( 15,2) =     2.0412090000D0

   GSS  ( 15) =    7.8016150D0
   GSP  ( 15) =    5.1869490D0
   GPP  ( 15) =    6.6184780D0
   GP2  ( 15) =    6.0620020D0
   HSP  ( 15) =    1.5428090D0

   ! . Sulfur
   SEPAR( 16)   = .TRUE.
   USS  ( 16)   =   -43.9063660000D0
   UPP  ( 16)   =   -43.4613480000D0
   BETAS( 16)   =    -2.9539120000D0
   BETAP( 16)   =    -8.5077790000D0
   ZS   ( 16)   =     1.0120020000D0
   ZP   ( 16)   =     1.8769990000D0
   ALP  ( 16)   =     2.5397510000D0
   EISOL( 16)   =  -166.3365540000D0
   DD   ( 16)   =     1.0069890000D0
   QQ   ( 16)   =     0.8914870000D0
   AM   ( 16)   =     0.3294621530D0
   AD   ( 16)   =     0.7025708472D0
   AQ   ( 16)   =     0.6628345989D0
   FN1  ( 16,1) =    -0.3306920000D0
   FN2  ( 16,1) =     6.0000000000D0
   FN3  ( 16,1) =     0.8238370000D0
   FN1  ( 16,2) =     0.0241710000D0
   FN2  ( 16,2) =     6.0000000000D0
   FN3  ( 16,2) =     2.0177560000D0
   PDDGC( 16,1) =     0.1204340000D0
   PDDGC( 16,2) =    -0.0026630000D0
   PDDGE( 16,1) =     0.6728700000D0
   PDDGE( 16,2) =     2.0323400000D0

   GSS  ( 16) =    8.9646670D0
   GSP  ( 16) =    6.7859360D0
   GPP  ( 16) =    9.9681640D0
   GP2  ( 16) =    7.9702470D0
   HSP  ( 16) =    4.0418360D0

   ! . Chlorine
   SEPAR(17)  = .TRUE.
   USS  ( 17) =  -95.0944340D0
   UPP  ( 17) =  -53.9216510D0
   BETAS( 17) =  -26.9131290D0
   BETAP( 17) =  -14.9911780D0
   ZS   ( 17) =    2.5482680D0
   ZP   ( 17) =    2.2846240D0
   ZD   ( 17) =    1.0000000D0
   ALP  ( 17) =    2.4976170D0
   EISOL( 17) = -305.7152010D0
   GSS  ( 17) =   16.0136010D0
   GSP  ( 17) =    8.0481150D0
   GPP  ( 17) =    7.5222150D0
   GP2  ( 17) =    7.5041540D0
   HSP  ( 17) =    3.4811530D0
   DD   ( 17) =    0.8275610D0
   QQ   ( 17) =    0.7324270D0
   AM   ( 17) =    0.5885190D0
   AD   ( 17) =    0.7182216D0
   AQ   ( 17) =    0.2174760D0

   FN1( 17,1) = -0.1122220D0 ; FN2( 17,1) = 5.9637190D0 ; FN3( 17,1) = 1.0277190D0
   FN1( 17,2) = -0.0130610D0 ; FN2( 17,2) = 1.9995560D0 ; FN3( 17,2) = 2.2863770D0
   
   PDDGC(  17,1) =    -0.0165520000D0
   PDDGC(  17,2) =    -0.0166460000D0
   PDDGE(  17,1) =     1.7276900000D0
   PDDGE(  17,2) =     1.7846550000D0
   END SUBROUTINE PARAMETERS_PDDG


   SUBROUTINE PARAMETERS_PM3
   ! . Hydrogen
   SEPAR(1)   = .TRUE.
   USS  (  1) =  -13.0733210D0
   BETAS(  1) =   -5.6265120D0
   ZS   (  1) =    0.9678070D0
   ALP  (  1) =    3.3563860D0
   EISOL(  1) =  -13.0733210D0
   GSS  (  1) =   14.7942080D0
   AM   (  1) =    0.5437048D0
   AD   (  1) =    0.5437048D0
   AQ   (  1) =    0.5437048D0

   FN1(  1,1) =  1.1287500D0 ; FN2(  1,1) = 5.0962820D0 ; FN3(  1,1) = 1.5374650D0
   FN1(  1,2) = -1.0603290D0 ; FN2(  1,2) = 6.0037880D0 ; FN3(  1,2) = 1.5701890D0

   ! . Beryllium
   SEPAR(4)   = .TRUE.
   USS  (  4) =  -17.2647520D0
   UPP  (  4) =  -11.3042430D0
   BETAS(  4) =   -3.9620530D0
   BETAP(  4) =   -2.7806840D0
   ZS   (  4) =    0.8774390D0
   ZP   (  4) =    1.5087550D0
   ALP  (  4) =    1.5935360D0
   EISOL(  4) =  -25.5166530D0
   GSS  (  4) =    9.0128510D0
   GSP  (  4) =    6.5761990D0
   GPP  (  4) =    6.0571820D0
   GP2  (  4) =    9.0052190D0
   HSP  (  4) =    0.5446790D0
   DD   (  4) =    1.0090531D0
   QQ   (  4) =    0.8117586D0
   AM   (  4) =    0.3312330D0
   AD   (  4) =    0.2908996D0
   AQ   (  4) =    0.3530008D0

   FN1(  4,1) =  1.6315720D0 ; FN2(  4,1) = 2.6729620D0 ; FN3(  4,1) = 1.7916860D0
   FN1(  4,2) = -2.1109590D0 ; FN2(  4,2) = 1.9685940D0 ; FN3(  4,2) = 1.7558710D0

   ! . Carbon
   SEPAR(6)   = .TRUE.
   USS  (  6) =  -47.2703200D0
   UPP  (  6) =  -36.2669180D0
   BETAS(  6) =  -11.9100150D0
   BETAP(  6) =   -9.8027550D0
   ZS   (  6) =    1.5650850D0
   ZP   (  6) =    1.8423450D0
   ALP  (  6) =    2.7078070D0
   EISOL(  6) = -111.2299170D0
   GSS  (  6) =   11.2007080D0
   GSP  (  6) =   10.2650270D0
   GPP  (  6) =   10.7962920D0
   GP2  (  6) =    9.0425660D0
   HSP  (  6) =    2.2909800D0
   DD   (  6) =    0.8332396D0
   QQ   (  6) =    0.6647750D0
   AM   (  6) =    0.4116394D0
   AD   (  6) =    0.5885862D0
   AQ   (  6) =    0.7647667D0

   FN1(  6,1) = 0.0501070D0 ; FN2(  6,1) = 6.0031650D0 ; FN3(  6,1) = 1.6422140D0
   FN1(  6,2) = 0.0507330D0 ; FN2(  6,2) = 6.0029790D0 ; FN3(  6,2) = 0.8924880D0

   ! . Nitrogen
   SEPAR(7)  = .TRUE.
   USS  (  7) =  -49.3356720D0
   UPP  (  7) =  -47.5097360D0
   BETAS(  7) =  -14.0625210D0
   BETAP(  7) =  -20.0438480D0
   ZS   (  7) =    2.0280940D0
   ZP   (  7) =    2.3137280D0
   ALP  (  7) =    2.8305450D0
   EISOL(  7) = -157.6137755D0
   GSS  (  7) =   11.9047870D0
   GSP  (  7) =    7.3485650D0
   GPP  (  7) =   11.7546720D0
   GP2  (  7) =   10.8072770D0
   HSP  (  7) =    1.1367130D0
   DD   (  7) =    0.6577006D0
   QQ   (  7) =    0.5293383D0
   AM   (  7) =    0.4375151D0
   AD   (  7) =    0.5030995D0
   AQ   (  7) =    0.7364933D0

   FN1(  7,1) =  1.5016740D0 ; FN2(  7,1) = 5.9011480D0 ; FN3(  7,1) = 1.7107400D0
   FN1(  7,2) = -1.5057720D0 ; FN2(  7,2) = 6.0046580D0 ; FN3(  7,2) = 1.7161490D0

   ! . Oxygen
   SEPAR(8)   = .TRUE.
   USS  (  8) =  -86.9930020D0
   UPP  (  8) =  -71.8795800D0
   BETAS(  8) =  -45.2026510D0
   BETAP(  8) =  -24.7525150D0
   ZS   (  8) =    3.7965440D0
   ZP   (  8) =    2.3894020D0
   ALP  (  8) =    3.2171020D0
   EISOL(  8) = -289.3422065D0
   GSS  (  8) =   15.7557600D0
   GSP  (  8) =   10.6211600D0
   GPP  (  8) =   13.6540160D0
   GP2  (  8) =   12.4060950D0
   HSP  (  8) =    0.5938830D0
   DD   (  8) =    0.4086173D0
   QQ   (  8) =    0.5125738D0
   AM   (  8) =    0.5790430D0
   AD   (  8) =    0.5299630D0
   AQ   (  8) =    0.8179630D0

   FN1(  8,1) = -1.1311280D0 ; FN2(  8,1) = 6.0024770D0 ; FN3(  8,1) = 1.6073110D0
   FN1(  8,2) =  1.1378910D0 ; FN2(  8,2) = 5.9505120D0 ; FN3(  8,2) = 1.5983950D0

   ! . Fluorine
   SEPAR(9)   = .TRUE.
   USS  (  9) = -110.4353030D0
   UPP  (  9) = -105.6850470D0
   BETAS(  9) =  -48.4059390D0
   BETAP(  9) =  -27.7446600D0
   ZS   (  9) =    4.7085550D0
   ZP   (  9) =    2.4911780D0
   ALP  (  9) =    3.3589210D0
   EISOL(  9) = -437.5171690D0
   GSS  (  9) =   10.4966670D0
   GSP  (  9) =   16.0736890D0
   GPP  (  9) =   14.8172560D0
   GP2  (  9) =   14.4183930D0
   HSP  (  9) =    0.7277630D0
   DD   (  9) =    0.3125302D0
   QQ   (  9) =    0.4916328D0
   AM   (  9) =    0.3857650D0
   AD   (  9) =    0.6768503D0
   AQ   (  9) =    0.6120047D0

   FN1(  9,1) = -0.0121660D0 ; FN2(  9,1) = 6.0235740D0 ; FN3(  9,1) = 1.8568590D0
   FN1(  9,2) = -0.0028520D0 ; FN2(  9,2) = 6.0037170D0 ; FN3(  9,2) = 2.6361580D0

   ! . Sodium
   SEPAR(11) = .TRUE.
   ALP  (11) = 1.681D0
   EISOL(11) = 0.000D0
   AM   (11) = 0.500D0

   ! . Magnesium
   SEPAR(12) = .TRUE.
   USS  ( 12) =  -14.6236880D0
   UPP  ( 12) =  -14.1734600D0
   BETAS( 12) =   -2.0716910D0
   BETAP( 12) =   -0.5695810D0
   ZS   ( 12) =    0.6985520D0
   ZP   ( 12) =    1.4834530D0
   ALP  ( 12) =    1.3291470D0
   EISOL( 12) =  -22.5530760D0
   GSS  ( 12) =    6.6943000D0
   GSP  ( 12) =    6.7939950D0
   GPP  ( 12) =    6.9104460D0
   GP2  ( 12) =    7.0908230D0
   HSP  ( 12) =    0.5433000D0
   DD   ( 12) =    1.1403950D0
   QQ   ( 12) =    1.1279899D0
   AM   ( 12) =    0.2460235D0
   AD   ( 12) =    0.2695751D0
   AQ   ( 12) =    0.2767522D0

   FN1( 12,1) =  2.1170500D0 ; FN2( 12,1) = 6.0094770D0 ; FN3( 12,1) = 2.0844060D0
   FN1( 12,2) = -2.5477670D0 ; FN2( 12,2) = 4.3953700D0 ; FN3( 12,2) = 2.0636740D0

   ! . Aluminium
   SEPAR(13)  = .TRUE.
   USS  ( 13) =  -24.8454040D0
   UPP  ( 13) =  -22.2641590D0
   BETAS( 13) =   -0.5943010D0
   BETAP( 13) =   -0.9565500D0
   ZS   ( 13) =    1.7028880D0
   ZP   ( 13) =    1.0736290D0
   ZD   ( 13) =    1.0000000D0
   ALP  ( 13) =    1.5217030D0
   EISOL( 13) =  -46.8647630D0
   GSS  ( 13) =    5.7767370D0
   GSP  ( 13) =   11.6598560D0
   GPP  ( 13) =    6.3477900D0
   GP2  ( 13) =    6.1210770D0
   HSP  ( 13) =    4.0062450D0
   DD   ( 13) =    1.2102799D0
   QQ   ( 13) =    1.5585645D0
   AM   ( 13) =    0.2123020D0
   AD   ( 13) =    0.6418584D0
   AQ   ( 13) =    0.2262838D0

   FN1( 13,1) = -0.4730900D0 ; FN2( 13,1) = 1.9158250D0 ; FN3( 13,1) = 1.4517280D0
   FN1( 13,2) = -0.1540510D0 ; FN2( 13,2) = 6.0050860D0 ; FN3( 13,2) = 2.5199970D0

   ! . Silicon
   SEPAR(14)  = .TRUE.
   USS  ( 14) =  -26.7634830D0
   UPP  ( 14) =  -22.8136350D0
   BETAS( 14) =   -2.8621450D0
   BETAP( 14) =   -3.9331480D0
   ZS   ( 14) =    1.6350750D0
   ZP   ( 14) =    1.3130880D0
   ZD   ( 14) =    1.0000000D0
   ALP  ( 14) =    2.1358090D0
   EISOL( 14) =  -67.7882140D0
   GSS  ( 14) =    5.0471960D0
   GSP  ( 14) =    5.9490570D0
   GPP  ( 14) =    6.7593670D0
   GP2  ( 14) =    5.1612970D0
   HSP  ( 14) =    0.9198320D0
   DD   ( 14) =    1.3144550D0
   QQ   ( 14) =    1.2743396D0
   AM   ( 14) =    0.1854905D0
   AD   ( 14) =    0.3060715D0
   AQ   ( 14) =    0.4877432D0

   FN1( 14,1) = -0.3906000D0 ; FN2( 14,1) = 6.0000540D0 ; FN3( 14,1) = 0.6322620D0
   FN1( 14,2) =  0.0572590D0 ; FN2( 14,2) = 6.0071830D0 ; FN3( 14,2) = 2.0199870D0

   ! . Phosphorus
   SEPAR(15)  = .TRUE.
   USS  ( 15) =  -40.4130960D0
   UPP  ( 15) =  -29.5930520D0
   BETAS( 15) =  -12.6158790D0
   BETAP( 15) =   -4.1600400D0
   ZS   ( 15) =    2.0175630D0
   ZP   ( 15) =    1.5047320D0
   ZD   ( 15) =    1.0000000D0
   ALP  ( 15) =    1.9405340D0
   EISOL( 15) = -117.9591740D0
   GSS  ( 15) =    7.8016150D0
   GSP  ( 15) =    5.1869490D0
   GPP  ( 15) =    6.6184780D0
   GP2  ( 15) =    6.0620020D0
   HSP  ( 15) =    1.5428090D0
   DD   ( 15) =    1.0644947D0
   QQ   ( 15) =    1.1120386D0
   AM   ( 15) =    0.2867187D0
   AD   ( 15) =    0.4309446D0
   AQ   ( 15) =    0.3732517D0

   FN1( 15,1) = -0.6114210D0 ; FN2( 15,1) = 1.9972720D0 ; FN3( 15,1) = 0.7946240D0
   FN1( 15,2) = -0.0939350D0 ; FN2( 15,2) = 1.9983600D0 ; FN3( 15,2) = 1.9106770D0

   ! . Sulphur
   SEPAR(16)  = .TRUE.
   USS  ( 16) =  -49.8953710D0
   UPP  ( 16) =  -44.3925830D0
   BETAS( 16) =   -8.8274650D0
   BETAP( 16) =   -8.0914150D0
   ZS   ( 16) =    1.8911850D0
   ZP   ( 16) =    1.6589720D0
   ZD   ( 16) =    1.0000000D0
   ALP  ( 16) =    2.2697060D0
   EISOL( 16) = -183.4537395D0
   GSS  ( 16) =    8.9646670D0
   GSP  ( 16) =    6.7859360D0
   GPP  ( 16) =    9.9681640D0
   GP2  ( 16) =    7.9702470D0
   HSP  ( 16) =    4.0418360D0
   DD   ( 16) =    1.1214313D0
   QQ   ( 16) =    1.0086488D0
   AM   ( 16) =    0.3294622D0
   AD   ( 16) =    0.6679118D0
   AQ   ( 16) =    0.6137472D0

   FN1( 16,1) = -0.3991910D0 ; FN2( 16,1) = 6.0006690D0 ; FN3( 16,1) = 0.9621230D0
   FN1( 16,2) = -0.0548990D0 ; FN2( 16,2) = 6.0018450D0 ; FN3( 16,2) = 1.5799440D0

   ! . Chlorine
   SEPAR(17)  = .TRUE.
   USS  ( 17) = -100.6267470D0
   UPP  ( 17) =  -53.6143960D0
   BETAS( 17) =  -27.5285600D0
   BETAP( 17) =  -11.5939220D0
   ZS   ( 17) =    2.2462100D0
   ZP   ( 17) =    2.1510100D0
   ZD   ( 17) =    1.0000000D0
   ALP  ( 17) =    2.5172960D0
   EISOL( 17) = -315.1949480D0
   GSS  ( 17) =   16.0136010D0
   GSP  ( 17) =    8.0481150D0
   GPP  ( 17) =    7.5222150D0
   GP2  ( 17) =    7.5041540D0
   HSP  ( 17) =    3.4811530D0
   DD   ( 17) =    0.9175856D0
   QQ   ( 17) =    0.7779230D0
   AM   ( 17) =    0.5885190D0
   AD   ( 17) =    0.6814522D0
   AQ   ( 17) =    0.3643694D0

   FN1( 17,1) = -0.1715910D0 ; FN2( 17,1) = 6.0008020D0 ; FN3( 17,1) = 1.0875020D0
   FN1( 17,2) = -0.0134580D0 ; FN2( 17,2) = 1.9666180D0 ; FN3( 17,2) = 2.2928910D0

   ! . Potassium
   SEPAR(19) = .TRUE.
   ALP  (19) = 1.400D0
   EISOL(19) = 0.000D0
   AM   (19) = 0.500D0

    ! . Zinc
    SEPAR(30)  = .TRUE.
    USS  ( 30) =  -18.5321980D0
    UPP  ( 30) =  -11.0474090D0
    BETAS( 30) =   -0.7155780D0
    BETAP( 30) =   -6.3518640D0
    ZS   ( 30) =    1.8199890D0
    ZP   ( 30) =    1.5069220D0
    ZD   ( 30) =    1.0000000D0
    ALP  ( 30) =    1.3501260D0
    EISOL( 30) =  -27.3872000D0
    GSS  ( 30) =    9.6771960D0
    GSP  ( 30) =    7.7362040D0
    GPP  ( 30) =    4.9801740D0
    GP2  ( 30) =    4.6696560D0
    HSP  ( 30) =    0.6004130D0
    DD   ( 30) =    1.5005758D0
    QQ   ( 30) =    1.4077174D0
    AM   ( 30) =    0.3556485D0
    AD   ( 30) =    0.2375689D0
    AQ   ( 30) =    0.2661069D0

    FN1( 30,1) = -0.1112340D0 ; FN2( 30,1) = 6.0014780D0 ; FN3( 30,1) = 1.5160320D0
    FN1( 30,2) = -0.1323700D0 ; FN2( 30,2) = 1.9958390D0 ; FN3( 30,2) = 2.5196420D0

   ! . Zinc (J. Comput. Chem., 25: 1677-1692, 2004, Merz)
   SEPAR(30)  = .TRUE.
   USS  ( 30) =  -16.9746360D0
   UPP  ( 30) =   -9.7941560D0
   BETAS( 30) =   -0.5705140D0
   BETAP( 30) =   -4.1245870D0
   ZS   ( 30) =    1.4341890D0
   ZP   ( 30) =    1.4545820D0
   ZD   ( 30) =    1.0000000D0
   ALP  ( 30) =    1.3602520D0
   EISOL( 30) =  -25.4527010D0
   GSS  ( 30) =    8.4965710D0
   GSP  ( 30) =    8.3019450D0
   GPP  ( 30) =    6.4853290D0
   GP2  ( 30) =    6.2358930D0
   HSP  ( 30) =    0.8023580D0
   DD   ( 30) =    1.7983380D0
   QQ   ( 30) =    1.4583710D0
   AM   ( 30) =    0.3122590D0
   AD   ( 30) =    0.2411480D0
   AQ   ( 30) =    0.2438070D0

   FN1( 30,1) = -0.2621700D0 ; FN2( 30,1) = 4.7309390D0 ; FN3( 30,1) = 1.8022450D0
   FN1( 30,2) = -0.1329170D0 ; FN2( 30,2) = 0.9599290D0 ; FN3( 30,2) = 2.3824630D0

   ! . Gallium
   SEPAR(31)  = .TRUE.
   USS  ( 31) =  -29.8555930D0
   UPP  ( 31) =  -21.8753710D0
   BETAS( 31) =   -4.9456180D0
   BETAP( 31) =   -0.4070530D0
   ZS   ( 31) =    1.8470400D0
   ZP   ( 31) =    0.8394110D0
   ALP  ( 31) =    1.6051150D0
   EISOL( 31) =  -57.3280250D0
   GSS  ( 31) =    8.4585540D0
   GSP  ( 31) =    8.9256190D0
   GPP  ( 31) =    5.0868550D0
   GP2  ( 31) =    4.9830450D0
   HSP  ( 31) =    2.0512600D0
   DD   ( 31) =    0.9776692D0
   QQ   ( 31) =    2.5271534D0
   AM   ( 31) =    0.3108620D0
   AD   ( 31) =    0.5129360D0
   AQ   ( 31) =    0.1546208D0

   FN1( 31,1) = -0.5601790D0 ; FN2( 31,1) = 5.6232730D0 ; FN3( 31,1) = 1.5317800D0
   FN1( 31,2) = -0.2727310D0 ; FN2( 31,2) = 1.9918430D0 ; FN3( 31,2) = 2.1838640D0

   ! . Germanium
   SEPAR(32) = .TRUE.
   USS  ( 32) =  -35.4671955D0
   UPP  ( 32) =  -31.5863583D0
   BETAS( 32) =   -5.3250024D0
   BETAP( 32) =   -2.2501567D0
   ZS   ( 32) =    2.2373526D0
   ZP   ( 32) =    1.5924319D0
   ALP  ( 32) =    1.9723370D0
   EISOL( 32) =  -84.0156006D0
   GSS  ( 32) =    5.3769635D0
   GSP  ( 32) =   10.2095293D0
   GPP  ( 32) =    7.6718647D0
   GP2  ( 32) =    6.9242663D0
   HSP  ( 32) =    1.3370204D0
   DD   ( 32) =    1.1920304D0
   QQ   ( 32) =    1.3321263D0
   AM   ( 32) =    0.1976098D0
   AD   ( 32) =    0.3798182D0
   AQ   ( 32) =    0.3620669D0

   FN1( 32,1) =  0.9631726D0 ; FN2( 32,1) = 6.0120134D0 ; FN3( 32,1) = 2.1633655D0
   FN1( 32,2) = -0.9593891D0 ; FN2( 32,2) = 5.7491802D0 ; FN3( 32,2) = 2.1693724D0

   ! . Arsenic
   SEPAR(33)  = .TRUE.
   USS  ( 33) =  -38.5074240D0
   UPP  ( 33) =  -35.1524150D0
   BETAS( 33) =   -8.2321650D0
   BETAP( 33) =   -5.0173860D0
   ZS   ( 33) =    2.6361770D0
   ZP   ( 33) =    1.7038890D0
   ALP  ( 33) =    1.7944770D0
   EISOL( 33) = -122.6326140D0
   GSS  ( 33) =    8.7890010D0
   GSP  ( 33) =    5.3979830D0
   GPP  ( 33) =    8.2872500D0
   GP2  ( 33) =    8.2103460D0
   HSP  ( 33) =    1.9510340D0
   DD   ( 33) =    0.9679655D0
   QQ   ( 33) =    1.2449874D0
   AM   ( 33) =    0.3230063D0
   AD   ( 33) =    0.5042239D0
   AQ   ( 33) =    0.2574219D0

   FN1( 33,1) = -0.4600950D0 ; FN2( 33,1) = 1.9831150D0 ; FN3( 33,1) = 1.0867930D0
   FN1( 33,2) = -0.0889960D0 ; FN2( 33,2) = 1.9929440D0 ; FN3( 33,2) = 2.1400580D0

   ! . Selenium
   SEPAR(34)  = .TRUE.
   USS  ( 34) =  -55.3781350D0
   UPP  ( 34) =  -49.8230760D0
   BETAS( 34) =   -6.1578220D0
   BETAP( 34) =   -5.4930390D0
   ZS   ( 34) =    2.8280510D0
   ZP   ( 34) =    1.7325360D0
   ALP  ( 34) =    3.0439570D0
   EISOL( 34) = -192.7748115D0
   GSS  ( 34) =    7.4325910D0
   GSP  ( 34) =   10.0604610D0
   GPP  ( 34) =    9.5683260D0
   GP2  ( 34) =    7.7242890D0
   HSP  ( 34) =    4.0165580D0
   DD   ( 34) =    0.8719813D0
   QQ   ( 34) =    1.2244019D0
   AM   ( 34) =    0.2731566D0
   AD   ( 34) =    0.7509697D0
   AQ   ( 34) =    0.5283737D0

   FN1( 34,1) = 0.0478730D0 ; FN2( 34,1) = 6.0074000D0 ; FN3( 34,1) = 2.0817170D0
   FN1( 34,2) = 0.1147200D0 ; FN2( 34,2) = 6.0086720D0 ; FN3( 34,2) = 1.5164230D0

   ! . Bromine
   SEPAR(35)  = .TRUE.
   USS  ( 35) = -116.6193110D0
   UPP  ( 35) =  -74.2271290D0
   BETAS( 35) =  -31.1713420D0
   BETAP( 35) =   -6.8140130D0
   ZS   ( 35) =    5.3484570D0
   ZP   ( 35) =    2.1275900D0
   ZD   ( 35) =    1.0000000D0
   ALP  ( 35) =    2.5118420D0
   EISOL( 35) = -352.5398970D0
   GSS  ( 35) =   15.9434250D0
   GSP  ( 35) =   16.0616800D0
   GPP  ( 35) =    8.2827630D0
   GP2  ( 35) =    7.8168490D0
   HSP  ( 35) =    0.5788690D0
   DD   ( 35) =    0.2759025D0
   QQ   ( 35) =    0.9970532D0
   AM   ( 35) =    0.5859399D0
   AD   ( 35) =    0.6755383D0
   AQ   ( 35) =    0.3823719D0

   FN1( 35,1) =  0.9604580D0 ; FN2( 35,1) = 5.9765080D0 ; FN3( 35,1) = 2.3216540D0
   FN1( 35,2) = -0.9549160D0 ; FN2( 35,2) = 5.9447030D0 ; FN3( 35,2) = 2.3281420D0

   ! . Cadmium
   SEPAR(48)  = .TRUE.
   USS  ( 48) =  -15.8285840D0
   UPP  ( 48) =    8.7497950D0
   BETAS( 48) =   -8.5819440D0
   BETAP( 48) =   -0.6010340D0
   ZS   ( 48) =    1.6793510D0
   ZP   ( 48) =    2.0664120D0
   ALP  ( 48) =    1.5253820D0
   EISOL( 48) =  -22.4502080D0
   GSS  ( 48) =    9.2069600D0
   GSP  ( 48) =    8.2315390D0
   GPP  ( 48) =    4.9481040D0
   GP2  ( 48) =    4.6696560D0
   HSP  ( 48) =    1.6562340D0
   DD   ( 48) =    1.5982681D0
   QQ   ( 48) =    1.2432402D0
   AM   ( 48) =    0.3383668D0
   AD   ( 48) =    0.3570290D0
   AQ   ( 48) =    0.2820582D0

   ! . Indium
   SEPAR(49)  = .TRUE.
   USS  ( 49) =  -26.1762050D0
   UPP  ( 49) =  -20.0058220D0
   BETAS( 49) =   -2.9933190D0
   BETAP( 49) =   -1.8289080D0
   ZS   ( 49) =    2.0161160D0
   ZP   ( 49) =    1.4453500D0
   ALP  ( 49) =    1.4183850D0
   EISOL( 49) =  -51.9750470D0
   GSS  ( 49) =    6.5549000D0
   GSP  ( 49) =    8.2298730D0
   GPP  ( 49) =    6.2992690D0
   GP2  ( 49) =    4.9842110D0
   HSP  ( 49) =    2.6314610D0
   DD   ( 49) =    1.5766241D0
   QQ   ( 49) =    1.7774563D0
   AM   ( 49) =    0.2409004D0
   AD   ( 49) =    0.4532655D0
   AQ   ( 49) =    0.3689812D0

   FN1( 49,1) = -0.3431380D0 ; FN2( 49,1) = 1.9940340D0 ; FN3( 49,1) = 1.6255160D0
   FN1( 49,2) = -0.1095320D0 ; FN2( 49,2) = 5.6832170D0 ; FN3( 49,2) = 2.8670090D0

   ! . Tin
   SEPAR(50) = .TRUE.
   USS  ( 50) =  -34.5501920D0
   UPP  ( 50) =  -25.8944190D0
   BETAS( 50) =   -2.7858020D0
   BETAP( 50) =   -2.0059990D0
   ZS   ( 50) =    2.3733280D0
   ZP   ( 50) =    1.6382330D0
   ALP  ( 50) =    1.6996500D0
   EISOL( 50) =  -78.8877790D0
   GSS  ( 50) =   10.1900330D0
   GSP  ( 50) =    7.2353270D0
   GPP  ( 50) =    5.6738100D0
   GP2  ( 50) =    5.1822140D0
   HSP  ( 50) =    1.0331570D0
   DD   ( 50) =    1.3120038D0
   QQ   ( 50) =    1.5681814D0
   AM   ( 50) =    0.3744959D0
   AD   ( 50) =    0.3218163D0
   AQ   ( 50) =    0.2832529D0

   FN1( 50,1) = -0.1503530D0 ; FN2( 50,1) = 6.0056940D0 ; FN3( 50,1) = 1.7046420D0
   FN1( 50,2) = -0.0444170D0 ; FN2( 50,2) = 2.2573810D0 ; FN3( 50,2) = 2.4698690D0

   ! . Antimony
   SEPAR(51) = .TRUE.
   USS  ( 51) =  -56.4321960D0
   UPP  ( 51) =  -29.4349540D0
   BETAS( 51) =  -14.7942170D0
   BETAP( 51) =   -2.8179480D0
   ZS   ( 51) =    2.3430390D0
   ZP   ( 51) =    1.8999920D0
   ALP  ( 51) =    2.0343010D0
   EISOL( 51) = -148.9382890D0
   GSS  ( 51) =    9.2382770D0
   GSP  ( 51) =    5.2776800D0
   GPP  ( 51) =    6.3500000D0
   GP2  ( 51) =    6.2500000D0
   HSP  ( 51) =    2.4244640D0
   DD   ( 51) =    1.4091903D0
   QQ   ( 51) =    1.3521354D0
   AM   ( 51) =    0.3395177D0
   AD   ( 51) =    0.4589010D0
   AQ   ( 51) =    0.2423472D0

   FN1( 51,1) =  3.0020280D0 ; FN2( 51,1) = 6.0053420D0 ; FN3( 51,1) = 0.8530600D0
   FN1( 51,2) = -0.0188920D0 ; FN2( 51,2) = 6.0114780D0 ; FN3( 51,2) = 2.7933110D0

   ! . Tellurium
   SEPAR(52) = .TRUE.
   USS  ( 52) =  -44.9380360D0
   UPP  ( 52) =  -46.3140990D0
   BETAS( 52) =   -2.6651460D0
   BETAP( 52) =   -3.8954300D0
   ZS   ( 52) =    4.1654920D0
   ZP   ( 52) =    1.6475550D0
   ALP  ( 52) =    2.4850190D0
   EISOL( 52) = -168.0945925D0
   GSS  ( 52) =   10.2550730D0
   GSP  ( 52) =    8.1691450D0
   GPP  ( 52) =    7.7775920D0
   GP2  ( 52) =    7.7551210D0
   HSP  ( 52) =    3.7724620D0
   DD   ( 52) =    0.3484177D0
   QQ   ( 52) =    1.5593085D0
   AM   ( 52) =    0.3768862D0
   AD   ( 52) =    1.1960743D0
   AQ   ( 52) =    0.2184786D0

   FN1( 52,1) =  0.0333910D0 ; FN2( 52,1) = 5.9563790D0 ; FN3( 52,1) = 2.2775750D0 
   FN1( 52,2) = -1.9218670D0 ; FN2( 52,2) = 4.9732190D0 ; FN3( 52,2) = 0.5242430D0

   ! . Iodine
   SEPAR(53)  = .TRUE.
   USS  ( 53) =  -96.4540370D0
   UPP  ( 53) =  -61.0915820D0
   BETAS( 53) =  -14.4942340D0
   BETAP( 53) =   -5.8947030D0
   ZS   ( 53) =    7.0010130D0
   ZP   ( 53) =    2.4543540D0
   ZD   ( 53) =    1.0000000D0
   ALP  ( 53) =    1.9901850D0
   EISOL( 53) = -288.3160860D0
   GSS  ( 53) =   13.6319430D0
   GSP  ( 53) =   14.9904060D0
   GPP  ( 53) =    7.2883300D0
   GP2  ( 53) =    5.9664070D0
   HSP  ( 53) =    2.6300350D0
   DD   ( 53) =    0.1581469D0
   QQ   ( 53) =    1.0467302D0
   AM   ( 53) =    0.5009902D0
   AD   ( 53) =    1.6699104D0
   AQ   ( 53) =    0.5153082D0

   FN1( 53,1) = -0.1314810D0 ; FN2( 53,1) = 5.2064170D0 ; FN3( 53,1) = 1.7488240D0
   FN1( 53,2) = -0.0368970D0 ; FN2( 53,2) = 6.0101170D0 ; FN3( 53,2) = 2.7103730D0

   ! . Mercury
   SEPAR(80)  = .TRUE.
   USS  ( 80) =  -17.7622290D0
   UPP  ( 80) =  -18.3307510D0
   BETAS( 80) =   -3.1013650D0
   BETAP( 80) =   -3.4640310D0
   ZS   ( 80) =    1.4768850D0
   ZP   ( 80) =    2.4799510D0
   ALP  ( 80) =    1.5293770D0
   EISOL( 80) =  -28.8997380D0
   GSS  ( 80) =    6.6247200D0
   GSP  ( 80) =   10.6392970D0
   GPP  ( 80) =   14.7092830D0
   GP2  ( 80) =   16.0007400D0
   HSP  ( 80) =    2.0363110D0
   DD   ( 80) =    1.2317811D0
   QQ   ( 80) =    1.2164033D0
   AM   ( 80) =    0.2434664D0
   AD   ( 80) =    0.4515472D0
   AQ   ( 80) =    0.2618394D0

   FN1( 80,1) =  1.0827200D0 ; FN2( 80,1) = 6.4965980D0 ; FN3( 80,1) = 1.1951460D0
   FN1( 80,2) = -0.0965530D0 ; FN2( 80,2) = 3.9262810D0 ; FN3( 80,2) = 2.6271600D0

   ! . Thallium
   SEPAR(81)  = .TRUE.
   USS  ( 81) =  -30.0531700D0
   UPP  ( 81) =  -26.9206370D0
   BETAS( 81) =   -1.0844950D0
   BETAP( 81) =   -7.9467990D0
   ZS   ( 81) =    6.8679210D0
   ZP   ( 81) =    1.9694450D0
   ALP  ( 81) =    1.3409510D0
   EISOL( 81) =  -56.6492050D0
   GSS  ( 81) =   10.4604120D0
   GSP  ( 81) =   11.2238830D0
   GPP  ( 81) =    4.9927850D0
   GP2  ( 81) =    8.9627270D0
   HSP  ( 81) =    2.5304060D0
   DD   ( 81) =    0.0781362D0
   QQ   ( 81) =    1.5317110D0
   AM   ( 81) =    0.3844326D0
   AD   ( 81) =    2.5741815D0
   AQ   ( 81) =    0.2213264D0

   FN1( 81,1) = -1.3613990D0 ; FN2( 81,1) = 3.5572260D0 ; FN3( 81,1) = 1.0928020D0
   FN1( 81,2) = -0.0454010D0 ; FN2( 81,2) = 2.3069950D0 ; FN3( 81,2) = 2.9650290D0

   ! . Lead
   SEPAR(82)  = .TRUE.
   USS  ( 82) =  -30.3227560D0
   UPP  ( 82) =  -24.4258340D0
   BETAS( 82) =   -6.1260240D0
   BETAP( 82) =   -1.3954300D0
   ZS   ( 82) =    3.1412890D0
   ZP   ( 82) =    1.8924180D0
   ALP  ( 82) =    1.6200450D0
   EISOL( 82) =  -73.4660775D0
   GSS  ( 82) =    7.0119920D0
   GSP  ( 82) =    6.7937820D0
   GPP  ( 82) =    5.1837800D0
   GP2  ( 82) =    5.0456510D0
   HSP  ( 82) =    1.5663020D0
   DD   ( 82) =    0.9866290D0
   QQ   ( 82) =    1.5940562D0
   AM   ( 82) =    0.2576991D0
   AD   ( 82) =    0.4527678D0
   AQ   ( 82) =    0.2150175D0

   FN1( 82,1) = -0.1225760D0 ; FN2( 82,1) = 6.0030620D0 ; FN3( 82,1) = 1.9015970D0
   FN1( 82,2) = -0.0566480D0 ; FN2( 82,2) = 4.7437050D0 ; FN3( 82,2) = 2.8618790D0

   ! . Bismuth
   SEPAR( 83) = .TRUE.
   USS  ( 83) =  -33.4959380D0
   UPP  ( 83) =  -35.5210260D0
   BETAS( 83) =   -5.6072830D0
   BETAP( 83) =   -5.8001520D0
   ZS   ( 83) =    4.9164510D0
   ZP   ( 83) =    1.9349350D0
   ALP  ( 83) =    1.8574310D0
   EISOL( 83) = -109.2774910D0
   GSS  ( 83) =    4.9894800D0
   GSP  ( 83) =    6.1033080D0
   GPP  ( 83) =    8.6960070D0
   GP2  ( 83) =    8.3354470D0
   HSP  ( 83) =    0.5991220D0
   DD   ( 83) =    0.2798609D0
   QQ   ( 83) =    1.5590294D0
   AM   ( 83) =    0.1833693D0
   AD   ( 83) =    0.6776013D0
   AQ   ( 83) =    0.2586520D0

   FN1( 83,1) = 2.5816930D0 ; FN2( 83,1) = 5.0940220D0 ; FN3( 83,1) = 0.4997870D0
   FN1( 83,2) = 0.0603200D0 ; FN2( 83,2) = 6.0015380D0 ; FN3( 83,2) = 2.4279700D0
   END SUBROUTINE PARAMETERS_PM3

   
   SUBROUTINE PARAMETERS_RM1
   ! . Hydrogen
   SEPAR ( 1) = .TRUE.
   ALP   ( 1) =      3.0683595D0
   EISOL ( 1) =    -11.9606770D0
   BETAS ( 1) =     -5.7654447D0
   ZS    ( 1) =      1.0826737D0
   AM    ( 1) =      0.5138998D0
   AD    ( 1) =      0.5138998D0
   AQ    ( 1) =      0.5138998D0
   USS   ( 1) =    -11.9606770D0
   GSS   ( 1) =     13.9832130D0

   FN1( 1, 1) =      0.1028888D0 ; FN2( 1, 1) =      5.9017227D0 ; FN3( 1, 1) =      1.1750118D0
   FN1( 1, 2) =      0.0645745D0 ; FN2( 1, 2) =      6.4178567D0 ; FN3( 1, 2) =      1.9384448D0
   FN1( 1, 3) =     -0.0356739D0 ; FN2( 1, 3) =      2.8047313D0 ; FN3( 1, 3) =      1.6365524D0

   ! . Carbon
   SEPAR ( 6) = .TRUE.
   ALP   ( 6) =      2.7928208D0
   EISOL ( 6) =   -117.8673444D0
   BETAS ( 6) =    -15.4593243D0
   BETAP ( 6) =     -8.2360864D0
   ZS    ( 6) =      1.8501880D0
   ZP    ( 6) =      1.7683009D0
   DD    ( 6) =      0.7967571D0
   QQ    ( 6) =      0.6926111D0
   AM    ( 6) =      0.4797179D0
   AD    ( 6) =      0.5097516D0
   AQ    ( 6) =      0.6614863D0
   USS   ( 6) =    -51.7255603D0
   UPP   ( 6) =    -39.4072894D0
   GSS   ( 6) =     13.0531244D0
   GSP   ( 6) =     11.3347939D0
   GPP   ( 6) =     10.9511374D0
   GP2   ( 6) =      9.7239510D0
   HSP   ( 6) =      1.5521513D0

   FN1( 6, 1) =      0.0746227D0 ; FN2( 6, 1) =      5.7392160D0 ; FN3( 6, 1) =      1.0439698D0
   FN1( 6, 2) =      0.0117705D0 ; FN2( 6, 2) =      6.9240173D0 ; FN3( 6, 2) =      1.6615957D0
   FN1( 6, 3) =      0.0372066D0 ; FN2( 6, 3) =      6.2615894D0 ; FN3( 6, 3) =      1.6315872D0
   FN1( 6, 4) =     -0.0027066D0 ; FN2( 6, 4) =      9.0000373D0 ; FN3( 6, 4) =      2.7955790D0

   ! . Nitrogen
   SEPAR ( 7) = .TRUE.
   ALP   ( 7) =      2.9642254D0
   EISOL ( 7) =   -205.0876419D0
   BETAS ( 7) =    -20.8712455D0
   BETAP ( 7) =    -16.6717185D0
   ZS    ( 7) =      2.3744716D0
   ZP    ( 7) =      1.9781257D0
   DD    ( 7) =      0.6495620D0
   QQ    ( 7) =      0.6191441D0
   AM    ( 7) =      0.4809762D0
   AD    ( 7) =      0.9706271D0
   AQ    ( 7) =      0.8023683D0
   USS   ( 7) =    -70.8512372D0
   UPP   ( 7) =    -57.9773092D0
   GSS   ( 7) =     13.0873623D0
   GSP   ( 7) =     13.2122683D0
   GPP   ( 7) =     13.6992432D0
   GP2   ( 7) =     11.9410395D0
   HSP   ( 7) =      5.0000085D0

   FN1( 7, 1) =      0.0607338D0 ; FN2( 7, 1) =      4.5889295D0 ; FN3( 7, 1) =      1.3787388D0
   FN1( 7, 2) =      0.0243856D0 ; FN2( 7, 2) =      4.6273052D0 ; FN3( 7, 2) =      2.0837070D0
   FN1( 7, 3) =     -0.0228343D0 ; FN2( 7, 3) =      2.0527466D0 ; FN3( 7, 3) =      1.8676382D0

   ! . Oxygen
   SEPAR ( 8) = .TRUE.
   ALP   ( 8) =      4.1719672D0
   EISOL ( 8) =   -312.0403540D0
   BETAS ( 8) =    -29.8510121D0
   BETAP ( 8) =    -29.1510131D0
   ZS    ( 8) =      3.1793691D0
   ZP    ( 8) =      2.5536191D0
   DD    ( 8) =      0.4886701D0
   QQ    ( 8) =      0.4796114D0
   AM    ( 8) =      0.5146059D0
   AD    ( 8) =      1.0065837D0
   AQ    ( 8) =      0.8953404D0
   USS   ( 8) =    -96.9494807D0
   UPP   ( 8) =    -77.8909298D0
   GSS   ( 8) =     14.0024279D0
   GSP   ( 8) =     14.9562504D0
   GPP   ( 8) =     14.1451514D0
   GP2   ( 8) =     12.7032550D0
   HSP   ( 8) =      3.9321716D0

   FN1( 8, 1) =      0.2309355D0 ; FN2( 8, 1) =      5.2182874D0 ; FN3( 8, 1) =      0.9036355D0
   FN1( 8, 2) =      0.0585987D0 ; FN2( 8, 2) =      7.4293293D0 ; FN3( 8, 2) =      1.5175461D0

   ! . Fluorine
   SEPAR ( 9) = .TRUE.
   ALP   ( 9) =      6.0000006D0
   EISOL ( 9) =   -484.5957024D0
   BETAS ( 9) =    -70.0000051D0
   BETAP ( 9) =    -32.6798271D0
   ZS    ( 9) =      4.4033791D0
   ZP    ( 9) =      2.6484156D0
   DD    ( 9) =      0.3488927D0
   QQ    ( 9) =      0.4624444D0
   AM    ( 9) =      0.6145135D0
   AD    ( 9) =      0.9225494D0
   AQ    ( 9) =      0.6236567D0
   USS   ( 9) =   -134.1836959D0
   UPP   ( 9) =   -107.8466092D0
   GSS   ( 9) =     16.7209132D0
   GSP   ( 9) =     16.7614263D0
   GPP   ( 9) =     15.2258103D0
   GP2   ( 9) =     14.8657868D0
   HSP   ( 9) =      1.9976617D0

   FN1( 9, 1) =      0.4030203D0 ; FN2( 9, 1) =      7.2044196D0 ; FN3( 9, 1) =      0.8165301D0
   FN1( 9, 2) =      0.0708583D0 ; FN2( 9, 2) =      9.0000156D0 ; FN3( 9, 2) =      1.4380238D0

   ! . Phosphorus
   SEPAR (15) = .TRUE.
   ALP   (15) =      1.9099329D0
   EISOL (15) =   -123.1797789D0
   BETAS (15) =     -6.1351497D0
   BETAP (15) =     -5.9444213D0
   ZS    (15) =      2.1224012D0
   ZP    (15) =      1.7432795D0
   DD    (15) =      1.0106955D0
   QQ    (15) =      0.9598690D0
   AM    (15) =      0.4072250D0
   AD    (15) =      0.3932223D0
   AQ    (15) =      0.3123478D0
   USS   (15) =    -41.8153318D0
   UPP   (15) =    -34.3834253D0
   GSS   (15) =     11.0805926D0
   GSP   (15) =      5.6833920D0
   GPP   (15) =      7.6041756D0
   GP2   (15) =      7.4026518D0
   HSP   (15) =      1.1618179D0

   FN1(15, 1) =     -0.4106347D0 ; FN2(15, 1) =      6.0875283D0 ; FN3(15, 1) =      1.3165026D0
   FN1(15, 2) =     -0.1629929D0 ; FN2(15, 2) =      7.0947260D0 ; FN3(15, 2) =      1.9072132D0
   FN1(15, 3) =     -0.0488712D0 ; FN2(15, 3) =      8.9997931D0 ; FN3(15, 3) =      2.6585778D0

   ! . Sulphur
   SEPAR (16) = .TRUE.
   ALP   (16) =      2.4401564D0
   EISOL (16) =   -185.3861382D0
   BETAS (16) =     -1.9591072D0
   BETAP (16) =     -8.7743065D0
   ZS    (16) =      2.1334431D0
   ZP    (16) =      1.8746065D0
   DD    (16) =      0.9936921D0
   QQ    (16) =      0.8926247D0
   AM    (16) =      0.4589594D0
   AD    (16) =      0.6930934D0
   AQ    (16) =      0.4959415D0
   USS   (16) =    -55.1677512D0
   UPP   (16) =    -46.5293042D0
   GSS   (16) =     12.4882841D0
   GSP   (16) =      8.5691057D0
   GPP   (16) =      8.5230117D0
   GP2   (16) =      7.6686330D0
   HSP   (16) =      3.8897893D0

   FN1(16, 1) =     -0.7460106D0 ; FN2(16, 1) =      4.8103800D0 ; FN3(16, 1) =      0.5938013D0
   FN1(16, 2) =     -0.0651929D0 ; FN2(16, 2) =      7.2076086D0 ; FN3(16, 2) =      1.2949201D0
   FN1(16, 3) =     -0.0065598D0 ; FN2(16, 3) =      9.0000018D0 ; FN3(16, 3) =      1.8006015D0

   ! . Chlorine
   SEPAR (17) = .TRUE.
   ALP   (17) =      3.6935883D0
   EISOL (17) =   -382.4700938D0
   BETAS (17) =    -19.9243043D0
   BETAP (17) =    -11.5293520D0
   ZS    (17) =      3.8649107D0
   ZP    (17) =      1.8959314D0
   DD    (17) =      0.4541788D0
   QQ    (17) =      0.8825847D0
   AM    (17) =      0.5645068D0
   AD    (17) =      0.7489737D0
   AQ    (17) =      0.7707320D0
   USS   (17) =   -118.4730692D0
   UPP   (17) =    -76.3533034D0
   GSS   (17) =     15.3602310D0
   GSP   (17) =     13.3067117D0
   GPP   (17) =     12.5650264D0
   GP2   (17) =      9.6639708D0
   HSP   (17) =      1.7648990D0

   FN1(17, 1) =      0.1294711D0 ; FN2(17, 1) =      2.9772442D0 ; FN3(17, 1) =      1.4674978D0
   FN1(17, 2) =      0.0028890D0 ; FN2(17, 2) =      7.0982759D0 ; FN3(17, 2) =      2.5000272D0

   ! . Bromine
   SEPAR (35) = .TRUE.
   ALP   (35) =      2.8671053D0
   EISOL (35) =   -357.1164272D0
   BETAS (35) =     -1.3413984D0
   BETAP (35) =     -8.2022599D0
   ZS    (35) =      5.7315721D0
   ZP    (35) =      2.0314758D0
   DD    (35) =      0.2099005D0
   QQ    (35) =      1.0442262D0
   AM    (35) =      0.6290199D0
   AD    (35) =      1.3165755D0
   AQ    (35) =      0.5863476D0
   USS   (35) =   -113.4839818D0
   UPP   (35) =    -76.1872002D0
   GSS   (35) =     17.1156307D0
   GSP   (35) =     15.6241925D0
   GPP   (35) =     10.7354629D0
   GP2   (35) =      8.8605620D0
   HSP   (35) =      2.2351276D0

   FN1(35, 1) =      0.9868994D0 ; FN2(35, 1) =      4.2848419D0 ; FN3(35, 1) =      2.0001970D0
   FN1(35, 2) =     -0.9273125D0 ; FN2(35, 2) =      4.5400591D0 ; FN3(35, 2) =      2.0161770D0

   ! . Iodine
   SEPAR (53) = .TRUE.
   ALP   (53) =      2.1415709D0
   EISOL (53) =   -248.4933241D0
   BETAS (53) =     -4.1931615D0
   BETAP (53) =     -4.4003841D0
   ZS    (53) =      2.5300375D0
   ZP    (53) =      2.3173868D0
   DD    (53) =      1.2963425D0
   QQ    (53) =      1.1085963D0
   AM    (53) =      0.7350144D0
   AD    (53) =      0.3717412D0
   AQ    (53) =      0.3512697D0
   USS   (53) =    -74.8999784D0
   UPP   (53) =    -51.4102380D0
   GSS   (53) =     19.9997413D0
   GSP   (53) =      7.6895767D0
   GPP   (53) =      7.3048834D0
   GP2   (53) =      6.8542461D0
   HSP   (53) =      1.4160294D0

   FN1(53, 1) =     -0.0814772D0 ; FN2(53, 1) =      1.5606507D0 ; FN3(53, 1) =      2.0000206D0
   FN1(53, 2) =      0.0591499D0 ; FN2(53, 2) =      5.7611127D0 ; FN3(53, 2) =      2.2048880D0
   END SUBROUTINE PARAMETERS_RM1

END MODULE MOPAC


subroutine qm3_mopac_setup( nqm, nmm, met, chg, mul, siz, dat, con, cof )
    use mopac 
    implicit none
    integer, intent( in ) :: nqm, nmm, met, chg, mul, siz
    real*8, intent( in ) :: con, cof
    real*8, dimension(0:siz-1), intent( in ) :: dat
    integer :: i, j

    naqm = nqm
    natm = nqm + nmm
    allocate( atmchg(1:natm), atmcrd(1:3,1:natm), atmder(1:3,1:natm) )
    atmchg = 0.0d0
    do i = 1, naqm
        atmchg(i) = dat(i)
    end do
    if( con > .0d0 .and. cof > .0d0 .and. cof > con ) then
        cut_on  = con
        cut_off = cof
    end if
    select case( met )
        case( 0 ) ; call mopac_setup( "MNDO", chg, mul )
        case( 1 ) ; call mopac_setup(  "AM1", chg, mul )
        case( 2 ) ; call mopac_setup(  "RM1", chg, mul )
        case( 3 ) ; call mopac_setup(  "PM3", chg, mul )
        case( 4 ) ; call mopac_setup( "PDDG", chg, mul )
    end select
end subroutine qm3_mopac_setup


subroutine qm3_mopac_calc( mxit, siz, dat )
    use mopac 
    implicit none
    integer, intent( in ) :: mxit, siz
    real*8, dimension(0:siz-1), intent( inout ) :: dat
    integer :: i, j, k
    real*8  :: nuc, scf
    real*8, dimension(1:naqm) :: chg

    ! 1 + 3 * nQM [QM_crd/grd] + nQM [QM_mul] + nMM [MM_chg] + 3 * nMM [MM_crd/grd]
    do i = 1, naqm
        j = 1 + 3 * ( i - 1 )
        atmcrd(1,i) = dat(j)
        atmcrd(2,i) = dat(j+1)
        atmcrd(3,i) = dat(j+2)
    end do
    do i = 1, natm - naqm
        k = naqm + i
        atmchg(k) = dat(4*naqm+i)
        j = 1 + 4 * naqm + ( natm - naqm ) + 3 * ( i - 1 )
        atmcrd(1,k) = dat(j)
        atmcrd(2,k) = dat(j+1)
        atmcrd(3,k) = dat(j+2)
    end do

    atmder = 0.0d0
    call mopac_integrals( nuc, atmder )
    call mopac_scf( scf, mxit, .false. )
    call mopac_gradients( atmder )
    call mopac_charges( chg )

    dat(0) = ev_to_kj * ( scf + nuc ) + atheat
    ! 1 + 3 * nQM [QM_crd/grd] + nQM [QM_mul] + nMM [MM_chg] + 3 * nMM [MM_crd/grd]
    do i = 1, naqm
        k        = 1 + 3 * ( i - 1 )
        dat(k)   = atmder(1,i)
        dat(k+1) = atmder(2,i)
        dat(k+2) = atmder(3,i)
        j        = 3 * naqm + i
        dat(j)   = chg(i)
    end do
    do i = 1, natm - naqm
        k = naqm + i
        j = 1 + 4 * naqm + ( natm - naqm ) + 3 * ( i - 1 )
        dat(j)   = atmder(1,k)
        dat(j+1) = atmder(2,k)
        dat(j+2) = atmder(3,k)
    end do
end subroutine qm3_mopac_calc


subroutine qm3_mopac_clean
    use mopac
    implicit none
    if( allocated( atmcrd ) ) deallocate( atmcrd )
    if( allocated( atmchg ) ) deallocate( atmchg )
    if( allocated( atmder ) ) deallocate( atmder )
    call mopac_data_initialize
end subroutine qm3_mopac_clean


subroutine qm3_mopac_density_write
    use mopac
    implicit none
    call density_write( "mopac.dens" )
end subroutine qm3_mopac_density_write


subroutine qm3_mopac_density_read
    use mopac
    implicit none
    call density_read( "mopac.dens" )
end subroutine qm3_mopac_density_read

