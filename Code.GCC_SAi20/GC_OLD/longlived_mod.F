!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: longlived_mod
!     
! !DESCRIPTION: Module LONGLIVED\_MOD contains routines and variables which
!  are related to long-lived species.
!\\   
!\\   
! !INTERFACE: 
!
      MODULE LONGLIVED_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN
      USE ERROR_MOD,  ONLY : DEBUG_MSG
      USE CMN_SIZE_MOD

#if !defined(ESMF_)
      ! NcdfUtil modules for netCDF I/O
      USE m_netcdf_io_define
      USE m_netcdf_io_create
      USE m_netcdf_io_write

      USE m_netcdf_io_open                    ! netCDF open
      USE m_netcdf_io_get_dimlen              ! netCDF dimension queries
      USE m_netcdf_io_read                    ! netCDF data reads
      USE m_netcdf_io_readattr
      USE m_netcdf_io_close                   ! netCDF close
#endif

      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
#if !defined(ESMF_)
#     include "netcdf.inc"
#endif!  14 Jan 2015 - R. Yantosca -  

      PRIVATE

!
! !PUBLIC DATA MEMBERS:
!
      INTEGER,            PUBLIC :: LLSYEAR ! Year for LLS emissions
      INTEGER,            PUBLIC :: OFFSETYEAR
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: EMISS_BASIC
      PUBLIC  :: LLS_SCALEFAC
      PUBLIC  :: INIT_LLS
      PUBLIC  :: GET_SCALED_CH3Br
      PUBLIC  :: GET_SCALED_CH4
      PUBLIC  :: CLEANUP_LLS
!
! PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: READ_SFC
      PRIVATE :: GET_RCP_NCDF
      PRIVATE :: GET_WMO_NCDF
      PRIVATE :: GET_LLS_MR
      PRIVATE :: SFCMR_INIT
      PRIVATE :: SFCMR_READ
      PRIVATE :: SFCMR_GET
      PRIVATE :: CHECK_LLSYEAR
!
! !REVISION HISTORY: 
!  09 Apr 2015 - S. D. Eastham - Transferred from UCX_MOD and added RCP
!                                functionality
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! MODULE PARAMETERS
      !
      ! MAX_EMIT        : Limit on tracer emission count
      !
      !=================================================================

      INTEGER,  PARAMETER           :: MAX_EMIT=50

#if defined(ESMF_)
      ! Never use NETCDF in ESMF environment (ckeller, 12/05/14). 
      LOGICAL, PARAMETER            :: UCXNETCDF = .FALSE.
#else
      LOGICAL, PARAMETER            :: UCXNETCDF = .TRUE.
#endif

      ! OCS and H2 surface boundary conditions
      REAL(fp),PARAMETER :: OCS_MR = 5.00e-10_fp ! OCS is fixed to 500 pptv
      REAL(fp),PARAMETER :: H2_MR  = 5.00e-7_fp  ! H2 is fixed to 500 ppbv
!
! PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! Scalars
      !
      ! N_EMIT             : Number of tracers to emit
      ! SFC_MONTH          : Current month
      ! OFFSETYEAR         : Number of year-ends since simulation start
      !
      ! Arrays
      !
      ! IDT_EMIT           : IDs of emitted tracers
      ! GRID_EMIT          : Surface CFC mixing ratios
      !
      ! Strings
      !
      ! WMO_FILE_NCDF      : Location of file containing WMO data
      ! RCP_FILE_NCDF      : Location of file containing RCP data
      !
      !=================================================================

      ! Scalars
      INTEGER              :: N_EMIT
      INTEGER              :: SFC_MONTH

      ! Arrays
      INTEGER,DIMENSION(MAX_EMIT)             :: IDT_EMIT
      REAL(fp),DIMENSION(:,:),ALLOCATABLE     :: GRID_EMIT
 
      ! Strings
      CHARACTER(LEN=255)                :: RCP_FILE_NCDF
      CHARACTER(LEN=255)                :: WMO_FILE_NCDF
      CHARACTER(LEN=255)                :: WMO_ASC_FILE_ROOT

      ! Logicals
      LOGICAL                           :: doPrint

      !=================================================================
      ! Derived type object for the surface boundary conditions.
      ! This type holds the surface boundary mixing ratios MR of a given
      ! species (NAME) and for the given dates (YR). The dates are in
      ! decimal years, i.e. 1990.04 refers to (mid of) January 1990. SCAL
      ! is the scale factor to be applied to the MR values, e.g. to 
      ! convert ppt to mixing ratios.
      ! The SFCMR vector will be allocated and filled during the UCX
      ! initialization call (SFCMR_INIT). Surface boundary conditions 
      ! for any given year and month can then be obtained by calling
      ! SFCMR_GET without the need to read the file anymore
      ! (ckeller, 05/12/2014).
      !=================================================================
      TYPE :: SfcMrTyp
         CHARACTER(LEN=31)      :: NAME   ! Species name
         INTEGER                :: ID     ! Species ID 
         INTEGER                :: NX     ! # of values
         REAL(fp)               :: SCAL   ! Scale factor
         REAL(fp), ALLOCATABLE  :: YR(:)  ! year
         REAL(fp), ALLOCATABLE  :: MR(:)  ! mixing ratio
      END TYPE SfcMrTyp

      ! Vector holding all boundary conditions
      TYPE(SfcMrTyp), POINTER   :: SFCMR(:) => NULL()

      ! Maximum number of elements for SFCMR
      INTEGER,        PARAMETER :: MAXSFCMR = 20

      ! Actual number of used elements in SFCMR
      INTEGER                   :: NSFCMR   = 0

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_sfc
!
! !DESCRIPTION: Subroutine READ\_SFC fills out the surface emissions for some
!  species introduced as part of the unified chemistry upgrade.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_SFC( am_I_Root, Input_Opt )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TRACERID_MOD,       ONLY : IDTHCFCX
      USE GRID_MOD,           ONLY : GET_YMID
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TRACERID_MOD,       ONLY : IDTCFCX, IDTCH4,  IDTOCS
      USE TRACERID_MOD,       ONLY : IDTH2,   IDTCH3Br,IDTHCFCX
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT VARIABLES:
!
!
! !REVISION HISTORY: 
!  04 Apr 2013 - S. D. Eastham - Initial version
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt object via the arg list
!  10 Apr 2015 - S. D. Eastham - Added RCP functionality
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)         :: YLAT
      REAL(fp)         :: MR90S,MR30S,MR30N,MR90N
      REAL(fp)         :: MR55S,MR55N
      INTEGER          :: N, TARG_LINE, J
      CHARACTER(LEN=20):: LOCAL_NAME
      INTEGER          :: NCYR, NCMO, TARG_IDT
      CHARACTER(LEN=14):: LUCNAME, UCNAME
      LOGICAL          :: LRCPLLS
 

      !=================================================================
      ! READ_SFC begins here!
      !=================================================================

      ! Re-zero GRID_EMIT
      GRID_EMIT = 0e+0_fp

      NCYR = LLSYEAR+OFFSETYEAR
      NCMO = GET_MONTH()
      DO N=1,N_EMIT
         TARG_IDT = IDT_EMIT(N)
         
         IF (TARG_IDT == IDTCFCX) THEN
            ! Add 3 CFCs
            GRID_EMIT(:,N) = 
     &         GET_LLS_MR('CFC-113',NCYR,NCMO,Input_Opt,am_I_Root) + 
     &         GET_LLS_MR('CFC-114',NCYR,NCMO,Input_Opt,am_I_Root) + 
     &         GET_LLS_MR('CFC-115',NCYR,NCMO,Input_Opt,am_I_Root)
         ELSEIF (TARG_IDT == IDTHCFCX) THEN 
            ! Add 2 HCFCs
            GRID_EMIT(:,N) = 
     &         GET_LLS_MR('HCFC-141b',NCYR,NCMO,Input_Opt,am_I_Root) + 
     &         GET_LLS_MR('HCFC-142b',NCYR,NCMO,Input_Opt,am_I_Root)
         ELSEIF (TARG_IDT == IDTCH4) THEN 
            ! Get the 4 major mixing ratios
            CALL GET_SCALED_CH4(MR90S,MR30S,MR30N,MR90N,
     &                          am_I_Root,Input_Opt)
            ! Figure out which box each is in
            DO J=1,JJPAR   
               ! Grid-box latitude index
               YLAT  = GET_YMID( 1,J,1 )
               IF ( YLAT < -30d0 ) THEN
                  GRID_EMIT(J,N) = MR90S
               ELSE IF ( YLAT >= -30d0 .and. YLAT < 0d0  ) THEN
                  GRID_EMIT(J,N) = MR30S
               ELSE IF ( YLAT >=   0d0 .and. YLAT < 30d0 ) THEN
                  GRID_EMIT(J,N) = MR30N
               ELSE
                  GRID_EMIT(J,N) = MR90N
               ENDIF
            ENDDO
         ELSEIF (TARG_IDT == IDTOCS) THEN 
            ! Simple boundary condition
            GRID_EMIT(:,N) = OCS_MR
         ELSEIF (TARG_IDT == IDTH2) THEN 
            ! Simple boundary condition
            GRID_EMIT(:,N) = H2_MR
         ELSEIF (TARG_IDT == IDTCH3Br) THEN
            ! 4 latitude bands
            CALL GET_SCALED_CH3Br(MR90S,MR55S,MR55N,MR90N,
     &                                am_I_Root,Input_Opt)
            DO J=1,JJPAR   
               ! Grid-box latitude index
               YLAT  = GET_YMID( 1,J,1 )
               IF ( YLAT < -55d0 ) THEN
                  GRID_EMIT(J,N) = MR90S
               ELSE IF ( YLAT >= -55d0 .and. YLAT < 0d0  ) THEN
                  GRID_EMIT(J,N) = MR55S
               ELSE IF ( YLAT >=   0d0 .and. YLAT < 55d0 ) THEN
                  GRID_EMIT(J,N) = MR55N
               ELSE
                  GRID_EMIT(J,N) = MR90N
               ENDIF
            ENDDO
         ELSE

            LUCNAME = Input_Opt%TRACER_NAME(TARG_IDT)
            ! Standard 1:1
            GRID_EMIT(:,N) = 
     &           GET_LLS_MR(TRIM(LUCNAME),NCYR,NCMO,Input_Opt,am_I_Root)
         ENDIF
      ENDDO

      END SUBROUTINE READ_SFC
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_scaled_ch3br
!
! !DESCRIPTION: Subroutine GET\_SCALED\_CH3BR retrieves zonal mean surface mixing
!  ratios for 4 regions (90S-55S,55S-0,0-55N,55N-90N)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_SCALED_CH3Br(MR90S,MR55S,MR55N,MR90N,am_I_Root,
     &                                       Input_Opt)
!
! !USES:
!
      USE TIME_MOD,              ONLY : GET_YEAR
      USE FUTURE_EMISSIONS_MOD,  ONLY : GET_FUTURE_YEAR
      USE GIGC_Input_Opt_Mod,    ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
!
! !OUTPUT VARIABLES:
!
      REAL(fp), INTENT(OUT)       :: MR90S,MR55S,MR55N,MR90N
!
! !REVISION HISTORY: 
!  18 Dec 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! CH3Br values ( from pptv )
      REAL(fp), PARAMETER :: gt55N     = 8.35e-12_fp
      REAL(fp), PARAMETER :: gt0_lt55N = 8.27e-12_fp
      REAL(fp), PARAMETER :: lt0_gt55S = 6.94e-12_fp
      REAL(fp), PARAMETER :: lt55S     = 6.522e-12_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: CH3Br_SCALE, LOWLATSUM, UCX_CH3Br
      INTEGER, SAVE         :: LAST_YEAR=-1
      INTEGER               :: CH3Br_YEAR
      REAL(fp), SAVE        :: LAST_MR(4)
      REAL(fp), SAVE        :: GC_CH3Br

      !=================================================================
      ! GET_SCALED_CH3Br begins here!
      !=================================================================

      ! Calculate RCP adjustment factor
      CALL CHECK_LLSYEAR(am_I_Root,Input_Opt)
      CH3Br_YEAR = LLSYEAR + OFFSETYEAR

      IF (CH3Br_YEAR.ne.LAST_YEAR) THEN
         IF (LAST_YEAR.lt.0) THEN
            ! Calculate mean adjustment
            ! This performs area weighting for 4 segments divided at
            ! +90, +55, 0, -55 and -90
            LOWLATSUM = gt0_lt55N + lt0_gt55S
            GC_CH3Br = (DSIN(55.d0)*(gt55N+lt55S-LOWLATSUM)+LOWLATSUM)/2.d0
         ENDIF
         UCX_CH3Br = GET_LLS_MR('CH3Br',CH3Br_YEAR,1,
     &                                       Input_Opt,am_I_Root)
         CH3Br_SCALE = UCX_CH3Br/GC_CH3Br

         LAST_MR(1) = lt55s*CH3Br_SCALE
         LAST_MR(2) = lt0_gt55S*CH3Br_SCALE
         LAST_MR(3) = gt0_lt55N*CH3Br_SCALE
         LAST_MR(4) = gt55N*CH3Br_SCALE

         ! Debug
         IF ( doPrint ) THEN
            WRITE(6,'(a,e10.4,a)') 'Scaled CH3Br N90: ',
     &                                              LAST_MR(1), ' (v/v)'
            WRITE(6,'(a,e10.4,a)') 'Scaled CH3Br N55: ',
     &                                              LAST_MR(2), ' (v/v)'
            WRITE(6,'(a,e10.4,a)') 'Scaled CH3Br S55: ',
     &                                              LAST_MR(3), ' (v/v)'
            WRITE(6,'(a,e10.4,a)') 'Scaled CH3Br S90: ',
     &                                              LAST_MR(4), ' (v/v)'
         ENDIF

         LAST_YEAR = CH3Br_YEAR
      ENDIF

      MR90S = LAST_MR(1)
      MR55S = LAST_MR(2)
      MR55N = LAST_MR(3)
      MR90N = LAST_MR(4)

      END SUBROUTINE GET_SCALED_CH3Br
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_scaled_ch4
!
! !DESCRIPTION: Subroutine GET\_SCALED\_CH4 retrieves zonal mean surface mixing
!  ratios for 4 equal-area regions (90S-30S,30S-0,0-30N,30N-90N)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_SCALED_CH4(MR90S,MR30S,MR30N,MR90N,am_I_Root,
     &                                     Input_Opt)
!
! !USES:
!
      USE TIME_MOD,              ONLY : GET_YEAR, GET_MONTH
      USE FUTURE_EMISSIONS_MOD,  ONLY : GET_FUTURE_YEAR
      USE GIGC_Input_Opt_Mod,    ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
!
! !OUTPUT VARIABLES:
!
      REAL(fp), INTENT(OUT)       :: MR90S,MR30S,MR30N,MR90N
!
! !REVISION HISTORY: 
!  18 Dec 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: BASE_YEAR, CH4_TARG_YEAR
      INTEGER, PARAMETER    :: REF_YEAR=2007
      REAL(fp)              :: CH4_SCALE
      INTEGER, SAVE         :: LAST_YEAR=-1
      REAL(fp), SAVE        :: LAST_MR(4)
      REAL(fp), SAVE        :: BASE_MR(4)
      INTEGER               :: CH4_YEAR   ! Target year
      LOGICAL               :: WriteCH4, CustomCH4

      !=================================================================
      ! GET_SCALED_CH4 begins here!
      !=================================================================

      ! For CH4 read/write
      WriteCH4 = Input_Opt%LWriteCH4
      CustomCH4= Input_Opt%LCustomCH4

      IF (CustomCH4) THEN
         ! Use the 2007 CH4 gradient and the given target value
         CH4_YEAR = 2007
         BASE_YEAR = 2007
         CH4_TARG_YEAR = 2007
      ELSE
         ! Adapted from GEOS-Chem routine
         CALL CHECK_LLSYEAR(am_I_Root,Input_Opt)
         CH4_YEAR = LLSYEAR + OFFSETYEAR
   
         ! Pick the CH4 concentration [ppbv] for the proper lat bin
         ! CH4 values are read in "chemdr.f"
         IF ( Input_Opt%LFUTURE ) THEN
            CH4_TARG_YEAR = GET_FUTURE_YEAR()
            BASE_YEAR = CH4_TARG_YEAR
         ELSE
            CH4_TARG_YEAR = LLSYEAR + OFFSETYEAR
            IF (Input_Opt%LRCPLLS) THEN
               BASE_YEAR = 2007
            ELSE
               BASE_YEAR = MIN(REF_YEAR,CH4_TARG_YEAR)
            ENDIF
         ENDIF
      ENDIF

      IF (LAST_YEAR .ne. CH4_TARG_YEAR) THEN
         IF (LAST_YEAR.lt.0) THEN
            ! Get CH4 [ppbv] in 4 latitude bins for 2007 (last year with
            ! gradient information)
            CALL GET_GLOBAL_CH4( BASE_YEAR,.TRUE., MR90S,
     &                              MR30S,   MR30N, MR90N,
     &                              am_I_Root, Input_Opt )
            BASE_MR(1) = MR90S
            BASE_MR(2) = MR30S
            BASE_MR(3) = MR30N
            BASE_MR(4) = MR90N
         ELSE
            MR90S = BASE_MR(1)
            MR30S = BASE_MR(2)
            MR30N = BASE_MR(3)
            MR90N = BASE_MR(4)
         ENDIF

         ! Each zone has the same area
         ! Convert from [ppbv] to [vv] and scale by year
         ! Using just one month's data to be consistent with
         ! base GEOS-Chem approach to CH4 (1-yr resolution)

         ! For CH4 read/write
         IF (CustomCH4) THEN
            CH4_SCALE = Input_Opt%TargCH4*1.d-9
         ELSE
            CH4_SCALE = GET_LLS_MR('CH4',LLSYEAR+OFFSETYEAR,1,
     &                 Input_Opt,am_I_Root)
         ENDIF
         CH4_SCALE = CH4_SCALE/(0.25d0*(MR90S+MR30S+MR30N+MR90N))

         MR90S = MR90S * CH4_SCALE
         MR30S = MR30S * CH4_SCALE
         MR30N = MR30N * CH4_SCALE
         MR90N = MR90N * CH4_SCALE

         ! Debug
         IF ( doPrint ) THEN
            WRITE(6,'(a,e10.4,a)') 'Scaled CH4 N90: ', MR90N, ' (v/v)'
            WRITE(6,'(a,e10.4,a)') 'Scaled CH4 N30: ', MR30N, ' (v/v)'
            WRITE(6,'(a,e10.4,a)') 'Scaled CH4 S30: ', MR30S, ' (v/v)'
            WRITE(6,'(a,e10.4,a)') 'Scaled CH4 S90: ', MR90S, ' (v/v)'
         ENDIF

         ! Save for later calls
         LAST_YEAR = CH4_TARG_YEAR
         LAST_MR(1) = MR90S
         LAST_MR(2) = MR30S
         LAST_MR(3) = MR30N
         LAST_MR(4) = MR90N
      ELSE
         ! Re-use data
         MR90S = LAST_MR(1)
         MR30S = LAST_MR(2)
         MR30N = LAST_MR(3)
         MR90N = LAST_MR(4)
      ENDIF

      END SUBROUTINE GET_SCALED_CH4
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lls_scalefac
!
! !DESCRIPTION: Function LLS\_SCALEFAC returns (RCP or WMO) the
!  ratio of the projected surface MR in a target year to that in a
!  reference year.
!\\
!\\
! !INTERFACE:
!
      FUNCTION LLS_SCALEFAC(MS_TRAC,MS_INITYR,MS_TARGMO,
     &                 Input_Opt,am_I_Root) RESULT(SCALEFAC)
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input options
      CHARACTER(*),INTENT(IN)     :: MS_TRAC
      INTEGER,INTENT(IN)          :: MS_INITYR  ! Baseline year
      INTEGER,INTENT(IN)          :: MS_TARGMO  ! Target month
      LOGICAL,INTENT(IN)          :: am_I_Root
!
! !OUTPUT VARIABLES:
!
      REAL*8                      :: SCALEFAC   ! Scaling factor
!
! !REMARKS:
! (1) A remark
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Initial version
!  27 Mar 2014 - S. D. Eastham - Revised for RCP scenarios
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                     :: MS_FUTRYR  ! Future year

      !=================================================================
      ! LLS_SCALEFAC begins here!
      !=================================================================
       
      MS_FUTRYR = LLSYEAR + OFFSETYEAR
      IF (UCXNETCDF) THEN
         SCALEFAC = 
     &    GET_LLS_MR(MS_TRAC,MS_FUTRYR,MS_TARGMO,Input_Opt,am_I_Root)/
     &    GET_LLS_MR(MS_TRAC,MS_INITYR,MS_TARGMO,Input_Opt,am_I_Root)
      ELSE
         SCALEFAC = SFCMR_GET(MS_TRAC,MS_FUTRYR,MS_TARGMO)
     &            / SFCMR_GET(MS_TRAC,MS_INITYR,MS_TARGMO)
      ENDIF

      END FUNCTION LLS_SCALEFAC
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_llsyear
!
! !DESCRIPTION: Subroutine CHECK\_LLSYEAR updates LLSYEAR and checks if a
!  new month has occurred, necessitating a new surface read-in.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_LLSYEAR( am_I_Root, Input_Opt )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,           ONLY : GET_MONTH
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
      LOGICAL, INTENT(IN)           :: am_I_Root
!
! !REMARKS:
! (1) A remark
! 
! !REVISION HISTORY: 
!  15 Jul 2014 - S. D. Eastham - Moved out of READ_SFC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255) :: DBGMSG
      !=================================================================
      ! CHECK_LLSYEAR begins here!
      !=================================================================

      ! If it's a new month, read in new fields
      IF (SFC_MONTH/=GET_MONTH()) THEN
         IF (Input_Opt%LFIXLLS) THEN
            OFFSETYEAR = 0
         ELSEIF (GET_MONTH().lt.SFC_MONTH) THEN
            OFFSETYEAR = OFFSETYEAR + 1
         ENDIF
         SFC_MONTH = GET_MONTH()
         IF ( doPrint ) THEN
            WRITE(DBGMSG,'(a,a,1x,I4)') ' ### LLS: ',
     &         'Using reference year', (LLSYEAR+OFFSETYEAR)
            CALL DEBUG_MSG( DBGMSG )
         ENDIF
            
         CALL READ_SFC( am_I_Root, Input_Opt )
      END IF

      END SUBROUTINE CHECK_LLSYEAR
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_basic
!
! !DESCRIPTION: Subroutine EMISS\_BASIC sets surface mixing ratios of N2O,
!  OCS and ozone-depleting substances covered by the Montreal protocol.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_BASIC( am_I_Root, Input_Opt, State_Met,
     &                        State_Chm  )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState

      ! For CH4 read/write
      USE TIME_MOD,           ONLY : GET_MONTH, GET_YEAR, GET_DAY
      USE TIME_MOD,           ONLY : GET_TS_EMIS, GET_TIME_AHEAD
      USE TIME_MOD,           ONLY : YMD_EXTRACT
      USE GRID_MOD,           ONLY : GET_AREA_M2, GET_YMID
      USE TRACERID_MOD,       ONLY : IDTCH4

!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REMARKS:
!  
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Initial version
!  11 Feb 2014 - R. Yantosca   - Rewrote DO loop to be more efficient
!  14 Feb 2014 - R. Yantosca   - Now parallelize main DO loop
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm 
!                                objects via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, N, N_EMITTED
      REAL(fp)           :: N2O_MULT,OCS_MULT,CFC_MULT,HCFC_MULT
      REAL(fp)           :: DT_YR
      CHARACTER(LEN=255) :: OUTMSG

      ! Local variables for quantities from Input_Opt
      LOGICAL            :: LPRT
      INTEGER            :: N_TRACERS
      REAL(fp)           :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      REAL(fp), POINTER  :: STT(:,:,:,:)

      ! For saving out CH4 data
      LOGICAL            :: CH4EX, CH4READ, CH4WRITE
      REAL(fp),SAVE      :: CH4MASS(4)
      REAL(fp)           :: CH4OLD, CH4NEW
      CHARACTER(LEN=255) :: CH4Root,CH4File
      INTEGER, SAVE      :: LASTMONTH=-1
      INTEGER, SAVE      :: LASTYEAR=-1
      INTEGER            :: CURRMONTH, CURRYEAR, K
      LOGICAL, SAVE      :: FIRST=.TRUE.
      REAL(fp),SAVE      :: TOTALAREA(4)
      REAL(fp)           :: DTSRCE
      INTEGER            :: DTSMIN
      REAL(fp),SAVE      :: TELAPSED=0d0
      LOGICAL            :: LastStepYr, LastStepMo
      REAL(fp)           :: WriteMass, WriteTime
      REAL(fp)           :: TLEFT
      REAL(fp),SAVE      :: TAUMoEND
      INTEGER            :: DStepEnd(2)
      INTEGER            :: YEnd,MEnd,DEnd,fId
      REAL(fp),SAVE      :: CH4Year(12,4)
      LOGICAL            :: NEWFILE
      LOGICAL,PARAMETER  :: WriteAll=.FALSE.
      REAL*8             :: YLAT

      ! For NetCDF I/O
      INTEGER            :: omode, idTime, vId, idBand
      INTEGER            :: var1D(1), ct1D(1), st1D(1)
      INTEGER            :: var2D(2), ct2D(2), st2D(2)

      !=================================================================
      ! EMISS_BASIC begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      LPRT      = Input_Opt%LPRT
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg]
      STT => State_Chm%Tracers

      ! If not correct hour block and day read in new fields
      ! If it's a new month, read in new fields
      CALL CHECK_LLSYEAR( am_I_Root, Input_Opt )

      CH4Read = Input_Opt%LReadCH4
      CH4Write = Input_Opt%LWriteCH4
      CH4Root = Input_Opt%CH4Root

      ! Emission timestep [min]
      DTSMIN = GET_TS_EMIS()
      ! Emission timestep [s]
      DTSRCE = DBLE(DTSMIN) * 60d0
!
      CURRYEAR = GET_YEAR()
      CURRMONTH = GET_MONTH()
!
      IF (FIRST) THEN
         CH4MASS(:) = 0d0
         TELAPSED = 0d0
         CH4YEAR = 0d0
         IF (CH4Read) THEN
            ! Calculate total global grid area
            TOTALAREA = 0d0
            DO J=1,JJPAR
            ! Grid-box latitude index
            YLAT  = GET_YMID( 1,J,1 )
            IF ( YLAT < -30d0 ) THEN
               K = 1
            ELSE IF ( YLAT >= -30d0 .and. YLAT < 0d0  ) THEN
               K = 2
            ELSE IF ( YLAT >=   0d0 .and. YLAT < 30d0 ) THEN
               K = 3
            ELSE
               K = 4
            ENDIF
            DO I=1,IIPAR
               TOTALAREA(K) = TOTALAREA(K) + GET_AREA_M2(I,J,1)
            ENDDO
            ENDDO
         ENDIF
      ENDIF

      NEWFILE = (FIRST.or.(CURRYEAR.ne.LASTYEAR).or.WriteAll)

      IF (NEWFILE) THEN
         ! Target file name
         WRITE(CH4File,'(A,A,I4,A)')
     &       TRIM(CH4Root),'/CH4Emissions.',CURRYEAR,'.nc'
      ENDIF

      IF (FIRST.or.(CURRMONTH.ne.LASTMONTH)) THEN
         IF (CH4Read) THEN
            ! Read in global CH4 emissions rate in kg/s
            IF (NewFile) THEN
               ! Target file name
               WRITE(CH4File,'(A,A,I4,A)')
     &            TRIM(CH4Root),'/CH4Emissions.',CURRYEAR,'.nc'
               CALL Ncop_Rd( fId, CH4File )
               CALL NcRd( CH4Year, fId, 'CH4', (/1,1/), (/12,4/) )
               CALL NcCl( fId )
            ENDIF
            DO I=1,4
                ! Change to kg/timestep/m2
                CH4MASS(I)=CH4YEAR(CURRMONTH,I) * DTSRCE / TOTALAREA(I)
            ENDDO
         ENDIF
      ENDIF

      CH4OLD = 0d0
      CH4NEW = 0d0
  
       ! Loop over emitted species
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N, N_EMITTED, CH4EX, CH4NEW, K )
!$OMP+PRIVATE( YLAT, CH4OLD                         )
!$OMP+SCHEDULE( DYNAMIC )
       DO N = 1, N_EMIT

          ! Get the tracer # of the emitted species
          N_EMITTED = IDT_EMIT(N)

          ! Convert units if the tracer is valied
          IF ( N_EMITTED /= 0 ) THEN
        
             CH4EX = ((N_EMITTED .eq. IDTCH4).and.
     &                (CH4Read .or. CH4Write) )
        
             ! Loop over surface grid boxes
             DO J = 1, JJPAR

             IF (CH4EX) THEN
                ! Grid-box latitude index
                YLAT  = GET_YMID( 1,J,1 )
                IF ( YLAT < -30d0 ) THEN
                   K = 1
                ELSE IF ( YLAT >= -30d0 .and. YLAT < 0d0  ) THEN
                   K = 2
                ELSE IF ( YLAT >=   0d0 .and. YLAT < 30d0 ) THEN
                   K = 3
                ELSE
                   K = 4
                ENDIF
                ! Emissions rate - need to normalize
                ! by gridbox area
                IF (CH4Read) CH4NEW = CH4MASS(K) * GET_AREA_M2(1,J,1)
             ELSE
                ! Be sure - cause an error if this turns up
                K = -1
             ENDIF

             DO I = 1, IIPAR
          
                ! CH4 is special
                IF (CH4EX) THEN
                    ! Which index do we need?
                    IF (CH4Read) THEN
                       STT(I,J,1,N_EMITTED) = STT(I,J,1,N_EMITTED) + 
     &                                          CH4NEW
                    ELSE
                       CH4OLD = STT(I,J,1,N_EMITTED)
                       STT(I,J,1,N_EMITTED) = GRID_EMIT(J,N)
     &                               * State_Met%AD(I,J,1)   
     &                               / TCVV(N_EMITTED)
                       CH4MASS(K) = CH4MASS(K) + 
     &                         (STT(I,J,1,N_EMITTED)-CH4OLD)
                    ENDIF
                ELSE
                    ! Set kg/box using vol/vol (see CONVERT_UNITS in DAO_MOD)
                    STT(I,J,1,N_EMITTED) = GRID_EMIT(J,N)
     &                               * State_Met%AD(I,J,1)   
     &                               / TCVV(N_EMITTED)
                ENDIF

             ENDDO
             ENDDO
          ENDIF
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      ! Write estimated CH4 emissions to file
      IF (CH4Write) THEN
         TELAPSED = TELAPSED + DTSRCE

         ! Is this the last timestep of a month?
         DStepEnd = GET_TIME_AHEAD(DTSMIN+1)
         CALL YMD_EXTRACT(DStepEnd(1),YEnd,MEnd,DEnd)
         LastStepMo = (MEnd.ne.CURRMONTH)

         IF (LastStepMo.or.WriteAll) THEN
            ! Write latest data
            DO K=1,4
                CH4Year(CURRMONTH,K) = CH4MASS(K) / TELAPSED ! kg/s
            ENDDO
            ! Do we need to reset array after writing?
            LASTSTEPYR = (YEnd.ne.CURRYEAR)
            ! Write to NetCDF?
            IF (WriteAll.or.LastStepYr) THEN
               ! Target file name
               WRITE(CH4File,'(A,A,I4,A)')
     &            TRIM(CH4Root),'/CH4Emissions.',CURRYEAR,'.nc'
               INQUIRE( FILE=CH4File, EXIST=CH4EX)
               IF (.not.CH4EX) THEN
               CALL NcCr_Wr( fId, CH4File )
               CALL NcSetFill( fId, NF_NOFILL, omode )

               CALL NcDef_Dimension( fId, 'time', 12, idTime )
               CALL NcDef_Dimension( fId, 'band', 4, idBand )
               CALL NcDef_Glob_Attributes( fId, 'Title','CH4 data' )
               CALL NcDef_Glob_Attributes( fId, 'Format','netCDF-3')
               CALL NcDef_Glob_Attributes( fId, 'valid_range',
     &                            (/ -1e15, +1e15 /))

               ! Time index array
               var1D = (/ idTime /)
               vId   = 0
               CALL NcDef_Variable(fId,'time',NF_INT,1,var1D,vId)
               CALL NcDef_Var_Attributes( fId, vId, 'long_name','time')
               CALL NcDef_Var_Attributes( fId, vId, 'units','months') 

               ! Latitude band index array
               var1D = (/ idBand /)
               vId   = vId+1
               CALL NcDef_Variable(fId,'band',NF_INT,1,var1D,vId)
               CALL NcDef_Var_Attributes( fId, vId, 'long_name',
     &                                            'latitude band')
               CALL NcDef_Var_Attributes( fId, vId, 'units','index') 

               ! CH4 emissions (kg/s)
               var2D = (/ idTime, idBand /)
               vId   = vId+1
               CALL NcDef_Variable(fId,'CH4',NF_FLOAT,2,var2D,vId)
               CALL NcDef_Var_Attributes( fId, vId, 'long_name','CH4')
               CALL NcDef_Var_Attributes( fId, vId, 'units','kg/s') 
               CALL NcDef_Var_Attributes( fId, vId, 'missing_value',   
     &                       1e15)
               CALL NcDef_Var_Attributes( fId, vId, '_FillValue',
     &                       1e15)
               CALL NcDef_Var_Attributes( fId, vId, 'valid_range',     
     &                     (/-1e15, +1e15/)  )

               CALL NcEnd_def( fId )
               st1D = (/ 1  /)
               ct1D = (/ 12 /)
               CALL NcWr( (/1,2,3,4,5,6,7,8,9,10,11,12/),
     &                              fId, 'time',st1D,ct1D)
               st1D = (/ 1 /)
               ct1D = (/ 4 /)
               CALL NcWr( (/1,2,3,4/),
     &                              fId, 'band',st1D,ct1D)
               ELSE
               CALL NcOp_Wr( fId, CH4File )
               ENDIF
               st2D = (/ 1,1 /)
               ct2D = (/ 12,4 /)
               CALL NcWr( CH4Year, fId, 'CH4', st2D, ct2D)
               CALL NcCl( fId )
               IF (LastStepYr) CH4Year = 0d0
            ENDIF
            IF (LastStepMo) THEN
               ! Zero counters
               CH4MASS = 0d0
               TELAPSED = 0d0
            ENDIF
         ENDIF
      ENDIF

      LASTMONTH = CURRMONTH
      LASTYEAR = CURRYEAR
      FIRST = .FALSE.

      IF ( doPrint ) CALL DEBUG_MSG( '### UCX: Emissions complete' )

      END SUBROUTINE EMISS_BASIC
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_lls_mr
!
! !DESCRIPTION: Function GET\_LLS\_MR is a wrapper function which
! retrieves either the WMO or RCP estimate for a surface mixing ratio
! depending on user input switches.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_LLS_MR (TRAC_NAME,TARG_YR_NC,TARG_MO_NC,
     &                  Input_Opt,am_I_Root) RESULT (MEAN_MIXRATIO)
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : TRANUC
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      REAL(fp)                      :: MEAN_MIXRATIO
      CHARACTER*(*),INTENT(IN)      :: TRAC_NAME
      INTEGER,INTENT(IN)            :: TARG_YR_NC
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
      INTEGER,INTENT(IN)            :: TARG_MO_NC
      LOGICAL,INTENT(IN)            :: am_I_Root
!
! !REVISION HISTORY: 
!  27 Mar 2014 - S. D. Eastham - Initial version
!  15 Jul 2014 - S. D. Eastham - Updated for v10-01c
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)             :: TARG_TRAC, DBGMSG

      !=================================================================
      ! GET_LLS_MR begins here!
      !=================================================================

      ! Get a full upper-case version of TRAC_NAME
      TARG_TRAC = TRAC_NAME
      CALL TRANUC( TARG_TRAC )

      ! Are we using RCP or WMO?
      ! In each case, translate cases which are named differently in
      ! GEOS-Chem and in the NetCDF
      IF (Input_Opt%LRCPLLS) THEN
         SELECT CASE ( TRIM( TARG_TRAC ) )
            CASE ( 'CFC11' )
               TARG_TRAC = 'CFC-11'
            CASE ( 'CFC12' )
               TARG_TRAC = 'CFC-12'
            CASE ( 'CFC113' )
               TARG_TRAC = 'CFC-113'
            CASE ( 'CFC114' )
               TARG_TRAC = 'CFC-114'
            CASE ( 'CFC115' )
               TARG_TRAC = 'CFC-115'
            CASE ( 'HCFC22' )
               TARG_TRAC = 'HCFC-22'
            CASE ( 'HCFC141B' )
               TARG_TRAC = 'HCFC-141b'
            CASE ( 'HCFC142B' )
               TARG_TRAC = 'HCFC-142b'
            CASE ( 'CCL4' )
               TARG_TRAC = 'CCl4'
            CASE ( 'H1211' )
               TARG_TRAC = 'H-1211'
            CASE ( 'H1301' )
               TARG_TRAC = 'H-1301'
            CASE ( 'H2402' )
               TARG_TRAC = 'H-2402'
            CASE ( 'C2BR2F4' )
               ! Alias
               TARG_TRAC = 'H-2402'
            CASE ( 'CH3BR' )
               TARG_TRAC = 'CH3Br'
            CASE ( 'CH3CL' )
               TARG_TRAC = 'CH3Cl'
            CASE ( 'CH3CCL3' )
               TARG_TRAC = 'CH3CCl3'
            CASE DEFAULT
               ! Restore case to be the same as before
               TARG_TRAC = TRAC_NAME
         END SELECT
         MEAN_MIXRATIO = GET_RCP_NCDF(TARG_TRAC,TARG_YR_NC,
     &                        TARG_MO_NC,Input_Opt,am_I_Root)
      ELSE
!         WRITE(6,'(a,a)') 'Normal   : ',TRIM(TRAC_NAME)
!         WRITE(6,'(a,a)') 'Uppercase: ',TRIM(TARG_TRAC)
         SELECT CASE ( TRIM( TARG_TRAC ) )
            CASE ( 'CFC-113' )
               TARG_TRAC = 'CFC113'
            CASE ( 'CFC-114' )
               TARG_TRAC = 'CFC114'
            CASE ( 'CFC-115' )
               TARG_TRAC = 'CFC115'
            CASE ( 'HCFC-141B' )
               TARG_TRAC = 'HCFC141b'
            CASE ( 'HCFC-142B' )
               TARG_TRAC = 'HCFC142b'
            CASE ( 'H24O2' )
               TARG_TRAC = 'C2BR2F4'
            CASE ( 'CCL4' )
               TARG_TRAC = 'CCl4'
            CASE ( 'CH3BR' )
               TARG_TRAC = 'CH3Br'
            CASE ( 'CH3CL' )
               TARG_TRAC = 'CH3Cl'
            CASE ( 'CH3CCL3' )
               TARG_TRAC = 'CH3CCl3'
            CASE ( 'H2402' )
               TARG_TRAC = 'C2BR2F4'
            CASE DEFAULT
               ! Restore case to be the same as before
               TARG_TRAC = TRAC_NAME
         END SELECT
         MEAN_MIXRATIO = GET_WMO_NCDF(TARG_TRAC,TARG_YR_NC,
     &          TARG_MO_NC,Input_Opt,am_I_Root)
      ENDIF

      IF ( doPrint ) THEN
         WRITE(DBGMSG,'(a,a,a,E14.6,a)') 'Reading ', TRIM(TRAC_NAME), 
     &     ' MR: ', MEAN_MIXRATIO, ' (v/v)'
         CALL DEBUG_MSG(TRIM(DBGMSG))
      ENDIF

      ! Return to calling program
      END FUNCTION GET_LLS_MR
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_rcp_ncdf
!
! !DESCRIPTION: Function GET\_RCP\_NCDF gets mean surface mixing ratio
!  of a given species based on the current RCP scenario.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_RCP_NCDF (TRAC_NAME,TARG_YR_NC,TARG_MO_NC,Input_Opt,
     &                          am_I_Root) RESULT (MEAN_MIXRATIO)
!
! !USES:
!
      USE FILE_MOD,           ONLY : IOERROR
      USE ERROR_MOD,          ONLY : ERROR_STOP, ALLOC_ERR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      REAL(fp)                       :: MEAN_MIXRATIO
      CHARACTER*(*),INTENT(IN)       :: TRAC_NAME
      INTEGER,INTENT(IN)             :: TARG_YR_NC
      INTEGER,INTENT(IN)             :: TARG_MO_NC 
      TYPE(OptInput), INTENT(IN)     :: Input_Opt   ! Input options
      LOGICAL,INTENT(IN)             :: am_I_Root
!
! !REVISION HISTORY: 
!  27 Mar 2014 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)             :: FILENAME, DBGMSG
      ! For NetCDF
      INTEGER                        :: fId
      INTEGER                        :: st1d(1), ct1d(1)   ! Start + count, for 1D arrays 
      REAL(fp)                       :: ncTemp(1)

      !=================================================================
      ! GET_RCP_NCDF begins here!
      !=================================================================

      ! Read direct from master file
      FILENAME = TRIM(RCP_FILE_NCDF)

      IF ( doPrint ) THEN
         WRITE(DBGMSG,'(a,a,a,a)') ' ### UCX: Reading ', 
     &        TRIM(TRAC_NAME), ' from ', TRIM( FILENAME )
         CALL DEBUG_MSG( TRIM(DBGMSG) )
      ENDIF

      CALL NcOp_Rd (fId,TRIM(FILENAME))

      ! Start and count indices - data is only yearly
      ! May attempt interpolation at a later date
      st1d(1) = TARG_YR_NC-1764
      ct1d(1) = 1

      ! Data is stored as v/v, so no adjustment necessary
      ! Seems to be having difficulties with overloaded function
      CALL NcRd( ncTemp, fId, TRIM(TRAC_NAME), st1d, ct1d )
      MEAN_MIXRATIO = ncTemp(1)
      CALL NcCl( fId )

      ! Return to calling program
      END FUNCTION GET_RCP_NCDF
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_wmo_ncdf
!
! !DESCRIPTION: Function GET\_WMO\_NCDF gets mean surface mixing ratio
!  of a given species based on forecasts from WMO estimates.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_WMO_NCDF (TRAC_NAME,TARG_YR_NC,TARG_MO_NC,Input_Opt,
     &                     am_I_Root)  RESULT (MEAN_MIXRATIO)
!
! !USES:
!
      USE FILE_MOD,           ONLY : IOERROR
      USE ERROR_MOD,          ONLY : ERROR_STOP, ALLOC_ERR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      REAL(fp)                       :: MEAN_MIXRATIO
      CHARACTER*(*),INTENT(IN)       :: TRAC_NAME
      INTEGER,INTENT(IN)             :: TARG_YR_NC
      INTEGER,INTENT(IN)             :: TARG_MO_NC
      TYPE(OptInput), INTENT(IN)     :: Input_Opt   ! Input options
      LOGICAL,INTENT(IN)             :: am_I_Root
!
! !REVISION HISTORY: 
!  04 Apr 2013 - S. D. Eastham - Initial version
!  10 Oct 2013 - S. D. Eastham - Split into two routines (GET and SET)
!  04 Mar 2014 - S. D. Eastham - Added NetCDF compatibility
!  27 Mar 2014 - S. D. Eastham - Restructured into separate module
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)             :: FILENAME, DBGMSG
      ! For NetCDF
      INTEGER                        :: fId, YY, MM
      INTEGER                        :: st1d(1), ct1d(1)   ! Start + count, for 1D arrays 
      REAL*8                         :: ncTemp(1)

      !=================================================================
      ! GET_WMO_NCDF begins here!
      !=================================================================

      IF (UCXNETCDF) THEN
#if !defined(ESMF_)
         ! Read direct from master file
         FILENAME = TRIM(WMO_FILE_NCDF)

         IF ( doPrint ) THEN
            WRITE(DBGMSG,'(a,a,a,a)') ' ### UCX: Reading ', 
     &           TRIM(TRAC_NAME), ' from ', TRIM( FILENAME )
            CALL DEBUG_MSG( TRIM(DBGMSG) )
         ENDIF

         ! By the end of this, need to have data in the form
         ! of an array size JJPARx51 [LAT x LEV]

         ! UCX_MR:      JJPARx51
         ! UCX_MR_IN:   19x51
         CALL NcOp_Rd (fId,TRIM(FILENAME))

         ! Start and count indices
         st1d(1) = (12*(TARG_YR_NC-1959)) + TARG_MO_NC
         ct1d(1) = 1

         ! Data is stored as v/v, so no adjustment necessary
         ! Seems to be having difficulties with overloaded function
         CALL NcRd( ncTemp, fId, TRIM(TRAC_NAME), st1d, ct1d )
         MEAN_MIXRATIO = ncTemp(1)
         CALL NcCl( fId )
#else
         ! This shouldn't happen
         CALL ERROR_STOP('Cannot use NetCDF with ESMF','LLS_MOD')
#endif
      ELSE
         YY = LLSYEAR+OFFSETYEAR
         MM = SFC_MONTH
         MEAN_MIXRATIO = SFCMR_GET(TRAC_NAME,TARG_YR_NC,TARG_MO_NC)
      ENDIF
      
!      IF ( doPrint ) THEN
!         WRITE(DBGMSG,'(a,E16.4)') ' ### UCX: Retrieved value: ', 
!     &        MEAN_MIXRATIO
!         CALL DEBUG_MSG( TRIM(DBGMSG) )
!      ENDIF

      ! Return to calling program
      END FUNCTION GET_WMO_NCDF
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sfcmr_init
!
! !DESCRIPTION: Subroutine SFCMR\_INIT initializes the surface mixing ratio
! derived type object. This reads in all the mixing ratios from ASCII files
! and saves them in the SFCMR object.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SFCMR_INIT( am_I_Root, Input_Opt )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE TRACERID_MOD,       ONLY : IDTCFCX, IDTCH4,  IDTOCS
      USE TRACERID_MOD,       ONLY : IDTH2,   IDTCH3Br,IDTHCFCX
      USE TRACERID_MOD,       ONLY : IDTH2402
!
! !INPUT PARAMETERS:
!
      LOGICAL,          INTENT(IN) :: am_I_Root          ! Root CPU? 
      TYPE(OptInput),   INTENT(IN) :: Input_Opt          ! Input options
!
! !OUTPUT VARIABLES:
!
!
! !REVISION HISTORY: 
!  05 Dec 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I, N, ID, AS
      CHARACTER(LEN=20) :: LOCAL_NAME

      !=================================================================
      ! SFCMR_INIT begins here!
      !=================================================================

      ! Initialize the SFCMR derived type
      ALLOCATE(SFCMR(MAXSFCMR), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SFCMR' )

      ! Fill SFCMR for every emitted species. Some emitted species
      ! have multiple source files, so create a separate SFCMR entry
      ! for all of those.
      DO I = 1, N_EMIT
         ID = IDT_EMIT(I)
         IF (ID.eq.IDTCFCX) THEN
             LOCAL_NAME = TRIM('CFC113')
            CALL SFCMR_READ( am_I_Root, Input_Opt, LOCAL_NAME, ID ) 
            LOCAL_NAME = TRIM('CFC114')
            CALL SFCMR_READ( am_I_Root, Input_Opt, LOCAL_NAME, ID ) 
            LOCAL_NAME = TRIM('CFC115')
            CALL SFCMR_READ( am_I_Root, Input_Opt, LOCAL_NAME, ID ) 
         ELSEIF (ID.eq.IDTHCFCX) THEN
            LOCAL_NAME = TRIM('HCFC141b')
            CALL SFCMR_READ( am_I_Root, Input_Opt, LOCAL_NAME, ID ) 
            LOCAL_NAME = TRIM('HCFC142b')
            CALL SFCMR_READ( am_I_Root, Input_Opt, LOCAL_NAME, ID ) 
         ELSEIF (ID.eq.IDTH2402) THEN
            LOCAL_NAME = TRIM('C2BR2F4')
            CALL SFCMR_READ( am_I_Root, Input_Opt, LOCAL_NAME, ID )
         ELSEIF (ID.eq.IDTOCS) THEN
            LOCAL_NAME = TRIM('OCS')
            CALL SFCMR_READ( am_I_Root,  Input_Opt, 
     &                       LOCAL_NAME, ID,        iMR=OCS_MR )
         ELSEIF (ID.eq.IDTH2) THEN
            LOCAL_NAME = TRIM('H2')
            CALL SFCMR_READ( am_I_Root,  Input_Opt, 
     &                       LOCAL_NAME, ID,        iMR=H2_MR )
         ELSE
               LOCAL_NAME = TRIM(Input_Opt%TRACER_NAME(ID))
               CALL SFCMR_READ( am_I_Root, Input_Opt, LOCAL_NAME, ID ) 
         ENDIF
      ENDDO

      END SUBROUTINE SFCMR_INIT
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sfcmr_read
!
! !DESCRIPTION: Subroutine SFCMR\_READ reads the surface mixing rations from
! ASCII file and stores the values in the next available SFCMR object.
! The ascii file is expected to be in the following format:
!C2CL3F3              MIXING RATIO BOUNDARY CONDITIONS FOR WMO-2006
!1.0E-12
!
! time      CFC113
!1959.04        0.99
!1959.12        1.00
!1959.21        1.01
!1959.29        1.02
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SFCMR_READ( am_I_Root, Input_Opt, iName, iID, iMR )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE FILE_MOD,           ONLY : IOERROR
!
! !INPUT PARAMETERS:
!
      LOGICAL,            INTENT(IN   ) :: am_I_Root          ! Root CPU? 
      TYPE(OptInput),     INTENT(IN   ) :: Input_Opt          ! Input options
      CHARACTER(LEN=*),   INTENT(IN   ) :: iName              ! species name 
      INTEGER,            INTENT(IN   ) :: iID                ! species ID 
      REAL(fp), OPTIONAL, INTENT(IN   ) :: iMR                ! fixed MR
!
! !REVISION HISTORY: 
!  05 Dec 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, N, NX, IOS, AS
      INTEGER            :: IU_FILE
      CHARACTER(LEN=30)  :: FM 
      ! Max. number of lines in ascii files. Monthly values btw 1950-2100 = 1800
      INTEGER, PARAMETER :: MAXLINES = 3000
      REAL(fp)           :: SCALEFAC, TMPYR, TMPMR
      REAL(fp)           :: YR(MAXLINES), MR(MAXLINES)
      CHARACTER(LEN=255) :: FILENAME, DBGMSG

      !=================================================================
      ! SFCMR_READ begins here!
      !=================================================================

      ! If a fixed MR is given, just use this value
      IF ( PRESENT( iMR ) ) THEN
         NX       = 1
         YR(1)    = 1900.0_fp
         MR(1)    = iMR
         SCALEFAC = 1.0_fp

      ! If no fixed MR is given, read MR and years from ASCII file
      ELSE

         ! Construct file name
         FILENAME = TRIM(WMO_ASC_FILE_ROOT)//TRIM(iName)//'_boundary'
        
         IF ( doPrint ) THEN
            WRITE(DBGMSG,'(a,a,a,a)') ' ### LLS: Reading ',
     &           TRIM(iNAME), ' from ', TRIM( FILENAME )
            CALL DEBUG_MSG( TRIM(DBGMSG) )
         ENDIF      
   
         ! Get a free LUN
         IU_FILE = findFreeLUN()      
   
         ! Initialize values
         IOS = 1
         NX  = 0  ! number of values
   
         ! Open file
         OPEN( IU_FILE, FILE=TRIM( FILENAME ), IOSTAT=IOS )
         IF ( IOS /= 0 ) THEN
            WRITE(6,*) 'LLS: FAILED TO FIND ', TRIM(FILENAME)
             CALL IOERROR( IOS, IU_FILE, 'LLS_MOD:SFCMR_READ')
         ENDIF      
   
         ! First line not important
         READ(IU_FILE,*,IOSTAT=IOS)
         IF ( IOS /= 0 ) THEN
            WRITE(6,*) 'LLS: READ ERROR IN ', TRIM(FILENAME)
             CALL IOERROR( IOS, IU_FILE, 'LLS_MOD:SFCMR_READ')
         ENDIF
         
         ! Second line contains scaling factor
         READ(IU_FILE,101,IOSTAT=IOS) SCALEFAC
   
101      FORMAT(E7.1)
   
         ! Third and fourth line also to be ignored
         DO I=1,2
            READ(IU_FILE,*,IOSTAT=IOS)
            IF ( IOS /= 0 ) THEN
               WRITE(6,*) 'LLS: READ ERROR IN ', TRIM(FILENAME)
                CALL IOERROR( IOS, IU_FILE, 'LLS_MOD:SFCMR_READ')
            ENDIF
         ENDDO
  
         ! C2BR2F4 has a different format.
         IF ( TRIM(iName) == 'C2BR2F4' ) THEN
            FM = "(F7.2,F8.2)"
         ELSE
            !FM = "(F7.2,F11.2)"
            FM = "(F7.2,F12.2)"
         ENDIF
   
         ! Read until end of file 
         DO
            ! Read values into temporary buffer
            READ(IU_FILE,FM,IOSTAT=IOS) TMPYR, TMPMR
   
            ! Error?
            IF ( IOS > 0 ) THEN
               WRITE(6,*) 'UCX: READ ERROR IN ', TRIM(FILENAME)
                CALL IOERROR( IOS, IU_FILE, 'UCX_MOD:SFCMR_READ')
            ENDIF
   
            ! End of file?
            IF ( IOS < 0 ) EXIT
   
            ! Increase counter and pass data to buffer
            NX = NX + 1
            IF ( NX > MAXLINES ) THEN
               DBGMSG = 'Beyond max. number of lines: '//TRIM(FILENAME)
               CALL ERROR_STOP( DBGMSG, 'UCX_MOD: SFCMR_READ' )
            ENDIF
            YR(NX) = TMPYR
            MR(NX) = TMPMR
         ENDDO

         ! Close file
         CLOSE(IU_FILE)

      ENDIF

      ! Get next SFCMR index in list
      NSFCMR = NSFCMR + 1
      IF ( NSFCMR > MAXSFCMR ) THEN
         DBGMSG = 'Max. number of SFCMR objects exceeded! '
         CALL ERROR_STOP( DBGMSG, 'UCX_MOD: SFCMR_READ' )
      ENDIF
      N = NSFCMR

      ! Fill SFCMR type object
      ALLOCATE( SFCMR(N)%YR(NX), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SFCMR(N)%YR'//TRIM(iName) )
      ALLOCATE( SFCMR(N)%MR(NX), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SFCMR(N)%MR'//TRIM(iName) )

      SFCMR(N)%NAME  = TRIM(iName)
      SFCMR(N)%ID    = iID 
      SFCMR(N)%NX    = NX
      SFCMR(N)%YR(:) = YR(1:NX)
      SFCMR(N)%MR(:) = MR(1:NX)
      SFCMR(N)%SCAL  = SCALEFAC

      END SUBROUTINE SFCMR_READ
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sfcmr_get
!
! !DESCRIPTION: Subroutine SFCMR\_GET returns the surface mixing ratio of
! the given species for the given year and month. The value of the closest
! available date in the past is taken, and no interpolation between dates
! is performed.
!\\
!\\
! !INTERFACE:
!
      FUNCTION SFCMR_GET( TRAC_NAME, YEAR, MONTH ) RESULT( MR ) 
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      CHARACTER(LEN=*), INTENT(IN) :: TRAC_NAME          ! species name
      INTEGER,          INTENT(IN) :: YEAR               ! desired year 
      INTEGER,          INTENT(IN) :: MONTH              ! desired month
!
! !INPUT/OUTPUT VARIABLES:
!
      REAL(fp)                     :: MR 
!
! !REVISION HISTORY: 
!  05 Dec 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, N, NX, IDX
      REAL(fp)            :: YRMT
      REAL(fp), PARAMETER :: DECMT = 1.0_fp / 12.0_fp
      REAL(fp), PARAMETER :: MIDMT = DECMT  /  2.0_fp
      CHARACTER(LEN=255)  :: MSG

      !=================================================================
      ! SFCMR_GET begins here!
      !=================================================================

      ! Init
      MR = 0.0_fp

      ! In the SFCMR structure, the date is stored as decimal year/month
      ! value referring to the mid-month, e.g. January 1990 is 1990.04, 
      ! July 1990 is 1990.54, etc.
      YRMT = REAL(YEAR,4) + ( (MONTH-1) * DECMT ) + MIDMT
  
      ! Get index to desired SFCMR object
      N = -1
      DO I = 1, NSFCMR
         IF ( TRIM(SFCMR(I)%NAME) == TRIM(TRAC_NAME) ) THEN
            N = I
            EXIT
         ENDIF
      ENDDO
      IF ( N < 0 ) RETURN

      ! Get time index to be used. Use closest available time stamp
      ! that is in the past.
      IDX = -1
      NX  = SFCMR(N)%NX

      ! Check for lower boundary
      IF ( YRMT <= SFCMR(N)%YR(1) ) THEN
         IDX = 1

      ! Check for upper boundary
      ELSEIF ( YRMT >= SFCMR(N)%YR(NX) ) THEN
         IDX = NX

      ! Otherwise walk through dates until we pass the date line
      ELSE
         DO I = 1, NX-1
            IF ( SFCMR(N)%YR(I+1) > YRMT ) THEN
               IDX = I
               EXIT
            ENDIF
         ENDDO
      ENDIF

      ! This should never happen:
      IF ( IDX < 0 ) RETURN

      ! Get mixing ratio
      MR = SFCMR(N)%MR(IDX) * SFCMR(N)%SCAL 

      IF ( doPrint ) THEN
         WRITE(MSG,'(a,E16.4)') ' ### UCX: Retrieved value: ', MR
         CALL DEBUG_MSG( TRIM(MSG) )
      ENDIF      

      END FUNCTION SFCMR_GET
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lls
!
! !DESCRIPTION: Subroutine INIT\_LLS initializes module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_LLS( am_I_Root, Input_Opt )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,           ONLY : GET_YEAR
      USE TRACERID_MOD
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN) :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
!
! !OUTPUT VARIABLES:
!
!
! !REVISION HISTORY: 
!  04 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N, AS
      CHARACTER(LEN=255) :: DBGMSG

      ! Local variables for quantities from Input_Opt
      LOGICAL :: LPRT
      LOGICAL :: LUCX
      LOGICAL :: LFUTURECFC
      LOGICAL :: LCH4SBC
      LOGICAL :: LCFCEMIS
      LOGICAL :: LCLEMIS
      LOGICAL :: LOCSEMIS
      LOGICAL :: LN2OEMIS
      LOGICAL :: LBREMIS
      LOGICAL :: CH4Read
!      LOGICAL :: LBASICEMIS
      CHARACTER(LEN=255) :: RCPNAME

      ! Return code
      INTEGER :: RC

      !=================================================================
      ! INIT_LLS begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      LPRT         = ( Input_Opt%LPRT .and. am_I_Root )
      LUCX         = Input_Opt%LUCX
      LFUTURECFC   = Input_Opt%LFUTURECFC
      LCH4SBC      = Input_Opt%LCH4SBC
      LCFCEMIS     = Input_Opt%LCFCEMIS
      LCLEMIS      = Input_Opt%LCLEMIS
      LOCSEMIS     = Input_Opt%LOCSEMIS
      LN2OEMIS     = Input_Opt%LN2OEMIS
      LBREMIS      = Input_Opt%LBREMIS
      RCPNAME      = Input_Opt%RCPNAME
!      LBASICEMIS   = Input_Opt%LBASICEMIS
      CH4Read      = Input_Opt%LReadCH4

      ! Set module variable doPrint for use in internal functions
      doPrint      = LPRT .AND. am_I_Root

      IF ( am_I_Root ) THEN
         WRITE( 6,'(a)') REPEAT( '=', 79 )
         WRITE( 6,'(a)') 'L O N G   L I V E D   S P E C I E S'
         WRITE( 6,'(a)') 'Routines written by SEBASTIAN D. EASTHAM'
         WRITE( 6,'(a)') REPEAT( '=', 79 )
      ENDIF

      ! --------------------------------------------------------------
      ! Input data sources
      ! --------------------------------------------------------------

      ! Determine folder paths from root folder
      ! Regridding of netCDF input data is currently not supported.
      IF (UCXNETCDF) THEN
         WRITE( WMO_FILE_NCDF,'(a,a)') 
     &        TRIM(Input_Opt%CHEM_INPUTS_DIR), 
     &        'UCX_201403/Init1D/Montreal.nc'
         IF (Input_Opt%LRCPLLS) THEN
            WRITE(RCP_FILE_NCDF,'(a,a,a,a)')
     &        TRIM(Input_Opt%CHEM_INPUTS_DIR), 
     &         '/UCX_201403/Init1D/',TRIM(RCPNAME),'.nc'
!            ! Overwrite LLSYEAR with RCPYEAR
!            LLSYEAR = RCPYEAR
            WRITE(DBGMSG,'(a,a)') 
     &        '### UCX: Setting long-lived species BCs from ',
     &        TRIM(RCP_FILE_NCDF)
         ELSE
            IF (.not.LFUTURECFC) THEN
               LLSYEAR = GET_YEAR()
            ENDIF
            RCP_FILE_NCDF=''
            WRITE(DBGMSG,'(a,a)') 
     &        '### UCX: Setting long-lived species BCs from ',
     &        TRIM(WMO_FILE_NCDF)
         ENDIF
      ELSE
         ! For ASCII input, use 2x25 grid for all other grids than 4x5. 
         WRITE( WMO_ASC_FILE_ROOT,'(a,a)') 
     &        TRIM(Input_Opt%CHEM_INPUTS_DIR), 
     &      'UCX_201403/Montreal/'
         WRITE(DBGMSG,'(a,a)') 
     &     '### UCX: Setting long-lived species BCs from ',
     &     TRIM(WMO_ASC_FILE_ROOT)
      ENDIF
  
      IF ( doPrint ) THEN
         CALL DEBUG_MSG( TRIM(DBGMSG) )
      ENDIF

      ! Set "current" month to impossible value
      SFC_MONTH = 0
      OFFSETYEAR = 0
      IF (.not.LFUTURECFC) THEN
         LLSYEAR = GET_YEAR()
      ENDIF

      ! Determine which tracers are to be emitted
      N_EMIT = 0
      IDT_EMIT = 0

      IF (LCH4SBC.or.CH4Read) THEN
         ! Simple CH4 fixed surface emissions
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCH4
      ENDIF 

      IF (LCFCEMIS) THEN
         ! R-11/12/113/114/115
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCFC11

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCFC12

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCFCX

         ! HCFC-22/123/141b/142b
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTHCFC22

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTHCFCX

         ! Halons (H-1301/2402/1211)
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH1301

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH1211

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH2402
      ENDIF

      IF (LCLEMIS) THEN
         ! CCl4
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCCl4

         ! CH3Cl
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCH3Cl

         ! CH3CCl3
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCH3CCl3
      ENDIF

      IF (LOCSEMIS) THEN
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTOCS
      ENDIF

      IF (LN2OEMIS) THEN
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTN2O
      ENDIF 

      IF (LBREMIS) THEN
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCH3Br
      ENDIF

      IF (IDTH2.gt.0) THEN
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH2
      ENDIF

      IF (N_EMIT.gt.0) THEN
         ALLOCATE( GRID_EMIT( JJPAR, N_EMIT ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'GRID_EMIT' )
         GRID_EMIT = 0e+0_fp
      ENDIF

      ! Clear other tracer names for safety
      DO N=(N_EMIT+1),MAX_EMIT
         IDT_EMIT(N) = 0
      ENDDO

      ! Initialize surface boundary conditions (ckeller, 12/05/14)
      ! This reads the surface boundary conditions for all emitted species
      ! and all available years to disk and stores them in the SFCMR object.
      ! Call SFCMR_GET can then be used to pick the surface boundary condition
      ! for a given year and month.
      ! This is advantagous in an MPI environment where we try to avoid file
      ! I/O during runtime as much as possible.
      IF ( Input_Opt%LBASICEMIS .AND. .NOT. UCXNETCDF ) THEN
         CALL SFCMR_INIT ( am_I_Root, Input_Opt ) 
      ENDIF

      END SUBROUTINE INIT_LLS
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_lls
!
! !DESCRIPTION: Subroutine CLEANUP\_LLS deallocates module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_LLS ( am_I_Root )
!
! !USES:
!
! 
! !INPUT ARGUMENTS:
!
      LOGICAL, INTENT(IN)     :: am_I_Root
!
! !REVISION HISTORY: 
!  09 Apr 2015 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      INTEGER                   :: N

      !=================================================================
      ! CLEANUP_LLS begins here!
      !=================================================================

      IF ( ALLOCATED( GRID_EMIT  ) ) DEALLOCATE( GRID_EMIT  )

      ! Cleanup SFCMR object
      IF ( ASSOCIATED(SFCMR) ) THEN
         DO N=1,NSFCMR
            IF ( ALLOCATED(SFCMR(N)%YR ) ) DEALLOCATE( SFCMR(N)%YR )
            IF ( ALLOCATED(SFCMR(N)%MR ) ) DEALLOCATE( SFCMR(N)%MR )
         ENDDO
         DEALLOCATE(SFCMR)
      ENDIF

      END SUBROUTINE CLEANUP_LLS
!EOC
      END MODULE LONGLIVED_MOD
