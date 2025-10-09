!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: ucx_mod
!     
! !DESCRIPTION: Module UCX\_MOD contains routines and variables which
!  are associated with the addition of full stratospheric chemistry to
!  GEOS-Chem (based on the NASA GMI implementation, forming the Unified
!  Chemistry eXtension (UCX).
!\\   
!\\   
! !INTERFACE: 
!
      MODULE UCX_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN
      USE ERROR_MOD,  ONLY : DEBUG_MSG
      USE CMN_SIZE_MOD

!#if !defined(ESMF_)
      ! NcdfUtil modules for netCDF I/O
      USE m_netcdf_io_open                    ! netCDF open
      USE m_netcdf_io_get_dimlen              ! netCDF dimension queries
      USE m_netcdf_io_read                    ! netCDF data reads
      USE m_netcdf_io_close                   ! netCDF close
!#endif

      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
!#if !defined(ESMF_)
#     include "netcdf.inc"
!#endif

      PRIVATE

!
! !PUBLIC DATA MEMBERS:
!
      CHARACTER(LEN=255), PUBLIC :: AVG_FILE_ROOT ! Directory for diurnal avg data
      CHARACTER(LEN=255), PUBLIC :: NOON_FILE_ROOT ! Directory for noontime data
                          PUBLIC :: T_STS ! Max temperature of STS formation (K)
      REAL(fp),           PUBLIC :: T_NAT_SUPERCOOL ! NAT supercooling (K)
      REAL(fp),           PUBLIC :: P_ICE_SUPERSAT ! Ice supersaturation (-)
                          PUBLIC :: KHETI_SLA ! Reaction cofactors for SLA
                          PUBLIC :: NDENS_AER ! See below
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: SET_INITIAL_MIXRATIOS
      PUBLIC  :: SET_H2O_TRAC
      PUBLIC  :: SET_CLOCK_TRAC
      PUBLIC  :: SETTLE_STRAT_AER
      PUBLIC  :: UCX_NOX
      PUBLIC  :: UCX_H2SO4PHOT
      PUBLIC  :: CALC_STRAT_AER
      PUBLIC  :: GET_STRAT_OPT
      PUBLIC  :: KG_STRAT_AER
      PUBLIC  :: RHO_STRAT_AER
      PUBLIC  :: CHECK_NAT
      PUBLIC  :: INIT_UCX
      PUBLIC  :: DIAGINIT_UCX
      PUBLIC  :: CLEANUP_UCX
      PUBLIC  :: READ_PSC_FILE
      PUBLIC  :: WRITE_STATE_PSC
!
! PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: APPLY_2DTRAC
      PRIVATE :: TERNARY
      PRIVATE :: CARSLAW_DENSITY
      PRIVATE :: CALC_H2SO4_GAS
      PRIVATE :: CALC_SLA_GAMMA
      PRIVATE :: MOLEC_SPEED
      PRIVATE :: NOXCOEFF_INIT
      PRIVATE :: GET_JJNOX 
!
! !REVISION HISTORY: 
!  26 Mar 2013 - S. D. Eastham - Initial version
!  04 Apr 2013 - S. D. Eastham - Rolled several routines into module
!  20 Feb 2014 - M. Sulprizio  - Removed "define.h", this is now obsolete
!  24 Nov 2014 - M. Yannetti   - Added PRECISION_MOD
!  05 Dec 2014 - C. Keller     - Added SFCMR object for reading surface mixing
!                                ratios to buffer (from ASCII file).
!                                Now map NOx coefficients onto grids other than
!                                4x5 or 2x25.
!  14 Jan 2014 - C. Keller     - Now read and diagnose STATE_PSC through HEMCO.
!  08 Apr 2015 - R. Yantosca   - Now make READ_PSC_FILE a PUBLIC routine
!  08 May 2015 - C. Keller     - Added WRITE_STATE_PSC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! MODULE PARAMETERS
      !
      ! UCX_NLEVS       : Number of levels in AER data
      ! UCX_NLAT        : Number of latitudes in AER data
      ! T_STS           : Maximum temperature of STS formation (K)
      ! I_SLA           : Index of liquid aerosols
      ! I_SPA           : Index of particulate PSCs
      ! AVOGADRO        : Avogadro's number (#/mole)
      ! R_UNIV          : Universal gas constant (J/mol K)
      ! INITMR_BASIS    : Year for which the initializing mixing ratios
      !                   were calculated (needed for future-scaling)
      ! UCXNETCDF       : Read data from NetCDF
      !
      !=================================================================

      INTEGER,  PARAMETER           :: UCX_NLEVS=51
      INTEGER,  PARAMETER           :: UCX_NLAT=19
      REAL(fp), PARAMETER           :: T_STS=240.0e+0_fp
      INTEGER,  PARAMETER           :: I_SLA=1
      INTEGER,  PARAMETER           :: I_SPA=2
      REAL(fp), PARAMETER           :: AVOGADRO=6.022e+23_fp
      REAL(fp), PARAMETER           :: R_UNIV = 8.31448e+0_fp
      INTEGER,  PARAMETER           :: INITMR_BASIS = 2005

#if defined(ESMF_)
      ! Never use NETCDF in ESMF environment (ckeller, 12/05/14). 
      LOGICAL, PARAMETER            :: UCXNETCDF = .FALSE.
#else
      LOGICAL, PARAMETER            :: UCXNETCDF = .TRUE.
#endif

!
! PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! Scalars
      !
      ! TRAC_IDX           : Tracer index for output
      ! N_EMIT             : Number of tracers to emit
      ! MAX_EMIT           : Limit on tracer emission count
      ! SFC_MONTH          : Current month
      ! OFFSETYEAR         : Number of year-ends since simulation start
      ! MONTH_STR_INT      : Current month
      ! SLA_VA             : SLA volume-area conversion
      ! SLA_RR             : SLA effective-liquid radius conversion
      ! SLA_VR             : SLA volume-effective radius conversion
      ! NATMW              : Molar mass of NAT (kg/kmol)
      ! ICEMEW             : Molar mass of ice (kg/kmol)
      ! DENSNAT            : Density of pure NAT (kg/m3)
      ! DENSICE            : Density of pure ice (kg/m3)
      ! ISR_ClNO3          : ClNO3 MW (inverse sqrt) (kg/kmol)^-0.5 
      ! ISR_BrNO3          : BrNO3 MW (inverse sqrt) (kg/kmol)^-0.5 
      ! ISR_N2O5           : N2O5 MW (inverse sqrt) (kg/kmol)^-0.5
      ! ISR_HOCl           : HOCl MW (inverse sqrt) (kg/kmol)^-0.5
      ! ISR_HOBr           : HOBr MW (inverse sqrt) (kg/kmol)^-0.5
      ! CLOCKMR            : Current base MR for the clock tracer
      ! NSTSWACCM          : LBS/STS number density (WACCM) (#/m3)
      ! NNATWACCM          : NAT number density (WACCM) (#/m3)
      ! SIGSTSWACCM        : LBS/STS size std dev (WACCM) (-)
      ! SIGNATWACCM        : NAT size std dev (WACCM) (-)
      !
      ! Arrays
      !
      ! UCX_MR             : AER 2D mixing ratios (v/v)
      ! UCX_PLEVS          : Pressure levels of 2D data (hPa)
      ! UCX_LATS           : Latitude edges of 2D data (deg)
      ! IDT_EMIT           : IDs of emitted tracers
      ! GRID_EMIT          : Surface CFC mixing ratios
      ! SFC_MONTH          : Month of last read-in
      ! RAD_AER            : Strat. aerosol radius (cm)
      ! KG_AER             : Aerosol mass (kg/box)
      ! SAD_AER            : Aerosol surface area density (cm2/cm3)
      ! NDENS_AER          : Aerosol number density (#/m3)
      ! RHO_AER            : Aerosol mass density (kg/m3 aerosol)
      ! STATE_PSC          : Gridbox PSC type (see Kirner et al)
      ! KHETI_SLA          : Liquid aerosol reaction cofactors
      ! AERFRAC            : Mass fraction of species in liquid aerosols
      ! AERFRACIND         : Indices of liquid aerosol species
      ! NOX_O              : Monthly mean noontime O3P/O1D for NOx calcs
      ! NOX_J              : Monthly mean noontime J-rates for NOx calcs
      ! SO4_TOPPHOT        : Photolysis rate at the top of the chemgrid (1/s)
      !
      ! Strings
      !
      ! TRAC_FILE          : Input filename
      ! MONTH_STR          : Current month
      ! TRAC_SET           : Tracer name to be overwritten
      ! TRAC_EMIT          ! Tracer names for emission
      !
      ! Logicals
      ! 
      ! TRAC_ADD            : Add to (rather than overwrite) tracer
      ! STRAT2DCH4          : Use 2D estimates of stratospheric CH4?
      ! 
      !=================================================================

      ! Scalars
      INTEGER              :: TRAC_IDX
      INTEGER              :: SFC_MONTH
      INTEGER              :: OFFSETYEAR
      INTEGER              :: MONTH_STR_INT
      REAL(fp)             :: SLA_VA
      REAL(fp)             :: SLA_RR
      REAL(fp)             :: SLA_VR
      REAL(fp), PARAMETER  :: NATMW   = 117.0
      REAL(fp), PARAMETER  :: ICEMW   = 18.0
      REAL(fp), PARAMETER  :: DENSNAT = 1626.e+0_fp
      REAL(fp), PARAMETER  :: DENSICE = 990.0e+0_fp
      REAL(fp), PARAMETER  :: ISR_ClNO3=1.e+0_fp/sqrt(97.46e+0_fp)
      REAL(fp), PARAMETER  :: ISR_BrNO3=1.e+0_fp/sqrt(141.9e+0_fp)
      REAL(fp), PARAMETER  :: ISR_N2O5 =1.e+0_fp/sqrt(108.0e+0_fp)
      REAL(fp), PARAMETER  :: ISR_HOCl =1.e+0_fp/sqrt(52.46e+0_fp)
      REAL(fp), PARAMETER  :: ISR_HOBr =1.e+0_fp/sqrt(96.91e+0_fp)
      REAL(fp)             :: CLOCKMR
      REAL(fp)             :: NWACCM(2)
      REAL(fp)             :: LSSQWACCM(2)
      REAL(fp)             :: LN_RMED_TO_REFF(2)
      REAL(fp)             :: LN_RMED_TO_RLIQ(2)
      REAL(fp)             :: LN_ACONV(2)
      REAL(fp)             :: LN_VCONV(2)
      REAL(fp)             :: CHENA1, CHENA2
                                            
      ! Arrays                              
      REAL(fp),DIMENSION(:,:),ALLOCATABLE        :: UCX_MR_IN
      REAL(fp),DIMENSION(:,:),ALLOCATABLE        :: UCX_MR
      REAL(fp),DIMENSION(:,:),ALLOCATABLE        :: UCX_REGRID
      REAL(fp),DIMENSION(:,:),ALLOCATABLE        :: UCX_REGRID_225
      REAL(fp),DIMENSION(:),ALLOCATABLE          :: UCX_PLEVS
      REAL(fp),DIMENSION(:),ALLOCATABLE       :: UCX_LATS
      REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: RAD_AER
      REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: KG_AER
      REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: SAD_AER
      REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: NDENS_AER
      REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: RHO_AER
      !INTEGER,DIMENSION(:,:,:),ALLOCATABLE  :: STATE_PSC
      REAL(f4),DIMENSION(:,:,:),POINTER       :: STATE_PSC => NULL()
      REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: KHETI_SLA
      REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: AERFRAC
      INTEGER,DIMENSION(:),ALLOCATABLE      :: AERFRACIND
      REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: NOX_O
      REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: NOX_J
      REAL(fp),DIMENSION(:,:),ALLOCATABLE     :: SO4_TOPPHOT
 
      ! Strings
      CHARACTER(LEN=255)                :: TRAC_FILE
      CHARACTER(LEN=2)                  :: MONTH_STR
      CHARACTER(LEN=255)                :: TRAC_SET 
      CHARACTER*20,DIMENSION(50)        :: TRAC_EMIT

      ! Logicals
      LOGICAL                           :: doPrint
      LOGICAL                           :: TRAC_ADD
      LOGICAL                           :: USE2DDATA
#if defined(ESMF_)
      LOGICAL,PARAMETER                 :: STRAT2DCH4=.FALSE.
#else
      LOGICAL,PARAMETER                 :: STRAT2DCH4=.TRUE.
#endif

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

      ! Maximum number of elements fo SFCMR
      INTEGER,        PARAMETER :: MAXSFCMR = 20

      ! Actual number of used elements in SFCMR
      INTEGER                   :: NSFCMR   = 0

      !=================================================================
      ! Variables to use NOx coefficients in ESMF / grid-independent
      ! envionment. The NOx coefficients are climatological 2D 
      ! (lat/lev/12 months) data that are currently available for 
      ! horizontal (latitude) resolutions of 2 and 4 degrees. For other
      ! resolutions, the horizontal data becomes mapped onto the
      ! simulation grid (see GET_JJNOX).
      ! Similar to the surface mixing ratio boundary conditions, we now
      ! read all the NOx coefficients during initialization to avoid
      ! additional I/O calls during run time (only if UCXNETCDF = false)
      ! (ckeller, 05/12/2014).
      !=================================================================
      REAL(fp), ALLOCATABLE, TARGET :: NOXCOEFF(:,:,:,:)
      REAL(fp), ALLOCATABLE         :: NOXLAT(:)
      INTEGER                     :: JJNOXCOEFF

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
! !IROUTINE: set_initial_mixratios
!
! !DESCRIPTION: Subroutine SET\_INITIAL\_MIXRATIOS is a public interface.
!  Concentrations of tracers are read from 2D estimates made using the
!  AER 2D model.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_INITIAL_MIXRATIOS( am_I_Root, Input_Opt, State_Met,
     &                                  State_Chm )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : ALLOC_ERR,ERROR_STOP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_YMID
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : GET_YEAR
      USE LONGLIVED_MOD,      ONLY : LLS_SCALEFAC
      USE LONGLIVED_MOD,      ONLY : GET_SCALED_CH3Br
      USE LONGLIVED_MOD,      ONLY : GET_SCALED_CH4
      USE TRACERID_MOD
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(OUT) :: State_Chm   ! Chemistry State object
!
! !REMARKS:
!  Concentrations of each of the following tracers are read from
!  2D estimates made using the AER 2D model:
!                                                                             .
!     (# ) TRC                 = Constituents         < Controlled by
! ----------------------------------------------------------------------------
!     (1 ) CH4                 = CH4                  < LSETCH4
!     (2 ) N2O                 = N2O                  < LSETN2O
!     (3 ) OCS                 = OCS                  < LSETOCS
!     (4 ) SO4                 = H2SO4                < LSETH2SO4
!     (5 ) CFC                 = CFC-113/114/115      < LSETCFC
!     (6 ) HCFC                = HCFC-22/141b/142b    < LSETCFC
!     (7 ) CFC11               = CFC-11               < LSETCFC
!     (8 ) CFC12               = CFC-12               < LSETCFC
!     (9 ) H1202               = Halon 1202           < LSETCFC
!     (10) H1211               = Halon 1211           < LSETCFC
!     (11) H1301               = Halon 1301           < LSETCFC
!     (12) H2402               = Halon 2402           < LSETCFC
!     (13) Cl2                 = Cl2                  < LSETCL
!     (14) ClOx                = Cl + ClO             < LSETCL
!     (15) CCl4                = CCl4                 < LSETCL
!     (16) CH3Cl               = CH3Cl                < LSETCL
!     (17) CH3CCl3             = CH3CCl3              < LSETCL
!     (18) HCl                 = HCl                  < LSETCL
!     (19) HOCl                = HOCl                 < LSETCL
!     (20) Cl2O2               = Cl2O2                < LSETCL
!     (21) ClNO2               = ClNO2                < LSETCL
!     (22) ClONO2              = ClNO3                < LSETCL
!     (23) OClO                = OClO                 < LSETCL
!     (24) ClOO                = ClOO                 < LSETCL
!     (25) BrCl                = BrCl                 < (LSETCL || LSETBR || LSETBRSTRAT)
!     (26) Br2                 = Br2                  < LSETBR || LSETBRSTRAT
!     (27) Br                  = Br                   < LSETBR || LSETBRSTRAT 
!     (28) BrO                 = BrO                  < LSETBR || LSETBRSTRAT
!     (29) HOBr                = HOBr                 < LSETBR || LSETBRSTRAT
!     (30) HBr                 = HBr                  < LSETBR || LSETBRSTRAT
!     (31) BrNO2               = BrNO2                < LSETBR || LSETBRSTRAT
!     (32) BrNO3               = BrNO3                < LSETBR || LSETBRSTRAT
!     (33) CHBr3               = CHBr3                < LSETBR || LSETBRSTRAT
!     (34) CH2Br2              = CH2Br2               < LSETBR || LSETBRSTRAT
!     (35) CH3Br               = CH3Br                < LSETBR || LSETBRSTRAT
!
! !REVISION HISTORY: 
!  26 Mar 2013 - S. D. Eastham - Initial version
!  07 Feb 2014 - R. Yantosca   - !$OMP DO loops now go in L-J-I order
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm 
!                                objects via the arg list
!  18 Jun 2014 - S. D. Eastham - Now read data for individual NOx species
!  24 Nov 2014 - C. Keller     - Update for stratospheric NOy (as suggested by 
!                                Seb Eastham)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: AS, N, I, J, L
      REAL(fp)           :: PLAST, PCURR, PNEXT
      REAL(fp)           :: YLAT
      CHARACTER(LEN=255) :: MSG
      REAL(fp)           :: C3090S, C0030S, C0030N, C3090N
      REAL(fp)           :: C5590S, C0055S, C0055N, C5590N
      INTEGER            :: CH4_YEAR, INIT_MO
      REAL(fp)           :: SCALEFAC

      ! Local variables for quantities from Input_Opt
      LOGICAL            :: LPRT
      LOGICAL            :: LSETCH4
      LOGICAL            :: LSETOCS
      LOGICAL            :: LSETCFC
      LOGICAL            :: LSETCL
      LOGICAL            :: LSETBR
      LOGICAL            :: LSETN2O
      LOGICAL            :: LSETH2SO4
      LOGICAL            :: LSETBRSTRAT
      LOGICAL            :: LSETNOYSTRAT
      INTEGER            :: N_TRACERS
      REAL(fp)           :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      REAL(fp), POINTER    :: STT(:,:,:,:)

      !=================================================================
      ! SET_INITIAL_MIXRATIOS begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      LPRT         = ( Input_Opt%LPRT .and. am_I_Root )
      LSETCH4      = Input_Opt%LSETCH4
      LSETOCS      = Input_Opt%LSETOCS
      LSETCFC      = Input_Opt%LSETCFC
      LSETCL       = Input_Opt%LSETCL
      LSETBR       = Input_Opt%LSETBR
      LSETN2O      = Input_Opt%LSETN2O
      LSETH2SO4    = Input_Opt%LSETH2SO4
      LSETBRSTRAT  = Input_Opt%LSETBRSTRAT
      LSETNOYSTRAT = Input_Opt%LSETNOYSTRAT
      N_TRACERS    = Input_Opt%N_TRACERS
      TCVV         = Input_Opt%TCVV(1:N_TRACERS)

      ! Set module variable

      ! Initialize GEOS-Chem tracer array [kg]
      STT => State_Chm%Tracers

      IF ( am_I_Root ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100   )
         WRITE( 6, 110   )
         WRITE( 6, 120   )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      ENDIF
   
      ! FORMAT strings
 100  FORMAT( 'T R A C E R   I N I T I A L I Z A T I O N'   )
 110  FORMAT( 'Routine written by SEBASTIAN D. EASTHAM' )
 120  FORMAT( 'Last modified: 03/26/13' )

      ! Get current month string
      MONTH_STR_INT = GET_MONTH()
      WRITE( MONTH_STR, '(i0.2)' ) MONTH_STR_INT

      IF (USE2DDATA) THEN

         ! Calculate future scaling factor
         INIT_MO = GET_MONTH()

         ALLOCATE( UCX_MR_IN( UCX_NLAT, UCX_NLEVS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'UCX_MR_IN' )
         UCX_MR_IN = 0e+0_fp

         ALLOCATE( UCX_MR( JJPAR, UCX_NLEVS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'UCX_MR' )
         UCX_MR = 0e+0_fp

         IF ( am_I_Root ) THEN
            WRITE(6,*) 'Setting initial values for:'
            WRITE(6,*) '|---------|---------|---------|'
         ENDIF

         ! Bromine species (except Halons)
         IF (LSETBR.or.LSETBRSTRAT) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|Br2      |Br       |BrO      |'
               WRITE(6,*) '|HOBr     |HBr      |BrNO2    |'
               WRITE(6,*) '|BrONO2   |CHBr3    |CH2Br2   |'
               WRITE(6,*) '|CH3Br    |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF

            ! Scaling data only available for CH3Br
            SCALEFAC = 1.e+0_fp

            TRAC_IDX = IDTBr2
            TRAC_SET = 'BR2'
            CALL APPLY_2DTRAC( .FALSE., LSETBRSTRAT, .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
         
            TRAC_IDX = IDTBr
            TRAC_SET = 'BR'
            CALL APPLY_2DTRAC( .FALSE., LSETBRSTRAT, .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
         
            TRAC_IDX = IDTBrO
            TRAC_SET = 'BRO'
            CALL APPLY_2DTRAC( .FALSE., LSETBRSTRAT, .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
         
            TRAC_IDX = IDTHOBr
            TRAC_SET = 'HOBR'
            CALL APPLY_2DTRAC( .FALSE., LSETBRSTRAT, .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
         
            TRAC_IDX = IDTHBr
            TRAC_SET = 'HBR'
            CALL APPLY_2DTRAC( .FALSE., LSETBRSTRAT, .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
     
            TRAC_IDX = IDTBrNO2
            TRAC_SET = 'BRNO2'
            CALL APPLY_2DTRAC( .FALSE., LSETBRSTRAT, .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
     
            TRAC_IDX = IDTBrNO3
            TRAC_SET = 'BRNO3'
            CALL APPLY_2DTRAC( .FALSE., LSETBRSTRAT, .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
     
            TRAC_IDX = IDTCHBr3 
            TRAC_SET = 'CHBR3'
            CALL APPLY_2DTRAC( .FALSE., LSETBRSTRAT, .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
     
            TRAC_IDX = IDTCH2Br2
            TRAC_SET = 'CH2BR2'
            CALL APPLY_2DTRAC( .FALSE., LSETBRSTRAT, .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
     
            ! Change this to match CH4-style input
            ! Retrieve 4 large-scale averages, scaled to target year
            CALL GET_SCALED_CH3Br(C5590S, C0055S, C0055N, C5590N,
     &                        am_I_Root, Input_Opt )

            ! Convert from [vv] to [kg CH3Br/mol air]
            C5590S = C5590S / TCVV(IDTCH3Br)
            C0055S = C0055S / TCVV(IDTCH3Br)
            C0055N = C0055N / TCVV(IDTCH3Br)
            C5590N = C5590N / TCVV(IDTCH3Br)

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, YLAT )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Grid-box latitude index
               YLAT  = GET_YMID( I,J,L )

               ! Set CH3Br according to latitude
               IF ( YLAT < -55e+0_fp ) THEN
                  STT(I,J,L,IDTCH3Br) = C5590S * State_Met%AD(I,J,L)
               ELSE IF ( YLAT >= -55e+0_fp .and. YLAT < 0e0_fp  ) THEN
                  STT(I,J,L,IDTCH3Br) = C0055S * State_Met%AD(I,J,L)
               ELSE IF ( YLAT >= 0e+0_fp .and. YLAT < 55e+0_fp ) THEN
                  STT(I,J,L,IDTCH3Br) = C0055N * State_Met%AD(I,J,L)
               ELSE
                  STT(I,J,L,IDTCH3Br) = C5590N * State_Met%AD(I,J,L)
               ENDIF
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

!            TRAC_IDX = IDTCH3Br
!            TRAC_SET  = 'CH3BR'
!            SCALEFAC = LLS_SCALEFAC('CH3Br',INITMR_BASIS,INIT_MO,
!     &                                     Input_Opt,am_I_Root)
!            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.,SCALEFAC,
!     &                         Input_Opt, State_Met, State_Chm )
         ENDIF

         ! Stratospheric NOy (update ckeller, 11/25/14)
         IF (LSETNOYSTRAT) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|NO       |NO2      |NO3      |'
               WRITE(6,*) '|---------|---------|---------|'
               WRITE(6,*) '|HNO2     |HNO3     |HNO4     |'
               WRITE(6,*) '|---------|---------|---------|'
               WRITE(6,*) '|N2O5     |PAN      |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
 
            ! No Montreal scaling
            SCALEFAC = 1.e+0_fp

            TRAC_IDX = IDTNO
            TRAC_SET = 'NO'
            CALL APPLY_2DTRAC( .FALSE., .TRUE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            TRAC_IDX = IDTNO2
            TRAC_SET = 'NO2'
            CALL APPLY_2DTRAC( .FALSE., .TRUE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
 
            ! Unstable product of ETO2 + NO2 gets lumped into NO2 for now
            TRAC_SET = 'ETHO2NO2'
            CALL APPLY_2DTRAC( .TRUE., .TRUE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            TRAC_IDX = IDTNO3
            TRAC_SET = 'NO3'
            CALL APPLY_2DTRAC( .FALSE., .TRUE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            TRAC_IDX = IDTN2O5
            TRAC_SET = 'N2O5'
            CALL APPLY_2DTRAC( .FALSE., .TRUE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            TRAC_IDX = IDTHNO2
            TRAC_SET = 'HONO'
            CALL APPLY_2DTRAC( .FALSE., .TRUE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            TRAC_IDX = IDTHNO3
            TRAC_SET = 'HNO3'
            CALL APPLY_2DTRAC( .FALSE., .TRUE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            TRAC_IDX = IDTHNO4
            TRAC_SET = 'HO2NO2'
            CALL APPLY_2DTRAC( .FALSE., .TRUE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            TRAC_IDX = IDTPAN
            TRAC_SET = 'PAN'
            CALL APPLY_2DTRAC( .FALSE., .TRUE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
         ENDIF

         ! Br-Cl (special case)
         IF (LSETCL.or.LSETBR.or.LSETBRSTRAT) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|BrCl     |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF

            ! No Montreal scaling
            SCALEFAC = 1.e+0_fp

            TRAC_IDX = IDTBrCl
            TRAC_SET  = 'BRCL'
            CALL APPLY_2DTRAC( .FALSE., LSETBRSTRAT, .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
     
         ENDIF

         ! Chlorine species (except CFCS, HCFCs and Halons)
         IF (LSETCL) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|Cl2      |ClO      |Cl2O2    |'
               WRITE(6,*) '|Cl       |ClONO2   |HCl      |'
               WRITE(6,*) '|HOCl     |OClO     |CH3Cl    |'
               WRITE(6,*) '|CH3CCl3  |CCl4     |ClOO     |'
               WRITE(6,*) '|         |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF

            ! Scale to CFC-11
            SCALEFAC = LLS_SCALEFAC('CFC11',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)

            TRAC_IDX = IDTCl
            TRAC_SET = 'CL'
            CALL APPLY_2DTRAC( .FALSE., .FALSE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
            TRAC_IDX = IDTClO
            TRAC_SET = 'CLO'
            CALL APPLY_2DTRAC( .FALSE., .FALSE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
            TRAC_IDX = IDTClNO2
            TRAC_SET = 'CLNO2'
            CALL APPLY_2DTRAC( .FALSE., .FALSE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
            TRAC_IDX = IDTClNO3
            TRAC_SET = 'CLNO3'
            CALL APPLY_2DTRAC( .FALSE., .FALSE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
            TRAC_IDX = IDTHOCl
            TRAC_SET = 'HOCL'
            CALL APPLY_2DTRAC( .FALSE., .FALSE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
            TRAC_IDX = IDTOClO
            TRAC_SET = 'OCLO'
            CALL APPLY_2DTRAC( .FALSE., .FALSE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
            TRAC_IDX = IDTClOO
            TRAC_SET = 'CLOO'
            CALL APPLY_2DTRAC( .FALSE., .FALSE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
     
            TRAC_IDX = IDTCl2
            TRAC_SET = 'CL2'
            CALL APPLY_2DTRAC( .FALSE., .FALSE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
            TRAC_IDX = IDTCl2O2
            TRAC_SET = 'CL2O2'
            CALL APPLY_2DTRAC( .FALSE., .FALSE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
     
            TRAC_IDX = IDTHCl
            TRAC_SET = 'HCL'
            CALL APPLY_2DTRAC( .FALSE., .FALSE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            ! Handle chloroalkanes      
            TRAC_IDX = IDTCH3Cl
            TRAC_SET = 'CH3CL'
            SCALEFAC = LLS_SCALEFAC('CH3Cl',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
     
            TRAC_IDX = IDTCH3CCl3
            TRAC_SET = 'CH3CCL3'
            SCALEFAC = LLS_SCALEFAC('CH3CCl3',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
     
            TRAC_IDX = IDTCCl4
            TRAC_SET = 'CCL4'
            SCALEFAC = LLS_SCALEFAC('CCl4',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
         ENDIf
         
         ! CFCs, HCFCs and Halons
         IF (LSETCFC) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|H-1211   |H-1301   |H-2402   |'
               WRITE(6,*) '|HCFC-123 |HCFC-141b|HCFC-142b|'
               WRITE(6,*) '|HCFC-22  |R-11     |R-12     |'
               WRITE(6,*) '|R-113    |R-114    |R-115    |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            ! H-1211, 1301 and 2402
            TRAC_IDX = IDTH1211
            TRAC_SET = 'CBRCLF2'
            SCALEFAC = LLS_SCALEFAC('H1211',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
            TRAC_IDX = IDTH1301
            TRAC_SET = 'CBRF3'
            SCALEFAC = LLS_SCALEFAC('H1301',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
            TRAC_IDX = IDTH2402
            TRAC_SET = 'C2BR2F4'
            SCALEFAC = LLS_SCALEFAC('C2BR2F4',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
     
            ! CFC  11, 12, and [113/114/115]
            TRAC_IDX = IDTCFC11
            TRAC_SET = 'CFCL3'
            SCALEFAC = LLS_SCALEFAC('CFC11',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            TRAC_IDX = IDTCFC12
            TRAC_SET = 'CF2CL2'
            SCALEFAC = LLS_SCALEFAC('CFC12',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            TRAC_IDX = IDTCFCX
            TRAC_SET = 'C2CL3F3'
            SCALEFAC = LLS_SCALEFAC('CFC113',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            TRAC_SET = 'C2CL2F4'
            SCALEFAC = LLS_SCALEFAC('CFC114',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
            TRAC_SET = 'C2CLF5'
            SCALEFAC = LLS_SCALEFAC('CFC115',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
     
            ! HCFC 22, 123, 141b and 142b
            TRAC_IDX = IDTHCFC22
            TRAC_SET = 'CHCLF2'
            SCALEFAC = LLS_SCALEFAC('HCFC22',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            TRAC_IDX = IDTHCFCX
            TRAC_SET = 'C2HF3CL2'
            ! Take scaling factor for HCFC22 (no data for 123)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
            TRAC_SET = 'C2H3FCL2'
            SCALEFAC = LLS_SCALEFAC('HCFC141b',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            TRAC_SET = 'C2H3F2CL'
            SCALEFAC = LLS_SCALEFAC('HCFC142b',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.,.FALSE.,SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
         ENDIF

         IF (LSETOCS) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|OCS      |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            TRAC_IDX = IDTOCS
            TRAC_SET = 'OCS'

            ! No Montreal scaling
            SCALEFAC = 1.e+0_fp

            CALL APPLY_2DTRAC( .FALSE., .FALSE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
         ENDIF

         IF (LSETH2SO4) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|H2SO4    |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            ! No Montreal scaling
            SCALEFAC = 1.e+0_fp

            TRAC_IDX = IDTSO4
            ! Actually H2SO4 aerosol
            TRAC_SET = 'SO4'
            CALL APPLY_2DTRAC( .FALSE., .TRUE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )

            ! Gaseous H2SO4
            IF (IDTSO4g > 0) THEN
               TRAC_IDX = IDTSO4g
               TRAC_SET = 'SO4g'
               ! Zero out the troposphere
               STT(:,:,:,TRAC_IDX) = 0d0
               CALL APPLY_2DTRAC( .FALSE., .TRUE., .FALSE., SCALEFAC,
     &                            Input_Opt, State_Met, State_Chm )
            ENDIF
         ENDIF
      ENDIF ! USE2DDATA

      IF (LSETN2O) THEN
         IF ( am_I_Root ) THEN
            WRITE(6,*) '|N2O      |         |         |'
            WRITE(6,*) '|---------|---------|---------|'
            IF ( LPRT ) CALL DEBUG_MSG( '### UCX: Forcing 2D N2O' )
         ENDIF

         TRAC_IDX = IDTN2O
         TRAC_SET = 'N2O'
         SCALEFAC = LLS_SCALEFAC('N2O',INITMR_BASIS,INIT_MO,
     &                Input_Opt, am_I_Root)
         CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC,
     &                      Input_Opt, State_Met, State_Chm )
      ENDIF

      ! CH4 taken from the existing GEOS-Chem input data,
      ! rather than from 2D data

      IF (LSETCH4) THEN
         IF ( am_I_Root ) THEN
            WRITE(6,*) '|CH4      |         |         |'
            WRITE(6,*) '|---------|---------|---------|'
         ENDIF

         ! Retrieve 4 large-scale averages, scaled to target year
         CALL GET_SCALED_CH4( C3090S, C0030S, C0030N, C3090N,
     &                     am_I_Root, Input_Opt )

         ! Convert from [vv] to [kg CH4/mol air]
         C3090S = C3090S / TCVV(IDTCH4)
         C0030S = C0030S / TCVV(IDTCH4)
         C0030N = C0030N / TCVV(IDTCH4)
         C3090N = C3090N / TCVV(IDTCH4)

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, YLAT )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Grid-box latitude index
            YLAT  = GET_YMID( I,J,L )

            ! Set CH4 according to latitude
            IF ( YLAT < -30e+0_fp ) THEN
               STT(I,J,L,IDTCH4) = C3090S * State_Met%AD(I,J,L)
            ELSE IF ( YLAT >= -30e+0_fp .and. YLAT < 0e+0_fp  ) THEN
               STT(I,J,L,IDTCH4) = C0030S * State_Met%AD(I,J,L)
            ELSE IF ( YLAT >=   0e+0_fp .and. YLAT < 30e+0_fp ) THEN
               STT(I,J,L,IDTCH4) = C0030N * State_Met%AD(I,J,L)
            ELSE
               STT(I,J,L,IDTCH4) = C3090N * State_Met%AD(I,J,L)
            ENDIF
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Set stratospheric CH4 from AER data?
         IF ( STRAT2DCH4 ) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) 'Setting stratospheric CH4 from 2D data.'
            ENDIF
            ! No Montreal scaling (already applied)
            SCALEFAC = 1.e+0_fp

            ! Code to set strat using AER 2D data
            TRAC_IDX = IDTCH4
            TRAC_SET = 'CH4'
            CALL APPLY_2DTRAC( .FALSE., .TRUE., .FALSE., SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
         ENDIF
      ENDIF

      ! Free pointer
      NULLIFY( STT )
 
      END SUBROUTINE SET_INITIAL_MIXRATIOS
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ucx_nox
!
! !DESCRIPTION: Subroutine UCX\_NOX calculates NOx and N2O loss  
!  rates above the chemistry grid, based on estimates of j-rates from 
!  a 2D model and simple photochemical assumptiones.              
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE UCX_NOX( Input_Opt, State_Met, State_Chm )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_CHEMGRID
      USE CMN_GCTM_MOD,       ONLY : PI, RE
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_YMID
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TIME_MOD,           ONLY : GET_DAY_OF_YEAR
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : ITS_A_LEAPYEAR
      USE TIME_MOD,           ONLY : GET_HOUR
      USE TIME_MOD,           ONLY : GET_LOCALTIME
      USE TIME_MOD,           ONLY : GET_MINUTE
      USE TRACERID_MOD,       ONLY : IDTNO, IDTN2O, IDTO3
      USE TRACERID_MOD,       ONLY : IDTNO2,IDTNO3, IDN
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY: 
!  16 Jul 2013 - S. D. Eastham - Initial version
!  14 Feb 2014 - R. Yantosca   - Now pull computation of ZMID out of main
!                                parallel DO loop.  Switch loop order to L-J-I.
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm
!                                objects via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! For indexing the O1D/O3P concentration array
      INTEGER, PARAMETER :: O3PIDX = 1
      INTEGER, PARAMETER :: O1DIDX = 2

      ! For indexing the J-rate coefficient array
      INTEGER, PARAMETER :: JNOIDX = 1
      INTEGER, PARAMETER :: JNO2IDX= 2
      INTEGER, PARAMETER :: JNO3IDX= 3
      INTEGER, PARAMETER :: JN2OIDX= 4

      ! Reaction rates and indices
      REAL(fp),DIMENSION(12):: RRATE
      INTEGER, PARAMETER :: k_JNO = 6
      INTEGER, PARAMETER :: k_JNO2 = 4
      INTEGER, PARAMETER :: k_JNO3 = 5
      INTEGER, PARAMETER :: k_JN2O = 12

      ! Intermediate variables
      REAL(fp) :: LOCALNOX, LOCALN2O, LOCALO3, LOCALO1D, LOCALO3P
      REAL(fp) :: FRACNO, FRACNO2, FRACNO3, FRACN
      REAL(fp) :: NOXRATE, N2ORATE, KGNOX, KGN2O
      REAL(fp) :: NO_ALPHA, NO_BETA, NO_GAMMA, NO_EPSILON
      REAL(fp) :: NO_QA, NO_QX, NO_QC, DNOX, DN2O
      REAL(fp) :: MESONOX_DELTA
      REAL(fp) :: MESON2O_DELTA

      ! Local air number density (molec/cm3) and box mass (kg)
      REAL(fp) :: NDAIR
      REAL(fp) :: XAIR

      ! Local temperature (K) and inverted T (1/K)
      REAL(fp) :: T3K, TINV

      ! Chemistry timestep in seconds
      REAL(f8) :: DTCHEM

      ! Timestep in hours
      REAL(f8) :: DTCHEM_HR

      ! Timing information
      INTEGER       :: NDAY
      INTEGER, SAVE :: LASTMONTH=0

      ! Grid indexing
      INTEGER :: I,J,L
      REAL(fp)  :: MIDLAT, ZDEL, ZBASE
      REAL(fp)  :: ZMID(IIPAR,JJPAR,LLPAR)

      ! Local daylight fraction
      REAL(fp)             :: DAYFRAC
      REAL(fp)             :: SOL_ME, SOL_FE, SOL_LS, SOL_DEC
      REAL(fp)             :: SOL_ECC2, SOL_B, SOL_C, SOL_TEST
      REAL(fp)             :: TSTART, TEND, TMIN, TMAX
      REAL(f8)             :: GMTSTART
      REAL(fp), PARAMETER  :: SOL_ECC=0.016
      INTEGER            :: NDYEAR, AS
      LOGICAL            :: CYCLEBOX
      CHARACTER(LEN=255) :: DBGMSG

      ! Local variables for quantities from Input_Opt
      LOGICAL            :: LPRT
      INTEGER            :: N_TRACERS
      REAL(fp)             :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      REAL(fp), POINTER    :: STT(:,:,:,:)
      REAL(fp), POINTER    :: Spc(:,:,:,:)

      ! Updated NOx chemistry
      Real(fp) :: OLD_N, OLD_NO, OLD_NO2, OLD_NO3, OLD_N2O
      Real(fp) :: NEW_N, NEW_NO, NEW_NO2, NEW_NO3, NEW_N2O

      !=================================================================
      ! UCX_NOX begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      LPRT      = Input_Opt%LPRT
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg]
      STT => State_Chm%Tracers

      ! Also get the GEOS-Chem species array [molec/cm3]
      Spc => State_Chm%Species

      ! Retrieve monthly mean data if necessary
      IF (LASTMONTH.ne.GET_MONTH()) THEN
         LASTMONTH = GET_MONTH()
         CALL GET_NOXCOEFF( LASTMONTH, Input_Opt, State_Met )
      ENDIF

      ! Get chemistry step length in seconds
      DTCHEM = GET_TS_CHEM()*60.e+0_f8
      DTCHEM_HR = DTCHEM/3600.e+0_f8

      IF (ITS_A_LEAPYEAR()) THEN
         NDYEAR = 366
      ELSE
         NDYEAR = 365
      ENDIF

      NDAY = GET_DAY_OF_YEAR()

      ! Calculate timestep start and end (GMT)
      GMTSTART = REAL(GET_HOUR()) + (REAL(GET_MINUTE())/60.e+0_f8)
      GMTSTART = GMTSTART - (DTCHEM_HR/2.e+0_f8)

      ! Eccentricity^2
      SOL_ECC2 = SOL_ECC*SOL_ECC

      ! Calculate solar declination
      SOL_ME = 2.e+0_fp*PI*(REAL(NDAY)/REAL(NDYEAR)-0.009589)
      SOL_FE = SOL_ME + (2.*SOL_ECC*SIN(SOL_ME)) + 
     &         (1.25e+0_fp*SOL_ECC2*SIN(2.e+0_fp*SOL_ME))
      SOL_LS = SOL_FE - (77.240e+0_fp*PI/180.e+0_fp)
      SOL_DEC = ASIN(SIN(23.44*PI/180.e+0_fp)*SIN(SOL_LS))

      ! Reset NOx/N2O mass counters
      MESONOX_DELTA = 0e+0_fp
      MESON2O_DELTA = 0e+0_fp

      ! First compute ZMID outside of main parallel loop
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, ZBASE, ZDEL )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Zero base height
         ZBASE = 0e+0_fp

         ! Compute the array of midpoint heights [m]
         DO L = 1, LLPAR
            ZDEL        = State_Met%BXHEIGHT(I,J,L)
            ZMID(I,J,L) = ZBASE + (ZDEL/2.e+0_fp)
            ZBASE       = ZBASE + ZDEL
         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Main parallel DO loop over lon, lat, alt
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,        J,        MIDLAT,   SOL_B         ) 
!$OMP+PRIVATE( SOL_C,    TSTART,   TEND                    )
!$OMP+PRIVATE( DAYFRAC,  L,        DN2O,     CYCLEBOX      )
!$OMP+PRIVATE( TMAX,     RRATE,    T3K,      TINV          )
!$OMP+PRIVATE( NDAIR,    XAIR,     KGNOX,    KGN2O         )
!$OMP+PRIVATE( LOCALNOX, LOCALN2O, LOCALO3,  LOCALO3P      )
!$OMP+PRIVATE( LOCALO1D, NO_ALPHA, NO_BETA,  NO_GAMMA      )
!$OMP+PRIVATE( FRACNO2,  FRACNO3,  FRACNO,   FRACN         )
!$OMP+PRIVATE( NOXRATE,  N2ORATE,  DNOX                    )
!$OMP+PRIVATE( OLD_N, OLD_NO, OLD_NO2, OLD_NO3, OLD_N2O    )
!$OMP+PRIVATE( NEW_N, NEW_NO, NEW_NO2, NEW_NO3, NEW_N2O    )
!$OMP+PRIVATE( SOL_TEST, TMIN                              )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Are we in the chemistry grid? 
         CYCLEBOX = ITS_IN_THE_CHEMGRID(I,J,L,State_Met)
         IF (CYCLEBOX) CYCLE

         ! Grid box longitude/latitude [degrees]
         MIDLAT = PI*GET_YMID( I, J, 1 )/180.e+0_fp
         DayFrac = Max(0.e+0_fp,State_Met%SUNCOSmid(I,J))

         ! Additional check
         If (NOX_J(I,J,L,JNO3IDX).lt.1.0e-20) Then
            DayFrac = 0.0e+0_fp
         End If

         ! Retrieve air mass (kg) and local temperature (K)
         T3K   = State_Met%T(I,J,L)
         TINV  = 1.e+0_fp/T3K

         ! Calculate air number density (molec/cm3)
         NDAIR = GET_PCENTER(I,J,L)*1.d-4*AVOGADRO/(R_UNIV*T3K)
         XAIR  = NDAIR / State_Met%AD(I,J,L)

         ! Get local concentrations in molec/cm3
         ! Get everything except atomic N from tracer array. N is not
         ! included in any tracer family, so we can safely use it
         ! directly. If it was included in (eg) NO, we would need to
         ! subtract OLD_N from OLD_NO to get the "NO" portion of OLD_NO

         ! The "species" array is not trustworthy in this version of GC
         !OLD_N    = Spc(I,J,L,IDN)
         OLD_N = 0.0e+0_fp

         ! Calculate local NOx, N2O and O3 number density
         OLD_NO   = STT(I,J,L,IDTNO ) * TCVV(IDTNO ) * XAIR
         OLD_NO2  = STT(I,J,L,IDTNO2) * TCVV(IDTNO2) * XAIR
         OLD_NO3  = STT(I,J,L,IDTNO3) * TCVV(IDTNO3) * XAIR
         OLD_N2O  = STT(I,J,L,IDTN2O) * TCVV(IDTN2O) * XAIR

         ! Get the totals in molec/cm3
         LOCALNOX = OLD_N + OLD_NO + OLD_NO2 + OLD_NO3
         LOCALN2O = OLD_N2O
 
         ! Retrieve O3 from model
         LOCALO3  = STT(I,J,L,IDTO3)*TCVV(IDTO3)*XAIR

         ! Store this for later
         kgN2O = STT(I,J,L,IDTN2O)

         ! Reset reaction rates
         RRATE = 0e+0_fp

         ! Get the rates which are needed for both day and night-time
         ! chemistry
         ! 2:  NO + O3 -> NO2 + O2
         RRATE(2)  = 3.0e-12_fp*exp(-1500.e+0_fp*TINV)
         ! 3:  NO2 + O3 -> NO3 + O2
         RRATE(3)  = 1.2e-13_fp*exp(-2450.e+0_fp*TINV)

         ! Is it daytime?
         If (DayFrac.lt.1.0e-10_fp) Then
            ! Night-time chemistry
            ! Assume N -> NO instantaneously
            ! Calculate NO -> NO2 -> NO3 as a slow conversion
            ! Calculation below is analytically exact
            NEW_N   = 0.0e+0_fp
            NEW_NO  = (OLD_NO + OLD_N)*dexp(-dtChem*RRATE(2)*localO3)
            NEW_NO2 = OLD_NO2*dexp(-dtChem*RRATE(3)*localO3) + 
     &              ((OLD_NO + OLD_N)/((RRATE(3)/RRATE(2))-1.0e+0_fp))
     &             *(dexp((-dtChem*RRATE(2)*localO3)) - 
     &               dexp((-dtChem*RRATE(3)*localO3)))
            NEW_NO3 = OLD_NO3 + (OLD_NO2 + OLD_NO + OLD_N) - 
     &                          (NEW_NO2 + NEW_NO + NEW_N)
            FRACN   = 0.0e+0_fp
            FRACNO  = NEW_NO /localNOx
            FRACNO2 = NEW_NO2/localNOx
            FRACNO3 = NEW_NO3/localNOx
         Else
            ! Day-time chemistry
            ! Calculate reaction rates
            ! 1:  NO2 + O -> NO + O2
            RRATE(1)  = 5.1e-12_fp*exp(210.e+0_fp*TINV)
            ! 4:  NO2 + hv -> NO + O1D
            RRATE(k_JNO2) = NOX_J(I,J,L,JNO2IDX)*DAYFRAC
            ! 5:  NO3 + hv -> NO2 + O
            RRATE(k_JNO3) = NOX_J(I,J,L,JNO3IDX)*DAYFRAC
            ! 6:  NO + hv -> N + O
            RRATE(k_JNO ) = NOX_J(I,J,L,JNOIDX)*DAYFRAC
            ! 7:  N + NO2 -> N2O + O
            RRATE(7)  = 5.8e-12_fp*exp(220.e+0_fp*TINV)
            ! 8:  N + NO -> N2 + O
            RRATE(8)  = 2.1e-11_fp*exp(100.e+0_fp*TINV)
            ! 9:  N + O2 -> NO + O
            RRATE(9)  = 1.5e-11_fp*exp(-3600.e+0_fp*TINV)
            ! 10:  N2O + O1D -> N2 + O2
            RRATE(10) = 4.63e-11_fp*exp(20.e+0_fp*TINV)
            ! 11:  N2O + O1D -> 2NO
            RRATE(11) = 7.25e-11_fp*exp(20.e+0_fp*TINV)
            ! 12:  N2O + hv -> N2 + O1D
            RRATE(k_JN2O) = NOX_J(I,J,L,JN2OIDX)*DAYFRAC

            ! Sanity check
            Where(RRate.lt.0.0e+0_fp) RRate = 0.0e+0_fp

            ! Retrieve local O3P/O1D mixing ratios and relevant
            ! j-rates from interpolated 2D arrays
            LOCALO3P = NOX_O(I,J,L,1)*NDAIR*DAYFRAC
            LOCALO1D = NOX_O(I,J,L,2)*NDAIR*DAYFRAC

            ! Partition NOx into N, NO, NO2 and NO3 based on PSSA
            ! Two cases: Daytime/nighttime
            NO_ALPHA =          RRATE(k_JNO)/
     &                     (RRATE(9)*0.21e+0_fp*NDAIR)
            NO_BETA  = (RRATE(k_JNO2)+(RRATE(1)*LOCALO3P))/
     &                        (RRATE(2)*LOCALO3)
            NO_GAMMA =        (RRATE(3)*LOCALO3)/
     &                            RRATE(k_JNO3)

            ! Calculate the partition fractions
            FRACNO2 = 1.e+0_fp/(1.e+0_fp+NO_GAMMA+(NO_BETA*
     &               (1.e+0_fp+NO_ALPHA)))
            FRACNO3 = NO_GAMMA * FRACNO2
            FRACNO  = NO_BETA  * FRACNO2
            FRACN   = NO_ALPHA * FRACNO

            ! Estimate net production rates for NOx and N2O
            NOXRATE = 2.e+0_fp*((RRATE(11)*LOCALN2O*LOCALO1D) -
     &       (((RRATE(7)*FRACN*FRACNO2)+(RRATE(8)*FRACN*FRACNO))
     &                                      *LOCALNOX*LOCALNOX))
            N2ORATE = (RRATE(7)*FRACN*FRACNO2*LOCALNOX*LOCALNOX) - 
     &      ((RRATE(k_JN2O)+((RRATE(10)+RRATE(11))*LOCALO1D))*LOCALN2O)
        
            ! Get kg NOx (pre-chemistry) on an NO mass basis
            kgNOx = LOCALNOX / (XAIR * TCVV(IDTNO))

            ! Calculate fractional loss and apply to tracer
            ! Explicit Euler method (fast)
            ! These are in kg(NO) and kg(N2O) respectively
            DNOX = NOXRATE * DTCHEM / (TCVV(IDTNO)*XAIR)
            DN2O = N2ORATE * DTCHEM / (TCVV(IDTN2O)*XAIR)
            IF ((DNOX*-1e+0_fp).gt.KGNOX) THEN
               DNOX = -1e+0_fp*KGNOX
               KGNOX = 0e+0_fp
               LOCALNOX = 0.0e+0_fp
            ELSE
               KGNOX = KGNOX + DNOX
               LOCALNOX = LOCALNOX + (NOXRATE*DTCHEM)
            ENDIF
            If ((DN2O*(-1e+0_fp)).gt.KGN2O) Then
               DN2O = -1.0e+0_fp*KGN2O
               KGN2O = 0e+0_fp
            Else
               KGN2O = KGN2O + DN2O
            End If
            NEW_N   = FRACN   * LOCALNOX
            NEW_NO  = FRACNO  * LOCALNOX
            NEW_NO2 = FRACNO2 * LOCALNOX
            NEW_NO3 = FRACNO3 * LOCALNOX
         End If ! DayFrac < 0

         ! For now, say N = NO
         fracNO   = fracN + fracNO
         fracN    = 0.0e+0_fp

         ! Convert back to kg (except N, which is not a tracer)
         !NEW_N    = localNOx * fracN  !/ (XAIR * TCVV(IDTN))
         NEW_NO   = localNOx * fracNO / (XAIR * TCVV(IDTNO))
         NEW_NO2  = localNOx * fracNO2/ (XAIR * TCVV(IDTNO2))
         NEW_NO3  = localNOx * fracNO3/ (XAIR * TCVV(IDTNO3))
         NEW_N2O  = kgN2O

         !Spc(I,J,L,IDN)    = NEW_N   ! molec/cm3
         STT(I,J,L,IDTNO)  = NEW_NO  ! kg
         STT(I,J,L,IDTNO2) = NEW_NO2 ! kg
         STT(I,J,L,IDTNO3) = NEW_NO3 ! kg
         STT(I,J,L,IDTN2O) = NEW_N2O ! kg

         ! NOTE: This is not 100% reliable
         ! To be certain, bracket this in an OMP CRITICAL clause
         ! However, this will result in a speed penalty!
         MESONOX_DELTA = MESONOX_DELTA + DNOX
         MESON2O_DELTA = MESON2O_DELTA + DN2O
         
      ENDDO ! J
      ENDDO ! I
      ENDDO ! L
!$OMP END PARALLEL DO

      IF (LPRT) THEN
         ! Print mean NOx tendency in mesosphere
         DBGMSG = ' ### UCX_NOX: Mesospheric NOx processed'
         CALL DEBUG_MSG(TRIM(DBGMSG))
         WRITE(DBGMSG,'(a,2(1x,F10.5),a)') ' ### Time (start/step):',
     &       GMTSTART, DTCHEM_HR, ' hours'
         CALL DEBUG_MSG(TRIM(DBGMSG))
         WRITE(DBGMSG,'(a,1x,E10.4,1x,a)') ' ### NOx direct delta: ',
     &       MESONOX_DELTA, ' kg'
         CALL DEBUG_MSG(TRIM(DBGMSG))
         WRITE(DBGMSG,'(a,1x,E10.4,1x,a)') ' ### N2O direct delta: ',
     &       MESON2O_DELTA, ' kg'
         CALL DEBUG_MSG(TRIM(DBGMSG))
      ENDIF

      ! Free pointers
      NULLIFY( STT )
      NULLIFY( Spc )

      END SUBROUTINE UCX_NOX
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_noxcoeff
!
! !DESCRIPTION: Subroutine GET\_NOXCOEFF reads in O1D and O3P mixing  
!  ratios along with NO, NO2, NO3 and N2O J-rates from 2D data, 
!  interpolating onto the 3D grid and storing in NOX\_O and NOX\_J.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_NOXCOEFF( TARG_MONTH, Input_Opt, State_Met )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_CHEMGRID
      USE CHEMGRID_MOD,       ONLY : GET_CHEMGRID_LEVEL
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE FILE_MOD,           ONLY : IOERROR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: TARG_MONTH
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !REMARKS:
!  At some later point we should attempt to rewrite the parallel DO loop so
!  that the loop order is L-J-I.  Not sure how easy that is. (bmy, 2/14/14)
!
! !REVISION HISTORY: 
!  26 Mar 2013 - S. D. Eastham - Initial version
!  14 Feb 2014 - R. Yantosca   - Parallelize main DO loop 
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt and State_Met objects via
!                                the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                               :: IOS
      INTEGER                               :: ILEV
      INTEGER                               :: I, J, L, ITRAC
      INTEGER                               :: JJNOX
      REAL(fp)                              :: PCENTER
      REAL(fp)                              :: CURRVAL
      INTEGER                               :: VERTCOUNT, AS
      LOGICAL                               :: FOUNDLEV,EXTRAP
      INTEGER                               :: IU_FILE, fId
      INTEGER                               :: st3d(3), ct3d(3) ! Start/count
      LOGICAL                               :: ISRATE
      CHARACTER(LEN=255)                    :: NOX_FILE
      CHARACTER(LEN=255)                    :: TARG_TRAC
      CHARACTER(LEN=255)                    :: DBGMSG
      REAL(fp), DIMENSION(:,:,:), POINTER   :: NOXDATA2D => NULL()
      REAL(fp), DIMENSION(:,:), ALLOCATABLE :: NOXD2D_IN
      INTEGER                               :: LSTART

      ! Local variables for quantities from Input_Opt
      LOGICAL                               :: LPRT

      !=================================================================
      ! GET_NOXCOEFF begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      LPRT = Input_Opt%LPRT

      ! Clear interpolated arrays
      NOX_O = 0e+0_fp
      NOX_J = 0e+0_fp

      IF (UCXNETCDF) THEN
         ! Allocate and zero 2D data array
         ALLOCATE( NOXDATA2D( JJPAR, 51, 6 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOXDATA2D' )
         NOXDATA2D = 0e+0_fp

         ! Allocate and zero 2D input array
         ALLOCATE( NOXD2D_IN( UCX_NLAT, UCX_NLEVS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOXD2D_IN' )
         NOXD2D_IN = 0e+0_fp
         
         NOX_FILE = TRIM(NOON_FILE_ROOT)

         IF (LPRT) THEN
            WRITE(DBGMSG,'(a,a,I02,a,a)') ' ### UCX: Reading ', 
     &           'mesospheric NOx coeff for month ', MONTH_STR_INT ,
     &           ' from ', TRIM( NOX_FILE )
            CALL DEBUG_MSG( TRIM(DBGMSG) )
         ENDIF

         CALL NcOp_Rd (fId,TRIM(NOX_FILE))

         ! Start and count indices
         st3d = (/ 1,        1,          TARG_MONTH /)
         ct3d = (/ UCX_NLAT, UCX_NLEVS,  1          /)

         DO ITRAC = 1,6
            SELECT CASE (ITRAC)
               CASE ( 1 )
                  TARG_TRAC = 'O'
               CASE ( 2 )
                  TARG_TRAC = 'O1D'
               CASE ( 3 )
                  TARG_TRAC = 'JNO'
               CASE ( 4 )
                  TARG_TRAC = 'JNO2'
               CASE ( 5 )
                  TARG_TRAC = 'JNO3'
               CASE ( 6 )
                  TARG_TRAC = 'JN2O'
            END SELECT
            CALL NcRd( NOXD2D_IN, fId, TRIM(TARG_TRAC), st3d, ct3d )

            ! Regrid each level from 19 to JJPAR latitudes
            ! Precalculated matrix for simple linear algebra
            ! Need to reverse levels - have layer 1 = TOA
            DO ILEV = 1,51
               NOXDATA2D(:,UCX_NLEVS+1-ILEV,ITRAC) = 
     &            MATMUL(UCX_REGRID,NOXD2D_IN(:,ILEV))
            ENDDO

         ENDDO

         DEALLOCATE( NOXD2D_IN )
         CALL NcCl( fId )
      ELSE

          ! All the coefficients are now stored in NOXCOEFF. 
          ! NOXDATA2D points to the desired month slice.
          ! (ckeller, 05/12/14)
          NOXDATA2D => NOXCOEFF(:,:,:,TARG_MONTH)

      ENDIF

      ! Scan through target array, element by element
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,      J,      L,       VERTCOUNT )
!$OMP+PRIVATE( EXTRAP, LSTART, PCENTER, FOUNDLEV  )
!$OMP+PRIVATE( ITRAC,  ISRATE, CURRVAL, JJNOX     )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! The data array has already been regridded, so no 
         ! re-indexing is necessary (SDE 2015-08-15)
         JJNOX = J

         ! Vertcount is the layer count for the input, where layer 1
         ! is at the *top* of the atmosphere
         VERTCOUNT = 51
         EXTRAP = .TRUE.
         LSTART = GET_CHEMGRID_LEVEL( I, J, State_Met)

         DO L = LSTART, LLPAR
            IF (ITS_IN_THE_CHEMGRID(I,J,L,State_Met)) CYCLE

            ! Pressure at center of cell
            PCENTER = GET_PCENTER(I,J,L)
            FOUNDLEV = (PCENTER.gt.UCX_PLEVS(VERTCOUNT))
            DO WHILE (.not. FOUNDLEV)
               IF (VERTCOUNT.eq.1) THEN
                 ! At top layer; use it anyway
                 FOUNDLEV = .TRUE.
                 EXTRAP = .TRUE.
               ELSE
                 VERTCOUNT = VERTCOUNT - 1
                 FOUNDLEV = (PCENTER.gt.UCX_PLEVS(VERTCOUNT))
                 EXTRAP = .FALSE.
               ENDIF
            ENDDO

            DO ITRAC=1,6
               ISRATE = (ITRAC.gt.2)
               ! Interpolate data
               IF (EXTRAP) THEN
                  ! Just take outside value if at edges
                  CURRVAL = NOXDATA2D(JJNOX,VERTCOUNT,ITRAC)
               ELSE
                  ! Interpolate by pressure
                  CURRVAL = (UCX_PLEVS(VERTCOUNT+1)-PCENTER)
                  CURRVAL = CURRVAL/(UCX_PLEVS(VERTCOUNT+1)-
     &                               UCX_PLEVS(VERTCOUNT))
                  CURRVAL = CURRVAL * 
     &               (NOXDATA2D(JJNOX,VERTCOUNT+1,ITRAC)-
     &                NOXDATA2D(JJNOX,VERTCOUNT,ITRAC))
                  CURRVAL = NOXDATA2D(JJNOX,VERTCOUNT+1,ITRAC) - CURRVAL
               ENDIF
               IF (.not.ISRATE) THEN
                  ! Reading in mixing ratios (v/v)
                  NOX_O(I,J,L,ITRAC) = CURRVAL
               ELSE
                  ! J-rate (no conversion necessary)
                  NOX_J(I,J,L,ITRAC-2) = CURRVAL
               ENDIF
            ENDDO ! ITRAC
         ENDDO ! L

      ENDDO ! I
      ENDDO ! J
!$OMP END PARALLEL DO

      ! Cleanup
      IF (UCXNETCDF) THEN
         IF ( ASSOCIATED(NOXDATA2D) ) DEALLOCATE( NOXDATA2D )
      ELSE 
         NOXDATA2D => NULL()
      ENDIF

      END SUBROUTINE GET_NOXCOEFF
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apply_2dtrac
!
! !DESCRIPTION: Subroutine APPLY\_2DTRAC reads in and applies 2D data 
!  from AER 2D model output to the 3D GEOS-Chem grid. Zonally-averaged
!  mixing ratios vary by month, but no interpolation is performed.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE APPLY_2DTRAC( TRAC_ADD, STRAT_ONLY, SFC_ONLY, SCALEFAC,
     &                         Input_Opt, State_Met, State_Chm )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_STRATMESO
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE FILE_MOD,           ONLY : IOERROR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: TRAC_ADD    ! Add (instead of overwrite)
      LOGICAL,        INTENT(IN)    :: STRAT_ONLY  ! Only write to stratosphere
      LOGICAL,        INTENT(IN)    :: SFC_ONLY    ! Only write to surface
      REAL(fp),       INTENT(IN)    :: SCALEFAC    ! Future scaling factor
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REMARKS:
!  At some later point we should attempt to rewrite the parallel DO loop so
!  that the loop order is L-J-I.  Not sure how easy that is. (bmy, 2/14/14)
!
! !REVISION HISTORY: 
!  26 Mar 2013 - S. D. Eastham - Initial version
!  14 Feb 2014 - R. Yantosca   - Now parallelize main DO loop
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm 
!                                objects via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: IOS
      INTEGER            :: ILEV
      INTEGER            :: I, J, L
      REAL(fp)           :: PCENTER
      REAL(fp)           :: TRACMR,TRACMASS
      INTEGER            :: VERTCOUNT
      LOGICAL            :: FOUNDLEV,EXTRAP
      REAL(fp)           :: TRACMULT
      LOGICAL            :: SET_CELL
      INTEGER            :: IU_FILE
      INTEGER            :: L_MAX
      CHARACTER(LEN=255) :: DBGMSG
      INTEGER :: fId
      INTEGER :: st3d(3), ct3d(3)   ! Start + count, for 3D arrays 
      REAL(fp)           :: TEMP_225(91,1)
      !=====SDE DEBUG
      LOGICAL,SAVE       :: FIRST=.TRUE.
      !=====SDE DEBUG

      ! Local variables for quantities from Input_Opt
      LOGICAL            :: LPRT
      INTEGER            :: N_TRACERS
      REAL(fp)           :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      REAL(fp), POINTER  :: STT(:,:,:,:)

      !=================================================================
      ! APPLY_2DTRAC begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      LPRT      = Input_Opt%LPRT
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV = Input_Opt%TCVV(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg]
      STT => State_Chm%Tracers

      ! Filename
      TRAC_FILE = TRIM(AVG_FILE_ROOT)//TRIM(TRAC_SET )
      TRAC_FILE = TRIM(TRAC_FILE)//'_'//TRIM(MONTH_STR)//'.dat'

      ! Error check
      IF ((STRAT_ONLY).and.(SFC_ONLY)) THEN
         CALL ERROR_STOP('Cannot have both SFC_ONLY and STRAT_ONLY!',
     &                           'APPLY_2DTRAC(UCX_mod.F)')
      ENDIF

      IF (UCXNETCDF) THEN
         ! Read from NetCDF
         TRAC_FILE = TRIM(AVG_FILE_ROOT)

         IF (LPRT) THEN
            WRITE(DBGMSG,'(a,a,a,I02,a,a)') ' ### UCX: Reading ', 
     &           TRIM(TRAC_SET), ' for month ', MONTH_STR_INT ,
     &           ' from ', TRIM( TRAC_FILE )
            CALL DEBUG_MSG( TRIM(DBGMSG) )
         ENDIF

         ! By the end of this, need to have data in the form
         ! of an array size JJPARx51 [LAT x LEV]

         ! UCX_MR:      JJPARx51
         ! UCX_MR_IN:   19x51
         CALL NcOp_Rd (fId,TRIM(TRAC_FILE))

         ! Start and count indices
         st3d = (/ 1,                1,  MONTH_STR_INT /)
         ct3d = (/ UCX_NLAT, UCX_NLEVS,  1             /)

         CALL NcRd( UCX_MR_IN, fId, TRIM(TRAC_SET), st3d, ct3d )
         CALL NcCl( fId )
 
         ! Regrid each level from 19 to JJPAR latitudes
         ! Precalculated matrix for simple linear algebra
         ! Need to reverse levels - have layer 1 = TOA
         DO ILEV = 1,UCX_NLEVS
            UCX_MR(:,UCX_NLEVS+1-ILEV) = 
     &         MATMUL(UCX_REGRID,UCX_MR_IN(:,ILEV))
         ENDDO
      ELSE
         TRAC_FILE = TRIM(AVG_FILE_ROOT)//TRIM(TRAC_SET )
         TRAC_FILE = TRIM(TRAC_FILE)//'_'//TRIM(MONTH_STR)//'.dat'

         ! Get a free LUN
         IU_FILE = findFreeLUN()

         IOS = 1
         OPEN( IU_FILE, FILE=TRIM( TRAC_FILE ),STATUS='OLD',IOSTAT=IOS) 
         IF ( IOS /= 0 ) THEN
            WRITE(6,*) 'UCX_MOD: Could not read ', TRIM( TRAC_FILE )
            CALL IOERROR( IOS, IU_FILE,'UCX_MOD:APPLY_2DTRAC')
         ENDIF
  
         IF (LPRT) THEN
            WRITE(DBGMSG,'(a,a)') ' ### UCX: Reading ', 
     &           TRIM( TRAC_FILE )
            CALL DEBUG_MSG( TRIM(DBGMSG) )
         ENDIF
 
         ! Read in mixing ratios
         ! SDE 2015-08-18: Now always read in at 2x2.5 and regrid
         DO ILEV = 1,51
            ! Reverse the level order to be compatible with the NetCDF
            READ(IU_FILE, 110, IOSTAT=IOS ) TEMP_225(:,1)
 
            UCX_MR(:,ILEV) = MATMUL(UCX_REGRID_225,TEMP_225(:,1))
            IF ( IOS /= 0 ) THEN
               WRITE(6,'(a,I4,a,1x,a)') 'UCX_MOD: Error reading line '
     &    , ILEV, ' in file ', TRIM( TRAC_FILE )
               CALL IOERROR( IOS, IU_FILE,'UCX_MOD:APPLY_2DTRAC')
            ENDIF
         ENDDO

110   FORMAT(91E10.3)

         CLOSE(IU_FILE)    
      ENDIF

      IF (SFC_ONLY) THEN
         L_MAX = 1
      ELSE
         L_MAX = LLPAR
      ENDIF

      ! Scan through STT, element by element
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,        J,       L,        VERTCOUNT, EXTRAP   )
!$OMP+PRIVATE( SET_CELL, PCENTER, FOUNDLEV, TRACMR,    TRACMASS )
      DO J = 1,JJPAR
         DO I = 1,IIPAR
            ! Vertcount is the layer count for the input, where layer 1
            ! is at the *top* of the atmosphere
            VERTCOUNT = 51
            EXTRAP = .TRUE.
            DO L = 1, L_MAX
               IF (STRAT_ONLY) THEN
                  ! Only want to set non-trop cells
                  SET_CELL = ITS_IN_THE_STRATMESO(I,J,L,State_Met)
               ELSE
                  SET_CELL = .TRUE.
               ENDIF
               ! Pressure at center of cell
               PCENTER = GET_PCENTER(I,J,L)
               FOUNDLEV = (PCENTER.gt.UCX_PLEVS(VERTCOUNT))
               DO WHILE (.not. FOUNDLEV)
                  IF (VERTCOUNT.eq.1) THEN
                    ! At top layer; use it anyway
                    FOUNDLEV = .TRUE.
                    EXTRAP = .TRUE.
                  ELSE
                    VERTCOUNT = VERTCOUNT - 1
                    FOUNDLEV = (PCENTER.gt.UCX_PLEVS(VERTCOUNT))
                    EXTRAP = .FALSE.
                  ENDIF
               ENDDO

               ! Interpolate mixing ratio
               IF (EXTRAP) THEN
                  ! Just take outside value if at edges
                  TRACMR = UCX_MR(J,VERTCOUNT)
               ELSE
                  ! Interpolate by pressure
                  TRACMR = (UCX_PLEVS(VERTCOUNT+1)-PCENTER)
                  TRACMR = TRACMR/(UCX_PLEVS(VERTCOUNT+1)-
     &                               UCX_PLEVS(VERTCOUNT))
                  TRACMR = TRACMR * 
     &               (UCX_MR(J,VERTCOUNT+1)-UCX_MR(J,VERTCOUNT))
                  TRACMR = UCX_MR(J,VERTCOUNT+1) - TRACMR
               ENDIF
               ! Send to STT (currently kg/box)
               IF (SET_CELL) THEN
                  ! Convert from [ppbv] to [kg]
                  TRACMASS = SCALEFAC * TRACMR * State_Met%AD(I,J,L) /
     &                       TCVV(TRAC_IDX)
                  IF (TRAC_ADD) THEN
                     STT(I,J,L,TRAC_IDX) = STT(I,J,L,TRAC_IDX)+TRACMASS
                  ELSE
                     STT(I,J,L,TRAC_IDX) = TRACMASS
                  ENDIF
                ENDIF
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
 
      ! Free pointer
      NULLIFY( STT )

      ! Return to calling program
      END SUBROUTINE APPLY_2DTRAC
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: settle_strat_aer
!
! !DESCRIPTION: Subroutine SETTLE\_STRAT\_AER performs gravitational settling
!  of stratospheric aerosols. It is copied largely from GRAV\_SETTLING in
!  sulfate\_mod.F. All of this is ignored if APM is active.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETTLE_STRAT_AER( am_I_Root, Input_Opt, State_Met, 
     &                             State_Chm, RC )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_STRATMESO
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_TROP
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_CHEMGRID
      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND44
      USE CMN_GCTM_MOD             ! g0
      USE CMN_FJX_MOD,        ONLY : RAA, IND999
      USE COMODE_MOD,         ONLY : JLOP,WERADIUS,ERADIUS
      USE ERROR_MOD,          ONLY : IT_IS_NAN,ERROR_STOP,SAFE_DIV
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_ErrCode_Mod
      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE TIME_MOD,           ONLY : GET_ELAPSED_SEC, GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTNIT,IDTBCPI,IDTH2O,IDTHNO3
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      INTEGER,        INTENT(INOUT) :: RC          ! Return code
!
! !REVISION HISTORY: 
!  11 Apr 2013 - S. D. Eastham - Initial version
!  07 Feb 2014 - R. Yantosca   - Add missing variables to !$OMP PRIVATE clause
!  07 Feb 2014 - R. Yantosca   - Cosmetic changes
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm 
!                                objects via the arg list
!  08 Apr 2015 - R. Yantosca   - Remove call to READ_PSC_FILE, this is
!                                now done from DO_CHEMISTRY (chemistry_mod.F)  
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Common
      INTEGER                :: I,      J,     L,        DTCHEM
      REAL(fp)               :: DELZ,   DELZ1
      REAL(fp)               :: P,      TEMP,  DP,       PDP     
      REAL(fp)               :: CONST,  SLIP,  VISC,     FAC1
      REAL(fp)               :: FAC2,   FLUX,  AREA_CM2, RHB
      REAL(fp)               :: RCM
      REAL(fp)               :: TOT1,   TOT2
      INTEGER                :: KLOOP, K
      LOGICAL                :: NATCOL, NATCELL
     
      ! Specific to each class
      REAL(fp)               :: RWET(2),CONST_V(2)
      REAL(fp)               :: RHO(2),RATIO_R(2),REFF(2)
      REAL(fp)               :: VTS(LLPAR,2)
 
      ! Used for old Seinfeld & Pandis slip factor calc
      REAL(fp)               :: sp_Lambda, sp_Num

      ! Parameters
#if defined( APM )
      REAL(fp), PARAMETER    :: BCDEN = 1800.e+0_fp
#else
      REAL(fp), PARAMETER    :: BCDEN = 1000.e+0_fp ! density (kg/m3)
#endif

      ! Indexing
      INTEGER, PARAMETER     :: IBC  = 1
      INTEGER, PARAMETER     :: ILIQ = 2
      INTEGER, PARAMETER     :: NSETTLE = 2
      INTEGER                :: IAERO
      LOGICAL                :: RUNCALC

      ! Partitioning
      REAL(fp)               :: PHASEMASS(3,2)
      REAL(fp)               :: SEDMASS
      INTEGER                :: IDTCURRENT

      ! NAT only
      REAL(fp)               :: VFALLMAX
      REAL(fp)               :: VNAT(LLPAR)
      INTEGER                :: MAXALT, MINALT
      REAL(fp)               :: BXMIN, SEDSTEP, TEMPREAL, P_0, P_ABOVE
      REAL(fp)               :: PSED0, PSEDABOVE, INVAIR_0, INVAIR_ABOVE
      REAL(fp)               :: XPSC_ABOVE, XPSC_0, XPSC_BELOW
      REAL(fp)               :: XNAT_ABOVE, XNAT_0, XNAT_BELOW
      REAL(fp)               :: XNO3_ABOVE, XNO3_0, XNO3_BELOW
      REAL(fp)               :: XICE_ABOVE, XICE_0, XICE_BELOW
      REAL(fp)               :: SEDPSC, SEDNAT, SEDICE, SEDQUANT
      REAL(fp)               :: SEDH2O, SEDNO3
      REAL(fp)               :: BELOWGRAD, ABOVEGRAD
      INTEGER                :: LOCALPROFILE, NUMSEDSTEPS, STARTPT, ISED

      ! Log-normal aerosols
      REAL(fp)               :: CHENM2R2, CHENM2A2R2
      REAL(fp)               :: RCALC, RHOCALC, LNSSQ

      ! SDE 2013-07-07: Debug variables
      REAL(fp)               :: D_XH2O, D_PPA,  D_TK,   D_XHNO3,D_XNAT
      REAL(fp)               :: D_ICE1, D_ICE0, D_H2O1, D_H2O0
      REAL(fp)               :: D_NAT1, D_NAT0, D_SEDQ, D_PSED, D_HNO3
      REAL(fp)               :: D_RMAX, D_VMAX, MAXLOSS, FRACLOSS

      ! Local variables for quantities from Input_Opt
      LOGICAL                :: LPRT
      LOGICAL                :: LWACCMAER
      LOGICAL                :: LGRAVSTRAT
      INTEGER                :: N_TRACERS
      REAL(fp)               :: XNUMOL(Input_Opt%N_TRACERS)
      REAL(fp)               :: TRACER_MW_G(Input_Opt%N_TRACERS)

      ! Pointers
      REAL(fp), POINTER      :: STT(:,:,:,:)

      !=================================================================
      ! SETTLE_STRAT_AER begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT
      LPRT        = Input_Opt%LPRT
      LGRAVSTRAT  = Input_Opt%LGRAVSTRAT
      LWACCMAER   = Input_Opt%LWACCMAER
      N_TRACERS   = Input_Opt%N_TRACERS
      XNUMOL      = Input_Opt%XNUMOL
      TRACER_MW_G = Input_Opt%TRACER_MW_G(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg]
      STT => State_Chm%Tracers

      ! Return if gravitational settling disabled
      IF (.not. LGRAVSTRAT) RETURN

      ! Return if it's the start of the run
      IF ( GET_ELAPSED_SEC() == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60e+0_fp

      ! For debugging
      MAXLOSS= 0e+0_fp
      D_XHNO3= 0e+0_fp
      D_XNAT = 0e+0_fp
      D_PPA  = 0e+0_fp
      D_TK   = 0e+0_fp
      D_XH2O = 0e+0_fp
      D_SEDQ = 0e+0_fp
      D_PSED = 0e+0_fp
      D_ICE0 = 0e+0_fp
      D_ICE1 = 0e+0_fp
      D_H2O0 = 0e+0_fp
      D_H2O1 = 0e+0_fp
      D_NAT0 = 0e+0_fp
      D_NAT1 = 0e+0_fp
      D_HNO3 = 0e+0_fp
      D_RMAX = 0e+0_fp
      D_VMAX = 0e+0_fp

      ! First settle liquid aerosols (SLA) using scheme found
      ! elsewhere in GEOS-Chem
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED  )
!$OMP+PRIVATE( J,            I,            L,          VTS        )
!$OMP+PRIVATE( VNAT,         NATCOL,       BXMIN,      MINALT     )
!$OMP+PRIVATE( MAXALT,       VFALLMAX,     TEMP,       P          )
!$OMP+PRIVATE( RUNCALC,      KLOOP,        RWET,       RHO        )
!$OMP+PRIVATE( SP_NUM,       SP_LAMBDA,    VISC,       IAERO      )
!$OMP+PRIVATE( DP,           PDP,          CONST,      SLIP       )
!$OMP+PRIVATE( DELZ,         CONST_V,      PHASEMASS,  K          )
!$OMP+PRIVATE( IDTCURRENT,   SEDMASS,      DELZ1,      SEDSTEP    )
!$OMP+PRIVATE( NUMSEDSTEPS,  TEMPREAL,     ISED,       STARTPT    )
!$OMP+PRIVATE( XNO3_0,       XNAT_0,       XICE_0,     XPSC_0     )
!$OMP+PRIVATE( INVAIR_0,     INVAIR_ABOVE, XNO3_ABOVE, XNAT_ABOVE )
!$OMP+PRIVATE( XICE_ABOVE,   XPSC_ABOVE,   P_ABOVE,    ABOVEGRAD  )
!$OMP+PRIVATE( XNO3_BELOW,   XNAT_BELOW,   XICE_BELOW, XPSC_BELOW )
!$OMP+PRIVATE( P_0,          PSED0,        PSEDABOVE,  BELOWGRAD  )
!$OMP+PRIVATE( LOCALPROFILE, SEDQUANT,     SEDNO3,     SEDICE     )
!$OMP+PRIVATE( SEDH2O,       SEDNAT,       SEDPSC,     FRACLOSS   )
!$OMP+PRIVATE( MAXLOSS,      D_ICE0,       D_ICE1,     D_H2O1     )
!$OMP+PRIVATE( D_H2O0,       D_NAT1,       D_NAT0,     D_HNO3     )
!$OMP+PRIVATE( D_RMAX,       D_VMAX,       D_SEDQ,     D_PSED     )
!$OMP+PRIVATE( D_XH2O,       D_XHNO3,      D_XNAT,     D_PPA      )
!$OMP+PRIVATE( D_TK,         CHENM2R2,     CHENM2A2R2, RHOCALC    )
!$OMP+PRIVATE( LNSSQ,        RCALC,        NATCELL                )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR       

         ! Initialize 
         DO L = 1, LLPAR
            VTS(L,:) = 0e+0_fp
            VNAT(L) = 0e+0_fp
         ENDDO
         NATCOL = .FALSE.

         ! Arbitrary limits
         BXMIN = 1.e+20_fp
         MINALT = LLPAR+1
         MAXALT = 0
         VFALLMAX = 0e+0_fp

         ! Loop over levels
         DO L = 1, LLPAR

            ! Temperature [K]
            TEMP    = State_Met%T(I,J,L)
 
            ! Pressure at center of the level [kPa]
            P       = GET_PCENTER(I,J,L) * 0.1e+0_fp

            RUNCALC = ITS_IN_THE_STRATMESO(I,J,L,State_Met)

            IF (RUNCALC) THEN
               ! Need to translate for BC radii
               IF (ITS_IN_THE_CHEMGRID(I,J,L,State_Met)) THEN
                  KLOOP = JLOP(I,J,L)  
                  RWET(IBC) = WERADIUS(KLOOP,2)*1.e-2_fp
               ELSE
                  ! Use defaults, assume dry (!)
                  RWET(IBC) = RAA(IND999,29)*1.e-6_fp
               ENDIF

               ! Taken from aerosol_mod (MSDENS(2))
               RHO(IBC) = BCDEN

               ! Get aerosol properties
               RWET(ILIQ) = RAD_AER(I,J,L,I_SLA)*1.e-2_fp
               RHO(ILIQ)  = RHO_AER(I,J,L,I_SLA)

               ! Do we need to sediment NAT?
               NATCELL = NDENS_AER(I,J,L,I_SPA).gt.TINY(0e+0_fp)
               IF (NATCELL) THEN
                  NATCOL = .TRUE.
                  BXMIN  = MIN(BXMIN, State_Met%BXHEIGHT(I,J,L))
                  MINALT = MIN(L,MINALT)
                  MAXALT = MAX(L,MAXALT)
               ENDIF
            ENDIF

            IF (.not.RUNCALC) THEN
               VTS(L,:) = 0e+0_fp
               VNAT(L) = 0e+0_fp
            ELSE
               ! Calculate common variables first
               sp_Num = P * 1.e+3_fp * 6.023e+23_fp / (8.314 * Temp)
               sp_Lambda = 1.e+0_fp/(1.41421*sp_Num*3.141592*
     &                    (3.7e-10_fp)**2)

               ! Viscosity [Pa*s] of air as a function of temperature 
               VISC = 1.458e-6_fp * (Temp)**(1.5e+0_fp) 
     &               / ( Temp + 110.4e+0_fp )

               IF (LWACCMAER) THEN
                  ! LBS and NAT are now handled as log-normal

                  ! BC is handled as monodisperse
                  ! Dp = particle diameter [um]
                  IF (RWET(IBC).gt.1.0e-12_fp) THEN
                     RCALC   = RWET(IBC)
                     RHOCALC = RHO(IBC)
                     CONST = 2.e+0_fp * RHOCALC * RCALC * RCALC
     &                      * g0 / 9.e+0_fp
                     ! Monodisperse slip correction
                     Slip = 1. + (sp_Lambda/RCALC) * (1.257 + 0.4*
     &                     exp(-1.1*RCALC / sp_Lambda))
                     ! Settling velocity [m/s]
                     VTS(L,IBC) = CONST * Slip / VISC 
                  ELSE
                     VTS(L,IBC) = 0.0e+0_fp
                  ENDIF

                  ! Now LBS - log-normally distributed
                  IF (RWET(ILIQ).gt.1.0e-12_fp) THEN
                     RCALC      = RWET(ILIQ)
                     RHOCALC    = RHO(ILIQ)
                     LNSSQ      = LSSQWACCM(I_SLA)
                     CONST = 2.e+0_fp*RHOCALC*g0*RCALC*RCALC/9.e+0_fp
                     ! SNAP slip calculation
                     CHENM2R2   = DEXP(2.0e+0_fp*LNSSQ)
                     CHENM2A2R2 = (RCALC ** ( - CHENA2)) *
     &                              DEXP((2.0e+0_fp - CHENA2) *
     &                              0.5e+0_fp*LNSSQ)
                     Slip = CHENM2R2 + (CHENA1*(sp_Lambda**CHENA2)*
     &                                  CHENM2A2R2)
                     VTS(L,ILIQ) = CONST * Slip / VISC
                  ELSE
                     VTS(L,ILIQ) = 0.0e+0_fp
                  ENDIF

                  ! Finally NAT/ice, also log-normally distributed
                  IF (NATCELL) THEN
                     RCALC      = RAD_AER(I,J,L,I_SPA) * 1.0e-2_fp
                     RHOCALC    = RHO_AER(I,J,L,I_SPA)
                     LNSSQ      = LSSQWACCM(I_SPA)
                     CONST = 2.e+0_fp*RHOCALC*g0*RCALC*RCALC/9.e+0_fp
                     ! SNAP slip calculation
                     CHENM2R2   = DEXP(2.0e+0_fp*LNSSQ)
                     CHENM2A2R2 = (RCALC ** ( - CHENA2)) *
     &                              DEXP((2.0e+0_fp - CHENA2) *
     &                              0.5e+0_fp*LNSSQ)
                     Slip = CHENM2R2 + (CHENA1*(sp_Lambda**CHENA2)*
     &                                  CHENM2A2R2)
                     ! Account for log-normal distribution with
                     ! Considine parameterization instead?
                     !Slip = 1. + (sp_Lambda/RCALC) * (1.257 + 0.4*
     &               !      exp(-1.1*RCALC / sp_Lambda))
                     !Slip = Slip * DEXP(8.0*LNSSQ)
                     VNAT(L) = CONST * Slip / VISC
                  ELSE
                     VNAT(L) = 0.0e+0_fp
                  ENDIF
               ELSE
                  DO IAERO=1,NSETTLE
                     IF (RWET(IAERO).le.TINY(0e+0_fp)) THEN
                        VTS(L,IAERO) = 0e+0_fp
                     ELSE
                        ! Dp = particle diameter [um]
                        DP = 2.e+0_fp * RWET(IAERO) * 1.e+6_fp        

                        ! PdP = P * dP [hPa * um]
                        !PDp     = P * Dp

                        ! Constant
                        CONST = 2.e+0_fp * RHO(IAERO) * RWET(IAERO)**2 
     &                         * g0 / 9.e+0_fp

                        !=========================================================
                        ! NOTE: Slip correction factor calculations following 
                        ! Seinfeld, pp464 which is thought to be more accurate 
                        ! but more computation required. (rjp, 1/24/02)
                        !
                        ! # air molecule number density
                        ! num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
                        !
                        ! # gas mean free path
                        ! lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
                        !
                        ! # Slip correction
                        ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp  
                        !     &     / (2. * lamda))) / Dp
                        !
                        ! NOTE: Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
                        ! which produces slip correction factore with small error
                        ! compared to the above with less computation.
                        !========================================================= 
                      
                        ! Slip correction factor (as function of P*dp)
                        ! Slip = 1.e+0_fp+(15.60e+0_fp + 7.0e+0_fp *
     &                  !           EXP(-0.059e+0_fp * PDp)) / PDp
                        ! Reverting to Seinfeld and Pandis
                        Slip = 1. + 2. * sp_Lambda * (1.257 + 0.4*
     &                        exp(-1.1*Dp / (2. * sp_Lambda))) / Dp
   
                        ! Settling velocity [m/s]
                        VTS(L,IAERO) = CONST * Slip / VISC 
                     ENDIF ! RWET
                  ENDDO ! IAERO
               
                  ! Now solid PSC particles
                  IF (NATCOL) THEN
                     ! sp_Num: Air molecule#/m3
                     VNAT(L) = CALC_FALLVEL(RHO_AER(I,J,L,I_SPA),
     &                             RAD_AER(I,J,L,I_SPA),TEMP,P)
                  ELSE
                     VNAT(L) = 0.e+0_fp
                  ENDIF
               ENDIF ! LWACCMAER

               ! Calculate maximum fall velocity
               IF (VNAT(L).gt.VFALLMAX) VFALLMAX = VNAT(L)
            ENDIF ! RUNCALC
         ENDDO

         ! First apply simpler SLA sedimentation scheme
         ! Handle model top condition
         L    = LLPAR
         DELZ = State_Met%BXHEIGHT(I,J,L)
         
         DO IAERO=1,NSETTLE
            CONST_V(IAERO) = 1.e+0_fp / (1.e+0_fp + DTCHEM
     &                      * VTS(L,IAERO) / DELZ)
         ENDDO

         ! Zero arrays
         PHASEMASS(:,:) = 0e+0_fp

         ! Only want to sediment fraction of species currently
         ! in the aerosol
         DO K = 1,7
            ! Only process transported tracers
            IDTCURRENT = AERFRACIND(K)
            IF (IDTCURRENT.ne.0) THEN
               ! Calculate local phase partitioning
               ! Total upper gridbox mass
               PHASEMASS(3,2) = STT(I,J,L,IDTCURRENT)
               ! Aerosol-phase upper gridbox mass
               PHASEMASS(2,2) = AERFRAC(I,J,L,K)*PHASEMASS(3,2)
               ! Gas-phase upper gridbox mass
               PHASEMASS(1,2) = PHASEMASS(3,2) - PHASEMASS(2,2)

               ! Calculate total sedimented mass
               SEDMASS = PHASEMASS(2,2) * (1.e+0_fp-CONST_V(ILIQ))

               ! Remove from upper gridbox
               PHASEMASS(2,2) = PHASEMASS(2,2) - SEDMASS
               PHASEMASS(3,2) = PHASEMASS(1,2) + PHASEMASS(2,2)

               ! Recalculate phase fractions
               IF (PHASEMASS(3,2).gt.TINY(1e+0_fp)) THEN
                  AERFRAC(I,J,L,K) = PHASEMASS(2,2)/PHASEMASS(3,2)
               ELSE
                  AERFRAC(I,J,L,K) = 0e+0_fp
               ENDIF

               ! Store result
               STT(I,J,L,IDTCURRENT) = PHASEMASS(3,2)
            ENDIF
         ENDDO
         STT(I,J,L,IDTBCPI)= STT(I,J,L,IDTBCPI)* CONST_V(IBC) 

         DO L = LLPAR-1,1,-1
            IF (ITS_IN_THE_TROP(I,J,L+1,State_Met)) CYCLE
            DELZ  = State_Met%BXHEIGHT(I,J,L)
            DELZ1 = State_Met%BXHEIGHT(I,J,L+1)

            DO K=1,7
               IDTCURRENT = AERFRACIND(K)
               IF (IDTCURRENT.ne.0) THEN
                  ! Total upper gridbox mass
                  PHASEMASS(3,2) = STT(I,J,L+1,IDTCURRENT)
                  ! Aerosol-phase upper gridbox mass
                  PHASEMASS(2,2) = AERFRAC(I,J,L+1,K)*PHASEMASS(3,2)
                  ! Gas-phase upper gridbox mass
                  PHASEMASS(1,2) = PHASEMASS(3,2) - PHASEMASS(2,2)
                  
                  ! Total lower gridbox mass
                  PHASEMASS(3,1) = STT(I,J,L,IDTCURRENT)
                  ! Aerosol-phase lower gridbox mass
                  PHASEMASS(2,1) = AERFRAC(I,J,L,K)*PHASEMASS(3,1)
                  ! Gas-phase lower gridbox mass
                  PHASEMASS(1,1) = PHASEMASS(3,1) - PHASEMASS(2,1)
                  
                  ! New lower gridbox mass
                  PHASEMASS(2,1) = 1.e+0_fp/(1.e+0_fp+DTCHEM
     &               * VTS(L,ILIQ) / DELZ)
     &               * (PHASEMASS(2,1)+DTCHEM*VTS(L+1,ILIQ)/DELZ1
     &               * PHASEMASS(2,2))

                  ! Calculate new total mass in lower gridbox
                  PHASEMASS(3,1) = PHASEMASS(2,1) + PHASEMASS(1,1)

                  ! Recalculate phase fraction
                  IF (PHASEMASS(3,1).gt.TINY(1e+0_fp)) THEN
                     AERFRAC(I,J,L,K) = PHASEMASS(2,1)/PHASEMASS(3,1)
                  ELSE
                     AERFRAC(I,J,L,K) = 0e+0_fp
                  ENDIF

                  ! Store result
                  STT(I,J,L,IDTCURRENT) = PHASEMASS(3,1)
               ENDIF
            ENDDO
            STT(I,J,L,IDTBCPI) = 1.e+0_fp/(1.e+0_fp+DTCHEM
     &              * VTS(L,IBC) / DELZ)
     &              * (STT(I,J,L,IDTBCPI)+DTCHEM*VTS(L+1,IBC)/DELZ1
     &              * STT(I,J,L+1,IDTBCPI)) 
         ENDDO

         ! Now perform trapezoidal scheme for particulates
         ! Calculate maximum allowable timestep (seconds)
         IF (VFALLMAX.gt.TINY(1e+0_fp)) THEN
            SEDSTEP = BXMIN/VFALLMAX
            IF(DTCHEM.le.SEDSTEP)THEN
               NUMSEDSTEPS = 1
               SEDSTEP = DTCHEM
            ELSE
               ! Will need to run iteratively
               ! Calculate minimum necessary number of steps, limiting to
               ! 10 steps if excessive
               TEMPREAL = DTCHEM/SEDSTEP 
               NUMSEDSTEPS = CEILING(TEMPREAL)
               NUMSEDSTEPS = MIN(10,NUMSEDSTEPS)
               SEDSTEP = DTCHEM/(NUMSEDSTEPS*1.e+0_fp)
               VFALLMAX = BXMIN/SEDSTEP
            ENDIF

            SEDSTEPLOOP: DO ISED=1,NUMSEDSTEPS
               STARTPT = MAX(1,MINALT-1)
               ! XPSC is the number of molecules tied up in solid particles
               ! per m3 in a grid box
               L = STARTPT
               !IF (STATE_PSC(I,J,L).lt.2) THEN
               IF (STATE_PSC(I,J,L) < 2.0e+0_f4 ) THEN
                  XNO3_0 = 0e+0_fp
                  XNAT_0 = 0e+0_fp
                  XICE_0 = 0e+0_fp
               ELSE
                  XNO3_0 = STT(I,J,L,IDTNIT) * AIRMW / 
     &                     (TRACER_MW_G(IDTNIT)*State_Met%AD(I,J,L))
                  XNAT_0 = XNO3_0! * 4.e+0_fp
                  XICE_0 = (KG_AER(I,J,L,I_SPA)-
     &                     ((NATMW/TRACER_MW_G(IDTNIT))*
     &                     STT(I,J,L,IDTNIT)))*AIRMW/(ICEMW*
     &                     State_Met%AD(I,J,L))
               ENDIF
               XPSC_0   = XNAT_0 + XICE_0
               P_0      = 100.0e+0_fp * GET_PCENTER(I,J,L)
               INVAIR_0 = AIRMW/State_Met%AD(I,J,L)
               IF (L .lt. (LLPAR-1)) THEN
                  INVAIR_ABOVE  = AIRMW/State_Met%AD(I,J,L+1)
                  !IF (STATE_PSC(I,J,L+1).lt.2) THEN
                  IF (STATE_PSC(I,J,L+1) < 2.0e+0_f4 ) THEN
                     XNO3_ABOVE = 0e+0_fp
                     XNAT_ABOVE = 0e+0_fp
                     XICE_ABOVE = 0e+0_fp
                  ELSE
                     XNO3_ABOVE = STT(I,J,L+1,IDTNIT) * INVAIR_ABOVE / 
     &                                   (TRACER_MW_G(IDTNIT))
                     ! NAT = HNO3.3H2O = 4 molecules
                     XNAT_ABOVE = XNO3_ABOVE! * 4.e+0_fp
!                     XICE_ABOVE = (KG_AER(I,J,L+1,I_SPA)-STT(I,J,L+1,IDTNIT))
!     &                                          *INVAIR_ABOVE/ICEMW
                     XICE_ABOVE = (KG_AER(I,J,L+1,I_SPA)-
     &                   ((NATMW/TRACER_MW_G(IDTNIT))*
     &                   STT(I,J,L+1,IDTNIT)))*INVAIR_ABOVE/ICEMW
                  ENDIF
                  XPSC_ABOVE = XNAT_ABOVE + XICE_ABOVE
                  P_ABOVE    = 100.0e+0_fp * GET_PCENTER(I,J,L+1)
                  PSEDABOVE  = VNAT(L+1)*SEDSTEP*1.e+2_fp*
     &                         (GET_PEDGE(I,J,L+1)-GET_PEDGE(I,J,L+2))/
     &                         State_Met%BXHEIGHT(I,J,L+1)
!                  PSEDABOVE=g0 * (1.0d3*AIRMW) * P_ABOVE * 
!     &               VNAT(L+1) / (Rd * T(I,J,L+1))*SEDSTEP
               ELSE
                  INVAIR_ABOVE = 1.0e+0_fp
                  XNO3_ABOVE = 0.0e+0_fp
                  XNAT_ABOVE = 0.0e+0_fp
                  XICE_ABOVE = 0.0e+0_fp
                  XPSC_ABOVE = 0.0e+0_fp
                  P_ABOVE = GET_PEDGE(I,J,L+2) * 100.0e+0_fp
                  PSEDABOVE = 0.0e+0_fp
               ENDIF
               ABOVEGRAD = (XPSC_ABOVE-XPSC_0)/(P_0-P_ABOVE)
               SED_LLOOP: DO L=STARTPT,MAXALT
                  ! Actually calculating sedimentation for the box above this
                  ! one. By the time we get around to processing a given
                  ! box, it is actually the (i-1)th box
                  ! We are therefore concerned with calculating
                  ! sedimentation *out* of the box above us and *into* this
                  ! one
                  XNO3_BELOW = XNO3_0
                  XNAT_BELOW = XNAT_0
                  XICE_BELOW = XICE_0
                  XPSC_BELOW = XPSC_0
                  XNO3_0 = XNO3_ABOVE
                  XNAT_0 = XNAT_ABOVE
                  XICE_0 = XICE_ABOVE
                  XPSC_0 = XPSC_ABOVE
                  P_0 = P_ABOVE
                  PSED0 = PSEDABOVE
                  INVAIR_0 = INVAIR_ABOVE
                  
                  IF (L.lt.(LLPAR-1)) THEN
                     INVAIR_ABOVE  = AIRMW/State_Met%AD(I,J,L+2)
                     !IF (STATE_PSC(I,J,L+2).ge.2) THEN
                     IF (STATE_PSC(I,J,L+2) >= 2.0e+0_f4 ) THEN
                        XNO3_ABOVE = STT(I,J,L+2,IDTNIT)*INVAIR_ABOVE/ 
     &                                TRACER_MW_G(IDTNIT)
                        XNAT_ABOVE = XNO3_ABOVE! * 4.e+0_fp
                        XICE_ABOVE = (KG_AER(I,J,L+2,I_SPA)-
     &                      ((NATMW/TRACER_MW_G(IDTNIT))*
     &                      STT(I,J,L+2,IDTNIT)))*INVAIR_ABOVE/ICEMW
                     ELSE
                        XNO3_ABOVE = 0e+0_fp
                        XNAT_ABOVE = 0e+0_fp
                        XICE_ABOVE = 0e+0_fp
                     ENDIF
                     XPSC_ABOVE = XNAT_ABOVE + XICE_ABOVE
                     P_ABOVE    = GET_PCENTER(I,J,L+2) * 100.0e+0_fp
                     PSEDABOVE  = VNAT(L+2)*SEDSTEP*1.e+2_fp*
     &                          (GET_PEDGE(I,J,L+2)-GET_PEDGE(I,J,L+3))/
     &                          State_Met%BXHEIGHT(I,J,L+2)
!                     PSEDABOVE = g0 * (1.0d3*AIRMW) * P_ABOVE *
!     &                  VNAT(L+2)/(Rd * T(I,J,L+2))*SEDSTEP
                  ELSE
                     INVAIR_ABOVE = 0.0e+0_fp
                     XNO3_ABOVE = 0.0e+0_fp
                     XNAT_ABOVE = 0.0e+0_fp
                     XICE_ABOVE = 0.0e+0_fp
                     XPSC_ABOVE = 0.0e+0_fp
                     P_ABOVE = GET_PEDGE(I,J,L+2) * 100.0e+0_fp
                     PSEDABOVE = 0.0e+0_fp
                  ENDIF
                  ! Note reversal of pressure values as pressure falls with
                  ! height. If a value is positive, PSC particle substance
                  ! mixing ratios are increasing with height, and vice versa
                  BELOWGRAD = ABOVEGRAD
                  ABOVEGRAD = (XPSC_ABOVE-XPSC_0)/(P_0-P_ABOVE)

                  ! NB: Order of following conditionals is important! Could
                  ! theoretically do without localprofile as a variable, but
                  ! should be sure to check for the possible (but phenomenally
                  ! unlikely) condition that the gradients are both zero
                  IF ((ABOVEGRAD .eq. 0e+0_fp) .and. 
     &                            (BELOWGRAD .eq. 0e+0_fp)) THEN
                     ! Extremelely unlikely case, but use standard scheme
                     ! for flat vertical profile
                     LOCALPROFILE = 0
                  ELSEIF ((ABOVEGRAD .le. 0e+0_fp) .and.
     &                     (BELOWGRAD .le. 0e+0_fp)) THEN
                     ! Consistently decreasing with height - above a peak
                     LOCALPROFILE = -1
                  ELSEIF ((ABOVEGRAD .ge. 0e+0_fp) .and. 
     &                     (BELOWGRAD .ge. 0e+0_fp)) THEN
                     ! Consistently increasing with height - below a peak
                     LOCALPROFILE = +1
                  ELSE
                     ! Local minmum or maximum
                     LOCALPROFILE = 0
                  ENDIF

                  IF (LOCALPROFILE.eq.0) THEN
                     ! Standard sedimentation scheme
                     SEDQUANT = XPSC_0 * PSED0
                     SEDNO3   = XNO3_0 * PSED0
                     SEDICE   = XICE_0 * PSED0
                  ELSEIF (LOCALPROFILE.eq.1) THEN
                     ! Currently below a peak
                     ! Use gradient of bottom and centre boxes
                     IF (XPSC_0 .le. (0.5e+0_fp * XPSC_ABOVE)) THEN
                        ! Interpret as nearing peak; use lower gradient
                        SEDQUANT = (XPSC_0+XPSC_BELOW)*0.5e+0_fp*PSED0
                        SEDNO3   = (XNO3_BELOW+XNO3_0)*0.50*PSED0
                        SEDICE   = (XICE_BELOW+XICE_0)*0.50*PSED0
                     ELSE
                        SEDQUANT = (XPSC_ABOVE+XPSC_0)*0.50*PSED0
                        SEDNO3   = (XNO3_ABOVE+XNO3_0)*0.50*PSED0
                        SEDICE   = (XICE_ABOVE+XICE_0)*0.50*PSED0
                     ENDIF
                  ELSE
                     ! Above a peak
                     SEDQUANT = (XPSC_0+XPSC_BELOW)*0.5e+0_fp*PSED0
                     SEDNO3   = (XNO3_0+XNO3_BELOW)*0.5e+0_fp*PSED0
                     SEDICE   = (XICE_0+XICE_BELOW)*0.5e+0_fp*PSED0
                  ENDIF
                  ! Divide sedimenting quantity by the pressure difference
                  ! across the box being sedimented from
                  ! Note conversion from hPa to Pa for denominator, ie
                  ! multiply numerator by (1/100)
                  IF (L .ne. LLPAR) THEN
                     SEDQUANT = (1.e-2_fp) * SEDQUANT /
     &                  (GET_PEDGE(I,J,L+1)-GET_PEDGE(I,J,L+2))
                     SEDNO3 = (1.e-2_fp) * SEDNO3 / 
     &                  (GET_PEDGE(I,J,L+1)-GET_PEDGE(I,J,L+2))
                     SEDICE = (1.e-2_fp) * SEDICE / 
     &                  (GET_PEDGE(I,J,L+1)-GET_PEDGE(I,J,L+2))
                  ELSE
                     ! This shouldn't be possible?
                     CALL ERROR_STOP('Unknown sedimentation error',
     &                               'UCX_mod.F')
                  ENDIF
            
                  ! Apply limits so that sedimented quantity is:
                  ! Greater than or equal to zero
                  ! Less than or equal to the total available quantity

                  ! Note that we are calculating sedimentation using
                  ! the total ice and NAT, but are only actually
                  ! transporting the local NO3
                  SEDQUANT = MAX(0.0e+0_fp,SEDQUANT)
                  SEDNO3 = MAX(0.0e+0_fp,SEDNO3)
                  SEDICE = MAX(0.0e+0_fp,SEDICE)
                  
                  ! Convert v/v to kg/box
                  SEDNO3 = SEDNO3*TRACER_MW_G(IDTNIT)/INVAIR_0
                  SEDICE = SEDICE*ICEMW/INVAIR_0
                  
                  SEDNO3 = MIN(SEDNO3,STT(I,J,L+1,IDTNIT))
                  STT(I,J,L,IDTNIT) = STT(I,J,L,IDTNIT) + SEDNO3
                  STT(I,J,L+1,IDTNIT)=STT(I,J,L+1,IDTNIT)-SEDNO3

                  ! Settle the ice out too
                  SEDH2O = MIN(SEDICE,STT(I,J,L+1,IDTH2O))
                  STT(I,J,L,IDTH2O) = STT(I,J,L,IDTH2O) + SEDH2O
                  STT(I,J,L+1,IDTH2O)=STT(I,J,L+1,IDTH2O)-SEDH2O

                  ! Now correct aerosol totals
                  SEDNAT = SEDNO3*NATMW/TRACER_MW_G(IDTNIT)
                  SEDPSC = SEDNAT + SEDICE
                  SEDPSC = MIN(SEDPSC,KG_AER(I,J,L+1,I_SPA))
                  KG_AER(I,J,L,I_SPA) = KG_AER(I,J,L,I_SPA) + SEDPSC
                  KG_AER(I,J,L+1,I_SPA)=KG_AER(I,J,L+1,I_SPA)-SEDPSC
                  
                  ! Store for debug
                  IF (LPRT) THEN
                     FRACLOSS = SAFE_DIV(SEDH2O,STT(I,J,L+1,IDTH2O),
     &                                                       0.0e+0_fp)
                     IF (FRACLOSS .gt. MAXLOSS) THEN
                        MAXLOSS = FRACLOSS
                        D_ICE0 = XICE_0*ICEMW/INVAIR_0
                        D_ICE1 = D_ICE0 - SEDICE
                        D_H2O1 = STT(I,J,L+1,IDTH2O)
                        D_H2O0 = D_H2O1 + SEDH2O
                        D_NAT1 = STT(I,J,L+1,IDTNIT)*NATMW/
     &                                      TRACER_MW_G(IDTNIT)
                        D_NAT0 = D_NAT1 + SEDNAT
                        D_HNO3 = STT(I,J,L+1,IDTHNO3)
                        D_RMAX = RAD_AER(I,J,L+1,I_SPA)
                        D_VMAX = VNAT(L+1)
                        D_SEDQ = SEDQUANT
                        D_PSED = PSED0
                        D_XH2O = D_H2O0*INVAIR_0/ICEMW
                        D_XHNO3= D_HNO3*INVAIR_0/TRACER_MW_G(IDTHNO3)
                        D_XNAT = D_NAT0*INVAIR_0/NATMW
                        D_PPA  = P_0
                        D_TK   = TEMP
                     ENDIF
                  ENDIF

               ENDDO SED_LLOOP
            ENDDO SEDSTEPLOOP
         ENDIF ! VFALLMAX > 0

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      IF (LPRT) THEN
         CALL DEBUG_MSG( '### UCX: after SETTLE_STRAT_AER' )
         WRITE(6,'(a)') 'Max ice loss statistics:'
         WRITE(6,'(a,1x,2E12.4)') 'Source ice (before/after) (kg): ',
     &                         D_ICE0, D_ICE1
         WRITE(6,'(a,1x,2E12.4)') 'Source H2O (before/after) (kg): ',
     &                         D_H2O0, D_H2O1
         WRITE(6,'(a,1x,1E12.4)') 'Source HNO3 (total) (kg)      : ',
     &                         D_HNO3
         WRITE(6,'(a,1x,2E12.4)') 'Source NAT (before/after) (kg): ',
     &                         D_NAT0, D_NAT1
         WRITE(6,'(a,1x,3E12.4)') 'Source MR (H2O-HNO3-NAT) (v/v): ',
     &                         D_XH2O, D_XHNO3, D_XNAT
         WRITE(6,'(a,1x,2E12.4)') 'Source T (K) and P (Pa)       : ',
     &                         D_TK,   D_PPA
         WRITE(6,'(a,1x,2E12.4)') 'Particle radius (cm) & V (m/s): ',
     &                         D_RMAX, D_VMAX
         WRITE(6,'(a,1x,2E12.4)') 'Sedquant and PSed variables   : ',
     &                         D_SEDQ, D_PSED
      ENDIF

      END SUBROUTINE SETTLE_STRAT_AER
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_h2so4_gas
!
! !DESCRIPTION: Subroutine CALC\_H2SO4\_GAS calculates the fraction of strat.
!  SO4 aerosol which can be considered to be gaseous H2SO4.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_H2SO4_GAS( Input_Opt, State_Met, State_Chm )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_TROP
      USE CMN_GCTM_MOD,       ONLY : AIRMW
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACERID_MOD,       ONLY : IDTSO4, IDTSO4g
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY: 
!  11 Apr 2013 - S. D. Eastham - Initial version
!  07 Feb 2014 - R. Yantosca   - Cosmetic changes
!  14 Feb 2014 - R. Yantosca   - Changed parallel DO loop order to L-J-I,
!                                which is the most efficient configuration
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm
!                                objects via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST=.TRUE.
      REAL(fp),PARAMETER :: GF_THRESHOLD = 0.0e+0_fp
      REAL(fp),PARAMETER :: GF_RANGE     = 1.0e-8_fp
      REAL(fp),PARAMETER :: GF_DELTAHBYR = 10156.e+0_fp
      REAL(fp),PARAMETER :: GF_T0        = 360.e+0_fp
      REAL(fp),PARAMETER :: GF_TC        = 905.e+0_fp
      REAL(fp),SAVE      :: GF_LOGP0
      REAL(fp),SAVE      :: GF_BFACTOR
      REAL(fp),SAVE      :: GF_ATMCONV
      REAL(fp),SAVE      :: GF_INVT0
      REAL(fp)           :: GF_INVT,GF_LOGPSULFATE,GF_CFACTOR
      REAL(fp)           :: GF_AFACTOR
      REAL(fp)           :: GF_PP,GF_PVAP,GF_DIFF

      INTEGER            :: I, J, L
      REAL(fp)           :: PCENTER, TCENTER, H2SO4SUM, INVAIR
      REAL(fp)           :: H2SO4KG, LIQFRAC, GASFRAC
      LOGICAL            :: ISLIQ, ISGAS

      ! Local variables for quantities from Input_Opt
      LOGICAL            :: LPRT
      INTEGER            :: N_TRACERS
      REAL(fp)           :: TRACER_MW_G(Input_Opt%N_TRACERS)

      ! Pointers
      REAL(fp), POINTER    :: STT(:,:,:,:)

      !=================================================================
      ! CALC_H2SO4_GAS begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      LPRT        = Input_Opt%LPRT
      N_TRACERS   = Input_Opt%N_TRACERS
      TRACER_MW_G = Input_Opt%TRACER_MW_G(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg]
      STT => State_Chm%Tracers

      IF (FIRST) THEN
         FIRST = .FALSE.
         ! Calculate H2SO4 gas phase prefactors
         GF_INVT0 = 1.e+0_fp/GF_T0
         GF_LOGP0 = (-1.e+0_fp*GF_DELTAHBYR*GF_INVT0) + 16.259e+0_fp
         GF_BFACTOR = 0.38e+0_fp/(GF_TC - GF_T0)
         GF_ATMCONV = LOG(101325.e+0_fp)
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J,          I,              L       )
!$OMP+PRIVATE( PCENTER,    TCENTER,        INVAIR  )
!$OMP+PRIVATE( H2SO4SUM,   GF_PP,          GF_INVT )
!$OMP+PRIVATE( GF_CFACTOR, GF_LOGPSULFATE, GF_PVAP )
!$OMP+PRIVATE( GF_DIFF,    H2SO4KG, LIQFRAC,GASFRAC)
!$OMP+PRIVATE( ISLIQ, ISGAS                        )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ! Only interested in low-pressure boxes
         PCENTER = GET_PCENTER(I,J,L)
         TCENTER = State_Met%T(I,J,L)
         INVAIR  = AIRMW / State_Met%AD(I,J,L)
         H2SO4KG = (STT(I,J,L,IDTSO4)+STT(I,J,L,IDTSO4g))
         IF (PCENTER.ge.1.e+2_fp) THEN
            ISLIQ = .TRUE.
            ISGAS = .FALSE.
         ELSEIF (ITS_IN_THE_TROP(I,J,L,State_Met)) THEN
            ! Don't want to interfere with tropospheric aerosols
            ISLIQ = .TRUE.
            ISGAS = .FALSE.
         ELSE
            H2SO4SUM = H2SO4KG*INVAIR/TRACER_MW_G(IDTSO4  )
            ! Use approximation from Kumala (1990)
            GF_PP = H2SO4SUM*PCENTER
            GF_INVT = 1./TCENTER
            GF_CFACTOR = 1.e+0_fp+(LOG(GF_T0*GF_INVT))-(GF_T0*GF_INVT)
            GF_LOGPSULFATE = GF_LOGP0 + (GF_DELTAHBYR*(GF_INVT0 -
     &               GF_INVT + (GF_BFACTOR*GF_CFACTOR)))
            GF_LOGPSULFATE = GF_LOGPSULFATE + GF_ATMCONV
            GF_PVAP = 1.e-2_fp * EXP(GF_LOGPSULFATE)
            GF_DIFF = (GF_PVAP+GF_THRESHOLD) - GF_PP
            IF (GF_DIFF .lt. 0) THEN
               ISLIQ = .TRUE.
               ISGAS = .FALSE.
            ELSEIF (GF_DIFF .lt. GF_RANGE) THEN
               GASFRAC = (GF_DIFF/GF_RANGE)
               LIQFRAC = 1.e+0_fp-GASFRAC
               ISLIQ = .FALSE.
               ISGAS = .FALSE.
            ELSE
               ISLIQ = .FALSE.
               ISGAS = .TRUE.
            ENDIF
         ENDIF
         ! Repartition
         IF (ISLIQ) THEN
            AERFRAC(I,J,L,1) = 1.0e+0_fp
            STT(I,J,L,IDTSO4)  = H2SO4KG
            STT(I,J,L,IDTSO4g) = 0.0e+0_fp
         ELSEIF (ISGAS) THEN
            AERFRAC(I,J,L,1) = 0.0e+0_fp
            STT(I,J,L,IDTSO4g) = H2SO4KG
            STT(I,J,L,IDTSO4)  = 0.0e+0_fp
         ELSE
            AERFRAC(I,J,L,1) = 1.0e+0_fp
            STT(I,J,L,IDTSO4g) = H2SO4KG*GASFRAC
            STT(I,J,L,IDTSO4)  = H2SO4KG*LIQFRAC
         ENDIF
      ENDDO
      ENDDO
      ENDDO

!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      IF ( LPRT ) CALL DEBUG_MSG( '### UCX: H2SO4 partitioned' )

      END SUBROUTINE CALC_H2SO4_GAS
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_strat_aer
!
! !DESCRIPTION: Subroutine CALC\_STRAT\_AER calculates aerosol properties    
!  stratosphere using the thermodynamic parameterization described in
!  Kirner et al. (`Simulation of polar stratospheric clouds in the
!  chemistry-climate-model EMAC via the submodel PSC', Geosci. Mod. Dev.,
!  4, 169-182, 2011).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_STRAT_AER( am_I_Root, Input_Opt,
     &                           State_Met, State_Chm, RC )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_STRATMESO
      USE CMN_GCTM_MOD,       ONLY : AIRMW, PI
      USE COMODE_MOD,         ONLY : CSPEC,JLOP
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE ERROR_MOD,          ONLY : IS_SAFE_DIV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      !USE GRID_MOD,           ONLY : GET_YEDGE
      USE GRID_MOD,           ONLY : GET_YMID
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TRACERID_MOD,       ONLY : IDTSO4, IDTH2O
      USE TRACERID_MOD,       ONLY : IDTHNO3, IDTNIT, IDTBrNO3, IDTHBr
      USE TRACERID_MOD,       ONLY : IDTHOBr, IDTHCl, IDTClNO3, IDTHOCl
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
! !OUTPUT PARAMETERS:
!      
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  13 Apr 2013 - S. D. Eastham - Initial version
!  06 Feb 2014 - R. Yantosca   - Add missing variables to OpenMP loop
!  14 Feb 2014 - R. Yantosca   - Make the L-loop the outermost DO loop
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm 
!                                objects via the arg list
!  06 Nov 2014 - R. Yantosca   - Now use State_Met%CLDF(I,J,L)
!  08 Apr 2015 - R. Yantosca   - Remove call to READ_PSC_FILE, this is
!                                now done from DO_CHEMISTRY (chemistry_mod.F)  
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Limits on PSC formation
      REAL(fp),PARAMETER      :: PSC_MAXLAT =  45.0e+0_fp
      REAL(fp),PARAMETER      :: PSC_MINLAT = -55.0e+0_fp
      REAL(fp),PARAMETER      :: PSC_PMAX   =  18.0e+3_fp
      REAL(fp),PARAMETER      :: PSC_PMIN   =   5.0e+2_fp

      ! Allow PSC formation outsize Kirner limits?
      LOGICAL,PARAMETER     :: PSC_FULL   =.FALSE.

      ! Saturation and nucleation properties of HNO3
      REAL(fp), PARAMETER     :: TSATHNO3_A = -2.7836e+0_fp
      REAL(fp), PARAMETER     :: TSATHNO3_B = -0.00088e+0_fp
      REAL(fp), PARAMETER     :: TSATHNO3_C = 38.9855e+0_fp
      REAL(fp), PARAMETER     :: TSATHNO3_D = -11397.0e+0_fp
      REAL(fp), PARAMETER     :: TSATHNO3_E = 0.009179e+0_fp

      ! Saturation and nucleation properties of water
      REAL(fp), PARAMETER     :: TSATH2O_A  = -2663.5e+0_fp
      REAL(fp), PARAMETER     :: TSATH2O_B  = 12.537e+0_fp

      ! Peak pressure at which NAT can form homogeneously
      REAL(fp), PARAMETER     :: P_MAXNAT = 1.40e+4_fp ! Pa

      ! Maximum temperature for PSC formation (K)
      REAL(fp), PARAMETER     :: T_MAX = 215.0e+0_fp

      ! Limits on NAT/ice formation
      REAL(fp), PARAMETER     :: MIN_RAD = 1.0e-7_fp ! m
      REAL(fp), PARAMETER     :: MAX_NDENS=42.0e+3_fp ! #/m3

      ! Local conditions
      REAL(fp)                :: TCENTER, PCENTER, PCENTER_PA, DENAIR
      REAL(fp)                :: INVAIR, TINV, TOFFSET

      ! Gridbox mixing ratios and partial pressures
      REAL(fp)                :: HNO3SUM, H2OSUM
      REAL(fp)                :: HNO3PP,  H2OPP
      REAL(fp)                :: PSATHNO3_SUPERCOOL
      REAL(fp)                :: PSATH2O_SUPERSAT
      REAL(fp)                :: PSATHNO3, PSATH2O
      REAL(fp)                :: H2SO4SUM
      REAL(fp)                :: ClNO3SUM, HClSUM, HOClSUM
      REAL(fp)                :: BrNO3SUM, HBrSUM, HOBrSUM

      ! Gridbox aerosol and phase data
      REAL(fp)                :: HNO3_BOX_G, HNO3_BOX_L, HNO3_BOX_S
      REAL(fp)                :: H2O_BOX_G,  H2O_BOX_L,  H2O_BOX_S
      REAL(fp)                :: H2SO4_BOX_G,H2SO4_BOX_L
      REAL(fp)                :: HCl_BOX_G,  HCl_BOX_L
      REAL(fp)                :: HOCl_BOX_G, HOCl_BOX_L
      REAL(fp)                :: HBr_BOX_G,  HBr_BOX_L
      REAL(fp)                :: HOBr_BOX_G, HOBr_BOX_L
      REAL(fp)                :: HNO3GASFRAC, HClGASFRAC, HOClGASFRAC
      REAL(fp)                :: HBrGASFRAC, HOBrGASFRAC
      REAL(fp)                :: VOL_NAT, VOL_ICE, VOL_SLA, VOL_TOT
      REAL(fp)                :: RAD_AER_BOX,RHO_AER_BOX
      REAL(fp)                :: KG_AER_BOX,NDENS_AER_BOX,SAD_AER_BOX
      REAL(fp)                :: KG_NAT, KG_ICE, KG_NO3

      ! SLA weight fractions
      REAL(fp)                :: W_H2SO4, W_H2O, W_HNO3
      REAL(fp)                :: W_HCl, W_HOCl, W_HBr, W_HOBr

      ! Reaction prefactors
      REAL(fp)                :: KHET_COMMON
      REAL(fp)                :: KHET_SPECIFIC

      ! Grid box location
      REAL(fp)                :: BOX_LAT_S, BOX_LAT_N, BOX_LAT
      LOGICAL                 :: IS_VALID, IS_POLAR, IS_STRAT

      ! Local properties
      REAL(fp), DIMENSION(11) :: GAMMA_BOX
      INTEGER               :: STATE_LOCAL

      ! Loop variables
      INTEGER               :: I, J, L, K, JLOOP

      ! SDE 2013-07-07: Debug variables
      REAL(fp)                 :: D_XH2O, D_PPA,  D_TK,   D_XHNO3,D_XNAT
      REAL(fp)                 :: D_MAXN, D_FRACN,D_PP,   D_PS

      ! Local variables for quantities from Input_Opt
      LOGICAL                :: prtDebug
      LOGICAL                :: LPRT
      LOGICAL                :: LHOMNUCNAT
      LOGICAL                :: LSOLIDPSC
      LOGICAL                :: LACTIVEH2O
      LOGICAL                :: LWACCMAER
      INTEGER                :: N_TRACERS
      REAL(fp)                 :: TRACER_MW_G(Input_Opt%N_TRACERS)

      ! Pointers
      REAL(fp), POINTER        :: STT(:,:,:,:)

      !=================================================================
      ! CALC_STRAT_AER begins here!
      !=================================================================

      ! Assume success
      RC       = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT
      LPRT        = Input_Opt%LPRT
      LHOMNUCNAT  = Input_Opt%LHOMNUCNAT
      LSOLIDPSC   = Input_Opt%LSOLIDPSC
      LACTIVEH2O  = Input_Opt%LACTIVEH2O
      N_TRACERS   = Input_Opt%N_TRACERS
      TRACER_MW_G = Input_Opt%TRACER_MW_G(1:N_TRACERS)
      LWACCMAER   = Input_Opt%LWACCMAER

      ! Initialize GEOS-Chem tracer array [kg]
      STT => State_Chm%Tracers

      ! Do we have to print debug output?
      prtDebug = ( LPRT .and. am_I_Root )

      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### UCX: start CALC_STRAT_AER' )
      ENDIF

      ! Partition H2SO4 before proceeding
      CALL CALC_H2SO4_GAS( Input_Opt, State_Met, State_Chm )

      D_XH2O = 0e+0_fp
      D_PPA  = 0e+0_fp
      D_TK   = 0e+0_fp
      D_XHNO3= 0e+0_fp
      D_XNAT = 0e+0_fp
      D_MAXN = 0e+0_fp
      D_PP   = 0e+0_fp
      D_PS   = 0e+0_fp
      
      ! Loop over latitude boxes first
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,            J,                  L             )
!$OMP+PRIVATE( K,            IS_POLAR,           IS_STRAT      )
!$OMP+PRIVATE( PCENTER,      PCENTER_PA,         DENAIR        )
!$OMP+PRIVATE( INVAIR,       PSATHNO3,           PSATH2O       )
!$OMP+PRIVATE( JLOOP,        VOL_NAT,            VOL_ICE       )
!$OMP+PRIVATE( VOL_SLA,      PSATHNO3_SUPERCOOL, TCENTER       )
!$OMP+PRIVATE( TINV,         IS_VALID                          )
!$OMP+PRIVATE( RAD_AER_BOX,  RHO_AER_BOX                       )
!$OMP+PRIVATE( KG_AER_BOX,   NDENS_AER_BOX,      SAD_AER_BOX   )
!$OMP+PRIVATE( KG_NAT,       KG_ICE,             KG_NO3        )
!$OMP+PRIVATE( GAMMA_BOX,    PSATH2O_SUPERSAT,   H2OSUM        )
!$OMP+PRIVATE( H2OPP,        H2O_BOX_S,          H2O_BOX_L     )
!$OMP+PRIVATE( H2O_BOX_G,    H2SO4SUM,           HNO3SUM       )
!$OMP+PRIVATE( HNO3PP,       HNO3_BOX_S,         HNO3_BOX_L    )
!$OMP+PRIVATE( HNO3_BOX_G,   BrNO3SUM,           HBrSUM        )
!$OMP+PRIVATE( HOBrSUM,      ClNO3SUM,           HClSUM        )
!$OMP+PRIVATE( HOClSUM,      STATE_LOCAL,        HBrGASFRAC    )
!$OMP+PRIVATE( HOBrGASFRAC,  HNO3GASFRAC,        HClGASFRAC    )
!$OMP+PRIVATE( HOClGASFRAC,  TOFFSET,            W_H2SO4       )
!$OMP+PRIVATE( W_H2O,        W_HCl,              W_HOCl        )
!$OMP+PRIVATE( W_HBr,        W_HOBr,             W_HNO3        )
!$OMP+PRIVATE( HCl_BOX_G,    HCl_BOX_L,          HOCl_BOX_G    )
!$OMP+PRIVATE( HOCl_BOX_L,   H2SO4_BOX_G,        HBr_BOX_G     )
!$OMP+PRIVATE( HBr_BOX_L,    HOBr_BOX_G,         HOBr_BOX_L    )
!$OMP+PRIVATE( H2SO4_BOX_L,  KHET_COMMON,        KHET_SPECIFIC )
!$OMP+PRIVATE( D_FRACN,      D_MAXN,             D_XH2O        )
!$OMP+PRIVATE( D_XHNO3,      D_XNAT,             D_PPA         )
!$OMP+PRIVATE( D_PP,         D_PS,               D_TK          )
!$OMP+PRIVATE( VOL_TOT,      BOX_LAT                           )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
         DO I = 1, IIPAR

               ! Now do IS_POLAR check for every grid box separately
               ! and using mid-point instead of edges. This is to be
               ! applicable to cube-sphere systems (ckeller, 1/16/15).
               BOX_LAT  = GET_YMID(I,J,1)
               IS_POLAR = ( (BOX_LAT <= PSC_MINLAT) .OR.
     &                      (BOX_LAT >= PSC_MAXLAT)     )

               ! Get local conditions
               PCENTER    = GET_PCENTER(I,J,L) 
               PCENTER_PA = PCENTER*1.e+2_fp
               TCENTER    = State_Met%T(I,J,L)
               TINV       = 1.e+0_fp/TCENTER

               ! Apply other limits from Kirner et al.
               IS_STRAT = ITS_IN_THE_STRATMESO(I,J,L,State_Met)
               IS_VALID = (IS_POLAR).and.(IS_STRAT).and.       
     &                        (.not.((PCENTER_PA.lt.PSC_PMIN).or.
     &                               (PCENTER_PA.gt.PSC_PMAX)))
               IS_VALID = (IS_VALID.or.PSC_FULL)

               ! Initialize variables
               VOL_NAT      = 0.0e+0_fp
               KG_NAT       = 0.0e+0_fp
               VOL_ICE      = 0.0e+0_fp
               KG_ICE       = 0.0e+0_fp
 
               RAD_AER_BOX  = 0.0e+0_fp
               RHO_AER_BOX  = 0.0e+0_fp
               KG_AER_BOX   = 0.0e+0_fp
               NDENS_AER_BOX= 0.0e+0_fp

               H2O_BOX_S    = 0.0e+0_fp
               H2O_BOX_L    = 0.0e+0_fp
               H2O_BOX_G    = 0.0e+0_fp
               HNO3_BOX_S   = 0.0e+0_fp
               HNO3_BOX_L   = 0.0e+0_fp
               HNO3_BOX_G   = 0.0e+0_fp
               H2SO4_BOX_L  = 0.0e+0_fp
               H2SO4_BOX_G  = 0.0e+0_fp
               HCl_BOX_L    = 0.0e+0_fp
               HCl_BOX_G    = 0.0e+0_fp
               HOCl_BOX_L   = 0.0e+0_fp
               HOCl_BOX_G   = 0.0e+0_fp
               HBr_BOX_L    = 0.0e+0_fp
               HBr_BOX_G    = 0.0e+0_fp
               HOBr_BOX_L   = 0.0e+0_fp
               HOBr_BOX_G   = 0.0e+0_fp

               !STATE_LOCAL  = STATE_PSC(I,J,L)
               STATE_LOCAL  = NINT(STATE_PSC(I,J,L))

               ! Calculate local air density
               DENAIR = AVOGADRO * PCENTER_PA / (TCENTER*R_UNIV)
               INVAIR = AIRMW / State_Met%AD(I,J,L)

               ! Get available NO3 mass
               KG_NO3 = (STT(I,J,L,IDTHNO3)*
     &                   TRACER_MW_G(IDTNIT)/TRACER_MW_G(IDTHNO3)) + 
     &                   STT(I,J,L,IDTNIT)

               ! Calculate HNO3 mixing ratio
               HNO3SUM = STT(I,J,L,IDTHNO3)*INVAIR/
     &                    TRACER_MW_G(IDTHNO3)
               HNO3SUM = HNO3SUM + STT(I,J,L,IDTNIT)*INVAIR/
     &                    TRACER_MW_G(IDTNIT)

               H2OSUM = STT(I,J,L,IDTH2O)*INVAIR/
     &                   TRACER_MW_G(IDTH2O)

               ! Calculate partial pressures (Pa)
               HNO3PP = PCENTER_PA * HNO3SUM
               H2OPP  = PCENTER_PA * H2OSUM

               IF (.not.IS_VALID) THEN
                  ! No PSCs (SSA only)
                  STATE_LOCAL = 0
               ELSEIF (TCENTER.gt.T_MAX) THEN
                  ! STS/SSA only
                  STATE_LOCAL = 1
               ELSE
                  ! Calculate saturation pressures
                  PSATHNO3 = (10.e+0_fp**(((TSATHNO3_A+(TCENTER*
     &            TSATHNO3_B))*LOG10(H2OPP*760.0e+0_fp/101325.0e+0_fp))+
     &            (TSATHNO3_C+(TSATHNO3_D*TINV)+(TSATHNO3_E*
     &            TCENTER))))*101325.0e+0_fp/760.0e+0_fp

                  PSATH2O = 10.0e+0_fp**((TSATH2O_A*TINV)+TSATH2O_B)

                  ! Supersaturation requirement for ice
                  PSATH2O_SUPERSAT = PSATH2O * P_ICE_SUPERSAT

                  ! If homogeneous NAT nucleation allowed, calculate
                  ! threshold saturation pressure
                  IF (LHOMNUCNAT) THEN
                     ! Calculate as if temperature is (T+T_NAT_SUPERCOOL)
                     TOFFSET = TCENTER + T_NAT_SUPERCOOL
                     PSATHNO3_SUPERCOOL = 10.e+0_fp**(((TSATHNO3_A+
     &                  (TOFFSET*TSATHNO3_B))*LOG10(H2OPP*760.0e+0_fp/
     &                  101325.0e+0_fp))+
     &                  (TSATHNO3_C+(TSATHNO3_D/TOFFSET)+(TSATHNO3_E*
     &                  TOFFSET)))*101325.0e+0_fp/760.0e+0_fp
                  ELSE
                     ! Make homogeneous nucleation impossible
                     TOFFSET = 280.0e+0_fp
                     PSATHNO3_SUPERCOOL = HNO3PP + 1e+0_fp
                  ENDIF
                  
                  IF (LACTIVEH2O) THEN
                     ! Only interested in sign
                     IF (STATE_LOCAL .gt. 1) THEN
                        H2O_BOX_S = H2OPP-PSATH2O
                     ELSE
                        H2O_BOX_S = H2OPP-PSATH2O_SUPERSAT
                     ENDIF
                  ELSE
                     ! Use local ice mass ratio from GEOS-5 data
                     ! Note that we are using this only for its sign!
                     H2O_BOX_S = State_Met%QI(I,J,L)
                  ENDIF

                  ! Ice exists/possible?
                  IF (H2O_BOX_S.gt.TINY(0e+0_fp)) THEN
                     STATE_LOCAL = 3
                     HNO3_BOX_S = HNO3PP-PSATHNO3
                  ELSE
                     ! If ice not possible could still have NAT
                     H2O_BOX_S = 0e+0_fp
                     ! 1. Homogeneous nucleation
                     IF ((LHOMNUCNAT).and.(STATE_LOCAL.eq.1)) THEN
                        HNO3_BOX_S = HNO3PP-PSATHNO3_SUPERCOOL
                     ENDIF
                     ! 2. Box formerly contained ice or NAT
                     IF (STATE_LOCAL.eq.2) THEN
                        HNO3_BOX_S = HNO3PP-PSATHNO3
                     ENDIF
                     IF (HNO3_BOX_S.gt.TINY(1e+0_fp)) THEN
                        STATE_LOCAL = 2
                     ELSE
                        STATE_LOCAL = 1
                     ENDIF
                  ENDIF
               ENDIF

               ! Store state
               !STATE_PSC(I,J,L) = STATE_LOCAL
               STATE_PSC(I,J,L) = REAL(STATE_LOCAL,f4)

               ! Only continue if we want online solid PSCs
               IF (LSOLIDPSC) THEN

                  IF (STATE_LOCAL.eq.3) THEN
                     ! Form ice PSCs
                     IF (LACTIVEH2O) THEN
                        H2O_BOX_S = (H2OPP-PSATH2O)/PCENTER_PA
                        H2O_BOX_S = MAX(0e+0_fp,H2O_BOX_S)
                        KG_ICE    = H2O_BOX_S*ICEMW*
     &                              State_Met%AD(I,J,L)/AIRMW
                     ELSE
                        H2O_BOX_S = State_Met%QI(I,J,L)   *
     &                              State_Met%CLDF(I,J,L) *
     &                              AIRMW / ICEMW
                        KG_ICE    = State_Met%QI(I,J,L)   *
     &                              State_Met%CLDF(I,J,L) *
     &                              State_Met%AD(I,J,L)
                     ENDIF
                     VOL_ICE = H2O_BOX_S * DENAIR * (1.e-3_fp) * ICEMW /
     &                             (DENSICE*AVOGADRO) ! m3 ice/m3 air 
                  ELSE
                     VOL_ICE = 0e+0_fp
                     H2O_BOX_S = 0e+0_fp
                     KG_ICE = 0e+0_fp
                  ENDIF

                  ! Calculate NAT if relevant
                  IF ((HNO3_BOX_S.gt.TINY(1e+0_fp)).and.
     &                               (STATE_LOCAL.ge.2)) THEN
                     HNO3_BOX_S = (HNO3PP-PSATHNO3)/PCENTER_PA
                     HNO3_BOX_S = MAX(0e+0_fp,HNO3_BOX_S)

                     ! If using WACCM, a maximum of 20% of the total
                     ! NO3 in the grid cell can go to NAT
                     IF (Input_Opt%LWACCMAER) THEN
                        HNO3_BOX_S = MIN(HNO3_BOX_S,0.20e+0_fp*HNO3SUM)
                     ENDIF

                     ! Calculate m3 NAT/m3 air
                     ! HNO3_BOX_S is the number of moles of HNO3
                     ! which will be frozen into HNO3.3H2O (NAT)
                     ! Therefore volume calculation must be done
                     ! with care!
                     VOL_NAT = HNO3_BOX_S * DENAIR * (1.e-3_fp) *
     &                         NATMW / (DENSNAT*AVOGADRO)
                     KG_NAT  = HNO3_BOX_S * NATMW *
     &                         State_Met%AD(I,J,L) / AIRMW

                     ! Debug
                     IF (HNO3SUM.gt.TINY(0e+0_fp)) THEN
                        D_FRACN = HNO3_BOX_S/HNO3SUM
                     ELSE
                        D_FRACN = 0e+0_fp
                     ENDIF
!
!                     IF (D_FRACN.gt.D_MAXN) THEN
!                        D_MAXN = D_FRACN
                     IF (HNO3_BOX_S.gt.D_XNAT) THEN
                        D_MAXN = D_FRACN
                        D_XH2O = H2OSUM
                        D_XHNO3= HNO3SUM
                        D_XNAT = HNO3_BOX_S
                        D_PPA  = PCENTER_PA
                        D_PP   = HNO3PP
                        D_PS   = PSATHNO3
                        D_TK   = TCENTER
                     ENDIF
                     
                  ELSE
                     HNO3_BOX_S = 0e+0_fp
                     VOL_NAT = 0e+0_fp
                     KG_NAT = 0e+0_fp
                  ENDIF

                  ! Calculate particle properties
                  IF (STATE_LOCAL.lt.2) THEN
                     ! Zero all!
                     KG_AER_BOX   = 0e+0_fp
                     RAD_AER_BOX  = 0e+0_fp
                     RHO_AER_BOX  = DENSICE
                     NDENS_AER_BOX= 0e+0_fp
                     SAD_AER_BOX  = 0e+0_fp
                  ELSE
                     VOL_TOT = VOL_NAT + VOL_ICE
                     KG_AER_BOX = KG_NAT + KG_ICE
                     IF (LWACCMAER) THEN
                        ! Number density is a given parameter
                        NDENS_AER_BOX = NWACCM(I_SPA)
                        ! Derive effective radius from volume and
                        ! known properties of log-normal aerosol
                        ! distributions
                        RAD_AER_BOX = (VOL_TOT*LN_VCONV(I_SPA))
     &                                     **(1.0e+0_fp/3.0e+0_fp)
                        SAD_AER_BOX = 1.0e-2_fp*(LN_ACONV(I_SPA)*
     &                                 RAD_AER_BOX*RAD_AER_BOX)
                     ELSE
                        ! Monodisperse
                        RAD_AER_BOX = MIN_RAD
                        NDENS_AER_BOX = (3.0e+0_fp*(VOL_TOT)/
     &                   (4.0e+0_fp*PI*(RAD_AER_BOX**3.0e+0_fp)))
                        IF (NDENS_AER_BOX.gt.MAX_NDENS) THEN
                         NDENS_AER_BOX = MAX_NDENS
                         RAD_AER_BOX = (3.0e+0_fp*(VOL_TOT)/
     &                    (4.0e+0_fp*PI*MAX_NDENS))**(1.e+0_fp/3.e+0_fp)
                        ENDIF
                        SAD_AER_BOX = 4.0e-2_fp * RAD_AER_BOX *
     &                     RAD_AER_BOX * NDENS_AER_BOX * PI
                     ENDIF

                     ! Prevent div-zero (ckeller, 1/16/15)
                     IF ( VOL_TOT > 0.0_fp ) THEN
                        RHO_AER_BOX = 
     &                     ((VOL_ICE*DENSICE)+(VOL_NAT*DENSNAT))/VOL_TOT
                     ELSE
                        RHO_AER_BOX = DENSICE ! Is that correct?
                     ENDIF

                  ENDIF
               ELSE
                  ! Solid PSCs not active
                  RAD_AER_BOX = 0e+0_fp
                  RHO_AER_BOX = 1000e+0_fp
                  KG_AER_BOX = 0e+0_fp
                  NDENS_AER_BOX = 0e+0_fp
                  SAD_AER_BOX = 0e+0_fp
                  HNO3_BOX_S = 0e+0_fp
                  H2O_BOX_S = 0e+0_fp
               ENDIF

               ! Store in outer arrays
               RAD_AER(I,J,L,I_SPA)  = RAD_AER_BOX*1.e+2_fp ! cm
               RHO_AER(I,J,L,I_SPA)  = RHO_AER_BOX ! kg/m3
               KG_AER(I,J,L,I_SPA)   = KG_AER_BOX  ! kg
               NDENS_AER(I,J,L,I_SPA)= NDENS_AER_BOX!#/m3
               SAD_AER(I,J,L,I_SPA)  = SAD_AER_BOX ! cm2/cm3

               ! Repartition NIT and HNO3 in strat/meso
               IF (LSOLIDPSC.and.IS_STRAT) THEN
                  ! Convert NAT from kg NAT to kg NO3
                  STT(I,J,L,IDTNIT) = KG_NAT*TRACER_MW_G(IDTNIT)/NATMW
                  ! Remove (kg NO3 as NAT) from total kg NO3
                  ! then convert to kg HNO3
                  STT(I,J,L,IDTHNO3) = (KG_NO3-STT(I,J,L,IDTNIT))
     &                  *TRACER_MW_G(IDTHNO3)/TRACER_MW_G(IDTNIT) 
               ENDIF

               ! Now start liquid aerosol consideration
               ! Start by assuming all non-solid H2O/HNO3 is gaseous
               HNO3_BOX_G = HNO3SUM - HNO3_BOX_S
               HNO3_BOX_L = 0e+0_fp
               H2O_BOX_G = H2OSUM - H2O_BOX_S
               H2O_BOX_L = 0e+0_fp

               ! Calculate mixing ratios of other relevant species
               H2SO4SUM = STT(I,J,L,IDTSO4  )*INVAIR/
     &                  TRACER_MW_G(IDTSO4  )
               BrNO3SUM = STT(I,J,L,IDTBrNO3)*INVAIR/
     &                  TRACER_MW_G(IDTBrNO3)
               ClNO3SUM = STT(I,J,L,IDTClNO3)*INVAIR/
     &                  TRACER_MW_G(IDTClNO3)
               HOClSUM  = STT(I,J,L,IDTHOCl )*INVAIR/
     &                  TRACER_MW_G(IDTHOCl )
               HClSUM   = STT(I,J,L,IDTHCl  )*INVAIR/
     &                  TRACER_MW_G(IDTHCl  )
               HOBrSUM  = STT(I,J,L,IDTHOBr )*INVAIR/
     &                  TRACER_MW_G(IDTHOBr )
               HBrSUM   = STT(I,J,L,IDTHBr  )*INVAIR/
     &                  TRACER_MW_G(IDTHBr  )

               ! H2SO4 gas fraction calculated earlier throughout grid
               ! Consider gaseoues H2SO4 to be unavailable for SLA
               H2SO4_BOX_L = H2SO4SUM * AERFRAC(I,J,L,1)
               H2SO4_BOX_G = H2SO4SUM - H2SO4_BOX_L

               ! Zero local properties
               RHO_AER_BOX   = 1000e+0_fp
               RAD_AER_BOX   = 0e+0_fp 
               KG_AER_BOX    = 0e+0_fp
               NDENS_AER_BOX = 0e+0_fp
               SAD_AER_BOX   = 0e+0_fp
               VOL_SLA       = 0e+0_fp
               W_H2O         = 0e+0_fp
               W_H2SO4       = 1e+0_fp

               IF (.not.IS_STRAT) THEN
                  ! Use JPL 10-06/Oslo CTM data, where available,
                  ! for conventional sulfates/H2SO4
                  GAMMA_BOX(1)  = 0.1e+0_fp
                  GAMMA_BOX(2)  = 0.0e+0_fp
                  GAMMA_BOX(3)  = 0.0e+0_fp
                  GAMMA_BOX(4)  = 0.0e+0_fp
                  GAMMA_BOX(5)  = 0.3e+0_fp
                  GAMMA_BOX(6)  = 0.4e+0_fp
                  GAMMA_BOX(7)  = 0.9e+0_fp
                  GAMMA_BOX(8)  = 0.0e+0_fp
                  GAMMA_BOX(9)  = 0.0e+0_fp
                  GAMMA_BOX(10) = 0.2e+0_fp
                  GAMMA_BOX(11) = 0.0e+0_fp
               ELSEIF (H2SO4_BOX_L.lt.1e-15_fp) THEN
                  ! No aerosol to speak of
                  DO K=1,11
                     GAMMA_BOX(K) = 0.0e+0_fp
                  ENDDO
               ELSE
                  IF (STATE_LOCAL.eq.0) THEN
                     ! Allow binary H2SO4.nH2O only
                     CALL TERNARY (PCENTER,TCENTER,H2OSUM,H2SO4_BOX_L,
     &                  0.e+0_fp   ,HClSUM,HOClSUM,HBrSUM,HOBrSUM,
     &                  W_H2SO4,W_H2O,W_HNO3,W_HCl,W_HOCl,W_HBr,W_HOBr,
     &                  HNO3GASFRAC,HClGASFRAC,HOClGASFRAC,
     &                  HBrGASFRAC,HOBrGASFRAC,VOL_SLA,RHO_AER_BOX)

                        ! For safety's sake, zero out HNO3 uptake
                        HNO3GASFRAC = 1.e+0_fp
                        W_H2O = W_H2O + W_HNO3
                        W_HNO3 = 0.e+0_fp
                        HNO3_BOX_G = HNO3SUM - HNO3_BOX_S
                        HNO3_BOX_L = 0.e+0_fp
                  ELSE
                     ! As per Buchholz, use only non-NAT HNO3 for STS
                     ! This is also the WACCM approach
                     HNO3_BOX_G = HNO3SUM - HNO3_BOX_S
                     CALL TERNARY (PCENTER,TCENTER,H2OSUM,H2SO4_BOX_L,
     &                  HNO3_BOX_G,HClSUM,HOClSUM,HBrSUM,HOBrSUM,
     &                  W_H2SO4,W_H2O,W_HNO3,W_HCl,W_HOCl,W_HBr,W_HOBr,
     &                  HNO3GASFRAC,HClGASFRAC,HOClGASFRAC,
     &                  HBrGASFRAC,HOBrGASFRAC,VOL_SLA,RHO_AER_BOX)
        
                        ! Partition HNO3 here for safety
                        HNO3_BOX_G = HNO3_BOX_G*HNO3GASFRAC
                        HNO3_BOX_L = HNO3SUM - (HNO3_BOX_G+HNO3_BOX_S)
                  ENDIF

                  ! Partition minor species
                  HCl_BOX_G  = HClSUM *HClGASFRAC
                  HCl_BOX_L  = HClSUM -HCl_BOX_G
                  HOCl_BOX_G = HOClSUM*HOClGASFRAC
                  HOCl_BOX_L = HOClSUM-HOCl_BOX_G
                  HBr_BOX_G  = HBrSUM *HBrGASFRAC
                  HBr_BOX_L  = HBrSUM -HBr_BOX_G
                  HOBr_BOX_G = HOBrSUM*HOBrGASFRAC
                  HOBr_BOX_L = HOBrSUM-HOBr_BOX_G

                  IF (VOL_SLA.gt.1.e-30_fp) THEN

                     ! Approximate particles as spherical
                     ! SAD:   cm2/cm3
                     ! RAD:   m
                     ! NDENS: #/m3
                     IF (LWACCMAER) THEN
                        ! Number density is a given parameter
                        NDENS_AER_BOX = NWACCM(I_SLA)
                        ! Derive effective radius from volume and
                        ! known properties of log-normal aerosol
                        ! distributions
                        RAD_AER_BOX = (VOL_SLA*LN_VCONV(I_SLA))
     &                                     **(1.0e+0_fp/3.0e+0_fp)
                        SAD_AER_BOX = 1.0e-2_fp*(LN_ACONV(I_SLA)*
     &                                 RAD_AER_BOX*RAD_AER_BOX)
                     ELSE
                        ! Calculate SLA parameters (Grainger 1995)
                        SAD_AER_BOX=SLA_VA*(VOL_SLA**0.751e+0_fp)
                        RAD_AER_BOX=SLA_VR*SLA_RR*(VOL_SLA**0.249e+0_fp) 
                        NDENS_AER_BOX = VOL_SLA*3.e+0_fp/
     &                       (4.e+0_fp*PI*(RAD_AER_BOX**3.e+0_fp))
                     ENDIF
                     KG_AER_BOX  = RHO_AER_BOX*VOL_SLA*             ! kg
     &                             State_Met%AIRVOL(I,J,L)


                     ! DENAIR in #/m3 - convert to #/cm3
                     ! RHO_AER_BOX in kg/m3 - convert to g/cm3
                     ! RAD_AER_BOX in m - convert to cm
                     CALL CALC_SLA_GAMMA(DENAIR*1.e-6_fp,
     &                   TCENTER,PCENTER,
     &                   W_H2SO4,H2OSUM,HClSUM,HBrSUM,HOBrSUM,
     &                   ClNO3SUM,BrNO3SUM,RHO_AER_BOX*1.e-3_fp,
     &                   RAD_AER_BOX*1.e+2_fp,GAMMA_BOX)
                  ELSE
                     ! Ignore SLA
                     NDENS_AER_BOX = 0.0e+0_fp
                     SAD_AER_BOX   = 0.0e+0_fp
                     RAD_AER_BOX   = 0.1e+0_fp
                     DO K=1,11
                        GAMMA_BOX(K) = 0.0e+0_fp
                     ENDDO
                  ENDIF
               ENDIF

               ! Store liquid fractions
               ! Liquid H2O is removed from the sum, then it is assumed
               ! that the pre-calculated solid H2O is taken out of this
               ! liquid total
               H2O_BOX_L = (98.09e+0_fp/18.02e+0_fp)*H2SO4_BOX_L
     &                    * (W_H2O/W_H2SO4)
               H2O_BOX_L = MAX(0e+0_fp,MIN(H2O_BOX_L-H2O_BOX_S,H2OSUM))
               H2O_BOX_G = MAX(0e+0_fp,H2O_BOX_G-(H2O_BOX_L+H2O_BOX_S))
  
               ! If very low number density, ignore settling
               AERFRAC(I,J,L,2) = 0e+0_fp 
               AERFRAC(I,J,L,3) = 0e+0_fp 
               AERFRAC(I,J,L,4) = 0e+0_fp 
               AERFRAC(I,J,L,5) = 0e+0_fp 
               AERFRAC(I,J,L,6) = 0e+0_fp 
               AERFRAC(I,J,L,7) = 0e+0_fp
 
               IF ((HNO3SUM.gt.1e+0_fp).and.
     &            (IS_SAFE_DIV(HNO3_BOX_L,HNO3SUM))) THEN
                  AERFRAC(I,J,L,2) = HNO3_BOX_L/HNO3SUM
               ENDIF

               IF ((HClSUM.gt.1e+0_fp).and.
     &            (IS_SAFE_DIV(HCl_BOX_L,HClSUM))) THEN
                  AERFRAC(I,J,L,3) = HCl_BOX_L/HClSUM
               ENDIF

               IF ((HOClSUM.gt.1e+0_fp).and.
     &            (IS_SAFE_DIV(HOCl_BOX_L,HOClSUM))) THEN
                  AERFRAC(I,J,L,4) = HOCl_BOX_L/HOClSUM
               ENDIF

               IF ((HBrSUM.gt.1e+0_fp).and.
     &            (IS_SAFE_DIV(HBr_BOX_L,HBrSUM))) THEN
                  AERFRAC(I,J,L,5) = HBr_BOX_L/HBrSUM
               ENDIF

               IF ((HOBrSUM.gt.1e+0_fp).and.
     &            (IS_SAFE_DIV(HOBr_BOX_L,HOBrSUM))) THEN
                  AERFRAC(I,J,L,6) = HOBr_BOX_L/HOBrSUM
               ENDIF

               IF ((H2OSUM.gt.1e+0_fp).and.
     &            (IS_SAFE_DIV(H2O_BOX_L,H2OSUM))) THEN
                  AERFRAC(I,J,L,7) = H2O_BOX_L/H2OSUM
               ENDIF

               ! Send properties to larger array
               ! Convert sticking coefficients into
               ! premultiplying factors ((Kirner)
               KHET_COMMON = 0.25e+0_fp*MOLEC_SPEED(TCENTER,1e+0_fp)

               ! N2O5 + H2O/HCl
               KHET_SPECIFIC= KHET_COMMON*ISR_N2O5
               KHETI_SLA(I,J,L,1)  = GAMMA_BOX(1 )*KHET_SPECIFIC
               KHETI_SLA(I,J,L,2)  = GAMMA_BOX(2 )*KHET_SPECIFIC

               ! ClNO3 + H2O/HCl/HBr
               KHET_SPECIFIC= KHET_COMMON*ISR_ClNO3
               KHETI_SLA(I,J,L,3)  = GAMMA_BOX(3 )*KHET_SPECIFIC
               KHETI_SLA(I,J,L,4)  = GAMMA_BOX(4 )*KHET_SPECIFIC
               KHETI_SLA(I,J,L,5)  = GAMMA_BOX(5 )*KHET_SPECIFIC
 
               ! BrNO3 + H2O/HCl
               KHET_SPECIFIC= KHET_COMMON*ISR_BrNO3
               KHETI_SLA(I,J,L,6)  = GAMMA_BOX(6 )*KHET_SPECIFIC
               KHETI_SLA(I,J,L,7)  = GAMMA_BOX(7 )*KHET_SPECIFIC

               ! HOCl + HCl/HBr
               KHET_SPECIFIC= KHET_COMMON*ISR_HOCl
               KHETI_SLA(I,J,L,8)  = GAMMA_BOX(8 )*KHET_SPECIFIC
               KHETI_SLA(I,J,L,9)  = GAMMA_BOX(9 )*KHET_SPECIFIC

               ! HOBr + HBr/HCl
               KHET_SPECIFIC= KHET_COMMON*ISR_HOBr
               KHETI_SLA(I,J,L,10) = GAMMA_BOX(10)*KHET_SPECIFIC
               KHETI_SLA(I,J,L,11) = GAMMA_BOX(11)*KHET_SPECIFIC

               RAD_AER(I,J,L,I_SLA)  = RAD_AER_BOX*1.e+2_fp ! cm
               RHO_AER(I,J,L,I_SLA)  = RHO_AER_BOX      ! kg/m3
               KG_AER(I,J,L,I_SLA)   = KG_AER_BOX       ! kg
               NDENS_AER(I,J,L,I_SLA)= NDENS_AER_BOX    ! #/m3
               SAD_AER(I,J,L,I_SLA)  = SAD_AER_BOX      ! cm2/cm3

         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### UCX: after CALC_STRAT_AER' )
         WRITE(6,'(a)') 'Strat. aerosol diagnostics summary:'
         WRITE(6,'(a,1x,E12.4)') 'Peak particle # density (#/m3): ',
     &                               MAXVAL(NDENS_AER(:,:,:,I_SPA))
         WRITE(6,'(a,1x,E12.4)') 'Peak particle density (kg/m3) : ',
     &                               MAXVAL(RHO_AER(:,:,:,I_SPA))
         WRITE(6,'(a,1x,E12.4)') 'Peak particle SAD (cm2/cm3)   : ',
     &                               MAXVAL(SAD_AER(:,:,:,I_SPA))
         WRITE(6,'(a,1x,E12.4)') 'Peak particle radius (cm)     : ',
     &                               MAXVAL(RAD_AER(:,:,:,I_SPA))
         WRITE(6,'(a,1x,E12.4)') 'Peak SSA/STS # density (#/m3) : ',
     &                               MAXVAL(NDENS_AER(:,:,:,I_SLA))
         WRITE(6,'(a,1x,E12.4)') 'Peak SSA/STS density (kg/m3)  : ',
     &                               MAXVAL(RHO_AER(:,:,:,I_SLA))
         WRITE(6,'(a,1x,E12.4)') 'Peak SSA/STS SAD (cm2/cm3)    : ',
     &                               MAXVAL(SAD_AER(:,:,:,I_SLA))
         WRITE(6,'(a,1x,E12.4)') 'Peak SSA/STS radius (cm)      : ',
     &                               MAXVAL(RAD_AER(:,:,:,I_SLA))
         WRITE(6,'(a,1x,E12.4)') 'NAT/NO3 fraction in max-NAT box  : ',
     &                               D_MAXN                       
         WRITE(6,'(a,3(1x,E12.4))') 'Local MR (H2O-HNO3-NAT) (v/v) : ',
     &                               D_XH2O,D_XHNO3,D_XNAT        
         WRITE(6,'(a,3(1x,E12.4))') 'Pressure (box-HNO3-sat) (Pa)  : ',
     &                               D_PPA,D_PP,D_PS        
         WRITE(6,'(a,1x,E12.4)') 'Local temperature (K)         : ',
     &                               D_TK                    
      ENDIF

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CALC_STRAT_AER
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kg_strat_aer
!
! !DESCRIPTION: Function KG\_STRAT\_AER returns the calculated mass of a  
!  stratospheric aerosol. The routine is essentially just an
!  interface to allow external routines to "see" the arrays.
!\\
!\\
! !INTERFACE:
!
      REAL(fp) FUNCTION KG_STRAT_AER (I,J,L,IAER)
!
! !INPUT PARAMETERS:
!
      INTEGER,INTENT(IN)          :: I,J,L      ! Grid indices
      INTEGER,INTENT(IN)          :: IAER       ! Aerosol index
                                                ! 1 = SSA (pure H2SO4)
                                                ! 2 = STS
                                                ! 3 = Solid PSC
!
! !REVISION HISTORY: 
!  18 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! KG_STRAT_AER begins here!
      !=================================================================

      KG_STRAT_AER = KG_AER(I,J,L,IAER)

      END FUNCTION KG_STRAT_AER
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rho_strat_aer
!
! !DESCRIPTION: Function RHO\_STRAT\_AER returns the calculated
!  stratospheric aerosol mass density.
!\\
!\\
! !INTERFACE:
!
      REAL(fp) FUNCTION RHO_STRAT_AER (I,J,L,IAER)
!
! !INPUT PARAMETERS:
!
      INTEGER,INTENT(IN)          :: I,J,L      ! Grid indices
      INTEGER,INTENT(IN)          :: IAER       ! Aerosol index:
                                                ! 1 = Liquid aerosol
                                                ! 2 = Solid aerosol
!
! !REVISION HISTORY: 
!  18 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! RHO_STRAT_AER begins here!
      !=================================================================

      RHO_STRAT_AER = RHO_AER(I,J,L,IAER)

      END FUNCTION RHO_STRAT_AER
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_strat_opt
!
! !DESCRIPTION: Subroutine GET\_STRAT\_OPT returns local optical properties 
!  for a given stratospheric aerosol. The routine is essentially just an
!  interface to allow external routines to "see" the arrays. However, local
!  aerosol radius is adjusted from liquid to effective radius for aerosol
!  optical depth calculations with liquid aerosols.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_STRAT_OPT (I,J,L,IAER,RAER,REFF,SAD,XSA)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER,  INTENT(IN)        :: I, J, L    ! Grid indices
      INTEGER,  INTENT(IN)        :: IAER       ! Aerosol index
                                                ! 1 = Liquid aerosols
                                                ! 2 = Solid PSC
!
! !OUTPUT VARIABLES:
!
      REAL(fp), INTENT(OUT) :: REFF       ! Effective radius (cm)
      REAL(fp), INTENT(OUT) :: RAER       ! Physical radius (cm)
      REAL(fp), INTENT(OUT) :: SAD        ! Surface area density (cm2/cm3)
      REAL(fp), INTENT(OUT) :: XSA        ! X-S area density (m2/m3)
!
! !REMARKS:
!  Seb Eastham writes: "I would edit GET_STRAT_OPT so that, when SAD is less
!  than some small value (say 1 nm2/cm3,  which is a vanishingly small surface
!  area), it returns SADSTRAT=XSASTRAT=0.d0 and RAER=REFF=0.1d0 for safety's
!  sake. I think that will prevent code blow-up later on."
! 
! !REVISION HISTORY: 
!  17 Apr 2013 - S. D. Eastham - Initial version
!  07 Apr 2015 - R. Yantosca   - Add error check to prevent div-by-zero errors
!                                in other areas of GEOS-Chem caused by low SAD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)       :: RMED

      !=================================================================
      ! GET_STRAT_OPT begins here!
      !=================================================================

      ! Surface area density [cm2/cm3]
      SAD  = SAD_AER(I,J,L,IAER)
      
      ! Add error check: threshold is 1e-14 cm2/cm3 (bmy, 4/7/15)
      IF ( SAD < 1e-14_fp ) THEN

         !--------------------------------------------------------------
         ! FOR SAFETY'S SAKE: Set outputs to "safe" values and 
         ! exit if very small SAD is encountered. (bmy, 4/7/15)
         !--------------------------------------------------------------
         SAD  = 0.0e+0_fp
         XSA  = 0.0e+0_fp
         RAER = 0.1e+0_fp
         REFF = 0.1e+0_fp
         
      ELSE 

         !--------------------------------------------------------------
         ! Otherwise, compute RAER, REFF, XSA normally
         !--------------------------------------------------------------

         ! Retrieve median radius
         RMED = RAD_AER(I,J,L,IAER)

         ! Use precomputed conversion factors
         ! Calculate effective (area-weighted mean) radius
         REFF = RMED * LN_RMED_TO_REFF(IAER)

         ! Calculate volume-weighted mean radius
         RAER = RMED * LN_RMED_TO_RLIQ(IAER)

         ! Standard log-normal distribution approach
         ! Note cm2/cm3 to m2/m3 = x100
         XSA = 0.25e+2_fp*SAD

      ENDIF

      END SUBROUTINE GET_STRAT_OPT
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_nat
!
! !DESCRIPTION: Subroutine CHECK\_NAT determines whether the solid PSC is 
!  composed of ice or NAT (needed for heterogeneous chemistry), or indeed
!  if there is any direct PSC calculation at all. This is important for
!  determining whether to use the JPP or Kirner scheme for ice cloud radii.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_NAT( am_I_Root, I, J, L, IS_NAT, IS_PSC, 
     &                      IS_STRAT, Input_Opt, State_Met, State_Chm )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_STRAT
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACERID_MOD,       ONLY : IDTNIT
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root  ! Root CPU?
      INTEGER,        INTENT(IN)  :: I,J,L      ! Grid indices
      TYPE(OptInput), INTENT(IN)  :: Input_Opt  ! Input options
      TYPE(MetState), INTENT(IN)  :: State_Met  ! Meteorology State object
      TYPE(ChmState), INTENT(IN)  :: State_Chm  ! Chemistry State object
!
! !OUTPUT VARIABLES:
!
      LOGICAL,        INTENT(OUT) :: IS_NAT     ! Is surface NAT?
      LOGICAL,        INTENT(OUT) :: IS_PSC     ! Are there solid PSCs?
      LOGICAL,        INTENT(OUT) :: IS_STRAT   ! Are we in the strat?
!
! !REVISION HISTORY: 
!  17 Apr 2013 - S. D. Eastham - Initial version
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm
!                                objects via the arg list
!  08 Apr 2015 - R. Yantosca   - Remove call to READ_PSC_FILE, this is
!                                now done from DO_CHEMISTRY (chemistry_mod.F) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! CHECK_NAT begins here!
      !=================================================================

      IS_STRAT = ITS_IN_THE_STRAT(I,J,L,State_Met)

      IS_PSC   = ( ( Input_Opt%LPSCCHEM ) .and.
!     &             ( STATE_PSC(I,J,L) .ge. 2 ) .and. ( IS_STRAT ) )
     &             ( STATE_PSC(I,J,L) >= 2.0_f4 ) .and. ( IS_STRAT ) )

      IS_NAT   = ( ( IS_PSC ) .and.
     &         ( State_Chm%Tracers(I,J,L,IDTNIT) .gt. TINY(1e+0_fp) ) )

      END SUBROUTINE CHECK_NAT
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ternary
!
! !DESCRIPTION: Subroutine TERNARY calculates the composition of SSA/STS 
!  aerosols using a paramaterization from Carslaw et al. "A Thermodynamic
!  Model of the System HCl-HNO3-H2SO4-H2O, Including Solubilities of HBr,
!  from <200 to 328 K". The bulk of this code was taken directly from the
!  Global Modeling Initiative implementation by David Considine.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TERNARY (PCENTER_IN,TCENTER_IN,H2OSUM_IN,H2SO4SUM,
     &                 HNO3SUM,HClSUM,HOClSUM,HBRSum,HOBrSUM,
     &                 W_H2SO4,W_H2O,W_HNO3,W_HCl,W_HOCl,W_HBr,W_HOBr,
     &                 HNO3GASFRAC,HClGASFRAC,HOClGASFRAC,
     &                 HBrGASFRAC,HOBrGASFRAC,SLA_VOL,SLA_RHO)
!
! !USES:
!
      USE CMN_GCTM_MOD,      ONLY : PI

      ! Temporary - for debug
      USE ERROR_MOD,     ONLY : IT_IS_NAN,ERROR_STOP     ! Test for NaN
      USE ERROR_MOD,     ONLY : SAFE_EXP,SAFE_DIV, DEBUG_MSG
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN)  :: PCENTER_IN   ! Pressure (hPa)
      REAL(fp), INTENT(IN)  :: TCENTER_IN   ! Temperature (K)
      REAL(fp), INTENT(IN)  :: H2OSUM_IN    ! Total H2O mixing ratio
      REAL(fp), INTENT(IN)  :: H2SO4SUM     ! Liquid H2SO4 mixing ratio
      REAL(fp), INTENT(IN)  :: HNO3SUM      ! Total HNO3 mixing ratio
      REAL(fp), INTENT(IN)  :: HClSUM       ! Total HCl mixing ratio
      REAL(fp), INTENT(IN)  :: HOClSUM      ! Total HOCl mixing ratio
      REAL(fp), INTENT(IN)  :: HBrSUM       ! Total HBr mixing ratio
      REAL(fp), INTENT(IN)  :: HOBrSUM      ! Total HOBr mixing ratio
!
! !OUTPUT VARIABLES:
!
      REAL(fp), INTENT(OUT) :: W_H2SO4      ! kg H2SO4/kg SLA
      REAL(fp), INTENT(OUT) :: W_H2O        ! kg H2O  /kg SLA
      REAL(fp), INTENT(OUT) :: W_HNO3       ! kg HNO3 /kg SLA
      REAL(fp), INTENT(OUT) :: W_HCl        ! kg HCl  /kg SLA
      REAL(fp), INTENT(OUT) :: W_HOCl       ! kg HOCl /kg SLA
      REAL(fp), INTENT(OUT) :: W_HBr        ! kg HBr  /kg SLA
      REAL(fp), INTENT(OUT) :: W_HOBr       ! kg HOBr /kg SLA
      REAL(fp), INTENT(OUT) :: HNO3GASFRAC  ! Gas fraction HNO3     
      REAL(fp), INTENT(OUT) :: HClGASFRAC   ! Gas fraction HCl     
      REAL(fp), INTENT(OUT) :: HOClGASFRAC  ! Gas fraction HOCl     
      REAL(fp), INTENT(OUT) :: HBrGASFRAC   ! Gas fraction HBr     
      REAL(fp), INTENT(OUT) :: HOBrGASFRAC  ! Gas fraction HOBr     
      REAL(fp), INTENT(OUT) :: SLA_VOL      ! Aerosol volume (m3/m3)
      REAL(fp), INTENT(OUT) :: SLA_RHO      ! Aer. mass density (kg/m3)
!
! !REVISION HISTORY: 
!  19 Apr 2013 - S. D. Eastham - Initial version
!  16 Apr 2015 - M. Yannetti   - Changed EXP to SAFE_EXP in H_HCL
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Derived inputs
      REAL(fp) :: H2OSUM
      REAL(fp) :: TCENTER
      REAL(fp) :: PCENTER

      ! Partial pressures
      REAL(fp) :: PATMH2O
      REAL(fp) :: PATMHNO3
      REAL(fp) :: PATMHCl
      REAL(fp) :: PATMHOCl
      REAL(fp) :: PATMHBr
      REAL(fp) :: PATMHOBr

      ! Molar densities (mol/m3)
      REAL(fp) :: MOLDENS_H2SO4

      ! Mass totals
      REAL(fp) :: M_H2SO4,M_HNO3
      REAL(fp) :: M_HCl,M_HOCl
      REAL(fp) :: M_HBr,M_HOBr

      ! Binary solutions denoted with BIN
      ! Mole fractions
      REAL(fp) :: X_H2SO4_BIN
      REAL(fp) :: X_HNO3_BIN

      ! Mass fractions
      REAL(fp) :: M_H2SO4_BIN
      REAL(fp) :: M_HNO3_BIN

      ! Effective Henry's Law coefficients
      REAL(fp) :: H_H2SO4_BIN
      REAL(fp) :: H_HNO3_BIN
      REAL(fp) :: H_HCL, H_HOCL
      REAL(fp) :: H_HBr, H_HOBr

      ! Frost point
      REAL(fp) :: T_ICE

      ! Equilibrium vapor pressure
      REAL(fp) :: PVAP_HNO3
      REAL(fp) :: PVAP_HCl
      REAL(fp) :: PVAP_HBr
      REAL(fp) :: PVAP_HOBr

      ! R in m3.atm/(mol K)
      REAL(fp), PARAMETER :: R_ATM = 8.205e-5_fp

      ! Transitional variables
      REAL(fp) :: DENSITY
      REAL(fp) :: TEMPERATURE
      REAL(fp) :: PR
      REAL(fp) :: TR
      REAL(fp) :: TT

      ! Coefficients (q, k) for calculation of H* for H2SO4 and HNO3
      REAL(fp),DIMENSION(10) :: QN,QS
      REAL(fp),DIMENSION(7)  :: KN,KS

      ! Derived parameters
      REAL(fp) :: A,B,C,PHI
      
      ! Debug message
      CHARACTER(LEN=255)   :: DBGMSG

! ======================================================================
        DATA QN/14.5734e+0_fp,0.0615994e+0_fp,-1.14895e+0_fp,
     &   0.691693e+0_fp,-0.098863e+0_fp,
     &   0.0051579e+0_fp,0.123472e+0_fp,-0.115574e+0_fp,
     &   0.0110113e+0_fp,0.0097914e+0_fp/
        DATA QS/14.4700e+0_fp,0.0638795e+0_fp,-3.29597e+0_fp,
     &   1.778224e+0_fp,-0.223244e+0_fp,
     &   0.0086486e+0_fp,0.536695e+0_fp,-0.335164e+0_fp,
     &   0.0265153e+0_fp,0.0157550e+0_fp/
        DATA KN/-39.136e+0_fp,6358.4e+0_fp,83.29e+0_fp,
     &   -17650.0e+0_fp,198.53e+0_fp,
     &   -11948.e+0_fp,-28.469e+0_fp/
        DATA KS/-21.661e+0_fp,2724.2e+0_fp,51.81e+0_fp,
     &   -15732.0e+0_fp,47.004e+0_fp,
     &   -6969.0e+0_fp,-4.6183e+0_fp/
! ======================================================================

      SAVE QN, QS, KN, KS

      !=================================================================
      ! TERNARY begins here!
      !=================================================================

      ! Routine only valid for certain limits
      H2OSUM = MAX(H2OSUM_IN,5.0e-7_fp)
      PCENTER = MAX(PCENTER_IN,5.0e+0_fp)
      TCENTER = TCENTER_IN

      ! Calculate partial pressure of H2O & HNO3
      ! PCENTER is in hPa
      PATMH2O  = H2OSUM  * PCENTER / 1013.0e+0_fp

      ! Carslaw only valid for 2e-5 < PPH2O < 2e-3 (hPa)
      PATMH2O = MAX(PATMH2O,1.9738465e-8_fp)
      PATMH2O = MIN(PATMH2O,1.9738465e-6_fp)

      PATMHNO3 = HNO3SUM * PCENTER / 1013.0e+0_fp
      PATMHCl  = HClSUM  * PCENTER / 1013.0e+0_fp
      PATMHOCl = HOClSUM * PCENTER / 1013.0e+0_fp
      PATMHBr  = HBrSUM  * PCENTER / 1013.0e+0_fp
      PATMHOBr = HOBrSUM * PCENTER / 1013.0e+0_fp

      ! Moles of H2SO4 per m3 air
      MOLDENS_H2SO4 = 100.e+0_fp*PCENTER*H2SO4SUM/(8.314e+0_fp*TCENTER)

      ! Nucleation temperature of ice
      T_ICE = 2668.70e+0_fp/
     &  (10.4310e+0_fp-(LOG(PATMH2O)+LOG(760.0e+0_fp))/LOG(10.0e+0_fp))

      ! Pressure relation
      PR = LOG(PATMH2O)+18.4e+0_fp
      
      ! Therefore if temperature lower, set to T_ICE-3
      IF (TCENTER .lt. (T_ICE-3.0e+0_fp)) THEN
         TCENTER = (T_ICE-3.0e+0_fp)
      ENDIF

      IF (TCENTER .lt. 185.0e+0_fp) THEN
         TCENTER = 185.0e+0_fp
      ENDIF

      ! ??
      TT = TCENTER * R_ATM * MOLDENS_H2SO4
      
      ! Temperature relation
      TR = 1.0e+4_fp/TCENTER-43.4782608e+0_fp

      ! Determine H2SO4/H2O pure solution concentration
      ! Mole fraction of H2SO4 in binary solution
      X_H2SO4_BIN = 1.0e+0_fp/(2.0e+0_fp*(KS(3)+KS(4)/TCENTER))*
     &   (-KS(1)-KS(2)/TCENTER-((KS(1)+KS(2)/TCENTER)**
     &   2.0e+0_fp-4.0e+0_fp*(KS(3)+KS(4)/TCENTER)*(KS(5)+KS(6)/
     &   TCENTER+KS(7)*LOG(TCENTER)-LOG(PATMH2O)))**0.5e+0_fp)

      ! Molality (mol H2SO4/kg H2O) in binary solution
      M_H2SO4_BIN = 55.51e+0_fp*X_H2SO4_BIN/(1.0e+0_fp-X_H2SO4_BIN)

      IF ((TCENTER.le.215.0e+0_fp).AND.(PATMHNO3.gt.0.0e+0_fp)) THEN
         ! Determine HNO3/H2SO4/H2O solution composition
         H_H2SO4_BIN = EXP(QS(1)+QS(2)*TR**2+(QS(3)+QS(4)*TR+
     &      QS(5)*TR**2+QS(6)*TR**3)*PR+(QS(7)+QS(8)*TR+
     &      QS(9)*TR**2)*PR**2+QS(10)*TR*PR**3)
         X_HNO3_BIN=1.0e+0_fp/(2.0e+0_fp*(KN(3)+KN(4)/TCENTER))*
     &      (-KN(1)-KN(2)/TCENTER-((KN(1)+KN(2)/TCENTER)**
     &      2-4.0e+0_fp*(KN(3)+KN(4)/TCENTER)*(KN(5)+
     &      KN(6)/TCENTER+KN(7)*LOG(TCENTER)-LOG(PATMH2O)
     &      ))**0.5e+0_fp)
         M_HNO3_BIN=55.51e+0_fp*X_HNO3_BIN/(1.0e+0_fp-X_HNO3_BIN)
         H_HNO3_BIN=EXP(QN(1)+QN(2)*TR**2+(QN(3)+QN(4)*TR+QN(5)*
     &      TR**2+QN(6)*TR**3)*PR+(QN(7)+QN(8)*TR+QN(9)*TR**2)*
     &      PR**2+QN(10)*TR*PR**3)
         A=(TT*H_HNO3_BIN*M_HNO3_BIN**2-TT*H_H2SO4_BIN*M_HNO3_BIN*
     &      M_H2SO4_BIN-2.0e+0_fp*M_HNO3_BIN**2*M_H2SO4_BIN+
     &      M_HNO3_BIN*M_H2SO4_BIN**2+H_HNO3_BIN*M_HNO3_BIN*
     &      M_H2SO4_BIN*PATMHNO3-H_H2SO4_BIN*M_H2SO4_BIN**2*
     &      PATMHNO3)/(M_HNO3_BIN**2-M_HNO3_BIN*M_H2SO4_BIN)
         B=M_H2SO4_BIN*(-2.0e+0_fp*TT*H_HNO3_BIN*M_HNO3_BIN+TT*
     &      H_H2SO4_BIN*M_H2SO4_BIN+M_HNO3_BIN*M_H2SO4_BIN-
     &      H_HNO3_BIN*M_H2SO4_BIN*PATMHNO3)/(M_HNO3_BIN-
     &      M_H2SO4_BIN)
         C=(TT*H_HNO3_BIN*M_HNO3_BIN*M_H2SO4_BIN**2)/
     &      (M_HNO3_BIN-M_H2SO4_BIN)
         PHI=ATAN(SQRT(4.0e+0_fp*(A**2-3.0e+0_fp*B)**3-(-2.0e+0_fp*A**3+
     &      9.0e+0_fp*A*B-27.0e+0_fp*C)**2)/(-2.0e+0_fp*A**3+9.0e+0_fp*A
     &      *B-27.0e+0_fp*C))
         IF (PHI.lt.0.e+0_fp) THEN
            PHI = PHI + PI
         ENDIF
         M_H2SO4=-1.0e+0_fp/3.0e+0_fp*(A+2.0e+0_fp*
     &      SQRT(A**2-3.0e+0_fp*B)*
     &      COS((PI+PHI)/3.0e+0_fp))
         M_HNO3=M_HNO3_BIN*(1.0e+0_fp-M_H2SO4/M_H2SO4_BIN)
         W_H2SO4 = M_H2SO4*0.098076e+0_fp/(1.0e+0_fp+M_H2SO4*
     &      0.098076e+0_fp+M_HNO3*0.063012e+0_fp)

         ! Check for low H2SO4
         IF (M_H2SO4 .lt. 0e+0_fp) THEN
             M_H2SO4 = 0.0e+0_fp
             M_HNO3 = M_HNO3_BIN
             W_H2SO4 = 0.0e+0_fp
         ENDIF

         PVAP_HNO3=M_HNO3/(H_HNO3_BIN*M_HNO3/(M_HNO3+
     &      M_H2SO4)+H_H2SO4_BIN*M_H2SO4/(M_HNO3+M_H2SO4))
         W_HNO3 = (M_HNO3*0.063012e+0_fp)/(1.0e+0_fp+M_H2SO4*
     &      0.098076e+0_fp+M_HNO3*0.063012e+0_fp)

         HNO3GASFRAC=(1.0e+0_fp-(PATMHNO3-PVAP_HNO3)/PATMHNO3)
         M_HNO3=MAX(M_HNO3,0.0e+0_fp)
         W_HNO3=MAX(W_HNO3,0.0e+0_fp) ! needed for small underflow (mdy 04/15)

      ELSE
         ! Solution is pure H2SO4/H2O
         M_H2SO4 = M_H2SO4_BIN
         M_HNO3 = 0.0e+0_fp
         W_H2SO4 = M_H2SO4_BIN*0.098076e+0_fp/(1.0e+0_fp+M_H2SO4_BIN*
     &      0.098076e+0_fp)
         W_HNO3 = 0.0e+0_fp
         PVAP_HNO3 = 0.0e+0_fp
         HNO3GASFRAC=1.0e+0_fp
      ENDIF
      ! Handle HCl (Luo et al., Vapor pressures of
      ! H2SO4/HNO3/HCl/HBr/H2O solutions to low stratospheric
      ! temperatures, 1995)
      IF (PATMHCL.gt.0.0e+0_fp) THEN
         H_HCL = EXP(-(21.0e+0_fp+46.610e+0_fp*W_HNO3+
     &      4.0690e+0_fp*W_H2SO4-
     &      4.8370e+0_fp*SQRT(W_HNO3)+2.1860e+0_fp*
     &      SQRT(W_H2SO4)-63.0e+0_fp*
     &      W_HNO3**2-40.170e+0_fp*W_HNO3*W_H2SO4-
     &      1.5710e+0_fp*W_H2SO4**2)-
     &      1.0e+0_fp/TCENTER*(-7437.0e+0_fp-8327.80e+0_fp*
     &      W_HNO3+1300.90e+0_fp*
     &      W_H2SO4+1087.20e+0_fp*SQRT(W_HNO3)-242.710e+0_fp*
     &      SQRT(W_H2SO4)+
     &      18749.0e+0_fp*W_HNO3**2+18500.0e+0_fp*W_HNO3*W_H2SO4+
     &      5632.0e+0_fp*W_H2SO4**2)-LOG(W_HNO3+0.610e+0_fp*W_H2SO4)-
     &      LOG(36.461e+0_fp/(1000.0e+0_fp+98.076e+0_fp*
     &      M_H2SO4+63.012e+0_fp*
     &      M_HNO3)))*1013.0e+0_fp
         M_HCl = (1.0e+0_fp/R_ATM/TCENTER*PATMHCL)/
     &      (MOLDENS_H2SO4/M_H2SO4 + 1.0e+0_fp/R_ATM/TCENTER/H_HCL)
         W_HCL=M_HCL*36.461e+0_fp/(1.0e+3_fp+98.076e+0_fp*M_H2SO4+
     &      63.012e+0_fp*M_HNO3)
         PVAP_HCl = M_HCl/H_HCl
         HClGASFRAC=1.0e+0_fp-(PATMHCL-PVAP_HCL)/PATMHCL
      ELSE
         W_HCl=0.0e+0_fp
         HClGASFRAC=1.0e+0_fp
      ENDIF

      ! Now HOCl
      IF (PATMHOCl.gt.0.0e+0_fp) THEN
         H_HOCl=EXP(6.49460e+0_fp-(-0.041070e+0_fp+54.56e+0_fp/TCENTER)*
     &      (M_H2SO4+M_HNO3)-5862.0e+0_fp*(1.0e+0_fp/298.15e+0_fp-
     &      1.0e+0_fp/TCENTER))
         M_HOCl=(1.0e+0_fp/R_ATM/TCENTER*PATMHOCl)/
     &      (MOLDENS_H2SO4/M_H2SO4 + 1.0e+0_fp/R_ATM/TCENTER/H_HOCL)
         W_HOCL=M_HOCL*52.46e+0_fp/(1.0e+3_fp+98.076e+0_fp*M_H2SO4+
     &      63.012e+0_fp*M_HNO3)
         ! Realistically expect no gas phase removal
         HOCLGASFRAC=1.0e+0_fp 
      ELSE
         W_HOCl=0.0e+0_fp
         HOClGASFRAC=1.0e+0_fp
      ENDIF
 
      ! Now HBr (Luo et al., Vapor pressures of
      ! H2SO4/HNO3/HCl/HBr/H2O solutions to low stratospheric
      ! temperatures, 1995)
      IF (PATMHBr.gt.0.0e+0_fp) THEN
         H_HBr = EXP(-(17.83e+0_fp+1.02e+0_fp*W_HNO3-1.08e+0_fp*W_H2SO4+
     &      3.9e+0_fp*SQRT(W_HNO3)+4.38e+0_fp*SQRT(W_H2SO4)-8.87e+0_fp*
     &      W_HNO3**2-17.0e+0_fp*W_HNO3*W_H2SO4+3.73e+0_fp*W_H2SO4**2)-
     &      1.0e+0_fp/TCENTER*(-8220.5e+0_fp-362.76e+0_fp*
     &      W_HNO3+658.93e+0_fp*
     &      W_H2SO4-914.0e+0_fp*SQRT(W_HNO3)-955.3e+0_fp*SQRT(W_H2SO4)+
     &      9976.6e+0_fp*W_HNO3**2+19778.5e+0_fp*W_HNO3*W_H2SO4+
     &      7680.0e+0_fp*W_H2SO4**2)-LOG(W_HNO3+0.410e+0_fp*W_H2SO4)-
     &      LOG(36.461e+0_fp/(1000.0e+0_fp+98.076e+0_fp*
     &      M_H2SO4+63.012e+0_fp*
     &      M_HNO3)))*1013.0e+0_fp
         M_HBr = (1.0e+0_fp/R_ATM/TCENTER*PATMHBr)/
     &      (MOLDENS_H2SO4/M_H2SO4 + 1.0e+0_fp/R_ATM/TCENTER/H_HBr)
         W_HBr=M_HBr*80.91e+0_fp/(1.0e+3_fp+98.076e+0_fp*M_H2SO4+
     &      63.012e+0_fp*M_HNO3)
         PVAP_HBr = M_HBr/H_HBr
         HBrGASFRAC=1.0e+0_fp-(PATMHBr-PVAP_HBr)/PATMHBr
      ELSE
         W_HBr=0.0e+0_fp
         HBrGASFRAC=1.0e+0_fp
      ENDIF

      ! Finally HOBr (Hanson and Ravishankara, Heterogeneous
      ! chemistry of Bromine species in sulfuric acid under
      ! stratospheric conditions, 1995)
      IF (PATMHOBr.gt.0.0e+0_fp) THEN
         ! Hanson and Ravishankara state that the volume-based
         ! Henry's Law coefficient for HOBr in H2SO4 is 10^6 M/atm.
         ! The molality-based Henry's law constant, H_HOBr, is
         ! therefore:
         H_HOBr = (1.0e+6_fp) * MOLDENS_H2SO4 / M_H2SO4
         M_HOBr = (1.0e+0_fp/R_ATM/TCENTER*PATMHOBr)/
     &      (MOLDENS_H2SO4/M_H2SO4 + 1.0e+0_fp/R_ATM/TCENTER/H_HOBr)
         W_HOBr = M_HOBr*96.911e+0_fp/(1.0e+3_fp+98.076e+0_fp*M_H2SO4+
     &      63.012e+0_fp*M_HNO3)
         PVAP_HOBr = M_HOBr/H_HOBr
         HOBrGASFRAC=1.0e+0_fp-(PATMHOBr-PVAP_HOBr)/PATMHOBr
      ELSE
         W_HOBr=0e+0_fp
         HOBrGASFRAC=1.0e+0_fp
      ENDIF

      ! Take W_H2O as remainder
      W_H2O = 1.e+0_fp-(W_H2SO4+W_HNO3+W_HCl+W_HOCl+W_HBr+W_HOBr)

      ! Aerosol mass density in kg/m3 aerosol
      SLA_RHO = CARSLAW_DENSITY(M_H2SO4,M_HNO3,TCENTER)

      ! Aerosol volume in m3/m3 air
      SLA_VOL = (MOLDENS_H2SO4*98.076e+0_fp/W_H2SO4/SLA_RHO)*1.e-3_fp

      END SUBROUTINE TERNARY
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: carslaw_density
!
! !DESCRIPTION: Function CARSLAW\_DENSITY determines the density of a
!  sol'n through a relationship from Carslaw et al.. Result is in kg/m3.
!\\
!\\
! !INTERFACE:
!
      REAL(fp) FUNCTION CARSLAW_DENSITY(CS,CN,T)
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: CS         ! H2SO4 molality (mol H2SO4/kg solvent)
      REAL(fp), INTENT(IN) :: CN         ! HNO3 molality (mol HNO3/kg solvent)
      REAL(fp), INTENT(IN) :: T          ! Temperature (K)
!
! !REVISION HISTORY: 
!  19 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)                :: DENSS,DENSN

      !=================================================================
      ! CARSLAW_DENSITY begins here!
      !=================================================================

        DENSS=1000.0e+0_fp+123.64e+0_fp*CS-5.6e-4_fp*CS*T**2  
     &       -29.54e+0_fp*CS**1.5e+0_fp + 1.814e-4_fp*CS**1.5e+0_fp
     &       *T**2+2.343e+0_fp*CS**2  -1.487e-3_fp*CS**2*T  
     &       -1.324e-5_fp*CS**2*T**2

        DENSN=1000.0e+0_fp+85.107e+0_fp*CN-5.043e-4_fp*CN*T**2  
     &       -18.96e+0_fp*CN**1.5e+0_fp + 1.427e-4_fp*CN**1.5e+0_fp
     &       *T**2+1.458e+0_fp*CN**2  -1.198e-3_fp*CN**2*T  
     &       -9.703e-6_fp*CN**2*T**2

        CARSLAW_DENSITY=1.0e+0_fp/((1.0e+0_fp/DENSS*CS/(CS+CN)  
     &  +1.0e+0_fp/DENSN*CN/(CS+CN)))
      RETURN

      END FUNCTION CARSLAW_DENSITY
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_fallvel
!
! !DESCRIPTION: Function CALC\_FALLVEL calculates the terminal velocity of a
!  solid particle.
!\\
!\\
! !INTERFACE:
!
      FUNCTION CALC_FALLVEL(DENSITY,RADIUS,TCENTER,PCENTER) RESULT(VEL)
!
! !USES:
!
      USE CMN_GCTM_MOD,    ONLY : PI, g0
      USE ERROR_MOD,       ONLY : ERROR_STOP
!
! !INPUT PARAMETERS:
!
      REAL(fp),INTENT(IN)   :: RADIUS  ! Particle radius (cm)
      REAL(fp),INTENT(IN)   :: DENSITY ! Particle density (kg/m3)
      REAL(fp),INTENT(IN)   :: TCENTER ! Local temperature (K)
      REAL(fp),INTENT(IN)   :: PCENTER ! Local pressure (kPa)
!
! !OUTPUT VARIABLES:
!
      REAL(fp)              :: VEL ! Fall velocity (m/s)
!
! !REMARKS:
! (1) A remark
! 
! !REVISION HISTORY: 
!  11 Aug 2012 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Vy ! Intermediate velocity (m/s)
      REAL(fp),PARAMETER    :: eta=6.45e-8_fp ! Constant (kg/(msK))
      REAL(fp)              :: val_x ! Dimensionless variable
      REAL(fp),DIMENSION(3) :: alpha ! Auxiliary variables
      REAL(fp)              :: PR ! Pressure times radius

      !=================================================================
      ! CALC_FALLVEL begins here!
      !=================================================================

      DATA ALPHA/1.49e-5_fp,5.02e-6_fp,2.64e-5_fp/

      ! Sanity check
      IF ((RADIUS.le.0.e+0_fp).or.(DENSITY.le.0.e+0_fp)) THEN
         VEL=0.e+0_fp
      ELSE
         ! PCENTER (kPa -> Pa) = *1.d3
         ! RADIUS  (cm  -> m ) = *1.d-2
         ! Therefore multiply PR by 10
         PR = PCENTER * RADIUS * 10e+0_fp
         VAL_X = -1.0e+0_fp*PR/(ALPHA(3)*TCENTER)
         VAL_X = ALPHA(2)*TCENTER*EXP(VAL_X)/PR
         VAL_X = 1.0e+0_fp + VAL_X + (ALPHA(1)*TCENTER/PR)
         Vy = g0*DENSITY*RADIUS*RADIUS*(1.e-4_fp)/(4.5*ETA*TCENTER)
         VEL = 0.893e+0_fp * Vy * VAL_X
      ENDIF

      ! Velocities should be of the order of 0.1 m/s
      IF (VEL.gt.10.e+0_fp) THEN
         CALL ERROR_STOP(' Excessive fall velocity? ', 
     &                   ' CALC_FALLVEL, UCX_mod')
      ENDIF

      END FUNCTION CALC_FALLVEL
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cacl_sla_gamma
!
! !DESCRIPTION: Subroutine CALC\_SLA\_GAMMA calculates 11 different sticking 
!  coefficients on the surface of local stratospheric liquid aerosols, 
!  relevant to each of the 11 reactions listed in Kirner's paper.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_SLA_GAMMA( NDENS, T, P, WT_FRC, H2OSUM, HClSUM,
     &                           HBrSUM, HOBrSUM, ClNO3SUM, BrNO3SUM,
     &                           RHO, ARAD, RXNGAMMA )
!
! !USES:
!
      USE CMN_GCTM_MOD,       ONLY : PI
      ! Temporary - for debug
      USE ERROR_MOD,     ONLY : IT_IS_NAN,ERROR_STOP     ! Test for NaN
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN)  :: NDENS    ! Air number density (molec/cm3)
      REAL(fp), INTENT(IN)  :: T        ! Temperature (K)
      REAL(fp), INTENT(IN)  :: P        ! Pressure (hPa)
      REAL(fp), INTENT(IN)  :: WT_FRC   ! Weight fraction of H2SO4 (kg/kg)
      REAL(fp), INTENT(IN)  :: H2OSUM   ! H2O mixing ratio
      REAL(fp), INTENT(IN)  :: HClSUM   ! HCl mixing ratio
      REAL(fp), INTENT(IN)  :: HBrSUM   ! HBr mixing ratio
      REAL(fp), INTENT(IN)  :: HOBrSUM  ! HOBr mixing ratio
      REAL(fp), INTENT(IN)  :: ClNO3SUM ! ClNO3 mixing ratio
      REAL(fp), INTENT(IN)  :: BrNO3SUM ! BrNO3 mixing ratio
      REAL(fp), INTENT(IN)  :: RHO      ! STS density (g/cm3)
      REAL(fp), INTENT(IN)  :: ARAD     ! SLA radius (cm)
!
! !OUTPUT VARIABLES:
!
      REAL(fp), INTENT(OUT) :: RXNGAMMA(11) ! Premultiplying factors
!
! !REVISION HISTORY: 
!  10 Oct 2012 - S. D. Eastham - Initial version
!  15 Apr 2015 - M. Yannetti   - Adjustment to f8 from fp due to overflow concerns
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(f8)               :: WT      ! Weight percentage H2SO4 (100*kg/kg)
      REAL(f8)               :: H2OPP   ! Partial pressure of H2O (hPa)
      REAL(f8)               :: HClPP   ! Partial pressure of HCl (atm)
      REAL(f8)               :: HBrPP   ! Partial pressure of HBr (atm)
      REAL(f8)               :: HOBrPP  ! Partial pressure of HOBr (atm)
      REAL(f8)               :: ClNO3PP ! Partial p. of ClONO2 (atm)
      REAL(f8)               :: BrNO3PP ! Partial p. of BrONO2 (atm)
      REAL(f8)               :: PSATH2O ! Water vapor sat. pressure (hPa)
      REAL(f8)               :: ACTH2O  ! Activity of water
      REAL(f8)               :: MOLAL   ! Molality of H2SO4 (mol H2SO4/kg solvent)
      REAL(f8), DIMENSION(3) :: Z       ! Parameters for H2SO4 sol'n
      REAL(f8)               :: M_H2SO4 ! Mass of H2SO4

      ! HOBr parameters
      REAL(f8)               :: c_HOBr
      REAL(f8)               :: SHOBr
      REAL(f8)               :: HHOBr
      REAL(f8)               :: DHOBr
      REAL(f8)               :: kHOBr_HCl
      REAL(f8)               :: GHOBrrxn
      REAL(f8)               :: lHOBr
      REAL(f8)               :: fHOBr
      REAL(f8)               :: gHOBr_HCl

      ! HOCl parameters
      REAL(f8)               :: c_HOCl
      REAL(f8)               :: SHOCl
      REAL(f8)               :: HHOCl
      REAL(f8)               :: DHOCl
      REAL(f8)               :: kHOCl_HCl
      REAL(f8)               :: GHOClrxn
      REAL(f8)               :: lHOCl
      REAL(f8)               :: fHOCl
      REAL(f8)               :: gHOCl_HCl

      ! ClNO3 parameters
      REAL(f8)               :: c_ClNO3
      REAL(f8)               :: SClNO3
      REAL(f8)               :: HClNO3
      REAL(f8)               :: DClNO3
      REAL(f8)               :: GClNO3rxn
      REAL(f8)               :: lClNO3
      REAL(f8)               :: fClNO3
      REAL(f8)               :: gClNO3
      REAL(f8)               :: gClNO3_HCl
      REAL(f8)               :: gClNO3_H2O
 
      ! N2O5 parameters
      REAL(f8), DIMENSION(3) :: AK

      ! Other parameters
      REAL(f8) :: kH2O,kH,khdr,GbH2O,HHCl,MHCl,kHCl,GbHCl,Gs,FHCl,Gsp
      REAL(f8) :: GbHClp, Gb, khydr, kII, k_dl

      ! Interim variables
      REAL(f8)               :: X,A,H,T_THRESHOLD,aH
      REAL(f8), PARAMETER    :: MAX_T_DIFF = 6.0e+0_fp

      ! Control whether to run calculations
      LOGICAL              :: HClOK, HOBrOK

      ! Debug variables
      INTEGER              :: I
      CHARACTER(LEN=255)   :: DBGMSG
      !=================================================================
      ! CALC_SLA_GAMMA begins here!
      !=================================================================

      PSATH2O = EXP(18.452406985e+0_f8-3505.1578807e+0_f8
     &     /T-330918.55082e+0_f8/(T*T)
     &     +12725068.262e+0_f8/(T*T*T))       ! Saturation pressure of H2O
      H2OPP = H2OSUM * P                  ! Partial pressure of H2O
      ACTH2O = MAX((H2OPP/PSATH2O),1.0e+0_f8) ! Water activity

      ! Calculate molality of solution
      !WT = MIN(100.0e+0_fp,100.0e+0_fp*WT_FRC) ! Convert from fraction to %
      WT = 100.0e+0_f8*WT_FRC ! Convert from fraction to %
      MOLAL = 1000.0e+0_f8 * (WT/98.0e+0_f8/(100.0-WT))

      ! Parameters for H2SO4 solution
      !----------------------------------------------------------
      ! The solution density is calculated earlier, including
      ! contributions from HNO3. This code treats it as a binary
      ! solution - so far this is just a kludge. Need to update
      ! all this code to acknowledge the presence of at least
      ! HNO3 (e.g. X is still calculated based on pure H2O
      ! solvent!)
      !----------------------------------------------------------
      !Z(1) =   0.12364e+0_fp-5.6d-7*T*T
      !Z(2) =   -0.02954e+0_fp+1.814d-7*T*T
      !Z(3) =   2.343d-3-1.487d-6*T-1.324d-8*T*T
      !RHO  =   1.0e+0_fp+Z(1)*MOLAL+Z(2)*MOLAL**1.5+Z(3)*MOLAL*MOLAL
      !----------------------------------------------------------
      M_H2SO4 =   RHO*WT/9.8 ! Molality (mol H2SO4/kg solvent)
      X    =   WT/(WT+(100.-WT)*98./18.)
      A    =   169.5+5.18*WT-0.0825*WT*WT+3.27e-3_f8*WT*WT*WT
      T_THRESHOLD = 144.11+0.166*WT-0.015*WT*WT+2.18e-4_f8*WT*WT*WT
      IF ((T-T_THRESHOLD).gt.MAX_T_DIFF) THEN
         H = A*T**(-1.43)*EXP(448./(T-T_THRESHOLD))
      ELSE
         H = A*T**(-1.43)*EXP(448./MAX_T_DIFF)
      ENDIF

      aH   =   EXP(60.51-0.095*WT+0.0077*WT*WT-1.61e-5_f8*WT*WT*WT
     &     -(1.76+2.52e-4_f8*WT*WT)*SQRT(T) +
     &     (-805.89+253.05*WT**0.076)/SQRT(T))

      HClPP = HClSUM*P/1013.25e+0_f8 ! Note atm, not hPa
      ClNO3PP = ClNO3SUM*P/1013.25e+0_f8
      BrNO3PP = BrNO3SUM*P/1013.25e+0_f8
      HOBrPP = HOBrSUM*P/1013.25 ! Note atm, not hPa
 
      ! Should we bother running calculations?
      HClOK = (HClPP .gt. 1.e-30_f8)
      HOBrOK = (HOBrPP .gt. 1.e-30_f8)

      ! Reaction 1. N2O5 + H2O (hydrolysis of N2O5)
      AK(1)=-25.5265-0.133188*WT+0.0093084*WT**2-9.0194E-5_f8*WT**3
      AK(2)=9283.76+115.345*WT-5.19258*WT**2+0.0483464*WT**3
      AK(3)=-851801-22191.2*WT+766.916*WT**2-6.85427*WT**3
      RXNGAMMA(1)=exp(AK(1)+AK(2)/T+AK(3)/T**2)
      
      ! Reaction 2. N2O5 + HCl
      ! JPL 10-06 suggests near-zero gamma
      RXNGAMMA(2) = TINY(1e+0_f8)

      ! Reactions 3/4. ClNO3 + H2O/HCl
      ! Now od only if HCl concentrations are large enough
      ! (to avoid div-by-zero errors), ckeller, 2/10/15.
      IF (HClOK) THEN
         c_ClNO3     =  1474.e+0_f8*SQRT(T)
         SClNO3      =  0.306e+0_f8+24.e+0_f8/T
         HClNO3      =  1.6e-6_f8*EXP(4710.e+0_f8/T)
     &               *  EXP(-SClNO3*M_H2SO4)
         DClNO3      =  5e-8_f8*T/h
         kH2O        =  1.95e+10_f8*EXP(-2800.e+0_f8/T)
         kH          =  1.22e+12_f8*EXP(-6200.e+0_f8/T)
         khydr       =  kH2O*ACTH2O + kH*aH*ACTH2O
         GbH2O       =  4*HClNO3*0.082*T*SQRT(DClNO3*khydr)
     &               /  c_ClNO3
         HHCl        =  (0.094e+0_f8-0.61e+0_f8*X+1.2e+0_f8*X*X)
     &               *  EXP(-8.68+(8515-10718*X**0.7)/T)
         MHCl        =  HHCl *HClPP
         kHCl        =  7.9e+11_f8*aH*DClNO3*MHCl
         lClNO3      =  SQRT(DClNO3/(khydr+kHCl))
         if (lClNO3.gt.(1.e+5_f8*arad)) then
            ! Limiting rate
            fClNO3   =  arad/(3.e+0_f8*lClNO3)
         else
            fClNO3   =  1.e+0_f8/tanh(arad/lClNO3)- lClNO3/arad
         endif
         GClNO3rxn   =  fClNO3*GbH2O *SQRT(1.e+0_f8+kHCl/khydr)
         GbHCl       =  GClNO3rxn* kHCl/(kHCl+ khydr)
         Gs          =  66.12e+0_f8*EXP(-1374.e+0_f8/T)*HClNO3*MHCl
         FHCl        =  1.e+0_f8/(1.e+0_f8+0.612e+0_f8*(Gs+GbHCl)*
     &                ClNO3PP/ HClPP)
         Gsp         =  FHCl*Gs
         GbHClp      =  FHCl*GbHCl
         Gb          =  GbHClp  + GClNO3rxn* khydr/( kHCl+ khydr)
         gClNO3      =  1.e+0_f8/(1.e+0_f8+1.e+0_f8/(Gsp + Gb))
         gClNO3_HCl  =  gClNO3 *(Gsp + GbHClp)/(Gsp + Gb)
         gClNO3_H2O  =  gClNO3 - gClNO3_HCl
   
         !IF (HClOK) THEN
         RXNGAMMA(3) =  gClNO3_H2O
         RXNGAMMA(4) =  gClNO3_HCl
      ELSE
         RXNGAMMA(3) = TINY(1e+0_f8)
         RXNGAMMA(4) = TINY(1e+0_f8)
      ENDIF

      ! Reaction 5. ClNO3 + HBr
      ! Not present in JPL 10-06 for H2SO4
      RXNGAMMA(5) = TINY(1e+0_f8)

      ! Reaction 6. BrNO3 + H2O
!      RXNGAMMA(6) = 1.0/(1.0/0.88+exp(-17.832+0.245*WT))
      RXNGAMMA(6) = 1.e+0_f8/(1.e+0_f8/0.80e+0_f8+1.e+0_f8/
     &             (exp(29.2e+0_f8-0.4e+0_f8*WT )+0.11))
      
      ! Reaction 7. BrNO3 + HCl
      RXNGAMMA(7) = 0.9e+0_f8 ! JPL 10-06

      ! Reaction 8. HOCl + HCl
      IF (HClOK) THEN
         c_HOCl    =  MOLEC_SPEED(T,52.46e+0_fp) ! input variable, fp
         SHOCl     =  0.0776e+0_f8+59.18e+0_f8/T
         HHOCl     =  1.91e-6_f8*EXP(5862.4e+0_f8/T)*EXP(-SHOCl*M_H2SO4)
         DHOCl     =  6.4e-8_f8*T/H
         kHOCl_HCl =  1.25e+9_f8*aH*DHOCl*MHCl
         lHOCl     =  SQRT(DHOCl/kHOCl_HCl)
         if (lHOCl.gt.(1.e+5_f8*arad)) then
            ! Limiting rate
            fHOCl = arad/(3.e+0_f8*lHOCl)
         else
            fHOCl     =  1.e+0_f8/tanh(arad/lHOCl)- lHOCl/arad
         endif
         GHOClrxn  =  4.e+0_f8*HHOCl*0.082e+0_f8*T*
     &                sqrt(DHOCl*kHOCl_HCl)/c_HOCl
         IF (fHOCl.eq.0.) THEN
            gHOCl_HCl =  TINY(1e+0_f8)
         ELSE
            gHOCl_HCl =  1.e+0_f8/(1.e+0_f8+1.e+0_f8/
     &                 (fHOCl*GHOClrxn*FHCl))
         ENDIF
      ELSE
         gHOCl_HCl = TINY(1e+0_f8)
      ENDIF
   
      RXNGAMMA(8) = gHOCl_HCl

      ! Reaction 9. HOCl + HBr
      ! Not yet implemented for STS; JPL 10-06 suggests complex
      ! relationship, not yet sufficiently well understood or
      ! parameterized for the purposes of simulation. Ignore for now
      RXNGAMMA(9) = TINY(1e+0_f8)

      ! Reaction 10. HOBr + HCl
      IF ((HClOK).and.(HOBrOK)) THEN
         c_HOBr    =  MOLEC_SPEED(T,96.91e+0_fp) ! input variable, fp
         SHOBr     =  0.0776e+0_f8+59.18e+0_f8/T
!         HHOBr     =  30.D0
         HHOBR     = exp(-9.86e+0_f8+5427.e+0_f8/T)
         DHOBr     =  1.E-8_f8
         kII       = exp(154.e+0_f8-1.63e+0_f8*WT)*exp(-(3.85e+4_f8-
     &               478.e+0_f8*WT)/T)
         k_dl      = 7.5e+14_f8*(DHOBr*arad*1.e7_f8)
         IF (kII.gt.k_dl) kII=k_dl
         kHOBr_HCl =  kII*HHOBr*HOBrPP
!         IF (kHOBr_HCl.eq.0.) THEN   ! catch for zero (mdy, 04/15)
!           kHOBr_HCl = TINY(1e+0_fp)
!         ENDIF
         GHOBrrxn  =  4.e+0_f8*HHCl*0.082e+0_f8*T*sqrt(DHOBr*kHOBr_HCl)/
     &              c_HOBr
         lHOBr     =  sqrt(DHOBr/kHOBr_HCl)
         if (lHOBr.gt.(1.e+3_f8*arad)) then
            ! Limiting rate
            fHOBr = arad/(3.e+0_f8*lHOBr)
         else
            fHOBr     =  1.e+0_f8/tanh(arad/lHOBr)- lHOBr/arad
         endif
!         IF ((fHOBr*GHOBrrxn).eq.0.) THEN
!           gHOBr_HCl =  1.e+0_fp/(1.e+0_fp+1.e+0_fp/TINY(1e+0_fp))
!         ELSE
           gHOBr_HCl =  1.e+0_f8/(1.e+0_f8+1.e+0_f8/(fHOBr*GHOBrrxn))
!         ENDIF
         RXNGAMMA(10) = gHOBr_HCl
      ELSE
         RXNGAMMA(10) = TINY(1e+0_f8)
      ENDIF

      ! Reaction 11. HOBr + HBr
      ! Data from JPL limited; ignore for now
      RXNGAMMA(11) = TINY(1e+0_f8)

      ! SDE 2013-10-18: DEBUG
      DO I=1,11
         IF (IT_IS_NAN(RXNGAMMA(I))) THEN
            WRITE(DBGMSG,'(a,I2)') 'RXNGAMMA NaN: ', I
            CALL DEBUG_MSG( TRIM(DBGMSG) )
            WRITE(DBGMSG,'(a,E10.4)') 'NDENS: ', NDENS
            CALL DEBUG_MSG( TRIM(DBGMSG) )
            WRITE(DBGMSG,'(a,E10.4)') 'T: ', T
            CALL DEBUG_MSG( TRIM(DBGMSG) )
            WRITE(DBGMSG,'(a,E10.4)') 'P: ', P
            CALL DEBUG_MSG( TRIM(DBGMSG) )
            WRITE(DBGMSG,'(a,E10.4)') 'WT_FRC: ', WT_FRC
            CALL DEBUG_MSG( TRIM(DBGMSG) )
            WRITE(DBGMSG,'(a,E10.4)') 'H2OSUM: ', H2OSUM
            CALL DEBUG_MSG( TRIM(DBGMSG) )
            WRITE(DBGMSG,'(a,E10.4)') 'HClSUM: ', HClSUM
            CALL DEBUG_MSG( TRIM(DBGMSG) )
            WRITE(DBGMSG,'(a,E10.4)') 'HBrSUM: ', HBrSUM
            CALL DEBUG_MSG( TRIM(DBGMSG) )
            WRITE(DBGMSG,'(a,E10.4)') 'HOBrSUM: ', HOBrSUM
            CALL DEBUG_MSG( TRIM(DBGMSG) )
            WRITE(DBGMSG,'(a,E10.4)') 'ClNO3SUM: ', ClNO3SUM
            CALL DEBUG_MSG( TRIM(DBGMSG) )
            WRITE(DBGMSG,'(a,E10.4)') 'BrNO3SUM: ', BrNO3SUM
            CALL DEBUG_MSG( TRIM(DBGMSG) )
            WRITE(DBGMSG,'(a,E10.4)') 'RHO: ', RHO
            CALL DEBUG_MSG( TRIM(DBGMSG) )
            WRITE(DBGMSG,'(a,E10.4)') 'ARAD: ', ARAD
            CALL DEBUG_MSG( TRIM(DBGMSG) )
            CALL ERROR_STOP('BAD GAMMA','UCX_mod')
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE CALC_SLA_GAMMA
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: molec_speed
!
! !DESCRIPTION: Function MOLEC\_SPEED calculates the mean velocity of gas
!  phase particles based on temperature and molecular mass.
!\\
!\\
! !INTERFACE:
!
      REAL(fp) FUNCTION MOLEC_SPEED(T,MOLMASS)
!
! !USES:
!
      USE CMN_GCTM_MOD, ONLY : PI
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN)          :: T       ! Temperature (K)
      REAL(fp), INTENT(IN)          :: MOLMASS ! Molecular mass (g/mol)
!
! !REVISION HISTORY: 
!  10 Oct 2012 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! MOLEC_SPEED begins here!
      !=================================================================

      MOLEC_SPEED=SQRT(8.0e+0_fp*8.31451e+7_fp*T/(PI*MOLMASS))

      END FUNCTION MOLEC_SPEED
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_psc_file
!
! !DESCRIPTION: Subroutine READ\_PSC\_FILE initializes PSC state information
!  from a checkpoint file (binary punch file format).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_PSC_FILE( am_I_Root, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE HCO_ERROR_MOD
      USE HCO_RESTART_MOD,    ONLY : HCO_RestartGet
      USE HCO_STATE_MOD,      ONLY : HCO_State
      USE HCO_CLOCK_MOD,      ONLY : HcoClock_First 
      USE HCO_CLOCK_MOD,      ONLY : HcoClock_Rewind
      USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoState
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER,        INTENT(INOUT) :: RC
! 
! !REVISION HISTORY: 
!  28 Apr 2013 - S. D. Eastham - Initial version, based on READ_CSPEC_FILE
!  14 Feb 2014 - R. Yantosca   - Reorder DO loops for efficiency
!  16 Apr 2014 - M. Sulprizio  - Now get PSC restart file path from Input_Opt
!  23 Jul 2014 - R. Yantosca   - Remove reference to obsolete CMN_mod.F
!  14 Jan 2015 - C. Keller     - Now read from HEMCO
!   4 Mar 2015 - R. Yantosca   - Declare pointer args to HCO_GetPtr w/ REAL(f4)
!  25 Mar 2015 - C. Keller     - Now use HEMCO restart module
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!     
      ! Scalars
      LOGICAL            :: IT_EXISTS 

      ! Strings
      CHARACTER(LEN=255) :: MSG
     
      ! Pointer to HEMCO state obj.
      TYPE(HCO_State), POINTER  :: HcoState => NULL()

      !=================================================================
      ! READ_PSC_FILE begins here!
      !=================================================================

      ! Assume success until otherwise 
      RC = GIGC_SUCCESS 
 
      ! Only read on first call or after rewind
      IF ( .NOT. HcoClock_First ( .FALSE. ) .AND. 
     &     .NOT. HcoClock_Rewind( .FALSE. )      ) THEN
         RETURN
      ENDIF
 
      ! Get HEMCO state obj.
      CALL GetHcoState( HcoState )

      ! Get restart variable
      CALL HCO_RestartGet( am_I_Root, HcoState, 'STATE_PSC',  
     &                     STATE_PSC, RC, FILLED=IT_EXISTS ) 
      IF ( RC /= HCO_SUCCESS ) THEN
         MSG = 'Cannot get PSC restart variable'
         CALL ERROR_STOP ( MSG, 'READ_PSC_FILE (ucx_mod.F)' )
      ENDIF

      ! Cleanup
      HcoState => NULL()

      ! On first call, prompt some messages

      ! If restart exists, point to this field. Eventually deallocate
      ! pointer first.
      IF ( IT_EXISTS ) THEN
         IF (am_I_Root ) THEN
            WRITE(6,*) '    - UCX: Initialize PSC from restart' 
         ENDIF
      ELSE
         IF (am_I_Root ) THEN
         WRITE(6,*) 
     &      '    - UCX: PSC restart not found, initialize PSC-free'
         ENDIF
      ENDIF

      END SUBROUTINE READ_PSC_FILE
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_clock_trac
!
! !DESCRIPTION: Subroutine SET\_CLOCK\_TRAC sets the clock tracer mixing 
!  ratio within the bottom five grid levels, increasing by a fixed rate of
!  1.0 pptv/day
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CLOCK_TRAC( Input_Opt, STEPLEN, State_Chm )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE TIME_MOD,           ONLY : GET_JD, GET_NHMS, GET_NYMD
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TRACERID_MOD,       ONLY : IDTCLOCK
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)  :: STEPLEN     ! Step length (min)
      TYPE(OptInput), INTENT(IN)  :: Input_Opt  ! Input options
!
! !OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(OUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY: 
!  16 Oct 2013 - S. D. Eastham - Initial version
!  14 Feb 2014 - R. Yantosca   - Reorder DO loop for efficiency
!  21 Feb 2014 - M. Sulprizio  - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I,J,L,LCLOCK, NHMS, NYMD
      REAL(f8)          :: JDNow
      REAL(fp), PARAMETER :: INCRATE = 1.0e-12_fp ! vv/day increase
      
      !=================================================================
      ! SET_CLOCK_TRAC begins here!
      !=================================================================

      ! Get current date
      NHMS = GET_NHMS()
      NYMD = GET_NYMD()
      JDNOW = GET_JD(NYMD,NHMS)

      CLOCKMR = MAX(0.0e0_fp,(JDNOW - Input_Opt%JDClock) * INCRATE)

      If (Input_Opt%LPRT) Then
         Write(*,'(a,F16.3,a)') ' --> Setting clock tracer to ', 
     &      CLOCKMR*1e12, ' pptv'
      End If

      ! Figure out CLOCKMR based on date

      IF (CLOCKMR.lt.TINY(1e+0_fp)) THEN
         ! Reset clock tracer everywhere
         LCLOCK = LLPAR
      ELSE
         ! Only set bottom 5 levels
         LCLOCK = 5
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L=1,LCLOCK
      DO J=1,JJPAR
      DO I=1,IIPAR
         State_Chm%Tracers(I,J,L,IDTCLOCK) = CLOCKMR
      ENDDO ! L
      ENDDO ! I
      ENDDO ! J
!$OMP END PARALLEL DO

      ! Increment clock mixing ratio for next step
      !CLOCKMR = CLOCKMR + (INCRATE*STEPLEN/(60e+0_fp*24e+0_fp))

      END SUBROUTINE SET_CLOCK_TRAC
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_h2o_trac
!
! !DESCRIPTION: Subroutine SET\_H2O\_TRAC sets the H2O tracer throughout 
!  the selected domain (either troposphere only or the full grid).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_H2O_TRAC ( SETSTRAT, Input_Opt, State_Met,
     &                          State_Chm )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_TROP
      USE CMN_SIZE_MOD
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TRACERID_MOD,       ONLY : IDTH2O
      USE GRID_MOD,           ONLY : GET_YEDGE
      USE TIME_MOD,           ONLY : GET_DAY_OF_YEAR, GET_HOUR
      USE TIME_MOD,           ONLY : ITS_A_LEAPYEAR, GET_MINUTE
      USE CMN_GCTM_MOD,       ONLY : PI
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: SETSTRAT    ! Set strat H2O?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Initial version
!  14 Feb 2014 - R. Yantosca   - Reordered DO loop for efficiency
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm 
!                                objects via the arg list
!  24 Jun 2015 - S. D. Eastham - New boundary condition options
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I,J,L

      REAL(fp)            :: ESAT, SHMB, PRES, TEMP
      REAL(fp), PARAMETER :: A =  23.5518e+0_fp
      REAL(fp), PARAMETER :: B =  2937.4e+0_fp
      REAL(fp), PARAMETER :: C = -4.9283e+0_fp

      LOGICAL           :: READ_SPHU

      ! Local variables for quantities from Input_Opt
      LOGICAL           :: LACTIVEH2O

      ! CPT variables
      LOGICAL             :: LCPTH2O, LCPTLOCAL
      REAL(fp)            :: PMIN, H2OKG, H2OMR, YMAX, LATMAX, TMIN
      INTEGER             :: LCOLD
      REAL(fp)            :: TCVV(Input_Opt%N_TRACERS)
      INTEGER             :: N_TRACERS

      INTEGER,PARAMETER   :: CPTOPT = 3
      LOGICAL,PARAMETER   :: CPT_SAT = .FALSE.
      REAL(fp)            :: H2OBC, CPTMIN, CPTMAX
      REAL(fp)            :: YRFRAC, CPTRANGE, CPTMEAN

      CHARACTER(LEN=255)  :: CPTMSG

      ! Pointers
      REAL(fp), POINTER   :: STT(:,:,:,:)

      !=================================================================
      ! SET_H2O_TRAC begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      LACTIVEH2O = Input_Opt%LACTIVEH2O

      LCPTH2O     = ((Input_Opt%LCPTH2O.and.LACTIVEH2O)
     &                .and.(.not.SETSTRAT))
      CPTMIN      = Input_Opt%CPTH2OMIN
      CPTMAX      = Input_Opt%CPTH2OMAX
      N_TRACERS   = Input_Opt%N_TRACERS
      TCVV        = Input_Opt%TCVV(1:N_TRACERS)

      ! Fixed parameters
      PMIN        = 7.5e+1_fp
      LATMAX      = 30.e+0_fp
      YRFRAC      = 0e+0_fp

      ! Derived parameters
      CPTMEAN = (CPTMIN+CPTMAX)/2.e+0_fp
      CPTRANGE = (CPTMAX-CPTMIN)/2.e+0_fp

      ! Initialize GEOS-Chem tracer array [kg]
      STT => State_Chm%Tracers

      SELECT CASE (CPTOPT)
          CASE (0)
              ! Fixed to 3 ppmv
              H2OBC = 3.e-6_fp
          CASE (1)
              ! Fixed to 3.5 ppmv
              H2OBC = 3.5e-6_fp
          CASE (2)
              ! 150 day offset to get minimum between Dec/Mar
              ! Match Mccormack and Siskind (2002)
              YRFRAC = DBLE(GET_DAY_OF_YEAR())-155.e+0_fp
              ! Just to ensure a smooth oscillation...
              YRFRAC = YRFRAC + ((1.e+0_fp/24.e+0_fp)*
     &              (DBLE(GET_HOUR()) + (DBLE(GET_MINUTE())/60.e+0_fp)))
              IF (ITS_A_LEAPYEAR()) THEN
                  YRFRAC = YRFRAC/366.e+0_fp
              ELSE
                  YRFRAC = YRFRAC/365.e+0_fp
              ENDIF
              H2OBC = (3.5e+0_fp + SIN(2.e+0_fp*PI*YRFRAC))*1.e-6_fp
          CASE (3)
              ! Use input.geos with annual variation
              YRFRAC = DBLE(GET_DAY_OF_YEAR())-155.e+0_fp
              ! Just to ensure a smooth oscillation...
              YRFRAC = YRFRAC + ((1.e+0_fp/24.e+0_fp)*
     &              (DBLE(GET_HOUR()) + (DBLE(GET_MINUTE())/60.e+0_fp)))
              IF (ITS_A_LEAPYEAR()) THEN
                  YRFRAC = YRFRAC/366.e+0_fp
              ELSE
                  YRFRAC = YRFRAC/365.e+0_fp
              ENDIF
              H2OBC = (CPTMEAN + CPTRANGE*SIN(2.e+0_fp*PI*YRFRAC))*1.e-6_fp
          CASE DEFAULT
              H2OBC = 3.e-6_fp
      END SELECT

      IF (Input_Opt%LPRT) THEN
            WRITE(CPTMSG,'(a,I0.3,a,2(I0.2,a),F7.3)') 
     &        ' ### Time: ',GET_DAY_OF_YEAR(),
     &        ' (', GET_HOUR(),':',GET_MINUTE(),
     &        ' - fraction ', YRFRAC
            CALL DEBUG_MSG( TRIM(CPTMSG) )
            WRITE(CPTMSG,'(a,I0.3,a,F7.3,a)') 
     &        ' ### Setting CPT H2O for day ',GET_DAY_OF_YEAR(), ' to ',
     &        H2OBC*1.e+6_fp, ' ppmv'
            CALL DEBUG_MSG( TRIM(CPTMSG) )
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PRES, TEMP, ESAT, SHMB, READ_SPHU, LCOLD )
!$OMP+PRIVATE( YMAX, LCPTLOCAL, TMIN, H2OMR, H2OKG               )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      IF (LCPTH2O) THEN
         YMAX = MAX(ABS(GET_YEDGE(I,J,1)),ABS(GET_YEDGE(I,J+1,1)))
         IF (YMAX.LE.LATMAX) THEN
            LCOLD = -1
            LCPTLOCAL = .TRUE.
            L = 0
            TMIN = 300.e+0_fp
            PRES = 1.e+3_fp
            DO WHILE (PRES.gt.PMIN)
               L = L+1
               TEMP = State_Met%T(I,J,L)
               ! Is this the minimum temperature below 10 hPa?
               IF (TEMP.lt.TMIN) THEN
                  TMIN = TEMP
                  LCOLD = L
               ENDIF
               PRES = GET_PCENTER(I,J,L)
            ENDDO
         ELSE
            ! Just use tropopause definition
            LCPTLOCAL = .FALSE.
         ENDIF
      ELSE
         LCPTLOCAL = .FALSE.
      ENDIF
      ! L now has to be the inner loop
      DO L = 1, LLPAR
         ! Use hygropause parameterization?
         IF (LCPTLOCAL) THEN
             ! Find cold point
             READ_SPHU = (L.LT.LCOLD)
         ELSE
             READ_SPHU = ( ITS_IN_THE_TROP( I, J, L, State_Met ) .or.
     &                 SETSTRAT .or. ( .not. LACTIVEH2O )  )
         ENDIF

         IF ( READ_SPHU ) THEN
            STT(I,J,L,IDTH2O) = State_Met%SPHU(I,J,L) * 1.e-3_fp *
     &                          State_Met%AD(I,J,L)
         ELSE

            ! Pressure at midpoint of box (I,J,L)
            PRES = GET_PCENTER(I,J,L)

            ! Temperature at grid box (I,J,L)
            TEMP = State_Met%T(I,J,L)

            ! Saturation water vapor pressure in mbar 
            ! (from NASA GTE PEM-Tropics handbook)
            ESAT = ( 10.e+0_fp**( A - ( B / TEMP ) ) ) * ( TEMP**C )
            
            IF (LCPTLOCAL.and.(L.eq.LCOLD)) THEN
               IF (CPT_SAT) THEN
                  ! Prevent supersaturation (assume frozen out)
                   H2OMR = MIN(H2OBC,ESAT/PRES)
               ELSE
                   ! Set to fixed MR
                   H2OMR = H2OBC
               ENDIF

               H2OKG = H2OMR * State_Met%AD(I,J,L) / TCVV(IDTH2O)
               STT(I,J,L,IDTH2O) = H2OKG
            ENDIF

            ! Set specific humidity to transported H2O
            State_Met%SPHU(I,J,L) = STT(I,J,L,IDTH2O) * 1.e+3_fp /
     &                              State_Met%AD(I,J,L)

            ! Specific humidity in mb
            SHMB = State_Met%SPHU(I,J,L) * 1.6072e-3_fp * PRES
               
            ! Relative humidity as a percentage
            State_Met%RH(I,J,L) = ( SHMB / ESAT ) * 100e+0_fp 
         ENDIF
      ENDDO ! L
      ENDDO ! I
      ENDDO ! J
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE SET_H2O_TRAC
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ucx_h2so4phot
!
! !DESCRIPTION: Subroutine UCX\_H2SO4PHOT propagates the calculated H2SO4
!  photolysis (J) rate at the top of the chemistry grid through to the top
!  of the transport grid, approximating H2SO4 photolysis in the mesosphere.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE UCX_H2SO4PHOT( Input_Opt, State_Met, State_Chm )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : GET_CHEMGRID_LEVEL
      USE COMODE_LOOP_MOD,    ONLY : IRM, NKSO4PHOT, NCSGAS
      USE COMODE_LOOP_MOD,    ONLY : DEFPRAT, NAMEGAS, NRATES
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTSO2,IDTSO4g
      USE FAST_JX_MOD,        ONLY : FJXFUNC
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object

!
! !REMARKS:
! (1) A remark
! 
! !REVISION HISTORY: 
!  17 Aug 2013 - S. D. Eastham - Initial version
!  14 Feb 2014 - R. Yantosca   - Reorder DO loops for efficiency
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm
!                                objects via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: I, J, L, IJWINDOW
      REAL(fp)                :: GMU,SO4_IN,PHOTDELTA,DTCHEM
      REAL(fp)                :: SO4_DELTA
      LOGICAL               :: DAYCOLUMN
      INTEGER               :: LMINPHOT
      INTEGER,SAVE          :: LOCAL_NKSO4PHOT
      CHARACTER(fp),SAVE      :: LOCAL_NAME
      INTEGER,SAVE          :: LOCAL_IPHOT
      INTEGER,SAVE          :: IFNC,IBRCH
      REAL(fp),SAVE           :: RELWT
      LOGICAL,SAVE          :: FIRST=.TRUE.
      INTEGER               :: ICS

      ! Local variables for quantities from Input_Opt
      INTEGER               :: N_TRACERS
      REAL(fp)                :: TRACER_MW_G(Input_Opt%N_TRACERS)

      ! Pointers
      REAL(fp), POINTER       :: STT(:,:,:,:)

      !=================================================================
      ! UCX_H2SO4PHOT begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      N_TRACERS  = Input_Opt%N_TRACERS
      TRACER_MW_G = Input_Opt%TRACER_MW_G(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg]
      STT => State_Chm%Tracers

      IF (FIRST) THEN
         FIRST = .FALSE.
         DO ICS = 1,NCSGAS
            IF (NKSO4PHOT(ICS).ne.0) THEN
               LOCAL_NKSO4PHOT = NKSO4PHOT(ICS)
               LOCAL_NAME = NAMEGAS(IRM(1,LOCAL_NKSO4PHOT,ICS))
               LOCAL_IPHOT = LOCAL_NKSO4PHOT - NRATES(ICS)
               IFNC = DEFPRAT(LOCAL_NKSO4PHOT,ICS) + 0.01e+0_fp
               IBRCH = 10.e+0_fp*(DEFPRAT(LOCAL_NKSO4PHOT,ICS)-IFNC) + 0.5e+0_fp
               RELWT = TRACER_MW_G(IDTSO2)/TRACER_MW_G(IDTSO4g)
            ENDIF
         ENDDO 
      ENDIF

      ! Allow for the possibility of variable timestep
      DTCHEM = GET_TS_CHEM()*60.e+0_fp

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, IJWINDOW, DAYCOLUMN, GMU )
!$OMP+PRIVATE( LMINPHOT, PHOTDELTA, SO4_IN, SO4_DELTA )
      DO J=1,JJPAR
      DO I=1,IIPAR

         IJWINDOW  = (J-1)*IIPAR + I
         GMU       = State_Met%SUNCOSmid(I,J)
         DAYCOLUMN = (GMU.gt.0e+0_fp)

         IF (DAYCOLUMN) THEN
            LMINPHOT  = GET_CHEMGRID_LEVEL( I, J, State_Met)

            ! Retrieve photolysis rate as a fraction of gaseous SO4
            PHOTDELTA = FJXFUNC( I, J, LMINPHOT, LOCAL_IPHOT,
     &                           IBRCH, LOCAL_NAME) * DTCHEM
            PHOTDELTA = MIN(1.e+0_fp,PHOTDELTA)

            DO L=LMINPHOT+1,LLPAR
               ! Apply photolysis to SO4
               ! First retrieve gaseous fraction
               SO4_IN = STT(I,J,L,IDTSO4g)
               SO4_DELTA = PHOTDELTA*SO4_IN
               ! Remove from SO4
               STT(I,J,L,IDTSO4g) = STT(I,J,L,IDTSO4g) - SO4_DELTA
               ! Add to SO2. Note change in molar mass
               STT(I,J,L,IDTSO2) = STT(I,J,L,IDTSO2) + (SO4_DELTA*RELWT)
            ENDDO

         ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE UCX_H2SO4PHOT
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: noxcoeff_init
!
! !DESCRIPTION: Subroutine NOXCOEFF\_INIT initializes the NOX 2D interpolation
! values. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NOXCOEFF_INIT ( am_I_Root, Input_Opt )
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
      INTEGER            :: I, AS, IOS
      INTEGER            :: IMON, ITRAC, ILEV
      INTEGER            :: IU_FILE
      CHARACTER(LEN=255) :: NOX_FILE
      CHARACTER(LEN=255) :: TARG_TRAC
      CHARACTER(LEN=255) :: DBGMSG
      LOGICAL            :: LPRT
      REAL(fp)           :: TEMP_225(91,1)

      !=================================================================
      ! NOXCOEFF_INIT begins here!
      !=================================================================

      IF (.not.UCXNETCDF) THEN
         ! Copy fields from INPUT_OPT
         LPRT = Input_Opt%LPRT

         ! Number of latitude levels of NOXCOEFF array. Use 2x25
         ! and map NOx coeffs onto simulation grid when calling 
         ! GET_NOXCOEFF.
         JJNOXCOEFF = 91

         ! Fill NOx latitudes
         ALLOCATE(NOXLAT(JJNOXCOEFF+1), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOXLAT' )

         ! Fill manually
         NOXLAT(1) = -90.5e+0_fp
         NOXLAT(2) = -89.0e+0_fp
         DO I = 3,JJNOXCOEFF
            NOXLAT(I) = NOXLAT(I-1) + 2e+0_fp
         ENDDO

         ! Overshoot to make sure that a latitude of 90.0 will be properl
         ! matched onto JJNOXCOEFF. 
         NOXLAT(JJNOXCOEFF+1) = 90.5e+0_fp 

         ! Initialize the NOXCOEFF array. This array holds monthly NOx 
         ! coefficients on 51 levels and for 6 species.
         ALLOCATE(NOXCOEFF(JJNOXCOEFF,UCX_NLEVS,6,12), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOXCOEFF' )
         NOXCOEFF = 0.0e+0_fp

         ! Fill array
         DO IMON  = 1,12 
         DO ITRAC = 1,6
            SELECT CASE (ITRAC)
               CASE ( 1 )
                  TARG_TRAC = 'O'
               CASE ( 2 )
                  TARG_TRAC = 'O1D'
               CASE ( 3 )
                  TARG_TRAC = 'JNO'
               CASE ( 4 )
                  TARG_TRAC = 'JNO2'
               CASE ( 5 )
                  TARG_TRAC = 'JNO3'
               CASE ( 6 )
                  TARG_TRAC = 'JN2O'
            END SELECT
            WRITE(NOX_FILE,'(a,a,a,I0.2,a)') TRIM(NOON_FILE_ROOT),
     &       TRIM(TARG_TRAC), '_', IMON, '.dat'

            ! Get a free LUN
            IU_FILE = findFreeLUN()

            IOS = 1
            OPEN( IU_FILE,FILE=TRIM(NOX_FILE),STATUS='OLD',IOSTAT=IOS) 
            IF ( IOS /= 0 ) THEN
               WRITE(6,*) 'UCX_MOD: Could not read ', TRIM(NOX_FILE)
               CALL IOERROR( IOS, IU_FILE,'UCX_MOD:NOXCOEFF_INIT')
            ENDIF
  
            IF (LPRT) THEN
               WRITE(DBGMSG,'(a,a)') ' ### UCX: Reading ', 
     &           TRIM( NOX_FILE )
               CALL DEBUG_MSG( TRIM(DBGMSG) )
            ENDIF
 
            ! Read in data
            DO ILEV = 1,UCX_NLEVS
               READ(IU_FILE, 110, IOSTAT=IOS ) TEMP_225(:,1)
               ! Regrid from 2x2.5 onto target grid
               NOXCOEFF(:,ILEV,ITRAC,IMON) = 
     &            MATMUL(UCX_REGRID_225,TEMP_225(:,1))
               IF ( IOS /= 0 ) THEN
                  WRITE(6,'(a,a,I4,a,1x,a)') 'UCX_MOD: Error reading '
     &    , 'line ', ILEV, ' in file ', TRIM( NOX_FILE )
                  CALL IOERROR( IOS, IU_FILE,'UCX_MOD:NOXCOEFF_INIT')
               ENDIF
            ENDDO

110   FORMAT(91E10.3)

            CLOSE(IU_FILE)

         ENDDO !ITRAC
         ENDDO !IMON
      ENDIF

      END SUBROUTINE NOXCOEFF_INIT
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_jjnox 
!
! !DESCRIPTION: Subroutine GET\_JJNOX maps grid box at location IISIM, JJSIM of  
! the simulation grid onto the latitude grid of the NOXCOEFF array. JJNOX can
! differ from JJSIM if it's not a 4x5 or 2x25 simulation.
!\\
!\\
! This routine simply returns the index of the NOx latitude vector that covers
! the latitude value of interest. No grid box weighting, etc. is performed.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_JJNOX( IISIM, JJSIM ) RESULT ( JJNOX )
!
! !USES:
!
      USE GRID_MOD,           ONLY : GET_YMID
!
! !INPUT PARAMETERS:
!
      INTEGER,          INTENT(IN) :: IISIM    ! Latitude index on simulation grid 
      INTEGER,          INTENT(IN) :: JJSIM    ! Latitude index on simulation grid 
!
! !OUTPUT VARIABLES:
!
      INTEGER                      :: JJNOX    ! Latitude index on NOXCOEFF grid
!
! !REVISION HISTORY: 
!  05 Dec 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I
      REAL(fp)              :: LAT 

      !=================================================================
      ! GET_JJNOX begins here!
      !=================================================================

      ! Nothing to do for 'standard' grids
#if defined( GRID4x5 ) || defined( GRID2x25 )
      JJNOX = JJSIM
#else

      ! Init
      JJNOX = -1

      ! Get latitude in degrees north on simulation grid
      LAT = GET_YMID ( IISIM, JJSIM, 1 )

      ! Loop over all latitudes of the NOx grid until we reach the grid
      ! box where the simulation latitude sits in. 
      I = 0
      DO WHILE (JJNOX < 0)
         I = I + 1
         IF ( LAT < NOXLAT(I+1) ) THEN
            JJNOX = I
         ENDIF
         IF (I.GT.JJNOXCOEFF) THEN
            JJNOX = JJNOXCOEFF * 2
         ENDIF
      ENDDO

      !=====SDE DEBUG
      WRITE(6,'(a,2(I2,x),F10.3,x,I3)') 'SDE DEBUG: II JJ LAT JJNOX: ',
     &                              IISIM, JJSIM, LAT, JJNOX
      !=====SDE DEBUG
#endif

      END FUNCTION GET_JJNOX
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_ucx
!
! !DESCRIPTION: Subroutine INIT\_UCX initializes module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_UCX( am_I_Root, Input_Opt )
      !SUBROUTINE INIT_UCX( am_I_Root, Input_Opt, YYYYMMDD, HHMMSS )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE ERROR_MOD,          ONLY : IS_SAFE_DIV
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GRID_MOD,           ONLY : GET_YEDGE
      USE TIME_MOD,           ONLY : GET_YEAR
      USE TRACERID_MOD
      USE CMN_GCTM_MOD,       ONLY : PI
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!      INTEGER,        INTENT(IN) :: YYYYMMDD    ! YYYY/MM/DD GMT date
!      INTEGER,        INTENT(IN) :: HHMMSS      ! hh:mm:ss   GMT time
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
      CHARACTER(LEN=255) :: DBGMSG, GRIDSPEC
      INTEGER            :: JIN, JOUT
      REAL(fp)             :: JMIN_IN, JMAX_IN, JDIF_IN
      REAL(fp)             :: JMIN_OUT,JMAX_OUT,JDIF_OUT
      REAL(fp)             :: JMIN_TMP,JMAX_TMP,JDIF_TMP
      REAL(fp)             :: JRATIO
      REAL(fp)             :: DEG_SUM

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
      LOGICAL :: LSETOCS
      LOGICAL :: LSETCFC
      LOGICAL :: LSETCL
      LOGICAL :: LSETBR
      LOGICAL :: LSETH2SO4
      LOGICAL :: LSETBRSTRAT
      LOGICAL :: LSETNOYSTRAT
      LOGICAL :: LSETN2O
!      LOGICAL :: LBASICEMIS

      ! Return code
      INTEGER :: RC

      !=================================================================
      ! INIT_UCX begins here!
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
      LSETOCS      = Input_Opt%LSETOCS
      LSETCFC      = Input_Opt%LSETCFC
      LSETCL       = Input_Opt%LSETCL
      LSETBR       = Input_Opt%LSETBR
      LSETH2SO4    = Input_Opt%LSETH2SO4
      LSETBRSTRAT  = Input_Opt%LSETBRSTRAT
      LSETNOYSTRAT = Input_Opt%LSETNOYSTRAT
      LSETN2O      = Input_Opt%LSETN2O
!      LBASICEMIS   = Input_Opt%LBASICEMIS

      ! Set module variable doPrint for use in internal functions
      doPrint      = LPRT .AND. am_I_Root

      IF ( am_I_Root ) THEN
         WRITE( 6,'(a)') REPEAT( '=', 79 )
         WRITE( 6,'(a)') 'U N I F I E D   C H E M I S T R Y'
         WRITE( 6,'(a)') 'Routines written by SEBASTIAN D. EASTHAM'
         WRITE( 6,'(a)') REPEAT( '=', 79 )
      ENDIF

      ! Don't use BASICEMIS switch in check below (ckeller, 05/12/14)
      USE2DDATA = ( LSETOCS     .or. LSETCFC      .or. LSETCL     .or.
     &              LSETBR      .or. LSETH2SO4    .or. STRAT2DCH4 .or.
     &              LSETBRSTRAT .or. LSETNOYSTRAT .or. LSETN2O        )

      ! --------------------------------------------------------------
      ! Input data sources
      ! --------------------------------------------------------------
      
      ! Determine folder paths from root folder
      IF (UCXNETCDF) THEN
         WRITE( AVG_FILE_ROOT,'(a,a)') 
     &        TRIM(Input_Opt%CHEM_INPUTS_DIR),
     &      'UCX_201403/Init2D/DiurnalAvg.nc'
         WRITE( NOON_FILE_ROOT,'(a,a)') 
     &        TRIM(Input_Opt%CHEM_INPUTS_DIR), 
     &        'UCX_201403/Init2D/Noontime.nc'

      ! For ASCII input, use 2x25 grid for all other grids than 4x5. 
      ! This is ok for the NOx coeffs which can be regridded on the fly
      ! from 2x25 onto any other grid. This won't work for the 2D 
      ! boundary conditions, but those have been checked in the logical 
      ! check above (USE2DDATA). 
      ELSE
         ! SDE 2015-08-18: Now always regrid 2x2.5
         GRIDSPEC = 'Grid2x25/InitCFC_'
         WRITE( AVG_FILE_ROOT,'(a,a,a)') 
     &        TRIM(Input_Opt%CHEM_INPUTS_DIR),
     &      'UCX_201403/DiurnalAvg/', TRIM(GRIDSPEC)
         WRITE(   NOON_FILE_ROOT,'(a,a,a)') 
     &        TRIM(Input_Opt%CHEM_INPUTS_DIR), 
     &        'UCX_201403/NoonTime/', TRIM(GRIDSPEC)
      ENDIF
  
      IF ( LPRT ) THEN
         WRITE(DBGMSG,'(a,a)') '### UCX: Reading O1D/O3P from ',
     &     TRIM(NOON_FILE_ROOT)
         CALL DEBUG_MSG( TRIM(DBGMSG) )
         WRITE(DBGMSG,'(a,a)') 
     &     '### UCX: Initializing long-lived species from ',
     &     TRIM(AVG_FILE_ROOT)
         CALL DEBUG_MSG( TRIM(DBGMSG) )
         WRITE(DBGMSG,'(a,a)') 
      ENDIF

      ! Allocate arrays of input pressure levels and lat edges
      ALLOCATE( UCX_PLEVS( UCX_NLEVS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'UCX_PLEVS' )

      ALLOCATE( UCX_LATS( UCX_NLAT+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'UCX_LATS' )

      ! Set input pressure levels (hPa)
      UCX_PLEVS = (/ 0.2200e+00_fp, 0.2600e+00_fp, 
     &               0.3100e+00_fp, 0.3600e+00_fp,
     &               0.4300e+00_fp, 0.5100e+00_fp, 
     &               0.6000e+00_fp, 0.7100e+00_fp,
     &               0.8400e+00_fp, 0.9900e+00_fp, 
     &               1.1700e+00_fp, 1.3800e+00_fp,
     &               1.6300e+00_fp, 1.9300e+00_fp, 
     &               2.2800e+00_fp, 2.6900e+00_fp,
     &               3.1800e+00_fp, 3.7600e+00_fp, 
     &               4.4400e+00_fp, 5.2500e+00_fp,
     &               6.2000e+00_fp, 7.3200e+00_fp, 
     &               8.6500e+00_fp, 1.0220e+01_fp,
     &               1.2070e+01_fp, 1.4260e+01_fp, 
     &               1.6850e+01_fp, 1.9910e+01_fp,
     &               2.3520e+01_fp, 2.7780e+01_fp, 
     &               3.2820e+01_fp, 3.8770e+01_fp,
     &               4.5810e+01_fp, 5.4110e+01_fp, 
     &               6.3930e+01_fp, 7.5220e+01_fp,
     &               8.9220e+01_fp, 1.0540e+02_fp, 
     &               1.2451e+02_fp, 1.4710e+02_fp,
     &               1.7377e+02_fp, 2.0529e+02_fp, 
     &               2.4252e+02_fp, 2.8650e+02_fp,
     &               3.3847e+02_fp, 3.9985e+02_fp, 
     &               4.7237e+02_fp, 5.5804e+02_fp,
     &               6.5924e+02_fp, 7.7880e+02_fp, 
     &               9.2004e+02_fp /)

      ! Set input latitude edges (degrees)
      UCX_LATS(1)  = -90.0e+0_fp
      UCX_LATS(2)  = (-90.0e+0_fp) + (9.5e+0_fp/2.0e+0_fp)
      UCX_LATS(20) = 90.0e+0_fp
      DO N=2,18
         UCX_LATS(N+1) = UCX_LATS(N) + 9.5e+0_fp
      ENDDO

      ! Calculate conversion factors for SLA
      ! Factor to convert volume (m3 SLA/m3 air) to
      ! surface area density (cm2 SLA/cm3 air)
      SLA_VA = (8.406e-8_fp)*(10.e+0_fp**(12.e+0_fp*0.751e+0_fp))

      ! Factor to convert effective radius to 
      ! liquid radius (unitless)
      SLA_RR = EXP(-0.173e+0_fp)

      ! Factor to convert volume (m3/m3) to effective
      ! radius (m)
      SLA_VR = (0.357e-6_fp)*(10.e+0_fp**(12.e+0_fp*0.249))

      ! Conversion factors (unitless) - median -> effective radius
      ! Reff = Rmed * LN_RMED_TO_REFF(I_S?A)
      IF (Input_Opt%LWACCMAER) THEN
         ! Use WACCM standard parameters for now
         NWACCM(I_SLA)    = 10.00e+6_fp
         NWACCM(I_SPA)    = 00.01e+6_fp
         LSSQWACCM(I_SLA) = (DLOG(1.6e+0_fp))**2.0e+0_fp
         LSSQWACCM(I_SPA) = (DLOG(1.6e+0_fp))**2.0e+0_fp

         ! Set up parameters for the Chen et al SNAP scheme
         CHENA1 = 1.43089d0
         CHENA2 = 1.0295d0
         DO N=1,2
            ! (LN_RMED_TO_REFF * Rmed) = Effective (A-weighted) radius
            LN_RMED_TO_REFF(N)= DEXP(2.5e+0_fp * LSSQWACCM(N))
            ! (LN_RMED_TO_RLIQ * Rmed) = Volume-weighted radius
            LN_RMED_TO_RLIQ(N)= DEXP(3.5e+0_fp * LSSQWACCM(N))
            ! (LN_VCONV * V0)^(1/3)    = Median radius
            LN_VCONV(N)       = 0.75e+0_fp*DEXP(-4.5e+0_fp*LSSQWACCM(N))
     &                          /(PI * NWACCM(N))
            ! LN_ACONV * (Rmed^2)      = SAD
            LN_ACONV(N)       = 4.0e+0_fp*PI*NWACCM(N)*
     &                          DEXP(2.0e+0_fp*LSSQWACCM(N))
         ENDDO
      ELSE
         ! Grainger SLA is basically log-normal with sigma ~ 1.3
         LN_RMED_TO_REFF(I_SLA) = 1.0e+0_fp/SLA_RR
         LN_RMED_TO_REFF(I_SPA) = 1.0e+0_fp
         LN_RMED_TO_RLIQ(I_SLA) = 1.0e+0_fp
         LN_RMED_TO_RLIQ(I_SPA) = 1.0e+0_fp
         ! Ignore the others
         LN_ACONV = 0.0e+0_fp
         LN_VCONV = 0.0e+0_fp
      ENDIF

      ! Set clock tracer mixing ratio to zero
      CLOCKMR = 0e+0_fp

      ! Initialize NOx coefficient arrays
      ALLOCATE( NOX_O( IIPAR, JJPAR, LLPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOX_O' )
      NOX_O = 0e+0_fp

      ALLOCATE( NOX_J( IIPAR, JJPAR, LLPAR, 4 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOX_J' )
      NOX_J = 0e+0_fp

      ! Initialize PSC variables
      ALLOCATE( RAD_AER( IIPAR, JJPAR, LLPAR, NSTRATAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RAD_AER' )
      RAD_AER = 0e+0_fp
      
      ALLOCATE( KG_AER( IIPAR, JJPAR, LLPAR, NSTRATAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'KG_AER' )
      KG_AER = 0e+0_fp
      
      ALLOCATE( SAD_AER( IIPAR, JJPAR, LLPAR, NSTRATAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SAD_AER' )
      SAD_AER = 0e+0_fp
      
      ALLOCATE( NDENS_AER( IIPAR, JJPAR, LLPAR, NSTRATAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NDENS_AER' )
      NDENS_AER = 0e+0_fp
      
      ALLOCATE( RHO_AER( IIPAR, JJPAR, LLPAR, NSTRATAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RHO_AER' )
      RHO_AER = 0e+0_fp
      
      ALLOCATE( STATE_PSC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'STATE_PSC' )
      !STATE_PSC = 0e+0_fp
      STATE_PSC     = 0.0e+0_f4
      
      ! Mass fraction of tracers contained in liquid aerosol
      ! Indices: 1 - SO4
      !          2 - HNO3
      !          3 - HCl
      !          4 - HOCl
      !          5 - HBr
      !          6 - HOBr
      !          7 - H2O
      ALLOCATE( AERFRAC( IIPAR, JJPAR, LLPAR, 7 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AERFRAC' )
      AERFRAC = 0e+0_fp

      ALLOCATE( AERFRACIND( 7 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AERFRACIND' )
      AERFRACIND(1) = IDTSO4
      AERFRACIND(2) = IDTHNO3
      AERFRACIND(3) = IDTHCl
      AERFRACIND(4) = IDTHOCl
      AERFRACIND(5) = IDTHBr
      AERFRACIND(6) = IDTHOBr
      AERFRACIND(7) = IDTH2O
    
      ! Sticking coefficients for PSC reactions on SLA 
      ALLOCATE( KHETI_SLA( IIPAR, JJPAR, LLPAR, 11 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'KHETI_SLA' )
      KHETI_SLA = 0e+0_fp
     
      ! H2SO4 photolysis rate at the top of the chemgrid
      ALLOCATE( SO4_TOPPHOT( IIPAR,JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4_TOPPHOT' )
      SO4_TOPPHOT = 0.e+0_fp

      ALLOCATE( UCX_REGRID( JGLOB, UCX_NLAT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'UCX_REGRID' )
      UCX_REGRID = 0e+0_fp

      IF (.not.UCXNETCDF) THEN
         ALLOCATE( UCX_REGRID_225( JGLOB, 91 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'UCX_REGRID_225' )
         UCX_REGRID_225 = 0e+0_fp
      ENDIF

      ! Calculate the scaling matrices
      ! Note cosine (area-weighted)
      JMAX_OUT = GET_YEDGE(1,1,1)
      DO JOUT=1,JJPAR
         JMIN_OUT = JMAX_OUT
         JMAX_OUT = GET_YEDGE(1,JOUT+1,1)
         JDIF_OUT = SIND(JMAX_OUT)-SIND(JMIN_OUT)
         DEG_SUM = 0e+0_fp
         DO JIN=1,UCX_NLAT
            JMIN_IN = UCX_LATS(JIN)
            JMAX_IN = UCX_LATS(JIN+1)
            IF ((JMAX_OUT.ge.JMIN_IN).and.(JMIN_OUT.le.JMAX_IN)) THEN
               JMAX_TMP = MIN(JMAX_IN,JMAX_OUT)
               JMIN_TMP = MAX(JMIN_IN,JMIN_OUT)
               JDIF_TMP = SIND(JMAX_TMP)-SIND(JMIN_TMP)
               IF (IS_SAFE_DIV(JDIF_TMP,JDIF_OUT)) THEN
                  JRATIO = JDIF_TMP/JDIF_OUT
                  UCX_REGRID(JOUT,JIN) = JRATIO
                  DEG_SUM = DEG_SUM + JRATIO
               ENDIF
            ENDIF
         ENDDO
         ! Normalize
         IF (DEG_SUM.gt.0e+0_fp) THEN
            DO JIN=1,UCX_NLAT
               UCX_REGRID(JOUT,JIN) = UCX_REGRID(JOUT,JIN)/DEG_SUM
            ENDDO
         ELSE
            UCX_REGRID(JOUT,:) = 0e+0_fp
         ENDIF

         IF (.not.UCXNETCDF) THEN
            DEG_SUM = 0e+0_fp
            DO JIN=1,91
               ! Deal with half-polar input grid
               IF (JIN.eq.1) THEN
                  JMIN_IN = -90.0e+0_fp
                  JMAX_IN = -89.0e+0_fp
               ELSEIF (JIN.eq.91) THEN
                  JMIN_IN = 89.0e+0_fp
                  JMAX_IN = 90.0e+0_fp
               ELSE
                  JMIN_IN = -91.0e+0_fp + (JIN * 2.0e+0_fp)
                  JMAX_IN = JMIN_IN + 2.0e+0_fp
               ENDIF
               IF ((JMAX_OUT.ge.JMIN_IN).and.(JMIN_OUT.le.JMAX_IN)) THEN
                  JMAX_TMP = MIN(JMAX_IN,JMAX_OUT)
                  JMIN_TMP = MAX(JMIN_IN,JMIN_OUT)
                  JDIF_TMP = SIND(JMAX_TMP)-SIND(JMIN_TMP)
                  IF (IS_SAFE_DIV(JDIF_TMP,JDIF_OUT)) THEN
                     JRATIO = JDIF_TMP/JDIF_OUT
                     UCX_REGRID_225(JOUT,JIN) = JRATIO
                     DEG_SUM = DEG_SUM + JRATIO
                  ENDIF
               ENDIF
            ENDDO
            ! Normalize
            IF (DEG_SUM.gt.0e+0_fp) THEN
               DO JIN=1,UCX_NLAT
                  UCX_REGRID_225(JOUT,JIN) = 
     &                      UCX_REGRID_225(JOUT,JIN)/DEG_SUM
               ENDDO
            ELSE
               UCX_REGRID_225(JOUT,:) = 0e+0_fp
            ENDIF
         ENDIF
      ENDDO

      ! Initialize NOXCOEFF arrays
      IF (.not.UCXNETCDF) THEN
         CALL NOXCOEFF_INIT( am_I_Root, Input_Opt ) 
      ENDIF

      ! Initialize diagnostics
      CALL DIAGINIT_UCX ( am_I_Root, Input_Opt, RC )

      END SUBROUTINE INIT_UCX
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_ucx
!
! !DESCRIPTION: Subroutine DIAGINIT\_UCX initializes diagnostics containers
! for the UCX module. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAGINIT_UCX( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE HCO_RESTART_MOD,    ONLY : HCO_RestartDefine
      USE HCO_STATE_MOD,      ONLY : HCO_STATE
      USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoState
      USE HCO_ERROR_MOD
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY: 
!  14 Jan 2015 - C. Keller     - Initial version
!  25 Mar 2015 - C. Keller     - Now use HEMCO restart module
!EOP
!------------------------------------------------------------------------------
!BOC
      TYPE(HCO_State), POINTER :: HcoState => NULL()
      CHARACTER(LEN=255)       :: MSG
      CHARACTER(LEN=255)       :: LOC = 'DIAGINIT_UCX (ucx_mod.F)' 

      !=================================================================
      ! DIAGINIT_UCX begins here!
      !=================================================================

      ! Define restart quantity    
      CALL GetHcoState( HcoState ) 
      CALL HCO_RestartDefine ( am_I_Root, HcoState, 
     &                       'STATE_PSC', STATE_PSC, 'count', RC )
      HcoState => NULL() 

      IF ( RC /= HCO_SUCCESS ) THEN
         MSG = 'Cannot define STATE_PSC restart'
         CALL ERROR_STOP ( MSG, LOC ) 
      ENDIF

      ! Return w/ success
      RC = GIGC_SUCCESS

      END SUBROUTINE DIAGINIT_UCX
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_state_psc
!
! !DESCRIPTION: Subroutine WRITE\_STATE\_PSC saves the STATE\_PSC array into
! diagnostics. This is only of relevance in an ESMF environment, where 
! STATE\_PSC is not automatically written into the HEMCO restart diagnostics
! but needs to be passed explicitly to the internal state. This should be
! done on every time step to make sure that replay simulations will always
! have the most current STATE\_PSC values in the internal state (important
! for checkpointing!).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_STATE_PSC ( am_I_Root, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE HCO_ERROR_MOD
      USE HCO_RESTART_MOD,  ONLY : HCO_RestartWrite
      USE HCO_STATE_MOD,    ONLY : HCO_STATE
      USE HCOI_GC_MAIN_MOD, ONLY : GetHcoState
! 
! !INPUT ARGUMENTS:
!
      LOGICAL, INTENT(IN)     :: am_I_Root
! 
! !INPUT/OUTPUT ARGUMENTS:
!
      INTEGER, INTENT(INOUT)  :: RC
!
! !REVISION HISTORY: 
!  08 May 2015 - C. Keller     - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      TYPE(HCO_STATE), POINTER  :: HcoState => NULL()

      !=================================================================
      ! WRITE_STATE_PSC begins here!
      !=================================================================

      ! Eventually write PSC_STATE to restart. 
      CALL GetHcoState ( HcoState )

      CALL HCO_RestartWrite( am_I_Root,   HcoState, 
     &                       'STATE_PSC', STATE_PSC, RC )
      HcoState => NULL()
      IF ( RC /= HCO_SUCCESS ) THEN
         RC = GIGC_FAILURE
         RETURN
      ENDIF

      ! Return w/ success
      RC = GIGC_SUCCESS

      END SUBROUTINE WRITE_STATE_PSC 
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_ucx
!
! !DESCRIPTION: Subroutine CLEANUP\_UCX deallocates module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_UCX ( am_I_Root )
! 
! !INPUT ARGUMENTS:
!
      LOGICAL, INTENT(IN)     :: am_I_Root
!
! !REVISION HISTORY: 
!  04 Apr 2013 - S. D. Eastham - Initial version
!  13 Apr 2013 - S. D. Eastham - Added PSC arrays
!EOP
!------------------------------------------------------------------------------
!BOC
      INTEGER :: N

      !=================================================================
      ! CLEANUP_UCX begins here!
      !=================================================================

      IF ( ALLOCATED( RAD_AER    ) ) DEALLOCATE( RAD_AER    )
      IF ( ALLOCATED( SAD_AER    ) ) DEALLOCATE( SAD_AER    )
      IF ( ALLOCATED( KG_AER     ) ) DEALLOCATE( KG_AER     )
      IF ( ALLOCATED( RHO_AER    ) ) DEALLOCATE( RHO_AER    )
      IF ( ASSOCIATED( STATE_PSC ) ) DEALLOCATE( STATE_PSC  )
      IF ( ALLOCATED( NDENS_AER  ) ) DEALLOCATE( NDENS_AER  )
      IF ( ALLOCATED( AERFRAC    ) ) DEALLOCATE( AERFRAC    )
      IF ( ALLOCATED( AERFRACIND ) ) DEALLOCATE( AERFRACIND )
      IF ( ALLOCATED( KHETI_SLA  ) ) DEALLOCATE( KHETI_SLA  )
      IF ( ALLOCATED( UCX_MR     ) ) DEALLOCATE( UCX_MR     )
      IF ( ALLOCATED( UCX_MR_IN  ) ) DEALLOCATE( UCX_MR_IN  )
      IF ( ALLOCATED( UCX_REGRID ) ) DEALLOCATE( UCX_REGRID )
      IF ( ALLOCATED( UCX_REGRID_225 ) ) DEALLOCATE( UCX_REGRID_225 )
      IF ( ALLOCATED( UCX_PLEVS  ) ) DEALLOCATE( UCX_PLEVS  )
      IF ( ALLOCATED( UCX_LATS   ) ) DEALLOCATE( UCX_LATS   )
      IF ( ALLOCATED( NOX_O      ) ) DEALLOCATE( NOX_O      )
      IF ( ALLOCATED( NOX_J      ) ) DEALLOCATE( NOX_J      )
      IF ( ALLOCATED( SO4_TOPPHOT) ) DEALLOCATE( SO4_TOPPHOT)

      ! Cleanup the NOx coeff arrays 
      IF ( ALLOCATED( NOXCOEFF ) ) DEALLOCATE( NOXCOEFF )
      IF ( ALLOCATED( NOXLAT   ) ) DEALLOCATE( NOXLAT   )

      END SUBROUTINE CLEANUP_UCX
!EOC
      END MODULE UCX_MOD
