!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag_moss_flux
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
      MODULE DIAG_MASS_FLUX_QUS
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: INIT_MASS_FLUX
      PUBLIC  :: CLEANUP_MASS_FLUX
      PUBLIC  :: WRITE_DIAG_MASS_FLUX
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !PUBLIC DATA MEMBERS:
      REAL*8, PUBLIC, ALLOCATABLE :: DIAG_PU(:,:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: DIAG_PU_TR(:,:,:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: DIAG_PV(:,:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: DIAG_PV_TR(:,:,:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: DIAG_SD(:,:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: DIAG_SD_TR(:,:,:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: DIAG_STE(:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: DIAG_STE_TR(:,:,:)
      REAL*8, PUBLIC :: DIAG_MFLUX_CNT, DIAG_MFLUX_TR_CNT
      LOGICAL, PUBLIC, PARAMETER :: DIAG_MFLUX = .TRUE.
!
! !REMARKS:
!
!  References:
!  ============================================================================
!  (1 )
! !REVISION HISTORY:
!  1 Feb 2011 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Scalars
      INTEGER, SAVE        :: fID
      INTEGER, SAVE        :: NRXNS
      INTEGER, SAVE        :: tindex
      REAL*8, SAVE         :: tstart

      ! Arrays

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diag_mass_flux
!
! !DESCRIPTION: Subroutine WRITE\_DIAG\_MASS\_FLUX
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_DIAG_MASS_FLUX( Input_Opt )
!
! !USES:
      USE TIME_MOD, ONLY   : GET_TS_CHEM, GET_TAU
      USE GIGC_Input_Opt_Mod, ONLY : OptInput

      USE m_netcdf_io_create
      USE m_netcdf_io_write

      USE CMN_SIZE_MOD

! !REVISION HISTORY: 
!  1 Feb 2011 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      REAL*8  :: dT
      INTEGER :: varID, I, J, K, N
      INTEGER :: N_TRACERS
      LOGICAL, SAVE :: FIRST
      INTEGER, SAVE :: timeIND

      IF ( FIRST ) timeIND = 0
      timeIND = timeIND + 1
      first = .false.
      N_TRACERS = Input_Opt%N_TRACERS

      ! Scale to get average over the diagnostic period
      IF ( DIAG_MFLUX_CNT .gt. 0 ) THEN
         DIAG_PU = DIAG_PU / DIAG_MFLUX_CNT
         DIAG_PV = DIAG_PV / DIAG_MFLUX_CNT
         DIAG_SD = DIAG_SD / DIAG_MFLUX_CNT
         DIAG_STE = DIAG_STE / DIAG_MFLUX_CNT
      ENDIF

      ! Write to file
      call NcWr( DIAG_PU, fID, 'pu',
     &           (/     1,     1,     1, tindex /),    ! Start
     &           (/ IIPAR, JJPAR, LLPAR,      1 /)  )  ! Count

      call NcWr( DIAG_PV, fID, 'pv',
     &           (/     1,     1,     1, tindex /),    ! Start
     &           (/ IIPAR, JJPAR, LLPAR,      1 /)  )  ! Count

      call NcWr( DIAG_SD, fID, 'sd',
     &           (/     1,     1,       1, tindex /),    ! Start
     &           (/ IIPAR, JJPAR, LLPAR-1,      1 /)  )  ! Count

      call NcWr( DIAG_STE, fID, 'ste_flux',
     &           (/     1,     1, tindex /),    ! Start
     &           (/ IIPAR, JJPAR,      1 /)  )  ! Count


      dT = ( GET_TAU() - tstart ) * 3600d0 ! Seconds

      ! Write to file and convert from [mb m2] to [kg s-1]
      call NcWr( DIAG_PU_TR*1d2/9.8d0/dT, fID, 'pu_tracer',
     &           (/     1,     1,     1,         1, tindex /),    ! Start
     &           (/ IIPAR, JJPAR, LLPAR, N_TRACERS,      1 /)  )  ! Count

      call NcWr( DIAG_PV_TR*1d2/9.8d0/dT, fID, 'pv_tracer',
     &           (/     1,     1,     1,         1, tindex /),    ! Start
     &           (/ IIPAR, JJPAR, LLPAR, N_TRACERS,      1 /)  )  ! Count

      call NcWr( DIAG_SD_TR*1d2/9.8d0/dT, fID, 'sd_tracer',
     &           (/     1,     1,       1,         1, tindex /),    ! Start
     &           (/ IIPAR, JJPAR, LLPAR-1, N_TRACERS,      1 /)  )  ! Count

      call NcWr( DIAG_STE_TR*1d2/9.8d0/dT, fID, 'ste_tracer_flux',
     &           (/     1,     1,         1, tindex /),    ! Start
     &           (/ IIPAR, JJPAR, N_TRACERS,      1 /)  )  ! Count

      ! Push to file
      call NcDo_Sync( fID )

      ! Reset diagnostic array and counters
      DIAG_PU        = 0d0
      DIAG_PV        = 0d0
      DIAG_SD        = 0d0
      DIAG_STE       = 0d0
      DIAG_MFLUX_CNT = 0d0

      DIAG_PU_TR        = 0d0
      DIAG_PV_TR        = 0d0
      DIAG_SD_TR        = 0d0
      DIAG_STE_TR       = 0d0
      DIAG_MFLUX_TR_CNT = 0d0

      tstart = GET_TAU()

      ! Advance time index
      tindex = tindex + 1

      END SUBROUTINE WRITE_DIAG_MASS_FLUX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_strat_chem
!
! !DESCRIPTION: Subroutine INIT\_STRAT\_CHEM allocates all module arrays.  
!  It also opens the necessary rate files.
!\\
!\\
! !INTERFACE:
!      
      SUBROUTINE INIT_MASS_FLUX( Input_Opt )
!
! !USES:
!
      USE GRID_MOD,    ONLY : GET_YMID, GET_XMID
      USE BPCH2_MOD,   ONLY : GET_NAME_EXT, GET_RES_EXT
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE TIME_MOD,    ONLY : GET_TS_CHEM, EXPAND_DATE
      USE TIME_MOD,    ONLY : GET_NYMDb, GET_NHMSb, GET_TAU
      USE GIGC_Input_Opt_Mod, ONLY : OptInput

      USE m_netcdf_io_define
      USE m_netcdf_io_create
      USE m_netcdf_io_write
      USE m_netcdf_io_close

      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD ! NTRATES

      IMPLICIT NONE

      ! netCDF include files

# include "netcdf.inc"

! 
! !REVISION HISTORY:
!  1 Feb 2011 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS
      CHARACTER(LEN=255) :: FILENAME, FILENAMEOUT
      INTEGER            :: N, NN, F, I, J, L
      INTEGER            :: varID, varID2, varID3, varID4, varID5,
     &                      dim1, dim2, dim3, dim4, dim5

      INTEGER :: omode, vID, var1(1), var2(2), var3(3), var4(4), var5(5)
      INTEGER :: idLon, idLat, idLev, idLevE, idLevI, idtime, idTr
      INTEGER :: strlen_14_dimid, strlen_256_dimid      

      INTEGER :: levels(LLPAR)
      REAL*8  :: ymid(JJPAR), xmid(IIPAR)

      INTEGER            :: NK, IREAC1
      LOGICAL            :: IS_DRYDEP_RXN, EQRXN, PMRXN
      CHARACTER(LEN=256) :: rxn_name(NTRATES(NCS)), rxn_string, rxn_prev
      CHARACTER(LEN=10)  :: rxn_coeff
      CHARACTER(LEN=14)  :: gasname
      LOGICAL            :: LLINOZ
      CHARACTER(LEN=255) :: DATA_DIR
      INTEGER            :: N_TRACERS

      LLINOZ = Input_Opt%LLINOZ
      DATA_DIR = TRIM(Input_Opt%DATA_DIR)
      N_TRACERS = Input_Opt%N_TRACERS

      ! Initialize counter
      DIAG_MFLUX_CNT = 0d0

      ! Allocate arrays
      ALLOCATE( DIAG_PU(IIPAR,JJPAR,LLPAR),STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DIAG_PU' )
      DIAG_PU = 0d0

      ALLOCATE( DIAG_PV(IIPAR,JJPAR,LLPAR),STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DIAG_PV' )
      DIAG_PV = 0d0

      ALLOCATE( DIAG_SD(IIPAR,JJPAR,LLPAR-1),STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DIAG_SD' )
      DIAG_SD = 0d0

      ALLOCATE( DIAG_STE(IIPAR,JJPAR),STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DIAG_STE' )
      DIAG_STE = 0d0

      ALLOCATE( DIAG_PU_TR(IIPAR,JJPAR,LLPAR,N_TRACERS),STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DIAG_PU_TR' )
      DIAG_PU_TR = 0d0

      ALLOCATE( DIAG_PV_TR(IIPAR,JJPAR,LLPAR,N_TRACERS),STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DIAG_PV_TR' )
      DIAG_PV_TR = 0d0

      ALLOCATE( DIAG_SD_TR(IIPAR,JJPAR,LLPAR-1,N_TRACERS),STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DIAG_SD_TR' )
      DIAG_SD_TR = 0d0

      ALLOCATE( DIAG_STE_TR(IIPAR,JJPAR,N_TRACERS),STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DIAG_STE_TR' )
      DIAG_STE_TR = 0d0

      tstart = GET_TAU()

      IF ( DIAG_MFLUX ) THEN        ! Do no prod/loss diagnostic

         !%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Create the netCDF file
         !%%%%%%%%%%%%%%%%%%%%%%%%%
         filenameout = 'mflux.YYYY.MM.DD.nc'      
         call expand_date( filenameout, get_nymdb(), get_nhmsb() )
         CALL NcCr_Wr( fId, trim(filenameout) )
         
         ! Turn filling off
         CALL NcSetFill( fId, NF_NOFILL, omode )
         
         !%%%%%%%%%%%%%%%%%%%%%%%%
         ! Define the dimensions
         !%%%%%%%%%%%%%%%%%%%%%%%%
         
         CALL NcDef_Dimension( fId, 'lon', IIPAR, idLon )
         CALL NcDef_Dimension( fId, 'lat', JJPAR, idLat )
         CALL NcDef_Dimension( fId, 'lev', LLPAR, idLev )
         CALL NcDef_Dimension( fId, 'leve', LLPAR+1, idLevE )
         CALL NcDef_Dimension( fId, 'levi', LLPAR-1, idLevI ) ! Level interfaces
         CALL NcDef_Dimension( fId, 'tracer', N_TRACERS, idTr )
         CALL NcDef_Dimension( fId, 'time', 0, idtime ) ! 0 makes it unlimited

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Set the global attributes
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         CALL NcDef_Glob_Attributes( fId, 'title',
     &      'Archived horizonal and vertical QUS mass flux diagnostic' )

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Define variables and their attributes
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! Define longitude variable
         vId  = 0
         var1 = (/ idLon /)
         CALL NcDef_Variable( fId, 'lon', NF_DOUBLE, 1, var1, vId )
         CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'Longitude' )
         CALL NcDef_Var_Attributes( fId, vId, 'units', 'degrees_east' )
  
         ! Define latitude variable
         vId  = vId + 1
         var1 = (/ idLat /)
         CALL NcDef_Variable( fId, 'lat', NF_DOUBLE, 1, var1, vId )
         CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'Latitude' )
         CALL NcDef_Var_Attributes( fId, vId, 'units', 'degrees_north' )

         ! Define level variable
         vId  = vId + 1
         var1 = (/ idLev /)
         CALL NcDef_Variable( fId, 'lev', NF_INT, 1, var1, vId )
         CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'Level' )
         CALL NcDef_Var_Attributes( fId, vId, 
     &                              'units', 'level from surface' )

         ! Define level variable
         vId  = vId + 1
         var1 = (/ idLevI /)
         CALL NcDef_Variable( fId, 'levi', NF_INT, 1, var1, vId )
         CALL NcDef_Var_Attributes( fId, vId, 'long_name', 
     &                              'Level interface' )
         CALL NcDef_Var_Attributes( fId, vId, 
     &                       'units', 'level interface from surface' )

         ! Define east-west variable
         vId  = vId + 1
         var4 = (/ idLon, idLat, idLev, idtime /)
         CALL NcDef_Variable( fId,'pu', NF_DOUBLE, 4, var4, vId )
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'long_name', 'West-East QUS mass flux' )
         CALL NcDef_Var_Attributes( fId, vId, 
     &    'units', 'mb m^2 s^-1' )

         ! Define north-south variable
         vId  = vId + 1
         var4 = (/ idLon, idLat, idLev, idtime /)
         CALL NcDef_Variable( fId,'pv', NF_DOUBLE, 4, var4, vId )
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'long_name', 'South-North QUS mass flux' )
         CALL NcDef_Var_Attributes( fId, vId, 
     &    'units', 'mb m^2 s^-1' )

         ! Define up-down variable
         vId  = vId + 1
         var4 = (/ idLon, idLat, idLevI, idtime /)
         CALL NcDef_Variable( fId,'sd', NF_DOUBLE, 4, var4, vId )
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'long_name', 'Up-down QUS mass flux' )
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'units', 'mb m^2 s^-1' )

         ! Cross-tropopause mass flux
         vId  = vId + 1
         var3 = (/ idLon, idLat, idtime /)
         CALL NcDef_Variable( fId,'ste_flux', NF_DOUBLE, 3, var3, vId )
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'long_name', 'Up-down mass flux at the tropopause' )
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'units', 'mb m^2 s^-1' )

         ! Define east-west tracer flux variable
         vId  = vId + 1
         var5 = (/ idLon, idLat, idLev, idTr, idtime /)
         CALL NcDef_Variable( fId,'pu_tracer', NF_DOUBLE, 5, var5, vId )
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'long_name', 'West-East QUS tracer mass flux' )
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'units', 'kg s^-1' )
         
         ! Define north-south tracer flux variable
         vId  = vId + 1
         var5 = (/ idLon, idLat, idLev, idTr, idtime /)
         CALL NcDef_Variable( fId,'pv_tracer', NF_DOUBLE, 5, var5, vId )
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'long_name', 'South-North QUS tracer mass flux' )
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'units', 'kg s^-1' )

         ! Define up-down tracer flux variable
         vId  = vId + 1
         var5 = (/ idLon, idLat, idLevI, idTr, idtime /)
         CALL NcDef_Variable( fId,'sd_tracer', NF_DOUBLE, 5, var5, vId )
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'long_name', 'Up-down QUS tracer mass flux' )
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'units', 'kg s^-1' )

         ! Cross-tropopause tracer mass flux
         vId  = vId + 1
         var4 = (/ idLon, idLat, idTr, idtime /)
         CALL NcDef_Variable( fId,'ste_tracer_flux', 
     &        NF_DOUBLE, 4, var4, vId )
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'long_name', 'Up-down tracer mass flux at the tropopause')
         CALL NcDef_Var_Attributes( fId, vId, 
     &        'units', 'kg s^-1' )

         !%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! End definition section
         !%%%%%%%%%%%%%%%%%%%%%%%%%%
         CALL NcEnd_def( fID )

         !%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Populate dimension data
         !%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! Latitude
         DO J=1,JJPAR 
            YMID(J) = GET_YMID(1,J,1)
         ENDDO
         CALL NcWr( ymid, fId, 'lat', (/1/), (/JJPAR/) )

         ! Longitude
         DO I=1,IIPAR 
            XMID(I) = GET_XMID(I,1,1)
         ENDDO
         CALL NcWr( xmid, fId, 'lon', (/1/), (/IIPAR/) )

         ! Level and interface
         DO L=1,LLPAR
            Levels(L) = L
         ENDDO
         CALL NcWr( levels, fId, 'lev', (/1/), (/LLPAR/) )
         CALL NcWr( levels(1:(LLPAR-1)), fID, 'levi',(/1/),(/LLPAR-1/))

         CALL NcDo_Sync( fID )
         
         tindex = 1

      ENDIF
      
      END SUBROUTINE INIT_MASS_FLUX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_mass_flux
!
! !DESCRIPTION: Subroutine CLEANUP\_MASS\_FLUX deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_MASS_FLUX( Input_Opt )
!
! !USES:

      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE m_netcdf_io_close
! !REVISION HISTORY: 
!  1 Feb 2011 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      
      IF ( ALLOCATED( DIAG_PU ) ) DEALLOCATE( DIAG_PU )
      IF ( ALLOCATED( DIAG_PV ) ) DEALLOCATE( DIAG_PV )
      IF ( ALLOCATED( DIAG_SD ) ) DEALLOCATE( DIAG_SD )
      IF ( ALLOCATED( DIAG_STE ) ) DEALLOCATE( DIAG_STE )
      IF ( ALLOCATED( DIAG_PU_TR ) ) DEALLOCATE( DIAG_PU_TR )
      IF ( ALLOCATED( DIAG_PV_TR ) ) DEALLOCATE( DIAG_PV_TR )
      IF ( ALLOCATED( DIAG_SD_TR ) ) DEALLOCATE( DIAG_SD_TR )
      IF ( ALLOCATED( DIAG_STE_TR ) ) DEALLOCATE( DIAG_STE_TR )
      
      IF ( Input_Opt%LTRAN .and. DIAG_MFLUX ) CALL NcCl( fID )

      END SUBROUTINE CLEANUP_MASS_FLUX
!EOC
      END MODULE DIAG_MASS_FLUX_QUS
