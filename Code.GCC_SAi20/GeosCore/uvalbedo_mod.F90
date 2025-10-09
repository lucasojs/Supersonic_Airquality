!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: uvalbedo_mod.F90
!
! !DESCRIPTION: Module UVALBEDO\_MOD contains variables and routines for 
!  reading the UV Albedo data.  This data is required by the FAST-JX photolysis
!  module.  UV albedo data will now be obtained from the HEMCO data structure.
!\\
!\\
! !INTERFACE:
!
MODULE UValbedo_Mod
!
! !USES:
!
  USE Precision_Mod    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Get_UValbedo
  PUBLIC :: Cleanup_UValbedo
!
! !PRIVATE TYPES:
!
  ! Arrays                              
  REAL(f4),ALLOCATABLE :: AVG_ALBEDO(:,:)
!
! !REMARKS:
!  References:
!  ============================================================================
!  Herman, J.R and Celarier, E.A., "Earth surface reflectivity climatology
!    at 340-380 nm from TOMS data", __J. Geophys. Res__, Vol. 102, D23, 
!    pp. 28003-28011, Dec 20, 1997.
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!   19 Apr 2002 - R. Yantosca - Initial version
!  (1 ) Now read uvalbedo file directly from DATA_DIR/uvalbedo_200111
!        subdirectory.  (bmy, 4/2/02)
!  (2 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections. (bmy, 5/28/02)
!  (3 ) Now references "error_mod.f" (bmy, 10/15/02)
!  (4 ) Minor modification in READ_UVALBEDO (bmy, 3/14/03)
!  (5 ) Now references "directory_mod.f" (bmy, 7/20/04)
!  (6 ) Bug fix for GCAP grid in READ_UVALBEDO (bmy, 8/16/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  24 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  12 Jan 2015 - R. Yantosca - Remove CLEANUP_UVALBEDO routine
!  04 Mar 2015 - R. Yantosca - UV albedo data now comes via HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_uvalbedo
!
! !DESCRIPTION: Copies the UV Albedo data from the HEMCO data structure
!  into the State\_Met derived type object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_UValbedo( am_I_root, Input_Opt, State_Met, RC )
!
! !USES:
!
    USE Error_Mod,          ONLY : Error_Stop
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr 
#if defined( GISS )
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR
    USE TIME_MOD,           ONLY : GET_MONTH
    USE ERROR_MOD,          ONLY : ALLOC_ERR

    USE m_netcdf_io_open
    USE m_netcdf_io_close
    USE m_netcdf_io_read
#endif

!
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    ! Scalars
    LOGICAL :: FND

    ! Integers
    INTEGER :: MONTH, fID

#if defined( GISS )
    ! Don't bother with pointer
    REAL(f4)            :: Q2D(IIPAR,JJPAR)
    INTEGER             :: ST3D(3), CT3D(3), AS, I, J
    LOGICAL,SAVE        :: FIRST=.TRUE.
    INTEGER,ALLOCATABLE :: AVG_COUNT(:,:)
#else
    ! Pointers
    REAL(f4), POINTER   :: Ptr2D(:,:) => NULL()
#endif

    !=======================================================================
    ! READ_UVALBEDO begins here!
    !=======================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Skip unless we are doing a fullchem or aerosol-only simulation
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) .and. &
         ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
       RETURN
    ENDIF

#if defined( GISS )
    ! Assume that UV albedo is similar to VIS albedo (!)
    CALL NcOp_Rd( fID, TRIM( Input_Opt%DATA_DIR ) // 'BGRID/' // &
         TRIM( Input_Opt%RUNID ) // '.acc.ymonmean.nc' )
    CT3D = (/IIPAR,JJPAR,1/)

    ! If it's the first month, get the annual average albedo
    ! Some areas have NaN albedo values in some months, so we need to
    ! store the annual average value for those times
    IF (FIRST) THEN
       ! Get average of non-NaN values
       ! Ignore differences in month lengths
       ALLOCATE( AVG_COUNT( IIPAR, JJPAR ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'AVG_COUNT' )
       AVG_COUNT = 0
       ALLOCATE( AVG_ALBEDO( IIPAR, JJPAR ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'AVG_ALBEDO' )
       AVG_ALBEDO = 0.0e+0_fp
       DO MONTH=1,12
          ST3D = (/1,1,MONTH/)
          CALL NcRd( Q2D, fID, 'grnd_alb', ST3D, CT3D )
          ! Test for NaNs
          WHERE (Q2D.gt.0.0e+0_f4)
             AVG_ALBEDO = AVG_ALBEDO + Q2D
             AVG_COUNT = AVG_COUNT + 1
          ENDWHERE
       ENDDO
       AVG_ALBEDO = AVG_ALBEDO / (1.0e+0_f4 * AVG_COUNT)
       DEALLOCATE( AVG_COUNT )
       FIRST=.FALSE.
    ENDIF

    MONTH = GET_MONTH()
    ST3D = (/1,1,MONTH/)
    CALL NcRd( Q2D, fID, 'grnd_alb', ST3D, CT3D )
    CALL NcCl( fID )

    ! Some entries are NaN
    WHERE( Q2D .lt. 0.0e+0_f4 ) Q2D = AVG_ALBEDO

    ! Convert from % to fraction
    Q2D = Q2D / 100.0e+0_f4

    ! Add to State_Met
    State_Met%UVALBEDO = Q2D(:,:) 
#else
    ! Nullify pointer
    Ptr2D => NULL()

    ! Get the pointer to the UV albedo data in the HEMCO data structure
    CALL HCO_GetPtr( am_I_Root, 'UV_ALBEDO', Ptr2D, RC, FOUND=FND )

      ! Stop with error message
    IF ( RC /= GIGC_SUCCESS .or. ( .not. FND ) ) THEN
       CALL ERROR_STOP ( 'Could not find UV_ALBEDO in HEMCO data list!', & 
                         'READ_UVALBEDO (uvalbedo_mod.F90)' )
    ENDIF

    ! Add to State_Met
    State_Met%UVALBEDO = Ptr2D(:,:) 

    ! Free the pointer
    Ptr2d => NULL()
#endif

  END SUBROUTINE Get_UValbedo
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_uvalbedo
!
! !DESCRIPTION: Subroutine CLEANUP\_UVALBEDO deallocates all allocated
!  global module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_UValbedo( am_I_Root )
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN) :: am_I_Root   ! Are we on the root CPU?
!
! !REVISION HISTORY:'
!  18 Aug 2015 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( ALLOCATED( AVG_ALBEDO ) ) DEALLOCATE( AVG_ALBEDO )

  END SUBROUTINE Cleanup_UValbedo
!EOC
END MODULE UValbedo_Mod
