!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: global_co2_mod.F90
!
! !DESCRIPTION: Module GLOBAL_CO2_MOD contains variables and routines
!  for reading global CO2.  This data is obtained from a NetCDF file.
!\\
!\\
! !INTERFACE:
!
MODULE Global_CO2_Mod
!
! !USES:
!
  USE Precision_Mod    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
#     include "netcdf.inc"
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Get_Global_CO2
  PUBLIC :: Cleanup_Global_CO2
!
  REAL(f4), ALLOCATABLE :: ALLCO2(:)
  INTEGER, ALLOCATABLE  :: ALLYEARS(:)
! !REMARKS:
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
! !IROUTINE: get_global_co2
!
! !DESCRIPTION: Retrieves the current global CO2 from a NetCDF file 
!  into the State\_Met derived type object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Global_CO2( am_I_root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE Error_Mod,          ONLY : Error_Stop
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE TIME_MOD,           ONLY : GET_YEAR

    USE m_netcdf_io_open                    ! netCDF open
    USE m_netcdf_io_get_dimlen              ! netCDF dimension queries
    USE m_netcdf_io_read                    ! netCDF data reads
    USE m_netcdf_io_close                   ! netCDF close
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  06 Jan 2015 - R. Yantosca - Initial version
!  26 May 2015 - S. D. Eastam - Modified for contrails
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    ! Scalars
    REAL(fp), PARAMETER :: STANDARD_CO2 = 3.90e-4_fp

    ! Logicals
    INTEGER, SAVE :: LAST_YEAR=-1
    INTEGER       :: IYEAR, fID

    ! Do this the old-fashioned way
    REAL(f4)      :: CURRCO2
    INTEGER, SAVE :: NYR_CO2
    INTEGER       :: CURR_YEAR
    INTEGER       :: st1d(1), ct1d(1)

    !=======================================================================
    ! GET_GLOBAL_CO2 begins here!
    !=======================================================================

    IF (.not.Input_Opt%LVARCO2) THEN
        State_Chm%GlobalCO2 = STANDARD_CO2
        RETURN
    ELSE

        IF (LAST_YEAR.lt.0) THEN
            ! Read in all data
            CALL NcOp_Rd (fId,TRIM(Input_Opt%CO2FILE))
            CALL Ncget_Dimlen( fID, 'time', NYR_CO2 )
            
            ! Start and count indices
            st1d = (/ 1       /)
            ct1d = (/ NYR_CO2 /)
            CALL NcRd( ALLCO2,   fId, 'ANNUAL_CO2', st1d, ct1d )
            CALL NcRd( ALLYEARS, fId, 'YEAR',       st1d, ct1d )
            CALL NcCl( fId )
        ENDIF

        CURR_YEAR = GET_YEAR()

        IF (CURR_YEAR.ne.LAST_YEAR) THEN

            IYEAR = 1
            DO WHILE (ALLYEARS(IYEAR).ne.CURR_YEAR)
                IYEAR = IYEAR + 1
            ENDDO

            CURRCO2 = ALLCO2(IYEAR)
    
            ! Add to State_Chm
            State_Chm%GlobalCO2 = CURRCO2

            LAST_YEAR = CURR_YEAR

        ENDIF
    ENDIF

  END SUBROUTINE Get_Global_CO2

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
  SUBROUTINE CLEANUP_GLOBAL_CO2 ( )
! 
! !REVISION HISTORY: 
!  29 Aug 2015 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

  !=================================================================
  ! CLEANUP_CO2 begins here!
  !=================================================================

  ! Cleanup the CO2 arrays 
  IF ( ALLOCATED( ALLYEARS ) ) DEALLOCATE( ALLYEARS )
  IF ( ALLOCATED( ALLCO2   ) ) DEALLOCATE( ALLCO2   )

  END SUBROUTINE Cleanup_Global_CO2
!EOC
END MODULE Global_CO2_Mod
