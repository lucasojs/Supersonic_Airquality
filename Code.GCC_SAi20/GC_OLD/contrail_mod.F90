!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: contrail_mod.F90
!
! !DESCRIPTION: Module CONTRAIL\_MOD contains variables and routines for 
!  reading and storing contrail data.  This data is  obtained from the 
!  HEMCO data structure.
!\\
!\\
! !INTERFACE:
!
MODULE Contrail_Mod
!
! !USES:
!
  USE Precision_Mod    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Get_Contrail_Grids
!
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
! !IROUTINE: get_contrail
!
! !DESCRIPTION: Copies the contrail data from the HEMCO data structure
!  into the State\_Met derived type object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Contrail_Grids( am_I_root, Input_Opt, State_Met, &
                                  iYear, iMonth, iDay, RC )
!
! !USES:
!
    USE Error_Mod,          ONLY : Error_Stop
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr 
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    INTEGER,        INTENT(IN)    :: iYear       ! Target year
    INTEGER,        INTENT(IN)    :: iMonth      ! Target month
    INTEGER,        INTENT(IN)    :: iDay        ! Target day
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
!  26 May 2015 - S. D. Eastam - Modified for contrails
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    ! Scalars
    LOGICAL :: FND

    ! Pointers
    REAL(f4), POINTER :: Ptr3D(:,:,:) => NULL()

    ! Time index
    INTEGER           :: tIdx
    LOGICAL           :: LCONTRAILS
    LOGICAL, SAVE     :: FIRST = .TRUE.

    !=======================================================================
    ! GET_CONTRAIL_GRIDS begins here!
    !=======================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Skip unless we are doing a fullchem simulation and using contrails
    LCONTRAILS = Input_Opt%ITS_A_FULLCHEM_SIM.and.Input_Opt%LCONTRAILS
    IF (.not.LCONTRAILS) THEN
       RETURN
    ENDIF

    !! Determine time index
    !tIdx = iYear + iMonth + iDay

    ! Nullify pointer
    Ptr3D => NULL()

    ! Get the pointer to the contrail data in the HEMCO data structure
    CALL HCO_GetPtr( am_I_Root, 'CON_OD_TOT', Ptr3D, RC, FOUND=FND )

    If ((RC == GIGC_SUCCESS) .and. (.not. FND)) Then
       ! Assume user wants these zeroed out
       IF (FIRST) THEN
          Write(*,*) ' ---   CONTRAIL DATA MISSING   --- '
          Write(*,*) ' --- SETTING CONTRAILS TO ZERO --- '
          FIRST=.FALSE.
       End If
       State_Met%CON_OD_TOT = 0.0d0
       State_Met%CON_SAD    = 0.0d0
       State_Met%CON_RAREA  = 1.0d-6
       State_Met%CON_RVOL   = 1.0d-6
       State_Met%CON_VOL    = 0.0d0
       State_Met%CON_NDENS  = 0.0d0
       State_Met%CON_FRC    = 0.0d0
       RETURN
    End If

    ! Stop with error message
    IF ( RC /= GIGC_SUCCESS .or. ( .not. FND ) ) THEN
       CALL ERROR_STOP ( 'Could not find CON_OD_TOT in HEMCO data list!', & 
                         'GET_CONTRAIL_GRIDS (contrail_mod.F90)' )
    ENDIF

    ! Add to State_Met
    State_Met%CON_OD_TOT = Ptr3D(:,:,:) 

    ! Get the pointer to the contrail data in the HEMCO data structure
    CALL HCO_GetPtr( am_I_Root, 'CON_SAD', Ptr3D, RC, FOUND=FND )

      ! Stop with error message
    IF ( RC /= GIGC_SUCCESS .or. ( .not. FND ) ) THEN
       CALL ERROR_STOP ( 'Could not find CON_SAD in HEMCO data list!', & 
                         'GET_CONTRAIL_GRIDS (contrail_mod.F90)' )
    ENDIF

    ! Add to State_Met
    State_Met%CON_SAD = Ptr3D(:,:,:) 

    ! Get the pointer to the contrail data in the HEMCO data structure
    CALL HCO_GetPtr( am_I_Root, 'CON_RAREA', Ptr3D, RC, FOUND=FND )

      ! Stop with error message
    IF ( RC /= GIGC_SUCCESS .or. ( .not. FND ) ) THEN
       CALL ERROR_STOP ( 'Could not find CON_RAREA in HEMCO data list!', & 
                         'GET_CONTRAIL_GRIDS (contrail_mod.F90)' )
    ENDIF

    ! Add to State_Met
    State_Met%CON_RAREA = Ptr3D(:,:,:) 

    ! Get the pointer to the contrail data in the HEMCO data structure
    CALL HCO_GetPtr( am_I_Root, 'CON_RVOL', Ptr3D, RC, FOUND=FND )

      ! Stop with error message
    IF ( RC /= GIGC_SUCCESS .or. ( .not. FND ) ) THEN
       CALL ERROR_STOP ( 'Could not find CON_RVOL in HEMCO data list!', & 
                         'GET_CONTRAIL_GRIDS (contrail_mod.F90)' )
    ENDIF

    ! Add to State_Met
    State_Met%CON_RVOL = Ptr3D(:,:,:) 

    ! Get the pointer to the contrail data in the HEMCO data structure
    CALL HCO_GetPtr( am_I_Root, 'CON_VOL', Ptr3D, RC, FOUND=FND )

      ! Stop with error message
    IF ( RC /= GIGC_SUCCESS .or. ( .not. FND ) ) THEN
       CALL ERROR_STOP ( 'Could not find CON_VOL in HEMCO data list!', & 
                         'GET_CONTRAIL_GRIDS (contrail_mod.F90)' )
    ENDIF

    ! Add to State_Met
    State_Met%CON_VOL = Ptr3D(:,:,:) 

    ! Get the pointer to the contrail data in the HEMCO data structure
    CALL HCO_GetPtr( am_I_Root, 'CON_NDENS', Ptr3D, RC, FOUND=FND )

      ! Stop with error message
    IF ( RC /= GIGC_SUCCESS .or. ( .not. FND ) ) THEN
       CALL ERROR_STOP ( 'Could not find CON_NDENS in HEMCO data list!', & 
                         'GET_CONTRAIL_GRIDS (contrail_mod.F90)' )
    ENDIF

    ! Add to State_Met
    State_Met%CON_NDENS = Ptr3D(:,:,:) 

    ! Get the pointer to the contrail data in the HEMCO data structure
    CALL HCO_GetPtr( am_I_Root, 'CON_FRC', Ptr3D, RC, FOUND=FND )

      ! Stop with error message
    IF ( RC /= GIGC_SUCCESS .or. ( .not. FND ) ) THEN
       CALL ERROR_STOP ( 'Could not find CON_FRC in HEMCO data list!', & 
                         'GET_CONTRAIL_GRIDS (contrail_mod.F90)' )
    ENDIF

    ! Add to State_Met
    State_Met%CON_FRC = Ptr3D(:,:,:) 

    ! Free the pointer
    Ptr3d => NULL()
    FIRST=.FALSE.

  END SUBROUTINE Get_Contrail_Grids
!EOC
END MODULE Contrail_Mod
