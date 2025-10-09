!@sum  RES_M24T Resolution info for 24 layer, 4x5 non-strat model
!@+    (M20a with 4 more layers between 984mb and 150mb - NO GWdrag)
!@auth Original Development Team
!@ver  1.0

      MODULE RESOLUTION
!@sum  RESOLUTION contains horiz/vert resolution variables
!@auth Original Development Team
!@ver  1.0
      IMPLICIT NONE
      SAVE
!@var IM,JM longitudinal and latitudinal number of grid boxes
!@var LM number of vertical levels
!@var LS1 Layers LS1->LM: constant pressure levels, L<LS1: sigma levels
      INTEGER, PARAMETER :: IM=72,JM=46,LM=24, LS1=17

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
      REAL*8, PARAMETER :: PSF=984.d0, PTOP = 150.d0, PMTOP = .3d0
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
      REAL*8, PARAMETER :: PSFMPT=PSF-PTOP, PSTRAT=PTOP-PMTOP

!@var PLbot pressure levels at bottom of layers (mb)
      REAL*8, PARAMETER, DIMENSION(LM+1) :: PLbot = (/
     t     PSF, 974.d0, 959.d0, 934.d0, 900.d0, 857.d0, ! Pbot L=1,..
     t  804.d0, 744.d0, 677.d0, 607.d0, 535.d0, 463.d0, !      L=...
     t  393.d0, 325.d0, 260.d0, 200.d0,                 !      L=...
     1    PTOP,                                         !      L=LS1
     s  110.d0,  80.d0,   55.d0,   35.d0,    20.d0,     !      L=...
     s  10.d0,    3.d0,   PMTOP /)                      !      L=..,LM+1

C**** KEP depends on whether stratos. EP flux diagnostics are calculated
C**** If dummy EPFLUX is used set KEP=0, otherwise KEP=21
!@param KEP number of lat/height E-P flux diagnostics
      INTEGER, PARAMETER :: KEP=0

C**** Based on model top, determine how much of stratosphere is resolved
C****         PMTOP >= 10 mb,    ISTRAT = 0
C**** 1 mb <= PMTOP <  10 mb,    ISTRAT = 1
C****         PMTOP <   1 mb,    ISTRAT = 2
      INTEGER, PARAMETER :: ISTRAT = 2

      END MODULE RESOLUTION

C**** The vertical resolution also determines whether
C**** stratospheric wave drag will be applied or not.
C**** Hence also included here are some dummy routines for non-strat
C**** models.

      SUBROUTINE DUMMY_STRAT
!@sum DUMMY dummy routines for non-stratospheric models
C**** Dummy routines in place of STRATDYN
      ENTRY init_GWDRAG
      ENTRY GWDRAG
      ENTRY VDIFF
      ENTRY io_strat
      ENTRY alloc_strat_com
C**** Dummy routines in place of STRAT_DIAG (EP flux calculations)
C**** Note that KEP=0 is set to zero above for the dummy versions.
      ENTRY EPFLUX
      ENTRY EPFLXI
      ENTRY EPFLXP
C****
      RETURN
      END SUBROUTINE DUMMY_STRAT
