!@sum  RES_M20S Resolution info for 20 layer, 4x5 strat model
!@+    (M18b with 2 more layers between 934mb and 150mb)
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
      INTEGER, PARAMETER :: IM=72,JM=46,LM=20, LS1=13

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
      REAL*8, PARAMETER :: PSF=984.d0, PTOP = 150.d0, PMTOP = .3d0
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
      REAL*8, PARAMETER :: PSFMPT=PSF-PTOP, PSTRAT=PTOP-PMTOP

!@var PLbot pressure levels at bottom of layers (mb)
      REAL*8, PARAMETER, DIMENSION(LM+1) :: PLbot = (/
     t     PSF, 964.d0, 934.d0, 890.d0, 830.d0, 750.d0, ! Pbot L=1,..
     t  650.d0, 545.d0, 445.d0, 355.d0, 275.d0, 205.d0, !      L=...
     1    PTOP,                                         !      L=LS1
     s  110.d0,  80.d0,   55.d0,   35.d0,    20.d0,     !      L=...
     s  10.d0,    3.d0,   PMTOP /)                      !      L=..,LM+1

C**** KEP depends on whether stratos. EP flux diagnostics are calculated
C**** If dummy EPFLUX is used set KEP=0, otherwise KEP=21
!@param KEP number of lat/height E-P flux diagnostics
      INTEGER, PARAMETER :: KEP=21

C**** Based on model top, determine how much of stratosphere is resolved
C****         PMTOP >= 10 mb,    ISTRAT = 0
C**** 1 mb <= PMTOP <  10 mb,    ISTRAT = 1
C****         PMTOP <   1 mb,    ISTRAT = 2
      INTEGER, PARAMETER :: ISTRAT = 2

      END MODULE RESOLUTION

C**** The vertical resolution also determines whether
C**** stratospheric wave drag will be applied or not.
C**** This resolution is for a stratospheric model and so must be used
C**** in conjunction with the strat. modules
