!@sum  RES_F53 Resolution info for fine 53 layer, 2x2.5 strat model
!@+    (top at .002 mb, GWDRAG)
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
      INTEGER, PARAMETER :: IM=144,JM=90,LM=53, LS1=27

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
      REAL*8, PARAMETER :: PSF=984.d0, PTOP = 154.d0, PMTOP = 2.0576514d-3
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
      REAL*8, PARAMETER :: PSFMPT=PSF-PTOP, PSTRAT=PTOP-PMTOP

!@var PLbot pressure levels at bottom of layers (mb)
      REAL*8, PARAMETER, DIMENSION(LM+1) :: PLbot = (/
     t     PSF,   954d0,  924d0,  894d0,  864d0,  834d0,  804d0,  774d0,
     t     744d0, 714d0,  680d0,  636d0,  591d0,  545d0,  500d0,  457d0,
     t     417d0, 381d0,  348d0,  318d0,  290d0,  264d0,  240d0,  217d0,
     t     195d0, 174d0,                               ! Pbot L=1,26
     1     PTOP,                                       !      L=LS1
     s     135d0, 117d0,  100d0,   84d0,   69d0,   55d0,   42d0,   30d0,
     s      19d0,  10d0, 7.1043d0, 5.04711d0, 3.58563d0, 2.54734d0,
     s     1.80971d0, 1.28567d0, 0.913382d0, 0.561994d0, 0.316003d0,
     s     0.177995d0, 0.1d0, 0.0562045d0, 0.0315979d0, 0.0178033d0, 
     s     0.01d0, 0.00518932d0, PMTOP /)              !      L=..,LM+1

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

