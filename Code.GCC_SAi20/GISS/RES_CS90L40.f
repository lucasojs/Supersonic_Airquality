!@sum  RES_CS90L40 Resolution info for 40 layer, CS90 cubed sphere model
!@ver  1.0

      MODULE RESOLUTION
!@sum  RESOLUTION contains horiz/vert resolution variables
!@ver  1.0
      IMPLICIT NONE
      SAVE
!@var IM,JM number of grid boxes on eahc cube face
!@var LM number of vertical levels
!@var LS1 Layers LS1->LM: constant pressure levels, L<LS1: sigma levels
      INTEGER, PARAMETER :: IM=90,JM=90,LM=40, LS1=24

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
      REAL*8, PARAMETER :: PSF=984.d0, PTOP = 150.d0, PMTOP = .1d0
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
      REAL*8, PARAMETER :: PSFMPT=PSF-PTOP, PSTRAT=PTOP-PMTOP

!@var PLbot pressure levels at bottom of layers (mb)
      REAL*8, PARAMETER, DIMENSION(LM+1) :: PLbot = (/
     *     PSF,   964d0, 942d0, 917d0, 890d0, 860d0, 825d0, ! Pbot L=1,..
     *     785d0, 740d0, 692d0, 642d0, 591d0, 539d0, 489d0, !      L=...
     *     441d0, 396d0, 354d0, 316d0, 282d0, 251d0, 223d0,
     *     197d0, 173d0,
     *     PTOP,                                            !      L=LS1
     *     128d0, 108d0,  90d0,  73d0,  57d0,  43d0,  31d0, !      L=...
     *     20d0,   10d0,  5.62d0,  3.16d0,  1.78d0,   1.d0,
     *     .562d0,  .316d0,  .178d0,  PMTOP /)              !      L=..,LM+1

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
