!@sum  RES_F23 Resolution info for 23 layer, 2x2.5 strat model
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
      INTEGER, PARAMETER :: IM=144,JM=90,LM=23, LS1=12

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
      REAL*8, PARAMETER :: PSF=984.d0, PTOP = 150.d0, PMTOP=2.0576514d-3
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
      REAL*8, PARAMETER :: PSFMPT=PSF-PTOP, PSTRAT=PTOP-PMTOP

!@var PLbot pressure levels at bottom of layers (mb)
      REAL*8, PARAMETER, DIMENSION(LM+1) :: PLbot = (/
     t     PSF,   960.d0,  929.d0,  884.d0,  819.d0,    ! Pbot L=1,5
     t  710.d0,   570.d0,  425.d0,  314.d0,  245.d0,    !      L=6,10
     t  192.d0,                                         !      L=11
     1    PTOP,                                         !      L=LS1
     s  117.d0,   86.2d0,  56.2d0,  31.6d0,  17.8d0,    !      L=13-17
     s  10.0d0,   4.63d0,  1.46d0,  0.460995258d0,      !      L=18-21
     s  0.1449994968d0,    3.12002802d-2,    PMTOP /)   !      L=..,LM+1

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
