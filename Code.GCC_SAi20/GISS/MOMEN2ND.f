      module MOMENTS
      implicit none
      private

      public ADVECV, moment_enq_order

      contains
!      subroutine init_MOM
!!@sum  init_MOM sets an order dependent coefficient for AVRX
!      USE DYNAMICS, only : XAVRX
!      XAVRX = 1. ! for second order scheme, byrt2 for 4th order scheme
!      CALL AVRX0
!
!      RETURN
!      END subroutine init_MOM

      subroutine moment_enq_order(order)
!@sum moment_enq_order returns order of the scheme
      implicit none
      integer, intent(out) :: order
      order = 2
      return
      end subroutine moment_enq_order

      SUBROUTINE ADVECV (PA,UT,VT,PB,U,V,P,DT1)
!@sum  ADVECV Advects momentum (incl. coriolis) using mass fluxes
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,imh,jm,lm,ls1,mrch,dsig,psfmpt,modd5k
     &     ,do_polefix
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, GRID,NORTH,SOUTH,GET
      USE DOMAIN_DECOMP_1D, only : haveLatitude
      USE GEOM, only : fcor,dxyv,dxyn,dxys,dxv,ravpn,ravps
     &     ,sini=>siniv,cosi=>cosiv,acor,polwt
      USE DYNAMICS, only : pu,pv,pit,sd,spa,dut,dvt,conv
      USE DYNAMICS, only : t_advecv
c      USE DIAG, only : diagcd
      IMPLICIT NONE
      REAL*8 U(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     * V(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     * P(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)
      REAL*8 UT(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     * VT(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     * PA(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     * PB(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     * FD(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      REAL*8, SAVE :: SMASS(JM)

      INTEGER,SAVE :: IFIRST = 1
      INTEGER I,J,IP1,IM1,L  !@var I,J,IP1,IM1,L  loop variables
      REAL*8 VMASS,RVMASS,ALPH,PDT4,DT1,DT2,DT4,DT8,DT12,DT24
     *     ,FLUX,FLUXU,FLUXV
      REAL*8 FLUXU_N_S,FLUXV_N_S
      REAL*8 FLUXU_SW_NE,FLUXV_SW_NE
      REAL*8 FLUXU_SE_NW,FLUXV_SE_NW

      REAL*8   VMASS2(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     *         ASDU(IM,GRID%J_STRT_SKP:GRID%J_STOP_HALO,LM-1)
c pole fix variables
      integer :: hemi,jpo,jns,jv,jvs,jvn,jj,ipole
      real*8 :: utmp,vtmp,wts
      real*8, dimension(im) :: dmt
      real*8, dimension(im,2) :: usv,vsv,usv0,vsv0
C****
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO=J_0H,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      !IF(MODD5K.LT.MRCH) CALL DIAG5F (U,V)
      IF(IFIRST.EQ.1) THEN
        IFIRST=0
        DO 10 J=J_0S,J_1
 10     SMASS(J)=PSFMPT*DXYV(J)
      END IF

      DT2=DT1/2.
      DT4=DT1/4.
      DT8=DT1/8.
      DT12=DT1/12.
      DT24=DT1/24.
C****
C**** SCALE UT AND VT WHICH MAY THEN BE PHYSICALLY INTERPRETED AS
C**** MOMENTUM COMPONENTS
C****
C     I=IM
C     DO 120 J=2,JM
C     DO 120 IP1=1,IM
C     VMASS=.5*((PA(I,J-1)+PA(IP1,J-1))*DXYN(J-1)
C    *  +(PA(I,J)+PA(IP1,J))*DXYS(J))
C     DO 110 L=1,LS1-1
C     UT(I,J,L)=UT(I,J,L)*VMASS*DSIG(L)
C 110 VT(I,J,L)=VT(I,J,L)*VMASS*DSIG(L)
C 120 I=IP1
C     DO 150 L=LS1,LM
C     DO 150 J=2,JM
C     VMASS=SMASS(J)*DSIG(L)
C     DO 150 I=1,IM
C     UT(I,J,L)=UT(I,J,L)*VMASS
C 150 VT(I,J,L)=VT(I,J,L)*VMASS
C     DUT=0.
C     DVT=0.
C
      CALL HALO_UPDATE(grid, PA, FROM=SOUTH)
      DO 110 J=J_0S,J_1
      I=IM
      DO 110 IP1=1,IM
         VMASS2(I,J)=.5*((PA(I,J-1)+PA(IP1,J-1))*DXYN(J-1)
     *                 +(PA(I,J)+PA(IP1,J))*DXYS(J))
         I=IP1
  110 CONTINUE
!$OMP  PARALLEL DO PRIVATE(J,L,VMASS)
      DO L=1,LM
        IF(L.LT.LS1) THEN  !  DO L=1,LS1-1
          DO J=J_0S,J_1
            DUT(:,J,L)=0.0
            DVT(:,J,L)=0.0
            UT(:,J,L)=UT(:,J,L)*VMASS2(:,J)*DSIG(L)
            VT(:,J,L)=VT(:,J,L)*VMASS2(:,J)*DSIG(L)
          END DO
        ELSE               !  DO L=LS1,LM
          DO J=J_0S,J_1
            VMASS=SMASS(J)*DSIG(L)
            DUT(:,J,L)=0.0
            DVT(:,J,L)=0.0
            UT(:,J,L)=UT(:,J,L)*VMASS
            VT(:,J,L)=VT(:,J,L)*VMASS
          END DO
        END IF
      END DO
!$OMP  END PARALLEL DO
C****
C**** BEGINNING OF LAYER LOOP
C****
      CALL HALO_UPDATE(GRID,U  ,FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(GRID,V  ,FROM=SOUTH+NORTH)
CAOO no need to communicate, local compute      CALL HALO_UPDATE(GRID,DUT,FROM)
CAOO no need to communicate, local compute      CALL HALO_UPDATE(GRID,DVT,FROM)

!$OMP  PARALLEL DO PRIVATE(I,IP1,J,L,FLUX,FLUXU,FLUXV,
!$OMP+   FLUXU_N_S,FLUXV_N_S,FLUXU_SW_NE, FLUXV_SW_NE,
!$OMP+   FLUXU_SE_NW, FLUXV_SE_NW, IPOLE,JV,JVS,JVN,WTS,USV0,VSV0) 
      DO 300 L=1,LM

c
c interpolate polar velocities to the appropriate latitude
c
      do ipole=1,2
      if((haveLatitude(grid,J=2) .or. haveLatitude(grid,J=3)) .and. 
     &        ipole.eq.1) then
         jv = 2 ! why not staggered grid
         jvs = 2 ! jvs is the southernmost velocity row
         jvn = jvs + 1 ! jvs is the northernmost velocity row
         wts = polwt
      else if(haveLatitude(grid,J=JM) .and. ipole.eq.2) then
         jv = JM ! why not staggered grid
         jvs = jv - 1
         jvn = jvs + 1
         wts = 1.-polwt
      else
         cycle
      endif
      usv0(:,ipole) = u(:,jv,l)
      vsv0(:,ipole) = v(:,jv,l)
      u(:,jv,l) = wts*u(:,jvs,l) + (1.-wts)*u(:,jvn,l)
      v(:,jv,l) = wts*v(:,jvs,l) + (1.-wts)*v(:,jvn,l)
      enddo
C****
C**** HORIZONTAL ADVECTION OF MOMENTUM
C****
      I=IM
      DO 230 IP1=1,IM
C**** CONTRIBUTION FROM THE WEST-EAST MASS FLUX
      DO 210 J=J_0STG,J_1STG
      FLUX=DT12*(PU(IP1,J,L)+PU(IP1,J-1,L)+PU(I,J,L)+PU(I,J-1,L))
      FLUXU=FLUX*(U(IP1,J,L)+U(I,J,L))
      DUT(IP1,J,L)=DUT(IP1,J,L)+FLUXU
      DUT(I,J,L)  =DUT(I,J,L)  -FLUXU
      FLUXV=FLUX*(V(IP1,J,L)+V(I,J,L))
      DVT(IP1,J,L)=DVT(IP1,J,L)+FLUXV
  210 DVT(I,J,L)  =DVT(I,J,L)  -FLUXV

      IF (haveLatitude(grid, J=1)) THEN ! No southern contribution off boundary.
        FLUXU_N_S  =0.
        FLUXV_N_S  =0.
        FLUXU_SW_NE=0.
        FLUXV_SW_NE=0.
        FLUXU_SE_NW=0.
        FLUXV_SE_NW=0.
      ELSE ! from lower boundary
        J=J_0-1
        FLUX=DT12*(PV(I,J,L)+PV(IP1,J,L)+PV(I,J+1,L)+PV(IP1,J+1,L))
        FLUXU_N_S=FLUX*(U(I,J,L)+U(I,J+1,L))
        FLUXV_N_S=FLUX*(V(I,J,L)+V(I,J+1,L))
        FLUX=DT24*(PU(IP1,J,L)+PU(I,J,L)+PV(IP1,J,L)+PV(IP1,J+1,L))
        FLUXU_SW_NE=FLUX*(U(IP1,J+1,L)+U(I,J,L))
        FLUXV_SW_NE=FLUX*(V(IP1,J+1,L)+V(I,J,L))
        FLUX=DT24*(-PU(IP1,J,L)-PU(I,J,L)+PV(IP1,J,L)+PV(IP1,J+1,L))
        FLUXU_SE_NW=FLUX*(U(I,J+1,L)+U(IP1,J,L))
        FLUXV_SE_NW=FLUX*(V(I,J+1,L)+V(IP1,J,L))
      END IF

      DO J=J_0S,J_1S

        DUT(I,J,L)   = DUT(I,J,L)   + FLUXU_N_S
        DVT(I,J,L)   = DVT(I,J,L)   + FLUXV_N_S
        DUT(IP1,J,L) = DUT(IP1,J,L) + FLUXU_SW_NE
        DVT(IP1,J,L) = DVT(IP1,J,L) + FLUXV_SW_NE
        DUT(I,J,L)   = DUT(I,J,L)   + FLUXU_SE_NW
        DVT(I,J,L)   = DVT(I,J,L)   + FLUXV_SE_NW

C**** CONTRIBUTION FROM THE SOUTH-NORTH MASS FLUX

        FLUX=DT12*(PV(I,J,L)+PV(IP1,J,L)+PV(I,J+1,L)+PV(IP1,J+1,L))
        FLUXU_N_S=FLUX*(U(I,J,L)+U(I,J+1,L))
        FLUXV_N_S=FLUX*(V(I,J,L)+V(I,J+1,L))
      
        DUT(I,J,L)=DUT(I,J,L)-FLUXU_N_S
        DVT(I,J,L)=DVT(I,J,L)-FLUXV_N_S

C**** CONTRIBUTION FROM THE SOUTHWEST-NORTHEAST MASS FLUX
        FLUX=DT24*(PU(IP1,J,L)+PU(I,J,L)+PV(IP1,J,L)+PV(IP1,J+1,L))
        FLUXU_SW_NE=FLUX*(U(IP1,J+1,L)+U(I,J,L))
        FLUXV_SW_NE=FLUX*(V(IP1,J+1,L)+V(I,J,L))

        DUT(I,J,L) = DUT(I,J,L) - FLUXU_SW_NE
        DVT(I,J,L) = DVT(I,J,L) - FLUXV_SW_NE

C**** CONTRIBUTION FROM THE SOUTHEAST-NORTHWEST MASS FLUX
        FLUX=DT24*(-PU(IP1,J,L)-PU(I,J,L)+PV(IP1,J,L)+PV(IP1,J+1,L))
        FLUXU_SE_NW=FLUX*(U(I,J+1,L)+U(IP1,J,L))
        FLUXV_SE_NW=FLUX*(V(I,J+1,L)+V(IP1,J,L))

        DUT(IP1,J,L) = DUT(IP1,J,L) - FLUXU_SE_NW
        DVT(IP1,J,L) = DVT(IP1,J,L) - FLUXV_SE_NW

      END DO

      IF (HAVE_NORTH_POLE) THEN
        J=JM
        DUT(I,J,L) = DUT(I,J,L) + FLUXU_N_S
        DVT(I,J,L) = DVT(I,J,L) + FLUXV_N_S

        DUT(IP1,J,L) = DUT(IP1,J,L) + FLUXU_SW_NE
        DVT(IP1,J,L) = DVT(IP1,J,L) + FLUXV_SW_NE

        DUT(I,J,L) = DUT(I,J,L) + FLUXU_SE_NW
        DVT(I,J,L) = DVT(I,J,L) + FLUXV_SE_NW
      END IF

  230 I=IP1

c restore uninterpolated values of u,v at the pole
      do ipole=1,2
      if((haveLatitude(grid,J=2) .or. haveLatitude(grid,J=3)) .and. 
     &        ipole.eq.1) then
         jv = 2 ! why not staggered grid
      else if(haveLatitude(grid,J=JM) .and. ipole.eq.2) then
         jv = JM ! why not staggered grid
      else
         cycle
      endif
      u(:,jv,l) = usv0(:,ipole)
      v(:,jv,l) = vsv0(:,ipole)
      enddo

  300 CONTINUE
!$OMP  END PARALLEL DO

      if(do_polefix.eq.1) then
c Horizontal advection for the polar row is performed upon
c x-y momentum rather than spherical-coordinate momentum,
c eliminating the need for the metric term (and for the latter
c to be exactly consistent with the advection scheme).
c Southwest-Northeast and Southeast-Northwest corner fluxes
c are discarded since they produce erroneous tendencies at the pole
c in models with a polar half-box.
c Explicit cross-polar advection will be ignored until issues with
c corner fluxes and spherical geometry can be resolved.
       
        do ipole=1,2
          if(have_south_pole .and. ipole.eq.1) then
            hemi = -1
            jpo = 1
            jns = jpo + 1
            jv = 2           ! why not staggered grid
            jvs = 2          ! jvs is the southernmost velocity row
            jvn = jvs + 1       ! jvs is the northernmost velocity row
            wts = polwt
          else if(have_north_pole .and. ipole.eq.2) then
            hemi = +1
            jpo = JM
            jns = jpo - 1
            jv = JM            ! why not staggered grid
            jvs = jv - 1
            jvn = jvs + 1
            wts = 1.-polwt
          else
            cycle
          endif
c loop over layers
      do l=1,lm
c
c Copy u,v into temporary storage and transform u,v to x-y coordinates
c
      do j=jvs,jvn
         jj = j-jvs+1
         do i=1,im
            usv(i,jj) = u(i,j,l)
            vsv(i,jj) = v(i,j,l)
            u(i,j,l) = cosi(i)*usv(i,jj)-hemi*sini(i)*vsv(i,jj)
            v(i,j,l) = cosi(i)*vsv(i,jj)+hemi*sini(i)*usv(i,jj)
         enddo
      enddo
c
c interpolate polar velocities to the appropriate latitude
c
      u(:,jv,l) = wts*u(:,jvs,l) + (1.-wts)*u(:,jvn,l)
      v(:,jv,l) = wts*v(:,jvs,l) + (1.-wts)*v(:,jvn,l)

c
c Compute advective tendencies of xy momentum in the polar rows
c
      dmt(:) = 0.
      dut(:,jv,l) = 0.
      dvt(:,jv,l) = 0.
c
      i = im
      do ip1=1,im
C**** CONTRIBUTION FROM THE WEST-EAST MASS FLUX
      FLUX=DT8*(PU(IP1,JPO,L)+PU(I,JPO,L)+PU(IP1,JNS,L)+PU(I,JNS,L))
      FLUXU=FLUX*(U(IP1,JV,L)+U(I,JV,L))
      DUT(IP1,JV,L)=DUT(IP1,JV,L)+FLUXU
      DUT(I,JV,L)  =DUT(I,JV,L)  -FLUXU
      FLUXV=FLUX*(V(IP1,JV,L)+V(I,JV,L))
      DVT(IP1,JV,L)=DVT(IP1,JV,L)+FLUXV
      DVT(I,JV,L)  =DVT(I,JV,L)  -FLUXV
      DMT(IP1)=DMT(IP1)+FLUX+FLUX
      DMT(I)  =DMT(I)  -FLUX-FLUX
C**** CONTRIBUTION FROM THE SOUTH-NORTH MASS FLUX
      FLUX=DT8*(PV(I,JVS,L)+PV(IP1,JVS,L)+PV(I,JVN,L)+PV(IP1,JVN,L))
      FLUX=FLUX*HEMI
      DUT(I,JV,L)=DUT(I,JV,L)+FLUX*(U(I,JVS,L)+U(I,JVN,L))
      DVT(I,JV,L)=DVT(I,JV,L)+FLUX*(V(I,JVS,L)+V(I,JVN,L))
      DMT(I)=DMT(I)+FLUX+FLUX
      i = ip1
      enddo ! i
c
c correct for the too-large dxyv in the polar row
c and convert dut,dvt from xy to polar coordinates
c
      do i=1,im
         dut(i,jv,l) = dut(i,jv,l) +
     &        (acor-1d0)*(dut(i,jv,l)-dmt(i)*u(i,jv,l))
         dvt(i,jv,l) = dvt(i,jv,l) +
     &        (acor-1d0)*(dvt(i,jv,l)-dmt(i)*v(i,jv,l))
         utmp = dut(i,jv,l)
         vtmp = dvt(i,jv,l)
         dut(i,jv,l) = cosi(i)*utmp+hemi*sini(i)*vtmp
         dvt(i,jv,l) = cosi(i)*vtmp-hemi*sini(i)*utmp
      enddo
c
c get the untransformed u,v back from storage space
c
      do j=jvs,jvn
         jj = j-jvs+1
         do i=1,im
            u(i,j,l) = usv(i,jj)
            v(i,j,l) = vsv(i,jj)
         enddo
      enddo
      enddo ! loop over layers
      enddo ! loop over poles
      endif

C****
C**** VERTICAL ADVECTION OF MOMENTUM
C****
C     DO 310 L=1,LM-1
C     DO 310 J=2,JM
C     I=IM
C     DO 310 IP1=1,IM
C     SDU=DT2*((SD(I,JJ(J-1),L)+SD(IP1,JJ(J-1),L))*RAVPN(J-1)+
C    *  (SD(I,JJ(J),L)+SD(IP1,JJ(J),L))*RAVPS(J))
C     DUT(I,J,L)  =DUT(I,J,L)  +SDU*(U(I,J,L)+U(I,J,L+1))
C     DUT(I,J,L+1)=DUT(I,J,L+1)-SDU*(U(I,J,L)+U(I,J,L+1))
C     DVT(I,J,L)  =DVT(I,J,L)  +SDU*(V(I,J,L)+V(I,J,L+1))
C     DVT(I,J,L+1)=DVT(I,J,L+1)-SDU*(V(I,J,L)+V(I,J,L+1))
C 310 I=IP1
!!! MUST USE CONV for HALO update not SD which is aliased to it.
      CALL HALO_UPDATE(GRID,SD,FROM=SOUTH)
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO L=1,LM-1
      DO J=J_0S,J_1
         DO I=1,IM-1
           ASDU(I,J,L)=DT2*(
     *                (SD(I,J-1,L)+SD(I+1,J-1,L))
     *                 *RAVPN(J-1)+
     *                (SD(I,J,L)+SD(I+1,J,L))
     *                 *RAVPS(J) )
         END DO
           ASDU(IM,J,L)=DT2*(
     *                (SD(IM,J-1,L)+SD(1,J-1,L))
     *                  *RAVPN(J-1)+
     *                (SD(IM,J,L)+SD(1,J,L))
     *                  *RAVPS(J) )
      END DO
      END DO
!$OMP  END PARALLEL DO

      L=1
      DO J=J_0S,J_1
         DUT(:,J,L)  =DUT(:,J,L)  +ASDU(:,J,L)  *(U(:,J,L)+U(:,J,L+1))
         DVT(:,J,L)  =DVT(:,J,L)  +ASDU(:,J,L)  *(V(:,J,L)+V(:,J,L+1))
      END DO
!$OMP  PARALLEL DO PRIVATE(J,L)
      DO L=2,LM-1
        DO J=J_0S,J_1
         DUT(:,J,L)  =DUT(:,J,L)  -ASDU(:,J,L-1)*(U(:,J,L-1)+U(:,J,L))
         DUT(:,J,L)  =DUT(:,J,L)  +ASDU(:,J,L)  *(U(:,J,L)+U(:,J,L+1))
         DVT(:,J,L)  =DVT(:,J,L)  -ASDU(:,J,L-1)*(V(:,J,L-1)+V(:,J,L))
         DVT(:,J,L)  =DVT(:,J,L)  +ASDU(:,J,L)  *(V(:,J,L)+V(:,J,L+1))
        END DO
      END DO
!$OMP  END PARALLEL DO
      L=LM
      DO J=J_0S,J_1
         DUT(:,J,L)=DUT(:,J,L)-ASDU(:,J,L-1)*(U(:,J,L-1)+U(:,J,L))
         DVT(:,J,L)=DVT(:,J,L)-ASDU(:,J,L-1)*(V(:,J,L-1)+V(:,J,L))
      END DO
!C**** CALL DIAGNOSTICS
!         IF(MODD5K.LT.MRCH) CALL DIAG5D (4,MRCH,DUT,DVT)
!         IF(MRCH.GT.0) CALL DIAGCD (grid,1,U,V,DUT,DVT,DT1)!,PIT)
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO L=1,LM
      DO J=J_0S,J_1
      DO I=1,IM
        UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)
        VT(I,J,L)=VT(I,J,L)+DVT(I,J,L)
        DUT(I,J,L)=0.
        DVT(I,J,L)=0.
      END DO
      END DO
      END DO
!$OMP  END PARALLEL DO
C****
C**** CORIOLIS FORCE
C****
        CALL HALO_UPDATE(GRID,P ,FROM=SOUTH+NORTH)
c$$$!$OMP  PARALLEL DO PRIVATE(I,IM1,J,L,FD,PDT4,ALPH)
      DO L=1,LM
        IM1=IM
        DO I=1,IM
C         FD(I,1)=FCOR(1)*2.  -.5*(SPA(IM1,1,L)+SPA(I,1,L))*DXV(2)
C         FD(I,JM)=FCOR(JM)*2.+.5*(SPA(IM1,JM,L)+SPA(I,JM,L))*DXV(JM)
C****     Set the Coriolis term to zero at the Poles:
          IF(haveLatitude(grid,J=1))
     *    FD(I,1)= -.5*(SPA(IM1,1,L)+
     *                 SPA(I,1,L))*
     *                 DXV(2)
          IF(haveLatitude(grid,J=JM))
     *    FD(I,JM)=  .5*(SPA(IM1,JM,L)+
     *                 SPA(I,JM,L))*
     *                 DXV(JM)
          IM1=I
        END DO

        DO J=J_0S,J_1S
          IM1=IM
          DO I=1,IM
            FD(I,J)=FCOR(J)+.25*(SPA(IM1,J,L)+SPA(I,J,L))
     *             *(DXV(J)-DXV(J+1))
            IM1=I
          END DO
        END DO
        CALL HALO_UPDATE(GRID,FD,FROM=SOUTH)
        DO J=J_0S,J_1
          IM1=IM
          DO I=1,IM
            PDT4=DT8*(P(I,J-1,L)+P(I,J,L))
            IF(L.GE.LS1) PDT4=DT4*PSFMPT
            ALPH=PDT4*(FD(I,J)+FD(I,J-1))*DSIG(L)
            DUT(I,J,L)=DUT(I,J,L)+ALPH*V(I,J,L)
            DUT(IM1,J,L)=DUT(IM1,J,L)+ALPH*V(IM1,J,L)
            DVT(I,J,L)=DVT(I,J,L)-ALPH*U(I,J,L)
            DVT(IM1,J,L)=DVT(IM1,J,L)-ALPH*U(IM1,J,L)
            IM1=I
          END DO
        END DO
      END DO
c$$$!$OMP  END PARALLEL DO

      if(do_polefix.eq.1) then
c apply the full coriolis force at the pole and ignore the metric term
c which has already been included in advective form
         do ipole=1,2
            if(haveLatitude(grid,J=2) .and. ipole.eq.1) then
               jpo = 1
               jns = jpo + 1
               j = 2
            else if(haveLatitude(grid,J=JM) .and. ipole.eq.2) then
               jpo = JM
               jns = jpo - 1
               j = JM 
            else
               cycle
            endif
            do l=1,lm
               dut(:,j,l) = 0.
               dvt(:,j,l) = 0.
               im1=im
               do i=1,im
                  pdt4=dt8*(p(i,jpo,l)+p(i,jns,l))
                  if(l.ge.ls1) pdt4=dt4*psfmpt
                  alph=pdt4*(2.*fcor(jpo) + fcor(jns))*dsig(l)
                  dut(i  ,j,l)=dut(i  ,j,l)+alph*v(i  ,j,l)
                  dut(im1,j,l)=dut(im1,j,l)+alph*v(im1,j,l)
                  dvt(i  ,j,l)=dvt(i  ,j,l)-alph*u(i  ,j,l)
                  dvt(im1,j,l)=dvt(im1,j,l)-alph*u(im1,j,l)
                  im1=i
               enddo
            enddo
         enddo
      endif

!C**** CALL DIAGNOSTICS
!         IF(MODD5K.LT.MRCH) CALL DIAG5D (5,MRCH,DUT,DVT)
!         IF(MRCH.GT.0) CALL DIAGCD (grid,2,U,V,DUT,DVT,DT1)
C****
C**** ADD CORIOLIS FORCE INCREMENTS TO UT AND VT
C**** AND UNDO SCALING PERFORMED AT BEGINNING OF ADVECV
C****
      CALL HALO_UPDATE(GRID,PB,FROM=SOUTH)
      DO J=J_0S,J_1
        I=IM
        DO IP1=1,IM
          VMASS2(I,J)=.5*((PB(I,J-1)+PB(IP1,J-1))*DXYN(J-1)
     *                 +(PB(I,J)+PB(IP1,J))*DXYS(J))
          I=IP1
        END DO
      END DO
!$OMP  PARALLEL DO PRIVATE(J,L,RVMASS)
      DO L=1,LM
        IF(L.LT.LS1) THEN  !  DO L=1,LS1-1
          DO J=J_0S,J_1
            VT(:,J,L)=(VT(:,J,L)+DVT(:,J,L))/(VMASS2(:,J)*DSIG(L))
            UT(:,J,L)=(UT(:,J,L)+DUT(:,J,L))/(VMASS2(:,J)*DSIG(L))
            DUT(:,J,L)=0.0
            DVT(:,J,L)=0.0
          END DO
        ELSE               !  DO L=LS1,LM
          DO J=J_0S,J_1
            RVMASS=1./(SMASS(J)*DSIG(L))
            UT(:,J,L)=(UT(:,J,L)+DUT(:,J,L))*RVMASS
            VT(:,J,L)=(VT(:,J,L)+DVT(:,J,L))*RVMASS
            DUT(:,J,L)=0.0
            DVT(:,J,L)=0.0
          END DO
        END IF
      END DO
!$OMP  END PARALLEL DO
C
      RETURN
      END SUBROUTINE ADVECV

      end module MOMENTS
