#include "rundeck_opts.h"

      module ATMDYN
      implicit none
      private

      public init_ATMDYN, DYNAM
     &     ,FILTER,AVRX
!     &     ,COMPUTE_DYNAM_AIJ_DIAGNOSTICS
     &     ,AFLUX !, COMPUTE_MASS_FLUX_DIAGS
     &     ,SDRAG,shap1,shap1d
     &     ,fcuva,fcuvb

C**** Variables used in DIAG5 calculations
!@var FCUVA,FCUVB fourier coefficients for velocities
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: FCUVA,FCUVB

      contains

      SUBROUTINE init_ATMDYN
      use domain_decomp_1d, only : grid
      use model_com, only : imh,lm
      CALL AVRX
      ALLOCATE( FCUVA(0:IMH, grid%j_strt_halo:grid%j_stop_halo, LM, 2),
     &          FCUVB(0:IMH, grid%j_strt_halo:grid%j_stop_halo, LM, 2))
      end SUBROUTINE init_ATMDYN


      SUBROUTINE DYNAM
!@sum  DYNAM Integrate dynamic terms
!@auth Original development team
!@ver  1.0
      USE CONSTANT, only : by3,sha,mb2kg,rgas,bygrav
      USE MODEL_COM, only : im,jm,lm,u,v,t,p,q,wm,NIdyn,dt,MODD5K
     *     ,NSTEP,NDA5K,ndaa,mrch,ls1,byim,QUVfilter,DTsrc,USE_UNR_DRAG
      USE GEOM, only : dyv,dxv,dxyp,areag,bydxyp
      USE SOMTQ_COM, only : tmom,mz
      USE DYNAMICS, only : ptold,pu,pv,sd,phi,dut,dvt
     &    ,pua,pva,sda,ps,mb,pk,pmid,pedn
     &    ,cos_limit
      !USE DIAG_COM, only : aij => aij_loc,ij_fmv,ij_fgzv
      USE DOMAIN_DECOMP_1D, only : grid, GET
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, GLOBALSUM
      USE DOMAIN_DECOMP_1D, only : NORTH, SOUTH
      USE DOMAIN_DECOMP_1D, only : haveLatitude
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
      USE MOMENTS, only : advecv

      IMPLICIT NONE

      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     PA, PB, PC, FPEU, FPEV, AM1, AM2
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     &     UT,VT,TT,TZ,TZT,MA,
     &     UX,VX,PIJL
     &    ,UNRDRAG_x,UNRDRAG_y

      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     PUSUM_filter

      REAL*8 DTFS,DTLF, DAMSUM
      INTEGER I,J,L,IP1,IM1   !@var I,J,L,IP1,IM1  loop variables
      INTEGER NS, NSOLD,MODDA    !? ,NIdynO

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

!?    NIdynO=MOD(NIdyn,2)   ! NIdyn odd is currently not an option
      DTFS=DT*2./3.
      DTLF=2.*DT
      NS=0
      NSOLD=0                            ! strat

!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         PUA(:,:,L) = 0.
         PVA(:,:,L) = 0.
         SDA(:,:,L) = 0.
      ENDDO
!$OMP  END PARALLEL DO
C**** Leap-frog re-initialization: IF (NS.LT.NIdyn)
  300 CONTINUE
c     UX(:,:,:)  = U(:,:,:)
c     UT(:,:,:)  = U(:,:,:)
c     VX(:,:,:)  = V(:,:,:)
c     VT(:,:,:)  = V(:,:,:)
!     copy z-moment of temperature into contiguous memory
c     tz(:,:,:) = tmom(mz,:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         UX(:,:,L)  = U(:,:,L)
         UT(:,:,L)  = U(:,:,L)
         VX(:,:,L)  = V(:,:,L)
         VT(:,:,L)  = V(:,:,L)
         TZ(:,:,L)  = TMOM(MZ,:,:,L)
      ENDDO
!$OMP  END PARALLEL DO
C
#ifdef NUDGE_ON
      CALL NUDGE_PREP
#endif
      PA(:,:) = P(:,:)
      PB(:,:) = P(:,:)
      PC(:,:) = P(:,:)

      PUSUM_filter = 0.

C**** INITIAL FORWARD STEP, QX = Q + .667*DT*F(Q)
      MRCH=0
C     CALL DYNAM (UX,VX,TX,PX,Q,U,V,T,P,Q,DTFS)
      CALL CALC_PIJL(LM,P,PIJL)
      CALL AFLUX (U,V,PIJL,PUSUM_filter)
      CALL ADVECM (P,PB,DTFS)
      CALL GWDRAG (PB,UX,VX,U,V,T,TZ,DTFS,.true.)   ! strat
#ifdef NUDGE_ON
      CALL NUDGE (UX,VX,DTFS)
#endif
      CALL VDIFF (PB,UX,VX,U,V,T,DTFS)       ! strat
      CALL ADVECV (P,UX,VX,PB,U,V,Pijl,DTFS)  !P->pijl
      CALL PGF (UX,VX,PB,U,V,T,TZ,Pijl,DTFS)
c      if (QUVfilter) CALL FLTRUV(UX,VX,U,V)
      call isotropuv(ux,vx,COS_LIMIT)
C**** INITIAL BACKWARD STEP IS ODD, QT = Q + DT*F(QX)
      MRCH=-1
C     CALL DYNAM (UT,VT,TT,PT,QT,UX,VX,TX,PX,Q,DT)
      CALL CALC_PIJL(LS1-1,PB,PIJL)
      CALL AFLUX (UX,VX,PIJL,PUSUM_filter)
      CALL ADVECM (P,PA,DT)
#ifdef NUDGE_ON
      CALL NUDGE  (UT,VT,DT)
#endif
      CALL GWDRAG (PA,UT,VT,UX,VX,T,TZ,DT,.false.)   ! strat
      CALL VDIFF (PA,UT,VT,UX,VX,T,DT)       ! strat
      CALL ADVECV (P,UT,VT,PA,UX,VX,Pijl,DT)   !PB->pijl
      CALL PGF (UT,VT,PA,UX,VX,T,TZ,Pijl,DT)

c      if (QUVfilter) CALL FLTRUV(UT,VT,UX,VX)
      call isotropuv(ut,vt,COS_LIMIT)
      GO TO 360
C**** ODD LEAP FROG STEP, QT = QT + 2*DT*F(Q)
  340 MRCH=-2
C     CALL DYNAM (UT,VT,TT,PT,QT,U,V,T,P,Q,DTLF)
      CALL CALC_PIJL(LS1-1,P,PIJL)
#ifdef MASSFLUX_SLPFIL
      CALL PFILTER(pa,PUSUM_filter)
#endif
      CALL AFLUX (U,V,PIJL,PUSUM_filter)
      CALL ADVECM (PA,PB,DTLF)
#ifdef NUDGE_ON
      CALL NUDGE  (UT,VT,DTLF)
#endif
      CALL GWDRAG (PB,UT,VT,U,V,T,TZ,DTLF,.false.)   ! strat
      CALL VDIFF (PB,UT,VT,U,V,T,DTLF)       ! strat
      CALL ADVECV (PA,UT,VT,PB,U,V,Pijl,DTLF)   !P->pijl
      CALL PGF (UT,VT,PB,U,V,T,TZ,Pijl,DTLF)
c      if (QUVfilter) CALL FLTRUV(UT,VT,U,V)
      call isotropuv(ut,vt,COS_LIMIT)
      PA(:,:) = PB(:,:)     ! LOAD PB TO PA
C**** EVEN LEAP FROG STEP, Q = Q + 2*DT*F(QT)
  360 NS=NS+2
         MODD5K=MOD(NSTEP+NS-NIdyn+NDA5K*NIdyn+2,NDA5K*NIdyn+2)
      MRCH=2
C     CALL DYNAM (U,V,T,P,Q,UT,VT,TT,PT,QT,DTLF)
      CALL CALC_PIJL(LS1-1,PA,PIJL)
#ifdef MASSFLUX_SLPFIL
      CALL PFILTER(pc,PUSUM_filter)
#endif
      CALL AFLUX (UT,VT,PIJL,PUSUM_filter)
      CALL ADVECM (PC,P,DTLF)
#ifdef NUDGE_ON
      CALL NUDGE  (U,V,DTLF)
#endif
      CALL GWDRAG (P,U,V,UT,VT,T,TZ,DTLF,.false.)   ! strat
      CALL ADVECV (PC,U,V,P,UT,VT,Pijl,DTLF)     !PA->pijl
         MODDA=MOD(NSTEP+NS-NIdyn+NDAA*NIdyn+2,NDAA*NIdyn+2)  ! strat
!         IF(MODDA.LT.MRCH) CALL DIAGA0   ! strat
C**** ACCUMULATE MASS FLUXES FOR TRACERS and Q
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
        PUA(:,:,L)=PUA(:,:,L)+PU(:,:,L)
        PVA(:,:,L)=PVA(:,:,L)+PV(:,:,L)
        IF (L.LE.LM-1) SDA(:,:,L)=SDA(:,:,L)+SD(:,:,L)
      END DO
!$OMP  END PARALLEL DO

C**** ADVECT Q AND T
CCC   TT(:,:,:) = T(:,:,:)
CCC   TZT(:,:,:)= TZ(:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         TT(:,:,L)  = T(:,:,L)
         TZT(:,:,L) = TZ(:,:,L)
      ENDDO
!$OMP  END PARALLEL DO
      call calc_amp(pc,ma)
      CALL AADVT (MA,T,TMOM, SD,PU,PV, DTLF,.FALSE.,FPEU,FPEV)
!     save z-moment of temperature in contiguous memory for later
CCC   tz(:,:,:) = tmom(mz,:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         TZ(:,:,L) = TMOM(MZ,:,:,L)
      ENDDO
!$OMP  END PARALLEL DO
      CALL VDIFF (P,U,V,UT,VT,T,DTLF)          ! strat
      PC(:,:)    = .5*(P(:,:)+PC(:,:))
CCC   TT(:,:,:)  = .5*(T(:,:,:)+TT(:,:,:))
CCC   TZT(:,:,:) = .5*(TZ(:,:,:)+TZT(:,:,:))
!$OMP PARALLEL DO PRIVATE (L)
      DO L=1,LM
         TT(:,:,L)  = .5*(T(:,:,L)+TT(:,:,L))
         TZT(:,:,L) = .5*(TZ(:,:,L)+TZT(:,:,L))
      ENDDO
!$OMP END PARALLEL DO

      CALL CALC_PIJL(LS1-1,PC,PIJL)
c      CALL CALC_PIJL(LS1-1,PA,PIJL) ! true leapfrog
      CALL PGF (U,V,P,UT,VT,TT,TZT,Pijl,DTLF)    !PC->pijl

      !call compute_mass_flux_diags(PHI, PU, PV, dt)

      CALL CALC_AMPK(LS1-1)
c      if (QUVfilter) CALL FLTRUV(U,V,UT,VT)
      call isotropuv(u,v,COS_LIMIT)
      PC(:,:) = P(:,:)      ! LOAD P TO PC
      if (USE_UNR_DRAG==0) CALL SDRAG (DTLF)
         IF (MOD(NSTEP+NS-NIdyn+NDAA*NIdyn+2,NDAA*NIdyn+2).LT.MRCH) THEN
!           CALL DIAGA
!           CALL DIAGB
           CALL EPFLUX (U,V,T,P)
         ENDIF
C**** Restart after 8 steps due to divergence of solutions
      IF (NS-NSOLD.LT.8 .AND. NS.LT.NIdyn) GO TO 340
      NSOLD=NS
      IF (NS.LT.NIdyn) GO TO 300

      if (USE_UNR_DRAG==1) then
      CALL UNRDRAG (P,U,V,T,TZ,UNRDRAG_x,UNRDRAG_y)
      U(:,:,:) = U(:,:,:) + UNRDRAG_x(:,:,:) * DTsrc
      V(:,:,:) = V(:,:,:) + UNRDRAG_y(:,:,:) * DTsrc
      end if

      PUA = PUA * DTLF
      PVA = PVA * DTLF
      SDA(:,:,1:LM-1) = SDA(:,:,1:LM-1) * DTLF

c apply east-west filter to U and V once per physics timestep
      CALL FLTRUV(U,V,UT,VT)
c apply north-south filter to U and V once per physics timestep
      call conserv_amb_ext(u,am1) ! calculate ang. mom. before filter
      call fltry2(u,1d0) ! 2nd arg could be set using DT_YUfilter
      call fltry2(v,1d0) ! 2nd arg could be set using DT_YVfilter
      call conserv_amb_ext(u,am2) ! calculate ang. mom. after filter
      am2(:,j_0stg:j_1stg) = am1(:,j_0stg:j_1stg)-am2(:,j_0stg:j_1stg)
      if(have_south_pole) am2(:,1) = 0.
      call globalsum(grid,am2,damsum,all=.true.)
      call add_am_as_solidbody_rotation(u,damsum) ! maintain global ang. mom.

      RETURN
      END SUBROUTINE DYNAM

!      Subroutine compute_mass_flux_diags(PHI, PU, PV, dt)
!      use MODEL_COM, only: IM, LM
!      use DOMAIN_DECOMP_1D, only: grid, halo_update, SOUTH, get
!      use DIAG_COM, only: AIJ => AIJ_loc, IJ_FGZU, IJ_FGZV
!
!      real*8, intent(inout) :: PHI(:,grid%J_STRT_HALO:,:)
!      real*8, intent(in) :: PU(:,grid%J_STRT_HALO:,:)
!      real*8, intent(in) :: PV(:,grid%J_STRT_HALO:,:)
!      real*8, intent(in) :: dt
!
!
!      integer :: J_0S, J_1S
!      integer :: J_0STG, J_1STG
!      integer :: I, IP1, L, J
!
!      call get(grid, J_STRT_STGR=J_0STG,J_STOP_STGR=J_1STG,
!     &               J_STRT_SKP =J_0S,  J_STOP_SKP =J_1S)
!
!      CALL HALO_UPDATE(grid, PHI, FROM=SOUTH)
!!$OMP  PARALLEL DO PRIVATE (J,L,I,IP1)
!      DO J=J_0S,J_1S ! eastward transports
!      DO L=1,LM
!         I=IM
!         DO IP1=1,IM
!            AIJ(I,J,IJ_FGZU)=AIJ(I,J,IJ_FGZU)+
!     &           (PHI(I,J,L)+PHI(IP1,J,L))*PU(I,J,L)*DT ! use DT=DTLF/2
!            I=IP1
!         END DO
!      END DO
!      END DO
!!$OMP  END PARALLEL DO
!!$OMP  PARALLEL DO PRIVATE (J,L,I)
!      DO J=J_0STG,J_1STG ! northward transports
!      DO L=1,LM
!         DO I=1,IM
!            AIJ(I,J,IJ_FGZV)=AIJ(I,J,IJ_FGZV)+
!     &           (PHI(I,J-1,L)+PHI(I,J,L))*PV(I,J,L)*DT ! use DT=DTLF/2
!         END DO
!      END DO
!      END DO
!!$OMP  END PARALLEL DO
!
!      end subroutine compute_mass_flux_diags
!
!      subroutine COMPUTE_DYNAM_AIJ_DIAGNOSTICS(
!     &    PUA, PVA, dt)
!      use CONSTANT,      only: BY3
!      use DOMAIN_DECOMP_1D, only: grid, get, halo_update, SOUTH
!      USE DOMAIN_DECOMP_1D, only : haveLatitude
!      use DIAG_COM, only: AIJ => AIJ_loc,
!     &     IJ_FGZU, IJ_FGZV, IJ_FMV, IJ_FMU
!      use MODEL_COM, only: IM,JM,LM
!
!      real*8, intent(in) :: PUA(:,grid%J_STRT_HALO:,:)
!      real*8, intent(in) :: PVA(:,grid%J_STRT_HALO:,:)
!      real*8, intent(in) :: dt
!
!      integer :: I, IP1, J, L
!      integer :: J_0STG, J_1STG, J_0S, J_1S
!      logical :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE
!      real*8 :: dtlf
!
!      dtlf = 2.*dt
!
!      call get(grid, J_STRT_STGR=J_0STG,J_STOP_STGR=J_1STG,
!     &               J_STRT_SKP =J_0S,  J_STOP_SKP =J_1S,
!     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE,
!     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)
!
!!$OMP  PARALLEL DO PRIVATE (J,L)
!      do j=J_0STG,J_1STG
!      do L=1,LM
!         AIJ(:,J,IJ_FMV)  = AIJ(:,J,IJ_FMV )+PVA(:,J,L)*DTLF
!      enddo
!      enddo
!!$OMP  END PARALLEL DO
!
!!$OMP  PARALLEL DO PRIVATE (J,L)
!      do j=J_0S,J_1S
!      do l=1,lm
!         AIJ(:,J,IJ_FMU) = AIJ(:,J,IJ_FMU)+PUA(:,J,L)*DTLF
!      enddo
!      enddo
!!$OMP  END PARALLEL DO
!
!      if (haveLatitude(grid, J=1)) then
!         do l=1,lm
!            AIJ(:,1,IJ_FMU)  = AIJ(:, 1,IJ_FMU )+PUA(:, 1,L)*DTLF*BY3
!         enddo
!      endif
!      if(haveLatitude(grid, J=JM)) then
!         do l=1,lm
!            AIJ(:,JM,IJ_FMU) = AIJ(:,JM,IJ_FMU )+PUA(:,JM,L)*DTLF*BY3
!         enddo
!      endif
!      end subroutine COMPUTE_DYNAM_AIJ_DIAGNOSTICS

      SUBROUTINE AFLUX (U,V,PIJL,PUSUM_filter)
!@sum  AFLUX Calculates horizontal/vertical air mass fluxes
!@+    Input: U,V velocities, PIJL pressure
!@+    Output: PIT  pressure tendency (mb m^2/s)
!@+            SD   sigma dot (mb m^2/s)
!@+            PU,PV horizontal mass fluxes (mb m^2/s)
!@+            CONV  horizontal mass convergence (mb m^2/s)
!@+            SPA
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,imh,jm,lm,ls1,dsig,bydsig,byim
     &     ,zatmo,sige,do_polefix
      USE GEOM, only : dyp,dxv,polwt,imaxj
      USE DYNAMICS, only : pit,sd,conv,pu,pv,spa
      USE DOMAIN_DECOMP_1D, only : grid, GET
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : NORTH, SOUTH
      USE DOMAIN_DECOMP_1D, only : haveLatitude
      IMPLICIT NONE
C**** CONSTANT PRESSURE AT L=LS1 AND ABOVE, PU,PV CONTAIN DSIG
!@var U,V input velocities (m/s)
!@var PIJL input 3-D pressure field (mb) (no DSIG)
      REAL*8, INTENT(INOUT),
     &  DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LM) :: U,V
      REAL*8, INTENT(INOUT),
     &  DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LM) :: PIJL
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     PUSUM_filter
      REAL*8, DIMENSION(IM) :: DUMMYS,DUMMYN
      INTEGER I,J,L,IP1,IM1,IPOLE
      REAL*8 PUS,PUN,PVS,PVN,PBS,PBN
      REAL*8 WT,DXDSIG,DYDSIG,PVSA(LM),PVNA(LM),xx,twoby3
      real*8, dimension(im,2,LM) :: usv0,vsv0
      integer :: jvs,jvn,jv
      real*8 :: wts
      real*8, dimension(lm) :: dsig_or_zero
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      dsig_or_zero(1:ls1-1) = dsig(1:ls1-1)
      dsig_or_zero(ls1:lm) = 0.

C****
C**** BEGINNING OF LAYER LOOP
C****
      CALL HALO_UPDATE(grid, U,FROM=NORTH+SOUTH)
      CALL HALO_UPDATE(grid, V,FROM=NORTH+SOUTH)

c
c interpolate polar velocities to the appropriate latitude
c
      do ipole=1,2
        if (haveLatitude(grid, J=2) .and. ipole == 1) then
          jv  = 2
          jvs = 2          ! jvs is the southernmost velocity row
          jvn = jvs+1      ! jvs is the northernmost velocity row
          wts = polwt
        else if(haveLatitude(grid,JM) .and. ipole == 2) then
          jv = JM
          jvs = jv - 1
          jvn = jvs + 1
          wts = 1.-polwt
        else
          cycle
        endif

        do L = 1, LM
          usv0(:,ipole,l) = u(:,jv,l)
          vsv0(:,ipole,l) = v(:,jv,l)
          u(:,jv,l) = wts*u(:,jvs,l) + (1.-wts)*u(:,jvn,l)
          v(:,jv,l) = wts*v(:,jvs,l) + (1.-wts)*v(:,jvn,l)
        end do
      enddo

C****
C**** COMPUTATION OF MASS FLUXES     P,T  PU     PRIMARY GRID ROW
C**** ARAKAWA'S SCHEME B             PV   U,V    SECONDARY GRID ROW
C****
C**** COMPUTE PU, THE WEST-EAST MASS FLUX, AT NON-POLAR POINTS
      do L = 1, LM
      DO 2154 J=J_0S,J_1S
      DO 2154 I=1,IM
 2154 SPA(I,J,L)=U(I,J,L)+U(I,J+1,L)
      CALL AVRX (SPA(1,J_0H,L))
      I=IM
      DO 2166 J=J_0S,J_1S
      DYDSIG = 0.25D0*DYP(J)*DSIG(L)
      DO 2165 IP1=1,IM
      PU(I,J,L)=DYDSIG*SPA(I,J,L)*(PIJL(I,J,L)+PIJL(IP1,J,L))
     &       +pusum_filter(i,j)*dsig_or_zero(l)
 2165 I=IP1
 2166 CONTINUE
      end do
C**** COMPUTE PV, THE SOUTH-NORTH MASS FLUX
      CALL HALO_UPDATE(grid, PIJL, FROM=SOUTH)
      do L = 1, LM
      IM1=IM
      DO 2172 J=J_0STG, J_1STG
      DXDSIG = 0.25D0*DXV(J)*DSIG(L)
      DO 2170 I=1,IM
      PV(I,J,L)=DXDSIG*(V(I,J,L)+V(IM1,J,L))*(PIJL(I,J,L)+PIJL(I,J-1,L))

 2170 IM1=I
 2172 CONTINUE
      end do

c restore uninterpolated values of u,v at the pole
      do ipole=1,2
        if (haveLatitude(grid, J=2) .and. ipole == 1) then
          jv = 2           ! why not staggered grid
        else if(haveLatitude(grid, J=JM) .and. ipole.eq.2) then
          jv = JM            ! why not staggered grid
        else
          cycle
        endif
        do L = 1, LM
          u(:,jv,l) = usv0(:,ipole,l)
          v(:,jv,l) = vsv0(:,ipole,l)
        end do
      enddo

C**** COMPUTE PU*3 AT THE POLES
      IF (haveLatitude(grid, J=1)) Then
        do L = 1, LM
        PUS=0.
        PVS=0.
        DO I=1,IM
          PUS=PUS+U(I,2,L)
          PVS=PVS+PV(I,2,L)
        END DO
        PUS=.25*DYP(2)*PUS*PIJL(1,1,L)*BYIM
        PVS=PVS*BYIM
        PVSA(L)=PVS
        DUMMYS(1)=0.
        DO I=2,IM
          DUMMYS(I)=DUMMYS(I-1)+(PV(I,2,L) -PVS)*BYDSIG(L)
        END DO
        PBS=0.
        PBN=0.
        DO I=1,IM
          PBS=PBS+DUMMYS(I)
        END DO
        PBS=PBS*BYIM
        DO I=1,IM
          SPA(I,1,L)=4.*(PBS-DUMMYS(I)+PUS)/(DYP(2)*PIJL(1,1,L))
          PU(I,1,L)=3.*(PBS-DUMMYS(I)+PUS)*DSIG(L)
        END DO
        end do
      END IF

      IF (haveLatitude(grid, J=JM)) THEN
        do L = 1,LM
        PUN=0.
        PVN=0.
        DO I=1,IM
          PUN=PUN+U(I,JM,L)
          PVN=PVN+PV(I,JM,L)
        END DO
        PUN=.25*DYP(JM-1)*PUN*PIJL(1,JM,L)*BYIM
        PVN=PVN*BYIM
        PVNA(L)=PVN
        DUMMYN(1)=0.
        DO I=2,IM
          DUMMYN(I)=DUMMYN(I-1)+(PV(I,JM,L)-PVN)*BYDSIG(L)
        END DO
        PBN=0.
        DO I=1,IM
          PBN=PBN+DUMMYN(I)
        END DO
        PBN=PBN*BYIM
        DO  I=1,IM
          SPA(I,JM,L)=4.*(DUMMYN(I)-PBN+PUN)/(DYP(JM-1)*PIJL(1,JM,L))
          PU(I,JM,L)=3.*(DUMMYN(I)-PBN+PUN)*DSIG(L)
        END DO
        END DO
      END IF
C****
C**** CONTINUITY EQUATION
C****
C**** COMPUTE CONV, THE HORIZONTAL MASS CONVERGENCE
c     CALL HALO_UPDATE(grid, PV, FROM=NORTH)
c     DO 1510 J=J_0S,J_1S
c     IM1=IM
c     DO 1510 I=1,IM
c     CONV(I,J,L)=(PU(IM1,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))
c1510 IM1=I
c     IF (HAVE_SOUTH_POLE) CONV(1,1,L)=-PVS
c     IF (HAVE_NORTH_POLE) CONV(1,JM,L)=PVN

      if(do_polefix.eq.1) then
c To maintain consistency with subroutine ADVECV,
c adjust pu at the pole if no corner fluxes are used there
c in ADVECV.
      do ipole=1,2
        if(haveLatitude(grid,J=1) .and. ipole.eq.1) then
          j = 1
        else if(haveLatitude(grid,J=JM) .and. ipole.eq.2) then
          j = JM
        else
          cycle
        endif
        twoby3 = 2d0/3d0
        do l=1,lm
           pu(:,j,l) = pu(:,j,l)*twoby3
        enddo
      enddo
      endif
C****
C**** END OF HORIZONTAL ADVECTION LAYER LOOP
C****
c
c modify uphill air mass fluxes around steep topography

      !DO J=90,1,-1
      !WRITE(6,'(A7,144(F8.0))') 'ZATMO: ', ZATMO(1:IM,J) 
      !END DO
      !CALL FLUSH(6)
      !STOP

      do 2015 j=J_0S,J_1S
      i = im
      do 2010 ip1=1,im
         xx = zatmo(ip1,j)-zatmo(i,j)
         if(xx.eq.0.0)  go to 2007
         DO 2005 L=1,LS1-1
ccc         if((zatmo(ip1,j)-zatmo(i,j))*pu(i,j,l).gt.0.) then
            if(xx*pu(i,j,l).gt.0.) then
               if(pu(i,j,l).gt.0.) then
                  wt = (pijl(ip1,j,l)/pijl(i,j,l)-sige(l+1))/dsig(l)
               else
                  wt = (pijl(i,j,l)/pijl(ip1,j,l)-sige(l+1))/dsig(l)
               endif
               if(wt.le.0.) then
                  pu(i,j,l+1) = pu(i,j,l+1) + pu(i,j,l)
                  pu(i,j,l) = 0.
               else
                  go to 2007
               endif
            endif
 2005    CONTINUE
 2007    CONTINUE
         i = ip1
 2010 CONTINUE
 2015 CONTINUE

ccc   do j=2,jm
c**** Exceptional J loop boundaries
      do 2035 j=max(J_0S,3), J_1S
      do 2035 i=1,im
         xx = zatmo(i,j)-zatmo(i,j-1)
         if(xx.eq.0.0)  go to 2035
         DO 2020 L=1,LS1-1
ccc         if((zatmo(i,j)-zatmo(i,j-1))*pv(i,j,l).gt.0.) then
            if(xx*pv(i,j,l).gt.0.) then
               if(pv(i,j,l).gt.0.) then
                  wt = (pijl(i,j,l)/pijl(i,j-1,l)-sige(l+1))/dsig(l)
               else
                  wt = (pijl(i,j-1,l)/pijl(i,j,l)-sige(l+1))/dsig(l)
               endif
               if(wt.le.0.) then
                  pv(i,j,l+1) = pv(i,j,l+1) + pv(i,j,l)
                  pv(i,j,l) = 0.
               else
                  go to 2035
               endif
            endif
 2020    CONTINUE
 2035 CONTINUE
C
C     Now Really Do  CONTINUITY EQUATION
C
C     COMPUTE CONV, THE HORIZONTAL MASS CONVERGENCE
C
      CALL HALO_UPDATE(GRID,PU ,FROM=SOUTH+NORTH) ! full halos needed later
      CALL HALO_UPDATE(GRID,PV ,FROM=SOUTH+NORTH)
!$OMP  PARALLEL DO PRIVATE (I,J,L,IM1)
      DO 2400 L=1,LM
      DO 1510 J=J_0S,J_1S
      IM1=IM
      DO 1510 I=1,IM
      CONV(I,J,L)=(PU(IM1,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))
 1510 IM1=I
      IF (haveLatitude(grid,J=1)) CONV(1,1,L)=-PVSA(L)
      If (haveLatitude(grid,J=JM)) CONV(1,JM,L)=PVNA(L)
 2400 CONTINUE
!$OMP  END PARALLEL DO
C
C**** COMPUTE PIT, THE PRESSURE TENDENCY
      PIT(:,J_0:J_1) = CONV(:,J_0:J_1,1)
      SD(:,J_0:J_1,1:LM-1) = CONV(:,J_0:J_1,2:LM)
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO 2420 J=J_0,J_1
         DO 2410 I=1,IMAXJ(J)
         DO 2410 L=LM-1,1,-1
           PIT(I,J) = PIT(I,J) + SD(I,J,L)
 2410    CONTINUE
 2420 CONTINUE
!$OMP  END PARALLEL DO
C**** COMPUTE SD, SIGMA DOT                                             -------
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO 2435 J=J_0,J_1
         DO 2430 I=1,IMAXJ(J)
         DO 2430 L=LM-2,LS1-1,-1
            SD(I,J,L)=SD(I,J,L+1)+SD(I,J,L)
 2430    CONTINUE
 2435 CONTINUE
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO 2440 J=J_0,J_1
         DO 2438 I=1,IMAXJ(J)
         DO 2438 L=LS1-2,1,-1
           SD(I,J,L)=SD(I,J,L+1)+SD(I,J,L)-
     &           DSIG(L+1)*PIT(I,J)
 2438    CONTINUE
 2440 CONTINUE
!$OMP  END PARALLEL DO
      DO 2450 L=1,LM-1
      DO 2450 I=2,IM
        IF (haveLatitude(grid,J=1))
     *     SD(I,1,L)=SD(1,1,L)
 2450   IF (haveLatitude(grid,J=JM))
     *       SD(I,JM,L)=SD(1,JM,L)

        ! Recopy into CONV to support prior usage
      CONV(:,J_0:J_1,1) = PIT(:,J_0:J_1)
      CONV(:,J_0:J_1,2:LM) = SD(:,J_0:J_1,1:LM-1)

      RETURN
      END SUBROUTINE AFLUX


      SUBROUTINE ADVECM (P,PA,DT1)
!@sum  ADVECM Calculates updated column pressures using mass fluxes
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,mrch,zatmo,u,v,t,q,ptop
      USE GEOM, only : bydxyp,imaxj
      USE DYNAMICS, only : pit
      USE DOMAIN_DECOMP_1D, only : grid, GET
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, GLOBALSUM
      USE DOMAIN_DECOMP_1D, only : NORTH, SOUTH
      USE DOMAIN_DECOMP_1D, only : haveLatitude
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: P(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8, INTENT(OUT) :: PA(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8, INTENT(IN) :: DT1
      INTEGER I,J,L  !@var I,J,L  loop variables
      INTEGER IM1 ! @var IM1 = I - 1
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      INTEGER :: n_exception, n_exception_all

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1, J_STRT_HALO = J_0H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE )

C**** COMPUTE PA, THE NEW SURFACE PRESSURE
      ! 1st pass count warning/termination events
      ! This avoids the need for 2 halo fills during normal
      ! execution.
      n_exception = 0
      outer_loop:  DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          PA(I,J)=P(I,J)+(DT1*PIT(I,J)*BYDXYP(J))
          IF (PA(I,J)+PTOP.GT.1160. .or. PA(I,J)+PTOP.LT.350.) THEN
            n_exception = n_exception + 1
            IF (PA(I,J)+PTOP.lt.250. .or. PA(I,J)+PTOP.GT.1200.)
     *           Exit outer_loop
          End If
        END DO
      END DO outer_loop

      Call GLOBALSUM(grid, n_exception, n_exception_all, all=.true.)
      IF (n_exception_all > 0) Then ! need halos
        CALL HALO_UPDATE(grid, U, FROM=NORTH)
        CALL HALO_UPDATE(grid, V, FROM=NORTH)
      END IF

      IF (n_exception > 0)  Then ! 2nd pass report problems
        Do J = J_0, J_1
          DO I = 1, IMAXJ(J)
            IF (PA(I,J)+PTOP.GT.1160. .or. PA(I,J)+PTOP.LT.350.) THEN
              IM1 = 1 + MOD(IM+I-2,IM)
              WRITE (6,990) I,J,MRCH,P(I,J),PA(I,J),ZATMO(I,J),DT1,
     *             (U(IM1,J,L),U(I,J,L),U(IM1,J+1,L),U(I,J+1,L),
     *             V(IM1,J,L),V(I,J,L),V(IM1,J+1,L),V(I,J+1,L),
     *             T(I,J,L),Q(I,J,L),L=1,LM)
              write(6,*) "Pressure diagnostic error"
              IF (PA(I,J)+PTOP.lt.250. .or. PA(I,J)+PTOP.GT.1200.)
     &          call stop_model('ADVECM: Pressure diagnostic error',11)
            END IF
          END DO
        END DO
      END IF

      IF (haveLatitude(grid, J=1)) PA(2:IM, 1)=PA(1,1)
      IF (haveLatitude(grid, J=JM)) PA(2:IM,JM)=PA(1,JM)

C****
      RETURN
  990 FORMAT (/'0PRESSURE DIAGNOSTIC     I,J,MRCH,P,PA=',3I4,2F10.2/
     *  '     ZATMO=',F10.3,' DT=',F6.1/
     *  '0    U(I-1,J)     U(I,J)   U(I-1,J+1)    U(I,J+1)    V(I-1,J)',
     *   '     V(I,J)   V(I-1,J+1)    V(I,J+1)     T(I,J)     Q(I,J)'/
     *  (1X,9F12.3,F12.6))
      END SUBROUTINE ADVECM


      SUBROUTINE PGF (UT,VT,PB,U,V,T,SZ,P,DT1)
!@sum  PGF Adds pressure gradient forces to momentum
!@auth Original development team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,kapa,bykapa,bykapap1,bykapap2
      USE MODEL_COM, only : im,jm,lm,ls1,mrch,dsig,psfmpt,sige,ptop
     *     ,zatmo,sig,modd5k,bydsig
     &     ,do_polefix
      USE GEOM, only : imaxj,dxyv,dxv,dyv,dxyp,dyp,dxp,acor,acor2
      USE DYNAMICS, only : gz,pu,pit,phi,spa,dut,dvt
c      USE DIAG, only : diagcd
      USE DOMAIN_DECOMP_1D, Only : grid, GET
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : NORTH, SOUTH
      USE DOMAIN_DECOMP_1D, only : haveLatitude
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM):: U,V,T
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO):: FD,RFDUX
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *  UT, VT, QT, P, SZ
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *  PB

      REAL*8 PKE(LS1:LM+1)
      REAL*8 DT4,DT1
      REAL*8 PIJ,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP,DP,P0,X
     *     ,BYDP
      REAL*8 TZBYDP,FLUX,FDNP,FDSP,RFDU,PHIDN,FACTOR
      INTEGER I,J,L,IM1,IP1,IPOLE  !@var I,J,IP1,IM1,L,IPOLE loop variab.
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)
C****
      DT4=DT1/4.
      DO L=LS1,LM+1
        PKE(L)=(PSFMPT*SIGE(L)+PTOP)**KAPA
      END DO
C****
C**** VERTICAL DIFFERENCING
C****
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=LS1,LM
      SPA(:,:,L)=0.
      END DO
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO PRIVATE(I,J,L,DP,P0,PIJ,PHIDN,TZBYDP,X,
!$OMP*             BYDP,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP)
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
        PIJ=P(I,J,1)
        PDN=PIJ+PTOP
        PKDN=PDN**KAPA
        PHIDN=ZATMO(I,J)
C**** LOOP OVER THE LAYERS
        DO L=1,LM
          PKPDN=PKDN*PDN
          PKPPDN=PKPDN*PDN
          IF(L.GE.LS1) THEN
            DP=DSIG(L)*PSFMPT
            BYDP=1./DP
            P0=SIG(L)*PSFMPT+PTOP
            TZBYDP=2.*SZ(I,J,L)*BYDP
            X=T(I,J,L)+TZBYDP*P0
            PUP=SIGE(L+1)*PSFMPT+PTOP
            PKUP=PKE(L+1)
            PKPUP=PKUP*PUP
            PKPPUP=PKPUP*PUP
          ELSE
            DP=DSIG(L)*PIJ
            BYDP=1./DP
            P0=SIG(L)*PIJ+PTOP
            TZBYDP=2.*SZ(I,J,L)*BYDP
            X=T(I,J,L)+TZBYDP*P0
            PUP=SIGE(L+1)*PIJ+PTOP
            PKUP=PUP**KAPA
            PKPUP=PKUP*PUP
            PKPPUP=PKPUP*PUP
C****   CALCULATE SPA, MASS WEIGHTED THROUGHOUT THE LAYER
            SPA(I,J,L)=RGAS*((X+TZBYDP*PTOP)*(PKPDN-PKPUP)*BYKAPAP1
     *      -X*PTOP*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2)
     *      *BYDP
          END IF
C**** CALCULATE PHI, MASS WEIGHTED THROUGHOUT THE LAYER
          PHI(I,J,L)=PHIDN+RGAS*(X*PKDN*BYKAPA-TZBYDP*PKPDN*BYKAPAP1
     *      -(X*(PKPDN-PKPUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2)
     *      *BYDP*BYKAPAP1)
C**** CALULATE PHI AT LAYER TOP (EQUAL TO BOTTOM OF NEXT LAYER)
          PHIDN=PHIDN+RGAS*(X*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPDN-PKPUP)
     *     *BYKAPAP1)
          PDN=PUP
          PKDN=PKUP
        END DO
      END DO
      END DO
!$OMP END PARALLEL DO
C**** SET POLAR VALUES FROM THOSE AT I=1
      IF (haveLatitude(grid, J=1)) THEN
        DO L=1,LM
          SPA(2:IM,1,L)=SPA(1,1,L)
          PHI(2:IM,1,L)=PHI(1,1,L)
        END DO
      END IF
      IF (haveLatitude(grid, J=JM)) THEN
        DO L=1,LM
          SPA(2:IM,JM,L)=SPA(1,JM,L)
          PHI(2:IM,JM,L)=PHI(1,JM,L)
        END DO
      END IF

!$OMP  PARALLEL DO PRIVATE(L)
      DO L=1,LM
        GZ(:,:,L)=PHI(:,:,L)
      END DO
!$OMP END PARALLEL DO
C****
C**** PRESSURE GRADIENT FORCE
C****
C**** NORTH-SOUTH DERIVATIVE AFFECTS THE V-COMPONENT OF MOMENTUM
C
      CALL HALO_UPDATE(grid, P,   FROM=SOUTH)
      CALL HALO_UPDATE(grid, PHI, FROM=SOUTH)
      CALL HALO_UPDATE(grid, SPA, FROM=SOUTH)
!$OMP  PARALLEL DO PRIVATE(I,IM1,J,L,FACTOR,FLUX)
      DO 3236 L=1,LM
      DO 3236 J=J_0STG,J_1STG
      FACTOR = DT4*DXV(J)*DSIG(L)
      IM1=IM
      DO 3234 I=1,IM
      FLUX=    ((P(I,J,L)+P(I,J-1,L))*(PHI(I,J,L)-PHI(I,J-1,L))+
     *  (SPA(I,J,L)+SPA(I,J-1,L))*(P(I,J,L)-P(I,J-1,L)))*FACTOR
      DVT(I,J,L)  =DVT(I,J,L)  -FLUX
      DVT(IM1,J,L)=DVT(IM1,J,L)-FLUX
 3234 IM1=I
 3236 CONTINUE
!$OMP  END PARALLEL DO
C
C**** SMOOTHED EAST-WEST DERIVATIVE AFFECTS THE U-COMPONENT
C
C Although PU appears to require a halo update, the halos
C of PHI, SPA, and P enable implementation without the additional halo.
C
!$OMP  PARALLEL DO PRIVATE(I,IP1,J,L,FACTOR)
      DO L=1,LM
        IF (haveLatitude(grid, J=1)) PU(:,1,L)=0.
        IF (haveLatitude(grid, J=JM)) PU(:,JM,L)=0.
        I=IM

        DO J=Max(2,J_0STG-1),J_1STG
          DO IP1=1,IM
            PU(I,J,L)=(P(IP1,J,L)+P(I,J,L))*(PHI(IP1,J,L)-PHI(I,J,L))+
     *           (SPA(IP1,J,L)+SPA(I,J,L))*(P(IP1,J,L)-P(I,J,L))
            I=IP1
          END DO
        END DO

        CALL AVRX (PU(1,J_0H,L),jrange=(/MAX(2,J_0H),MIN(JM-1,J_1H)/))

        DO J=J_0STG,J_1STG
          FACTOR = -DT4*DYV(J)*DSIG(L)
          DO I=1,IM
            DUT(I,J,L)=DUT(I,J,L)+FACTOR*(PU(I,J,L)+PU(I,J-1,L))
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO

c correct for erroneous dxyv at the poles
      if(do_polefix.eq.1) then
         do ipole=1,2
            if(haveLatitude(grid,J=2) .and. ipole.eq.1) then
               j = 2
            else if(haveLatitude(grid,J=JM) .and. ipole.eq.2) then
               j = JM
            else
               cycle
            endif
            dut(:,j,:) = dut(:,j,:)*acor
            dvt(:,j,:) = dvt(:,j,:)*acor2
         enddo
      endif
C
C**** CALL DIAGNOSTICS
      IF(MRCH.GT.0) THEN
!         IF(MODD5K.LT.MRCH) CALL DIAG5D (6,MRCH,DUT,DVT)
!         CALL DIAGCD (grid,3,U,V,DUT,DVT,DT1)
      ENDIF
C****
C****
C**** UNDO SCALING PERFORMED AT BEGINNING OF DYNAM
C****
      DO 3410 J=J_0STG,J_1STG
      DO 3410 I=1,IM
 3410 FD(I,J)=PB(I,J)*DXYP(J)
      IF (haveLatitude(grid, J=1)) THEN
        FDSP=PB(1, 1)*DXYP( 1)
        FDSP=FDSP+FDSP
        DO I=1,IM
          FD(I, 1)=FDSP
        END DO
      END IF
      IF (haveLatitude(grid, J=JM)) THEN
        FDNP=PB(1,JM)*DXYP(JM)
        FDNP=FDNP+FDNP
        DO I=1,IM
          FD(I,JM)=FDNP
        END DO
      END IF
C
      CALL HALO_UPDATE(grid, FD, FROM=SOUTH)
!$OMP  PARALLEL DO PRIVATE(I,IP1,J)
      DO 3530 J=J_0STG,J_1STG
      I=IM
      DO 3525 IP1=1,IM
      RFDUX(I,J)=4./(FD(I,J)+FD(IP1,J)+FD(I,J-1)+FD(IP1,J-1))
 3525 I = IP1
 3530 CONTINUE
!$OMP  END PARALLEL DO
C
!$OMP  PARALLEL DO PRIVATE(I,J,L,RFDU)
      DO 3550 L=1,LM
      DO 3550 J=J_0STG,J_1STG
      RFDU=1./(PSFMPT*DXYV(J)*DSIG(L))
      DO 3540 I=1,IM
      IF(L.LT.LS1) RFDU=RFDUX(I,J)*BYDSIG(L)
      VT(I,J,L)=VT(I,J,L)+DVT(I,J,L)*RFDU
      UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)*RFDU
 3540 CONTINUE
 3550 CONTINUE
!$OMP  END PARALLEL DO
C
      RETURN
      END SUBROUTINE PGF

      SUBROUTINE AVRX(X,jrange)
!@sum  AVRX Smoothes zonal mass flux and geopotential near the poles
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,imh
      USE GEOM, only : dlon,dxp,dyp,bydyp
      !USE DYNAMICS, only : xAVRX
C**** THIS VERSION OF AVRX DOES SO BY TRUNCATING THE FOURIER SERIES.
      USE DOMAIN_DECOMP_1D, Only : grid, GET
      USE MOMENTS, only : moment_enq_order
      USE constant, only : byrt2
      IMPLICIT NONE
      REAL*8, INTENT(INOUT), optional ::
     &     X(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
      Integer, Intent(In), optional :: jrange(2)
      REAL*8, ALLOCATABLE, SAVE  :: DRAT(:)
      REAL*8, SAVE ::  BYSN(IMH)
      REAL*8, DIMENSION(0:IMH) :: AN,BN
CCC   INTEGER, SAVE :: NMIN(grid%J_STRT_HALO:grid%J_STOP_HALO)
CCC   INTEGER, SAVE :: IFIRST = 1
      INTEGER, ALLOCATABLE, SAVE :: NMIN(:)
      INTEGER J,N
      LOGICAL, SAVE :: init = .false.
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H, J0, J1
      REAL*8, SAVE :: xAVRX
      INTEGER order

      if ( present(X) ) goto 1000

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_HALO = J_0H, J_STOP_HALO = J_1H,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)
C
      call moment_enq_order(order)
      if ( order==4 ) then
        xAVRX = byrt2
      else if ( order==2 ) then
        xAVRX = 1.d0
      else
        call stop_model("unsupported scheme order in AVRX0",255)
      endif
C
      IF (.NOT. init) THEN
        init = .true.
C       CALL FFT0(IM)
        j0 = MAX(1,J_0H)
        j1 = MIN(JM,J_1H)
        ALLOCATE(DRAT(j0:j1), NMIN(j0:j1))
        DO N=1,IMH
          BYSN(N)=xAVRX/SIN(.5*DLON*N)
        END DO
        DO J=j0,j1
          DRAT(J) = DXP(J)*BYDYP(3)
          DO N=IMH,1,-1
            IF(BYSN(N)*DRAT(J) .GT.1.) THEN
              NMIN(J) = N+1
              EXIT
            ENDIF
          END DO
        END DO
      END IF
      RETURN
C****
!!!      ENTRY AVRX (X)
 1000 continue
C****

      If (Present(jrange)) Then
        j0 = jrange(1)
        j1 = jrange(2)
      Else
        CALL GET(grid, J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)
        j0=J_0S
        j1=J_1S
      End If

      DO J=j0,j1
        IF (DRAT(J).GT.1) CYCLE
        CALL FFT (X(1,J),AN,BN)
        DO N=NMIN(J),IMH-1
          AN(N)=BYSN(N)*DRAT(J) * AN(N)
          BN(N)=BYSN(N)*DRAT(J) * BN(N)
        END DO
        AN(IMH) = BYSN(IMH)*DRAT(J) * AN(IMH)
        CALL FFTI(AN,BN,X(1,J))
      END DO

      RETURN
      END SUBROUTINE AVRX

      SUBROUTINE FILTER
!@sum  FILTER Performs 8-th order shapiro filter in zonal direction
!@auth Original development team
!@ver  1.0
!@calls SHAP1D
C****
C**** MFILTR=1  SMOOTH P USING SEA LEVEL PRESSURE FILTER
C****        2  SMOOTH T USING TROPOSPHERIC STRATIFICATION OF TEMPER
C****        3  SMOOTH P AND T
C****
      USE CONSTANT, only : bygrav,kapa,sha,mb2kg
      USE MODEL_COM, only : im,jm,lm,ls1,t,p,q,wm,mfiltr,zatmo,ptop
     *     ,byim,sig,itime,psf,pmtop
      USE GEOM, only : areag,dxyp
      USE SOMTQ_COM, only : tmom,qmom
      USE DYNAMICS, only : pk, COS_LIMIT
#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm,trname,ITIME_TR0,trm,trmom
#endif
      USE PBLCOM, only : tsavg
      USE DOMAIN_DECOMP_1D, Only : grid, GET, GLOBALSUM
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) :: X,Y
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *        POLD, PRAT
      REAL*8 PSUMO,PSUMN,PDIF,AKAP,PS,ZS
      REAL*8, EXTERNAL :: SLP
      INTEGER I,J,L,N  !@var I,J,L  loop variables
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO) :: KEJ,PEJ
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      REAL*8 initialTotalEnergy, finalTotalEnergy
      real*8 getTotalEnergy ! external for now

#ifdef MASSFLUX_SLPFIL
      ! SLP filter already performed in mass-flux form
      return
#endif

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)

      IF (MOD(MFILTR,2).NE.1) GO TO 200
C**** Initialise total energy (J/m^2)
      initialTotalEnergy = getTotalEnergy()
C****
C**** SEA LEVEL PRESSURE FILTER ON P
C****
!$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,PS,ZS)
      DO J=J_0S,J_1S
        DO I=1,IM
          POLD(I,J)=P(I,J)      ! Save old pressure
          PS=P(I,J)+PTOP
          ZS=ZATMO(I,J)*BYGRAV
          X(I,J)=SLP(PS,TSAVG(I,J),ZS)
          Y(I,J)=X(I,J)/PS
        END DO
      END DO
!$OMP  END PARALLEL DO
      CALL SHAP1D (8,X)
      call isotropslp(x,COS_LIMIT)
!$OMP  PARALLEL DO PRIVATE(I,J,PSUMO,PSUMN,PDIF)
      DO J=J_0S,J_1S
        PSUMO=0.
        PSUMN=0.
        DO I=1,IM
          PSUMO=PSUMO+P(I,J)
          P(I,J)=X(I,J)/Y(I,J)-PTOP
C**** reduce large variations (mainly due to topography)
          P(I,J)=MIN(MAX(P(I,J),0.99d0*POLD(I,J)),1.01d0*POLD(I,J))
          PSUMN=PSUMN+P(I,J)
        END DO
        PDIF=(PSUMN-PSUMO)*BYIM
        DO I=1,IM
          P(I,J)=P(I,J)-PDIF
        END DO
      END DO
!$OMP  END PARALLEL DO
C**** Scale mixing ratios (incl moments) to conserve mass/heat
!$OMP  PARALLEL DO PRIVATE(I,J)
      DO J=J_0S,J_1S
        DO I=1,IM
          PRAT(I,J)=POLD(I,J)/P(I,J)
        END DO
      END DO
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO PRIVATE (I,J,L)
      DO L=1,LS1-1
      DO J=J_0S,J_1S
      DO I=1,IM
c adjust pot. temp. to maintain unchanged absolute temp.
        T(I,J,L)= T(I,J,L)*
     &       ((POLD(I,J)*SIG(L)+PTOP)/(P(I,J)*SIG(L)+PTOP))**KAPA
        Q(I,J,L)= Q(I,J,L)*PRAT(I,J)
        WM(I,J,L)=WM(I,J,L)*PRAT(I,J)
        QMOM(:,I,J,L)=QMOM(:,I,J,L)*PRAT(I,J)
      END DO
      END DO
      END DO
!$OMP  END PARALLEL DO
#ifdef TRACERS_ON
C**** In general, only an air tracer is affected by the filter
C**** This fix conserves tracer concentration, BUT NOT MASS!
C****   It is not wanted for most tracers.
C**** Thus, this code slows model and should be removed if not wanted
C**** Instead of if(trname...) could use n=n_air and avoid the IF-test
C**** But if n_air=0 this will cause problems...
      do n=1,ntm
      if (trname(n).ne.'Air' .and. trname(n).ne.'CO2n') cycle
!     if (itime.lt.itime_tr0(n)) cycle   !probably not needed
!$OMP  PARALLEL DO PRIVATE (I,J,L)
      DO L=1,LS1-1
        DO J=J_0S,J_1S
          DO I=1,IM
             trm(I,J,L,n)=  trm(I,J,L,n)/PRAT(I,J)
             trmom(:,I,J,L,n)=trmom(:,I,J,L,n)/PRAT(I,J)
      end do; end do; end do
!$OMP  END PARALLEL DO
      end do
#endif
      CALL CALC_AMPK(LS1-1)

C**** This fix adjusts thermal energy to conserve total energy TE=KE+PE
      finalTotalEnergy = getTotalEnergy()
      call addEnergyAsDiffuseHeat(finalTotalEnergy - initialTotalEnergy)

  200 IF (MFILTR.LT.2) RETURN
C****
C**** TEMPERATURE STRATIFICATION FILTER ON T
C****
      AKAP=KAPA-.205d0    ! what is this number?
!$OMP  PARALLEL DO PRIVATE (J,L,X,Y)
      DO L=1,LM
        IF(L.LT.LS1) THEN
          DO J=J_0S,J_1S
            Y(:,J)=(SIG(L)*P(:,J)+PTOP)**AKAP
            X(:,J)=T(:,J,L)*Y(:,J)
          END DO
          CALL SHAP1D (8,X)
          DO J=J_0S,J_1S
            T(:,J,L)=X(:,J)/Y(:,J)
          END DO
        ELSE
          DO J=J_0S,J_1S
            X(:,J)=T(:,J,L)
          END DO
          CALL SHAP1D (8,X)
          DO J=J_0S,J_1S
            T(:,J,L)=X(:,J)
          END DO
        END IF
      END DO
!$OMP  END PARALLEL DO
C
      RETURN
      END SUBROUTINE FILTER

      SUBROUTINE PFILTER(P,PUSUM_filter)
!@sum
!@+    SEA LEVEL PRESSURE FILTER ON P in zonal direction
!@+    in mass-flux form
!@+
!@auth Original development team
!@ver  1.0
!@calls SHAP1D
      USE CONSTANT, only : bygrav
      USE MODEL_COM, only : im,jm,zatmo,ptop,byim,dtsrc
      USE GEOM, only : dxyp
      USE DYNAMICS, only : COS_LIMIT
      USE PBLCOM, only : tsavg
      USE DOMAIN_DECOMP_1D, Only : grid, get, halo_update
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     P,PUSUM_filter
c
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     X,Y,POLD,PNEW
      REAL*8 PSUMO,PSUMN,PDIF,PS,ZS
      REAL*8, EXTERNAL :: SLP
      INTEGER I,J,L,N  !@var I,J,L  loop variables

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)

      DO J=J_0S,J_1S
        DO I=1,IM
          POLD(I,J)=P(I,J)      ! Save old pressure
          PS=P(I,J)+PTOP
          ZS=ZATMO(I,J)*BYGRAV
          X(I,J)=SLP(PS,TSAVG(I,J),ZS)
          Y(I,J)=X(I,J)/PS
        END DO
      END DO
      CALL SHAP1D (8,X)
      call isotropslp(x,COS_LIMIT)
      DO J=J_0S,J_1S
        PSUMO=0.
        PSUMN=0.
        DO I=1,IM
          PSUMO=PSUMO+P(I,J)
          PNEW(I,J)=X(I,J)/Y(I,J)-PTOP
C**** reduce large variations (mainly due to topography)
          PNEW(I,J)=
     &         MIN(MAX(PNEW(I,J),0.99d0*POLD(I,J)),1.01d0*POLD(I,J))
          PSUMN=PSUMN+PNEW(I,J)
        END DO
        PDIF=(PSUMN-PSUMO)*BYIM
        DO I=1,IM
          PNEW(I,J)=PNEW(I,J)-PDIF
        END DO
! dp(i) = x(i,j)-y(i,j) = (pusum(i-1,j) - pusum(i,j))*dt/dxyp(j)
        pusum_filter(1,j) = 0.
        do i=2,im
          pusum_filter(i,j) = pusum_filter(i-1,j)
     &         - (pnew(i,j)-pold(i,j))
        enddo
        ! Subtract zonal mean
        pusum_filter(:,j) = pusum_filter(:,j)
     &       -sum(pusum_filter(:,j))*byim
        ! Convert to units of pu: mb*m2/s
        ! Using physics timestep, since original FILTER routine
        ! is called once per physics timestep.
        pusum_filter(:,j) = pusum_filter(:,j)*dxyp(j)/DTsrc
      enddo

      call halo_update(grid,pusum_filter)

      RETURN
      END SUBROUTINE PFILTER

      subroutine fltry2(q3d,strength)
!@sum  fltry2 noise reduction filter for a velocity-type field
!@sum  at secondary latitudes
!@ver  1.0
      use resolution, only : im,jm,lm
      use domain_decomp_1d, only : get,grid,halo_update,north,south
      implicit none
      integer, parameter :: nshap=8
      real*8, parameter :: by4ton=1./(4.**nshap)
      real*8 :: dt
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lm) :: q3d
      real*8 :: strength
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lm) :: yn
      real*8 yvby4ton
      integer i,j,l,jj,n,j_0stg,j_1stg,j_1f
      real*8, dimension(im) :: yjm1,yj
      logical :: have_south_pole,have_north_pole
c      real*8 :: by4ton
c      by4ton=1./(4.**nshap)

      call get(grid, j_strt_stgr = j_0stg, j_stop_stgr = j_1stg,
     &         have_south_pole = have_south_pole,
     &         have_north_pole = have_north_pole)


      yvby4ton = min(strength,1d0)*by4ton*((-1)**(nshap))

      if(have_north_pole) then
        j_1f=jm-1
      else
        j_1f=j_1stg
      endif

      do l=1,lm
        do j=j_0stg,j_1stg
          yn(:,j,l)=q3d(:,j,l)
        enddo
      enddo
      do n=1,nshap
        call halo_update(grid, yn, from=north+south)
        do l=1,lm
          if(have_south_pole) then ! pole-crossing conditions
            yjm1(1:im/2)    = -yn(im/2+1:im,2,l)
            yjm1(im/2+1:im) = -yn(1:im/2,2,l)
          else
            yjm1(:)   = yn(:,j_0stg-1,l)
          endif
          do j=j_0stg,j_1f
            do i=1,im
              yj(i)   = yn(i,j,l)
              yn(i,j,l) = yjm1(i)-yj(i)-yj(i)+yn(i,j+1,l)
              yjm1(i) = yj(i)
            enddo
          enddo
          if(have_north_pole) then ! pole-crossing conditions
            j=jm
            do i=1,im/2
              yj(i)   = yn(i,j,l)
              yn(i,j,l) = yjm1(i)-yj(i)-yj(i)-yn(i+im/2,j,l)
            enddo
            do i=im/2+1,im
              yj(i)   = yn(i,j,l)
              yn(i,j,l) = yjm1(i)-yj(i)-yj(i)-yj(i-im/2)
            enddo
          endif
        enddo                 ! l
      enddo                   ! nshap
      do l=1,lm
        do j=j_0stg,j_1stg
          q3d(:,j,l) = q3d(:,j,l) -yn(:,j,l)*yvby4ton
        enddo
      enddo
      return
      end subroutine fltry2

      SUBROUTINE FLTRUV(U,V,UT,VT)
!@sum  FLTRUV Filters 2 gridpoint noise from the velocity fields
!@auth Original development team
!@ver  1.0
      USE CONSTANT, only : sha
      USE MODEL_COM, only : im,jm,lm,byim,mrch,dt,t,ang_uv
     *  ,DT_XUfilter,DT_XVfilter,DT_YVfilter,DT_YUfilter
     &  ,do_polefix
      USE GEOM, only : dxyn,dxys
      USE DYNAMICS, only : pdsig,pk, COS_LIMIT
c      USE DIAG, only : diagcd
C**********************************************************************
C**** FILTERING IS DONE IN X-DIRECTION WITH A 8TH ORDER SHAPIRO
C**** FILTER. THE EFFECT OF THE FILTER IS THAT OF DISSIPATION AT
C**** THE SMALLEST SCALES.
C**********************************************************************
      USE DOMAIN_DECOMP_1D, only : grid, GET
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, HALO_UPDATE_COLUMN
      USE DOMAIN_DECOMP_1D, only : NORTH, SOUTH
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM),
     *     INTENT(INOUT) :: U,V
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM),
     *     INTENT(IN) :: UT,VT
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *     DUT,DVT,USAVE,VSAVE
      REAL*8 X(IM),YV(max(2*JM,IM)),DP(IM)
      REAL*8 XUby4toN,XVby4toN,YVby4toN,YUby4toN
      REAL*8 :: DT1=0.
      INTEGER I,J,K,L,N,IP1  !@var I,J,L,N  loop variables
      REAL*8 YV2,YVJ,YVJM1,X1,XI,XIM1
      INTEGER, PARAMETER :: NSHAP=8  ! NSHAP MUST BE EVEN
      REAL*8, PARAMETER :: BY16=1./16., by4toN=1./(4.**NSHAP)
      REAL*8 angm,dpt,D2V,D2U
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      INTEGER :: II
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
C****
      USAVE=U ; VSAVE=V
      if (DT_XUfilter.gt.0.) then
        XUby4toN = (DT/DT_XUfilter)*by4toN
      else
        XUby4toN = 0.
      end if
      if (DT_XVfilter.gt.0.) then
        XVby4toN = (DT/DT_XVfilter)*by4toN
      else
        XVby4toN = 0.
      end if
C****
C**** Filtering in east-west direction
C****
!$OMP  PARALLEL DO PRIVATE (I,J,L,N,X,X1,XI,XIM1)
      DO 350 L=1,LM
C**** Filter U component of velocity
      DO 240 J=J_0STG,J_1STG
      DO 210 I=1,IM
  210 X(I) = U(I,J,L)
      DO 230 N=1,NSHAP
      X1   = X(1)
      XIM1 = X(IM)
      DO 220 I=1,IM-1
      XI   = X(I)
      X(I) = XIM1-XI-XI+X(I+1)
  220 XIM1 = XI
  230 X(IM)= XIM1-X(IM)-X(IM)+X1
      DO 240 I=1,IM
  240 U(I,J,L) = U(I,J,L) - X(I)*XUby4toN
C**** Filter V component of velocity
      DO 340 J=J_0STG,J_1STG
      DO 310 I=1,IM
  310 X(I) = V(I,J,L)
      DO 330 N=1,NSHAP
      X1   = X(1)
      XIM1 = X(IM)
      DO 320 I=1,IM-1
      XI   = X(I)
      X(I) = XIM1-XI-XI+X(I+1)
  320 XIM1 = XI
  330 X(IM)= XIM1-X(IM)-X(IM)+X1
      DO 340 I=1,IM
  340 V(I,J,L) = V(I,J,L) - X(I)*XVby4toN
  350 CONTINUE
!$OMP  END PARALLEL DO

c This routine is now called by subroutine DYNAM
c      if(do_polefix.eq.1) then
c         call isotropuv(u,v,COS_LIMIT)
c      endif

C**** Conserve angular momentum along latitudes
c***  The following halo is not needed because PDSIG halo is up to date
c***      CALL HALO_UPDATE_COLUMN(grid, PDSIG, FROM=SOUTH)
!$OMP  PARALLEL DO PRIVATE (I,IP1,J,L,DP,ANGM,DPT)
      DO L=1,LM
        DO J=J_0STG,J_1STG
          ANGM=0.
          DPT=0.
          I=IM
          DO IP1=1,IM
            DP(I)=0.5*((PDSIG(L,IP1,J-1)+PDSIG(L,I,J-1))*DXYN(J-1)
     *           +(PDSIG(L,IP1,J  )+PDSIG(L,I,J  ))*DXYS(J  ))
            ANGM=ANGM-DP(I)*(U(I,J,L)-USAVE(I,J,L))
            DPT=DPT+DP(I)
            I=IP1
          END DO
          DO I=1,IM
            if (ang_uv.eq.1) U(I,J,L)=U(I,J,L)+ANGM/DPT
            DUT(I,J,L)=(U(I,J,L)-USAVE(I,J,L))*DP(I)
            DVT(I,J,L)=(V(I,J,L)-VSAVE(I,J,L))*DP(I)
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO

C**** Call diagnostics only for even time step
      IF (MRCH.eq.2) THEN
!        CALL DIAGCD(grid,5,UT,VT,DUT,DVT,DT1)
      END IF

      RETURN
      END SUBROUTINE FLTRUV

      subroutine isotropslp(slp,coscut)
      use MODEL_COM, only : im,jm,dt
      USE DOMAIN_DECOMP_1D, Only : GET,grid
      use GEOM, only : cosp,dxp
      implicit none
      real*8, parameter :: k=1d3
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: slp
      real*8 :: coscut,fac
      integer :: ipole,j,jcut,jinc,jp
      integer :: j_0s, j_1s, hemi

      CALL GET(grid, J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)

      Do j = j_0s, j_1s
        If (far_from_pole(j, cosp, coscut)) Cycle
        fac = k*dt/(dxp(j)*dxp(j))
        Call shap1(slp(1,j),im,fac)
      enddo

      return
      end subroutine isotropslp

      subroutine isotropuv(u,v,coscut)
!@sum  isotropuv isotropizes the velocity field in the near-polar row(s)
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : im,imh,jm,lm,dt
      USE DOMAIN_DECOMP_1D, Only : GET, grid
      USE GEOM, only : cosv,dxv,cosi=>cosiv,sini=>siniv
      implicit none
      real*8, parameter :: klo=1d3,khi=1d7
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *  U, V
      real*8 :: coscut,fac,k
      real*8, dimension(im) :: ua,va
      real*8, dimension(0:imh) :: an,bn
      integer :: i,j,l,hemi,jp,jinc,jcut,ipole
      Integer :: J_0STG, J_1STG

      CALL GET(grid, J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG)

      do l=1,lm
        do j = J_0STG, J_1STG
          hemi = Hemisphere(j)
          If (far_from_pole(j, cosv, coscut)) Cycle

c compute xy velocities
          do i=1,im
            ua(i) = cosi(i)*u(i,j,l)-hemi*sini(i)*v(i,j,l)
            va(i) = cosi(i)*v(i,j,l)+hemi*sini(i)*u(i,j,l)
          enddo
c filter the xy velocities

          k = maxval(abs(u(:,j,l)))*2.*dt/dxv(j)
          if(k.lt.0.5) then
            k = klo
          else if(k.gt.1d0) then
            k = khi
          else
            k = klo + 2d0*(k-0.5)*(khi-klo)
          endif
          fac = k*dt/(dxv(j)*dxv(j))
          call shap1(ua,im,fac)
          call shap1(va,im,fac)
          if(at_pole(j)) then   ! really strong filtering right at the pole
            call fft(ua,an,bn)
            an(2:imh) = 0.
            bn(2:imh) = 0.
            call ffti(an,bn,ua)
            call fft(va,an,bn)
            an(2:imh) = 0.
            bn(2:imh) = 0.
            call ffti(an,bn,va)
          endif
c convert xy velocities back to polar coordinates
          do i=1,im
            u(i,j,l) = cosi(i)*ua(i)+hemi*sini(i)*va(i)
            v(i,j,l) = cosi(i)*va(i)-hemi*sini(i)*ua(i)
          enddo
        enddo                   ! j
      enddo                     ! l

      return
      end subroutine isotropuv

      Integer Function Hemisphere(j)
      Use GEOM, only: FJEQ
      Integer :: j

      If (J < FJEQ) Then
        hemisphere = -1
      Else
        hemisphere = +1
      End If
      End Function Hemisphere

      ! Detect whether at the pole on staggered grid
      Logical Function at_pole(j)
      Use model_com, only : jm
      Integer :: j
      If (j == jm .or. j == 2) Then
        at_pole = .true.
      else
        at_pole = .false.
      end if
      End Function at_pole

      Logical Function far_from_pole(j, cosj, coscut)
      Use MODEL_COM, only: JM
        Integer :: j
        Real*8 :: cosj(JM), coscut

        far_from_pole = (cosj(j) >= coscut)

      End Function far_from_pole

      subroutine shap1(x,im,fac)
      implicit none
      integer :: im
      real*8, dimension(im) :: x
      real*8 :: fac,facby4,x1,xim1,xi
      integer :: i,n,nn
      n = int(fac) + 1
      facby4 = fac*.25d0/n
      do nn=1,n
      x1 = x(1)
      xim1 = x(im)
      do i=1,im-1
         xi = x(i)
         x(i) = x(i) + facby4*(xim1-xi-xi+x(i+1))
         xim1 = xi
      enddo
      i = im
      x(i) = x(i) + facby4*(xim1-x(i)-x(i)+x1)
      enddo
      return
      end subroutine shap1

      SUBROUTINE SHAP1D (NORDER,X)
!@sum  SHAP1D Smoothes in zonal direction use n-th order shapiro filter
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only :im,jm
      USE DOMAIN_DECOMP_1D, Only : grid, GET
      IMPLICIT NONE
!@var NORDER order of shapiro filter (must be even)
      INTEGER, INTENT(IN) :: NORDER
      REAL*8, INTENT(INOUT),
     *        DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) :: X
      REAL*8, DIMENSION(IM)::XS
      REAL*8 by4toN,XS1,XSIM1,XSI
      INTEGER I,J,N   !@var I,J,N  loop variables
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)

      by4toN=1./4.**NORDER
      DO J=J_0S,J_1S
        XS(:) = X(:,J)
        DO N=1,NORDER
          XS1=XS(1)
          XSIM1=XS(IM)
          DO I=1,IM-1
            XSI=XS(I)
            XS(I)=XSIM1-XSI-XSI+XS(I+1)
            XSIM1=XSI
          END DO
          XS(IM)=XSIM1-XS(IM)-XS(IM)+XS1
        END DO
        X(:,J)=X(:,J)-XS(:)*by4toN
      END DO
      RETURN
      END SUBROUTINE SHAP1D

      SUBROUTINE SDRAG(DT1)
!@sum  SDRAG puts a drag on the winds in the top layers of atmosphere
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,sha
      USE MODEL_COM, only : im,jm,lm,ls1,u,v,t,q,x_sdrag,csdragl,lsdrag
     *     ,lpsdrag,ang_sdrag,itime,Wc_Jdrag,wmax,vsdragl,byim
      USE GEOM, only : cosv,imaxj,kmaxj,idij,idjj,rapj,dxyv,dxyn,dxys
     *     ,rapvs,rapvn
      !USE DIAG_COM, only : ajl=>ajl_loc,jl_dudtsdrg
      USE DYNAMICS, only : pk,pdsig,pedn
c      USE DIAG, only : diagcd
      USE DOMAIN_DECOMP_1D, only : grid, GET
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, HALO_UPDATE_COLUMN
      USE DOMAIN_DECOMP_1D, only : NORTH, SOUTH
      IMPLICIT NONE

!@var DT1 time step (s)
      REAL*8, INTENT(IN) :: DT1
!@var L(P)SDRAG lowest level at which SDRAG_lin is applied (near poles)
C**** SDRAG_const is applied above PTOP (150 mb) and below the SDRAG_lin
C**** regime (but not above P_CSDRAG)
      REAL*8 WL,TL,RHO,CDN,X,DP,DPL(LM),du,dps
!@var DUT,DVT change in momentum (mb m^3/s)
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *        DUT,DVT
      INTEGER I,J,L,IP1,K,Lmax
      logical cd_lin
!@var ang_mom is the sum of angular momentun at layers LS1 to LM
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)    ::
     *        ang_mom, sum_airm
!@var wmaxp =.75*wmax,the imposed limit for stratospheric winds (m/s)
      real*8 wmaxp,wmaxj,xjud
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      ang_mom=0. ;  sum_airm=0. ; dut=0.
C*
      DUT=0. ; DVT=0.
      wmaxp = wmax*3.d0/4.d0
c***  The following halo is not needed because PDSIG halo is up to date
c***      CALL HALO_UPDATE_COLUMN(grid, PDSIG, FROM=SOUTH)
      DO L=LS1,LM
      DO J=J_0STG, J_1STG
      cd_lin=.false.
      IF( L.ge.LSDRAG .or.
     *   (L.ge.LPSDRAG.and.COSV(J).LE..15) ) cd_lin=.true.
      wmaxj=wmax
      if(COSV(J).LE..15) wmaxj=wmaxp
      I=IM
      DO IP1=1,IM
        TL=T(I,J,L)*PK(L,I,J)   ! not quite correct - should be on UV grid
C**** check T to make sure it stayed within physical bounds
        if (TL.lt.100..or.TL.gt.373.) then
          write(99,'(a,i8,3i4,a,3f10.2)')
     *    ' SDRAG:',itime,i,j,l,'  T,U,V=',TL,U(I,J,L),V(I,J,L)
          call stop_model('Stopped in ATMDYN::SDRAG',11)
        end if
        RHO=PEDN(L+1,I,J)/(RGAS*TL)   ! not quite correct - should be on UV grid
        WL=SQRT(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
        xjud=1.
        if(Wc_JDRAG.gt.0.) xjud=(Wc_JDRAG/(Wc_JDRAG+min(WL,wmaxj)))**2
C**** WL is restricted to Wmax by adjusting X, if necessary;
C**** the following is equivalent to first reducing (U,V), if necessary,
C**** then finding the drag and applying it to the reduced winds
                    CDN=CSDRAGl(l)*xjud
        IF (cd_lin) CDN=(X_SDRAG(1)+X_SDRAG(2)*min(WL,wmaxj))*xjud
        DPS= (PDSIG(L,IP1,J-1)+PDSIG(L,I,J-1))*RAPVN(J-1)+
     *       (PDSIG(L,IP1,J  )+PDSIG(L,I,J  ))*RAPVS(J)
        X=DT1*RHO*CDN*min(WL,wmaxj)*GRAV*VSDRAGL(L)/DPS
        if (wl.gt.wmaxj) X = 1. - (1.-X)*wmaxj/wl
C**** adjust diags for possible difference between DT1 and DTSRC
c        call inc_ajl(i,j,l,JL_DUDTSDRG,-U(I,J,L)*X) ! for a-grid only
        !ajl(j,l,jl_dudtsdrg) = ajl(j,l,jl_dudtsdrg) -u(i,j,l)*x*byim
        DP=DPS*DXYV(J)
        ang_mom(i,j) = ang_mom(i,j)+U(I,J,L)*X*DP
        DUT(I,J,L)=-X*U(I,J,L)*DP
        DVT(I,J,L)=-X*V(I,J,L)*DP
        U(I,J,L)=U(I,J,L)*(1.-X)
        V(I,J,L)=V(I,J,L)*(1.-X)
        I=IP1
      END DO
      END DO
      END DO

C*
C***  Add the lost angular momentum uniformly back in if ang_sdrag>0
C***  only below 150mb if ang_sdrag=1, into whole column if ang_sdrag>1
C*
      if (ang_sdrag.gt.0) then
        lmax=ls1-1
        if (ang_sdrag.gt.1) lmax=lm
        do j = J_0STG,J_1STG
        I=IM
        do ip1 = 1,im
          do l = 1,lmax
            DPL(L)=0.5*((PDSIG(L,IP1,J-1)+PDSIG(L,I,J-1))*DXYN(J-1)
     *        +(PDSIG(L,IP1,J  )+PDSIG(L,I,J  ))*DXYS(J  ))
            sum_airm(i,j) = sum_airm(i,j)+DPL(L)
          end do
C*
          do l = 1,lmax
            du = ang_mom(i,j)/sum_airm(i,j)
            DUT(I,J,L) = DUT(I,J,L) + du*dpl(l)
c            call inc_ajl(i,j,l,JL_DUDTSDRG,du) ! for a-grid only
            !ajl(j,l,jl_dudtsdrg) = ajl(j,l,jl_dudtsdrg) +du*byim
            U(I,J,L)=U(I,J,L) + du
          end do
          I=IP1
        end do
        end do
      end if

C**** conservation diagnostic
C**** (technically we should use U,V from before but this is ok)
!      CALL DIAGCD (grid,4,U,V,DUT,DVT,DT1)

      RETURN
      END SUBROUTINE SDRAG

      end module ATMDYN

      subroutine add_am_as_solidbody_rotation(u,dam)
      use constant, only : radius,mb2kg
      use model_com, only : im,jm,lm,fim,p,pstrat
      use geom, only : cosv,dxyn,dxys
      use domain_decomp_1d, only : get, south, halo_update, grid
      use domain_decomp_1d, only : globalsum
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lm) :: u
      real*8 :: dam
      integer :: j,l
      real*8 :: u0,xintsum
      real*8, dimension(grid%j_strt_halo:grid%j_stop_halo) ::
     &     psumj,xintj
      integer :: j_0stg, j_1stg, j_0, j_1
      logical :: have_south_pole, have_north_pole

      call get(grid, j_strt=j_0, j_stop=j_1,
     &               j_strt_stgr=j_0stg, j_stop_stgr=j_1stg,
     &               have_south_pole=have_south_pole,
     &               have_north_pole=have_north_pole)

c      call halo_update(grid, p, from=south) ! already done
      do j=j_0stg-1,j_1
        psumj(j) = sum(p(:,j))+fim*pstrat
      enddo
      do j=j_0stg,j_1stg
        xintj(j) = cosv(j)*cosv(j)*
     &       (psumj(j-1)*dxyn(j-1)+psumj(j)*dxys(j))
      enddo
      if(have_south_pole) xintj(1)=0.
      call globalsum(grid,xintj,xintsum,all=.true.)
      u0 = dam/(radius*mb2kg*xintsum)
      do l=1,lm
      do j=j_0stg,j_1stg
        u(:,j,l) = u(:,j,l) + u0*cosv(j)
      enddo
      enddo
      return
      end subroutine add_am_as_solidbody_rotation

      SUBROUTINE conserv_AMB_ext(U,AM)
      USE CONSTANT, only : omega,radius,mb2kg
      USE MODEL_COM, only : im,jm,lm,fim,ls1,dsig,p,psfmpt,pstrat
      USE GEOM, only : cosv,dxyn,dxys,dxyv,byaxyp
      USE DOMAIN_DECOMP_1D, only : GET, SOUTH, HALO_UPDATE, GRID
      USE DOMAIN_DECOMP_1D, only : CHECKSUM
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: U
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: AM
      INTEGER :: I,IP1,J,L
      REAL*8 :: PSJ,PSIJ,UE,UEDMS,FACJ

      INTEGER :: J_0S, J_1S, J_0STG, J_1STG, J_0, J_1, I_0, I_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     *               I_STRT=I_0, I_STOP=I_1,
     *               J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)


C****
C**** ANGULAR MOMENTUM ON B GRID
C****
      CALL HALO_UPDATE(grid, P, FROM=SOUTH)

      DO J=J_0STG,J_1STG
      PSJ=(2.*PSFMPT*DXYV(J))
      UE=RADIUS*OMEGA*COSV(J)
      UEDMS=2.*UE*PSTRAT*DXYV(J)
      FACJ=.5*COSV(J)*RADIUS*mb2kg
      I=IM
      DO IP1=1,IM
        PSIJ=(P(I,J-1)+P(IP1,J-1))*DXYN(J-1)+(P(I,J)+P(IP1,J))*DXYS(J)
        AM(I,J)=0.
        DO L=1,LS1-1
          AM(I,J)=AM(I,J)+U(I,J,L)*DSIG(L)
        END DO
        AM(I,J)=AM(I,J)*PSIJ
        DO L=LS1,LM
          AM(I,J)=AM(I,J)+U(I,J,L)*PSJ*DSIG(L)
        END DO
        AM(I,J)=(UEDMS+UE*PSIJ+AM(I,J))*FACJ
        I=IP1
      END DO
      END DO

      RETURN
C****
      END SUBROUTINE conserv_AMB_ext

      SUBROUTINE conserv_AM(AM)
!@sum  conserv_AM calculates A-grid column-sum atmospheric angular momentum,
!@sum  per unit area
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,u
      USE GEOM, only : byaxyp
      USE DOMAIN_DECOMP_1D, only : GET, GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: AM
      INTEGER :: I,J
      INTEGER :: J_0, J_1, I_0, I_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     *               I_STRT=I_0, I_STOP=I_1)

C****
C**** ANGULAR MOMENTUM ON B GRID
C****
      call conserv_AMB_ext(U,AM)

c move to A grid
      call regrid_btoa_ext(am)

c scale by area
      DO J=J_0,J_1
        DO I=I_0,I_1
          am(I,J)=am(I,J)*BYAXYP(I,J)
        END DO
      END DO

      RETURN
C****
      END SUBROUTINE conserv_AM

      SUBROUTINE conserv_KE(RKE)
!@sum  conserv_KE calculates A-grid column-sum atmospheric kinetic energy,
!@sum  (J/m2)
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg
      USE MODEL_COM, only : im,jm,lm,fim,dsig,ls1,p,u,v,psfmpt
      USE GEOM, only : dxyn,dxys,dxyv,byaxyp
      USE DOMAIN_DECOMP_1D, only : GET, CHECKSUM, HALO_UPDATE, GRID
      USE DOMAIN_DECOMP_1D, only : SOUTH
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: RKE
      INTEGER :: I,IP1,J,L
      INTEGER :: J_0STG,J_1STG, J_0, J_1, I_0, I_1
      REAL*8 :: PSJ,PSIJ

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     *               I_STRT=I_0, I_STOP=I_1,
     *               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG)

C****
C**** KINETIC ENERGY ON B GRID
C****

      CALL HALO_UPDATE(grid, P, FROM=SOUTH)
      DO J=J_0STG,J_1STG
      PSJ=(2.*PSFMPT*DXYV(J))
      I=IM
      DO IP1=1,IM
        PSIJ=(P(I,J-1)+P(IP1,J-1))*DXYN(J-1)+(P(I,J)+P(IP1,J))*DXYS(J)
        RKE(I,J)=0.
        DO L=1,LS1-1
          RKE(I,J)=RKE(I,J)+
     &         (U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))*DSIG(L)
        END DO
        RKE(I,J)=RKE(I,J)*PSIJ
        DO L=LS1,LM
          RKE(I,J)=RKE(I,J)+
     &         (U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))*DSIG(L)*PSJ
        END DO
        RKE(I,J)=0.25*RKE(I,J)*mb2kg
        I=IP1
      END DO
      END DO

c move to A grid
      call regrid_btoa_ext(rke)

c scale by area
      DO J=J_0,J_1
        DO I=I_0,I_1
          rke(I,J)=rke(I,J)*BYAXYP(I,J)
        END DO
      END DO

      RETURN
C****
      END SUBROUTINE conserv_KE

      SUBROUTINE calc_kea_3d(kea)
!@sum  calc_kea_3d calculates square of wind speed on the A grid
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,byim,u,v
c      USE GEOM, only : ravps,ravpn
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, GRID
      USE DOMAIN_DECOMP_1D, only : NORTH
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: KEA
      INTEGER :: I,J,L
      DO L=1,LM
      DO J=GRID%J_STRT_STGR,GRID%J_STOP_STGR
      DO I=1,IM
        KEA(I,J,L)=.5*(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
      ENDDO
      ENDDO
      ENDDO
      call regrid_btoa_3d(kea)
      RETURN

      END SUBROUTINE calc_kea_3d

      subroutine recalc_agrid_uv
!@sum recalc_agrid_uv Computes u_a,v_a from u and v
!@var u x-component at secondary grids (B_grid)
!@var v y-component at secondary grids (B_grid)
!@var u_a x-component at primary grids (A_grid)
!@var v_a y-component at primary grids (A_grid)
!@auth Ye Cheng
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,u,v
      USE DYNAMICS, only : ua=>ualij,va=>valij
      USE DOMAIN_DECOMP_1D, only : grid,get,NORTH, HALO_UPDATE_COLUMN
      USE DOMAIN_DECOMP_1D, only : halo_update
      USE GEOM, only : imaxj,idij,idjj,kmaxj,rapj,cosiv,siniv
      implicit none
      real*8, dimension(im) :: ra
      integer, dimension(im) :: idj
      real*8 :: HEMI,u_t,v_t,rak,ck,sk,uk,vk
      integer :: i,j,l,k,idik,idjk,kmax

      integer :: J_0S,J_1S
      logical :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      call get(grid, J_STRT_SKP=J_0S,   J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE    )
!     polar boxes

C**** Update halos of U and V
      CALL HALO_UPDATE(grid,u, from=NORTH)
      CALL HALO_UPDATE(grid,v, from=NORTH)


      if (HAVE_SOUTH_POLE) then
        J=1
        KMAX=KMAXJ(J)
        HEMI=-1.
!$OMP  PARALLEL DO PRIVATE (I,L,u_t,v_t,K,IDIK,IDJK,RAK,ck,sk,uk,vk)
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJJ(K,J)
              RAK=RAPJ(K,J)
              ck=cosiv(k)
              sk=siniv(k)
              uk=u(idik,idjk,L)
              vk=v(idik,idjk,L)
              u_t=u_t+rak*(uk*ck-hemi*vk*sk)
              v_t=v_t+rak*(vk*ck+hemi*uk*sk)
            END DO
            ua(l,i,j)=u_t
            va(l,i,j)=v_t
          END DO
        END DO
!$OMP  END PARALLEL DO
      end if              !south pole
!
      if (HAVE_NORTH_POLE) then
        J=JM
        KMAX=KMAXJ(J)
        HEMI=1.
!$OMP  PARALLEL DO PRIVATE (I,L,u_t,v_t,K,IDIK,IDJK,RAK,ck,sk,uk,vk)
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJJ(K,J)
              RAK=RAPJ(K,J)
              ck=cosiv(k)
              sk=siniv(k)
              uk=u(idik,idjk,L)
              vk=v(idik,idjk,L)
              u_t=u_t+rak*(uk*ck-hemi*vk*sk)
              v_t=v_t+rak*(vk*ck+hemi*uk*sk)
            END DO
            ua(l,i,j)=u_t
            va(l,i,j)=v_t
          END DO
        END DO
!$OMP  END PARALLEL DO
      end if                !north pole

!     non polar boxes
C**** Update halos of u and v. (Needed bcs. IDJJ(3:4,J_1S)=J_1S+1)
C     ---> done by calling routine...

c      CALL HALO_UPDATE(grid, u, FROM=NORTH)
c      CALL HALO_UPDATE(grid, v, FROM=NORTH)
!$OMP  PARALLEL DO PRIVATE (J,I,L,u_t,v_t,K,KMAX,IDJ,RA,IDIK,IDJK,RAK)
!$OMP*    SCHEDULE(DYNAMIC,2)
      DO J=J_0S,J_1S
        KMAX=KMAXJ(J)
        DO K=1,KMAX
          IDJ(K)=IDJJ(K,J)
          RA(K)=RAPJ(K,J)
        END DO
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJ(K)
              RAK=RA(K)
              u_t=u_t+u(IDIK,IDJK,L)*RAK
              v_t=v_t+v(IDIK,IDJK,L)*RAK
            END DO
            ua(l,i,j)=u_t
            va(l,i,j)=v_t
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO
C****
      return
      end subroutine recalc_agrid_uv

      subroutine regrid_atov_1d(u_a,v_a,uv1d)
      USE MODEL_COM, only : im,jm
      USE DOMAIN_DECOMP_1D, only : grid,halo_update,SOUTH
      USE GEOM, only : rapvs,rapvn,cosiv,siniv
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo)  ::
     &          u_a,v_a
      real*8, dimension(2*im*(1+grid%j_stop_stgr-grid%j_strt_stgr)),
     &        intent(out) :: uv1d
      real*8 :: hemi
      integer :: i,ip1,j,n
      real*8, dimension(im) :: usouth,vsouth,unorth,vnorth
      integer :: j_0stg, j_1stg
      j_0stg = grid%j_strt_stgr
      j_1stg = grid%j_stop_stgr
      CALL HALO_UPDATE(grid,U_A,from=SOUTH)
      CALL HALO_UPDATE(grid,V_A,from=SOUTH)
      j=j_0stg-1
      if(j.eq.1) then
        hemi = -1.
        usouth(:)=2.*(u_a(1,j)*cosiv(:)+v_a(1,j)*siniv(:)*hemi)
        vsouth(:)=2.*(v_a(1,j)*cosiv(:)-u_a(1,j)*siniv(:)*hemi)
      else
        i=im
        do ip1=1,im
          usouth(i)=(u_a(i,j)+u_a(ip1,j))
          vsouth(i)=(v_a(i,j)+v_a(ip1,j))
          i=ip1
        enddo
      endif
      n = 0
      do j=j_0stg,j_1stg
        if(j.lt.jm) then
          i=im
          do ip1=1,im
            unorth(i)=(u_a(i,j)+u_a(ip1,j))
            vnorth(i)=(v_a(i,j)+v_a(ip1,j))
            i=ip1
          enddo
        else
          hemi = +1.
          unorth(:)=2.*(u_a(1,j)*cosiv(:)+v_a(1,j)*siniv(:)*hemi)
          vnorth(:)=2.*(v_a(1,j)*cosiv(:)-u_a(1,j)*siniv(:)*hemi)
        endif
        do i=1,im
          n = n + 1
          uv1d(n) = rapvn(j-1)*usouth(i)+rapvs(j)*unorth(i)
          n = n + 1
          uv1d(n) = rapvn(j-1)*vsouth(i)+rapvs(j)*vnorth(i)
          usouth(i) = unorth(i)
          vsouth(i) = vnorth(i)
        enddo
      enddo
      return
      end subroutine regrid_atov_1d

      subroutine get_nuv(nuv)
      use model_com, only : im
      USE DOMAIN_DECOMP_1D, only : GRID
      implicit none
      integer :: nuv
      nuv = 2*im*(1+grid%j_stop_stgr-grid%j_strt_stgr)
      return
      end subroutine get_nuv

      subroutine get_vpkey_of_n(n,vpkey)
      implicit none
      integer :: n,vpkey
      vpkey = 1+(n-1)/2
      return
      end subroutine get_vpkey_of_n

      subroutine get_regrid_info_for_n(n,ilist,jlist,wts,nnbr)
      use model_com, only : im
      use geom, only : rapvn,rapvs
      implicit none
      integer :: n
      integer, dimension(4) :: ilist,jlist
      real*8, dimension(4) :: wts
      integer :: nnbr
      integer :: iv,jv,ivp1
      call get_ivjv_of_n(n,iv,jv)
      nnbr = 4
      ivp1 = iv+1 - im*(iv/im)
      ilist(1:4) = (/ iv, ivp1, iv, ivp1 /)
      jlist(1:4) = (/ jv-1, jv-1, jv, jv /)
      wts(1:4) = (/ rapvn(jv-1), rapvn(jv-1), rapvs(jv), rapvs(jv) /)
      return
      end subroutine get_regrid_info_for_n

      subroutine get_uv_of_n(n,uv)
      use model_com, only : im,lm,u,v
      use domain_decomp_1d, only : am_i_root
      implicit none
      integer :: n
      real*8, dimension(lm) :: uv
      integer :: iv,jv
      call get_ivjv_of_n(n,iv,jv)
      if(mod(n,2).eq.1) then
        uv(1:lm) = u(iv,jv,1:lm)
      else
        uv(1:lm) = v(iv,jv,1:lm)
      endif
      return
      end subroutine get_uv_of_n

      subroutine store_uv_of_n(n,uv)
      use model_com, only : im,lm,u,v
      implicit none
      integer :: n
      real*8, dimension(lm) :: uv
      integer :: iv,jv
      call get_ivjv_of_n(n,iv,jv)
      if(mod(n,2).eq.1) then
        u(iv,jv,1:lm) = uv(1:lm)
      else
        v(iv,jv,1:lm) = uv(1:lm)
      endif
      return
      end subroutine store_uv_of_n

      subroutine get_ivjv_of_n(n,iv,jv)
      use model_com, only : im
      USE DOMAIN_DECOMP_1D, only : GRID
      implicit none
      integer :: n
      integer :: iv,jv
      integer :: nv,njm1
      nv = 1+(n-1)/2
      njm1 = (nv-1)/im
      jv = grid%j_strt_stgr + njm1
      iv = nv - njm1*im
      return
      end subroutine get_ivjv_of_n

      subroutine replicate_uv_to_agrid(ur,vr,k,ursp,vrsp,urnp,vrnp)
      USE MODEL_COM, only : im,jm,lm,u,v
      USE DOMAIN_DECOMP_1D, only : GRID,
     &     hasSouthPole, hasNorthPole
      implicit none
      integer :: k
      REAL*8, DIMENSION(k,LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     UR,VR
      real*8, dimension(im,lm) :: ursp,vrsp,urnp,vrnp
      integer :: i,j,l
      integer :: J_0S,J_1S
      if(k.ne.4)
     &     call stop_model('incorrect k in replicate_uv_to_agrid',255)
      J_0S = GRID%J_STRT_SKP
      J_1S = GRID%J_STOP_SKP
      do j=j_0s,j_1s
      do i=2,im
      do l=1,lm
        ur(1,l,i,j) = u(i-1,j  ,l)
        vr(1,l,i,j) = v(i-1,j  ,l)
        ur(2,l,i,j) = u(i  ,j  ,l)
        vr(2,l,i,j) = v(i  ,j  ,l)
        ur(3,l,i,j) = u(i-1,j+1,l)
        vr(3,l,i,j) = v(i-1,j+1,l)
        ur(4,l,i,j) = u(i  ,j+1,l)
        vr(4,l,i,j) = v(i  ,j+1,l)
      enddo ! l
      enddo ! i
      i = 1
      do l=1,lm
        ur(1,l,i,j) = u(im ,j  ,l)
        vr(1,l,i,j) = v(im ,j  ,l)
        ur(2,l,i,j) = u(i  ,j  ,l)
        vr(2,l,i,j) = v(i  ,j  ,l)
        ur(3,l,i,j) = u(im ,j+1,l)
        vr(3,l,i,j) = v(im ,j+1,l)
        ur(4,l,i,j) = u(i  ,j+1,l)
        vr(4,l,i,j) = v(i  ,j+1,l)
      enddo ! l
      enddo ! j
      if(hasSouthPole(grid)) then
        ursp(:,:) = u(:,2,:)
        vrsp(:,:) = v(:,2,:)
      endif
      if(hasNorthPole(grid)) then
        urnp(:,:) = u(:,jm,:)
        vrnp(:,:) = v(:,jm,:)
      endif
      return
      end subroutine replicate_uv_to_agrid

      subroutine avg_replicated_duv_to_vgrid(du,dv,k,
     &     dusp,dvsp,dunp,dvnp)
      USE MODEL_COM, only : im,jm,lm,u,v
      USE DOMAIN_DECOMP_1D, only : GRID, HALO_UPDATE_BLOCK,SOUTH,
     &     hasSouthPole, hasNorthPole
      implicit none
      integer :: k
      REAL*8, DIMENSION(k,LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     DU,DV
      real*8, dimension(im,lm) :: dusp,dvsp,dunp,dvnp
      integer :: i,j,l
      integer :: J_0STG,J_1STG
      if(k.ne.4) call stop_model(
     &     'incorrect k in avg_replicated_duv_to_vgrid',255)
      J_0STG = GRID%J_STRT_STGR
      J_1STG = GRID%J_STOP_STGR
      CALL HALO_UPDATE_BLOCK(GRID,DU,FROM=SOUTH)
      CALL HALO_UPDATE_BLOCK(GRID,DV,FROM=SOUTH)
c
c copy circumpolar data into the appropriate spots in du,dv
c
      if(hasSouthPole(grid)) then
        j=1
        do i=2,im
        do l=1,lm
          du(3,l,i,j) = dusp(i-1,l)
          du(4,l,i,j) = dusp(i  ,l)
          dv(3,l,i,j) = dvsp(i-1,l)
          dv(4,l,i,j) = dvsp(i  ,l)
        enddo
        enddo
        i=1
        do l=1,lm
          du(3,l,i,j) = dusp(im ,l)
          du(4,l,i,j) = dusp(i  ,l)
          dv(3,l,i,j) = dvsp(im ,l)
          dv(4,l,i,j) = dvsp(i  ,l)
        enddo
#ifndef ALT_CLDMIX_UV
c compensate for the factor of 2 in ravj(1).  change ravj(1) later.
        du(3:4,:,:,j) = du(3:4,:,:,j)*.5
        dv(3:4,:,:,j) = dv(3:4,:,:,j)*.5
#endif
      endif
      if(hasNorthPole(grid)) then
        j=jm
        do i=2,im
        do l=1,lm
          du(1,l,i,j) = dunp(i-1,l)
          du(2,l,i,j) = dunp(i  ,l)
          dv(1,l,i,j) = dvnp(i-1,l)
          dv(2,l,i,j) = dvnp(i  ,l)
        enddo
        enddo
        i=1
        do l=1,lm
          du(1,l,i,j) = dunp(im ,l)
          du(2,l,i,j) = dunp(i  ,l)
          dv(1,l,i,j) = dvnp(im ,l)
          dv(2,l,i,j) = dvnp(i  ,l)
        enddo
#ifndef ALT_CLDMIX_UV
c compensate for the factor of 2 in ravj(jm).  change ravj(jm) later.
        du(1:2,:,:,j) = du(1:2,:,:,j)*.5
        dv(1:2,:,:,j) = dv(1:2,:,:,j)*.5
#endif
      endif
c
c now do the averaging
c
      do j=j_0stg,j_1stg
      do i=1,im-1
      do l=1,lm
        u(i,j,l)=u(i,j,l)+
#ifdef ALT_CLDMIX_UV
     &       .25*
#endif
     &       (du(4,l,i,j-1)+du(3,l,i+1,j-1)+du(2,l,i,j)+du(1,l,i+1,j))
        v(i,j,l)=v(i,j,l)+
#ifdef ALT_CLDMIX_UV
     &       .25*
#endif
     &       (dv(4,l,i,j-1)+dv(3,l,i+1,j-1)+dv(2,l,i,j)+dv(1,l,i+1,j))
      enddo ! l
      enddo ! i
      i = im
      do l=1,lm
        u(i,j,l)=u(i,j,l)+
#ifdef ALT_CLDMIX_UV
     &       .25*
#endif
     &       (du(4,l,i,j-1)+du(3,l,1,j-1)+du(2,l,i,j)+du(1,l,1,j))
        v(i,j,l)=v(i,j,l)+
#ifdef ALT_CLDMIX_UV
     &       .25*
#endif
     &       (dv(4,l,i,j-1)+dv(3,l,1,j-1)+dv(2,l,i,j)+dv(1,l,1,j))
      enddo ! l
      enddo ! j
      return
      end subroutine avg_replicated_duv_to_vgrid

      SUBROUTINE regrid_btoa_3d(x)
      USE MODEL_COM, only : im,jm,lm,byim
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, GRID
      USE DOMAIN_DECOMP_1D, only : NORTH,
     &     hasSouthPole, hasNorthPole
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: X
      INTEGER :: I,IM1,J,L
      REAL*8 :: XIM1J,XIJ
      call halo_update(grid,x,from=north)
      DO L=1,LM
      if(hasSouthPole(grid)) then
        x(:,1,l) = sum(x(:,2,l))*byim
      endif
      DO J=GRID%J_STRT_SKP,GRID%J_STOP_SKP
      IM1=IM
      XIM1J = x(im1,j,l)
      DO I=1,IM
        XIJ = x(i,j,l)
        X(I,J,L)=.25*(XIM1J+XIJ+
     &       x(im1,j+1,l)+x(i,j+1,l))
        XIM1J = XIJ
        IM1=I
      ENDDO
      ENDDO
      if(hasNorthPole(grid)) then
        x(:,jm,l) = sum(x(:,jm,l))*byim
      endif
      ENDDO
      RETURN
      END SUBROUTINE regrid_btoa_3d

      subroutine regrid_btoa_ext(x)
c regrids scalar x_bgrid*dxyv -> x_agrid*dxyp
      USE MODEL_COM, only : im,jm,byim
      USE GEOM, only : rapvs,rapvn,dxyp,dxyv
      USE DOMAIN_DECOMP_1D, only : GET, HALO_UPDATE, GRID, NORTH,
     &     hasSouthPole, hasNorthPole
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: X
      INTEGER :: I,IM1,J
      INTEGER :: J_0S,J_1S
      REAL*8 :: XIM1J,XIJ
      CALL GET(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)
      call halo_update(grid,x,from=north)
      if(hasSouthPole(grid)) then
        X(:,1) = SUM(X(:,2))*BYIM*(DXYP(1)/DXYV(2))
      endif
      DO J=J_0S,J_1S
      IM1=IM
      XIM1J = X(IM1,J)
      DO I=1,IM
        XIJ = X(I,J)
c        X(I,J) = .25*(XIM1J+X(I,J)+X(IM1,J+1)+X(I,J+1))
        X(I,J) = (
     &       (XIM1J+X(I,J))*RAPVS(J)
     &      +(X(IM1,J+1)+X(I,J+1))*RAPVN(J) )
        XIM1J = XIJ
        IM1 = I
      ENDDO
      ENDDO
      if(hasNorthPole(grid)) then
        X(:,JM) = SUM(X(:,JM))*BYIM*(DXYP(JM)/DXYV(JM))
      endif
      return
      end subroutine regrid_btoa_ext

!c      module DIAG
!c      contains
!      SUBROUTINE DIAGCD (grid,M,UX,VX,DUT,DVT,DT1)!,PIT)
!!@sum  DIAGCD Keeps track of the conservation properties of angular
!!@+    momentum and kinetic energy inside dynamics routines
!!@auth Gary Russell
!!@ver  1.0
!      USE CONSTANT, only : omega,mb2kg, radius
!      USE MODEL_COM, only : im,jm,lm,fim,byim,mdiag,mdyn
!      USE GEOM, only : cosv, ravpn,ravps,bydxyp
!      !USE DIAG_COM, only : consrv=>consrv_loc
!      USE DYNAMICS, only : PIT
!      USE DOMAIN_DECOMP_1D, only : GET, CHECKSUM, HALO_UPDATE, DIST_GRID
!      USE DOMAIN_DECOMP_1D, only : SOUTH, NORTH
!      USE GETTIME_MOD
!      IMPLICIT NONE
!C****
!C**** THE PARAMETER M INDICATES WHEN DIAGCD IS BEING CALLED
!C**** M=1  AFTER ADVECTION IN DYNAMICS
!C****   2  AFTER CORIOLIS FORCE IN DYNAMICS
!C****   3  AFTER PRESSURE GRADIENT FORCE IN DYNAMICS
!C****   4  AFTER STRATOS DRAG IN DYNAMICS
!C****   5  AFTER FLTRUV IN DYNAMICS
!C****   6  AFTER GRAVITY WAVE DRAG IN DYNAMICS
!C****
!      TYPE (DIST_GRID), INTENT(IN) :: grid
!!@var M index denoting from where DIAGCD is called
!      INTEGER, INTENT(IN) :: M
!!@var DT1 current time step
!      REAL*8, INTENT(IN) :: DT1
!!@var UX,VX current velocities
!      REAL*8, INTENT(IN),
!     &        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
!     &        UX,VX
!!@var DUT,DVT current momentum changes
!      REAL*8, INTENT(IN),
!     &        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
!     &        DUT,DVT
!!@var PIT current pressure tendency
!!      REAL*8, INTENT(IN), OPTIONAL,
!!     &        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: PIT
!      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: PI
!     &     ,DAMB,DKEB
!      INTEGER :: I,J,L,N,IP1
!      LOGICAL dopit
!      REAL*8 :: DUTI,DUTIL,RKEI,RKEIL,begin
!      INTEGER, DIMENSION(6) ::
!     *     NAMOFM=(/2,3,4,5,6,7/), NKEOFM=(/14,15,16,17,18,19/)
!
!      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H
!      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
!
!      CALL GETTIME(BEGIN)
!
!      CALL GET(grid, J_STRT=J_0,         J_STOP=J_1,
!     &               J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
!     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
!     &               J_STRT_HALO=J_0H,
!     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
!     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
!
!C****
!C**** PRESSURE TENDENCY FOR CHANGE BY ADVECTION
!C****
!      IF (M.eq.1) THEN
!        dopit=.true.
!        IF(HAVE_SOUTH_POLE) PI(1)=FIM*PIT(1,1)
!        IF(HAVE_NORTH_POLE) PI(JM)=FIM*PIT(1,JM)
!        DO J=J_0S,J_1S
!          PI(J)=SUM(PIT(:,J))
!        END DO
!      ELSE
!        PI=0.
!        dopit=.false.
!      END IF
!C****
!C**** CHANGE OF ANGULAR MOMENTUM AND KINETIC ENERGY BY VARIOUS
!C**** PROCESSES IN DYNAMICS
!C****
!C****
!
!      CALL HALO_UPDATE(grid, PI, FROM=SOUTH)
!
!!$OMP PARALLEL DO PRIVATE (J,L,I,DUTIL,RKEIL,DUTI,RKEI,N)
!      DO J=J_0STG,J_1STG
!        DUTIL=0.
!        RKEIL=0.
!        DO L=1,LM
!          DUTI=0.
!          RKEI=0.
!          DO I=1,IM
!            DUTI=DUTI+DUT(I,J,L)
!            RKEI=RKEI+(UX(I,J,L)*DUT(I,J,L)+VX(I,J,L)*DVT(I,J,L))
!          END DO
!          DUTIL=DUTIL+DUTI
!          RKEIL=RKEIL+RKEI
!        END DO
!        if (dopit) DUTIL=DUTIL+2.*DT1*RADIUS*OMEGA*COSV(J)*
!     *       (PI(J-1)*RAVPN(J-1)+PI(J)*RAVPS(J))
!        DAMB(J)=DUTIL*COSV(J)*RADIUS*mb2kg
!        DKEB(J)=RKEIL*mb2kg
!      END DO
!!$OMP END PARALLEL DO
!C****
!
!c
!c regrid to primary latitudes
!c
!      call regrid_to_primary_1d(damb)
!      N=NAMOFM(M)
!      DO J=J_0,J_1
!        CONSRV(J,N)=CONSRV(J,N)+DAMB(J)*BYDXYP(J)*BYIM
!      ENDDO
!      call regrid_to_primary_1d(dkeb)
!      N=NKEOFM(M)
!      DO J=J_0,J_1
!        CONSRV(J,N)=CONSRV(J,N)+DKEB(J)*BYDXYP(J)*BYIM
!      ENDDO
!      CALL TIMEOUT(BEGIN,MDIAG,MDYN)
!      RETURN
!      END SUBROUTINE DIAGCD
c      end module DIAG

      subroutine regrid_to_primary_1d(x)
      USE MODEL_COM, only : jm
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, GRID, NORTH,
     &     hasSouthPole, hasNorthPole

      implicit none
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: X
      integer :: j
      CALL HALO_UPDATE(grid, X, FROM=NORTH)
      if(hasSouthPole(grid)) X(1)=0.
      DO J=GRID%J_STRT,GRID%J_STOP_SKP
        X(J)=.5*(X(J)+X(J+1))
      ENDDO
      if(hasNorthPole(grid)) X(JM)=.5*X(JM)
      return
      end subroutine regrid_to_primary_1d

!      SUBROUTINE DIAG5D (M5,NDT,DUT,DVT)
!      USE MODEL_COM, only : im,imh,jm,lm,fim,
!     &     DSIG,JEQ,LS1,MDIAG,MDYN
!      !USE DIAG_COM, only : speca,nspher,klayer
!      USE ATMDYN, only : FCUVA,FCUVB
!      USE DOMAIN_DECOMP_1D, only : GRID,GET,GLOBALSUM, WRITE_PARALLEL
!      USE GETTIME_MOD
!      IMPLICIT NONE
!
!      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
!     &        DUT,DVT
!
!c      REAL*8, DIMENSION(0:IMH,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM,2) :: FCUVA,FCUVB
!c      COMMON/WORK7/FCUVA,FCUVB
!
!      INTEGER :: M5,NDT
!
!      REAL*8, DIMENSION(IMH+1) :: X
!      REAL*8, DIMENSION(0:IMH) :: FA,FB
!      REAL*8, DIMENSION(IMH+1,NSPHER) :: KE
!      REAL*8, DIMENSION
!     &  (IMH+1,GRID%J_STRT_HALO:GRID%J_STOP_HALO,NSPHER) :: KE_part
!      REAL*8 BEGIN
!
!      INTEGER :: J,J45N,KUV,KSPHER,L,MKE,N,NM
!      INTEGER :: J_0STG,J_1STG
!
!      CALL GET(GRID,J_STRT_STGR=J_0STG,J_STOP_STGR=J_1STG)
!
!      NM=1+IM/2
!      J45N=2.+.75*(JM-1.)
!      MKE=M5
!
!      GO TO (810,810,810,100,100,  100,810),M5
!C****  810 WRITE (6,910) M5
!  810 CALL WRITE_PARALLEL(M5, UNIT=6, format=
!     & "('0INCORRECT VALUE OF M5 WHEN CALLING DIAG5D.  M5=',I5)")
!C****  910 FORMAT ('0INCORRECT VALUE OF M5 WHEN CALLING DIAG5D.  M5=',I5)
!      call stop_model('INCORRECT VALUE OF M5 WHEN CALLING DIAG5D',255)
!C****
!C**** KINETIC ENERGY
!C****
!C**** TRANSFER RATES FOR KINETIC ENERGY IN THE DYNAMICS
!  100 CALL GETTIME(BEGIN)
!      KE(:,:)=0.
!      KE_part(:,:,:)=0.
!
!      DO L=1,LM
!        DO J=J_0STG,J_1STG
!          KSPHER=KLAYER(L)
!          IF (J > JEQ) KSPHER= KSPHER+1
!          DO KUV=1,2 ! loop over u,v
!            IF(KUV.EQ.1) CALL FFT(DUT(1,J,L),FA,FB)
!            IF(KUV.EQ.2) CALL FFT(DVT(1,J,L),FA,FB)
!            DO N=1,NM
!              X(N)=.5*FIM*
!     &          (FA(N-1)*FCUVA(N-1,J,L,KUV)+FB(N-1)*FCUVB(N-1,J,L,KUV))
!            ENDDO
!            X(1)=X(1)+X(1)
!            X(NM)=X(NM)+X(NM)
!            IF (J.NE.JEQ) KE_part(:,J,KSPHER)=KE_part(:,J,KSPHER)+
!     &                                        X(:)*DSIG(L)
!            IF (J.EQ.J45N) THEN     ! 45 N
!               KE_part(:,J,KSPHER+2)=KE_part(:,J,KSPHER+2)+X(:)*DSIG(L)
!            ELSE IF (J.EQ.JEQ) THEN ! EQUATOR
!              DO N=1,NM
!                KE_part(N,J,KSPHER+2)=KE_part(N,J,KSPHER+2)+
!     &                                X(N)*DSIG(L)
!                KE_part(N,J,KSPHER  )=KE_part(N,J,KSPHER  )+
!     &                                .5D0*X(N)*DSIG(L)       ! CONTRIB TO SH
!                KE_part(N,J,KSPHER+1)=KE_part(N,J,KSPHER+1)+
!     &                                .5D0*X(N)*DSIG(L)       ! CONTRIB TO NH
!              ENDDO
!              IF (KUV.EQ.2) KSPHER=KSPHER+1
!            ENDIF
!          ENDDO
!        ENDDO
!      ENDDO
!
!      CALL GLOBALSUM(grid, KE_part(1:NM,:,1:NSPHER), KE(1:NM,1:NSPHER),
!     &   ALL=.TRUE.)
!
!      DO 180 KSPHER=1,NSPHER
!      DO 180 N=1,NM
!  180 SPECA(N,MKE,KSPHER)=SPECA(N,MKE,KSPHER)+KE(N,KSPHER)/NDT
!C****
!      CALL TIMEOUT(BEGIN,MDIAG,MDYN)
!      RETURN
!      END SUBROUTINE DIAG5D
!
!      SUBROUTINE DIAG5F(UX,VX)
!C**** FOURIER COEFFICIENTS FOR CURRENT WIND FIELD
!C****
!      USE MODEL_COM, only : im,imh,jm,lm,
!     &     IDACC,MDIAG,MDYN
!      !USE DIAG_COM, only : ia_d5f
!      USE ATMDYN, only : FCUVA,FCUVB
!      USE DOMAIN_DECOMP_1D, only : GRID,GET
!      USE GETTIME_MOD
!      IMPLICIT NONE
!
!      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
!     &        UX,VX
!c      REAL*8, DIMENSION(0:IMH,JM,LM,2) :: FCUVA,FCUVB
!c      COMMON/WORK7/FCUVA,FCUVB
!      INTEGER :: J,L
!      INTEGER :: J_0STG, J_1STG
!      REAL*8 BEGIN
!
!      CALL GET(GRID, J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG)
!      CALL GETTIME(BEGIN)
!      IDACC(ia_d5f)=IDACC(ia_d5f)+1
!      DO L=1,LM
!         DO J=J_0STG,J_1STG
!            CALL FFT(UX(1,J,L),FCUVA(0,J,L,1),FCUVB(0,J,L,1))
!            CALL FFT(VX(1,J,L),FCUVA(0,J,L,2),FCUVB(0,J,L,2))
!         ENDDO
!      ENDDO
!      CALL TIMEOUT(BEGIN,MDIAG,MDYN)
!
!      RETURN
!      END SUBROUTINE DIAG5F

c      module ATMDYN_QDYNAM
c      USE ATMDYN
c      implicit none
c      private
c      public QDYNAM
c      contains
      SUBROUTINE QDYNAM( am_i_root )
!@sum  QDYNAM is the driver to integrate dynamic terms by the method
!@+          of pre-computing Courant limits using mean fluxes
!@+    It replaces CALL AADVT (MA,Q,QMOM, SD,PU,PV, DTLF,.TRUE.,
!@auth J. Lerner
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,q,dt,byim
      USE SOMTQ_COM, only : qmom
      !USE DIAG_COM, only: agc=>agc_loc
      !USE GCDIAG, only : jl_totntlh,jl_zmfntlh,jl_totvtlh,jl_zmfvtlh
      USE DYNAMICS, only: ps,mb,ma,sda
      USE DYNAMICS, only: mu=>pua, mv=>pva, mw=>sda
      USE TRACER_ADV, only:
     *    AADVQ,AADVQ0,sbf,sbm,sfbm,scf,scm,sfcm,ncyc
      USE DOMAIN_DECOMP_1D, only : grid, GET, halo_update, south, north

![eml
#ifdef CACHED_SUBDD
      USE subdd_mod, only : massfluxV1, massfluxW1, massfluxU1,
     &                      nmassflux, do_mass_flux
#endif
!eml]

      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: am_i_root

      REAL*8 DTLF,byncyc,byma
      INTEGER I,J,L   !@var I,J,L loop variables

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG)


      DTLF=2.*DT
      CALL CALC_AMP(PS,MB)
      CALL HALO_UPDATE(grid, MB, FROM=SOUTH+NORTH) ! for convenience later

![eml
#ifdef CACHED_SUBDD
      if(do_mass_flux) then
         massfluxU1(:,:,:)=MU(:,:,:) 
         massfluxV1(:,:,:)=MV(:,:,:)
         massfluxW1(:,:,:)=MW(:,:,1:LM-1)
         !if ( am_I_root ) then
         !   WRITE(6,*) 'CALL Accumulating Mass Flux', DTLF
         !endif
      endif
#endif
!eml]

      !IF ( am_I_root ) THEN
      !   WRITE(6,*) 'DTLF', DTLF
      !   WRITE(6,'(A6,144(E10.1))') 'MU', MU(:,:,1)
      !   WRITE(6,'(A6,144(E10.1))') 'MV', MV(:,:,1)
      !   WRITE(6,'(A6,144(E10.1))') 'MW', MW(:,:,1)
      !   CALL FLUSH(6)
      !END IF

      CALL AADVQ0 (1._8)  ! uses the fluxes pua,pva,sda from DYNAM
C****
C**** convert from concentration to mass units
C****
!$OMP PARALLEL DO PRIVATE (L,J,I)
      DO L=1,LM
      DO J=J_0,J_1
      DO I=1,IM
        Q(I,J,L)=Q(I,J,L)*MB(I,J,L)
        QMOM(:,I,J,L)=QMOM(:,I,J,L)*MB(I,J,L)
      enddo; enddo; enddo
!$OMP END PARALLEL DO
C**** ADVECT
        sfbm = 0.; sbm = 0.; sbf = 0.
        sfcm = 0.; scm = 0.; scf = 0.
      CALL AADVQ (Q,QMOM, .TRUE. ,'q       ')
        byncyc = 1./ncyc
!        AGC(:,:,jl_totntlh) = AGC(:,:,jl_totntlh) + sbf(:,:)
!        AGC(:,:,jl_zmfntlh) = AGC(:,:,jl_zmfntlh)
!     &    + sbm(:,:)*sfbm(:,:)*byim*byncyc
!        AGC(:,:,jl_totvtlh) = AGC(:,:,jl_totvtlh) + scf(:,:)
!        AGC(:,:,jl_zmfvtlh)  = AGC(:,:,jl_zmfvtlh)
!     &    + scm(:,:)*sfcm(:,:)*byim*byncyc
C****
C**** convert from mass to concentration units (using updated MA)
C****
!$OMP PARALLEL DO PRIVATE (L,I,J,BYMA)
      DO L=1,LM
      DO J=J_0,J_1
      DO I=1,IM
        BYMA = 1.D0/MA(I,J,L)
        Q(I,J,L)=Q(I,J,L)*BYMA
        QMOM(:,I,J,L)=QMOM(:,I,J,L)*BYMA
      enddo; enddo; enddo
!$OMP END PARALLEL DO

#ifndef TRACERS_ON
c Unscale the vertical mass flux accumulation for use by column physics.
c Switch the sign convention back to "positive downward".
      SDA(:,:,:) = -SDA(:,:,:)*NCYC
#else
c TRDYNAM will do the unscaling
#endif

      RETURN
      END SUBROUTINE QDYNAM
c      end module ATMDYN_QDYNAM

#ifdef TRACERS_ON
      SUBROUTINE TrDYNAM
!@sum  TrDYNAM is the driver to integrate tracer dynamic terms
!@auth J. Lerner
!@ver  1.0
      USE MODEL_COM, only: im,jm,lm,itime,dt,byim
      USE TRACER_COM, only: itime_tr0,trm,trmom,trname,t_qlimit,ntm
      USE TRACER_ADV
#ifndef SKIP_TRACER_DIAGS
      USE TRDIAG_COM, only: TAJLN=>TAJLN_loc, TAIJN=>TAIJN_LOC,
     *     jlnt_nt_tot,jlnt_nt_mm,jlnt_vt_tot,jlnt_vt_mm,
     *     tij_uflx,tij_vflx
#endif
      USE DYNAMICS, only : sda
      IMPLICIT NONE
      REAL*8 DTLF,byncyc
      INTEGER N

C**** uses the fluxes pua,pva,sda from DYNAM and QDYNAM
      DO N=1,NTM
        IF (itime.LT.itime_tr0(N)) cycle
        sfbm = 0.; sbm = 0.; sbf = 0.
        sfcm = 0.; scm = 0.; scf = 0.
        safv = 0.; sbfv = 0.

        CALL AADVQ (TRM(:,:,:,n),TrMOM(:,:,:,:,n),t_qlimit(n),trname(n))

C**** Flux diagnostics
#ifndef SKIP_TRACER_DIAGS
        byncyc = 1./ncyc
        TAJLN(:,:,jlnt_nt_tot,n) = TAJLN(:,:,jlnt_nt_tot,n) + sbf(:,:)
        TAJLN(:,:,jlnt_nt_mm, n) = TAJLN(:,:,jlnt_nt_mm, n)
     &    + sbm(:,:)*sfbm(:,:)*byim*byncyc
        TAJLN(:,:,jlnt_vt_tot,n) = TAJLN(:,:,jlnt_vt_tot,n) + scf(:,:)
        TAJLN(:,:,jlnt_vt_mm, n) = TAJLN(:,:,jlnt_vt_mm, n)
     &    + scm(:,:)*sfcm(:,:)*byim*byncyc

#ifdef TRACERS_WATER
#ifndef TRACERS_TOMAS 
C**** vertically integrated atmospheric fluxes
        TAIJN(:,:,tij_uflx,n) = TAIJN(:,:,tij_uflx,n) + safv(:,:)
        TAIJN(:,:,tij_vflx,n) = TAIJN(:,:,tij_vflx,n) + sbfv(:,:)
#endif
#endif
#endif

      ENDDO

c Unscale the vertical mass flux accumulation for use by column physics.
c Switch the sign convention back to "positive downward".
      SDA(:,:,:) = -SDA(:,:,:)*NCYC

      RETURN
      END SUBROUTINE TrDYNAM
#endif

      module UNRDRAG_COM
      !@sum  UNRDRAG_COM model variables for (alternative) gravity wave drag
      !@auth Tiehan Zhou / Marvin A. Geller
      !@ver  1.0
      USE MODEL_COM, only: JDPERY
      USE RESOLUTION, only: IM, JM
      implicit none
      save
      !@var r8: kind parameter of real*8
            integer, parameter :: r8 = selected_real_kind(12)
      !@var Z4var: Subgrid-scale orographic variances multiplied by 4 at uv grids
      !@+          of 2.5 x 2 resolution. It is real*4 rather than real*8.
      !@+          Its unit is m^2.
            real :: Z4var(IM, 2:JM)
      !@var Eke_by2: It is a tunable parameter from Eq.(3.1b) in McFarlane (JAS, 1987).
      !@+            Its unit is m^(-1).
            real(r8) :: Eke_by2 = 5.5E-6_r8
      !
      !@  Following parameters/variables are related to calculating nonorogrphic drag.
      !@+ The parameterization is described in Alexander and Dunkerton (JAS, 1999).
      !
      !   A. Parameters related to gravity wave source (assumed Gaussian shape in C):
      !
      !@var flag = 1 for B1 ( peak flux at c0 = 0 )
      !@+   flag = 0 for B2 ( peak flux at ci = 0 )
            integer, parameter :: flag = 0
      !@var Bt: sum of |momentum flux| for all +/-c (kg/m/s^2)
            real(r8) :: Bt(JM,JDPERY)
      !@var N_Kh: number of horizontal wavenumbers
            integer, parameter :: N_Kh = 1
      !@var Bm: amplitude for the spectrum (m^2/s^2) ~ u'w'
            real(r8), parameter :: Bm(N_Kh) = (/0.01_r8/)
      !@var Cw: half-width for the spectrum in C (m/s)
            real(r8), parameter :: Cw(N_Kh) = (/10.0_r8/)
      !@var [C_inf, C_sup]: the range of phase velocities for the Gaussian exp(-(C/Cw)**2)
            real(r8), parameter :: C_inf(N_Kh) = -3.0_r8 * Cw
            real(r8), parameter :: C_sup(N_Kh) =  3.0_r8 * Cw
      !@var N_C: number of C samples in source spectrum
            integer, parameter :: N_C = 100
      !@var dc: spectral resolution (m/s)
            real(r8), parameter :: dc(N_Kh) = (C_sup - C_inf)/(N_C - 1)
      !@var C: horizontal phase speed grid
            real(r8) :: C(N_C, N_Kh)
      !@var IZ0: vertical grid index of GW source (function of latitude)
            integer::IZ0(2:JM)
      !@var Wavelenth: wavelength of included gravity wave
      !@+              Unit of Wavelenth is km.
            real(r8), parameter :: Wavelenth(N_Kh) = (/100.0_r8/)
      !
      !    B. Other parameters and variables:
      !
      !@var N_Az: number of azimuths which is even.
            integer, parameter :: N_Az = 4
      !@var Kh: horizontal wave number grid
            real(r8) :: Kh(N_Kh)
      !@var Ah1, Ah2: used to compute components of velocity in azimuthal directions
            real(r8) :: Ah1(N_Az), Ah2(N_Az)
      !@param aLn2: ln(2.0)
            real(r8), parameter :: aLn2 = 0.69314718055994529_r8
      !@var L_min:
            integer :: L_min
      end module UNRDRAG_COM
      subroutine UNRDRAG (PB,U,V,T,SZ,UNRDRAG_x,UNRDRAG_y)
      !@sum  UNRDRAG is the driver for (alternative) gravity wave drag
      !@auth Tiehan Zhou / Marvin A. Geller
      !@ver  1.0
      USE UNRDRAG_COM
      USE CONSTANT, only : grav, bygrav, kapa, rgas
      USE GEOM, only: RAPVS, RAPVN
      USE MODEL_COM, only: IM, JM, LM, LS1
     *         , SIG, DSIG, SIGE, PSFMPT, PTOP, JDAY
      USE DOMAIN_DECOMP_1D, Only : GRID, GET
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : NORTH, SOUTH
      USE DOMAIN_DECOMP_1D, only : haveLatitude

      implicit none
      real(r8), dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *                   U, V, T, SZ, UNRDRAG_x, UNRDRAG_y
      real(r8), dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: PB
      intent(inout) :: PB, T, SZ
      intent(in) :: U, V
      intent(out) :: UNRDRAG_x, UNRDRAG_y
      real(r8), parameter :: byrgas = 1.0_r8/rgas
      real(r8), parameter :: dkapa = 1.0_r8 - kapa
      real(r8), parameter :: g_sq = grav * grav
      real(r8), parameter :: byPSFMPT = 1.0_r8/PSFMPT
      real(r8), dimension(LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     *                   T2, SZ2
      real(r8), dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: PB2
      real(r8), dimension(LM) :: dp, P_mid
      real(r8), dimension(LM+1) :: P_edge
      real(r8), dimension(LM) :: uc, vc, rho, bvf_sq, drag_x, drag_y
      real(r8), dimension(2:LM) :: ue, ve, rhoe, bvfe
      real(r8), dimension(2:LM) :: hb, U_Comp
      real(r8), dimension(LM) :: GW, GWF_X, GWF_Y
      real(r8) :: bvf_sqe
      !@var Eps: intermittency factor
      real(r8) :: Eps
      integer :: I, J, L, IP1
      integer :: Ikh, IC, IAZ
      integer :: IC0, MC
      real :: h_4sq
      real(r8) :: byPB2, by_dp_sum
      real(r8) :: Bsum
      real(r8) :: Ugw_S, Vgw_S, U_Comp_S
      real(r8) :: SGN, x
      real(r8) :: B(N_C,N_Az,N_kh)
      real(r8) :: C_actual(N_C)
      !
      !Extract domain decomposition info
      !
      integer :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S, J_0H, J_1H
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call GET(grid, J_STRT = J_0, J_STOP = J_1,             !&
     &         J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,   !&
     &         J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,     !&
     &         J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,     !&
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,            !&
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      if (HAVE_SOUTH_POLE) then
      do I = 2, IM
      PB(I,1) = PB(1,1)
      end do
      do L = 1, LM
      do I = 2, IM
      T(I,1,L) = T(1,1,L)
      SZ(I,1,L) = SZ(1,1,L)
      end do
      end do
      end if

      if (HAVE_NORTH_POLE) then
      do I = 2, IM
      PB(I,JM) = PB(1,JM)
      end do
      do L = 1, LM
      do I = 2, IM
      T(I,JM,L) = T(1,JM,L)
      SZ(I,JM,L) = SZ(1,JM,L)
      end do
      end do
      end if

      call HALO_UPDATE(GRID, T,  from=SOUTH)
      call HALO_UPDATE(GRID, SZ, from=SOUTH)
      call HALO_UPDATE(GRID, PB, from=SOUTH)

      do L= LS1, LM
      dp(L) = PSFMPT * DSIG(L)
      P_mid(L) = SIG(L) * PSFMPT + PTOP
      P_edge(L) = SIGE(L) * PSFMPT + PTOP
      end do
      P_edge(LM+1) = SIGE(LM+1) * PSFMPT + PTOP

      do J = J_0S, J_1
         I = IM
         do IP1 = 1, IM
            PB2(I,J) = (PB(I,J-1) + PB(IP1,J-1))*RAPVN(J-1)     !&
     &                  + (PB(I,J) + PB(IP1,J))*RAPVS(J)
            do L= 1, LM
            T2(L,I,J) = 0.25_r8 *
     *          (T(I,J-1,L) + T(IP1,J-1,L) + T(I,J,L) + T(IP1,J,L))
            SZ2(L,I,J) = 0.25_r8 *
     *       (SZ(I,J-1,L) + SZ(IP1,J-1,L) + SZ(I,J,L) + SZ(IP1,J,L))
            end do
         I = IP1
         end do
      end do

      Latitude:  do J = J_0STG, J_1STG
      Longitude: do I = 1, IM
      byPB2 = 1.0_r8/PB2(I,J)
      h_4sq = Z4var(I,J)
      !
      !Following loop calculates dp(1:LS1-1), P_mid(1:LS1-1), P_edge(1:LS1-1)
      !
      do L= 1, LS1-1
      dp(L) = PB2(I,J) * DSIG(L)
      P_mid(L) = SIG(L) * PB2(I,J) + PTOP
      P_edge(L) = SIGE(L) * PB2(I,J) + PTOP
      end do
      !
      !Following loop calculates rho(:) ,bvf_sq(:), uc(:), vc(:) at the middle levels.
      !
      do L= 1, LM
      rho(L) = P_mid(L)**dkapa / T2(L,I,J) * byrgas
      bvf_sq(L) = 2.0_r8*SZ2(L,I,J) / (dp(L)*T2(L,I,J)) * g_sq*rho(L)
      uc(L) = U(I,J,L)
      vc(L) = V(I,J,L)
      end do
      !
      !Following loop calculates rhoe(:) ,bvf_sqe(:), ue(:), ve(:) at the edge levels.
      !
      do L= 2, LM
      by_dp_sum = 1.0_r8 / (dp(L) + dp(L-1))
      rhoe(L) = (rho(L) * dp(L-1) + rho(L-1) * dp(L)) * by_dp_sum
      bvf_sqe = (bvf_sq(L) * dp(L-1) + bvf_sq(L-1) * dp(L)) * by_dp_sum
      if (bvf_sqe <= 0.0_r8) then
          bvfe(L) = -1.0
      else
          bvfe(L) = sqrt (bvf_sqe)
      end if
      ue(L) = (uc(L) * dp(L-1) + uc(L-1) * dp(L)) * by_dp_sum
      ve(L) = (vc(L) * dp(L-1) + vc(L-1) * dp(L)) * by_dp_sum
      hb(L) = P_edge(L) / rhoe(L) * bygrav
      end do

      if (h_4sq <= 0.0) then 
          !For efficiency, orogrphic variances over oceans were set to a negative value. 
          drag_x(:) = 0.0_r8
          drag_y(:) = 0.0_r8
      else
          call orographic_drag (ue, ve, rhoe, bvfe, h_4sq, Eke_by2 
     *                           , drag_x, drag_y)
      end if
          UNRDRAG_x(I,J,1:LS1-1) = drag_x(1:LS1-1) * byPB2
          UNRDRAG_y(I,J,1:LS1-1) = drag_y(1:LS1-1) * byPB2
          UNRDRAG_x(I,J,LS1:LM) = drag_x(LS1:LM) * byPSFMPT
          UNRDRAG_y(I,J,LS1:LM) = drag_y(LS1:LM) * byPSFMPT

      !...Calculating Eps
      Eps = 0.0_r8
      do Ikh = 1, N_Kh
         Bsum = 0.0_r8
         do IC = 1, N_C
            Bsum = Bsum + Bm(Ikh) * exp(-(C(IC,Ikh)/Cw(Ikh))**2 * aLn2)
         end do
         Eps = Eps + Bsum * dc(Ikh) * ( rhoe(IZ0(J)) * 100.0_r8 )
                                                       !!!100.0_r8 arises from the units of P and rho.
      end do
      Eps = Bt(J,JDAY) / Eps
      !...Calculating source spectra (function of azimuth, horizontal wave number)
      do IAZ = 1, N_Az
         Ugw_S = ue(IZ0(J))
         Vgw_S = ve(IZ0(J))
         U_Comp_S = Ugw_S * Ah1(IAZ) + Vgw_S * Ah2(IAZ)
         do Ikh = 1, N_Kh
            do IC = 1, N_C
               x = C(IC,Ikh)*(1-flag) + (C(IC,Ikh) - U_Comp_S)*flag
               if ( x == 0.0_r8 ) then
                  SGN = 0.0_r8
               else
                  SGN = sign(1.0_r8, x)
               end if
               B(IC,IAZ,Ikh) = SGN * Bm(Ikh) *
     *            exp( -(C(IC,Ikh)/Cw(Ikh))**2 * aLn2 ) * rhoe(IZ0(J))
            end do
         end do
      end do

      GWF_X(:) = 0.0_r8
      GWF_Y(:) = 0.0_r8
      !...Loop over azimuth
      do IAZ = 1, N_Az
         U_Comp(:) = ue(:) * Ah1(IAZ) + ve(:) * Ah2(IAZ)
         !...Loop over horizontal wave number
         do Ikh = 1, N_kh
            IC = N_C
            IC0 = 1
            MC = 0
            do while ( B(IC,IAZ,Ikh) > 0.0_r8 )
               MC = MC + 1
               IC0 = IC
               IC = IC - 1
               if ( IC == 0) exit
            end do
         !...if flag == 0 then parameterize spectrum with c_actual = c - ucomp(iz0)
         C_actual(:) = C(:,Ikh) * real(flag, r8)  !&
     &          + ( C(:,Ikh) + U_Comp(IZ0(J)) ) * real(1 - flag, r8)
         call nonorographic_drag (C_actual(IC0), dc(Ikh), B(IC,IAZ,Ikh)
     &           , Eps, Kh(Ikh), hb, rhoe, U_Comp, bvfe, MC, IZ0(J), GW)
         GWF_X(:) = GWF_X(:) + GW(:) * Ah1(IAZ)
         GWF_Y(:) = GWF_Y(:) + GW(:) * Ah2(IAZ)
         end do      !horizontal wavenumber grid
      end do      !azimuth grid
      do L = L_min, LM
      if (L < LS1) then
         UNRDRAG_x(I,J,L) = UNRDRAG_x(I,J,L) + GWF_X(L) * byPB2
         UNRDRAG_y(I,J,L) = UNRDRAG_y(I,J,L) + GWF_Y(L) * byPB2
      else
         UNRDRAG_x(I,J,L) = UNRDRAG_x(I,J,L) + GWF_X(L) * byPSFMPT
         UNRDRAG_y(I,J,L) = UNRDRAG_y(I,J,L) + GWF_Y(L) * byPSFMPT
      end if
      end do

      end do Longitude
      end do Latitude
      end subroutine UNRDRAG

      subroutine init_UNRDRAG
      !@sum  init_UNRDRAG initializes parameters for (alternative) gravity wave drag
      !@auth Tiehan Zhou / Marvin A. Geller
      !@ver  1.0
      USE RESOLUTION, only: JM, LM, PLbot
      USE CONSTANT, only : pi, twopi
      USE GEOM, only: LAT_DG
      USE MODEL_COM, only: JDPERY
      USE UNRDRAG_COM, only: Z4var, Bt
      USE UNRDRAG_COM, only: r8, N_C, C_inf, dc, C, IZ0, N_Kh, Wavelenth
      USE UNRDRAG_COM, only: Kh, Ah1, Ah2, N_Az, aLn2, L_min
      !USE FILEMANAGER, only: openunit, closeunit
      implicit none
      integer :: iu_Z4var, I, IAZ, J, IT
      real(r8) :: x, Phi
      real(r8) :: Bt_Smax, Bt_Nmax
      real(r8) :: Bt_Tmax
      character(Len=80) :: Title
      !call openunit("Z4var", iu_Z4var, .true., .true.)
      !read(iu_Z4var) Title, Z4var
      !call closeunit(iu_Z4var)

      Bt_Smax = 6.0_r8 * 0.001_r8
      Bt_Nmax = 0.5_r8 * 0.001_r8
      do IT = 1, JDPERY
         x = cos ( twopi * real(IT - 16, r8) / real(JDPERY, r8) )
         do J = 1, JM
            if ( LAT_DG(J,2) <= 1.0E-8 .and. x <= 0.0_r8 ) then
               Bt(J,IT) = -Bt_Smax *
     *          exp(-((LAT_DG(J,2) + 60.0_r8)/15.0_r8)**2 * aLn2 ) * x
            elseif ( LAT_DG(J,2) > 1.0E-8 .and. x >= 0.0_r8 ) then
               Bt(J,IT) =  Bt_Nmax *
     *          exp(-((LAT_DG(J,2) - 60.0_r8)/15.0_r8)**2 * aLn2 ) * x
            else
               Bt(J,IT) = 0.0_r8
            end if
         end do
      end do

      Bt_Tmax = 0.5_r8 * 0.001_r8
      do IT = 1, JDPERY
         x = cos ( twopi * real(IT - 16, r8) / real(JDPERY, r8) )
         Phi = -10.0_r8 * x
         do J = 1, JM
            Bt(J,IT) = Bt(J,IT) + Bt_Tmax *
     *          exp(-( (LAT_DG(J,2) - Phi)/5.0_r8 )**2 * aLn2 ) *
     *          0.25_r8 * ( 3.0_r8 - x )
         end do
      end do

      Bt = Bt + 1.0_r8 * 0.001_r8

      do I = 1, N_C
         C(I, :) = C_inf(:) + real(I - 1, r8) * dc(:)
      end do
      do I = 1, N_Kh
      Kh(I) = twopi / (Wavelenth(I) * 1000.0_r8)
                            !!!Factor 1000.0 arises from the unit of Wavelenth.
      end do
      do IAZ = 1, N_Az
      x = twopi / real(N_Az, r8) * real(IAZ - 1, r8)
      Ah1(IAZ) = cos(x)
      Ah2(IAZ) = sin(x)
      end do
      I = 1
      do while ( PLbot(I) >= 100.0_r8 )
         I = I + 1
         if ( I == LM + 2 ) exit
      end do
         IZ0(:) = I - 1
      L_min = minval(IZ0)
      end subroutine init_UNRDRAG

      subroutine orographic_drag (u,v,rho, bvf,h_4sq,coef,drag_x,drag_y)
      !@sum   orographic_drag
      !@auth Tiehan Zhou / Marvin A. Geller
      !@ver  1.0
      USE CONSTANT, only : grav
      USE MODEL_COM, only: LM, byDSIG
      USE UNRDRAG_COM, only : r8
      implicit none
      real(r8), intent(in) :: coef
      real(r8), dimension(2:LM), intent(in) :: u, v, rho, bvf
      real(r8), dimension(LM), intent(out) :: drag_x, drag_y
      real, intent(in) :: h_4sq
      !@param  byFc_sq: inverse Froude number squared
      real(r8), parameter :: byFc_sq = 0.5_r8
      real(r8), parameter :: Fc_sq = 1.0_r8/byFc_sq
      real(r8) :: flux(2:LM+1)
      real(r8) :: wind(2:LM)
      real(r8) :: he_sq, u0_sq, u0, flux_temp
      real(r8) :: Fr_sq, const, Fr_sq_min
      real(r8) :: drag, proj_x, proj_y
      integer :: L
      do L = 1, LM
      flux(L+1) = 0.0_r8
      drag_x(L) = 0.0_r8
      drag_y(L) = 0.0_r8
      end do
      Fr_sq_min = Fc_sq
      u0_sq = u(2) * u(2) + v(2) * v(2)
      u0 = sqrt(u0_sq)
      if (u0 < 1.0E-8_r8) then
          return
      else
          wind(2) = u0
          do L = 3, LM
          wind(L) = (u(L) * u(2) + v(L) * v(2)) / u0
          end do
      end if

      if (bvf(2) <= 0.0_r8) then
          return
      else
          he_sq = min(h_4sq, byFc_sq * u0_sq / ( bvf(2) * bvf(2) ) )
      end if

      const = he_sq * rho(2) * bvf(2) * wind(2)
      flux(2) = -coef * const
      flux_temp = flux(2) * byFc_sq

      do L = 3, LM
      if (wind(L) < 1.0E-8_r8) then
          flux(L) = 0.0_r8
          exit
      elseif (bvf(L) <= 0.0_r8) then
          flux(L) = flux(L-1)
      else
          Fr_sq = rho(L) * wind(L)**3 / ( const * bvf(L) )
          if (Fr_sq >= Fr_sq_min) then
          flux(L) = flux(L-1)
          else
          flux(L) = flux_temp * Fr_sq
          Fr_sq_min = Fr_sq
          end if
      end if
      end do

      proj_x = u(2) / u0
      proj_y = v(2) / u0

      do L = 2, LM
      drag = -grav * (flux(L+1) - flux(L)) * byDSIG(L)
      drag_x(L) = drag * proj_x
      drag_y(L) = drag * proj_y
      end do
      end subroutine orographic_drag

      subroutine nonorographic_drag (c,dc,b,eps,kh,hb,rho
     *                                ,u,bf,nc,iz0, gwfrc)
      !@sum   nonorographic_drag
      !@auth Tiehan Zhou / Marvin A. Geller
      !@ver  1.0
      USE CONSTANT, only : grav, by3
      USE MODEL_COM, only: LM, byDSIG
      USE UNRDRAG_COM, only : r8
      !===============================================
      !...AD parameterizaion with arbitrary tabulated
      !...momentum flux (phase speed) spectrum: rho*u'*v' = b(c)*dc
      !...computes force gwfrc in a single azimuthal direction
      !
      !...Input:
      !...c(1:nc) phase speed (m/s)
      !...b(1:nc) mom flux density (kg/m**2/s)
      !...kh=horizontal wavenumber (1/m)
      !...eps=intermittency factor
      !...rho(2:LM)=density profile (kg/m**3)
      !...u(2:LM)=wind profile (m/s)
      !...bf(2:LM)=bouyancy frequency (1/s)
      !...iz0=source grid level
      !
      !...gwfrc(1:LM)=output force
      !

      implicit none

      !...arguements
      integer, intent(in) :: nc, iz0
      real(r8), intent(in), dimension(nc) :: c, b                  !global variables
      real(r8), intent(in), dimension(2:LM) :: hb, rho, u, bf   !global variables
      real(r8), intent(out) :: gwfrc(LM)                        !global variables
      real(r8), intent(in) :: kh, dc                            !global variables
      real(r8), intent(in) :: eps                               !intermittency factor
      !...The following are global variables
      real(r8) :: crfl(iz0:LM+1)      !phase speed of turning point at z
      integer :: n_rfl(iz0:LM+1)
      real(r8) :: dc3(iz0:LM+1), total_flux(iz0:LM+1)
      real(r8) :: k2, alp2, cm
      !!!...The following are local variables
      integer :: i, j, number, n1, n2, n3, n_end
      real(r8)::unsat              !unsat > 0 if unsaturated at z
      real(r8)::const, dc1, dc2, s1, s2, ratio, dc_end
      !------------------------------------------------------------------
      gwfrc(:) = 0.0_r8  !necessary when the subroutine 'unsaturated' returns quickly.

      if ( nc == 0 ) return

      k2 = kh * kh
      total_flux(:) = 0.0_r8

      !find the turning point, which must decrease with altitude
      cm = 2000.0
      do i = iz0, LM
      alp2 = 1.0_r8 / (4.0_r8 * hb(i) * hb(i) )
      if (bf(i) < 0.0_r8) then
         crfl(i) = u(i)
      else
         crfl(i) = u(i) + bf(i) / sqrt( k2 + alp2 )
      end if
      crfl(i) = min ( crfl(i), cm )
      cm = crfl(i)
      end do
      n_rfl(:) = int( ( crfl(:) - c(1) ) / dc ) + 1
      dc3(:) = crfl(:) - ( c(1) + dc * real ( n_rfl(:) - 1 ) )
      i = iz0
      if ( crfl(i) < c(nc) .and. c(1) < crfl(i) ) then
         n_end = n_rfl(i)
         dc_end = dc3(i)
      else if ( crfl(i) >= c(nc) ) then
         n_end = nc - 1
         dc_end = dc
      else
         return
      end if

      !...find unsaturated waves in the spectrum at z(i)
      if (bf(i) > 0.0_r8) then
         const = kh * rho(i) / bf(i) / 2.0_r8
         number = 0
         SOURCE: do j = 1, n_end
            unsat = const * ( c(j) - u(i) )**3 - b(j)
            if ( number == 0 .and. unsat > 0.0_r8 ) then
               number = 1
               if ( j == 1 ) then
                  n1 = j
                  dc1 = 0.0_r8
               else
                  n1 = j -1
                  s2 = const * ( c(n1) - u(i) )**3 - b(n1)
                  s1 = const * ( c(n1+1) - u(i) )**3 - b(n1+1)
                  ratio = s2 / ( s2 - s1 )
                  dc1 = ratio * dc
               end if
            else if ( number == 1 .and. unsat <= 0.0_r8 ) then
               number = 2
               n2 = j - 1
               s2 = const * ( c(n2) - u(i) )**3 - b(n2)
               s1 = const * ( c(n2+1) - u(i) )**3 - b(n2+1)
               ratio = s2 / ( s2 - s1 )
               dc2 = ratio * dc
            else if( ( number == 2 ) .and. ( unsat > 0.0_r8 ) ) then
               number = 3
               n3 = j - 1
               exit SOURCE
            end if
         end do SOURCE
      else
         number = 1
         n1 = 1
         dc1 = 0.0_r8
      end if

      if ( number == 0 ) return

      BAND_1_2: if ( number == 1 ) then
         n2 = n_end
         dc2 = dc_end
         call unsaturated ( n1, dc1, n2, dc2, i )
      else if ( number >= 2 ) then BAND_1_2
         call unsaturated ( n1, dc1, n2, dc2, i )
         BAND_2:if ( number == 3 ) then
            n1 = n3
            s2 = const * ( c(n1) - u(i) )**3 - b(n1)
            s1 = const * ( c(n1+1) - u(i) )**3 - b(n1+1)
            ratio = s2 / ( s2 - s1 )
            dc1 = ratio * dc
            n2 = n_end
            dc2 = dc_end
            call unsaturated ( n1, dc1, n2, dc2, i )
         end if BAND_2
      end if BAND_1_2

      total_flux(LM) = total_flux(LM-2) * by3
      total_flux(LM-1) = total_flux(LM-2) * 2.0_r8 * by3

      do i = iz0, LM
         gwfrc(i) = -grav * ( total_flux(i+1) - total_flux(i) ) *
     *                                          byDSIG(i) * eps
      end do

      contains

      recursive subroutine unsaturated ( nz1, dcz1, nz2, dcz2, level)
      implicit none
      integer, intent(in)::nz1, nz2, level
      real(r8), intent(in)::dcz1, dcz2
      integer::number, n1, n2, n3, n_end, n_start
      integer::i, j
      real(r8)::const, unsat, dc1, dc2, s1, s2, ratio, flux_rfl, dc_end

      i = level
      do j = nz1, nz2 - 1
         total_flux(i) = total_flux(i) + b(j) * dc
      end do
      total_flux(i) = total_flux(i) - b(nz1) * dcz1 + b(nz2) * dcz2

      if (level == LM) return

      i = level + 1
      if ( c(nz1) + dcz1 < crfl(i) .and. crfl(i) < c(nz2) + dcz2 ) then
         n_end = n_rfl(i)
         dc_end = dc3(i)
      else if ( crfl(i) >= c(nz2) + dcz2 ) then
         n_end = nz2
         dc_end = dcz2
      else
         n_end = nz1
         dc_end = dcz1
      end if

      flux_rfl = 0.0_r8
      do j = n_end, nz2 - 1
      flux_rfl =  flux_rfl + b(j) * dc
      end do
      flux_rfl = flux_rfl - b(n_end) * dc_end + b(nz2) * dcz2
      total_flux(:level) = total_flux(:level) - flux_rfl

      if ( crfl(i) <= c(nz1) + dcz1 ) return

      if (bf(i) > 0.0_r8) then
         const = kh * rho(i) / bf(i) / 2.0_r8
         number = 0
         if ( nz1 == 1 .and. dcz1 == 0.0_r8 ) then
            n_start = nz1
         else
            n_start = nz1 + 1
         end if
         find_bands: do j = n_start, n_end
            unsat = const * ( c(j) - u(i) )**3 - b(j)
            upward: if ( number == 0 .and. unsat > 0.0_r8 ) then
               number = 1
               if ( j == 1 ) then
                  n1 = j
                  dc1 = 0.0_r8
               else
                  n1 = j - 1
                  s2 = const * ( c(n1) - u(i) )**3 - b(n1)
                  s1 = const * ( c(n1+1) - u(i) )**3 - b(n1+1)
                  ratio = s2 / ( s2 - s1 )
                  dc1 = ratio * dc
                  if ( n1 == nz1 .and. s1 * s2 <= 0.0_r8 ) then
                     dc1 = max ( dc1, dcz1 )
                  else if ( n1 == nz1 .and. s1 * s2 > 0.0_r8 ) then
                     dc1 = dcz1
                  end if
               end if
            else if ( number == 1 .and. unsat <= 0.0_r8 ) then upward
               number = 2
               n2 = j - 1
               s2 = const * ( c(n2) -u(i) )**3 - b(n2)
               s1 = const * ( c(n2+1) - u(i) )**3 - b(n2+1)
               ratio = s2 / ( s2 - s1 )
               dc2 = ratio * dc
               if ( n2 == n_end ) dc2 = min ( dc2, dc_end )
            else if ( number == 2 .and. unsat > 0.0_r8 ) then upward
               number = 3
               n3 = j - 1
               exit find_bands
            end if upward
         end do find_bands
      else
         number = 1
         n1 = nz1
         dc1 = dcz1
      end if

      if ( number == 0 ) return

      outgoing_1_2: if ( number == 1 ) then
         n2 = n_end
         dc2 = dc_end
         call unsaturated ( n1, dc1, n2, dc2, i )
      else if ( number >= 2 ) then outgoing_1_2
         call unsaturated ( n1, dc1, n2, dc2, i )
         outgoing_2: if ( number == 3 ) then
            n1 = n3
            s2 = const * ( c(n1) - u(i) )**3 - b(n1)
            s1 = const * ( c(n1+1) - u(i) )**3 - b(n1+1)
            ratio = s2 / ( s2 - s1 )
            dc1 = ratio * dc
            n2 = n_end
            dc2 = dc_end
            call unsaturated ( n1, dc1, n2, dc2, i )
         end if outgoing_2
      end if outgoing_1_2

      end subroutine unsaturated
      end subroutine nonorographic_drag
