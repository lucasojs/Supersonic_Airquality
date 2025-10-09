#include "rundeck_opts.h"

!==============================================================================
! LTM: MODEL_COM.f
!==============================================================================

      MODULE MODEL_COM
!@sum  MODEL_COM Main model variables, independent of resolution
!@auth Original Development Team
!@ver  1.0
      USE RESOLUTION, only : im,jm,lm,ls1,kep,istrat,
     *     psf,pmtop,ptop,psfmpt,pstrat,plbot
!      USE DOMAIN_DECOMP_ATM, only : grid
!AOO      USE ESMF_CUSTOM_MOD, ONLY: FIELD
c$$$#ifdef USE_FVCORE
c$$$      USE ESMF_MOD, only: esmf_clock
c$$$#endif
      IMPLICIT NONE
! just to make all compilers happy (should check later)
#if ! defined(COMPILER_G95)
      SAVE
#endif

!@param IMH half the number of latitudinal boxes
      INTEGER, PARAMETER :: IMH=IM/2
!@param IVNP,IVSP V at north/south pole is stored in U(IVNP,JM)/U(IVSP,1)
      INTEGER, PARAMETER :: IVNP = IM/4 , IVSP = 3*IM/4
!@param FIM,BYIM real values related to number of long. grid boxes
      REAL*8, PARAMETER :: FIM=IM, BYIM=1./FIM
!@param JEQ grid box zone around or immediately north of the equator
      INTEGER, PARAMETER :: JEQ=1+JM/2

      CHARACTER*132 XLABEL !@var XLABEL=runID+brief description of run
      INTEGER :: LRUNID    !@var Run name stored in XLABEL(1:LRUNID)
      INTEGER :: LHEAD=15  !@var length of crucial beg of module_headers

!@var LM_REQ Extra number of radiative equilibrium layers
      INTEGER, PARAMETER :: LM_REQ=3
!@var REQ_FAC/REQ_FAC_M factors for REQ layer pressures
      REAL*8, PARAMETER, DIMENSION(LM_REQ-1) ::
     *     REQ_FAC=(/ .5d0, .2d0 /)               ! edge
      REAL*8, PARAMETER, DIMENSION(LM_REQ) ::
     *     REQ_FAC_M=(/ .75d0, .35d0, .1d0 /),    ! mid-points
     *     REQ_FAC_D=(/ .5d0,  .3d0,  .2d0 /)     ! delta

!**** Vertical resolution dependent variables (set in INPUT)
!@var SIGE sigma levels at layer interfaces (1)
!!!!  Note:   sige(1)=1,  sige(ls1)=0,  sige(lm+1)=-pstrat/psfmpt
      REAL*8, DIMENSION(LM+1) :: SIGE
!@var SIG,DSIG,byDSIG mid point, depth, 1/depth of sigma levels (1)
      REAL*8, DIMENSION(LM) ::
     &     SIG,    ! = (sige(1:lm)+sige(2:lm+1))*0.5d0,
     &     DSIG,   ! =  sige(1:lm)-sige(2:lm+1),
     &     byDSIG  ! =  1./DSIG

!@var PL00, PMIDL00, PDSIGL00, AML00 press (mb), mid-pressure (mb),
!@+        mass (kg/m2) for mean profile
!@var PEDNL00 edge pressure for mean profile (mb)
      REAL*8, DIMENSION(LM+LM_REQ) :: PL00, PMIDL00, PDSIGL00, AML00,
     *     BYAML00
      REAL*8, DIMENSION(LM+LM_REQ+1) :: PEDNL00

!**** Model control parameters:
!@dbparam KOCEAN: if 0 => specified, if 1 => predicted ocean
!@dbparam MFILTR: if 1 => SLP, if 2 => T, if 3 => SLP&T is filtered
      integer :: KOCEAN = 1, MFILTR = 1
!@dbparam COUPLED_CHEM: if 0 => uncoupled, if 1 => coupled
      integer :: COUPLED_CHEM = 0
!@var ij_debug: if i > 0, print out some extra info on bad ij box
      integer, dimension(2) :: ij_debug = (/ 0 , 1 /)

!**** Default simulation year. If set to zero, transient run.
!@dbparam master_yr year of simulation. This value will define aero_yr,
!@+       aer_int_yr, albsn_yr, crops_yr, ghg_yr, o3_yr, s0_yr, volc_yr,
!@+       variable_orb_par,orb_par_year_bp,and PI_run, unless any of these
!@+       are specifically defined.
      INTEGER ::  master_yr = 1951

!**** Controls on FLTRUV (momentum/velocity filter)
!@dbparam DT_XUfilter dU is multiplied by dt/DT_XUfilter in E-W
!@dbparam DT_XVfilter dV is multiplied by dt/DT_XVfilter in E-W
!@dbparam DT_YUfilter dU is multiplied by dt/DT_YUfilter in N-S
!@dbparam DT_YVfilter dV is multiplied by dt/DT_YVfilter in N-S
      REAL*8 :: DT_XUfilter=0. ! U-filter is NOT used in E-W direction
      REAL*8 :: DT_XVfilter=0. ! V-filter is NOT used in E-W direction
      REAL*8 :: DT_YUfilter=0. ! U-filter is NOT used in N-S direction
      REAL*8 :: DT_YVfilter=0. ! V-filter is NOT used in N-S direction
!     Target Coordinates for SCM
      INTEGER*4 :: I_TARG,J_TARG   !TWP I=125,J=39  set targets in parameter list
      INTEGER*4 :: NSTEPSCM=0      !Time step counter for SCM
!@var QUVfilter: True if any of DT_[XY][UV]filter are not=0
      LOGICAL :: QUVfilter
!@dbparam ang_uv =1 to conserve ang mom in UVfilter
      INTEGER :: ang_uv = 1 ! UV-filter conserves ang mom

!@dbparam X_SDRAG.  SDRAG ~X_SDRAG(1)+X_SDRAG(2)*wind_magnitude
      REAL*8, DIMENSION(2) :: X_SDRAG = (/2.5D-4,2.5D-5/)
!@dbparam C_SDRAG.  SDRAG=C_SDRAG (const.) above PTOP
      REAL*8 :: C_SDRAG = 2.5D-5
      REAL*8, DIMENSION(LS1:LM) :: CSDRAGL
!@dbparam P_CSDRAG pressure level above which const.drag is increased
      REAL*8 :: P_CSDRAG=0.
!@dbparam P(P)_SDRAG pressure level above which SDRAG is applied (mb)
      REAL*8 :: P_SDRAG=0., PP_SDRAG = 1.d0 ! (PP_... near poles)
!@var L(P)SDRAG level above which SDRAG is applied (near pole)
      INTEGER :: LSDRAG=LM, LPSDRAG=LM  ! non-polar, polar limit
!@var ANG_SDRAG if =1: ang.momentum lost by SDRAG is added in below PTOP
      INTEGER :: ANG_SDRAG=1  ! default: SDRAG does conserve ang.mom
!@dbparam Wc_JDRAG critical velocity for J.Hansen/Judith Perlwitz drag
      REAL*8 :: Wc_JDRAG=30.d0  !  if 0.: no JDRAG-feature in Sdrag
!@dbparam wmax imposed limit for stratospheric winds (m/s) in SDRAG
      real*8 :: wmax=200.d0
!@dbparam do_polefix if =1 : u,v tendencies are corrected near the pole
      INTEGER :: do_polefix=1     ! default is to enable corrections
!@dbparam VSDRAGL layer dependent tuning factor for stratospheric drag
!@+   (not =1 e.g. if used with explicit grav.wave drag scheme)
      real*8 :: VSDRAGL(LS1:LM) = 1d0
!@dbparam  USE_UNR_DRAG   if 1 =>SDRAG is turned off and GWD is applied.
!@+    if 0 => SDRAG is kept intact and alternative GWD is not employed.
      INTEGER :: USE_UNR_DRAG=0  ! default: SDRAG is kept intact.

!**** Diagnostic control parameters
!@dbparam KCOPY: if 1 => acc, if 2 => +rsf, if 3 => +od are saved
!@dbparam NMONAV number of months in a diagnostic accuml. period
!@dbparam Kvflxo if 1 => vert.fluxes into ocean are saved daily
!@dbparam Kradia if -1 save data for, if 1|2 do   inst|adj forcing run
!@dbparam NIPRNT number of instantaneous initial printouts
      integer :: KCOPY=2, NMONAV=1, Kvflxo=0, Kradia=0,iu_rad, NIPRNT=1

C**** (Simplified) Calendar Related Terms
!@param JDperY,JMperY    number of days,months per year
!@var   JDendOfM(0:12)   last Julian day in month
!@var   JDmidOfM(0:13)   middle Julian day in month
      integer, PARAMETER :: JDPERY = 365, JMPERY = 12
      integer :: JDendOfM(0:JMPERY) = (
     *     /0,31,59,90,120,151,181,212,243,273,304,334,365/)
      integer :: JDmidOfM(0:JMPERY+1) = (
     *     /-15,16,45,75,106,136,167,197,228,259,289,320,350,381/)

!@var AMON,AMONTH(0:12)  (3-4 letter) names for current,all months
!@var AMON0  (3-4 letter) name of first month of the current acc-period
      CHARACTER*4 :: AMON='none',AMON0='none', AMONTH(0:12) = (/'IC  ',
     *  'JAN ','FEB ','MAR ','APR ','MAY ','JUNE',
     *  'JULY','AUG ','SEP ','OCT ','NOV ','DEC '/)

!@var NDAY and IYEAR1 relate CALENDAR TIME and INTERNAL TIME Itime :
!@var NDAY number of Internal Time Units per day (1 ITU = DTsrc sec)
!@nlparam IYEAR1  year 1 of internal clock (Itime=0 to 365*NDAY)
      INTEGER :: NDAY,IYEAR1=-1   !@var relate internal to calendar time

!@var ITIME current time in ITUs (1 ITU = DTsrc sec, currently 1 hour)
!@var JDAY,JMON,JDATE,JYEAR,JHOUR current Julian day,month,day,year,hour
      INTEGER :: Itime,JDAY,JMON,JDATE,JYEAR,JHOUR
!@var ItimeI,ItimeE   time at start,end of run
!@var Itime0          time at start of current accumulation period
!@var JMON0,JDATE0,JYEAR0,JHOUR0 date-info about Itime0 (beg.of acc.per)
      INTEGER :: ItimeI,ItimeE,   Itime0,JMON0,JDATE0,JYEAR0,JHOUR0

!@var ESMF clock required for some interfaces
c$$$#ifdef USE_FVCORE
c$$$      Type (ESMF_CLOCK) :: clock
c$$$#endif
!@dbparam DTSRC source time step (s)   = 1 ITU
      REAL*8 :: DTsrc = 3600.
!@dbparam DT (atmospheric) dynamics time step (s)
      REAL*8 :: DT    =  450.         ! DT = DTdyn_atm

C**** Time step related multipliers:  N... NI...
C**** general rule:   DTxxx = Nxxx*DTsrc  and  DTxxx = DTsrc/NIxxx
C**** except that the time steps related to NDAa, NDA5k, NDAsf are
C**** slightly larger, to sample all points within the cycle

!@var     NIdyn:  DT atm_dyn  =  DTsrc/NIdyn     (NIdyn=DTsrc/DT)
!@dbparam NIsurf: DT_Surface  =  DTsrc/NIsurf
!@dbparam NRad:   DT_Rad      =  NRad*DTsrc
!@dbparam NFILTR: DT_filter   =  NFILTR*DTsrc
      INTEGER :: NIdyn, NIsurf = 2, NRad = 5 , NFILTR = 1

!@dbparam Ndisk:  DT_saversf    =  Ndisk *DTsrc fort.1/fort.2 saves
!@dbparam Nssw:   DT_checkSsw   =  Nssw  *DTsrc
      INTEGER :: NDisk = 24, Nssw = 1

!@dbparam NDAA:   DT_DiagA    =  NDAA*DTsrc + 2*DT(dyn)
!@dbparam NDA5k:  DT_Diag5k   =  NDA5k*DTsrc + 2*DT(dyn) SpAnal KE
!@dbparam NDA5d:  DT_Diag5d   =  NDA5d*DTsrc     Consrv  SpAnal dyn
!@dbparam NDA5s:  DT_Diag5s   =  NDA5s*DTsrc     Consrv  SpAnal src
!@dbparam NDASf:  DT_DiagSrfc =  NDASf*DTsrc + DTsrc/NIsurf
!@dbparam NDA4:   DT_Diag4    =  NDA4 *DTsrc   Energy history
      INTEGER :: NDAa=7, NDA5d=1, NDA5k=7, NDA5s=1, NDASf=1, NDA4=24

!**** Accounting variables
!@dbparam IRAND last seed used by rand.number generator
!@var KDISK next rsf (fort.)1 or 2 to be written to
!@var NSTEP number of dynamics steps since start of run
!@var MRCH  flags position in dynamics cycle (>0 fw, <0 bw step)
      INTEGER :: IRAND=123456789, KDISK=1, NSTEP,MRCH
!@param rsf_file_name names of restart files
      CHARACTER(6), PARAMETER :: rsf_file_name(2)=(/'fort.1','fort.2'/)
!@var MODRD,MODD5K,MODD5S: if MODxxx=0 do xxx, else skip xxx
      INTEGER :: MODRD, MODD5K, MODD5S
!@var MDYN,MCNDS,MRAD,MSURF,MDIAG,MELSE timing-indices
      INTEGER  MDYN,MCNDS,MRAD,MSURF,MDIAG,MELSE
!@param NSAMPL number of diagnostic sampling schemes
      INTEGER, PARAMETER :: NSAMPL = 12
!@var IDACC(NSAMPL) counters for diagn. accumulations
      INTEGER, DIMENSION(NSAMPL) :: IDACC


!**** IO read/write flags used by the io_xyz routines
!@param IOWRITE Flag used for writing normal restart files
!@param IOWRITE_SINGLE Flag used for saving diags in single precision
!@param IOWRITE_MON Flag used for saving restart part only (no diags)
!@param IOREAD Flag used for reading in (composite) restart files
!@param IOREADNT Flag used for reading in restart files (w/o tracers)
!@param IRSFIC Flag used for reading in restart part to start NEW run
!@param IRSFICNT Flag used for reading restart (w/o tracers) for NEW run
!@param IRSFICNO Flag used for reading restart (w/o ocean) for NEW run
!@param IRERUN Flag used for reading in restart part to extend OLD run
      INTEGER, PARAMETER :: ioread=1,ioread_single=2,
     *     irerun=3,irsfic=4,irsficnt=5,ioreadnt=6,irsficno=7,
     *     ioread_nodiag=8,
     *     iowrite=-1,iowrite_single=-2,iowrite_mon=-3

!**** Main model prognostic variables
!@var U,V east-west, and north-south velocities (m/s)
!@var T potential temperature (referenced to 1 mb) (K)
!@var Q specific humidity (kg water vapor/kg air)
!@var WM cloud liquid water amount (kg water/kg air)
#ifdef BLK_2MOM
!@var WMICE cloud ice amount (kg water/kg air)
#endif
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: U
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: V
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: T
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: Q
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: WM
#ifdef BLK_2MOM
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: WMICE
#endif

!**** Boundary condition arrays:
!@var ZATMO,HLAKE Topography arrays: elevation (m), lake sill depth (m)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: ZATMO
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: HLAKE
!@var Fxx fraction of gridbox of type xx (land,ocean,...)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FLAND
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FOCEAN
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FLICE
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FLAKE0
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FEARTH0

!@var WFCS water field capacity of first ground layer (kg/m2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: WFCS

!@var P surface pressure (hecto-Pascals - PTOP)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: P

C**** Define surface types (mostly used for weighting diagnostics)
!@param NTYPE number of different surface types
      INTEGER, PARAMETER :: NTYPE=6   ! orig = 3
!@var FTYPE fractions of each surface type
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)   :: FTYPE

!@param ITxx indices of various types (used only when it matters)
      INTEGER, PARAMETER :: ITOCEAN=1, ITOICE=2, ITEARTH=3,
     *                      ITLANDI=4, ITLAKE=5, ITLKICE=6

C**** Variables specific for stratosphere and/or strat diagnostics
!@var DO_GWDRAG when true, prints Gravity Wave diagnostics
      LOGICAL :: DO_GWDRAG = .false.
!@var iDO_GWDRAG number if AIJ Gravity wave diagnostics
      INTEGER :: iDO_GWDRAG = 0

!@dbparam UOdrag parameter that decides whether ocean.ice velocities
!@+   feed into drag calculation in surface (default = 0)
      INTEGER :: UOdrag = 0

!@nlparam QCHECK TRUE for running diagnostic checks
      LOGICAL :: QCHECK = .FALSE.

!@var stop_on TRUE stops the model (set with "kill -15 PID)
      LOGICAL :: stop_on = .FALSE.

      END MODULE MODEL_COM

      module PARAM
      implicit none
      private
      save

      public :: sync_param

      interface sync_param
        module procedure sync_iparam, sync_rparam, sync_cparam
        module procedure sync_aiparam, sync_arparam, sync_acparam
      end interface
      
      contains

      subroutine sync_iparam( name, value )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(inout) :: value
      end subroutine sync_iparam

      subroutine sync_aiparam( name, value, np )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(in) :: np
      integer, intent(inout) :: value(np)
      end subroutine sync_aiparam

      subroutine sync_rparam( name, value )
      implicit none
      character*(*), intent(in) :: name
      real*8, intent(inout) :: value
      end subroutine sync_rparam

      subroutine sync_arparam( name, value, np )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(in) :: np
      real*8, intent(inout) :: value(np)
      end subroutine sync_arparam

      subroutine sync_cparam( name, value )
      implicit none
      character*(*), intent(in) :: name
      character*(*), intent(inout) :: value
      end subroutine sync_cparam

      subroutine sync_acparam( name, value, np )
      implicit none
      character*(*), intent(in) :: name
      integer, intent(in) :: np
      character*(*), intent(inout) :: value(np)
      end subroutine sync_acparam

      end module PARAM

      subroutine stop_model( message, retcode )
      !@sum Aborts the execution of the program. Passes an error message and
      !@+ a return code to the calling script. Should be used instead of STOP
        !USE ERROR_MOD, ONLY : ERROR_STOP
        implicit none
      !@var message an error message (reason to stop)
        character*(*), intent (in) :: message
      !@var retcode return code to be passed to the calling script
        integer, intent(in) :: retcode

        !call error_stop( message, 'giss_wrapper' )        
        WRITE(6,*) message, 'giss_wrapper'
        CALL FLUSH(6)
        CALL EXIT(99999)
        STOP

      end subroutine stop_model

!==============================================================================
! LTM: shared/CONST.f
!==============================================================================

      MODULE CONSTANT
!@sum  CONSTANT definitions for physical constants and useful numbers
!@auth G. Schmidt
!@ver  1.0
      IMPLICIT NONE
      SAVE
C**** Conventions: 'by' implies reciprocal, 'rt' implies square root

C**** Numerical constants

      real*8,parameter :: pi = 3.1415926535897932d0 !@param pi    pi
      real*8,parameter :: twopi = 2d0*pi           !@param twopi 2*pi
      real*8,parameter :: radian = pi/180d0        !@param radian pi/180
!@param zero,one 0 and 1 for occasional use as arguments
      real*8,parameter :: zero = 0d0, one=1d0
!@param rt2,byrt2   sqrt(2), 1/sqrt(2)
      real*8,parameter :: rt2 = 1.4142135623730950d0
      real*8,parameter :: byrt2 = 1./rt2
!@param rt3,byrt3   sqrt(3), 1/sqrt(3)
      real*8,parameter :: rt3 = 1.7320508075688772d0
      real*8,parameter :: byrt3 = 1./rt3
!@param rt12,byrt12   sqrt(12), 1/sqrt(12)
      real*8,parameter :: rt12 = 3.4641016151377546d0
      real*8,parameter :: byrt12 = 1./rt12
      real*8,parameter :: by3 =1./3d0  !@param by3  1/3
      real*8,parameter :: by6 =1./6d0  !@param by6  1/6
      real*8,parameter :: by9 =1./9d0  !@param by9  1/9
      real*8,parameter :: by12=1./12d0 !@param by12 1/12
!@param undef Missing value
      real*8,parameter :: undef=-1.d30
!@param teeny  small positive value used in num/(den+teeny) to avoid 0/0
      real*8,parameter :: teeny=1.d-30
      integer*8,parameter :: intNaN=-1  ! i.e. = Z'FFFFFFFFFFFFFFFF'
!@param NaN NaN
#if (defined COMPILER_PGI || defined COMPILER_NAG)
      real*8,parameter :: NaN=1d30
#else
      real*8,parameter :: NaN=transfer(intNaN,1.d0)
#endif

C**** Physical constants

!@param stbo Stefan-Boltzmann constant (W/m^2 K^4)
      real*8,parameter :: stbo =5.67051d-8 !current best estimate

c**** Latent heats:
c**** Note that for energy conservation the efective latent heat at any
c**** temperature must follow these formulae (assuming a reference
c**** temperature of 0 Celcius, and constant specific heats).
c**** If specific heats vary as a function of temperature, the extra
c**** term becomes an integral
c**** lhe(T) = lhe(0) + (shv-shw) T (in C)
c**** lhm(T) = lhm(0) + (shw-shi) T (in C)
c**** lhs(T) = lhs(0) + (shv-shi) T (in C)
!@param lhe   latent heat of evap at 0 C (2.5008d6 J/kg)
      real*8,parameter :: lhe = 2.5d6
!@param lhm   latent heat of melt at 0 C (334590 J/kg)
      real*8,parameter :: lhm = 3.34d5
!@param bylhm  1/lhm
      real*8,parameter :: bylhm = 1./lhm
!@param lhs  latent heat of sublimation at 0 C (J/kg)
      real*8,parameter :: lhs = lhe+lhm

!@param rhow density of pure water (1000 kg/m^3)
      real*8,parameter :: rhow = 1d3
!@param rhows density of average sea water (1030 kg/m^3)
      real*8,parameter :: rhows = 1030d0
!@param byrhows recip. density of average sea water (1/1030 m^3/kg)
      real*8,parameter :: byrhows = 1d0/rhows
!@param rhoi density of pure ice (916.6 kg/m^3)
      real*8,parameter :: rhoi = 916.6d0
!@param byrhoi 1/rhoi (m^3/kg)
      real*8,parameter :: byrhoi = 1d0/rhoi

!@param tf freezing point of water at 1 atm (273.16 K)
      real*8,parameter :: tf = 273.16d0
!@param bytf 1/tf (K^-1)
      real*8,parameter :: bytf = 1d0/tf

!@param shw heat capacity of water (at 20 C) (4185 J/kg C)
      real*8,parameter :: shw  = 4185.
!@param byshw 1/shw
      real*8,parameter :: byshw = 1d0/shw

!@param shi heat capacity of pure ice (at 0 C) (2060 J/kg C)
      real*8,parameter :: shi  = 2060.
!@param byshi 1/shi
      real*8,parameter :: byshi = 1d0/shi

!@param fraction of O2 in the atmosphere (0-1)
      real*8,parameter :: pO2 = 0.209476d0

c**** RGAS = R/M_A = 1000* 8.314510 J/mol K /28.9655 g/mol
c**** For values of CO2 much larger than present day (> 4x conc)
c**** the molecular weight of dry air M_A could change.
c**** Assume that M_O2 = 31.9988 and M_CO2 = 44.00995
c**** and current percentages 20.946% and 0.0350% (US Stand. Atm.)
c**** Assuming CO2 displaces other gases equally M_A=28.9602 + n*0.00527
c**** where n is multiple of present day CO2 conc (350 ppm)
c**** For 4xCO2  M_A = 28.9813  => rgas = 286.89
c**** For 10xCO2 M_A = 29.0129  => rgas = 286.58
!@param gasc  gas constant (8.314510 J/mol K)
      real*8,parameter :: gasc = 8.314510d0
!@param bygasc  1/gasc
      real*8,parameter :: bygasc = 1./gasc
!@param mair molecular weight of dry air (28.9655 g/mol)
      real*8,parameter :: mair = 28.9655d0
!@param rgas gas constant (287.05 J/K kg)
      real*8,parameter :: rgas = 1d3 * gasc / mair ! = 287.05...

!@param mwat molecular weight of water vapour
      real*8,parameter :: mwat = 18.015d0
!@param rvap  gas constant for water vapour (461.5 J/K kg)
c**** defined as R/M_W = 1000* 8.314510 J/mol K /18.015 g/mol
      real*8,parameter :: rvap = 1d3 * gasc / mwat ! = 461.5...

!@param mrat  mass ratio of air to water vapour (0.62197)
      real*8,parameter :: mrat = mwat/mair    ! = 0.62197....
!@param bymrat 1/mrat (1.6078)
      real*8,parameter :: bymrat = 1./mrat    ! = 1.6078....
!@param deltx coeff. of humidity in virtual temperature defn. (0.6078)
      real*8,parameter :: deltx = bymrat-1.   ! = 0.6078....

!@param srat ratio of specific heats at const. press. and vol. (=1.401)
      real*8,parameter :: srat = 1.401d0
!@param kapa ideal gas law exponent for dry air (.2862)
c**** kapa = (g-1)/g where g=1.401 = c_p/c_v
      real*8,parameter :: kapa = (srat - 1.)/srat  ! =.2862....
!@param bykapa,bykapap1,bykapap2 various useful reciprocals of kapa
      real*8,parameter :: bykapa = 1./kapa
      real*8,parameter :: bykapap1 = 1./(kapa+1.)
      real*8,parameter :: bykapap2 = 1./(kapa+2.)

!@param sha specific heat of dry air (const. pres.) (rgas/kapa J/kg C)
      real*8,parameter :: sha = rgas/kapa
!@param bysha 1/sha
      real*8,parameter :: bysha = 1./sha

!@param shv specific heat of water vapour (const. pres.) (J/kg C)
c**** shv is currently assumed to be zero to aid energy conservation in
c**** the atmosphere. Once the heat content associated with water
c**** vapour is included, this can be set to the standard value
c**** Literature values are 1911 (Arakawa), 1952 (Wallace and Hobbs)
c**** Smithsonian Met Tables = 4*rvap + delta = 1858--1869 ????
c     real*8,parameter :: shv = 4.*rvap  ????
      real*8,parameter :: shv = 0.

C**** air viscosity - temperature independent
!@var visc_air0 dynamic viscosity of air (kg/m s)
      real*8,parameter :: visc_air0 = 1.7d-5

!@var visc_air_kin0 kinematic viscosity of air (1 bar 15 deg C) (m^2/s)
      real*8,parameter :: visc_air_kin0 = 1.46d-5

!@var visc_wtr_kin kinematic viscosity of water (35 psu, 20 deg C) (m^2/s)
      real*8,parameter :: visc_wtr_kin = 1.05d-6

!@var avog Avogadro's constant (atmos/mole)
      real*8,parameter :: avog=6.023d23

C**** Astronomical constants

!@param sday  sec per day (s)
      real*8,parameter :: sday = 86400.
!@param syr  sec per year (s)
      real*8,parameter :: syr = sday*365.

!@param hrday  hours in a day (hrs)
      real*8,parameter :: hrday = sday/3600.

!@param omega earth's rotation rate (7.29 s^-1)
c      real*8,parameter :: omega = 7.2921151467d-5 ! NOVAS value
      real*8,parameter :: EDPERD = 1.
      real*8,parameter :: EDPERY = 365.
      real*8,parameter :: omega = TWOPI*(EDPERD+EDPERY)/
     *                            (EDPERD*EDPERY*SDAY)
!@param omega2 2*omega
      real*8,parameter :: omega2 = 2.*omega

!@param radius radius of the earth (6371000 m, IUGG)
      real*8,parameter :: radius = 6371000.
!@param areag surface area of the earth (m^2)
      real*8,parameter :: areag = 4.*pi*radius*radius

!@param grav gravitaional accelaration (9.80665 m/s^2)
c**** SI reference gravity (at 45 deg) = 9.80665
      real*8,parameter :: grav = 9.80665d0
!@param bygrav 1/grav
      real*8,parameter :: bygrav = 1d0/grav

C**** lapse rate related variables
!@param GAMD dry adiabatic lapse rate (=0.0098 K/m)
      real*8, parameter :: gamd = grav*kapa/rgas
!@param BMOIST moist adiabatic lapse rate (K/m)
      real*8, parameter :: bmoist = 0.0065d0
!@param BBYG moist adiabatic lapse rate divided by grav
      real*8, parameter :: bbyg = bmoist*bygrav
!@param GBYRB grav divided by rgas and bmoist
      real*8, parameter :: gbyrb = grav/(rgas*bmoist)

C**** Useful conversion factors

!@param kg2mb,mb2kg conversion from milli-bars to kg/m^2
      real*8,parameter :: kg2mb = 1d-2*grav, mb2kg = 1d2*bygrav
!@param kgpa2mm,mm2kgpa conversion from kg/m^2 water to mm
      real*8,parameter :: kgpa2mm = 1d0, mm2kgpa = 1d0

      CONTAINS

      real*8 function visc_air(T)
!@sum visc_air dynamic viscosity of air (function of T) (kg/m s)
!@auth Sutherland formula
      real*8, intent(in) :: T  ! temperature (K)
      real*8, parameter :: n0=1.827d-5, T0=291.15d0, C=120d0

      visc_air = n0*sqrt((T/T0)**3)*(T0+C)/(T+C)

      return
      end function

      real*8 function visc_air_kin(T)
!@sum visc_air_kin kinematic viscosity of air (function of T) (m2/s)
!@auth COARE formula - Andreas (1989) CRREL Rep. 89-11
      real*8, intent(in) :: T  ! temperature (K)
      real*8, parameter :: nu0=1.326d-5, a0=6.542d-3, b0=8.301d-6,
     *     c0=4.84d-9
      real*8 :: Tc  ! temperature in deg C

      Tc=T-tf
      visc_air_kin = nu0*(1.+Tc*(a0+Tc*(b0-c0*Tc)))   !m2/s

      return
      end function

      END MODULE CONSTANT

!==============================================================================
! LTM: ESMF_Interface/Hidden_mod.F90
!==============================================================================

      module Hidden_mod
      private
      public :: Hidden_type

      type Hidden_type
      integer :: mpi_comm
      integer :: numProcesses
      integer :: numAllProcesses
      integer :: mpi_tag
      
      logical :: hasSouthPole
      logical :: hasNorthPole
      logical :: hasEquator
      
      logical :: periodicBC
      integer :: log_unit
    
!     @var lookup_pet index of PET for a given J
      integer, dimension(:), pointer :: lookup_pet
      
      end type Hidden_type
      
      end module Hidden_mod

!==============================================================================
! LTM: ESMF_Interface/dd2d_utils.f
!==============================================================================

      module dd2d_utils
!@auth M. Kelley
!@ver 1.0
!@sum  dd2d_utils provides communication procedures for domain
!@+    decomposition using MPI.  Currently, 2D cubed-sphere
!@+    and 1D lat-lon layouts are enabled.  These routines
!@+    were developed to meet the requirements of modelE
!@+    not met by other packages.
!@+
!@usage Access is through the following types/interfaces.
!@+
!@+     Below, grid is an instance of the dist_grid derived type,
!@+     and the optional argument jdim indicates which array index
!@+     corresponds to the last distributed dimension (usually, j).
!@+     If jdim is not specified, it is assumed to be 2.
!@+
!@+     dist_grid: a derived type containing domain decomp info,
!@+     initialized by calling subroutine
!@+     init_dist_grid(
!@+     &     npx,npy,ntiles, ! x,y size of tile and # of tiles
!@+     &     is,ie,js,je,    ! start, end of local i,j domain
!@+     &     isd,ied,jsd,jed, ! start, end of local data domain (not used)
!@+     &     grid)
!@+     ntiles = 1 for a latlon layout, and 6 for a cubed sphere layout
!@+     NOTE: all arguments save for grid are inputs - this routine does
!@+           NOT decide the domain decomposition is,ie,js,je!
!@+
!@+     halo_update(grid,array,jdim)
!@+        fills the halo cells of array with the appropriate values
!@+        from neighboring PEs.  The depth of the halo is inferred
!@+        by comparing the dimension sizes of array to the start,end
!@+        info in grid.  For a latlon grid, no halos are performed in
!@+        the I direction.
!@+        For the cubed-sphere layout, the current implementation
!@+        requires npx=npy and an N x N layout of PEs on each face.
!@+
!@+     pack_data(grid,local_array,global_array,jdim)
!@+        gathers local arrays into a global array on root
!@+     unpack_data(grid,global_array,local_array,jdim)
!@+        scatters a global array on root into local arrays
!@+     If global_array has the same rank as local_array, the
!@+        gather/scatter operations are over the local tile.
!@+        Each tile has a root PE with a valid copy of global_array.
!@+     If the rank of global_array = 1 + rank of local_array, the
!@+        gather/scatter operations are over all tiles.  The last
!@+        index of global_array is the tile index in this case.
!@+
!@+     globalsum(grid,local_arr,arrsum,all)
!@+        calculates the global sum of distributed array local_arr
!@+        and stores the result in arrsum; the result is
!@+        independent of the number of PEs.  If all=.true., all
!@+        PEs receive the result. Currently only implemented for
!@+        real*8 local_arr(i,j).
!@+
!@+     pack_row(grid,local_array,row_array_1d,jdim)
!@+     unpack_row(grid,row_array_1d,local_array,jdim)
!@+        These routines are like pack_data/unpack_data, but only
!@+        gather/scatter to/from one PE for each "row" of PEs.
!@+        For now, the row array with the contents of all
!@+        PEs in the row is 1D since it is only used for I/O.
!@+        Dimensioned versions of row_array will be added soon.
!@+
!@+

c
c Implementation notes for halo_update
c
c For simplicity, halo_update is programmed using mpi_sendrecv
c and predefined MPI communicators for "east-west" and
c "north-south" directions. For the latlon grid, "north-south"
c corresponds to geographic N-S.
c For the cubed sphere 2D decomposition, these communicators
c trace periodic paths.  Given the following example 2x2 PE
c layout on each face,
c
c              18 19 22 23
c              16 17 20 21
c        10 11 14 15
c         8  9 12 13
c   2  3  6  7
c   0  1  4  5
c
c  the PE groupings for the first periodic direction are
c
c  0,2,3,6,8,10,11,14,16,18,19,22
c  1,4,5,7,9,12,13,15,17,20,21,23
c
c  and for the second periodic direction they are
c
c  2,10,18
c  0,1,3,8,9,11,16,17,19
c  4,6,7,12,14,15,20,22,23
c  5,13,21
c
c Most corner halos can be obtained by adjusting halo bounds
c during the sendrecv along the second periodic direction.
c However, for PEs on or next to the UL-LR diagonal of each cubed
c sphere tile, one extra call to mpi_sendrecv is needed.
c

      use Hidden_mod, only: Hidden_type

      implicit none
      private
      save

      public :: init_dist_grid
      public :: dist_grid
      type dist_grid
        type (Hidden_type) :: private
        integer :: npx ! number of i cells
        integer :: npy ! number of j cells
        integer :: is  ! first i index of computational domain
        integer :: ie  ! last  i index of computational domain
        integer :: js  ! first j index of computational domain
        integer :: je  ! last  j index of computational domain
        integer :: isd ! first i index of data domain
        integer :: ied ! last  i index of data domain
        integer :: jsd ! first j index of data domain
        integer :: jed ! last  j index of data domain
        logical :: am_i_globalroot ! am I the root for all processes
        integer :: ntiles ! number of tiles globally

c the following is for use in modelE only:
c
        logical :: have_domain ! whether this PE has any of the domain
        integer :: mpi_comm    ! MPI communicator for PEs having this domain
        integer :: npes_comm   ! number of PEs having this domain
        integer :: mpi_tag     ! for MPI book-keeping

        INTEGER :: NPES_USED
         ! Parameters for Global domain
        INTEGER :: IM_WORLD     ! Number of Longitudes
        INTEGER :: JM_WORLD     ! Number of latitudes
         ! Parameters for local domain
        INTEGER :: I_STRT       ! Begin local domain longitude index
        INTEGER :: I_STOP       ! End   local domain longitude index
        INTEGER :: J_STRT       ! Begin local domain latitude  index
        INTEGER :: J_STOP       ! End   local domain latitude  index
        INTEGER :: J_STRT_SKP   ! Begin local domain exclusive of S pole
        INTEGER :: J_STOP_SKP   ! End   local domain exclusive of N pole
        INTEGER :: ni_loc       ! for transpose
         ! Parameters for halo of local domain
        INTEGER :: I_STRT_HALO  ! Begin halo longitude index
        INTEGER :: I_STOP_HALO  ! End   halo longitude index
        INTEGER :: J_STRT_HALO  ! Begin halo latitude  index
        INTEGER :: J_STOP_HALO  ! End   halo latitude  index
         ! Parameters for staggered "B" grid
         ! Note that global staggered grid begins at "2".
        INTEGER :: J_STRT_STGR  ! Begin local staggered domain
        INTEGER :: J_STOP_STGR  ! End   local staggered domain
         ! Controls for special cases

        INTEGER, DIMENSION(:), POINTER :: DJ_MAP
        INTEGER :: DJ
        INTEGER :: log_unit     ! for debugging
         !@var lookup_pet index of PET for a given J
        INTEGER, DIMENSION(:), POINTER :: lookup_pet
        LOGICAL :: BC_PERIODIC

      end type dist_grid

      contains

      subroutine init_dist_grid(
     &     npx,npy,ntiles,
     &     is,ie,js,je,
     &     isd,ied,jsd,jed,
     &     grid)
      integer, intent(in) :: npx,npy,ntiles,is,ie,js,je,isd,ied,jsd,jed
      type(dist_grid) :: grid
      integer :: i,itile,ierr,iproc,n,ihem,ihem_sv,modrank,midp1
     &     ,group_world,group_mytile,group_intertile
     &     ,group_halo,halo_comm,group_row,row_comm,xpos,ypos
      integer, parameter :: ntiles_max=6
      integer :: nproc_comm
      integer, dimension(ntiles_max) ::
     &      tile_comms,tile_root_procs
      integer, dimension(:), allocatable :: proclist
      logical :: swap_ne,swap_sw

      grid%am_i_globalroot = .true. ! initialize to serial-mode default
      grid%npx = npx
      grid%npy = npy
      grid%ntiles = ntiles
      grid%is = is
      grid%ie = ie
      grid%js = js
      grid%je = je
      grid%isd = isd
      grid%ied = ied
      grid%jsd = jsd
      grid%jed = jed

      return
      end subroutine init_dist_grid

      end module dd2d_utils
      
!==============================================================================
! LTM: ESMF_Interface/dist_grid_mod.F90
!==============================================================================

      MODULE dist_grid_mod
!@sum  This module encapsulates lat-lon decomposition information for the
!@+    message passing implementation.
!@auth Tom Clune GSFC/SSSO/610.3

      use Hidden_mod
      use dd2d_utils, only : dist_grid, init_dist_grid

      IMPLICIT NONE

      SAVE
      PRIVATE ! Except for

!@var DIST_GRID derived type to provide decomposition info
!@+   public components are used to minimize overhead for accessing
!@+   routine components
      PUBLIC :: DIST_GRID
!@var  grid Default decomposition; globally accessible for convenience.
      PUBLIC :: grid
!@var INIT_APP set some parameters
      PUBLIC :: INIT_GRID
!@var GLOBALMAX determine max value across pes
      PUBLIC :: GLOBALMAX
!@var GET - extracts bounds information from DIST_GRID object
      PUBLIC :: GET
      
      INTEGER, PARAMETER :: HALO_WIDTH = 1
      integer ::  root
      
      type (DIST_GRID) :: GRID
!not used      TYPE (DIST_GRID) :: GRID_TRANS
 
      public :: haveLatitude
       
! Remaining variables are private to the module.
 
!@var NPES_WORLD number of total processes
      INTEGER :: NPES_WORLD
!@var NP_LON number of azimuthal processes.
      INTEGER :: NP_LON
!@var NP_LAT number of meridional     processes.
      INTEGER :: NP_LAT
!@var MY_PET index of _this_ PET (analagous to MPI rank)
      INTEGER :: my_pet
!@var RANK_LON index of _this_ process in azimuthal set.
      INTEGER :: RANK_LON
!@var RANK_LAT_RANK index of _this_ process in meridional set.
      INTEGER :: RANK_LAT
 
      INTEGER, PUBLIC :: CHECKSUM_UNIT
 
! temporary public during refactoring
      public :: isPeriodic
      public :: my_pet, npes_world
      public :: am_i_root

!@var GLOBALMAX Generic wrapper for Real/integer
      INTERFACE GLOBALMAX
        MODULE PROCEDURE GLOBALMAX_R
        MODULE PROCEDURE GLOBALMAX_I
        MODULE PROCEDURE GLOBALMAX_I_1D
      END INTERFACE
 
! Direction bits
      public :: NORTH, SOUTH
      Integer, Parameter :: NORTH = 2**0
      Integer, Parameter :: SOUTH = 2**1
 
      public :: hasSouthPole
      public :: hasNorthPole
      public :: hasPeriodicBC
      
      integer, parameter :: maxStrLen = 40
      public :: maxStrLen

      CONTAINS

      FUNCTION AM_I_ROOT() RESULT(AMI)
      LOGICAL :: AMI
      AMI = .TRUE.
      END FUNCTION

      SUBROUTINE INIT_GRID(grd_dum,IM,JM,LM,width,J_SCM,bc_periodic, 
     &                     CREATE_CAP,npes_max,excess_on_pe0)
!@sum  This routine initializes the ESMF grid as well as data
!@+    structures associated with modelE's domain decomposition
!@auth Tom Clune GSFC/SSSO/610.3
      !USE FILEMANAGER, Only : openunit
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(INOUT) :: grd_dum
      INTEGER, INTENT(IN) :: IM, JM,LM
      INTEGER, OPTIONAL, INTENT(IN) :: J_SCM ! single column model
      INTEGER, OPTIONAL :: width
      LOGICAL, OPTIONAL, INTENT(IN) :: bc_periodic
      LOGICAL, OPTIONAL, INTENT(IN) :: CREATE_CAP
      INTEGER, OPTIONAL, INTENT(IN) :: npes_max
      LOGICAL, OPTIONAL, INTENT(IN) :: excess_on_pe0
      integer, parameter :: numDims=2
      integer, dimension(numDims) :: grid_size
      integer             :: rc
      INTEGER :: J_EQUATOR
      INTEGER :: I0_DUM, I1_DUM
      INTEGER :: J0_DUM, J1_DUM
      INTEGER :: width_
      INTEGER :: pet
      INTEGER :: NTILES
      INTEGER :: AIbounds(4)

      INTEGER :: p
      integer :: npes_used
      integer, dimension(:), allocatable :: pelist
      integer :: group_world,group_used,ierr
      integer :: newCommunicator

      logical :: useCubedSphere

      useCubedSphere = .false.

      if (useCubedSphere) then
        grid_size = (/IM, JM*6/)
!!$$        grd_dum%localSubdomain = newDomain(IM, 6*JM)
      else
        grid_size = (/IM, JM/)
!!$$        grd_dum%localSubdomain = newDomain(IM, JM)
      end if

      grd_dum%IM_WORLD      = IM
      grd_dum%JM_WORLD      = JM

      npes_used = 1
      grd_dum%npes_used = npes_used

      I0_DUM = 1
      I1_DUM = IM
      J0_DUM = 1
      J1_DUM = JM

      if (present(J_SCM)) then
         J0_DUM = J_SCM
         J1_DUM = J_SCM
      end if


      width_ = HALO_WIDTH
      If (Present(width)) width_=width

! Wrapped grid
      grd_dum%I_STRT        = I0_DUM
      grd_dum%I_STOP        = I1_DUM
      grd_dum%I_STRT_HALO   = MAX( 1, I0_DUM-width_)
      grd_dum%I_STOP_HALO   = MIN(IM, I1_DUM+width_)
      grd_dum%ni_loc = (RANK_LAT+1)*IM/NPES_used - RANK_LAT*IM/NPES_used

      grd_dum%j_strt        = J0_DUM
      grd_dum%J_STOP        = J1_DUM
      grd_dum%HAVE_DOMAIN   = J0_DUM <= JM

      grd_dum%J_STRT_SKP = max (   2, J0_DUM)
      grd_dum%J_STOP_SKP = min (JM-1, J1_DUM)

! I guess we don't need HALO in SCM mode...
!grd_dum%J_STRT_HALO = MAX(1,  grd_dum % J_STRT - 1)
!grd_dum%J_STOP_HALO = MIN(JM, grd_dum % J_STOP + 1)
      grd_dum%J_STRT_HALO = MAX(1,  grd_dum % J_STRT)
      grd_dum%J_STOP_HALO = MIN(JM, grd_dum % J_STOP)


      grd_dum%J_STRT_STGR   = max(2,J0_DUM)
      grd_dum%J_STOP_STGR   = J1_DUM

      grd_dum%private%hasSouthPole = J0_DUM == 1  !(RANK_LAT == 0)
      grd_dum%private%hasNorthPole = J1_DUM == JM .and. J0_DUM <= JM !(RANK_LAT == NP_LAT - 1) &

      J_EQUATOR = JM/2
      grd_dum%private%hasEquator =  (J0_DUM <= J_EQUATOR) .AND. 
     &                              (J1_DUM >= J_EQUATOR)


! need to initialize the dd2d version of dist_grid for I/O
           call init_dist_grid( 
     &     grd_dum%IM_WORLD,grd_dum%JM_WORLD,1,  
     &     grd_dum%I_STRT,grd_dum%I_STOP, 
     &     grd_dum%j_strt,grd_dum%J_STOP, 
     &     grd_dum%I_STRT_HALO,grd_dum%I_STOP_HALO, 
     &     grd_dum%J_STRT_HALO,grd_dum%J_STOP_HALO,grd_dum)


      if (present(J_SCM)) then
! assume J_SCM is in "general position"
        grd_dum%private%hasSouthPole = .false.
        grd_dum%private%hasNorthPole = .false.
        grd_dum%private%hasEquator    = .false.
      endif

      grd_dum % private%PERIODICBC = isPeriodic(bc_periodic)

! assumption: decomposition along "east-west" direction
! not used for lat-lon grids
      if(grd_dum % private%PERIODICBC) then
        grd_dum%private%hasSouthPole = .false.
        grd_dum%private%hasNorthPole = .false.
        grd_dum%private%hasEquator    = .false.
      endif

! set lookup table PET(J)
      Allocate(grd_dum%private%lookup_pet(1:JM))
      grd_dum%private%lookup_pet(:) = 0

      END SUBROUTINE INIT_GRID

      SUBROUTINE GET(grd_dum, I_STRT, I_STOP, 
     &                        I_STRT_HALO, I_STOP_HALO, 
     &                        J_STRT, J_STOP, J_STRT_HALO, J_STOP_HALO, 
     &                        J_STRT_SKP, J_STOP_SKP, 
     &                        J_STRT_STGR, J_STOP_STGR, 
     &                        have_south_pole, have_north_pole)
!@sum  This routine provides general access to various values that characterize the local
!@+ domain bounds on this process.
!@auth Tom Clune GSFC/SSSO/610.3
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      INTEGER, OPTIONAL :: I_STRT, I_STOP
      INTEGER, OPTIONAL :: I_STRT_HALO, I_STOP_HALO
      INTEGER, OPTIONAL :: J_STRT, J_STOP
      INTEGER, OPTIONAL :: J_STRT_HALO, J_STOP_HALO
      INTEGER, OPTIONAL :: J_STRT_SKP, J_STOP_SKP
      INTEGER, OPTIONAL :: J_STRT_STGR, J_STOP_STGR
      LOGICAL, OPTIONAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      IF (PRESENT(I_STRT)) I_STRT = grd_dum%I_STRT
      IF (PRESENT(I_STOP)) I_STOP = grd_dum%I_STOP

      IF (PRESENT(I_STRT_HALO)) I_STRT_HALO = grd_dum%I_STRT_HALO
      IF (PRESENT(I_STOP_HALO)) I_STOP_HALO = grd_dum%I_STOP_HALO

      IF (PRESENT(J_STRT)) J_STRT = grd_dum%j_strt
      IF (PRESENT(J_STOP)) J_STOP = grd_dum%J_STOP

      IF (PRESENT(J_STRT_HALO)) J_STRT_HALO = grd_dum%J_STRT_HALO
      IF (PRESENT(J_STOP_HALO)) J_STOP_HALO = grd_dum%J_STOP_HALO

      IF (PRESENT(J_STRT_SKP)) J_STRT_SKP = grd_dum%J_STRT_SKP
      IF (PRESENT(J_STOP_SKP)) J_STOP_SKP = grd_dum%J_STOP_SKP

      IF (PRESENT(J_STRT_STGR)) J_STRT_STGR = grd_dum%J_STRT_STGR
      IF (PRESENT(J_STOP_STGR)) J_STOP_STGR = grd_dum%J_STOP_STGR

      IF (PRESENT(HAVE_SOUTH_POLE)) 
     &     HAVE_SOUTH_POLE= grd_dum%private%hasSouthPole
      IF (PRESENT(HAVE_NORTH_POLE)) 
     &     HAVE_NORTH_POLE= grd_dum%private%hasNorthPole

      END SUBROUTINE GET

      logical function haveLatitude(grd_dum, j)
!@sum returns true if latitude j is in the local domain
!@auth Tom Clune GSFC/SSSO/610.3
      type (DIST_GRID), intent(in) :: grd_dum
      integer, intent(in) :: j

      haveLatitude = (j >= grd_dum%j_strt .and. j <= grd_dum%J_STOP)

      end function haveLatitude

      logical function isPeriodic(override)
!@sum Helper function to handle optional arguments related to 
!@+ periodic boundaries
!@auth Tom Clune GSFC/SSSO/610.3
        logical, optional, intent(in) :: override

        isPeriodic = .false.
        if (present(override)) isPeriodic = override

      end function isPeriodic

      logical function hasSouthPole(this)
!@sum  returns true if GRID has a south pole
!@auth Tom Clune GSFC/SSSO/610.3
        type (dist_grid), intent(in) :: this
        hasSouthPole = this%private%hasSouthPole
      end function hasSouthPole

      logical function hasNorthPole(this)
!@sum  returns true if GRID has a north pole
!@auth Tom Clune GSFC/SSSO/610.3
        type (dist_grid), intent(in) :: this
        hasNorthPole = this%private%hasNorthPole
      end function hasNorthPole

      logical function hasPeriodicBC(this)
!@sum  returns true if GRID has periodic boundary conditions
!@auth Tom Clune GSFC/SSSO/610.3
        type (dist_grid), intent(in) :: this
        hasPeriodicBC = this%private%periodicBC
      end function hasPeriodicBC

      SUBROUTINE GLOBALMAX_R(grd_dum, val, val_max)
!@sum determine max value of real variable across processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      REAL*8,            INTENT(IN)  :: val
      REAL*8,            INTENT(OUT) :: val_max

      INTEGER  :: ierr

      val_max = val

      END SUBROUTINE

      SUBROUTINE GLOBALMAX_I(grd_dum, val, val_max)
!@sum  determine max value of integer variable across processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      INTEGER,            INTENT(IN)  :: val
      INTEGER,            INTENT(OUT) :: val_max

      INTEGER  :: ierr

      val_max = val

      END SUBROUTINE

      SUBROUTINE GLOBALMAX_I_1D(grd_dum, val, val_max)
!@sum determine max value of integer 1-D array across processes
!@auth Tom Clune GSFC/SSSO/610.3
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      INTEGER,            INTENT(IN)  :: val(:)
      INTEGER,            INTENT(OUT) :: val_max(:)

      INTEGER  :: n,ierr

      val_max(:) = val(:)

      END SUBROUTINE

      end module dist_grid_mod

!==============================================================================
! LTM: ESMF_Interface/Halo_mod.F90
!==============================================================================

#define FILL(N) IAND(USABLE_FROM,N)==N
  
      module Halo_mod
!@sum  Collection of routines to perform halo updates.  
!@auth Tom Clune GSFC/SSSO/610.3
      use dist_grid_mod
      implicit none
      private
      
!@var halo_update Update data in halo for local domain using data from
!@+   neighbouring processes
      public :: halo_update     ! Communicate overlapping portions of subdomains
      public :: halo_updateJ    ! jx
      public :: halo_update_column ! K, I, J
      public :: halo_update_block ! K, L, I, J
      public :: halo_update_mask ! K, L, I, J
      
      public :: NORTH, SOUTH
      integer, parameter :: ALL = NORTH + SOUTH ! no east/west for now
      
!@var halo_update generic wrapper for 2D and 3D routines
      interface halo_update
      module procedure halo_update_1d ! j
      module procedure halo_update_2d ! i,j
      module procedure halo_update_3d ! i,j,k
      module procedure halo_update_2dint ! i,j,k
      end interface
      
      interface halo_updatej
      module procedure halo_updatej_2d
      end interface
      
      interface halo_update_column
      module procedure halo_update_column_2d ! m,j
      module procedure halo_update_column_3d ! m,i,j
      module procedure int_halo_update_column_3d ! m,i,j
      module procedure halo_update_column_4d ! m,i,j,k
      module procedure halo_update_column_7d ! m1,m2,m3,m4,i,j
      end interface
      
      interface halo_update_block
      module procedure halo_update_block_4d ! k,l,i,j
      end interface
            
      contains
      
      subroutine halo_update_1d(grid, arr, from)
      type (dist_grid),   intent(in)   :: grid
      real*8,   intent(inout) :: arr(grid%j_strt_halo:)
      integer, optional, intent(in)    :: from
      
      call sendrecv(grid, arr, shape(arr), 1, from
     &     ,hasperiodicbc(grid))
      
      end subroutine halo_update_1d
      
      subroutine halo_update_2d(grid, arr, from)
      type (dist_grid),   intent(in)   :: grid
      real*8,   intent(inout) :: 
     &     arr(grid%i_strt_halo:,grid%j_strt_halo:)
      integer, optional, intent(in)    :: from
      
      call sendrecv(grid, arr, shape(arr), 2, from, hasperiodicbc(grid))
      end subroutine halo_update_2d
      
      subroutine halo_update_2dint(grid, arr, from)
      type (dist_grid),   intent(in)    :: grid
      integer,   intent(inout) :: 
     &     arr(grid%i_strt_halo:,grid%j_strt_halo:)
      integer, optional, intent(in)    :: from
      
      call sendrecv_int
     &     (grid, arr, shape(arr), 2, from, hasperiodicbc(grid))
      end subroutine halo_update_2dint
      
      subroutine halo_updatej_2d(grid, arr, from)
      type (dist_grid),   intent(in)   :: grid
      real*8,   intent(inout) :: arr(grid%j_strt_halo:,:)
      integer, optional, intent(in)    :: from
      
      call sendrecv(grid, arr, shape(arr), 1, from, hasperiodicbc(grid))
      
      end subroutine halo_updatej_2d
      
      subroutine halo_update_3d(grid, arr, from, jdim)
      type (dist_grid),   intent(in)    :: grid
      real*8,   intent(inout) :: arr(:,:,:)
      integer, optional, intent(in)    :: from
      integer, optional, intent(in)    :: jdim
      
      integer :: jd
      
      if(present(jdim)) then
         jd = jdim
      else
         jd = 2
      endif
      call sendrecv
     &     (grid, arr, shape(arr), jd, from, hasperiodicbc(grid))
      
      end subroutine halo_update_3d
      
      subroutine halo_update_mask(grid, sbufs, sbufn, rbufs, rbufn)
      type (dist_grid),   intent(in)    :: grid
      real*8,   intent(in) :: sbufn(:)
      real*8,   intent(in) :: sbufs(:)
      real*8,   intent(out) :: rbufn(:)
      real*8,   intent(out) :: rbufs(:)
! NOOP
      end subroutine halo_update_mask
      
      subroutine halo_update_column_2d(grid, arr, from)
      type (dist_grid),   intent(in)   :: grid
      real*8,   intent(inout) :: arr(:,grid%j_strt_halo:)
      integer, optional, intent(in)    :: from
      
      call sendrecv(grid, arr, shape(arr), 2, from, hasperiodicbc(grid))
      
      end subroutine halo_update_column_2d
      
      subroutine halo_update_column_3d(grid, arr, from)
      type (dist_grid),   intent(in)    :: grid
      real*8,   intent(inout) :: 
     &     arr(:,grid%i_strt_halo:,grid%j_strt_halo:)
      integer, optional, intent(in)    :: from
      
      call sendrecv(grid, arr, shape(arr), 3, from, hasperiodicbc(grid))
      
      end subroutine halo_update_column_3d
      
      subroutine int_halo_update_column_3d(grid, arr, from)
      type (dist_grid),   intent(in)   :: grid
      integer,  intent(inout) ::
     &     arr(:,grid%i_strt_halo:,grid%j_strt_halo:)
      integer, optional, intent(in)    :: from
      
      real*8 :: foo(size(arr,1),size(arr,2),size(arr,3))
      
! cheat
      foo = arr
      call halo_update_column(grid, foo, from)
      arr = nint(foo)
      
      end subroutine int_halo_update_column_3d
      
      subroutine halo_update_column_4d(grid, arr, from)
      type (dist_grid),   intent(in)    :: grid
      real*8,   intent(inout) ::
     &     arr(:,grid%i_strt_halo:,grid%j_strt_halo:,:)
      integer, optional, intent(in)    :: from
      
      call sendrecv(grid, arr, shape(arr), 3, from)
      
      end subroutine halo_update_column_4d
      
      subroutine halo_update_block_4d(grid, arr, from)
      type (dist_grid),   intent(in)   :: grid
      real*8,   intent(inout) ::
     &     arr(:,:,grid%i_strt_halo:,grid%j_strt_halo:)
      integer, optional, intent(in)    :: from
      
      call sendrecv(grid, arr, shape(arr), 4, from)
      
      end subroutine halo_update_block_4d
      
      subroutine halo_update_column_7d(grid, arr, from)
      type (dist_grid),   intent(in)   :: grid
      real*8, intent(inout) ::
     &     arr(:,:,:,:,grid%i_strt_halo:,grid%j_strt_halo:)
      integer, optional, intent(in)    :: from
      
      call sendrecv(grid, arr, shape(arr), 6, from, hasperiodicbc(grid))
      
      end subroutine halo_update_column_7d
      
      subroutine sendrecv
     &     (grid, arr, shp, dist_idx, from, bc_periodic_)
      Type (Dist_Grid) :: grid
      Real(Kind=8) :: arr(*)
      Integer :: shp(:)
      Integer :: dist_idx
      Integer, optional :: from
      Logical, optional :: bc_periodic_
      end subroutine sendrecv
      
      subroutine sendrecv_int
     &     (grid, arr, shp, dist_idx, from, bc_periodic_)
      Type (Dist_Grid) :: grid
      Integer :: arr(*)
      Integer :: shp(:)
      Integer :: dist_idx
      Integer, optional :: from
      Logical, optional :: bc_periodic_
      end subroutine sendrecv_int
      
      ! Helper function to handle optional arguments related to halo directions
      integer function usableFrom(fromDirection)
      integer, optional, intent(in) :: fromDirection
      usableFrom = ALL
      if (present(fromDirection)) usableFrom = fromDirection
      end function usableFrom
      
      end module Halo_mod

#include "rundeck_opts.h"

!==============================================================================
! LTM: ESMF_Interface/GatherScatter_mod.F90
!==============================================================================

      module GatherScatter_mod
      !@sum Collection of routines to support MPI gather-scatter functionality
      !@+   for various data types and array sizes.
      !@auth Tom Clune GSFC/SSSO/610.3
        use dist_grid_mod
        implicit none
        private
      
        public :: pack_data
        public :: pack_block
        public :: pack_dataj
        public :: pack_j
      
        public :: unpack_data
        public :: unpack_block
        public :: unpack_dataj
        public :: unpack_j
      
        public :: pack_column
        public :: unpack_column
      
        public :: gatherReal8
        public :: scatterReal8
      
      !@var PACK Generic routine to pack a global array
      !@+   with the data from the corresponding distributed array.
        interface pack_data
          module procedure pack_1d       ! (i)
          module procedure ipack_1d      ! (i)
          module procedure pack_2d       ! (i,j)
          module procedure ipack_2d      ! (i,j)
          module procedure lpack_2d      ! (i,j)
          module procedure pack_3d       ! (i,j,l)
          module procedure ipack_3d      ! (i,j,l)
          module procedure pack_4d       ! (i,j,l,m)
          module procedure pack_5d       ! (i,j,l,m,n)
        end interface
      
        interface pack_dataj
          module procedure packj_2d     ! (j,k)
          module procedure ipackj_2d    ! (j,k)
          module procedure packj_3d     ! (j,k,l)
          module procedure packj_4d     ! (j,k,l,m)
        end interface
      
      !@var UNPACK Generic routine to unpack into a distributed
      !@+   array the data from the corresponding global array.
        interface unpack_data
          module procedure unpack_1d      ! (i)
          module procedure unpack_2d      ! (i,j)
          module procedure iunpack_2d     ! (i,j)
          module procedure lunpack_2d     ! (i,j)
          module procedure unpack_3d      ! (i,j,l)
          module procedure iunpack_3d     ! (i,j,l)
          module procedure unpack_4d      ! (i,j,l,m)
          module procedure unpack_5d      ! (i,j,l,m,n)
        end interface
      
        interface unpack_dataj
          module procedure unpackj_2d     ! (j,k)
          module procedure unpackj_3d     ! (j,k,l)
          module procedure unpackj_4d     ! (j,k,l,m)
        end interface
      
      !@var PACK_COLUMN Generic routine to pack  a global array
      !@+   with the data from the corresponding distributed array.
        interface pack_column
          module procedure pack_column_1d  ! (k,  j  )
          module procedure pack_column_2d  ! (k,i,j  )
          module procedure pack_column_i2d ! (k,i,j  )
          module procedure pack_column_3d  ! (k,i,j,l)
        end interface
      
      !@var UNPACK_COLUMN Generic routine to unpack into a distributed
      !@+   array the data from the corresponding global array.
        interface unpack_column
          module procedure unpack_column_1d  ! (k,  j  )
          module procedure unpack_column_2d  ! (k,i,j  )
          module procedure iunpack_column_2d ! (k,i,j  )
          module procedure unpack_column_3d  ! (k,i,j,l)
        end interface
      
      !@var PACK_BLOCK  Generic routine to pack  a global array
      !@+   with the data from the corresponding distributed array.
        interface pack_block
          module procedure ipack_block_2d    ! (k,l,i,j  )
          module procedure  pack_block_2d    ! (k,l,i,j  )
          module procedure  pack_block_3d    ! (k,l,m,i,j)
        end interface
      
      !@var UNPACK_BLOCK  Generic routine to unpack into a distributed
      !@+   array the data from the corresponding global array.
        interface unpack_block
          module procedure iunpack_block_2d    ! (k,l,i,j  )
          module procedure  unpack_block_2d    ! (k,l,i,j  )
          module procedure  unpack_block_3d    ! (k,l,m,i,j)
        end interface
      
      !@var PACK_J Generic routine to pack  a global array
      !@+   with the data from the corresponding distributed array.
        interface pack_j
          module procedure unpackj_2d     ! (j,k)
          module procedure unpackj_3d     ! (j,k,l)
          module procedure unpackj_4d     ! (j,k,l,m)
        end interface
      !@var UNPACK Generic routine to unpack into a distributed
      !@+   array the data from the corresponding global array.
        interface unpack_j
          module procedure unpackj_2d     ! (j,k)
          module procedure unpackj_3d     ! (j,k,l)
          module procedure unpackj_4d     ! (j,k,l,m)
        end interface
      
        interface localCopy
          module procedure localCopyReal8
          module procedure localCopyInteger
          module Procedure localCopyLogical
        end interface
      
      contains
      
        subroutine 
     &     gatherReal8(grid, arr_loc, arr_glob, shp, dist_idx, all)
          type (dist_grid) :: grid
          real(kind=8) :: arr_loc(*)
          real(kind=8) :: arr_glob(*)
          integer, intent(in) :: shp(:)
          integer, intent(in) :: dist_idx
          logical, optional, intent(in) :: all
      
          integer :: nn
      
          nn = product(shp)
          arr_glob(1:nn) = arr_loc(1:nn)
      
        end subroutine gatherReal8
      
        subroutine 
     &     gatherInteger(grid, arr_loc, arr_glob, shp, dist_idx, all)
          type (dist_grid) :: grid
          integer :: arr_loc(*)
          integer :: arr_glob(*)
          integer, intent(in) :: shp(:)
          integer, intent(in) :: dist_idx
          logical, optional, intent(in) :: all
      
          integer :: nn
      
          nn = product(shp)
          arr_glob(1:nn) = arr_loc(1:nn)
      
        end subroutine gatherInteger
      
        subroutine 
     &     gatherLogical(grid, arr_loc, arr_glob, shp, dist_idx, all)
          type (dist_grid) :: grid
          logical :: arr_loc(*)
          logical :: arr_glob(*)
          integer, intent(in) :: shp(:)
          integer, intent(in) :: dist_idx
          logical, optional, intent(in) :: all
      
          integer :: nn
      
          nn = product(shp)
          arr_glob(1:nn) = arr_loc(1:nn)
      
        end subroutine gatherLogical
        
        subroutine 
     &      scatterReal8(grid, arr_glob, arr_loc, shp, dist_idx, local)
          type (dist_grid) :: grid
          real(kind=8) :: arr_loc(*)
          real(kind=8) :: arr_glob(*)
          integer, intent(in) :: shp(:)
          integer, intent(in) :: dist_idx
          logical, optional, intent(in) :: local ! unused in serial mode
      
          call localCopyReal8 
     &     (grid, arr_glob, arr_loc,  
     &     getGlobalShape(grid, shp, dist_idx))
      
        end subroutine scatterReal8
      
        subroutine 
     &      scatterInteger(grid, arr_glob, arr_loc, 
     &      shp, dist_idx, local)
          type (dist_grid) :: grid
          integer :: arr_loc(*)
          integer :: arr_glob(*)
          integer, intent(in) :: shp(:)
          integer, intent(in) :: dist_idx
          logical, optional, intent(in) :: local ! unused in serial mode
      
          call localCopyInteger 
     &     (grid, arr_glob, arr_loc, 
     &      getGlobalShape(grid, shp, dist_idx))
      
        end subroutine scatterInteger
      
        subroutine 
     &      scatterLogical(grid, arr_glob, arr_loc,  
     &                     shp, dist_idx, local)
          type (dist_grid) :: grid
          logical :: arr_loc(*)
          logical :: arr_glob(*)
          integer, intent(in) :: shp(:)
          integer, intent(in) :: dist_idx
          logical, optional, intent(in) :: local ! unused in serial mode
      
          call localCopyLogical 
     &     (grid, arr_glob, arr_loc, 
     &      getGlobalShape(grid, shp, dist_idx))
      
        end subroutine scatterLogical
      
        subroutine 
     &      localCopyReal8(grid, arr_glob, arr_loc, globalShape)
          type (dist_grid), intent(in) :: grid
          integer, intent(in) :: globalShape(3)
          real(kind=8) :: 
     &      arr_glob(globalShape(1),globalShape(2),globalShape(3))
          real(kind=8) :: 
     &      arr_loc(globalShape(1), 
     &     grid%j_strt_halo:grid%j_stop_halo,globalShape(3))
            
          integer :: j_0, j_1
            
          call get(grid, j_strt=j_0, j_stop=j_1)
          arr_loc(:,j_0:j_1,:) = arr_glob(:,j_0:j_1,:)
          
        end subroutine localCopyReal8
      
        subroutine 
     &      localCopyInteger(grid, arr_glob, arr_loc, globalShape)
          type (dist_grid), intent(in) :: grid
          integer, intent(in) :: globalShape(3)
          integer ::  
     &     arr_glob(globalShape(1),globalShape(2),globalShape(3))
          integer ::  
     &     arr_loc(globalShape(1), 
     &     grid%j_strt_halo:grid%j_stop_halo,globalShape(3))
            
          integer :: j_0, j_1
            
          call get(grid, j_strt=j_0, j_stop=j_1)
          arr_loc(:,j_0:j_1,:) = arr_glob(:,j_0:j_1,:)
          
        end subroutine localCopyInteger
      
        subroutine 
     &      localCopyLogical(grid, arr_glob, arr_loc, globalShape)
          type (dist_grid), intent(in) :: grid
          integer, intent(in) :: globalShape(3)
          logical :: 
     &      arr_glob(globalShape(1),globalShape(2),globalShape(3))
          logical :: 
     &         arr_loc(globalShape(1), 
     &         grid%j_strt_halo:grid%j_stop_halo,globalShape(3))
            
          integer :: j_0, j_1
            
          call get(grid, j_strt=j_0, j_stop=j_1)
          arr_loc(:,j_0:j_1,:) = arr_glob(:,j_0:j_1,:)
          
        end subroutine localCopyLogical
          
        function 
     &    getGlobalShape(grid, localShape, dist_idx) result(globalShape)
          type (dist_grid), intent(in) :: grid
          integer, intent(in) :: localShape(:)
          integer, intent(in) :: dist_idx
      
          integer :: globalShape(3)
      
          globalShape(1) = product(localShape(1:dist_idx-1))
          globalShape(2) = grid%jm_world
          globalShape(3) = product(localShape(dist_idx+1:))
      
        end function getGlobalShape
      
        subroutine pack_1d(grd_dum,ARR,ARR_GLOB)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr(grd_dum%j_strt_halo:)
          real*8, intent(out) :: arr_glob(grd_dum%jm_world)
      
          call gatherReal8 
     &     (grd_dum, arr, arr_glob, shape(arr), 1, all=.false.)
      
        end subroutine pack_1d
      
        subroutine ipack_1d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          integer, intent(in) :: arr(grd_dum%j_strt_halo:)
          integer, intent(out) :: arr_glob(grd_dum%jm_world)
      
          call gatherInteger 
     &     (grd_dum, arr, arr_glob, shape(arr), 1, all=.false.)
      
        end subroutine ipack_1d
      
        subroutine pack_2d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
      
          real*8, intent(in) ::  
     &     arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
          real*8, intent(inout) :: arr_glob(:,:)
      
          call gatherReal8 
     &     (grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)
      
        end subroutine pack_2d
      
        subroutine ipack_2d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          integer, intent(in) :: 
     &      arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
          integer, intent(inout) :: arr_glob(:,:)
      
          call gatherInteger 
     &     (grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)
      
        end subroutine ipack_2d
      
        subroutine lpack_2d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          logical, intent(in) :: 
     &      arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
          logical, intent(inout) :: arr_glob(:,:)
      
          call gatherLogical 
     &     (grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)
      
        end subroutine lpack_2d
      
        subroutine pack_3d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) ::  
     &     arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
          real*8, intent(out) :: arr_glob(:,:,:)
      
          call gatherReal8 
     &     (grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)
      
        end subroutine pack_3d
      
        subroutine ipack_3d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          integer, intent(in) :: 
     &      arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
          integer, intent(out) :: arr_glob(:,:,:)
      
          call gatherInteger 
     &     (grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)
      
        end subroutine ipack_3d
      
        subroutine pack_4d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: 
     &      arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:,:)
          real*8, intent(out) :: arr_glob(:,:,:,:)
      
          call gatherReal8 
     &     (grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)
      
        end subroutine pack_4d
      
        subroutine pack_5d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) ::  
     &     arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:,:,:)
          real*8, intent(out) :: arr_glob(:,:,:,:,:)
      
          call gatherReal8 
     &     (grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)
      
        end subroutine pack_5d
      
        subroutine packj_2d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr(grd_dum%j_strt_halo:,:)
          real*8, intent(inout) :: arr_glob(:,:)
      
          call gatherReal8 
     &     (grd_dum, arr, arr_glob, shape(arr), 1, all=.false.)
      
        end subroutine packj_2d
      
        subroutine ipackj_2d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          integer, intent(in) :: arr(grd_dum%j_strt_halo:,:)
          integer, intent(inout) :: arr_glob(:,:)
      
          call gatherInteger 
     &     (grd_dum, arr, arr_glob, shape(arr), 1, all=.false.)
      
        end subroutine ipackj_2d
      
        subroutine packj_3d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr(grd_dum%j_strt_halo:,:,:)
          real*8, intent(inout) :: arr_glob(:,:,:)
      
          call gatherReal8 
     &     (grd_dum, arr, arr_glob, shape(arr), 1, all=.false.)
      
        end subroutine packj_3d
      
        subroutine packj_4d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr(grd_dum%j_strt_halo:,:,:,:)
          real*8, intent(inout) :: arr_glob(:,:,:,:)
      
          call gatherReal8 
     &     (grd_dum, arr, arr_glob, shape(arr), 1, all=.false.)
      
        end subroutine packj_4d
      
        subroutine unpack_1d(grd_dum,arr_glob,arr)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr_glob(grd_dum%jm_world)
          real*8, intent(out) :: arr(grd_dum%j_strt_halo:)
      
          call scatterReal8(grd_dum, arr_glob, arr, shape(arr), 1)
      
        end subroutine unpack_1d
      
        subroutine unpack_2d(grd_dum,arr_glob,arr,local)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr_glob(:,:)
          real*8, intent(out) :: arr(:,grd_dum%j_strt_halo:)
          logical, optional :: local
      
          call scatterReal8 
     &     (grd_dum, arr_glob, arr, shape(arr), 2, local=local)
      
        end subroutine unpack_2d
      
        subroutine iunpack_2d(grd_dum, arr_glob, arr, local)
          type (dist_grid),  intent(in) :: grd_dum
          integer, intent(in) :: arr_glob(:,:)
          integer, intent(out) :: arr(:,grd_dum%j_strt_halo:)
          logical, optional :: local
      
          call scatterInteger 
     &     (grd_dum, arr_glob, arr, shape(arr), 2, local=local)
      
        end subroutine iunpack_2d
      
        subroutine lunpack_2d(grd_dum, arr_glob, arr, local)
          type (dist_grid),  intent(in) :: grd_dum
          logical, intent(in) :: arr_glob(:,:)
          logical, intent(out) :: arr(:,grd_dum%j_strt_halo:)
          logical, optional :: local
      
          call scatterLogical 
     &     (grd_dum, arr_glob, arr, shape(arr), 2, local=local)
      
        end subroutine lunpack_2d
      
        subroutine unpack_3d(grd_dum, arr_glob, arr, local)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr_glob(:,:,:)
          real*8, intent(out) :: arr(:,grd_dum%j_strt_halo:,:)
          logical, optional :: local
      
          call scatterReal8 
     &     (grd_dum, arr_glob, arr, shape(arr), 2, local=local)
      
        end subroutine unpack_3d
      
        subroutine iunpack_3d(grd_dum,arr_glob,arr,local)
          type (dist_grid),  intent(in) :: grd_dum
          integer, intent(in) :: arr_glob(:,:,:)
          integer, intent(out) :: arr(:,grd_dum%j_strt_halo:,:)
          logical, optional :: local
      
          call scatterInteger 
     &     (grd_dum, arr_glob, arr, shape(arr), 2, local=local)
      
        end subroutine iunpack_3d
      
        subroutine unpack_4d(grd_dum,arr_glob,arr,local)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr_glob(:,:,:,:)
          real*8, intent(out) :: arr(:,grd_dum%j_strt_halo:,:,:)
          logical, optional :: local
          
          call scatterReal8 
     &     (grd_dum, arr_glob, arr, shape(arr), 2, local=local)
      
        end subroutine unpack_4d
        
        subroutine unpack_5d(grd_dum,arr_glob,arr,local)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr_glob(:,:,:,:,:)
          real*8, intent(out) :: arr(:,grd_dum%j_strt_halo:,:,:,:)
          logical, optional :: local
          
          call scatterReal8 
     &     (grd_dum, arr_glob, arr, shape(arr), 2, local=local)
          
        end subroutine unpack_5d
        
        subroutine unpackj_2d(grd_dum,arr_glob,arr,local)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr_glob(:,:)
          real*8, intent(out) :: arr(grd_dum%j_strt_halo:,:)
          logical, optional :: local
          
          call scatterReal8 
     &     (grd_dum, arr_glob, arr, shape(arr), 1, local=local)
          
        end subroutine unpackj_2d
        
        subroutine unpackj_3d(grd_dum,arr_glob,arr,local)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr_glob(:,:,:)
          real*8, intent(out) :: arr(grd_dum%j_strt_halo:,:,:)
          logical, optional :: local
          
          call scatterReal8 
     &     (grd_dum, arr_glob, arr, shape(arr), 1, local=local)
          
        end subroutine unpackj_3d
        
        subroutine unpackj_4d(grd_dum,arr_glob,arr,local)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr_glob(:,:,:,:)
          real*8, intent(out) :: arr(grd_dum%j_strt_halo:,:,:,:)
          logical, optional :: local
      
          call scatterReal8 
     &     (grd_dum, arr_glob, arr, shape(arr), 1, local=local)
          
        end subroutine unpackj_4d
        
        subroutine pack_column_1d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr(:,grd_dum%j_strt_halo:)
          real*8, intent(out) :: arr_glob(:,:)
          
          call gatherReal8 
     &     (grd_dum, arr, arr_glob, shape(arr), 2, all=.false.)
      
        end subroutine pack_column_1d
        
        subroutine pack_column_2d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr(:,:,grd_dum%j_strt_halo:)
          real*8, intent(inout) :: arr_glob(:,:,:)
          
          call gatherReal8 
     &     (grd_dum, arr, arr_glob, shape(arr), 3, all=.false.)
          
        end subroutine pack_column_2d
        
        subroutine pack_column_i2d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          integer, intent(in) :: arr(:,:,grd_dum%j_strt_halo:)
          integer, intent(inout) :: arr_glob(:,:,:)
          integer :: i, k
          
          call gatherInteger 
     &     (grd_dum, arr, arr_glob, shape(arr), 3, all=.false.)
          
        end subroutine pack_column_i2d
        
        subroutine pack_column_3d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr(:,:,grd_dum%j_strt_halo:,:)
          real*8, intent(out) :: arr_glob(:,:,:,:)
          
          call gatherReal8 
     &     (grd_dum, arr, arr_glob, shape(arr), 3, all=.false.)
          
        end subroutine pack_column_3d
        
        subroutine pack_block_3d(grd_dum,arr,arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr(:,:,:,:,grd_dum%j_strt_halo:)
          real*8, intent(out) :: arr_glob(:,:,:,:,:)
          
          call gatherReal8 
     &     (grd_dum, arr, arr_glob, shape(arr), 5, all=.false.)
          
        end subroutine pack_block_3d
        
        subroutine unpack_column_1d(grd_dum,arr_glob,arr,local)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr_glob(:,:)
          real*8, intent(out) :: arr(:,grd_dum%j_strt_halo:)
          logical, optional :: local
      
          call scatterReal8 
     &     (grd_dum, arr_glob, arr, shape(arr), 2, local=local)
      
        end subroutine unpack_column_1d
        
        subroutine unpack_column_2d(grd_dum,arr_glob,arr,local)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr_glob(:,:,:)
          real*8, intent(out) :: arr(:,:,grd_dum%j_strt_halo:)
          logical, optional, intent(in) :: local
          
          call scatterReal8 
     &     (grd_dum , arr_glob, arr, shape(arr), 3, local=local)
      
        end subroutine unpack_column_2d
        
        subroutine iunpack_column_2d(grd_dum,arr_glob,arr,local)
          type (dist_grid),  intent(in) :: grd_dum
          integer, intent(in) :: arr_glob(:,:,:)
          integer, intent(out) :: arr(:,:,grd_dum%j_strt_halo:)
          logical, optional :: local
          
          call scatterInteger 
     &     (grd_dum, arr_glob, arr, shape(arr), 3, local=local)
      
        end subroutine iunpack_column_2d
        
        subroutine unpack_column_3d(grd_dum,arr_glob,arr,local)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr_glob(:,:,:,:)
          real*8, intent(out) :: arr(:,:,grd_dum%j_strt_halo:,:)
          logical, optional, intent(in) :: local
          
          call scatterReal8 
     &     (grd_dum , arr_glob, arr, shape(arr), 3, local=local)
          
        end subroutine unpack_column_3d
      
        subroutine unpack_block_3d(grd_dum, arr_glob, arr, local)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr_glob(:,:,:,:,:)
          real*8, intent(out) :: arr(:,:,:,:,grd_dum%j_strt_halo:)
          logical, optional, intent(in) :: local
          
          call scatterReal8 
     &     (grd_dum , arr_glob, arr, shape(arr), 5, local=local)
          
        end subroutine unpack_block_3d
        
        subroutine ipack_block_2d(grd_dum, arr, arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          integer, intent(in) :: arr(:,:,:,grd_dum%j_strt_halo:)
          integer, intent(inout) :: arr_glob(:,:,:,:)
      
          call gatherInteger(grd_dum, arr, arr_glob, shape(arr), 4)
      
        end subroutine ipack_block_2d
        
        subroutine pack_block_2d(grd_dum, arr, arr_glob)
          type (dist_grid),  intent(in) :: grd_dum
          real*8 , intent(in) :: arr(:,:,:,grd_dum%j_strt_halo:)
          real*8 , intent(inout) :: arr_glob(:,:,:,:)
          
          call gatherReal8 
     &     (grd_dum, arr, arr_glob, shape(arr), 4, all=.false.)
          
        end subroutine pack_block_2d
        
        subroutine iunpack_block_2d(grd_dum,arr_glob,arr,local)
          type (dist_grid),  intent(in) :: grd_dum
          integer, intent(in) :: arr_glob(:,:,:,:)
          integer, intent(out) :: arr(:,:,:,grd_dum%j_strt_halo:)
          logical, optional :: local
          
          call scatterInteger 
     &     (grd_dum, arr_glob, arr, shape(arr), 4, local=local)
      
        end subroutine iunpack_block_2d
        
        subroutine unpack_block_2d(grd_dum,arr_glob,arr,local)
          type (dist_grid),  intent(in) :: grd_dum
          real*8, intent(in) :: arr_glob(:,:,:,:)
          real*8, intent(out) :: arr(:,:,:,grd_dum%j_strt_halo:)
          logical, optional, intent(in) :: local
      
          call scatterReal8 
     &     (grd_dum , arr_glob, arr, shape(arr), 4, local=local)
      
        end subroutine unpack_block_2d
        
      end module GatherScatter_mod

!==============================================================================
! LTM: ESMF_Interface/GlobalSum_mod.F90
!==============================================================================

      module GlobalSum_mod
!@sum Collection of routines to support MPI reduce operations for various
!@+   data types and array sizes
!@auth Tom Clune GSFC/SSSO/610.3
      !use MpiSupport_mod, only: ROOT_PROCESS
      use dist_grid_mod
      use GatherScatter_mod, only: gatherReal8
      implicit none
      private
      
!@var GLOBALSUM output a bit-reproducible global-hemisphere-zonal sum for an array
      public :: globalsum
      
      interface globalsum
      module procedure globalsum_int_reduce
      module procedure globalsum_j
      module procedure globalsum_ij
      module procedure globalsum_ijk
      module procedure globalsum_ijk_ik
      module procedure globalsum_other_ijk
      module procedure globalsum_other_ijk_ik ! not used?
      module procedure globalsum_jk
      module procedure globalsum_xxxj_xxx
      module procedure globalsum_xxxij_xxx
      end interface
      
!@var checksum output a bit-reproducible checksum for an array
      public :: checksum        ! Communicate overlapping portions of subdomains
      public :: checksumJ       ! Communicate overlapping portions of subdomains
      public :: checksum_column ! (k, i, j)
      
      interface checksum
      module procedure checksum_1d
      module procedure checksum_2d
      module procedure checksum_3d
      end interface
      
      interface checksumj
      module procedure checksumj_2d
      module procedure checksumj_3d
      end interface
      
      interface checksum_column
      module procedure checksum_column_2d
      module procedure checksum_column_3d
      module procedure int_checksum_column_3d
      module procedure checksum_column_4d
      end interface
            
      contains
      
!---------------------------
      subroutine globalsum_int_reduce(grd_dum, ivar, isum, all)
!@sum Output a bit-reproducible global-hemisphere-zonal sum for an integer array
      type(dist_grid), intent(in) :: grd_dum
      integer, intent(in)  :: ivar
      integer, intent(out) :: isum
      logical, optional, intent(in) :: all
      
      logical :: all_
      integer :: ierr
      
      all_ = .false.
      if (present(all)) all_ = all
      
      isum = ivar
      
      end subroutine globalsum_int_reduce
      
      subroutine globalSum_J
     &     (grd_dum, arr, gsum, hsum, istag, iskip, all, jband)
!@sum Output a bit-reproducible global-hemisphere-zonal sum for a real 1-D array
      type (dist_grid),  intent(IN) :: grd_dum
      real*8,            intent(IN) :: arr(grd_dum%j_strt_halo:)
      real*8,            intent(OUT):: gsum
      real*8, optional,  intent(OUT):: hsum(2)
      integer,optional,  intent(IN) :: istag
      integer,optional,  intent(IN) :: iskip
      logical,optional,  intent(IN) :: all
      integer, optional, intent(IN) :: jband(2)
      
      integer :: i_0, i_1, j_0, j_1, IM, JM, J, ierr
      real*8  :: garr(grd_dum%jm_world)
      logical :: istag_, iskip_
      
      integer :: JSTRT, JSTOP
! now local
      
      call get(grd_dum, J_STRT=j_0, J_STOP=j_1)
      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      
      istag_ = .false.
      if (present(istag)) then
         if (istag == 1) istag_ = .true.
      end if
      
      iskip_ = .false.
      if (present(iskip)) then
         if (iskip == 1) iskip_ = .true.
      end if
      
      JSTRT=1
      JSTOP=JM
      
      if (present(jband)) then
         JSTRT=jband(1)
         JSTOP=jband(2)
      end if
      
      call gatherReal8(grd_dum, arr, garr, shape(arr), 1, .false.)
      
      if (AM_I_ROOT()) then
         if (istag_) then
            gsum = sum(garr(2:JM),1)
         elseif (iskip_) then
            gsum = sum(garr(2:JM-1),1)
         else
            gsum = sum(garr(JSTRT:JSTOP),1)
         endif
         if (present(hsum)) then
            if (istag_) then
               hsum(1)   = sum( garr(2     :JM/2),1   )
               hsum(2)   = sum( garr(2+JM/2:JM  ),1   )
               hsum(1)   = hsum(1) + 0.5*garr(1+JM/2)
               hsum(2)   = hsum(2) + 0.5*garr(1+JM/2)
            else
               hsum(1)   = Sum( garr(1     :JM/2),1   )
               hsum(2)   = sum( garr(1+JM/2:JM  ),1   )
            endif
         endif
      endif
      
      end subroutine globalSum_J
      
      subroutine globalSum_IJ
     &     (grd_dum, arr, gsum, hsum, zsum, istag, all, iskip)
!@sum Output a bit-reproducible global-hemisphere-zonal sum for a real 2-D array
      type (dist_grid),  intent(IN) :: grd_dum
      real*8,            intent(IN) :: arr(grd_dum%i_strt_halo:, 
     &                                     grd_dum%j_strt_halo:)
      real*8,            intent(OUT):: gsum
      real*8, OPTIONAL,  INTENT(OUT):: hsum(2)
      real*8, OPTIONAL,  INTENT(OUT):: zsum(grd_dum%j_strt: 
     &                                      grd_dum%j_stop)
      integer,OPTIONAL,  INTENT(IN) :: istag
      integer,OPTIONAL,  INTENT(IN) :: iskip
      
      integer :: i_0, i_1, j_0, j_1, IM, JM, J,J_0STG, ierr, J_0S,J_1S
      real*8, allocatable  :: zon(:)
      real*8  :: garr(grd_dum%jm_world)
      logical :: istag_,iskip_
      logical,OPTIONAL,  INTENT(IN) :: all
      
! now local
      
      call get(grd_dum, j_strt=j_0, j_stop=j_1)
      allocate(zon(j_0:j_1))
      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0STG  = grd_dum%j_strt_stgr
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      j_0S = grd_dum%j_strt_skp
      j_1S  = grd_dum%j_stop_skp
      
      istag_ = .false.
      if (present(istag)) Then
         if (istag == 1) istag_ = .true.
      end if
      
      iskip_ = .false.
      if (present(iskip)) Then
         if (iskip == 1) iskip_ = .true.
      end if
      
      if (istag_) Then
         zon(J_0STG:J_1)  = sum(arr(:,j_0STG:j_1),1)
      Else if (iskip_) Then
         zon(J_0S:J_1S)  = sum(arr(:,j_0S:j_1S),1)
      Else
         zon  = sum(arr(:,j_0:j_1),1)
      end if
      
      Call gatherReal8(grd_dum, zon, garr, shape(zon), 1, .false.)
      
      if (AM_I_ROOT()) then
         if (istag_) then
            gsum = sum(garr(2:JM),1)
         Else if(iskip_) then
            gsum = sum(garr(2:JM-1),1)
         Else
            gsum = sum(garr(1:JM),1)
         endif
         if (present(hsum)) then
            if (istag_) then
               hsum(1)   = Sum( garr(2     :JM/2),1   )
               hsum(2)   = Sum( garr(2+JM/2:JM  ),1   )
               hsum(1)   = hsum(1) + 0.5*garr(1+JM/2)
               hsum(2)   = hsum(2) + 0.5*garr(1+JM/2)
            Else
               hsum(1)   = Sum( garr(1     :JM/2),1   )
               hsum(2)   = Sum( garr(1+JM/2:JM  ),1   )
            endif
         endif
      endif
      if (present(zsum)) zsum = zon
      
      end subroutine globalSum_IJ
      
      subroutine globalSum_IJK
     &     (grd_dum, arr, gsum, hsum, zsum, istag, iskip)
!@sum Output a bit-reproducible global-hemisphere-zonal sum for a real 3-D array
      type (dist_grid),   INTENT(IN) :: grd_dum
      real*8,             INTENT(IN) :: arr(grd_dum%i_strt_halo:, 
     &                                      grd_dum%j_strt_halo:,:)
      real*8,             INTENT(OUT):: gsum(size(arr,3))
      real*8, OPTIONAL,   INTENT(OUT):: hsum(2,size(arr,3))
      real*8, OPTIONAL,   INTENT(OUT):: zsum(grd_dum%j_strt: 
     &                                       grd_dum%j_stop, 
     &                                       size(arr,3))
      integer,OPTIONAL,   INTENT(IN) :: istag,iskip
      
      integer :: k
      integer :: i_0, i_1, j_0, j_1, j_0S, j_1S, IM, JM, j_0STG
      real*8  :: zon(grd_dum%j_strt:grd_dum%j_stop,size(arr,3))
      real*8  :: garr(grd_dum%jm_world,size(arr,3))
      logical :: istag_,iskip_
! now local
      
      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_0S = grd_dum%j_strt_skp
      j_0STG  = grd_dum%j_strt_stgr
      j_1  = grd_dum%j_stop
      j_1S = grd_dum%j_stop_skp
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      
      istag_ = .false.
      if (present(istag)) Then
         if (istag == 1) istag_ = .true.
      end if
      
      iskip_ = .false.
      if (present(iskip)) Then
         if (iskip == 1) iskip_ = .true.
      end if
      
      if (istag_) Then
         zon(j_0STG:j_1,:)  = sum(arr(:,j_0STG:j_1,:),1)
      Else if (iskip_) Then
         zon(j_0S:j_1S,:)  = sum(arr(:,j_0S:j_1S,:),1)
      Else
         zon  = sum(arr(:,j_0:j_1,:),1)
      end if
      
      Call gatherReal8(grd_dum, zon, garr, shape(zon), 1, .false.)
      
      if (AM_I_ROOT()) then
         if (istag_) then
            gsum = sum(garr(2:JM,:),1)
         Else if(iskip_) then
            gsum = sum(garr(2:JM-1,:),1)
         Else
            gsum = sum(garr(1:JM,:),1)
         endif
         if (present(hsum)) then
            if (istag_) then
               hsum(1,:)   = Sum( garr(2     :JM/2,:),1   )
               hsum(2,:)   = Sum( garr(2+JM/2:JM  ,:),1   )
               hsum(1,:)   = hsum(1,:) + 0.5*garr(1+JM/2,:)
               hsum(2,:)   = hsum(2,:) + 0.5*garr(1+JM/2,:)
            Else
               hsum(1,:)   = Sum( garr(1     :JM/2,:),1   )
               hsum(2,:)   = Sum( garr(1+JM/2:JM  ,:),1   )
            endif
         endif
      endif
      if (present(zsum)) zsum = zon
      
      end subroutine globalSum_IJK
      
      subroutine globalSum_OTHER_IJK(grd_dum, arr, gsum, jband, all)
!@sum Output a bit-reproducible global-hemisphere-zonal sum for a real 3-D array
      type (dist_grid),   intent(IN) :: grd_dum
      real*8,             intent(IN) :: arr(:,grd_dum%j_strt_halo:)
      real*8,             intent(OUT):: gsum(size(arr,1))
      integer,           intent(IN) :: jband(2)
      logical,optional,   intent(IN) :: all
      
      integer :: k
      integer :: i_0, i_1, j_0, j_1, IM, JM, jb1, jb2
      logical :: all_
      integer :: ierr
! now local
      
      all_ = .false.
      if (present(all)) all_ = all
      
      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM = SIZE(arr,1)
      JM   = grd_dum%JM_WORLD
      
      jb1 = jband(1)
      jb2 = jband(2)
      
      gsum = sum(arr(:, max(jb1,j_0):min(jb2,j_1) ),2)

      end subroutine globalSum_OTHER_IJK
      
      subroutine globalSum_OTHER_IJK_IK(grd_dum, arr, gsum, jband, all)
      type (dist_grid),   intent(IN) :: grd_dum
      real*8,             intent(IN) :: arr(:,grd_dum%j_strt_halo:,:)
      real*8,             intent(OUT):: gsum(size(arr,1), size(arr,3))
      integer,            intent(IN) :: jband(2)
      logical,optional,   intent(IN) :: all
      
      integer :: k
      integer :: i_0, i_1, j_0, j_1, IM, JM, jb1, jb2
      logical :: all_
      integer :: ierr
! now local
      
      all_ = .false.
      if (present(all)) all_ = all
      
      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM = SIZE(arr,1)
      JM   = grd_dum%JM_WORLD
      
      jb1 = jband(1)
      jb2 = jband(2)
      
      gsum = sum(arr(:, max(jb1,j_0):min(jb2,j_1) ,:),2)

      end subroutine globalSum_OTHER_IJK_IK
      
      subroutine globalSum_IJK_IK(grd_dum, arr, gsum, all)
!@sum Output a bit-reproducible global-hemisphere-zonal sum for a real array
      type (dist_grid),   intent(IN) :: grd_dum
      real*8,             intent(IN) :: arr(:,grd_dum%j_strt_halo:,:)
      real*8,             intent(OUT):: gsum(size(arr,1), size(arr,3))
      logical,optional,   intent(IN) :: all
      
      integer :: i_0, i_1, j_0, j_1, IM, JM, LM
      logical :: all_
      integer :: ierr

! now local
      
      all_ = .false.
      if (present(all)) all_ = all
      
      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM = SIZE(arr,1)
      JM   = grd_dum%JM_WORLD
      LM =  size(arr,3)
      
      gsum = sum(arr(:,j_0:j_1,:),2)

      end subroutine globalSum_IJK_IK
      
      subroutine globalSum_JK(grd_dum, arr, gsum, hsum, istag, all)
!@sum Output a bit-reproducible global-hemisphere-zonal sum for a real array
      type (dist_grid),  intent(IN) :: grd_dum
      real*8,            intent(IN) :: arr(grd_dum%j_strt_halo:,:)
      real*8,            intent(OUT):: gsum(size(arr,2))
      real*8, optional,  intent(OUT):: hsum(2,size(arr,2))
      integer,optional,  intent(IN) :: istag
      logical,optional,  intent(IN) :: all
      
      integer :: k
      integer :: ierr
      integer :: i_0, i_1, j_0, j_1, IM, JM
      real*8  :: garr(grd_dum%jm_world,size(arr,2))
      logical :: istag_
      
! now local
      
      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      istag_ = .false.
      if (present(istag)) then
         if (istag == 1) istag_ = .true.
      end if
      
      call gatherReal8(grd_dum, arr, garr, shape(arr), 1, .false.)
      
      if (AM_I_ROOT()) then
         if (istag_) then
            gsum = sum(garr(2:JM,:),1)
         else
            gsum = sum(garr(1:JM,:),1)
         endif
         if (present(hsum)) then
            if (istag_) then
               hsum(1,:)   = sum( garr(2     :JM/2,:),1   )
               hsum(2,:)   = sum( garr(2+JM/2:JM  ,:),1   )
               hsum(1,:)   = hsum(1,:) + 0.5*garr(1+JM/2,:)
               hsum(2,:)   = hsum(2,:) + 0.5*garr(1+JM/2,:)
            else
               hsum(1,:)   = sum( garr(1     :JM/2,:),1   )
               hsum(2,:)   = sum( garr(1+JM/2:JM  ,:),1   )
            endif
         endif
      endif
      
      end subroutine globalSum_JK
      
      subroutine globalSum_XXXJ_XXX(grd_dum, arr, gsum, all)
!@sum Output a bit-reproducible global-hemisphere-zonal sum for a real array
      type (dist_grid), intent(IN) :: grd_dum
      real*8, intent(IN) :: arr(:,:,:,grd_dum%j_strt_halo:)
      real*8, intent(Out) :: gsum(:,:,:)
      logical, optional, intent(IN) :: all
      
      integer :: k
      integer :: i_0, i_1, j_0, j_1, IM, JM
      real*8  :: garr(size(arr,1),size(arr,2),size(arr,3), 
     &     grd_dum%jm_world)
      logical :: istag_,all_
      
      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      
      gsum = SUM(arr(:,:,:,J_0:J_1),4)
      
      end subroutine globalSum_XXXJ_XXX
      
      subroutine globalSum_XXXIJ_XXX(grd_dum, arr, gsum, all)
!@sum Output a bit-reproducible global-hemisphere-zonal sum for a real array
      type (dist_grid), INTENT(IN) :: grd_dum
      Real*8, INTENT(IN) :: 
     &     arr(:,:,:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      Real*8, INTENT(Out) :: gsum(:,:,:)
      Logical, Optional, INTENT(IN) :: all
      
      integer :: k
      integer :: i_0, i_1, j_0, j_1, IM, JM
      real*8  :: garr(size(arr,1),size(arr,2),size(arr,3), 
     &     grd_dum%jm_world)
      real*8  :: larr(size(arr,1),size(arr,2),size(arr,3), 
     &     grd_dum%j_strt_halo:grd_dum%j_stop_halo)
      logical :: istag_,all_
      
      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      
      gsum = SUM(SUM(arr(:,:,:,I_0:I_1,J_0:J_1),4),4)
      
      end subroutine globalSum_XXXIJ_XXX
      
      subroutine checksum_1D(grd_dum, arr, line, file, unit, STGR, SKIP)
!@sum Output global-hemisphere-zonal checksum for a real 1-D array
      TYPE (dist_grid),   INTENT(IN) :: grd_dum
      real*8,            INTENT(IN) :: 
     &                arr(grd_dum%j_strt_halo:)
      integer,           INTENT(IN) :: line
      character(LEN=*),  INTENT(IN) :: file
      integer, OPTIONAL, INTENT(IN) :: unit
      logical, OPTIONAL, INTENT(IN) :: stgr
      logical, OPTIONAL, INTENT(IN) :: skip
      
      
      integer :: unit_
      real*8  :: asum, L1norm
      real*8  :: t_arr(grd_dum%j_strt_halo:grd_dum%j_stop_halo)
      integer :: J_0, J_1
      integer :: stgr_, skip_
      
      end subroutine checksum_1D
      
      subroutine checksum_2D(grd_dum, arr, line, file, unit, stgr, 
     &     skip)
!@sum Output global-hemisphere-zonal checksum for a real 2-D array
      TYPE (dist_grid),   INTENT(IN) :: grd_dum
      real*8,            INTENT(IN) :: 
     &                arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      integer,           INTENT(IN) :: line
      character(LEN=*),  INTENT(IN) :: file
      integer, OPTIONAL, INTENT(IN) :: unit
      logical, OPTIONAL, INTENT(IN) :: stgr
      logical, OPTIONAL, INTENT(IN) :: skip
      
      integer :: unit_
      real*8  :: asum, L1norm
      real*8  :: asum_glob, L1norm_glob
      real*8 :: 
     &  t_arr(size(arr,1),grd_dum%j_strt_halo:grd_dum%j_stop_halo)
      integer :: J_0, J_1, I,J
      integer :: stgr_,skip_
      
      J_0 = grd_dum%J_STRT
      J_1 = grd_dum%J_STOP
      
      unit_ = checksum_UNIT     ! default
      if (present(unit)) unit_ = unit
      
      stgr_ = 0
      if (present(stgr)) THEN
         if (stgr) stgr_=1
      end if
      
      skip_ = 0
      if (present(skip)) THEN
         if (skip) skip_=1
      end if
      
      t_arr(:,J_0:J_1) = arr(:,J_0:J_1)
      Call globalSum(grd_dum,      t_arr,asum, istag=stgr_,iskip=skip_)
      t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
      Call globalSum(grd_dum,      t_arr,L1norm,istag=stgr_,iskip=skip_)
      
      if (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') 
     &     file,line, asum, L1norm
      
      
      end subroutine checksum_2D
      
      subroutine checksum_3D(grd_dum, arr, line, file, unit, stgr, skip)
!@sum Output global-hemisphere-zonal checksum for a real 3-D array
      TYPE (dist_grid),   INTENT(IN) :: grd_dum
      real*8,            INTENT(IN) :: 
     &                arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      integer,           INTENT(IN) :: line
      character(LEN=*),  INTENT(IN) :: file
      integer, OPTIONAL, INTENT(IN) :: unit
      logical, OPTIONAL, INTENT(IN) :: stgr
      logical, OPTIONAL, INTENT(IN) :: skip
      
      
      integer :: unit_
      integer :: k
      real*8, DIMENSION(Size(arr,3))  :: asum, L1norm
      
      real*8 :: 
     &  t_arr(size(arr,1),grd_dum%j_strt_halo:grd_dum%j_stop_halo)
      integer :: J_0, J_1
      Integer :: stgr_,skip_
      
      
      J_0 = grd_dum%J_STRT
      J_1 = grd_dum%J_STOP
      
      unit_ = checksum_UNIT     ! default
      if (present(unit)) unit_ = unit
      
      stgr_ = 0
      if (present(stgr)) THEN
         if (stgr) stgr_=1
      end if
      
      skip_ = 0
      if (present(skip)) THEN
         if (skip) skip_=1
      end if
      
      Do k = 1, Size(arr, 3)
         t_arr(:,J_0:J_1) = arr(:,J_0:J_1,k)
         Call globalSum
     &        (grd_dum, t_arr, asum(k), istag=stgr_,iskip=skip_)
         t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
         Call globalSum
     &        (grd_dum, t_arr,L1norm(k),istag=stgr_,iskip=skip_)
      end Do
      if (AM_I_ROOT()) Then
         Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') 
     &        file,line, Sum(asum), Sum(L1norm)
      end if
      
      end subroutine checksum_3D
      
      subroutine checksumj_2D(grd_dum, arr, line, file, unit, stgr, 
     &     skip)
!@sum Output global-hemisphere-zonal checksum for a real 2-D array
      TYPE (dist_grid),   INTENT(IN) :: grd_dum
      real*8,            INTENT(IN) :: 
     &                arr(grd_dum%j_strt_halo:,:)
      integer,           INTENT(IN) :: line
      character(LEN=*),  INTENT(IN) :: file
      integer, OPTIONAL, INTENT(IN) :: unit
      logical, OPTIONAL, INTENT(IN) :: stgr
      logical, OPTIONAL, INTENT(IN) :: skip
      
      integer :: unit_
      real*8  :: asum, L1norm
      real*8  :: asum_glob, L1norm_glob
      real*8 :: 
     &  t_arr(size(arr,2),grd_dum%j_strt_halo:grd_dum%j_stop_halo)
      integer :: J_0, J_1
      integer :: stgr_,skip_
      
      J_0 = grd_dum%J_STRT
      J_1 = grd_dum%J_STOP
      
      unit_ = checksum_UNIT     ! default
      if (present(unit)) unit_ = unit
      
      stgr_ = 0
      if (present(stgr)) THEN
         if (stgr) stgr_=1
      end if
      
      skip_ = 0
      if (present(skip)) THEN
         if (skip) skip_=1
      end if
      
      t_arr(:,J_0:J_1) = Transpose(arr(J_0:J_1,:))
      Call globalSum(grd_dum,      t_arr,asum, istag=stgr_,iskip=skip_)
      t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
      Call globalSum(grd_dum,      t_arr,L1norm,istag=stgr_,iskip=skip_)
      
      if (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') 
     &     file,line, asum, L1norm
      
      
      end subroutine checksumj_2D
      
      subroutine checksumj_3D(grd_dum, arr, line, file, unit, stgr)
!@sum Output global-hemisphere-zonal checksum for a real 3-D array
      TYPE (dist_grid),   INTENT(IN) :: grd_dum
      real*8,            INTENT(IN) :: 
     &                arr(grd_dum%j_strt_halo:,:,:)
      integer,           INTENT(IN) :: line
      character(LEN=*),  INTENT(IN) :: file
      integer, OPTIONAL, INTENT(IN) :: unit
      logical, OPTIONAL, INTENT(IN) :: stgr
      
      
      integer :: unit_
      integer :: k
      real*8, DIMENSION(Size(arr,3))  :: asum, L1norm
      
      real*8 :: 
     &  t_arr(size(arr,2),grd_dum%j_strt_halo:grd_dum%j_stop_halo)
      integer :: J_0, J_1
      Integer :: stgr_
      
      
      J_0 = grd_dum%J_STRT
      J_1 = grd_dum%J_STOP
      
      unit_ = checksum_UNIT     ! default
      if (present(unit)) unit_ = unit
      
      stgr_ = 0
      if (present(stgr)) THEN
         if (stgr) stgr_=1
      end if
      
      Do k = 1, Size(arr, 3)
         t_arr(:,J_0:J_1) = Transpose(arr(J_0:J_1,:,k))
         Call globalSum(grd_dum, t_arr, asum(k), istag=stgr_)
         t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
         Call globalSum(grd_dum, t_arr, L1norm(k), istag=stgr_)
      end Do
      if (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') 
     &     file,line, Sum(asum), Sum(L1norm)
      
      
      end subroutine checksumj_3D
      
      subroutine checksum_COLUMN_2D(grd_dum,arr,line,file,unit,stgr)
!@sum Output global-hemisphere-zonal checksum for a real array
      TYPE (dist_grid),   INTENT(IN) :: grd_dum
      real*8,            INTENT(IN) :: 
     &                arr(:,grd_dum%j_strt_halo:)
      integer,           INTENT(IN) :: line
      character(LEN=*),  INTENT(IN) :: file
      integer, OPTIONAL, INTENT(IN) :: unit
      logical, OPTIONAL, INTENT(IN) :: stgr
      
      
      integer :: unit_
      integer :: k
      real*8, DIMENSION(SIZE(arr,1))  :: asum, L1norm
      integer :: stgr_
      
      unit_ = checksum_UNIT     ! default
      if (present(unit)) unit_ = unit
      stgr_ = 0
      if (present(stgr)) THEN
         if (stgr) stgr_=1
      end if
      
      Do k = 1, Size(arr, 1)
         Call globalSum(grd_dum,      arr(k,:), asum(k), istag=stgr_)
         Call globalSum(grd_dum, Abs(arr(k,:)), L1norm(k), istag=stgr_)
      end DO
      
      if (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') 
     &     file,line, Sum(asum), Sum(L1norm)
      
      end subroutine checksum_COLUMN_2D
      
      subroutine checksum_COLUMN_3D(grd_dum,arr,line,file,unit,stgr)
!@sum Output global-hemisphere-zonal checksum for a real array
      TYPE (dist_grid),   INTENT(IN) :: grd_dum
      real*8,            INTENT(IN) :: 
     &                arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      integer,           INTENT(IN) :: line
      character(LEN=*),  INTENT(IN) :: file
      integer, OPTIONAL, INTENT(IN) :: unit
      logical, OPTIONAL, INTENT(IN) :: stgr
      
      
      integer :: unit_
      integer :: k
      real*8, DIMENSION(Size(arr,1))  :: asum, L1norm
      integer :: stgr_
      
      unit_ = checksum_UNIT     ! default
      if (present(unit)) unit_ = unit
      stgr_ = 0
      if (present(stgr)) THEN
         if (stgr) stgr_=1
      end if
      
      Do k = 1, Size(arr, 1)
         Call globalSum(grd_dum,      arr(k,:,:), asum(k), istag=stgr_)
         Call globalSum
     &        (grd_dum, Abs(arr(k,:,:)), L1norm(k), istag=stgr_)
      end DO
      
      if (AM_I_ROOT()) Then
         Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') 
     &        file,line, Sum(asum), Sum(L1norm)
         CALL FLUSH(unit_)
      end if
      
      end subroutine checksum_COLUMN_3D
      
      subroutine INT_checksum_COLUMN_3D(grd_dum, arr, line, file, unit)
!@sum Output global-hemisphere-zonal checksum for an integer array
      TYPE (dist_grid),   INTENT(IN) :: grd_dum
      integer,           INTENT(IN) :: 
     &                arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      integer,           INTENT(IN) :: line
      character(LEN=*),  INTENT(IN) :: file
      integer, OPTIONAL, INTENT(IN) :: unit
      
      
      integer :: unit_
      integer :: k
      real*8, DIMENSION(Size(arr,1))  :: asum, L1norm
      
      Call checksum_COLUMN(grd_dum, Real(arr,KIND=KIND(1.0D+0)), line, 
     &     file, unit)
      
      end subroutine INT_checksum_COLUMN_3D
      
      subroutine checksum_COLUMN_4D(grd_dum,arr,line,file,unit,stgr, 
     &                              skip)
!@sum Output global-hemisphere-zonal checksum for an integer array
      TYPE (dist_grid),   INTENT(IN) :: grd_dum
      real*8,            INTENT(IN) :: 
     &                arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      integer,           INTENT(IN) :: line
      character(LEN=*),  INTENT(IN) :: file
      integer, OPTIONAL, INTENT(IN) :: unit
      logical, OPTIONAL, INTENT(IN) :: stgr
      logical, OPTIONAL, INTENT(IN) :: skip
      
      
      integer :: unit_
      integer :: i_0, i_1, j_0, j_1, k
      real*8, DIMENSION(Size(arr,1),Size(arr,4)) :: asum, L1norm
      Integer :: stgr_,skip_
      
      unit_ = checksum_UNIT     ! default
      if (present(unit)) unit_ = unit
      
      i_0 = grd_dum%i_strt
      i_1 = grd_dum%i_stop
      j_0 = grd_dum%j_strt
      j_1 = grd_dum%j_stop
      
      stgr_ = 0
      if (present(stgr)) THEN
         if (stgr) stgr_=1
      end if
      
      skip_ = 0
      if (present(skip)) THEN
         if (skip) skip_=1
      end if
      
      Do k = 1, Size(arr,1)
         CALL globalSum(grd_dum,    arr(k,:,:,:),   asum(k,:), 
     &        istag=stgr_, iskip=skip_)
         CALL globalSum(grd_dum,Abs(arr(k,:,:,:)),L1norm(k,:), 
     &       istag=stgr_, iskip=skip_)
      end Do
      
      if (AM_I_ROOT()) Then
         Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') 
     &       file,line, SUM(asum), SUM(L1norm)
      end if
      
      end subroutine checksum_COLUMN_4D
      
      end module GlobalSum_mod

!==============================================================================
! LTM: ESMF_Interface/DomainDecomposition_mod.F90
!==============================================================================

      module DOMAIN_DECOMP_1D
!     @sum  This module conveniently wraps (lumps together) various domain 
!     @+ decomposition modules.
!     @auth Tom Clune GSFC/SSSO/610.3
      use dist_grid_mod
      use Halo_mod
!      use SpecialIO_mod
      use GatherScatter_mod
      use GlobalSum_mod
      implicit none
      public

      contains
      
      end module DOMAIN_DECOMP_1D

      MODULE PBLCOM
!@sum  PBLCOM contains the arrays used by the Boundary Layer code
!@auth Greg Hartke/Ye Cheng
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
      IMPLICIT NONE
      SAVE
!@var TSAVG     COMPOSITE SURFACE AIR TEMPERATURE (K)
      REAL*8, allocatable, dimension(:,:) :: tsavg

      END MODULE PBLCOM

      SUBROUTINE ALLOC_PBL_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE PBLCOM
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID, GET
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_1H, I_0H, J_1H, J_0H
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      ALLOCATE(    
     *              tsavg(I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      END SUBROUTINE ALLOC_PBL_COM

!%! !==============================================================================
!%! ! LTM: ESMF_Interface/ESMF_MOD.F90
!%! !==============================================================================
!%! 
!%!       Module ESMF_MOD_private
!%!       Implicit None
!%!       Private
!%!       
!%! !     Derived Types
!%!       Public :: Grid
!%!       Public :: Array
!%!       Public :: Field
!%!       Public :: AxisIndex
!%!       Public :: RouteHandle
!%! 
!%! !     Methods
!%!       Public :: Field_SetDataPointer
!%!       Public :: Field_GetDataPointer
!%!       Public :: Field_GetArray
!%!       Public :: Field_SetDistRank, Field_GetDistRank
!%!       Public :: Field_GetRank, Field_GetShape
!%!       Public :: Field_Clone
!%!       Public :: Field_Destroy
!%!       Public :: Grid_Create
!%!       Public :: Grid_Destroy
!%!       Public :: Grid_Get
!%!       
!%! !     Constants
!%!       Public :: ESMF_KIND_R4
!%!       Public :: ESMF_KIND_R8
!%!       Integer, Parameter :: ESMF_KIND_R4 = Selected_Real_Kind(P=6)
!%!       Integer, Parameter :: ESMF_KIND_R8 = Selected_Real_Kind(P=14)
!%!       
!%! !     Constants
!%!       Public :: ESMF_FAILURE
!%!       Public :: ESMF_SUCCESS
!%!       Integer, Parameter :: ESMF_FAILURE = -1
!%!       Integer, Parameter :: ESMF_SUCCESS = 0
!%!       
!%!       Public :: TYPE_REAL
!%!       Public :: TYPE_INTEGER
!%!       Integer, Parameter :: TYPE_REAL = 1
!%!       Integer, Parameter :: TYPE_INTEGER = 2
!%!       
!%!       Public :: ESMF_CELL_NFACE
!%!       Public :: ESMF_CELL_SFACE
!%!       Public :: ESMF_CELL_CENTER
!%!       Public :: ESMF_CELL_SWCORNER
!%! 
!%!       Integer, Parameter :: ESMF_CELL_NFACE=1
!%!       Integer, Parameter :: ESMF_CELL_SFACE=2
!%!       Integer, Parameter :: ESMF_CELL_CENTER=3
!%!       Integer, Parameter :: ESMF_CELL_SWCORNER=4
!%! 
!%!       Public :: MAXSTR
!%!       Integer, Parameter :: MAXSTR = 100
!%!       Integer, Parameter :: N_DIMENSIONS = 3 ! dimensions of grid
!%! 
!%!       Public :: ESMF_DATATYPE, ESMF_DATAKIND
!%!       Public :: ESMF_DATA_REAL, ESMF_DATA_INTEGER
!%!       Public :: ESMF_R8, ESMF_R4
!%!       Public :: ESMF_I8, ESMF_I4
!%! 
!%!       Public :: Operator(==)
!%!       interface operator (.eq.)
!%!       module procedure ESMF_type_eq
!%!       module procedure ESMF_kind_eq
!%!       End interface
!%! 
!%!       Type ESMF_DATATYPE
!%!       Integer :: type
!%!       End Type ESMF_DATATYPE
!%! 
!%!       Type ESMF_DATAKIND
!%!       Integer :: kind
!%!       End Type ESMF_DATAKIND
!%! 
!%!       Type (ESMF_DATATYPE), Parameter :: ESMF_DATA_REAL=ESMF_DATATYPE(1)
!%!       Type (ESMF_DATATYPE),Parameter::ESMF_DATA_INTEGER=ESMF_DATATYPE(2)
!%!       Type (ESMF_DATAKIND), Parameter :: ESMF_R8=ESMF_DATAKIND(1)
!%!       Type (ESMF_DATAKIND), Parameter :: ESMF_R4=ESMF_DATAKIND(2)
!%!       Type (ESMF_DATAKIND), Parameter :: ESMF_I8=ESMF_DATAKIND(3)
!%!       Type (ESMF_DATAKIND), Parameter :: ESMF_I4=ESMF_DATAKIND(4)
!%!       
!%!       type AxisIndex
!%!       sequence
!%!       integer :: min
!%!       integer :: max
!%!       integer :: stride
!%!       end type AxisIndex
!%! 
!%!       Type Grid
!%!       Private
!%!       Integer :: IM
!%!       Integer :: JM
!%!       Integer :: LM
!%!       Type (AxisIndex) :: global(N_DIMENSIONS)
!%!       End Type Grid
!%! 
!%! !!$  Type HaloDirection
!%! !!$     Private
!%! !!$     Integer :: direction ! NORTH, SOUTH, EAST, WEST
!%! !!$  End Type HaloDirection
!%! !!$  Integer, Public, Parameter :: NORTH = 2**0
!%! !!$  Integer, Public, Parameter :: SOUTH = 2**1
!%!       Integer, Parameter :: MAX_RANK=5
!%!       Type Field
!%!       Private
!%! 
!%!       Logical :: active = .false.
!%!       Integer :: esmf_kind = -1
!%!       Integer :: esmf_type = -2
!%!       Integer :: rank, dist_rank
!%!       Integer :: data_shape(MAX_RANK)
!%!       Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_1D (:) => Null()
!%!       Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_2D (:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R4), Pointer :: 
!%!      &     ptr_Real4_3D (:,:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R4), Pointer :: 
!%!      &     ptr_Real4_4D (:,:,:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R4), Pointer :: 
!%!      &     ptr_Real4_5D (:,:,:,:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_1D (:) => Null()
!%!       Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_2D (:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R8), Pointer :: 
!%!      &     ptr_Real8_3D (:,:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R8), Pointer :: 
!%!      &     ptr_Real8_4D (:,:,:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R8), Pointer :: 
!%!      &     ptr_Real8_5D (:,:,:,:,:) => Null()
!%!       Integer , Pointer :: ptr_Integer_1D (:) => Null()
!%!       Integer , Pointer :: ptr_Integer_2D (:,:) => Null()
!%!       Integer , Pointer :: ptr_Integer_3D (:,:,:) => Null()
!%!       Integer , Pointer :: ptr_Integer_4D (:,:,:,:) => Null()
!%!       Integer , Pointer :: ptr_Integer_5D (:,:,:,:,:) => Null()
!%!       End Type Field
!%!       
!%!       Type Array
!%!       private
!%!       Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_1D (:) => Null()
!%!       Real (Kind=ESMF_KIND_R4), Pointer :: ptr_Real4_2D (:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R4), Pointer :: 
!%!      &     ptr_Real4_3D (:,:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R4), Pointer :: 
!%!      &     ptr_Real4_4D (:,:,:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R4), Pointer :: 
!%!      &     ptr_Real4_5D (:,:,:,:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_1D (:) => Null()
!%!       Real (Kind=ESMF_KIND_R8), Pointer :: ptr_Real8_2D (:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R8), Pointer :: 
!%!      &     ptr_Real8_3D (:,:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R8), Pointer :: 
!%!      &     ptr_Real8_4D (:,:,:,:) => Null()
!%!       Real (Kind=ESMF_KIND_R8), Pointer :: 
!%!      &     ptr_Real8_5D (:,:,:,:,:) => Null()
!%!       Integer , Pointer :: ptr_Integer_1D (:) => Null()
!%!       Integer , Pointer :: ptr_Integer_2D (:,:) => Null()
!%!       Integer , Pointer :: ptr_Integer_3D (:,:,:) => Null()
!%!       Integer , Pointer :: ptr_Integer_4D (:,:,:,:) => Null()
!%!       Integer , Pointer :: ptr_Integer_5D (:,:,:,:,:) => Null()
!%!       End Type Array
!%!       
!%!       Type RouteHandle
!%!       Private
!%!       Integer :: cannot_have_empty_types
!%!       End Type RouteHandle
!%!       
!%!       Interface Field_SetDataPointer
!%!       Module Procedure Field_SetDataPointer_Real4_1D
!%!       Module Procedure Field_SetDataPointer_Real4_2D
!%!       Module Procedure Field_SetDataPointer_Real4_3D
!%!       Module Procedure Field_SetDataPointer_Real4_4D
!%!       Module Procedure Field_SetDataPointer_Real4_5D
!%!       Module Procedure Field_SetDataPointer_Real8_1D
!%!       Module Procedure Field_SetDataPointer_Real8_2D
!%!       Module Procedure Field_SetDataPointer_Real8_3D
!%!       Module Procedure Field_SetDataPointer_Real8_4D
!%!       Module Procedure Field_SetDataPointer_Real8_5D
!%!       Module Procedure Field_SetDataPointer_Integer_1D
!%!       Module Procedure Field_SetDataPointer_Integer_2D
!%!       Module Procedure Field_SetDataPointer_Integer_3D
!%!       Module Procedure Field_SetDataPointer_Integer_4D
!%!       Module Procedure Field_SetDataPointer_Integer_5D
!%!       End Interface
!%!       
!%!       
!%!       
!%!       Interface Field_GetDataPointer
!%!       Module Procedure Field_GetDataPointer_Real4_1D
!%!       Module Procedure Field_GetDataPointer_Real4_2D
!%!       Module Procedure Field_GetDataPointer_Real4_3D
!%!       Module Procedure Field_GetDataPointer_Real4_4D
!%!       Module Procedure Field_GetDataPointer_Real4_5D
!%!       Module Procedure Field_GetDataPointer_Real8_1D
!%!       Module Procedure Field_GetDataPointer_Real8_2D
!%!       Module Procedure Field_GetDataPointer_Real8_3D
!%!       Module Procedure Field_GetDataPointer_Real8_4D
!%!       Module Procedure Field_GetDataPointer_Real8_5D
!%!       Module Procedure Field_GetDataPointer_Integer_1D
!%!       Module Procedure Field_GetDataPointer_Integer_2D
!%!       Module Procedure Field_GetDataPointer_Integer_3D
!%!       Module Procedure Field_GetDataPointer_Integer_4D
!%!       Module Procedure Field_GetDataPointer_Integer_5D
!%!       End Interface
!%!       
!%!       Contains
!%!       
!%! ! This routine is a simplified rendition of DataMap. It
!%! ! is used simply to specify which rank of a field is decomposed.
!%!       Subroutine Field_SetDistRank(aField, dist_rank)
!%!       Type (Field), Intent(InOut) :: aField
!%!       Integer, Intent(In) :: dist_rank
!%!       
!%!       aField%dist_rank = dist_rank
!%!       End Subroutine Field_SetDistRank
!%!       
!%! ! This routine is a simplified rendition of DataMap. It
!%! ! is used simply to determine which rank of a field is decomposed.
!%!       Function Field_GetDistRank(aField) Result(dist_rank)
!%!       Type (Field), Intent(In) :: aField
!%!       Integer :: dist_rank
!%!       
!%!       dist_rank = aField%dist_rank
!%!       End Function Field_GetDistRank
!%!       
!%!       Function Field_GetRank(aField) Result(rank)
!%!       Type (Field), Intent(In) :: aField
!%!       Integer :: rank
!%!       
!%!       rank = aField%rank
!%!       End Function Field_GetRank
!%!       
!%!       Subroutine Field_GetShape(aField, data_shape)
!%!       Type (Field), Intent(In) :: aField
!%!       Integer :: data_shape(aField%rank)
!%!       
!%!       data_shape = aField%data_shape(1:aField%rank)
!%!       
!%!       End Subroutine Field_GetShape
!%!       
!%!       Subroutine Field_SetDataPointer_Real4_1D (aField, dataPointer,
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Real (Kind=ESMF_KIND_R4), Pointer       :: dataPointer (:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 1 
!%!       aField% ptr_Real4_1D => dataPointer 
!%!       aField%data_shape(1:1) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Real4_1D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Real4_2D (aField, dataPointer, 
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Real (Kind=ESMF_KIND_R4), Pointer       :: dataPointer (:,:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 2 
!%!       aField% ptr_Real4_2D => dataPointer 
!%!       aField%data_shape(1:2) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Real4_2D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Real4_3D (aField, dataPointer,
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Real (Kind=ESMF_KIND_R4), Pointer       :: dataPointer (:,:,:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 3 
!%!       aField% ptr_Real4_3D => dataPointer 
!%!       aField%data_shape(1:3) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Real4_3D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Real4_4D (aField, dataPointer, 
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Real (Kind=ESMF_KIND_R4), Pointer       :: dataPointer (:,:,:,:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 4 
!%!       aField% ptr_Real4_4D => dataPointer 
!%!       aField%data_shape(1:4) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Real4_4D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Real4_5D (aField, dataPointer, 
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Real (Kind=ESMF_KIND_R4), Pointer       :: dataPointer (:,:,:,:,:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 5 
!%!       aField% ptr_Real4_5D => dataPointer 
!%!       aField%data_shape(1:5) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Real4_5D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Real8_1D (aField, dataPointer, 
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Real (Kind=ESMF_KIND_R8), Pointer       :: dataPointer (:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 1 
!%!       aField% ptr_Real8_1D => dataPointer 
!%!       aField%data_shape(1:1) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Real8_1D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Real8_2D (aField, dataPointer,
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Real (Kind=ESMF_KIND_R8), Pointer       :: dataPointer (:,:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 2 
!%!       aField% ptr_Real8_2D => dataPointer 
!%!       aField%data_shape(1:2) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Real8_2D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Real8_3D (aField, dataPointer,
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Real (Kind=ESMF_KIND_R8), Pointer       :: dataPointer (:,:,:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 3 
!%!       aField% ptr_Real8_3D => dataPointer 
!%!       aField%data_shape(1:3) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Real8_3D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Real8_4D (aField, dataPointer,
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Real (Kind=ESMF_KIND_R8), Pointer       :: dataPointer (:,:,:,:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 4 
!%!       aField% ptr_Real8_4D => dataPointer 
!%!       aField%data_shape(1:4) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Real8_4D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Real8_5D (aField, dataPointer, 
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Real (Kind=ESMF_KIND_R8), Pointer       :: dataPointer (:,:,:,:,:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 5 
!%!       aField% ptr_Real8_5D => dataPointer 
!%!       aField%data_shape(1:5) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Real8_5D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Integer_1D (aField, dataPointer, 
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Integer , Pointer       :: dataPointer (:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 1 
!%!       aField% ptr_Integer_1D => dataPointer 
!%!       aField%data_shape(1:1) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Integer_1D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Integer_2D (aField, dataPointer, 
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Integer , Pointer       :: dataPointer (:,:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 2 
!%!       aField% ptr_Integer_2D => dataPointer 
!%!       aField%data_shape(1:2) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Integer_2D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Integer_3D (aField, dataPointer, 
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Integer , Pointer       :: dataPointer (:,:,:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 3 
!%!       aField% ptr_Integer_3D => dataPointer 
!%!       aField%data_shape(1:3) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Integer_3D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Integer_4D (aField, dataPointer, 
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Integer , Pointer       :: dataPointer (:,:,:,:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 4 
!%!       aField% ptr_Integer_4D => dataPointer 
!%!       aField%data_shape(1:4) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Integer_4D 
!%!       
!%!       
!%!       Subroutine Field_SetDataPointer_Integer_5D (aField, dataPointer, 
!%!      &     haloWidth, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(InOut) :: aField 
!%!       Integer , Pointer       :: dataPointer (:,:,:,:,:) 
!%!       Integer, Optional, Intent(In) :: haloWidth ! not implemened (not needed?) 
!%!       Integer, Optional, Intent(out) :: rc 
!%!       
!%!       aField%rank = 5 
!%!       aField% ptr_Integer_5D => dataPointer 
!%!       aField%data_shape(1:5) = shape(dataPointer) 
!%!       
!%!       afield%active = .true. 
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_SetDataPointer_Integer_5D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Real4_1D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Real (Kind=ESMF_KIND_R4), Pointer           :: ptr (:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Real4_1D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Real4_1D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Real4_2D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Real (Kind=ESMF_KIND_R4), Pointer           :: ptr (:,:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Real4_2D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Real4_2D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Real4_3D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Real (Kind=ESMF_KIND_R4), Pointer           :: ptr (:,:,:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Real4_3D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Real4_3D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Real4_4D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Real (Kind=ESMF_KIND_R4), Pointer           :: ptr (:,:,:,:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Real4_4D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Real4_4D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Real4_5D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Real (Kind=ESMF_KIND_R4), Pointer           :: ptr (:,:,:,:,:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Real4_5D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Real4_5D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Real8_1D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Real (Kind=ESMF_KIND_R8), Pointer           :: ptr (:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Real8_1D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Real8_1D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Real8_2D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Real (Kind=ESMF_KIND_R8), Pointer           :: ptr (:,:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Real8_2D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Real8_2D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Real8_3D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Real (Kind=ESMF_KIND_R8), Pointer           :: ptr (:,:,:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Real8_3D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Real8_3D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Real8_4D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Real (Kind=ESMF_KIND_R8), Pointer           :: ptr (:,:,:,:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Real8_4D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Real8_4D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Real8_5D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Real (Kind=ESMF_KIND_R8), Pointer           :: ptr (:,:,:,:,:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Real8_5D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Real8_5D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Integer_1D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Integer , Pointer           :: ptr (:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Integer_1D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Integer_1D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Integer_2D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Integer , Pointer           :: ptr (:,:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Integer_2D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Integer_2D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Integer_3D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Integer , Pointer           :: ptr (:,:,:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Integer_3D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Integer_3D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Integer_4D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Integer , Pointer           :: ptr (:,:,:,:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Integer_4D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Integer_4D 
!%!       
!%!       
!%!       Subroutine Field_GetDataPointer_Integer_5D (aField, ptr, rc) 
!%!       Implicit None 
!%!       Type (Field), Intent(In) :: aField 
!%!       Integer , Pointer           :: ptr (:,:,:,:,:) 
!%!       integer, intent(out), optional :: rc  
!%!       
!%!       ptr => aField% ptr_Integer_5D 
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS 
!%!       
!%!       End Subroutine Field_GetDataPointer_Integer_5D 
!%!       
!%!       Function Grid_Create(im, jm, lm) Result(aGrid)
!%!       Implicit None
!%!       Integer :: im, jm, lm
!%!       Type (Grid) :: aGrid
!%!       
!%!       aGrid%im = im
!%!       aGrid%jm = jm
!%!       aGrid%lm = lm
!%!       
!%!       End Function Grid_Create
!%!       
!%!       Subroutine Grid_Destroy(agrid, rc)
!%!       Type (Grid) :: agrid
!%!       Integer, Optional, Intent(Out) :: rc
!%!       
!%!       If (present(rc)) rc = Esmf_Success
!%!       End Subroutine Grid_Destroy
!%!       
!%!       Subroutine Grid_Get(aGrid, ai)
!%!       Implicit None
!%!       Type (Grid) :: aGrid
!%!       Type (AxisIndex) :: ai(N_DIMENSIONS)
!%!       
!%! !!$$ ai%min = (/ 1, 1, 1 /)
!%! !!$$ ai%max = (/ aGrid%im, aGrid%jm, aGrid%lm /)
!%!       ai%min = (/ 1, 1, 1 /)
!%!       ai%max = (/ aGrid%im, aGrid%jm, aGrid%lm /)
!%!       
!%!       End Subroutine Grid_Get
!%!       
!%! !< Method: [[Field_HaloStore]] >>
!%! !< Method: [[Field_Halo]] >>
!%!       
!%!       Subroutine Field_GetArray(aField, anArray, rc)
!%!       Type (Field) :: aField
!%!       Type (Array) :: anArray
!%!       Integer, Optional, Intent(Out) :: rc
!%!       
!%! ! Not implemented
!%!       If (Present(rc)) rc = ESMF_FAILURE
!%!       
!%!       End Subroutine Field_GetArray
!%!       
!%!       Function Field_Clone(f1) Result(f2)
!%!       Type (Field), Intent(In), Target :: f1
!%!       Type (Field) :: f2
!%!       Integer, Pointer :: s(:)
!%!       
!%!       f2%active = f1%active
!%!       f2%esmf_kind = f1%esmf_kind
!%!       f2%esmf_type = f1%esmf_type
!%!       f2%rank = f1%rank
!%!       f2%dist_rank = f1%dist_rank
!%!       f2%data_shape = f1%data_shape
!%!       
!%!       s => f1%data_shape(:f1%rank)
!%!       
!%!       If (Associated(f1%ptr_Real4_1D)) Then 
!%!          Allocate(f2%ptr_Real4_1D(s(1))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Real4_2D)) Then 
!%!          Allocate(f2%ptr_Real4_2D(s(1),s(2))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Real4_3D)) Then 
!%!          Allocate(f2%ptr_Real4_3D(s(1),s(2),s(3))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Real4_4D)) Then 
!%!          Allocate(f2%ptr_Real4_4D(s(1),s(2),s(3),s(4))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Real4_5D)) Then 
!%!          Allocate(f2%ptr_Real4_5D(s(1),s(2),s(3),s(4),s(5))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Real8_1D)) Then 
!%!          Allocate(f2%ptr_Real8_1D(s(1))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Real8_2D)) Then 
!%!          Allocate(f2%ptr_Real8_2D(s(1),s(2))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Real8_3D)) Then 
!%!          Allocate(f2%ptr_Real8_3D(s(1),s(2),s(3))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Real8_4D)) Then 
!%!          Allocate(f2%ptr_Real8_4D(s(1),s(2),s(3),s(4))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Real8_5D)) Then 
!%!          Allocate(f2%ptr_Real8_5D(s(1),s(2),s(3),s(4),s(5))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Integer_1D)) Then 
!%!          Allocate(f2%ptr_Integer_1D(s(1))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Integer_2D)) Then 
!%!          Allocate(f2%ptr_Integer_2D(s(1),s(2))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Integer_3D)) Then 
!%!          Allocate(f2%ptr_Integer_3D(s(1),s(2),s(3))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Integer_4D)) Then 
!%!          Allocate(f2%ptr_Integer_4D(s(1),s(2),s(3),s(4))) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f1%ptr_Integer_5D)) Then 
!%!          Allocate(f2%ptr_Integer_5D(s(1),s(2),s(3),s(4),s(5))) 
!%!       End If 
!%!       
!%!       End Function Field_Clone
!%!       
!%!       Subroutine Field_Destroy(f, rc)
!%!       Type (Field), Intent(InOut), Target :: f
!%!       Integer, Optional, Intent(Out) :: rc
!%!       
!%!       f%active = .false.
!%!       f%rank = 0
!%!       f%dist_rank = 0
!%!       f%data_shape = 0
!%!       
!%!       If (Present(rc)) rc = ESMF_SUCCESS
!%!       
!%!       If (Associated(f%ptr_Real4_1D)) Then 
!%!          Deallocate(f%ptr_Real4_1D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Real4_2D)) Then 
!%!          Deallocate(f%ptr_Real4_2D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Real4_3D)) Then 
!%!          Deallocate(f%ptr_Real4_3D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Real4_4D)) Then 
!%!          Deallocate(f%ptr_Real4_4D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Real4_5D)) Then 
!%!          Deallocate(f%ptr_Real4_5D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Real8_1D)) Then 
!%!          Deallocate(f%ptr_Real8_1D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Real8_2D)) Then 
!%!          Deallocate(f%ptr_Real8_2D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Real8_3D)) Then 
!%!          Deallocate(f%ptr_Real8_3D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Real8_4D)) Then 
!%!          Deallocate(f%ptr_Real8_4D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Real8_5D)) Then 
!%!          Deallocate(f%ptr_Real8_5D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Integer_1D)) Then 
!%!          Deallocate(f%ptr_Integer_1D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Integer_2D)) Then 
!%!          Deallocate(f%ptr_Integer_2D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Integer_3D)) Then 
!%!          Deallocate(f%ptr_Integer_3D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Integer_4D)) Then 
!%!          Deallocate(f%ptr_Integer_4D) 
!%!       End If 
!%!       
!%!       
!%!       If (Associated(f%ptr_Integer_5D)) Then 
!%!          Deallocate(f%ptr_Integer_5D) 
!%!       End If 
!%!       
!%!       End Subroutine Field_Destroy
!%!       
!%!       Logical Function ESMF_type_eq(t1, t2)
!%!       Type (ESMF_DATATYPE), Intent(In) :: t1, t2
!%!       ESMF_type_eq = (t1%type == t2%type)
!%!       End Function ESMF_type_eq
!%!       
!%!       Logical Function ESMF_kind_eq(t1, t2)
!%!       Type (ESMF_DATAKIND), Intent(In) :: t1, t2
!%!       ESMF_kind_eq = (t1%kind == t2%kind)
!%!       End Function ESMF_kind_eq
!%!       
!%!       End Module ESMF_MOD_private
!%!       
!%!       Module ESMF_MOD
!%!       Use ESMF_MOD_private, Only: ESMF_Grid => Grid
!%!       Use ESMF_MOD_private, Only: ESMF_Field => Field
!%!       Use ESMF_MOD_private, Only: ESMF_Array => Array
!%!       Use ESMF_MOD_private, Only: ESMF_AxisIndex => AxisIndex
!%!       Use ESMF_MOD_private, Only: ESMF_KIND_R4
!%!       Use ESMF_MOD_private, Only: ESMF_KIND_R8
!%!       
!%!       Use ESMF_MOD_private, Only: ESMF_FAILURE, ESMF_SUCCESS
!%!       Use ESMF_MOD_private, Only: ESMF_CELL_NFACE, ESMF_CELL_SFACE
!%!       Use ESMF_MOD_private, Only: ESMF_CELL_CENTER, ESMF_CELL_SWCORNER
!%!       
!%!       Use ESMF_MOD_private, Only: ESMF_MAXSTR => MAXSTR
!%!       
!%!       Use ESMF_MOD_private, Only: 
!%!      &     ESMF_FieldGetDataPointer => Field_GetDataPointer
!%!       Use ESMF_MOD_private, Only: ESMF_FieldGetArray => Field_GetArray
!%!       Use ESMF_MOD_private, Only: ESMF_FieldDestroy => Field_Destroy
!%!       
!%!       Use ESMF_MOD_private, Only: ESMF_DATATYPE, ESMF_DATAKIND
!%!       Use ESMF_MOD_private, Only: ESMF_DATA_INTEGER, ESMF_DATA_REAL
!%!       Use ESMF_MOD_private, Only: ESMF_R8, ESMF_R4
!%!       Use ESMF_MOD_private, Only: ESMF_I8, ESMF_I4
!%!       
!%!       Use ESMF_Mod_private, Only: Operator(==)
!%!       
!%!       Use ESMF_Mod_private, Only: ESMF_GridCreate => Grid_Create
!%!       Use ESMF_Mod_private, Only: ESMF_GridDestroy => Grid_Destroy
!%!       Use ESMF_MOD_private, Only: 
!%!      &     ESMF_FieldSetDataPointer => Field_SetDataPointer
!%!       
!%!       Implicit None
!%!       Private
!%!       
!%!       Public :: ESMF_Grid
!%!       Public :: ESMF_Array
!%!       Public :: ESMF_Field
!%!       Public :: ESMF_AxisIndex
!%!       
!%!       Public :: ESMF_KIND_R4
!%!       Public :: ESMF_KIND_R8
!%!       
!%!       Public :: ESMF_FieldGetDataPointer
!%!       Public :: ESMF_FieldSetDataPointer
!%!       Public :: ESMF_FieldGetArray
!%!       Public :: ESMF_FieldDestroy
!%!       
!%!       Public :: ESMF_FAILURE, ESMF_SUCCESS
!%!       
!%!       Public :: ESMF_CELL_NFACE
!%!       Public :: ESMF_CELL_SFACE
!%!       Public :: ESMF_CELL_CENTER
!%!       Public :: ESMF_CELL_SWCORNER
!%!       
!%!       Public :: ESMF_MAXSTR
!%!       
!%!       Public :: ESMF_DATATYPE, ESMF_DATAKIND
!%!       Public :: ESMF_DATA_INTEGER, ESMF_DATA_REAL
!%!       Public :: ESMF_R8, ESMF_R4
!%!       Public :: ESMF_I8, ESMF_I4
!%!       
!%!       Public :: Operator(==)
!%!       
!%!       Public :: ESMF_Initialize, ESMF_Finalize
!%!       Public :: ESMF_GridCreate, ESMF_GridDestroy
!%!       
!%!       Contains
!%!       
!%!       Subroutine ESMF_Initialize(rc)
!%!       Integer, Optional, Intent(Out) :: rc
!%!       If (Present(rc)) rc = ESMF_SUCCESS
!%!       End Subroutine ESMF_Initialize
!%!       
!%!       Subroutine ESMF_Finalize(rc)
!%!       Integer, Optional, Intent(Out) :: rc
!%!       If (Present(rc)) rc = ESMF_SUCCESS
!%!       End Subroutine ESMF_Finalize
!%!       
!%!       End Module ESMF_MOD

!==============================================================================
! LTM: DOMAIN_DECOMPcs.f
!==============================================================================

!%!       module domain_decomp_atm
!%! !@auth M. Kelley
!%! !@ver 1.0
!%! !@sum  This module bundles the procedures necessary to a modelE
!%! !@+    cubed-sphere domain decomposition into a single module,
!%! !@+    and provides a version of init_grid specific to a cubed
!%! !@+    sphere atmosphere that runs the GSFC FVcubed dynamical core.
!%! 
!%! c get grid-independent procedures from domain_decomp_1d
!%!       use domain_decomp_1d, only : get, am_i_root, ! sumxpe,
!%! !     &     read_parallel,write_parallel,esmf_bcast,
!%! !     &     globalmax, setMpiCommunicator,
!%!      &     hasSouthPole, hasNorthPole
!%! 
!%! c get dist_grid, halo_update, globalsum, etc. from the dd2d_utils module
!%!       use dd2d_utils, only : dist_grid,init_dist_grid
!%! !     &                      ,halo_update !,globalsum
!%! c the following are also available, but use of global arrays is discouraged
!%! c     &                     ,pack_data,unpack_data
!%! 
!%! c get i/o routines for fortran binary sequential access
!%! !      use pario_fbsa  ! readt_parallel et al.
!%! 
!%! c get i/o routines for netcdf format
!%! !      use pario
!%! 
!%! c get procedures from the MPP package as they become compatible with
!%! c model E or vice versa
!%! c      use mpp_domains_mod
!%! 
!%!       USE ESMF_MOD, only : ESMF_Grid ! ESMF_VM
!%! 
!%!       implicit none
!%!       save
!%! 
!%! c declare an instance of dist_grid for the atmosphere
!%!       type(dist_grid) :: grid
!%! 
!%! c declare requisite instances of these types for the atmosphere
!%!       !Type (ESMF_VM) :: modelE_vm
!%!       TYPE (ESMF_Grid) :: ESMF_GRID_ATM
!%! 
!%!       contains
!%! 
!%!       SUBROUTINE INIT_GRID(grd_dum, IM,JM,LM,
!%!      &     width,J_SCM,bc_periodic,CREATE_CAP) ! optional arguments
!%!       USE ESMF_MOD
!%!       IMPLICIT NONE
!%!       TYPE (DIST_GRID), INTENT(INOUT) :: grd_dum
!%!       INTEGER, INTENT(IN) :: IM,JM,LM
!%!       INTEGER, OPTIONAL :: width
!%!       INTEGER, OPTIONAL, INTENT(IN) :: J_SCM       ! ignored here
!%!       LOGICAL, OPTIONAL, INTENT(IN) :: bc_periodic ! ignored here
!%!       LOGICAL, OPTIONAL, INTENT(IN) :: CREATE_CAP  ! always true for FVcubed
!%! c local variables:
!%!       integer :: rc,width_
!%!       TYPE(ESMF_Grid), external :: AppGridCreateF
!%!       !TYPE(ESMF_Config) :: cf
!%!       !Type (ESMF_DELayout)::layout
!%!       integer :: npes,npesx,npesy,  rank_glob,rank_tile,rank_i,rank_j
!%!       integer, parameter :: npesx_max=100
!%!       integer, dimension(npesx_max) :: istrt,istop
!%! c      integer :: ntiles, npx, npy, ng, isd, ied , jsd, jed ! for MPP
!%!       include 'mpif.h' ! temporary to see MPI_COMM_WORLD
!%! 
!%! c sanity checks
!%!       if(jm.ne.im)
!%!      &     call stop_model('cubed sphere init_grid: im != jm',255)
!%!       if(mod(im,2).ne.0)
!%!      &     call stop_model('cubed sphere init_grid: odd im',255)
!%! 
!%!       width_ = 1
!%!       If (Present(width)) width_=width
!%! 
!%!       grd_dum%private%PERIODICBC = .false.
!%!       grd_dum%private%hasSouthPole = .false.
!%!       grd_dum%private%hasNorthPole = .false.
!%!       grd_dum%private%hasEquator = .false.
!%! 
!%! !      Call ESMF_Initialize(vm=modelE_vm, rc=rc)
!%! !      Call ESMF_VMGet(modelE_vm, petCount=NPES, localPET=rank_glob)
!%! !
!%!       npesx = int(floor(sqrt(real(NPES/6))))
!%!       npesy = NPES / npesx  ! in ESMF, npesy = 6*npesx
!%! 
!%! c FVcubed uses a domain decomposition decided by the MPP library.
!%! c To _temporarily_ avoid build dependencies upon MPP, mimic
!%! c its procedure for grouping IM gridcells into NPESX divisions.
!%! c Store start/end indices for each division in istrt,istop.
!%!       call mimic_mpp_division_of_points(im,npesx,istrt,istop)
!%! 
!%!       rank_tile = mod(rank_glob,npes/6)
!%!       rank_i    = mod(rank_tile,npesx)
!%!       rank_j    = rank_tile/npesx
!%!       grd_dum%I_STRT        = istrt(1+rank_i)
!%!       grd_dum%I_STOP        = istop(1+rank_i)
!%!       grd_dum%J_STRT        = istrt(1+rank_j)
!%!       grd_dum%J_STOP        = istop(1+rank_j)
!%! 
!%! c
!%! c Once the MPP communication procedures start to be used, we will
!%! c initialize its domain2D object and retrieve I/J_STRT/STOP from it
!%! c
!%! c      npx=IM+1; npy=JM+1; ng = width_;
!%! c      ntiles=6
!%! c      my_pet = mpp_pe()
!%! c      NPES   = mpp_npes()
!%! c      root   = mpp_root_pe()
!%! c      call init_domain(grd_dum%domain,npx,npy,ntiles,ng,
!%! c     &                 grd_dum%bc_periodic)
!%! c      call mpp_get_compute_domain( grd_dum%domain,
!%! c     &                         I0_DUM, I1_DUM, J0_DUM, J1_DUM)
!%! c      call mpp_get_data_domain( grd_dum%domain, isd, ied , jsd, jed)
!%! c      write(*,*)'mpp-bounds',my_pet,I0_DUM, I1_DUM, J0_DUM, J1_DUM
!%! 
!%! 
!%!       grd_dum%IM_WORLD      = IM
!%!       grd_dum%JM_WORLD      = IM
!%! 
!%!       grd_dum%I_STRT_HALO   = grd_dum%I_STRT - width_
!%!       grd_dum%I_STOP_HALO   = grd_dum%I_STOP + width_
!%!       grd_dum%J_STRT_HALO   = grd_dum%J_STRT - width_
!%!       grd_dum%J_STOP_HALO   = grd_dum%J_STOP + width_
!%! 
!%!       grd_dum%J_STRT_SKP    = grd_dum%J_STRT 
!%!       grd_dum%J_STOP_SKP    = grd_dum%J_STOP
!%! 
!%! c set these to something huge to trigger segfaults if accidentlly used?
!%!       grd_dum%J_STRT_STGR   = grd_dum%J_STRT
!%!       grd_dum%J_STOP_STGR   = grd_dum%J_STOP
!%! 
!%! !     call setMpiCommunicator(grd_dum, MPI_COMM_WORLD)
!%! 
!%! c
!%! c initialize compoments of dist_grid specific to dd2d_utils routines
!%! c
!%!       call init_dist_grid(
!%!      &     grd_dum%IM_WORLD,grd_dum%JM_WORLD,6,
!%!      &     grd_dum%I_STRT,grd_dum%I_STOP,
!%!      &     grd_dum%J_STRT,grd_dum%J_STOP,
!%!      &     grd_dum%I_STRT_HALO,grd_dum%I_STOP_HALO,
!%!      &     grd_dum%J_STRT_HALO,grd_dum%J_STOP_HALO,grd_dum)
!%! 
!%!       grd_dum%have_domain = .true.
!%!       !grd_dum%mpi_comm = MPI_COMM_WORLD
!%!       !grd_dum%npes_comm = grd_dum%nproc
!%! 
!%! !c
!%! !c create the configuration file needed by the FVcubed core
!%! !c
!%! !      cf = load_cap_config('cap.rc',IM,JM*6,LM,npesx,npesy)
!%! !      print*, 'Started AppGridCreateF'
!%! !      ESMF_GRID_ATM = AppGridCreateF(cf, modelE_vm, rc)
!%! !c is this needed?
!%! !      call ESMF_GridGet(ESMF_GRID_ATM, delayout=layout, rc=rc)
!%! !      print*, 'Finished AppGridCreateF'
!%! 
!%!       return
!%!       END SUBROUTINE INIT_GRID
!%! 
!%!       subroutine mimic_mpp_division_of_points(npx,ndiv,is,ie)
!%! c Divides 1:npx gridpoints into ndiv intervals.
!%! c The start,end indices of interval m are stored in is(m),ie(m).
!%! c Intended only for even values of npx!
!%!       implicit none
!%!       integer :: npx,ndiv ! inputs
!%!       integer, dimension(ndiv) :: is,ie ! outputs
!%!       integer :: idiv,ni,jdiv
!%!       real*4 :: rhalf
!%!       rhalf = .5*real(ndiv)
!%!       is(1) = 1
!%!       ie(ndiv) = npx
!%!       do idiv=1,ndiv/2      ! loop over half of the intervals
!%!         jdiv = ndiv-idiv+1  ! and set the other half using symmetry
!%!         if(idiv.gt.1) then
!%!           is(idiv) = ie(idiv-1) + 1
!%!           ie(jdiv) = is(jdiv+1) - 1
!%!         endif
!%! ! the following formula for ni is based on the one in MPP:
!%! ! domain is sized by dividing the remaining points by remaining domains
!%! ! CEILING causes the larger intervals to lie at the edges
!%!         ni = CEILING( real(npx/2-is(idiv)+1)/(rhalf-real(idiv-1)) )
!%!         ie(idiv) = is(idiv) + ni - 1
!%!         is(jdiv) = ie(jdiv) - ni + 1
!%!       enddo
!%!       if(ndiv.gt.1 .and. mod(ndiv,2).eq.1) then
!%! c ndiv odd: the middle interval
!%!         idiv = 1+ndiv/2
!%!         is(idiv) = ie(idiv-1) + 1
!%!         ie(idiv) = is(idiv+1) - 1
!%!       endif
!%!       return
!%!       end subroutine mimic_mpp_division_of_points
!%! 
!%! !! ----------------------------------------------------------------------
!%! !      function load_cap_config(config_file,IM,JM,LM,NP_X,NP_Y)
!%! !     &     result( config )
!%! !! ----------------------------------------------------------------------
!%! !      use ESMF_mod
!%! !      !use FILEMANAGER    
!%! !      include 'mpif.h' ! temporary to see MPI_Barrier
!%! !      character(len=*), parameter :: Iam=
!%! !     &     "DOMAIN_DECOMP::load_cap_config"
!%! !      character(len=*), intent(in) :: config_file
!%! !      integer,          intent(in) :: IM,JM,LM,NP_X,NP_Y
!%! !      type (esmf_config)           :: config
!%! !    
!%! !      integer :: rc, iunit
!%! !    
!%! !      config = ESMF_ConfigCreate(rc=rc)
!%! !      !call openunit(config_file, iunit, qbin=.false., qold=.false.)
!%! !      !if(am_i_root()) then
!%! !      !  write(iunit,*)'IM:  ', IM
!%! !      !  write(iunit,*)'JM:  ', JM
!%! !      !  write(iunit,*)'LM:  ', LM
!%! !      !  write(iunit,*)'NX:  ', NP_X
!%! !      !  write(iunit,*)'NY:  ', NP_Y
!%! !      !endif
!%! !      !call closeUnit(iunit)
!%! !     
!%! !      call MPI_Barrier(mpi_comm_world, rc)
!%! !      call ESMF_ConfigLoadFile(config, config_file, rc=rc)
!%! !
!%! !      end function load_cap_config
!%! 
!%!       end module domain_decomp_atm

!==============================================================================
! LTM: ATMDYN_COM.f
!==============================================================================

      MODULE DYNAMICS
!@sum  DYNAMICS contains all the pressure and momentum related variables
!@auth Original development team
!@ver  1.0
      USE DOMAIN_DECOMP_1D, ONLY : grid
      USE RESOLUTION, ONLY : LM
      IMPLICIT NONE
      SAVE
C**** Some helpful arrays (arrays should be L first)
!@var  PLIJ  Surface pressure: P(I,J) or PSF-PTOP (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PLIJ
!@var  PDSIG  Surface pressure * DSIG(L) (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PDSIG
!@var  AM  Air mass of each box (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: AM     ! PLIJ*DSIG(L)*100/grav
!@var  BYAM  1/Air mass (m^2/kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: BYAM
!@var  PMID  Pressure at mid point of box (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PMID    ! SIG(L)*PLIJ+PTOP
!@var  PK   PMID**KAPA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PK
!@var  PEUP  Pressure at lower edge of box (incl. surface) (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PEDN  ! SIGE(L)*PLIJ+PTOP
!@var  PEK  PEUP**KAPA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PEK
!@var  SQRT  square root of P (used in diagnostics)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SQRTP
!@var  PTROPO  Pressure at mid point of tropopause level (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PTROPO
!@var  LTROPO  Tropopause layer
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LTROPO

C**** module should own dynam variables used by other routines
!@var PTOLD pressure at beginning of dynamic time step (for clouds)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)    :: PTOLD
!@var SD_CLOUDS vert. integrated horizontal convergence (for clouds)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SD_CLOUDS
!@var GZ geopotential height (for Clouds and Diagnostics)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GZ
!@var DPDX_BY_RHO,DPDY_BY_RHO (pressure gradients)/density at L=1
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDX_BY_RHO,DPDY_BY_RHO
!@var DPDX_BY_RHO_0,DPDY_BY_RHO_0 surface (pressure gradients)/density
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDX_BY_RHO_0
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDY_BY_RHO_0

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PU,PV
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), TARGET :: CONV
cgsfc      REAL*8, DIMENSION(IM,JM,LM-1) :: SD
      REAL*8, ALLOCATABLE :: SD(:,:,:)
!@var PIT  pressure tendency (mb m^2/s)
      REAL*8, ALLOCATABLE :: PIT(:,:)
cgsfc      EQUIVALENCE (SD(1,1,1),CONV(1,1,2))

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PHI,SPA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DUT,DVT
!!@var xAVRX scheme-depend. coefficient for AVRX: 1,byrt2 (2nd,4th order)
!      REAL*8 xAVRX

!@var PUA,PVA,SDA,PS save PU,PV,SD,P for hourly tracer advection
!@var MB Air mass array for tracers (before advection)
!@var MA Air mass array for tracers (updated during advection)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PUA,PVA,SDA,MB,MA
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PS

!@var WCP vertical mass flux in a constant-pressure vertical
!@+   coordinate whose pressure levels are the global means of
!@+   each layer.
!@var WCPsig: WCP interpolated to terrain-following coordinate
!@+   surfaces (sigma levels)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: WCP,WCPsig

!@var DKE change in KE due to dissipation (SURF/DC/MC) (m^2/s^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DKE
!@var KEA KE on the A grid (m^2/s^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: KEA ! ke on A grid
!@var UALIJ,VALIJ U,V on the A grid (m/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: UALIJ,VALIJ
!@var WSAVE vertical velocity (m/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: WSAVE
!@var SMASS local but "SAVE"d array ADVECV in MOMEN2ND made global
!@    here since its use does not go beyond ATMDYN that calls ADVECV
      REAL*8, ALLOCATABLE:: SMASS(:)

      INTEGER, DIMENSION(2) :: t_dyn_a, t_dyn_b, t_dyn_c, t_dyn_d
      INTEGER, DIMENSION(2) :: t_aflux, t_advecm, t_advecv
      INTEGER, DIMENSION(2) :: t_PGF, t_filter, t_fltruv, t_calc_pijl
      INTEGER, DIMENSION(2) :: t_sdrag, t_calc_trop

      REAL*8, PARAMETER :: COS_LIMIT = 0.15d0

      END MODULE DYNAMICS

      SUBROUTINE ALLOC_DYNAMICS(grid)
      USE DOMAIN_DECOMP_1d, ONLY : DIST_GRID
      USE RESOLUTION , ONLY : LM
      USE DYNAMICS, ONLY : PLIJ,PDSIG,AM,BYAM,PMID,PK,PEDN,PEK,
     $                     SD_CLOUDS,GZ,PU,PV,CONV,PHI,SPA,DUT,
     $                     DVT,PUA,PVA,SDA,MB,MA,DKE,WSAVE,SD,PIT,
     $                     SQRTP,PTROPO,LTROPO,PTOLD,DPDX_BY_RHO,
     $                     DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0,
     $                     PS,SMASS,KEA,UALIJ,VALIJ,WCP,WCPsig
      USE DYNAMICS, only : t_dyn_a, t_dyn_b, t_dyn_c, t_dyn_d
      USE DYNAMICS, only : t_aflux, t_advecm, t_advecv

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO
      
      ! K-I-J arrays
      ALLOCATE ( PLIJ(LM,I_0H:I_1H,J_0H:J_1H), 
     $           PDSIG(LM,I_0H:I_1H,J_0H:J_1H),
     $             AM(LM,I_0H:I_1H,J_0H:J_1H),  
     $           BYAM(LM,I_0H:I_1H,J_0H:J_1H),
     $           PMID(LM,I_0H:I_1H,J_0H:J_1H),    
     $           PK(LM,I_0H:I_1H,J_0H:J_1H),
     $           PEDN(LM+1,I_0H:I_1H,J_0H:J_1H), 
     $            PEK(LM+1,I_0H:I_1H,J_0H:J_1H),
     $   STAT = IER)
      ! I-J-K arrays
      ALLOCATE( SD_CLOUDS(I_0H:I_1H,J_0H:J_1H,LM),  
     $                 GZ(I_0H:I_1H,J_0H:J_1H,LM), 
     $                 PU(I_0H:I_1H,J_0H:J_1H,LM),  
     $                 PV(I_0H:I_1H,J_0H:J_1H,LM), 
     $               CONV(I_0H:I_1H,J_0H:J_1H,LM), 
     $               PHI(I_0H:I_1H,J_0H:J_1H,LM), 
     $                SPA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                DUT(I_0H:I_1H,J_0H:J_1H,LM), 
     $                DVT(I_0H:I_1H,J_0H:J_1H,LM), 
     $                PUA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                PVA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                SDA(I_0H:I_1H,J_0H:J_1H,LM),  
     $                MB(I_0H:I_1H,J_0H:J_1H,LM), 
     $                MA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                DKE(I_0H:I_1H,J_0H:J_1H,LM), 
     $                KEA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                UALIJ(LM,I_0H:I_1H,J_0H:J_1H), 
     $                VALIJ(LM,I_0H:I_1H,J_0H:J_1H), 
     $              WSAVE(I_0H:I_1H,J_0H:J_1H,LM-1), 
     $              PIT(I_0H:I_1H,J_0H:J_1H),
     $              SD(I_0H:I_1H,J_0H:J_1H,LM-1),
     $   STAT = IER)

      !hack to remove NaNs from "SD" - possible mistake somewhere
      CONV(:,:,:) = 0.d0

      ! I-J arrays
      ALLOCATE(  SQRTP(I_0H:I_1H,J_0H:J_1H), 
     $          PTROPO(I_0H:I_1H,J_0H:J_1H),
     $          LTROPO(I_0H:I_1H,J_0H:J_1H),  
     $          PTOLD(I_0H:I_1H,J_0H:J_1H),
     $          DPDX_BY_RHO(I_0H:I_1H,J_0H:J_1H), 
     $          DPDY_BY_RHO(I_0H:I_1H,J_0H:J_1H),
     $          DPDX_BY_RHO_0(I_0H:I_1H,J_0H:J_1H), 
     $          DPDY_BY_RHO_0(I_0H:I_1H,J_0H:J_1H),
     $          PS(I_0H:I_1H,J_0H:J_1H),
     $          WCP(I_0H:I_1H,J_0H:J_1H,LM),
     $          WCPsig(I_0H:I_1H,J_0H:J_1H,LM),
     $   STAT = IER)

      ! J arrays
      ALLOCATE(  SMASS(J_0H:J_1H), 
     $   STAT = IER)

! correct or wrong, but being static all arrays were initialized
! to zero by default. They have to be initialized to something now 
! to avoid floating point exceptions...
      DPDX_BY_RHO(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDY_BY_RHO(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDX_BY_RHO_0(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDY_BY_RHO_0(I_0H:I_1H,J_0H:J_1H) = 0.d0

      PU(I_0H:I_1H,J_0H:J_1H,1:LM) = 0.d0
      PV(I_0H:I_1H,J_0H:J_1H,1:LM) = 0.d0

      SD_CLOUDS(I_0H:I_1H,J_0H:J_1H,1:LM) = 0.d0
      PIT(:,:) = 0.d0

      END SUBROUTINE ALLOC_DYNAMICS

      SUBROUTINE ALLOC_MODEL_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      !USE DOMAIN_DECOMP_ATM, ONLY : GRID, DIST_GRID
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      USE RESOLUTION, ONLY : IM,JM,LM
      USE MODEL_COM, ONLY : NTYPE
      USE MODEL_COM, ONLY : ZATMO,HLAKE,FLAND,FOCEAN,FLICE,FLAKE0,
     *                      FEARTH0,WFCS,P,U,V,T,Q,WM,FTYPE
#ifdef BLK_2MOM
     *  ,WMICE
#endif
!AOO      USE ESMF_CUSTOM_MOD, ONLY: modelE_grid
!AOO      USE ESMF_CUSTOM_MOD, ONLY: ESMF_CELL_SFACE
!AOO      USE ESMF_CUSTOM_MOD, ONLY: ESMF_CELL_CENTER

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE(ZATMO(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(HLAKE(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(FLAND(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(FOCEAN(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(FLICE(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(FLAKE0(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(FEARTH0(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(WFCS(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(P(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(U(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)

      ALLOCATE(V(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)

      ALLOCATE(T(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)

      ALLOCATE(Q(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)

      ALLOCATE(WM(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)
#ifdef BLK_2MOM
      ALLOCATE(WMICE(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)
#endif


      ALLOCATE(FTYPE(NTYPE,I_0H:I_1H,J_0H:J_1H), STAT = IER)

! initialize non-initialized arrays
      FTYPE(:,:,:) = 0.d0

      END SUBROUTINE ALLOC_MODEL_COM

      SUBROUTINE CALC_VERT_AMP(P0,LMAX,PL,AM,PDSIG,PEDN,PMID)
!@sum  CALC_VERT_AMPK calculates air mass and pressure vertical arrays
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : bygrav
      USE MODEL_COM, only : lm,ls1,dsig,sig,sige,ptop,psfmpt,lm_req
     *     ,req_fac,req_fac_m,req_fac_d,pmtop
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: P0 !@var P0 surface pressure (-PTOP) (mb)
      INTEGER, INTENT(IN) :: LMAX !@var LMAX max level for calculation
!@var AM mass at each level (kg/m2)
!@var PDSIG pressure interval at each level (mb)
!@var PMID mid-point pressure (mb)
      REAL*8, INTENT(OUT), DIMENSION(LMAX) :: AM,PDSIG,PMID,PL
!@var PEDN edge pressure (top of box) (mb)
      REAL*8, INTENT(OUT), DIMENSION(LMAX+1) :: PEDN
      INTEGER :: L  !@var L  loop variables

C**** Calculate air mass, layer pressures
C**** Note that only layers LS1 and below vary as a function of surface
C**** pressure.
C**** Note Air mass is calculated in (kg/m^2)

      DO L=1,LS1-1
        PL(L)   = P0
        PDSIG(L)= P0*DSIG(L)
        PMID(L) = SIG(L)*P0+PTOP
        PEDN(L) = SIGE(L)*P0+PTOP
        AM  (L) = PDSIG(L)*1d2*BYGRAV
      END DO
      DO L=LS1,MIN(LMAX,LM)
        PL(L)   = PSFMPT
        PDSIG(L)= PSFMPT*DSIG(L)
        PMID(L) = SIG(L)*PSFMPT+PTOP
        PEDN(L) = SIGE(L)*PSFMPT+PTOP
        AM  (L) = PDSIG(L)*1d2*BYGRAV
      END DO
      IF (LMAX.ge.LM) PEDN(LM+1) = SIGE(LM+1)*PSFMPT+PTOP
C**** Rad. equ. layers if necessary (only PEDN,AM,PMID)
      IF (LMAX.eq.LM+LM_REQ) THEN
        PMID(LM+1:LM+LM_REQ) = REQ_FAC_M(1:LM_REQ)*PMTOP
          AM(LM+1:LM+LM_REQ) = REQ_FAC_D(1:LM_REQ)*PMTOP*1d2*BYGRAV
        PEDN(LM+2:LM+LM_REQ) = REQ_FAC(1:LM_REQ-1)*PEDN(LM+1)
        PEDN(LM+LM_REQ+1)=0.    ! 1d-5  ! why not zero?
      END IF

      RETURN
      END SUBROUTINE CALC_VERT_AMP

!==============================================================================
! LTM: GEOM_B.f
!==============================================================================

      MODULE GEOM
!@sum  GEOM contains spherical geometric variables and arrays
!@auth Original development team
!@ver  1.0 (B grid version)
!@cont GEOM_B
      USE CONSTANT, only : OMEGA,RADIUS,TWOPI,SDAY,radian,AREAG
      USE MODEL_COM, only : IM,JM,LM,FIM,BYIM
      IMPLICIT NONE
      PRIVATE

      public :: GEOM_B, areag
      public :: lonlat_to_ij
      public :: lonlat_to_tile
      public :: lon_to_I
      public :: lat_to_J

C**** The primary grid is the A grid (including both poles)
C**** The secondary grid is for the B grid velocities, located on the
C**** vertices of the A grid (Note: first velocity point is J=2)
C**** Polar boxes can have different latitudinal size and are treated
C**** as though they were 1/IM of their actual area
      SAVE
!@param  DLON grid spacing in longitude (deg)
      REAL*8, PUBLIC, PARAMETER :: DLON = TWOPI*BYIM
C**** For the wonderland model set DLON=DLON/3
c      REAL*8, PUBLIC, PARAMETER :: DLON=TWOPI/(IM*3)
!@var DLAT,DLAT_DG,DLATM grid spacing in latitude (rad,deg,minutes)
      REAL*8, PUBLIC  :: DLAT,DLAT_DG, DLATM
!@param FJEQ equatorial value of J
      REAL*8, PUBLIC, PARAMETER :: FJEQ=.5*(1+JM)
!@var  J1U index of southernmost latitude (currently 2, later 1)
      INTEGER, parameter :: J1U = 2
!@var  JRANGE_HEMI lowest,highest lat index for SH,NH for A,B grid
      INTEGER, PUBLIC, parameter, dimension(2,2,2) :: JRANGE_HEMI = 
     *     reshape(
     *  (/1,JM/2,  1+JM/2,JM,  J1U,J1U-1+JM/2, J1U-1+JM/2,JM+J1U-2/),
     *  (/2,2,2/))
!@var  LAT latitude of mid point of primary grid box (radians)
!@var  LATV latitude of southern edge of primary grid box (radians)
      REAL*8, PUBLIC, DIMENSION(JM) :: LAT,LATV
!@var  LAT_DG latitude of mid points of primary and sec. grid boxs (deg)
      REAL*8, PUBLIC, DIMENSION(JM,2) :: LAT_DG
!@var  LON longitude of mid points of primary grid box (radians)
!@var  LONV longitude of east edge of primary grid box (radians)
      REAL*8, PUBLIC, DIMENSION(IM) :: LON,LONV
!@var  LON_DG longitude of mid points of prim. and sec. grid boxes (deg)
      REAL*8, PUBLIC, DIMENSION(IM,2) :: LON_DG
!@var  DXYP,BYDXYP area of grid box (+inverse) (m^2)
C**** Note that this is not the exact area, but is what is required for
C**** some B-grid conservation quantities
      REAL*8, PUBLIC, DIMENSION(JM) :: DXYP,BYDXYP
      REAL*8, PUBLIC, DIMENSION(IM,JM) :: aDXYP
      REAL*8, PUBLIC, DIMENSION(:,:), ALLOCATABLE ::
     &     AXYP,BYAXYP,LAT2D,LON2D,LAT2D_DG,LON2D_DG,SINLAT2D,COSLAT2D
     &    ,ddx_ci,ddx_cj,ddy_ci,ddy_cj
!@var WTJ area weighting used in JLMAP, JKMAP (for hemispheric means)
      REAL*8, PUBLIC, DIMENSION(JM,2,2) :: WTJ

!@var DXYV,BYDXYV area of grid box around velocity point (recip.)(m^2)
      REAL*8, PUBLIC, DIMENSION(JM) :: DXYV,BYDXYV

!@var  DXP,DYP,BYDXP,BYDYP distance between points on primary grid
!@+     (+inverse)
      REAL*8, PUBLIC, DIMENSION(JM) :: DXP,DYP,BYDXP,BYDYP
!@var  DXV,DYV distance between velocity points (secondary grid)
      REAL*8, PUBLIC, DIMENSION(JM) :: DXV,DYV
!@var  DXYN,DXYS half box areas to the North,South of primary grid point
      REAL*8, PUBLIC, DIMENSION(JM) :: DXYS,DXYN
!@var  SINP sin of latitude at primary grid points
!@var  SINLATV,COSLATV sin(latv), cos(latv)
!@var  COSP, COSV cos of latitude at primary, secondary latitudes
      REAL*8, PUBLIC, DIMENSION(JM) :: SINP,SINLATV
      REAL*8, PUBLIC, DIMENSION(JM) :: COSP,COSV,COSLATV,DXLATV
!@var  RAPVS,RAPVN,RAVPS,RAVPN area scalings for primary and sec. grid
      REAL*8, PUBLIC, DIMENSION(JM) :: RAPVS,RAPVN,RAVPS,RAVPN

!@var SINIV,COSIV,SINIP,COSIP longitud. sin,cos for wind,pressure grid
      REAL*8, PUBLIC, DIMENSION(IM) :: SINIV,COSIV,SINIP,COSIP,SINU,COSU
!@var  RAVJ scaling for A grid U/V to B grid points (func. of lat. j)
!@var  RAPJ scaling for B grid -> A grid conversion (1/4,1/im at poles)
      REAL*8, PUBLIC, DIMENSION(IM,JM) :: RAPJ,RAVJ
!@var  IDJJ J index of adjacent U/V points for A grid (func. of lat. j)
      INTEGER, PUBLIC, DIMENSION(IM,JM) :: IDJJ
!@var  IDIJ I index of adjacent U/V points for A grid (func. of lat/lon)
      INTEGER, PUBLIC, DIMENSION(:,:,:), allocatable :: IDIJ
!@var  KMAXJ varying number of adjacent velocity points
      INTEGER, PUBLIC, DIMENSION(JM) :: KMAXJ
!@var  IMAXJ varying number of used longitudes
      INTEGER, PUBLIC, DIMENSION(JM) :: IMAXJ
!@var  FCOR latitudinally varying coriolis parameter
      REAL*8, PUBLIC, DIMENSION(JM) :: FCOR

!@var JG_U, JG_KE lat. grids on which U-wind and KE are defined
!@+   (1 for primary latitudes, 2 for secondary latitudes)
!@+   Information for diagnostics.
      integer, parameter :: jg_u=2, jg_ke=2

      real*8, public :: acor,acor2,polwt
!@var J_BUDG a mapping array that takes every grid point to the 
!@+   zonal mean budget array
      integer, public, allocatable, dimension(:,:) :: J_BUDG
!@var j_0b, j_1b are the min/max zonal budget latitudes for this processor
      integer, public :: j_0b, j_1b

      CONTAINS

      SUBROUTINE GEOM_B
!@sum  GEOM_B Calculate spherical geometry for B grid
!@auth Original development team (modifications by G. Schmidt)
!@ver  1.0 (B grid version)
      !use domain_decomp_atm, only : grid, hasSouthPole, hasNorthPole
      use domain_decomp_1d, only : grid, hasSouthPole, hasNorthPole
      IMPLICIT NONE
      REAL*8, PARAMETER :: EDPERD=1.,EDPERY = 365.

      INTEGER :: I,J,K,IM1  !@var I,J,K,IM1  loop variables
      INTEGER :: JVPO,JMHALF
      REAL*8  :: RAVPO,LAT1,COSP1,DXP1, SINV,SINVm1

      integer :: i_0h,i_1h,j_0h,j_1h,i_0,i_1,j_0,j_1,j_0s,j_1s

      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo
      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop
      j_0s = grid%j_strt_skp
      j_1s = grid%j_stop_skp

      allocate(idij(im,im,grid%j_strt_halo:grid%j_stop_halo))

      allocate(
     &       axyp(i_0h:i_1h,j_0h:j_1h)
     &    ,byaxyp(i_0h:i_1h,j_0h:j_1h)
     &    ,lat2d(i_0h:i_1h,j_0h:j_1h)
     &    ,lat2d_dg(i_0h:i_1h,j_0h:j_1h)
     &    ,lon2d(i_0h:i_1h,j_0h:j_1h)
     &    ,lon2d_dg(i_0h:i_1h,j_0h:j_1h)
     &    ,sinlat2d(i_0h:i_1h,j_0h:j_1h)
     &    ,coslat2d(i_0h:i_1h,j_0h:j_1h)
     &    ,ddx_ci(i_0h:i_1h,j_0h:j_1h)
     &    ,ddx_cj(i_0h:i_1h,j_0h:j_1h)
     &    ,ddy_ci(i_0h:i_1h,j_0h:j_1h)
     &    ,ddy_cj(i_0h:i_1h,j_0h:j_1h)
     &     )

C**** latitudinal spacing depends on whether you have even spacing or
C**** a partial box at the pole

      DLAT_DG=180./REAL(JM)                   ! even spacing (default)
      IF (JM.eq.46) DLAT_DG=180./REAL(JM-1)   ! 1/2 box at pole for 4x5
cc    IF (JM.eq.24) DLAT_DG=180./REAL(JM-1)   ! 1/2 box at pole, orig 8x10
      IF (JM.eq.24) DLAT_DG=180./REAL(JM-1.5) ! 1/4 box at pole, 'real' 8x10
      DLATM=60.*DLAT_DG
      DLAT=DLAT_DG*radian
      LAT(1)  = -.25*TWOPI
      LAT(JM) = -LAT(1)
      SINP(1)  = -1.
      SINP(JM) = 1.
      COSP(1)  = 0.
      COSP(JM) = 0.
      DXP(1)  = 0.
      DXP(JM) = 0.
      DO J=2,JM-1
        LAT(J)  = DLAT*(J-FJEQ)
        SINP(J) = SIN(LAT(J))
        COSP(J) = COS(LAT(J))
        DXP(J)  = RADIUS*DLON*COSP(J)
      END DO
      BYDXP(2:JM-1) = 1.D0/DXP(2:JM-1)
      LAT1    = DLAT*(1.-FJEQ)
      COSP1   = COS(LAT1)
      DXP1    = RADIUS*DLON*COSP1
      DO J=2,JM
        COSV(J) = .5*(COSP(J-1)+COSP(J))
        DXV(J)  = .5*(DXP(J-1)+DXP(J))
        DYV(J)  = RADIUS*(LAT(J)-LAT(J-1))
C**** The following corrections have no effect for half polar boxes
C**** but are important for full and quarter polar box cases.
        IF (J.eq.2) THEN
          polwt = cosv(j)
          COSV(J) = .5*(COSP1+COSP(J))
          DXV(J)  = .5*(DXP1+DXP(J))
        END IF
        IF (J.eq.JM) THEN
          COSV(J) = .5*(COSP(J-1)+COSP1)
          DXV(J)  = .5*(DXP(J-1)+DXP1)
        END IF
C****
      END DO
      DYP(1)  = RADIUS*(LAT(2)-LAT(1)-0.5*DLAT)
      DYP(JM) = RADIUS*(LAT(JM)-LAT(JM-1)-0.5*DLAT)

      SINV    = Sin (DLAT*(1+.5-FJEQ))
      DXYP(1) = RADIUS*RADIUS*DLON*(SINV+1)
      BYDXYP(1) = 1./DXYP(1)

      aDXYP(:,1) = DXYP(1) 

      SINVm1  = Sin (DLAT*(JM-.5-FJEQ))
      DXYP(JM)= RADIUS*RADIUS*DLON*(1-SINVm1)
      BYDXYP(JM) = 1./DXYP(JM)

      aDXYP(:,JM) = DXYP(JM) 

      DXYS(1)  = 0.
      DXYS(JM) = DXYP(JM)
      DXYN(1)  = DXYP(1)
      DXYN(JM) = 0.
      polwt = (cosv(3)-cosv(2))/(cosv(3)-polwt)
      DO J=2,JM-1
        DYP(J)  =  radius*dlat !.5*(DYV(J)+DYV(J+1))
        SINVm1  = Sin (DLAT*(J-.5-FJEQ))
        SINV    = Sin (DLAT*(J+.5-FJEQ))
        DXYP(J) = RADIUS*RADIUS*DLON*(SINV-SINVm1)

        aDXYP(:,J) = DXYP(J) 

        BYDXYP(J) = 1./DXYP(J)
        DXYS(J) = .5*DXYP(J)
        DXYN(J) = .5*DXYP(J)
      END DO
      do j=2,jm
        latv(j) = dlat*(j-.5-fjeq)
        coslatv(j) = cos(latv(j))
        sinlatv(j) = sin(latv(j))
        dxlatv(j) = radius*dlon*coslatv(j)
      enddo
      BYDYP(:) = 1.D0/DYP(:)
      RAVPS(1)  = 0.
      RAPVS(1)  = 0.
      RAVPN(JM) = 0.
      RAPVN(JM) = 0.
      DO J=2,JM
        DXYV(J) = DXYN(J-1)+DXYS(J)
        BYDXYV(J) = 1./DXYV(J)
        RAPVS(J)   = .5*DXYS(J)/DXYV(J)
        RAPVN(J-1) = .5*DXYN(J-1)/DXYV(J)
        RAVPS(J)   = .5*DXYS(J)/DXYP(J)
        RAVPN(J-1) = .5*DXYN(J-1)/DXYP(J-1)
      END DO
      acor = dxyv(2)/(.5*dxp(2)*dyv(2)) ! gridbox area correction factor
      acor2 = dxyv(2)/(dxv(2)*dyv(2))
C**** LONGITUDES (degrees); used in ILMAP
      LON_DG(1,1) = -180.+360./(2.*FLOAT(IM))
      LON_DG(1,2) = -180.+360./    FLOAT(IM)
      DO I=2,IM
        LON_DG(I,1) = LON_DG(I-1,1)+360./FLOAT(IM)
        LON_DG(I,2) = LON_DG(I-1,2)+360./FLOAT(IM)
      END DO
C**** LATITUDES (degrees); used extensively in the diagn. print routines
      LAT_DG(1,1:2)=-90.
      LAT_DG(JM,1)=90.
      DO J=2,JM-1
        LAT_DG(J,1)=DLAT_DG*(J-FJEQ)    ! primary (tracer) latitudes
      END DO
      DO J=2,JM
        LAT_DG(J,2)=DLAT_DG*(J-JM/2-1)  ! secondary (velocity) latitudes
      END DO
C**** WTJ: area weighting for JKMAP, JLMAP hemispheres
      JMHALF= JM/2
      DO J=1,JM
        WTJ(J,1,1)=1.
        WTJ(J,2,1)=2.*FIM*DXYP(J)/AREAG
      END DO
      DO J=2,JM
        WTJ(J,1,2)=1.
        WTJ(J,2,2)=2.*FIM*DXYV(J)/AREAG
      END DO
cgsfc      WTJ(JMHALF+1,1,2)=.5
cgsfc      WTJ(JMHALF+1,2,2)=WTJ(JMHALF+1,2,2)/2.
      WTJ(1,1,2)=0.
      WTJ(1,2,2)=0.
C**** CALCULATE CORIOLIS PARAMETER
c      OMEGA = TWOPI*(EDPERD+EDPERY)/(EDPERD*EDPERY*SDAY)
      FCOR(1)  = -OMEGA*DXV(2)*DXV(2)/DLON
      FCOR(JM) = OMEGA*DXV(JM)*DXV(JM)/DLON
      DO J=2,JM-1
        FCOR(J) = OMEGA*(DXV(J)*DXV(J)-DXV(J+1)*DXV(J+1))/DLON
      END DO

C**** Set indexes and scalings for the influence of A grid points on
C**** adjacent velocity points

C**** Calculate relative directions of polar box to nearby U,V points
      DO I=1,IM
        SINIV(I)=SIN((I-1)*DLON)
        COSIV(I)=COS((I-1)*TWOPI*BYIM) ! DLON)
        LON(I)=DLON*(I-.5)
        LONV(I)=DLON*I
        SINIP(I)=SIN(LON(I))
        COSIP(I)=COS(LON(I))
        SINU(I) =SIN (I*TWOPI/IM)
        COSU(I) =COS (I*TWOPI/IM)
      END DO
      SINU(IM) = 0
      COSU(IM) = 1

C**** Conditions at the poles
      DO J=1,JM,JM-1
        IF(J.EQ.1) THEN
          JVPO=2
          RAVPO=2.*RAPVN(1)
        ELSE
          JVPO=JM
          RAVPO=2.*RAPVS(JM)
        END IF
        KMAXJ(J)=IM
        IMAXJ(J)=1
        RAVJ(1:KMAXJ(J),J)=RAVPO
        RAPJ(1:KMAXJ(J),J)=BYIM
        IDJJ(1:KMAXJ(J),J)=JVPO
        if(j.ge.grid%j_strt_halo .and. j.le.grid%j_stop_halo) then
          DO K=1,KMAXJ(J)
            IDIJ(K,1:IM,J)=K
          END DO
        endif
      END DO
C**** Conditions at non-polar points
      DO J=2,JM-1
        KMAXJ(J)=4
        IMAXJ(J)=IM
        DO K=1,2
          RAVJ(K,J)=RAPVS(J)
          RAPJ(K,J)=RAVPS(J)    ! = .25
          IDJJ(K,J)=J
          RAVJ(K+2,J)=RAPVN(J)
          RAPJ(K+2,J)=RAVPN(J)  ! = .25
#ifdef SCM
          IDJJ(K+2,J)=J
#else
          IDJJ(K+2,J)=J+1
#endif
        END DO
        if(j.ge.grid%j_strt_halo .and. j.le.grid%j_stop_halo) then
          IM1=IM
          DO I=1,IM
#ifdef SCM
            IDIJ(1,I,J)=I
            IDIJ(2,I,J)=I
            IDIJ(3,I,J)=I
            IDIJ(4,I,J)=I
#else 
            IDIJ(1,I,J)=IM1
            IDIJ(2,I,J)=I
            IDIJ(3,I,J)=IM1
            IDIJ(4,I,J)=I
            IM1=I
#endif
          END DO
        endif
      END DO

c
c temporary
c
      do j=max(1,j_0h),min(jm,j_1h)
      do i=i_0h,i_1h
        axyp(i,j) = dxyp(j)
        byaxyp(i,j) = bydxyp(j)
        lat2d(i,j) = lat(j)
        lat2d_dg(i,j) = lat_dg(j,1)
        sinlat2d(i,j) = sin(lat(j))
        coslat2d(i,j) = cos(lat(j))
      enddo
      do i=i_0,i_1
        lon2d(i,j) = lon(i) ! i=1 is dlon/2
        lon2d_dg(i,j) = lon_dg(i,1) ! i=1 is -180 + dlon/2
      enddo
      enddo

      do j=j_0s,j_1s
      do i=i_0,i_1
        ddx_ci(i,j) =  .5/(radius*dlon*cosp(j))
        ddx_cj(i,j) = 0.
        ddy_ci(i,j) = 0.
        ddy_cj(i,j) =  .5/(radius*dlat)
c        dloni = lon2d(i+1,j)-lon2d(i-1,j)
c        dlati = lat2d(i+1,j)-lat2d(i-1,j)
c        dlonj = lon2d(i,j+1)-lon2d(i,j-1)
c        dlatj = lat2d(i,j+1)-lat2d(i,j-1)
c        if(dloni.lt.0.) dloni=dloni+twopi
c        if(dlonj.lt.0.) dlonj=dlonj+twopi
c        bydet = 1d0/(dloni*dlatj-dlonj*dlati)
c        ddx_ci(i,j) =  dlatj*bydet/(radius*coslat2d(i,j))
c        ddx_cj(i,j) = -dlonj*bydet/(radius*coslat2d(i,j))
c        ddy_ci(i,j) = -dlati*bydet/radius
c        ddy_cj(i,j) =  dloni*bydet/radius
      enddo
      enddo

C**** set up mapping arrays for budget/conserv diags
!      call set_j_budg   !called after lat2d_dg is initialized
!      call set_wtbudg() !sets area weights
 
      RETURN
      END SUBROUTINE GEOM_B

      subroutine lonlat_to_ij(ll,ij)
c converts lon,lat=ll(1:2) into model i,j=ij(1:2) (primary grid)
c this version is for the latlon grid.  ll are in degrees east.
      implicit none
      real*8, intent(in) :: ll(2)
      integer, intent(out) :: ij(2)
      real*8 :: dlon_dg
      dlon_dg = 360./dble(im)

      ij(1) = nint( .5*(im + 1) + ll(1)/dlon_dg )
      ij(2) = nint( .5*(jm + 1) + ll(2)/dlat_dg )
      if(ij(2)>jm) ij(2)=jm
      if(ij(2)<1 ) ij(2)=1
      return
      end subroutine lonlat_to_ij

      subroutine lonlat_to_tile(ll,tile)
c returns tile = 1, latlon grid has only one face 
      real*8, intent(in) :: ll(2)
      integer, intent(out) :: tile
      tile=1
      return
      end subroutine lonlat_to_tile


      function lon_to_I(deg)
c converts longitude into model i (primary grid)
c this version is for the latlon grid.  DEG is in degrees east.
      implicit none
      integer lon_to_I
      real*8, intent(in) :: DEG
      real*8 :: dlon_dg
      dlon_dg = 360./dble(im)
      lon_to_I = nint( .5*(im + 1) + DEG/dlon_dg )
      return
      end function lon_to_I

      function lat_to_J(deg)
c converts latitude into model j (primary grid)
c this version is for the latlon grid.
      implicit none
      integer lat_to_J
      real*8, intent(in) :: DEG
      lat_to_J = nint( .5*(jm + 1) + DEG/dlat_dg )
      if( lat_to_J>jm ) lat_to_J = jm
      if( lat_to_J<1 )  lat_to_J = 1
      return
      end function lat_to_J

      END MODULE GEOM

      module RAD_COSZ0
      implicit none
      save
      private
      public :: cosz_init,coszs,coszt,daily_cosz

!@var SIND,COSD sin,cos of solar declination angle
!@+   (these are local copies that are set by daily_cosz)
      real*8 :: sind,cosd
!@var SINJ,COSJ sines and cosines for zenith angle calculation
      REAL*8, ALLOCATABLE, DIMENSION(:) :: SINJ,COSJ

      contains

      subroutine daily_cosz(SIND_in,COSD_in, COSZ_day,DUSK)
c Resets parameters needed for zenith angle calculations, and calculates
c daily-average cosz and the hour of dusk
      !USE DOMAIN_DECOMP_ATM, ONLY: grid
      USE DOMAIN_DECOMP_1D, ONLY: grid
      USE MODEL_COM
      USE CONSTANT, only : PI
      implicit none
      REAL*8 SIND_in,COSD_in
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     COSZ_day, ! average of cos(zenith angle) for the current day
     &     DUSK      ! hour of dusk, radians from local noon
      INTEGER :: J
      REAL*8 :: SJSD,CJCD,CJCD_SDUSK
      SIND = SIND_in
      COSD = COSD_in
      do j=grid%j_strt,grid%j_stop
        SJSD=SINJ(J)*SIND
        CJCD=COSJ(J)*COSD
        IF(CJCD-SJSD.LT.0.) THEN     ! constant daylight
          DUSK(:,J) = PI
          COSZ_day(:,J) = SJSD
        ELSEIF(CJCD+SJSD.LT.0.) THEN ! constant darkness
          DUSK(:,J) = 0D0
          COSZ_day(:,J) = 0D0
        ELSE
          DUSK(:,J) = ACOS(-SJSD/CJCD)
          CJCD_SDUSK = SQRT((CJCD-SJSD)*(CJCD+SJSD))
          COSZ_day(:,J) = max(0d0,(SJSD*DUSK(1,J)+CJCD_SDUSK)/PI)
        ENDIF
      enddo
      return
      end subroutine daily_cosz

      subroutine cosz_init
      !use domain_decomp_atm, only : grid, hasSouthPole, hasNorthPole
      use domain_decomp_1D, only : grid, hasSouthPole, hasNorthPole
      use constant, only : twopi
      use model_com, only : jm
      use geom, only : dlat
      real*8 :: phis,cphis,sphis,phin,cphin,sphin,phim
      integer :: j_0,j_1,j_0s,j_1s,j_0h,j_1h,j

      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo
      j_0 = grid%j_strt
      j_1 = grid%j_stop
      j_0s = grid%j_strt_skp
      j_1s = grid%j_stop_skp
    
      allocate(sinj(j_0h:j_1h),cosj(j_0h:j_1h))

C**** COMPUTE THE AREA WEIGHTED LATITUDES AND THEIR SINES AND COSINES
      if (hasSouthPole(GRID)) then
        PHIS=-.25*TWOPI
        SPHIS=-1.
        CPHIS=0.
      else
        PHIS=DLAT*(J_0-1-.5*JM)
        SPHIS=SIN(PHIS)
        CPHIS=COS(PHIS)
      end if
      DO J=J_0,J_1S
        PHIN=DLAT*(J-.5*JM)
        SPHIN=SIN(PHIN)
        CPHIN=COS(PHIN)
        PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)
        SINJ(J)=SIN(PHIM)
        COSJ(J)=COS(PHIM)
        PHIS=PHIN
        SPHIS=SPHIN
        CPHIS=CPHIN
      END DO
      IF (hasNorthPole(GRID)) THEN
        PHIN=.25*TWOPI
        SPHIN=1.
        CPHIN=0.
        PHIM=( PHIN*SPHIN + CPHIN
     *        -PHIS*SPHIS - CPHIS)
     *        /(SPHIN - SPHIS)
        SINJ(JM)=SIN(PHIM)
        COSJ(JM)=COS(PHIM)
      END IF
      return
      end subroutine cosz_init

      subroutine COSZS (ROT1,ROT2,COSZ,COSZA)
      REAL*8 ROT1,ROT2
      REAL*8, DIMENSION(:,:) :: COSZ,COSZA
      call COSZT (ROT1,ROT2,COSZ,COSZA)
      end subroutine COSZS

      subroutine COSZT (ROT1,ROT2,COSZ,COSZA)
!@sum  COSZ0 calculates Earth's zenith angle, weighted by time/sunlight
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : twopi
      USE MODEL_COM
      USE GEOM, only : lon,sinip,cosip
c      USE RAD_COM, only : cosd,sind
      !USE DOMAIN_DECOMP_ATM, ONLY: grid, hasSouthPole, hasNorthPole
      USE DOMAIN_DECOMP_1D, ONLY: grid, hasSouthPole, hasNorthPole
      IMPLICIT NONE
      SAVE
      REAL*8 ROT1,ROT2
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     COSZ
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
     &     , optional :: COSZA
      REAL*8, DIMENSION(IM) :: LT1,LT2,SLT1,SLT2,S2LT1,S2LT2
C**** ZERO1 HAS TO EQUAL THE CUT-OFF VALUE FOR COSZ USED IN SOLAR
C**** COSZS WORKS CORRECTLY ONLY IF ZERO1 >> 1.D-3
      REAL*8, PARAMETER :: ZERO1=1.D-2
      INTEGER I,J
      REAL*8 S2DAWN,S2DUSK,ECOSZ,ECOSQZ,CLT1,CLT2,ZERO2,CDUSK,DUSK,DAWN
     *     ,SDUSK,SDAWN,CJCD,SJSD,SR1,CR1,SR2,CR2,DROT
      INTEGER :: J_0, J_1, J_0S, J_1S

C****
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      J_0S = grid%J_STRT_SKP
      J_1S = grid%J_STOP_SKP

      if ( present(COSZA) ) goto 777 ! COSZS
C****
!      ENTRY COSZT (ROT1,ROT2,COSZ)
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE WEIGHTED BY DAYTIME
C**** HOURS FROM ROT1 TO ROT2, GREENWICH MEAN TIME IN RADIANS.  ROT1
C**** MUST BE BETWEEN 0 AND 2*PI.  ROT2 MUST BE BETWEEN ROT1 AND
C**** ROT1+2*PI.  I=1 MUST LIE ON THE INTERNATIONAL DATE LINE.
C****
      DROT=ROT2-ROT1
C**** COMPUTE THE SINES AND COSINES OF THE INITIAL AND FINAL GMT'S
      SR1=SIN(ROT1)
      CR1=COS(ROT1)
      SR2=SIN(ROT2)
      CR2=COS(ROT2)
C**** COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO
C****   NOON) AND THEIR SINES AND COSINES
      DO I=1,IM
        LT1(I)=ROT1+LON(I)
        SLT1(I)=SR1*COSIP(I)+CR1*SINIP(I)
        LT2(I)=ROT2+LON(I)
        SLT2(I)=SR2*COSIP(I)+CR2*SINIP(I)
      END DO
C****
C**** CALCULATION FOR POLAR GRID BOXES
C****
      DO J=1,JM,JM-1
        IF(((J .EQ. 1) .AND. (hasSouthPole(grid))) .OR.
     *     ((J .EQ. JM) .AND. (hasNorthPole(grid)))) THEN
           SJSD=SINJ(J)*SIND
           CJCD=COSJ(J)*COSD
           IF (SJSD+CJCD.GT.ZERO1) THEN
              IF (SJSD-CJCD.LT.0.) THEN
C**** AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES
                 DUSK=ACOS(-SJSD/CJCD)
                 SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
                 DAWN=-DUSK
                 SDAWN=-SDUSK
                 COSZ(1,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/TWOPI
              ELSE
C**** CONSTANT DAYLIGHT NEAR THE POLES
                 COSZ(1,J)=SJSD
              END IF
           ELSE
C**** CONSTANT NIGHTIME NEAR THE POLES
              COSZ(1,J)=0.
           END IF
        END IF
      END DO
C****
C**** LOOP OVER NON-POLAR LATITUDES
C****
      DO 500 J=J_0S,J_1S
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.GT.ZERO1) THEN
      IF (SJSD-CJCD.LT.0.) THEN
C**** COMPUTE DAWN AND DUSK (AT LOCAL TIME) AND THEIR SINES
      DUSK=ACOS(-SJSD/CJCD)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      DAWN=-DUSK
      SDAWN=-SDUSK
C**** NEITHER CONSTANT DAYTIME NOR CONSTANT NIGHTIME AT THIS LATITUDE,
C**** LOOP OVER LONGITUDES
      ZERO2=ZERO1/CJCD
      DO 400 I=1,IM
C**** FORCE DUSK TO LIE BETWEEN LT1 AND LT1+2*PI
      IF (DUSK.GT.LT1(I)+ZERO2) GO TO 220
      DUSK=DUSK+TWOPI
      DAWN=DAWN+TWOPI
 220  IF (DAWN.LT.LT2(I)-ZERO2) GO TO 240
C**** CONTINUOUS NIGHTIME FROM INITIAL TO FINAL TIME
      COSZ(I,J)=0.
      GO TO 400
 240  IF (DAWN.GE.LT1(I)) GO TO 300
      IF (DUSK.LT.LT2(I)) GO TO 260
C**** CONTINUOUS DAYLIGHT FROM INITIAL TIME TO FINAL TIME
      COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/DROT
      GO TO 400
 260  IF (DAWN+TWOPI.LT.LT2(I)-ZERO2) GO TO 280
C**** DAYLIGHT AT INITIAL TIME AND NIGHT AT FINAL TIME
      COSZ(I,J)=(SJSD*(DUSK-LT1(I))+CJCD*(SDUSK-SLT1(I)))/DROT
      GO TO 400
C**** DAYLIGHT AT INITIAL AND FINAL TIMES WITH NIGHTIME IN BETWEEN
 280  COSZ(I,J)=(SJSD*(LT2(I)-DAWN-TWOPI+DUSK-LT1(I))+CJCD*
     *     (SLT2(I)-SDAWN+SDUSK-SLT1(I)))/DROT
      GO TO 400
 300  IF (DUSK.LT.LT2(I)) GO TO 320
C**** NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME
      COSZ(I,J)=(SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN))/DROT
      GO TO 400
C**** NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN
 320  COSZ(I,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/DROT
 400  CONTINUE
      ELSE
C**** CONSTANT DAYLIGHT AT THIS LATITUDE
         DO I=1,IM
          COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/DROT
        END DO
      END IF
      ELSE
C**** CONSTANT NIGHTIME AT THIS LATITUDE
        COSZ(1:IM,J)=0.
      END IF
 500  CONTINUE
      RETURN
C****
C****
!      ENTRY COSZS (ROT1,ROT2,COSZ,COSZA)
 777  continue
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE TWICE, FIRST WEIGHTED BY THE
C**** DAYTIME HOURS FROM ROT1 TO ROT2 AND SECONDLY WEIGHTED BY THE
C**** INCIDENT SUN LIGHT FROM ROT1 TO ROT2.  COSZT MUST HAVE BEEN
C**** CALLED JUST PREVIOUSLY.
C****
      DROT=ROT2-ROT1
C**** COMPUTE THE SINES AND COSINES OF THE INITIAL AND FINAL GMT'S
      SR1=SIN(ROT1)
      CR1=COS(ROT1)
      SR2=SIN(ROT2)
      CR2=COS(ROT2)
C**** COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO
C****   NOON) AND THEIR SINES AND COSINES
      DO I=1,IM
        LT1(I)=ROT1+LON(I)
        SLT1(I)=SR1*COSIP(I)+CR1*SINIP(I)
        CLT1=CR1*COSIP(I)-SR1*SINIP(I)
        S2LT1(I)=2.*SLT1(I)*CLT1
        LT2(I)=ROT2+LON(I)
        SLT2(I)=SR2*COSIP(I)+CR2*SINIP(I)
        CLT2=CR2*COSIP(I)-SR2*SINIP(I)
        S2LT2(I)=2.*SLT2(I)*CLT2
      END DO
C****
C**** CALCULATION FOR POLAR GRID BOXES
C****
      DO J=1,JM,JM-1
        IF(((J .EQ. 1) .AND. (hasSouthPole(grid))) .OR.
     *     ((J .EQ. JM) .AND. (hasNorthPole(grid)))) THEN
           SJSD=SINJ(J)*SIND
           CJCD=COSJ(J)*COSD
           IF (SJSD+CJCD.GT.ZERO1) THEN
              IF (SJSD-CJCD.LT.0.) THEN
C**** AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES
                 CDUSK=-SJSD/CJCD
                 DUSK=ACOS(CDUSK)
                 SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
                 S2DUSK=2.*SDUSK*CDUSK
                 DAWN=-DUSK
                 SDAWN=-SDUSK
                 S2DAWN=-S2DUSK
                 ECOSZ=SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN)
                 ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SDAWN)+
     *                .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))
                 COSZ(1,J)=ECOSZ/TWOPI
                 COSZA(1,J)=ECOSQZ/ECOSZ
              ELSE
C**** CONSTANT DAYLIGHT NEAR THE POLES
                 ECOSZ=SJSD*TWOPI
                 ECOSQZ=SJSD*ECOSZ+.5*CJCD*CJCD*TWOPI
                 COSZ(1,J)=ECOSZ/TWOPI
                 COSZA(1,J)=ECOSQZ/ECOSZ
              END IF
           ELSE
C**** CONSTANT NIGHTIME NEAR THE POLES
              COSZ(1,J)=0.
              COSZA(1,J)=0.
           END IF
        END IF
      END DO
C****
C**** LOOP OVER NON-POLAR LATITUDES
C****
      DO 900 J=J_0S,J_1S
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.GT.ZERO1) THEN
      IF (SJSD-CJCD.LT.0.) THEN
C**** COMPUTE DAWN AND DUSK (AT LOCAL TIME) AND THEIR SINES
      CDUSK=-SJSD/CJCD
      DUSK=ACOS(CDUSK)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      S2DUSK=2.*SDUSK*CDUSK
      DAWN=-DUSK
      SDAWN=-SDUSK
      S2DAWN=-S2DUSK
C**** NEITHER CONSTANT DAYTIME NOR CONSTANT NIGHTIME AT THIS LATITUDE,
C**** LOOP OVER LONGITUDES
      ZERO2=ZERO1/CJCD
      DO 800 I=1,IM
C**** FORCE DUSK TO LIE BETWEEN LT1 AND LT1+2*PI
      IF (DUSK.GT.LT1(I)+ZERO2) GO TO 620
      DUSK=DUSK+TWOPI
      DAWN=DAWN+TWOPI
  620 IF (DAWN.LT.LT2(I)-ZERO2) GO TO 640
C**** CONTINUOUS NIGHTIME FROM INITIAL TO FINAL TIME
      COSZ(I,J)=0.
      COSZA(I,J)=0.
      GO TO 800
  640 IF (DAWN.GE.LT1(I)) GO TO 700
      IF (DUSK.LT.LT2(I)) GO TO 660
C**** CONTINUOUS DAYLIGHT FROM INITIAL TIME TO FINAL TIME
      ECOSZ=SJSD*DROT+CJCD*(SLT2(I)-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SLT1(I))+
     *  .5*CJCD*(DROT+.5*(S2LT2(I)-S2LT1(I))))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
  660 IF (DAWN+TWOPI.LT.LT2(I)-ZERO2) GO TO 680
C**** DAYLIGHT AT INITIAL TIME AND NIGHT AT FINAL TIME
      ECOSZ=SJSD*(DUSK-LT1(I))+CJCD*(SDUSK-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SLT1(I))+
     *  .5*CJCD*(DUSK-LT1(I)+.5*(S2DUSK-S2LT1(I))))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
C**** DAYLIGHT AT INITIAL AND FINAL TIMES WITH NIGHTIME IN BETWEEN
  680 ECOSZ=SJSD*(DROT-DAWN-TWOPI+DUSK)+
     *  CJCD*(SLT2(I)-SDAWN+SDUSK-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SLT1(I)+SLT2(I)-SDAWN)+
     *  .5*CJCD*(DUSK+DROT-DAWN-TWOPI+
     *  .5*(S2DUSK-S2LT1(I)+S2LT2(I)-S2DAWN)))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
  700 IF (DUSK.GE.LT2(I)) THEN
C**** NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME
        ECOSZ=SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN)
        ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SDAWN)+
     *       .5*CJCD*(LT2(I)-DAWN+.5*(S2LT2(I)-S2DAWN)))
        COSZ(I,J)=ECOSZ/DROT
        COSZA(I,J)=ECOSQZ/ECOSZ
      ELSE
C**** NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN
        ECOSZ=SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN)
        ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SDAWN)+
     *       .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))
        COSZ(I,J)=ECOSZ/DROT
        COSZA(I,J)=ECOSQZ/ECOSZ
      END IF
  800 CONTINUE
      ELSE
C**** CONSTANT DAYLIGHT AT THIS LATITUDE
        DO I=1,IM
          ECOSZ=SJSD*DROT+CJCD*(SLT2(I)-SLT1(I))
          ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SLT1(I))+
     *         .5*CJCD*(DROT+.5*(S2LT2(I)-S2LT1(I))))
          COSZ(I,J)=ECOSZ/DROT
          COSZA(I,J)=ECOSQZ/ECOSZ
        END DO
      END IF
C**** CONSTANT NIGHTIME AT THIS LATITUDE
      ELSE
        COSZ(1:IM,J)=0.
        COSZA(1:IM,J)=0.
      END IF
  900 CONTINUE
      RETURN
      END subroutine COSZT

      end module RAD_COSZ0

!==============================================================================
! LTM: ATM_UTILS.f
!==============================================================================
#ifdef CUBED_SPHERE
      SUBROUTINE PGRAD_PBL
!@sum  PGRAD_PBL calculates surface/layer 1 pressure gradients for pbl
!@sum  This version works for a nonorthogonal grid
!@auth M. Kelley
!@ver  1.0
C**** As this is written, it must be called after the call to CALC_AMPK
C**** after DYNAM (since it uses pk/pmid). It would be better if it used
C**** SPA and PU directly from the dynamics. (Future work).
      USE CONSTANT, only : rgas
      USE MODEL_COM, only : t,p,zatmo,sig
      USE GEOM, only : ddx_ci,ddx_cj,ddy_ci,ddy_cj
      USE DYNAMICS, only : phi,dpdy_by_rho,dpdy_by_rho_0,dpdx_by_rho
     *     ,dpdx_by_rho_0,pmid,pk
      USE DOMAIN_DECOMP_1D, only : grid, GET, HALO_UPDATE
      IMPLICIT NONE
      REAL*8 :: by_rho1
      real*8 :: dpsi,dpsj,dg1i,dg1j,dgsi,dgsj
      real*8 :: dpsdx,dpsdy,dg1dx,dg1dy,dgsdx,dgsdy
      INTEGER I,J

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, I_0, I_1
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP


C**** (Pressure gradient)/density at first layer and surface
C**** to be used in the PBL, at the primary grids

      CALL HALO_UPDATE(grid, P)
      CALL HALO_UPDATE(grid, PHI(:,:,1))
c      CALL HALO_UPDATE(grid, ZATMO)

      DO J=J_0,J_1
      DO I=I_0,I_1
        dpsi = p(i+1,j)-p(i-1,j)
        dg1i = phi(i+1,j,1)-phi(i-1,j,1)
        dgsi = zatmo(i+1,j)-zatmo(i-1,j)
        dpsj = p(i,j+1)-p(i,j-1)
        dg1j = phi(i,j+1,1)-phi(i,j-1,1)
        dgsj = zatmo(i,j+1)-zatmo(i,j-1)
        dpsdx = dpsi*ddx_ci(i,j) + dpsj*ddx_cj(i,j)
        dpsdy = dpsi*ddy_ci(i,j) + dpsj*ddy_cj(i,j)
        dg1dx = dg1i*ddx_ci(i,j) + dg1j*ddx_cj(i,j)
        dg1dy = dg1i*ddy_ci(i,j) + dg1j*ddy_cj(i,j)
        dgsdx = dgsi*ddx_ci(i,j) + dgsj*ddx_cj(i,j)
        dgsdy = dgsi*ddy_ci(i,j) + dgsj*ddy_cj(i,j)
        by_rho1=(rgas*t(I,J,1)*pk(1,I,J))/(pmid(1,I,J))
        DPDX_BY_RHO(I,J)=dpsdx*sig(1)*by_rho1+dg1dx
        DPDX_BY_RHO_0(I,J)=dpsdx*by_rho1+dgsdx
        DPDY_BY_RHO(I,J)=dpsdy*sig(1)*by_rho1+dg1dy
        DPDY_BY_RHO_0(I,J)=dpsdy*by_rho1+dgsdy
      ENDDO
      ENDDO

      return
      END SUBROUTINE PGRAD_PBL

#else
      SUBROUTINE PGRAD_PBL
!@sum  PGRAD_PBL calculates surface/layer 1 pressure gradients for pbl
!@auth Ye Cheng
!@ver  1.0
C**** As this is written, it must be called after the call to CALC_AMPK
C**** after DYNAM (since it uses pk/pmid). It would be better if it used
C**** SPA and PU directly from the dynamics. (Future work).
      USE CONSTANT, only : rgas
      USE MODEL_COM, only : im,jm,t,p,zatmo,sig,byim
#ifdef CUBED_SPHERE
      USE GEOM, only : cosip,sinip
#else
      USE GEOM, only : bydyp,bydxp,cosip,sinip
#endif
      USE DYNAMICS, only : phi,dpdy_by_rho,dpdy_by_rho_0,dpdx_by_rho
     *     ,dpdx_by_rho_0,pmid,pk
      USE DOMAIN_DECOMP_1D, only : grid, GET
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : NORTH, SOUTH
      USE DOMAIN_DECOMP_1D, only : haveLatitude
      IMPLICIT NONE
      REAL*8 by_rho1,dpx1,dpy1,dpx0,dpy0,hemi
      INTEGER I,J,K,IP1,IM1,J1
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)


C**** (Pressure gradient)/density at first layer and surface
C**** to be used in the PBL, at the primary grids

      ! for dPdy/rho at non-pole grids
      CALL HALO_UPDATE(grid, P,   FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(grid, PHI, FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(grid, ZATMO, FROM=SOUTH+NORTH)

      DO I=1,IM
        DO J=J_0S,J_1S
          by_rho1=(rgas*t(I,J,1)*pk(1,I,J))/(100.*pmid(1,I,J))
#ifndef SCM
#ifndef CUBED_SPHERE  /* bydyp on cubed sphere? */
          DPDY_BY_RHO(I,J)=(100.*(P(I,J+1)-P(I,J-1))*SIG(1)*by_rho1
     2         +PHI(I,J+1,1)-PHI(I,J-1,1))*BYDYP(J)*.5d0
          DPDY_BY_RHO_0(I,J)=(100.*(P(I,J+1)-P(I,J-1))*by_rho1
     2         +ZATMO(I,J+1)-ZATMO(I,J-1))*BYDYP(J)*.5d0
#endif
#endif
        END DO
      END DO

      ! for dPdx/rho at non-pole grids

      DO J=J_0S,J_1S
        IM1=IM-1
        I=IM
        DO IP1=1,IM
          by_rho1=(rgas*t(I,J,1)*pk(1,I,J))/(100.*pmid(1,I,J))
#ifndef CUBED_SPHERE  /* bydxp on cubed sphere? */
          DPDX_BY_RHO(I,J)=(100.*(P(IP1,J)-P(IM1,J))*SIG(1)*by_rho1
     2         +PHI(IP1,J,1)-PHI(IM1,J,1))*BYDXP(J)*.5d0
          DPDX_BY_RHO_0(I,J)=(100.*(P(IP1,J)-P(IM1,J))*by_rho1
     2         +ZATMO(IP1,J)-ZATMO(IM1,J))*BYDXP(J)*.5d0
#endif
          IM1=I
          I=IP1
        END DO
      END DO

      ! at poles

      IF (haveLatitude(grid, J=1)) THEN
        hemi = -1.; J1 = 2
        dpx1=0. ; dpy1=0.
        dpx0=0. ; dpy0=0.
        DO K=1,IM
          dpx1=dpx1+(DPDX_BY_RHO(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO(K,J1)*SINIP(K))
          dpy1=dpy1+(DPDY_BY_RHO(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO(K,J1)*SINIP(K))
          dpx0=dpx0+(DPDX_BY_RHO_0(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO_0(K,J1)*SINIP(K))
          dpy0=dpy0+(DPDY_BY_RHO_0(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO_0(K,J1)*SINIP(K))
        END DO
        DPDX_BY_RHO(1,1)  =dpx1*BYIM
        DPDY_BY_RHO(1,1)  =dpy1*BYIM
        DPDX_BY_RHO_0(1,1)=dpx0*BYIM
        DPDY_BY_RHO_0(1,1)=dpy0*BYIM
      END IF

      If (haveLatitude(grid, J=JM)) THEN
          hemi= 1.; J1=JM-1
        dpx1=0. ; dpy1=0.
        dpx0=0. ; dpy0=0.
        DO K=1,IM
          dpx1=dpx1+(DPDX_BY_RHO(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO(K,J1)*SINIP(K))
          dpy1=dpy1+(DPDY_BY_RHO(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO(K,J1)*SINIP(K))
          dpx0=dpx0+(DPDX_BY_RHO_0(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO_0(K,J1)*SINIP(K))
          dpy0=dpy0+(DPDY_BY_RHO_0(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO_0(K,J1)*SINIP(K))
        END DO
        DPDX_BY_RHO(1,JM)  =dpx1*BYIM
        DPDY_BY_RHO(1,JM)  =dpy1*BYIM
        DPDX_BY_RHO_0(1,JM)=dpx0*BYIM
        DPDY_BY_RHO_0(1,JM)=dpy0*BYIM
      END IF

      END SUBROUTINE PGRAD_PBL
#endif

      SUBROUTINE CALC_PIJL(lmax,p,pijl)
!@sum  CALC_PIJL Fills in P as 3-D
!@auth Jean Lerner
!@ver  1.0
      USE MODEL_COM, only : lm,ls1,psfmpt
C****
      USE DOMAIN_DECOMP_1D, Only : grid, GET
      implicit none
      REAL*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: p
      REAL*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: pijl
      integer :: l,lmax

      do l=1,ls1-1
        pijl(:,:,l) = p(:,:)
      enddo
      do l=ls1,lmax
        pijl(:,:,l) = PSFMPT
      enddo
      return
      end subroutine calc_pijl

      SUBROUTINE CALC_AMPK(LMAX)
!@sum  CALC_AMPK calculate air mass and pressure arrays
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : bygrav,kapa
      USE MODEL_COM, only : im,jm,lm,ls1,p
      USE DYNAMICS, only : plij,pdsig,pmid,pk,pedn,pek,sqrtp,am,byam
      USE DOMAIN_DECOMP_1D, Only : grid, GET, HALO_UPDATE
      IMPLICIT NONE

      INTEGER :: I,J,L  !@var I,J,L  loop variables
      INTEGER, INTENT(IN) :: LMAX !@var LMAX max. level for update
      REAL*8, DIMENSION(LMAX) :: PL,AML,PDSIGL,PMIDL
      REAL*8, DIMENSION(LMAX+1) :: PEDNL
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H, I_0H, I_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &               J_STRT_HALO= J_0H, J_STOP_HALO= J_1H,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

C**** Calculate air mass, layer pressures, P**K, and sqrt(P)
C**** Note that only layers LS1 and below vary as a function of surface
C**** pressure. Routine should be called with LMAX=LM at start, and
C**** subsequentaly with LMAX=LS1-1
C**** Note Air mass is calculated in (kg/m^2)

C**** Fill in polar boxes
      IF (have_south_pole) P(2:IM,1) = P(1,1)
      IF (have_north_pole) P(2:IM,JM)= P(1,JM)
      Call HALO_UPDATE(grid, P)

!$OMP  PARALLEL DO DEFAULT(NONE)
!$OMP&    PRIVATE (I,J,L,PL,AML,PDSIGL,PEDNL,PMIDL)
!$OMP&    SHARED (J_0H, J_1H, I_0H, I_1H, LMAX, PLIJ, 
!$OMP&           PDSIG, PMID, PEDN, AM, PK, PEK, BYAM, SQRTP, P)
      DO J=J_0H,J_1H ! filling halo for P is faster than PDSIG
        DO I=I_0H,I_1H

          CALL CALC_VERT_AMP(P(I,J),LMAX,PL,AML,PDSIGL,PEDNL,PMIDL)

          DO L=1,MIN(LMAX,LM)
            PLIJ (L,I,J) = PL    (L)
            PDSIG(L,I,J) = PDSIGL(L)
            PMID (L,I,J) = PMIDL (L)
            PEDN (L,I,J) = PEDNL (L)
            AM   (L,I,J) = AML   (L)
            PK   (L,I,J) = PMIDL (L)**KAPA
            PEK  (L,I,J) = PEDNL (L)**KAPA
            BYAM (L,I,J) = 1./AM(L,I,J)
          END DO

          IF (LMAX.ge.LM) THEN
            PEDN(LM+1:LMAX+1,I,J) = PEDNL(LM+1:LMAX+1)
            PEK (LM+1:LMAX+1,I,J) = PEDN(LM+1:LMAX+1,I,J)**KAPA
          END IF
          SQRTP(I,J) = SQRT(P(I,J))
        END DO
      END DO
!$OMP  END PARALLEL DO

      RETURN
      END SUBROUTINE CALC_AMPK

      SUBROUTINE CALC_AMP(p,amp)
!@sum  CALC_AMP Calc. AMP: kg air*grav/100, incl. const. pressure strat
!@auth Jean Lerner/Max Kelley
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,ls1,dsig,psf,ptop
      USE GEOM, only : axyp
C****
      USE DIST_GRID_MOD, ONLY : grid, GET
      implicit none
      REAL*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: p
      REAL*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: amp
      integer :: j,l
c**** Extract domain decomposition info
      INTEGER :: I_0, I_1, J_0, J_1
      CALL GET(grid, I_STRT = I_0, I_STOP = I_1, J_STRT = J_0,
     &               J_STOP = J_1)


C
!$OMP  PARALLEL DO PRIVATE(J,L)
      DO L=1,LM
        IF(L.LT.LS1) THEN
ccc   do l=1,ls1-1
          do j=J_0,J_1
            amp(I_0:I_1,j,l) = p(I_0:I_1,j)*axyp(I_0:I_1,j)*dsig(l)
          enddo
ccc   enddo
        ELSE
ccc   do l=ls1,lm
          do j=J_0,J_1
            amp(I_0:I_1,j,l) = (psf-ptop)*axyp(I_0:I_1,j)*dsig(l)
          enddo
        END IF
ccc   enddo
      enddo
!$OMP  END PARALLEL DO
C
      return
C****
      end subroutine calc_amp

      SUBROUTINE CALC_TROP
!@sum  CALC_TROP (to calculate tropopause height and layer)
!@auth J. Lerner
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,t
      USE GEOM, only : imaxj
      !USE DIAG_COM, only : aij => aij_loc, ij_ptrop, ij_ttrop
      USE DYNAMICS, only : pk, pmid, PTROPO, LTROPO
      USE DOMAIN_DECOMP_1D, Only : grid, GET
      IMPLICIT NONE
      INTEGER I,J,L,IERR
      REAL*8, DIMENSION(LM) :: TL
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Find WMO Definition of Tropopause to Nearest L
!$OMP  PARALLEL DO PRIVATE (I,J,L,TL,IERR)
      do j=J_0,J_1
      do i=I_0,imaxj(j)
        do l=1,lm
          TL(L)=T(I,J,L)*PK(L,I,J)
        end do
        CALL TROPWMO(TL,PMID(1,I,J),PK(1,I,J),PTROPO(I,J),LTROPO(I,J)
     *       ,IERR)
        IF (IERR.gt.0) print*,"TROPWMO error: ",i,j
        !AIJ(I,J,IJ_PTROP)=AIJ(I,J,IJ_PTROP)+PTROPO(I,J)
        !AIJ(I,J,IJ_TTROP)=AIJ(I,J,IJ_TTROP)+TL(LTROPO(I,J))
      end do
      end do
!$OMP  END PARALLEL DO
      IF (have_south_pole) THEN
        PTROPO(2:IM,1) = PTROPO(1,1)
        LTROPO(2:IM,1) = LTROPO(1,1)
      END IF
      IF (have_north_pole) THEN
        PTROPO(2:IM,JM)= PTROPO(1,JM)
        LTROPO(2:IM,JM)= LTROPO(1,JM)
      END IF

      END SUBROUTINE CALC_TROP


      subroutine tropwmo(ptm1, papm1, pk, ptropo, ltropp,ierr)
!@sum  tropwmo calculates tropopause height according to WMO formula
!@auth D. Nodorp/T. Reichler/C. Land
!@+    GISS Modifications by Jean Lerner/Gavin Schmidt
!@ver  1.0
!@alg  WMO Tropopause Definition
!@+
!@+ From A Temperature Lapse Rate Definition of the Tropopause Based on
!@+ Ozone, J. M. Roe and W. H. Jasperson, 1981
!@+
!@+ In the following discussion the lapse rate is defined as -dT/dz.
!@+
!@+ The main features of the WMO tropopause definition are as follows:
!@+ * The first tropopause (i.e., the conventional tropopause) is
!@+   defined as the lowest level at which the lapse rate decreases to 2
!@+   K/km or less, and the average lapse rate from this level to any
!@+   level within the next higher 2 km does not exceed 2 K/km.
!@+ * If above the first tropopause the average lapse rate between any
!@+   level and all higher levels within 1 km exceed 3 K/km, then a
!@+   second tropopause is defined by the same criterion as under the
!@+   statement above. This tropopause may be either within or above the
!@+   1 km layer.
!@+ * A level otherwise satisfying the definition of tropopause, but
!@+   occuring at an altitude below that of the 500 mb level will not be
!@+   designated a tropopause unless it is the only level satisfying the
!@+   definition and the average lapse rate fails to exceed 3 K/km over
!@+   at least 1 km in any higher layer.
!@+ * (GISS failsafe) Some cases occur when the lapse rate never falls
!@+   below 2 K/km. In such cases the failsafe level is that where the
!@+   lapse rate first falls below 3 K/km. If this still doesn't work
!@+   (ever?), the level is set to the pressure level below 30mb.
!@+
      USE MODEL_COM, only : klev=>lm
      USE CONSTANT, only : zkappa=>kapa,zzkap=>bykapa,grav,rgas
      implicit none

      real*8, intent(in), dimension(klev) :: ptm1, papm1, pk
      real*8, intent(out) :: ptropo
      integer, intent(out) :: ltropp,ierr
      real*8, dimension(klev) :: zpmk, zpm, za, zb, ztm, zdtdz
!@param zgwmo min lapse rate (* -1) needed for trop. defn. (-K/km)
!@param zgwmo2 GISS failsafe minimum lapse rate (* -1) (-K/km)
!@param zdeltaz distance to check for lapse rate changes (km)
!@param zfaktor factor for caluclating height from pressure (-rgas/grav)
!@param zplimb min pressure at which to define tropopause (mb)
      real*8, parameter :: zgwmo  = -2d-3, zgwmo2=-3d-3,
     *     zdeltaz = 2000.0, zfaktor = -GRAV/RGAS, zplimb=500.
      real*8 zptph, zp2km, zag, zbg, zasum, zaquer, zptf
      integer iplimb,iplimt, jk, jj, kcount, ltset,l
      logical ldtdz
c****
c****  2. Calculate the height of the tropopause
c****  -----------------------------------------
      ltset = -999
      ierr=0
      iplimb=1
c**** set limits based on pressure
      do jk=2,klev-1
        if (papm1(jk-1).gt.600d0) then
          iplimb=jk
        else
          if (papm1(jk).lt.30d0) exit
        end if
      end do
      iplimt=jk
c****
c****  2.1 compute dt/dz
c****  -----------------
c****       ztm  lineare Interpolation in p**kappa
c****     gamma  dt/dp = a * kappa + papm1(jx,jk)**(kappa-1.)

      do jk=iplimb+1,iplimt       ! -1 ?????
        zpmk(jk)=0.5*(pk(jk-1)+pk(jk))

        zpm(jk)=zpmk(jk)**zzkap ! p mitte

        za(jk)=(ptm1(jk-1)-ptm1(jk))/(pk(jk-1)-pk(jk))
        zb(jk) = ptm1(jk)-(za(jk)*pk(jk))

        ztm(jk)=za(jk)*zpmk(jk)+zb(jk) ! T mitte
        zdtdz(jk)=zfaktor*zkappa*za(jk)*zpmk(jk)/ztm(jk)
      end do
c****
c****  2.2 First test: valid dt/dz ?
c****  -----------------------------
c****
      do 1000 jk=iplimb+1,iplimt-1

c**** GISS failsafe test
        if (zdtdz(jk).gt.zgwmo2.and.ltset.ne.1) then
          ltropp=jk
          ltset =1
        end if
c****
        if (zdtdz(jk).gt.zgwmo .and. ! dt/dz > -2K/km
     &       zpm(jk).le.zplimb) then ! zpm not too low
          ltropp = jk
          ltset = 1
c****
c****  2.3 dtdz is valid > something in German
c****  ----------------------------------------
c****    1.lineare in p^kappa (= Dieters neue Methode)

          zag = (zdtdz(jk)-zdtdz(jk+1))/
     &         (zpmk(jk)-zpmk(jk+1)) ! a-gamma
          zbg = zdtdz(jk+1) - zag*zpmk(jk+1) ! b-gamma
          if(((zgwmo-zbg)/zag).lt.0.) then
            zptf=0.
          else
            zptf=1.
          end if
          zptph = zptf*abs((zgwmo-zbg)/zag)**zzkap
          ldtdz=zdtdz(jk+1).lt.zgwmo
          if(.not.ldtdz) zptph=zpm(jk)
c****
c****  2.4 2nd test: dt/dz above 2km must not be lower than -2K/km
c****  -----------------------------------------------------------
c****
          zp2km = zptph + zdeltaz*zpm(jk)
     &         / ztm(jk)*zfaktor ! p at ptph + 2km
          zasum = 0.0           ! zdtdz above
          kcount = 0            ! number of levels above
c****
c****  2.5 Test until pm < p2km
c****  --------------------------
c****
          do jj=jk,iplimt-1
            if(zpm(jj).gt.zptph) cycle ! doesn't happen
            if(zpm(jj).lt.zp2km) goto 2000 ! ptropo valid
            zasum = zasum+zdtdz(jj)
            kcount = kcount+1
            zaquer = zasum/float(kcount) ! dt/dz mean
            if(zaquer.le.zgwmo) goto 1000 ! dt/dz above < 2K/1000
                                          ! discard it
          end do                ! test next level
          goto 2000
        endif
 1000 continue                  ! next level
 2000 continue

      if (ltset.eq.-999) then
        ltropp=iplimt-1  ! default = last level below 30mb
        print*,"In tropwmo ltropp not set, using default: ltropp ="
     *       ,ltropp
        write(6,'(12(I4,5F10.5,/))') (l,ptm1(l),papm1(l),pk(l),zdtdz(l)
     *       ,zpm(l),l=iplimb+1,iplimt-1)
        ierr=1
      end if
      ptropo = papm1(ltropp)
c****
      return
      end subroutine tropwmo

#ifndef CUBED_SPHERE
      module zonalmean_mod
      contains
      subroutine zonalmean_ij2ij(arr,arr_zonal)
c Computes zonal means of arr and stores the result in arr_zonal.
c Lat-lon version.
      use model_com, only : im
      use domain_decomp_1D, only : grid
      use geom, only : imaxj
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     arr,         ! input
     &     arr_zonal    ! output
      integer :: j
      do j=grid%j_strt,grid%j_stop
        arr_zonal(:,j)=sum(arr(1:imaxj(j),j))/real(imaxj(j),kind=8)
      enddo
      return
      end subroutine zonalmean_ij2ij
      end module zonalmean_mod
#endif

!If running SCM use dummy routines
#ifndef SCM


      function getTotalEnergy() result(totalEnergy)
!@sum  getTotalEnergy returns the sum of kinetic and potential energy.
!@auth Tom Clune (SIVO)
!@ver  1.0
      use GEOM, only: AXYP, AREAG
      use DOMAIN_DECOMP_1D, only: grid, GLOBALSUM, get
      REAL*8 :: totalEnergy
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     KEIJ,PEIJ,TEIJ
      INTEGER :: I,J
      integer :: I_0, I_1, J_0, J_1

      call get(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%i_strt
      I_1 = grid%i_stop

      call conserv_PE(PEIJ)
      call conserv_KE(KEIJ)
      DO J=J_0,J_1
      DO I=I_0,I_1
        TEIJ(I,J)= (KEIJ(I,J) + PEIJ(I,J))*AXYP(I,J)/AREAG
      ENDDO
      ENDDO

      CALL GLOBALSUM(grid, TEIJ, totalEnergy, ALL=.true.)

      end function getTotalEnergy

      subroutine addEnergyAsDiffuseHeat(deltaEnergy)
!@sum  addEnergyAsDiffuseHeat adds in energy increase as diffuse heat.
!@auth Tom Clune (SIVO)
!@ver  1.0
      use CONSTANT, only: sha, mb2kg
      use MODEL_COM, only: T, PSF, PMTOP, LM
      use DYNAMICS, only: PK
      use DOMAIN_DECOMP_1D, only: grid, get
      implicit none
      real*8, intent(in) :: deltaEnergy

      real*8 :: ediff
      integer :: l
      integer :: I_0, I_1, J_0, J_1

      call get(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      ediff = deltaEnergy / ((PSF-PMTOP)*SHA*mb2kg)

!$OMP  PARALLEL DO PRIVATE (L)
      do l=1,lm
        T(I_0:I_1,J_0:J_1,L)=T(I_0:I_1,J_0:J_1,L)
     &       -ediff/PK(L,I_0:I_1,J_0:J_1)
      end do
!$OMP  END PARALLEL DO

      end subroutine addEnergyAsDiffuseHeat

      SUBROUTINE DISSIP
!@sum DISSIP adds in dissipated KE (m^2/s^2) as heat locally
!@auth Gavin Schmidt
      USE MODEL_COM, only : t
      USE DYNAMICS, only : dke,kea,pk
      IMPLICIT NONE
C**** temporarily store latest KE in DKE array
      call calc_kea_3d(dke)
      dke(:,:,:) = dke(:,:,:) - kea(:,:,:)
      call addEnergyAsLocalHeat(DKE, T, PK)

      END SUBROUTINE DISSIP

C***** Add in dissipiated KE as heat locally
      subroutine addEnergyAsLocalHeat(deltaKE, T, PK)!, diagIndex)
!@sum  addEnergyAsLocalHeat adds in dissipated kinetic energy as heat locally.
!@sum  deltaKE is on the A grid (J/kg)
!@auth Tom Clune (SIVO)
!@ver  1.0
      use CONSTANT, only: SHA
      use GEOM, only: IMAXJ
      use MODEL_COM, only: LM
      use DOMAIN_DECOMP_1D, only: grid, get
      implicit none
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &     deltaKE,T
      real*8, dimension(lm,grid%i_strt_halo:grid%i_stop_halo,
     &                     grid%j_strt_halo:grid%j_stop_halo) :: PK
c      integer, optional, intent(in) :: diagIndex

      integer :: i, j, l
      real*8 :: ediff
      integer :: I_0, I_1, J_0, J_1

      call get(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%i_strt
      I_1 = grid%i_stop

!$OMP  PARALLEL DO PRIVATE(I,J,L,ediff)
      DO L=1,LM
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        ediff = deltaKE(I,J,L) / (SHA*PK(L,I,J))
        T(I,J,L)=T(I,J,L)-ediff
c        if (present(diagIndex)) then
c          CALL INC_AJL(I,J,L,diagIndex,-ediff)
c        end if
      END DO
      END DO
      END DO
!$OMP  END PARALLEL DO
      end subroutine addEnergyAsLocalHeat

C**** Calculate 3D vertical velocity (take SDA which has units
C**** mb*m2, needs division by physics time step)
C**** and convert to WSAVE, units of m/s):

      subroutine COMPUTE_WSAVE
      use CONSTANT, only: rgas, bygrav
      use DOMAIN_DECOMP_1D, only: grid, GET
      use GEOM, only: byaxyp
      use MODEL_COM, only: IM,JM,LM,DTsrc,T
      use DYNAMICS, only: wsave, sda,pk,pedn
      implicit none

      integer :: i, j, l
      integer :: I_0, I_1, J_0, J_1

      call get(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

!$OMP PARALLEL DO PRIVATE (l,i)
      do l=1,lm-1
        do j=J_0,J_1
        do i=I_0,I_1
         wsave(i,j,l)=sda(i,j,l)*byaxyp(i,j)*
     &   rgas*0.5*(T(i,j,l)*pk(l,i,j)+T(i,j,l+1)*
     &   pk(l+1,i,j))*bygrav/(DTsrc*pedn(l+1,i,j))
        end do
        end do
      end do
!$OMP END PARALLEL DO

      end subroutine COMPUTE_WSAVE

      SUBROUTINE COMPUTE_GZ(p,t,tz,gz)
!@sum  COMPUTE_GZ calculates geopotential on model levels.
!@auth Original development team
      USE model_com, only : im,jm,lm,ls1,ptop,psfmpt,dsig,sige,sig,
     &     zatmo
      USE DOMAIN_DECOMP_1D, only : GRID,GET
      USE CONSTANT, only : grav,rgas,kapa,bykapa,bykapap1,bykapap2
      USE GEOM, only : imaxj
      IMPLICIT NONE
      INTEGER I,J,L
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     p
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     &     t,tz,gz
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     phidn,pkdn
      REAL*8 PKE(LS1:LM+1)
      REAL*8 PDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP,DP,P0,
     &     BYDP,pkdnl,TZBYDP,dpk,dpkp,dpkpp,
     &     dphidt,dphidtz,dphimdt,dphimdtz

c**** Extract domain decomposition info
      INTEGER :: I_0, I_1, J_0, J_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid,
     &     I_STRT=I_0, I_STOP=I_1,
     &     J_STRT=J_0, J_STOP=J_1,
     &     HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      DO L=LS1,LM+1
        PKE(L)=(PSFMPT*SIGE(L)+PTOP)**KAPA
      END DO

      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          pkdn(i,j)=(p(i,j)+ptop)**kapa
          phidn(i,j)=zatmo(i,j)
        ENDDO
      ENDDO

      do l=1,ls1-1 ! sigma levels
        DO J=J_0,J_1
          DO I=I_0,IMAXJ(J)
            pdn=sige(l)*p(i,j)+ptop
            pup=sige(l+1)*p(i,j)+ptop
            pkpdn=pkdn(i,j)*pdn
            dp=dsig(l)*p(i,j)
            bydp=1./dp
            p0=sig(l)*p(i,j)+ptop
            pkup=pup**kapa
            pkpup=pkup*pup
            dpk = (pkdn(i,j)-pkup)*bykapa
            dpkp = (pkpdn-pkpup)*bykapap1
            dpkpp = (pkpdn*pdn-pkpup*pup)*bykapap2
            dphidt = dpk
            dphimdt = bykapa*(pkdn(i,j)-dpkp*bydp)
            dphidtz = dphidt*p0 -dpkp
            dphimdtz = dphimdt*p0 +bykapap1*(bydp*dpkpp-pkpdn)
            tzbydp = 2.*tz(i,j,l)*bydp
!**** CALCULATE PHI, MASS WEIGHTED THROUGHOUT THE LAYER
            gz(i,j,l) = phidn(i,j)
     &           +rgas*(dphimdt*t(i,j,l)+dphimdtz*tzbydp)
!**** CALCULATE PHI AT LAYER TOP (EQUAL TO BOTTOM OF NEXT LAYER)
            phidn(i,j) = phidn(i,j)
     &           +rgas*(dphidt*t(i,j,l)+dphidtz*tzbydp)
            pkdn(i,j) = pkup
          ENDDO
        ENDDO
        IF (have_south_pole) GZ(2:IM, 1,L)=GZ(1, 1,L)
        IF (have_north_pole) GZ(2:IM,JM,L)=GZ(1,JM,L)
      enddo

      DO L=LS1,LM ! constant-pressure levels
        pdn=sige(l)*psfmpt+ptop
        pup=sige(l+1)*psfmpt+ptop
        pkdnl=pke(l)
        pkup=pke(l+1)
        pkpdn=pkdnl*pdn
        dp=dsig(l)*psfmpt
        bydp=1./dp
        p0=sig(l)*psfmpt+ptop
        pkpup=pkup*pup
        dpk = (pkdnl-pkup)*bykapa
        dpkp = (pkpdn-pkpup)*bykapap1
        dpkpp = (pkpdn*pdn-pkpup*pup)*bykapap2
        dphidt = dpk
        dphimdt = bykapa*(pkdnl-dpkp*bydp)
        dphidtz = (dphidt*p0 -dpkp)*2.*bydp
        dphimdtz = (dphimdt*p0 +bykapap1*(bydp*dpkpp-pkpdn))*2.*bydp
        dphidt = dphidt*rgas
        dphimdt = dphimdt*rgas
        dphidtz = dphidtz*rgas
        dphimdtz = dphimdtz*rgas
        do j=j_0,j_1
          do i=i_0,imaxj(j)
            gz(i,j,l) = phidn(i,j)
     &           +dphimdt*t(i,j,l) + dphimdtz*tz(i,j,l)
            phidn(i,j) = phidn(i,j)
     &           +dphidt *t(i,j,l) + dphidtz *tz(i,j,l)
          enddo
        enddo
        IF (have_south_pole) GZ(2:IM, 1,L)=GZ(1, 1,L)
        IF (have_north_pole) GZ(2:IM,JM,L)=GZ(1,JM,L)
      ENDDO

      RETURN
      END SUBROUTINE compute_gz


!if running SCM end
#endif

      function nij_before_j0(j0)
#ifdef CUBED_SPHERE
      use resolution, only : im,jm
      use domain_decomp_1D, only : grid
#else
      use geom, only : imaxj
#endif
      implicit none
      integer :: nij_before_j0,j0
#ifdef CUBED_SPHERE
      nij_before_j0 = im*((grid%tile-1)*jm + (j0-1))
#else
      nij_before_j0 = SUM(IMAXJ(1:J0-1))
#endif
      return
      end function nij_before_j0

      function nij_after_j1(j1)
      use resolution, only : im,jm
#ifdef CUBED_SPHERE
      use domain_decomp_1d, only : grid
#else
      use geom, only : imaxj
#endif
      implicit none
      integer :: nij_after_j1,j1
#ifdef CUBED_SPHERE
      nij_after_j1 = im*((6-grid%tile)*jm + (jm-j1))
#else
      nij_after_j1 = SUM(IMAXJ(J1+1:JM))
#endif
      return
      end function nij_after_j1

      function nij_after_i1(i1)
      use resolution, only : im,jm
      implicit none
      integer :: nij_after_i1,i1
#ifdef CUBED_SPHERE
      nij_after_i1 = im-i1
#else
      nij_after_i1 = 0
#endif
      return
      end function nij_after_i1

      subroutine sample_latlon(imr,jmr,arr,arrout)
!@sum sample_latlon fills the domain-decomposed atmosphere-grid
!@+   array arrout using values sampled from the global-domain
!@+   input array arr that is on a latitude-longitude grid having
!@+   imr longitudes and jmr latitudes.  No interpolation is
!@+   performed; each output value is simply copied from the input
!@+   gridcell containing the longitude and latitude of the nominal
!@+   center point of the output gridcell.
      use domain_decomp_1d, only : grid
      use geom, only : lon2d,lat2d
      use constant, only : pi,twopi,radian
      implicit none
      integer :: imr,jmr
      real*8 arr(imr,jmr)
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: arrout
      integer :: i,j,ilon,jlat
      real*8 :: bydlon,bydlat

      bydlon = real(imr,kind=8)/twopi
      bydlat = real(jmr,kind=8)/pi
      if(jmr == 46) bydlat = 45d0/pi

      do j=grid%j_strt,grid%j_stop
      do i=grid%i_strt,grid%i_stop
        ilon = 1 + int( bydlon*lon2d(i,j) )
        jlat = int( bydlat*min(.5d0*pi - .001d0, abs(lat2d(i,j))) )
        if(lat2d(i,j) > 0.) then
          jlat = 1 + jlat + jmr/2
        else
          jlat = jmr/2 - jlat
        endif
        arrout(i,j) = arr(ilon,jlat)
      enddo
      enddo

      return
      end subroutine sample_latlon

      FUNCTION SLP(PS,TAS,ZS)
!@sum SLP estimates sea level pressure in the presence of topography
!@+   for better match to reanalyses.
      USE CONSTANT, only: bmoist, grav, rgas, by3
      IMPLICIT NONE
!@var PS surface pressure (mb)
!@var TAS surface temperature (K)
!@var ZS surface elevation (m)
      REAL*8, INTENT(IN) :: PS, TAS, ZS
      REAL*8 :: SLP, TSL, BETA, BZBYT, GBYRB, TASn

      IF (ZS.ne.0.) THEN
        TSL= TAS+BMOIST*ZS
        TASn=TAS
        BETA=BMOIST
        IF (TAS < 290.5 .and. TSL > 290.5) BETA= (290.5d0 - TAS)/ZS
        IF (TAS > 290.5 .and. TSL > 290.5) TASn = 0.5*(290.5d0 + TAS)
        IF (TAS < 255) TASn = 0.5*(255d0 + TAS)
        BZBYT=BETA*ZS/TASn
        GBYRB=GRAV/(RGAS*BETA)
        IF (BETA > 1d-6 ) THEN
          SLP=PS*(1.+BZBYT)**GBYRB
        ELSE
          SLP=PS*EXP((1.-0.5*BZBYT+BZBYT**(2.*by3))*GBYRB*BZBYT)
        END IF
      ELSE
        SLP=PS
      END IF
      RETURN
      END FUNCTION SLP

      SUBROUTINE conserv_PE(TPE)
!@sum  conserv_TPE calculates total atmospheric potential energy
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : sha,mb2kg
      USE MODEL_COM, only : im,jm,lm,t,p,ptop,zatmo
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pk,pdsig
      USE DOMAIN_DECOMP_1D, only : GET,GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: TPE
      INTEGER :: I,J,L
      INTEGER :: J_0,J_1,I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C****
C**** TOTAL POTENTIAL ENERGY (J/m^2)
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        TPE(I,J)=0.
        DO L=1,LM
          TPE(I,J)=TPE(I,J)+T(I,J,L)*PK(L,I,J)*PDSIG(L,I,J)
        ENDDO
        TPE(I,J)=(TPE(I,J)*SHA+ZATMO(I,J)*(P(I,J)+PTOP))*mb2kg
      ENDDO
      ENDDO
      IF(HAVE_SOUTH_POLE) TPE(2:im,1) =TPE(1,1)
      IF(HAVE_NORTH_POLE) TPE(2:im,JM)=TPE(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_PE
