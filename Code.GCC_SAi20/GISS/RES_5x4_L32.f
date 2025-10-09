C****   
!@sum  RES_5x4_L32 Resolution info for 4x5 horizontal
!!       and 32 vertical ocean layer
!@auth Gary Russel/Larissa Nazarenko
!@ver  1.0
C****   
C**** RES_5x4_L32.f    2008/08/08
C****
      Module OCEANRES
      Implicit None
C****
C**** Vertical resolution of ocean
C****
      Integer*4,Parameter ::
     *      IMO = 72,
     *      JMO = 46, 
     *      LMO = 32, !  maximum number of ocean layers in a column
     *  LMO_MIN =  2, !  minimum number of ocean layers in a column
     *    LSRPD =  4, !  deepest layer for penetrating solar radiation
     *    MAXGL =  6, !  maximum no. layers for depositing glac. melt (<200m)
     *   NOCEAN =  1  !  NDYNO must be multiple of 2*NOCEAN
     * , NORDER=4     !  order of Alternating Binomial Filter (must be even)
C*** 
C***
      Real*8,Parameter ::
     *    AKHMIN = 1.5d8    ! minimum horizontal diffusion
     *  , AKHFAC = 1d0      ! factor to multiply built-in scaling for k_diff
     *  , fkph   = 0.3d0    ! vertical diffusion coefficient (cm**2/sec)
     *  , fkpm   = 10.*fkph ! vertical viscosity coefficient (cm**2/sec)
     *  , OABFUX=0.         ! coef. for divergence filter in X dir. 
     *  , OABFVX=0.         ! coef. for vorticity filter in X dir.
     *  , by4tonv=0.        ! coef. for divergence filter in Y dir.
     *  , by4tonu=0.        ! coef. for vorticity filter in Y dir.
      Real*8,Parameter ::

     *     dZO(LMO) =   !  approximate thicknesses (m) of ocean layers
     *               (/ 12d0, 18d0, 26d0, 36d0,  48d0, 62d0, 78d0, 96d0,
     *               116d0,134d0,150d0,164d0, 176d0,186d0,194d0,200d0,
     *               204d0,206d0,206d0,206d0, 206d0,206d0,206d0,206d0,
     *               206d0,206d0,206d0,206d0, 206d0,206d0,206d0,206d0 /) 
C****
C**** Vertical resolution of ocean
C****
C****  L  dZO   ZOE   ZOC      L  dZO   ZOE   ZOC
C****  =  ===   ===   ===      =  ===   ===   ===
C****  1   12    12     6     17  204  1900  1798
C****  2   18    30    21     18  206  2100  2003
C****  3   26    56    43     19  206  2312  2209
C****  4   36    92    74     20  206  2518  2415
C****  5   48   140   116     21  206  2724  2621
C****  6   62   202   171     22  206  2930  2827
C****  7   78   280   241     23  206  3136  3033
C****  8   96   376   328     24  206  3342  3239
C****  9  116   492   434     25  206  3548  3445
C**** 10  134   626   559     26  206  3754  3651
C**** 11  150   776   701     27  206  3960  3857
C**** 12  164   940   858     28  206  4166  4063
C**** 13  176  1116  1028     29  206  4372  4269
C**** 14  186  1302  1209     30  206  4578  4475
C**** 15  194  1496  1399     31  206  4784  4681
C**** 16  200  1696  1596     32  206  4990  4887
C****
      EndModule OCEANRES
