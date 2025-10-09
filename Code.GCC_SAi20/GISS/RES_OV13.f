C****   
!@sum  RES_OV13 Resolution info for 13 ocean layer
!@auth Gary Russel/Larissa Nazarenko
!@ver  1.0
C****   
C**** RES_OV13.f    2008/07/03
C****
      Module OCEANRES
      Use RESOLUTION, Only: imo=>IM,jmo=>JM
      Use MODEL_COM,  Only: IVSP,IVNP
      Implicit None
C****
C**** Vertical resolution of ocean
C****
      Integer*4,Parameter ::
!     *      IMO = 72,
!     *      JMO = 46, 
     *      LMO = 13, !  maximum number of ocean layers in a column
     *  LMO_MIN =  2, !  minimum number of ocean layers in a column
     *    LSRPD =  3, !  deepest layer for penetrating solar radiation
     *    MAXGL =  5  !  maximum no. layers for depositing glac. melt
C*** 
C***
      Real*8,Parameter ::

     *     dZO(LMO)  =   !  approximate thicknesses (m) of ocean layers
     *                 (/ 12d0, 18d0, 27d0, 40.5d0, 60.75d0, 91.125d0, 
     *                   136.6875d0, 205.03125d0, 307.546875d0, 
     *                   461.3203125d0, 691.98046875d0, 
     *                   1037.970703125d0, 1556.9560546875d0 /) 

C****
C****  L  dZO    ZE   ZOC      L  dZO    ZE   ZOC
C****  =  ===    ==   ===      =  ===    ==   ===
C****  1   12    12     6      7  137   386   318
C****  2   18    30    21      8  205   591   488
C****  3   27    57    44      9  308   899   745
C****  4   41    98    77     10  461  1360  1129
C****  5   61   158   128     11  692  2052  1706
C****  6   91   249   204     12 1038  3090  2571
C****                         13 1557  4647  3868
C****

      EndModule OCEANRES
