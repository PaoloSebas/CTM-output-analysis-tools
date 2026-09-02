!---addFJX_H2CO.f  -- with user supplied subroutine that supplies X-section x Q-yield
!---generates fast-JX 18-bin X-sections (mprather modified file,
!---Sebastianelli 05/2024)

      implicit none

c defining indexes
      integer, parameter :: NB_ = 100      !--- in wave-bins.dat
      integer, parameter :: NS_ = 40000    !--- solar-p05nm-UNCI.dat
      integer, parameter :: NJ_ = 18       !--- will be used in CLOUD-J
c defining ARRAYS
      real*8 SRB(15,NJ_)   !Solar radiation in the bin
      real*8, dimension(NB_+1) :: WBIN
      real*8, dimension(NB_) :: FBIN, ABIN, XBIN
      real*8, dimension(NJ_) :: FFBIN,AABIN,XXBIN
      real*8 W(NS_),F(NS_)
c defining variables
      real*8 W1,W2,WW, XT
      real*8 XCH2O, QCH2O

      integer IBINJ(NS_)
      integer IJX(NB_), ITT
      integer NB,I,J,J1,J2,K,NB77,I778
      integer     TTT(2)
      character*1 ISXP
      character*6 TITLNEW
      character*6 TITLTBL(2)
      Character*20 TITL77

      open (1, file='wavel-bins.dat', status='OLD')
        SRB(:,:) = 0.d0
        read(1,'(i5)') NB
        if (NB .gt. NB_) stop
        read(1,'(5x,f8.3)') (WBIN(I), I=1,NB+1)
        read(1,*)
        read(1,*)
        read(1,'(2x,15f5.1)') ((SRB(I,J),I=1,15),J=1,8)
        read(1,*)
        read(1,'(5x,i5)') (IJX(I),I=16,NB)
      close (1)

      open (2, file='solar-p05nm-UCI.dat', status='OLD')
        read(2,*)
        read(2,*)
        read(2,'(f10.4,e10.3)') (W(J),F(J), J=1,NS_)
      close (2)

        IBINJ(:) = 0
       do I=1,NB
         W1 = WBIN(I)
         W2 = WBIN(I+1)
        do J=1,NS_
          if (W(J) .gt. W1) goto 11
        enddo
          J = NS_ + 1
   11     J1 = J
        do J=J1,NS_
          if (W(J) .gt. W2) goto 12
        enddo
          J = NS_ + 1
   12     J2 = J-1
        do J=J1,J2
          IBINJ(J) = I
        enddo
       enddo

!!!!this binning does not interpolate and is OK for large bins
!     it has 7% error in the very short wavel S-R bins of pratmo.

!not extrapolated, uses Q-yld at 300K min
        TTT(1) = 223.
        TTT(2) = 298.
        TITLTBL(1) = 'H2COa '
        TITLTBL(2) = 'H2COa '

! We have only two temperatures
        do K = 1,2
        XT = TTT(K)

!---now ready to do any flux-weighted means over the 77-pratmo bins
         FBIN(:) = 0.d0
         ABIN(:) = 0.d0
         I = 0

! primary high-resolution wavelength loop - generate input for pratmo reference J's

        do J = 1,NS_
         I = IBINJ(J)
         if (I .gt. 0) then
          WW = W(J)

           call X_CH2O (WW, XT, XCH2O)
           call Q_CH2O (WW, XT, QCH2O)

           FBIN(I) = FBIN(I) + F(J)
           ABIN(I) = ABIN(I) + F(J)*XCH2O*QCH2O


         endif
        enddo

        do I=1,NB
         if (FBIN(I) .gt. 0.d0) then
           NB77 = I
           ABIN(I) = ABIN(I)/FBIN(I)
         endif
        enddo

!---get Flsn
!---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
         FFBIN(:) = 0.d0
         AABIN(:) = 0.d0
         !XXBIN(:) = 0.d0
        do I=16,NB
         J = IJX(I)
         FFBIN(J) = FFBIN(J) + FBIN(I)
         AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)
         enddo

        do I=1,15
         do J=1,NJ_
           FFBIN(J) = FFBIN(J) + FBIN(I)*SRB(I,J)
           AABIN(J) = AABIN(J) + FBIN(I)*ABIN(I)*SRB(I,J)
         enddo
        enddo
        do J=1,NJ_
         if (FFBIN(J) .gt. 0.d0) AABIN(J) = AABIN(J)/FFBIN(J)
        enddo

! save UCI fast-JX v68 data bins for 'FJX_spec.dat'

        TITLNEW = TITLTBL(1)
        if (K .le. 2) then

         write(6,'(a6,a1,i3,1p,6e10.3/10x,6e10.3/10x,6e10.3)')
     &     TITLTBL(1),' ',TTT(K), AABIN

         endif

        enddo

      stop
      end

      subroutine X_CH2O (WW, TT, XCH2O)

c---JPL_2010 standard tables.
c---interpolates H2CO cross section vs. wavelength and temperature.
c---    H2CO is linear, limited to range 220 to 294K (THIS MUST BE
c---    CHECKED)

      implicit none
      real*8, intent(in) :: WW, TT
      real*8, intent(out):: XCH2O !

      integer IW,Z
      real*8 X1,X2,X1p,X2p,FT,FTX,FW

C----JPL 2010 Table 4D-3.
      real*4, parameter, dimension(1:107)  :: WNM =
     &[250, 251, 252, 253, 254, 255,
     & 256, 257, 258, 259, 260, 261,
     & 262, 263, 264, 265, 266, 267,
     & 268, 269, 270, 271, 272, 273,
     & 274, 275, 276, 277, 278, 279,
     & 280, 281, 282, 283, 284, 285,
     & 286, 287, 288, 289, 290, 291,
     & 292, 293, 294, 295, 296, 297,
     & 298, 299, 300, 301, 302, 303,
     & 304, 305, 306, 307, 308, 309,
     & 310, 311, 312, 313, 314, 315,
     & 316, 317, 318, 319, 320, 321,
     & 322, 323, 324, 325, 326, 327,
     & 328, 329, 330, 331, 332, 333,
     & 334, 335, 336, 337, 338, 339,
     & 340, 341, 342, 343, 344, 345,
     & 346, 347, 348, 349, 350, 351,
     & 352, 353, 354, 355, 356]


      real*4, parameter, dimension(1:107)  :: X223K =
     &[0.2550, 0.2030, 0.3380, 0.2880, 0.3410, 0.4480,
     & 0.6280, 0.4360, 0.3020, 0.6170, 0.6010, 0.6590,
     & 0.5970, 1.0860, 0.9430, 0.5220, 0.5350, 1.3670,
     & 1.2410, 0.9820, 0.9540, 1.9500, 1.4230, 0.8010,
     & 0.6500, 2.1700, 2.6050, 1.5560, 1.0190, 2.4690,
     & 2.3450, 1.5560, 0.9660, 0.7200, 4.3350, 4.0640,
     & 2.0780, 1.1440, 3.2050, 3.2290, 1.1570, 1.8490,
     & 0.7870, 3.1470, 7.2270, 4.0310, 2.4640, 1.3500,
     & 4.2570, 3.1590, 0.9310, 1.6490, 0.8690, 3.0470,
     & 7.2660, 4.7090, 4.2880, 1.7560, 1.3780, 3.2890,
     & 1.7220, 0.4600, 1.1990, 0.9130, 5.6450, 5.5410,
     & 2.5040, 5.8090, 3.1250, 0.9560, 1.1880, 1.6160,
     & 0.7200, 0.3240, 0.8640, 1.5230, 6.9120, 4.3290,
     & 1.1820, 3.1450, 3.8760, 1.3810, 0.3330, 0.2110,
     & 0.1610, 0.0970, 0.1240, 0.3810, 1.9080, 5.5490,
     & 3.1090, 0.9340, 0.4950, 1.9420, 1.2690, 0.4150,
     & 0.1150, 0.0460, 0.0810, 0.0420, 0.0410, 0.0880,
     & 0.7280, 2.3200, 1.6460, 0.6780, 0.1510]

c      real*4, parameter, dimension(1:89)  :: X223K =
c     &[0.2590, 0.2050, 0.3360, 0.2900, 0.3430, 0.4520,
c     & 0.6300, 0.4500, 0.3120, 0.6190, 0.6070, 0.6610,
c     & 0.6070, 1.0740, 0.9510, 0.5380, 0.5410, 1.3530,
c     & 1.2390, 0.9980, 0.9660, 1.9300, 1.4380, 0.8190,
c     & 0.6640, 2.1300, 2.5750, 1.5840, 1.0410, 2.4310,
c     & 2.3350, 1.5640, 0.9780, 0.7200, 4.2050, 4.0370,
c     & 2.1020, 1.1560, 3.1350, 3.2110, 1.1830, 1.8310,
c     & 0.8050, 3.0730, 7.0730, 4.0890, 2.4960, 1.3700,
c     & 4.1830, 3.1810, 0.9950, 1.6110, 0.8350, 2.9930,
c     & 7.1940, 4.7710, 4.2920, 1.8040, 1.3820, 3.2310,
c     & 1.7580, 0.4620, 1.1810, 0.8910, 5.6550, 5.5790,
c     & 2.5760, 5.7710, 3.1750, 0.9940, 1.1920, 1.5840,
c     & 0.7220, 0.3300, 0.8580, 1.5570, 6.8280, 4.4110,
c     & 1.2580, 3.0950, 3.8440, 1.4390, 0.3590, 0.2170,
c     & 0.1570, 0.0970, 0.1280, 0.3850, 1.9320]

      real*4, parameter, dimension(1:107)  :: X298K =
     &[0.2570, 0.2040, 0.3370, 0.2890, 0.3420, 0.4500,
     & 0.6290, 0.4430, 0.3070, 0.6180, 0.6040, 0.6600,
     & 0.6020, 1.0800, 0.9470, 0.5300, 0.5380, 1.3600,
     & 1.2400, 0.9900, 0.9600, 1.9400, 1.4300, 0.8100,
     & 0.6570, 2.1500, 2.5900, 1.5700, 1.0300, 2.4500,
     & 2.3400, 1.5600, 0.9720, 0.7200, 4.2700, 4.0500,
     & 2.0900, 1.1500, 3.1700, 3.2200, 1.1700, 1.8400,
     & 0.7960, 3.1100, 7.1500, 4.0600, 2.4800, 1.3600,
     & 4.2200, 3.1700, 0.9630, 1.6300, 0.8520, 3.0200,
     & 7.2300, 4.7400, 4.2900, 1.7800, 1.3800, 3.2600,
     & 1.7400, 0.4610, 1.1900, 0.9020, 5.6500, 5.5600,
     & 2.5400, 5.7900, 3.1500, 0.9750, 1.1900, 1.6000,
     & 0.7210, 0.3270, 0.8610, 1.5400, 6.8700, 4.3700,
     & 1.2200, 3.1200, 3.8600, 1.4100, 0.3460, 0.2140,
     & 0.1590, 0.0966, 0.1260, 0.3830, 1.9200, 5.5000,
     & 3.1500, 0.9780, 0.5040, 1.9200, 1.2700, 0.4360,
     & 0.1190, 0.0440, 0.0760, 0.0380, 0.0360, 0.0890,
     & 0.7310, 2.2800, 1.6500, 0.6960, 0.1480]

      XCH2O = 0.d0

      if (WW .lt. WNM(1) .or. WW .gt. WNM(107)) goto 2

      do Z=1,107
      if (WW .ge. WNM(Z)) IW=Z
      enddo

      FW = WW - WNM(IW)
      FW = max(0.d0, min(1.d0, FW))

      X1p = X223K(IW)
      X2p = X298K(IW)

      X1 = (X1p + FW*(X223K(IW+1)-X223K(IW)))*1.0d-20
      X2 = (X2p + FW*(X298K(IW+1)-X298K(IW)))*1.0d-20
      FTX = (TT - 223.d0) / (298.d0 - 223.d0)
      FT = max (0.d0, min (1.d0, FTX))
      XCH2O = X1 + FT*(X2-X1)

    2 continue
      return
      end


      subroutine Q_CH2O (WW, TT, QCH2O)

c---JPL_2010 standard tables.
c---interpolates H2CO quantum yield vs. wavelength

c---JPL 2006 recommendations
      implicit none
      real*8, intent(in) :: WW, TT
      real*8, intent(out):: QCH2O

      integer IW
      real*8 FW,Q1,FT

C----------JPL 2010 Table 4D-4 pg297.
      real*4, parameter, dimension(107)  :: WNQ =
     &[250.0, 251.0, 252.0, 253.0, 254.0, 255.0, 256.0,
     & 257.0, 258.0, 259.0, 260.0, 261.0, 262.0, 263.0,
     & 264.0, 265.0, 266.0, 267.0, 268.0, 269.0, 270.0, 271.0,
     & 272.0, 273.0, 274.0, 275.0, 276.0, 277.0, 278.0, 279.0,
     & 280.0, 281.0, 282.0, 283.0, 284.0, 285.0, 286.0, 287.0,
     & 288.0, 289.0, 290.0, 291.0, 292.0, 293.0, 294.0, 295.0,
     & 296.0, 297.0, 298.0, 299.0, 300.0, 301.0, 302.0, 303.0,
     & 304.0, 305.0, 306.0, 307.0, 308.0, 309.0, 310.0, 311.0,
     & 312.0, 313.0, 314.0, 315.0, 316.0, 317.0, 318.0, 319.0,
     & 320.0, 321.0, 322.0, 323.0, 324.0, 325.0, 326.0, 327.0,
     & 328.0, 329.0, 330.0, 331.0, 332.0, 333.0, 334.0, 335.0,
     & 336.0, 337.0, 338.0, 339.0, 340.0, 341.0, 342.0, 343.0,
     & 344.0, 345.0, 346.0, 347.0, 348.0, 349.0, 350.0, 351.0,
     & 352.0, 353.0, 354.0, 355.0, 356.0]

      real*4, parameter, dimension(107)  :: Q300 =
     &[0.000, 0.308, 0.307, 0.306, 0.305, 0.304, 0.304,
     & 0.303, 0.303, 0.304, 0.307, 0.312, 0.318, 0.325,
     & 0.333, 0.343, 0.354, 0.365, 0.377, 0.390, 0.404, 0.418,
     & 0.433, 0.448, 0.464, 0.479, 0.495, 0.512, 0.528, 0.544,
     & 0.560, 0.576, 0.591, 0.606, 0.620, 0.633, 0.645, 0.657,
     & 0.669, 0.680, 0.690, 0.700, 0.710, 0.718, 0.726, 0.734,
     & 0.740, 0.746, 0.751, 0.755, 0.758, 0.761, 0.762, 0.762,
     & 0.762, 0.760, 0.758, 0.754, 0.749, 0.744, 0.737, 0.729,
     & 0.720, 0.709, 0.698, 0.685, 0.671, 0.656, 0.639, 0.622,
     & 0.603, 0.583, 0.561, 0.539, 0.515, 0.489, 0.463, 0.435,
     & 0.406, 0.375, 0.343, 0.310, 0.276, 0.240, 0.203, 0.165,
     & 0.126, 0.085, 0.043, 0.000, 0.000, 0.000, 0.000, 0.000,
     & 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     & 0.000, 0.000, 0.000, 0.000, 0.000]

      if (WW .lt. WNQ(1) .or. WW .gt. WNQ(107)) goto 2

      IW = (WW - 249.d0)
      IW = max(1, min(107, IW))
      FW = WW - WNQ(IW)
      FW = max(0.d0, min(1.d0, FW))
      Q1 = Q300(IW) + FW*(Q300(IW+1)-Q300(IW))
      QCH2O = Q1

   2  continue
      return
      end


~
~
~
~
~
~
~
~
~
~
