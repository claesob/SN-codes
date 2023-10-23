      SUBROUTINE CROSSSECT
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameters.h'
      integer skl,skk      
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/SIK/SK(ncr,NE1:NE2)
      dimension ion(270)
      data ion/
     &2,     3,    28*0,
     &12,   13,   8,    9,    10,   11,   24*0,
     &18,   19,   20,   21,   22,   23,   24,   23*0,
     &14,   15,   16,   17,   7 ,   4 ,   5 ,   6 ,   22*0,
     &72,   73,   74,   75,   76,   77,   78,   79,   80,   81,
     &20*0,
     &39,   40,   41,   27*0,
     &25,   26,   27,   28,   29,   30,   31,   32,   33,   34,   
     &35,   36,   37,   38,   16*0,
     &56,   57,   58,   59,   60,   61,   62,   63,   64,   65,
     &66,   67,   68,   69,   70,   71,   14*0,
     &42,   43,   44,   45,   84,   85,   86,   87,   88,   89,
     &90,   91,   92,   93,   94,   15*0/
c      dimension nopr(13)
c      data nopr/3,10,11,23,24,5,6,80,81,37,38,70,71/
      do iz=1,100
         do j=jmin,jj
            SI(iz,j)=0.
            Sk(iz,j)=0.
         enddo
      enddo
      open(42,file='./ATDAT/reilman.cross3.dat',status='old')
      call enint
      call num_shell

      CALL CRVERN(1,1,1)
      CALL CRVERN(2,2,1)
      CALL CRVERN(6,6,3)
      CALL CRVERN(7,7,3)
      CALL CRVERN(8,8,3)
      CALL CRVERN(10,10,3)
      CALL CRVERN(11,11,4)
      CALL CRVERN(12,12,4)

      CALL CRVERN(13,4,5)
      CALL CRVERN(14,14,5)
      CALL CRVERN(16,16,5)
      CALL CRVERN(18,7,5)
      CALL CRVERN(20,3,7)
      CALL CRVERN(26,26,7)
      RETURN
      end




      SUBROUTINE SEAT(J,SI0R,ALR,SR,E0R,J0,SIG)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameters.h'
      REAL*4 SI0R,ALR,SR,E0R
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      SI0=DBLE(SI0R)
      AL=DBLE(ALR)
      S=DBLE(SR)
      E0=DBLE(E0R)
      SIG=0.
      IF(E1(J).LT.E0.AND.J0.LT.2) GOTO 100
      IF(J.LT.J0.AND.J0.GE.2) GOTO 100
      SIG=SI0*1.E-18*(AL*(E0/E1(J))**S+(1.-AL)*(E0/E1(J))**(S+1.))
  100 CONTINUE
      RETURN
      END

      SUBROUTINE SHELL(J,SI0R,SR,E0R,J0,SIG)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameters.h'
      REAL*4 SI0R,SR,E0R
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
        SI0=DBLE(SI0R)
        S=DBLE(SR)
        E0=DBLE(E0R)
      SIG=0.E0
      IF(E1(J).LT.E0.AND.J0.LT.2) GOTO 100
      IF(J.LT.J0.AND.J0.GE.2) GOTO 100
      SIG=SI0*1.D-18*(E0/E1(J))**S
  100 CONTINUE
      RETURN
      END 

      DOUBLE PRECISION FUNCTION CR(J,E0R,C1R,C2R,C3R,C4R,C5R,C6R,C7R,
     & S0R,JW)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameters.h'
      REAL*4 E0R,C1R,C2R,C3R,C4R,C5R,C6R,C7R,S0R
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/FRE/NINT,JMIN,JJ
      E0=DBLE(E0R)
      C1=DBLE(C1R)
      C2=DBLE(C2R)
      C3=DBLE(C3R)
      C4=DBLE(C4R)
      C5=DBLE(C5R)
      C6=DBLE(C6R)
      C7=DBLE(C7R)
      S0=DBLE(S0R)
      IF(JW.GT.2) J0=JW
      IF(JW.LE.2) J0=JW
      CR=0.
      IF(E1(J).LT.E0.AND.J0.LT.2) GOTO 100
      IF(J.LT.J0.AND.J0.GE.2) GOTO 100
      EL=LOG10(E1(J))
      S=C1+C2*EL+C4*EL**3+C6*EL**5+C3/EL**2+C5/EL**4+C7/EL**6
      CR=1.E-18*(E0/E1(J))**S0*10.**S
 100  CONTINUE
       RETURN
        END

      DOUBLE PRECISION FUNCTION CRB(J,E0R,C1R,C2R,C3R,C4R,C5R,C6R,C7R,
     & S0R,JW)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 E0R,C1R,C2R,C3R,C4R,C5R,C6R,C7R,S0R      
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/FRE/NINT,JMIN,JJ
      IF(JW.GT.2) J0=JW
      IF(JW.LE.2) J0=JW
      E0=DBLE(E0R)
      C1=DBLE(C1R)
      C2=DBLE(C2R)
      C3=DBLE(C3R)
      C4=DBLE(C4R)
      C5=DBLE(C5R)
      C6=DBLE(C6R)
      C7=DBLE(C7R)
      S0=DBLE(S0R)
      CRB=0.
      IF(E1(J).LT.E0.AND.J0.LT.2) GOTO 100
      IF(J.LT.J0.AND.J0.GE.2) GOTO 100
      EL=LOG10(E1(J))
      S=C1/EL**3+C2/EL**2+C3/EL+C4+C5*EL+C6*EL**2+C7*EL**3
      CRB=1.E-18*(E0/E1(J))**S0*10.**S
 100  CONTINUE
      RETURN
      END

      SUBROUTINE ENINT
C     ***********************************************************
C     *****
C     ENERGY INTERVALS IN EV
C        E(-1)    E(0)    E(1)    E(
C     *****
C     ***********************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      DATA WLEV/1.239854e4/,WLBAL/3704.9034d0/
c     INTERVALS FROM PASHEN LIMIT TO BALMER
      NPASH=30
      DEP=LOG10(13.598/3.4)/(NPASH+1)
      NBALM=60
      WLX=2500.
c!!!
      nbalm=100
      wlx=1.e4
      ELX=WLEV/WLX
      DEB=LOG10(13.598/ELX)/(NBALM+1)
      
      DWL1=100.
      WL=WLX
      DO J=-NBALM,jmin,-1
         IF(WL.GT.10000.) DWL1=0.020*WL
         IF(WL.GT.20000.) DWL1=0.040*WL
         wlP1=WL
         wl=WL+DWL1
         E(J)=WLEV/WL
         EP1=WLEV/WLP1
         E1(J)=SQRT(E(J)*EP1)
      enddo
300   CONTINUE
c     INTERVALS FROM BALMER LIMIT TO LYMAN
      DO J=1-NBALM,1
         E(J)=ELX*10.**(REAL(J+NBALM-1)*DEB)
         E1(J)=ELX*10.**((REAL(J)+NBALM-.5)*DEB)
      enddo
c!!!!!
c      DO 7857 J=JMIN,1
c      E(J)=3.4*10.**(REAL(J+13)*.04013733)
c      E1(J)=3.4*10.**((REAL(J)+13.5)*.04013733)
c 7857 CONTINUE
      N1=20
      DE1=0.2572/N1
      DO  J=2,2+N1-1
         E(J)=13.598*10.**((J-2)*DE1)
         E1(J)=13.598*10.**((J-1.5)*DE1)
      enddo
      J2=2+N1
      N2=25
      DE2=0.34502/N2
      DO J=J2,J2+N2-1
         E(J)=24.587*10.**((J-J2)*DE2)
         E1(J)=24.587*10.**((J-J2+0.5)*DE2)
      enddo
      J3=J2+N2
      N3=25
      DE3=0.4076629/N3
      DO 230 J=J3,J3+N3-1
      E(J)=54.416*10.**((J-J3)*DE3)
  230 E1(J)=54.416*10.**((J-J3+0.5)*DE3)
      J4=J3+N3
      DO 240 J=J4,J4+5
      E(J)=139.12*10.**((J-J4)*0.0526138)
  240 E1(J)=139.12*10.**((J-J4+0.5)*0.0526138)
      J5=J4+5-27
      E(28+J5)=280.
      E1(28+J5)=287.9
      E(29+J5)=296.
      E1(29+J5)=306.3
      E(30+J5)=317.
      E1(30+J5)=331.7
      E(31+J5)=347.
      E1(31+J5)=357.8
      E(32+J5)=369.
      E1(32+J5)=380.3
      E(33+J5)=392.
      E1(33+J5)=401.9
      E(34+J5)=412.
      E1(34+J5)=421.9
      E(35+J5)=432.
      E1(35+J5)=445.3
      E(36+J5)=459.
      E1(36+J5)=474.2
      E(37+J5)=490.
      E1(37+J5)=500.3
      E(38+J5)=511.
      E1(38+J5)=521.9
      E(39+J5)=533.
      E1(39+J5)=541.4
      E(40+J5)=550.
      E1(40+J5)=559.9
      E(41+J5)=570.
      E1(41+J5)=582.4
      E(42+J5)=595.
      E1(42+J5)=610.8
      E(43+J5)=627.
      E1(43+J5)=646.7
      E(44+J5)=667.
      E1(44+J5)=684.3
      E(45+J5)=702.
      E1(45+J5)=720.3
      E(46+J5)=739.
      E1(46+J5)=769.9
      E(47+J5)=802.
      E1(47+J5)=835.3
      DO 260 J=48+J5,65+J5
      E(J)=870.*10.**(REAL(J-48-J5)*0.0750)
  260 E1(J)=870.*10.**((REAL(J-J5)-47.5)*0.0750)
      DO 261 J=66+J5,JJ
      E(J)=1.6387E4*10.**(REAL(J-65-J5)*0.15)
  261 E1(J)=1.6387E4*10.**((REAL(J-J5)-64.5)*0.15)
2123  CONTINUE
c      WRITE(6,*)'  J      E     E1'
      DO J=JMIN,JJ
c      WRITE(6,98)J,E(J),E1(J)
      ENDDO
98    FORMAT(I5,2F10.3)
      RETURN
      END

      SUBROUTINE CROSS
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/TRES/ EL(100),EK(100)
      COMMON/SIK/SK(ncr,NE1:NE2)
      DIMENSION JLMIN(100),JKMIN(100)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      EK(1)=13.598
      EL(1)=13.598
      EK(2)=24.587
      EL(2)=24.587
      EK(3)=54.416
      EL(3)=54.416
      EK(4)=672.
      EL(4)=138.12
      EK(5)=739.
      EL(5)=739.
      EK(6)=870.
      EL(6)=870.
      EK(7)=627.
      EL(7)=113.9
      EK(8)=317.0
      EL(8)=47.9
      EK(9)=347.0
      EL(9)=64.5
      EK(10)=392.
      EL(10)=392.
      EK(11)=490.
      EL(11)=490.
      EL(12)=11.26
      EK(12)=280.
      EL(13)=24.38
      EK(13)=296.
      EL(14)=13.62
      EK(14)=533.
      EL(15)=35.12
      EK(15)=550.
      EL(16)=54.93
      EK(16)=570.
      EL(17)=77.41
      EK(17)=595.
      EL(18)=14.53
      EK(18)=395.
      EL(19)=29.60
      EK(19)=412.
      EL(20)=47.45
      EK(20)=432.
      EL(21)=77.47
      EK(21)=459.
      EL(22)=97.89
      EK(22)=496.
      EL(23)=552.06
      EK(23)=552.06
      EL(24)=667.03
      EK(24)=667.03
      EL(25)=8.15
      EK(25)=1870.
      EL(26)=16.35
      EK(26)=1880.
      EL(27)=33.49
      EK(27)=1910.
      EL(28)=45.14
      EK(28)=1930.
      EL(29)=166.8
      EK(29)=1950.
      EL(30)=205.1
      EK(30)=1980.
      EL(31)=246.5
      EK(31)=2010.
      EL(32)=303.2
      EK(32)=2050.
      EL(33)=351.
      EK(33)=2090
      EL(34)=401.
      EK(34)=2140.
      EL(35)=476.
      EK(35)=2210.
      EL(36)=523.
      EK(36)=2247.
      EL(37)=2438.
      EK(37)=2438.
      EL(38)=2673.
      EK(38)=2673.
      EL(39)=7.646
      EK(39)=1320.
      EL(40)=15.
      EK(40)=1340.
      EL(41)=80.
      EK(41)=1360.
      EL(42)=7.870
      EK(42)=6900.
      EL(43)=16.16
      EK(43)=6900.
      EL(44)=30.65
      EK(44)=6900.
      EL(45)=54.4
      EK(45)=6900.
      EL(46)=5.99
      EK(46)=1514.
      EL(47)=18.8
      EK(47)=1549.
      EL(48)=28.4
      EK(48)=1549.
      EL(49)=120.
      EK(49)=1531.
      EL(50)=6.113
      EK(50)=4000.
      EL(51)=11.871
      EK(51)=5470.
      EL(52)=50.91
      EK(52)=5470.
      EL(53)=5.139
      EK(53)=1648.
      EL(54)=47.286
      EK(54)=1648.
      EL(55)=71.64
      EK(55)=1648.
      EL(56)=10.360
CF      EK(56)=3494.
      EL(57)=23.33
CF      EK(57)=3494.
      EL(58)=34.83
CF      EK(58)=3494.
C
CK> S I - S XVI
C
C      EL(56)=10.287
      EK(56)=2447.606
C      EL(57)=20.733
      EK(57)=2460.713
C      EL(58)=32.498
      EK(58)=2476.220
      EL(59)=45.391
      EK(59)=2493.880
      EL(60)=71.727
      EK(60)=2513.518
      EL(61)=88.470
      EK(61)=2535.537
      EL(62)=278.414
      EK(62)=2560.151
      EL(63)=328.439
      EK(63)=2625.932
      EL(64)=381.067
      EK(64)=2696.165
      EL(65)=436.123
      EK(65)=2770.605
      EL(66)=493.420
      EK(66)=2849.185
      EL(67)=552.670
      EK(67)=2931.664
      EL(68)=645.228
      EK(68)=3017.675
      EL(69)=706.591
      EK(69)=3096.581
      EL(70)=3198.541
      EK(70)=3198.541
CK  EL,EK(71) is not correct!!!
      EL(71)=3494.
      EK(71)=3494.
CK<
CF      EL(72)=21.564
CF      EK(72)=1362.
CF      EL(73)=40.962
CF      EK(73)=1362.
C
CK>  NE I - NE X
C
      EL(72)=21.564
      EK(72)=857.042
      EL(73)=40.419
      EK(73)=881.953
      EL(74)=64.155
      EK(74)=912.009
      EL(75)=90.843
      EK(75)=946.896
      EL(76)=120.216
      EK(76)=986.283
      EL(77)=151.997
      EK(77)=1029.952
      EL(78)=203.879
      EK(78)=1077.539
      EL(79)=239.976
      EK(79)=1123.254
      EL(80)=1184.674
      EK(80)=1184.674
CK  EL,EK (81) is not correct!!!
      EL(81)=1184.674
      EK(81)=1184.674
CK<
C
C    AR I - AR II
C
      EL(82)=15.759
      EK(82)=4426.
      EL(83)=27.629
      EL(83)=4426.
C
CK>  FE V - FE XV
C
      EL(84)=73.302
      EK(84)=7084.4
      EL(85)=96.853
      EK(85)=7111.9
      EL(86)=122.347
      EK(86)=7142.1
      EL(87)=149.639
      EK(87)=7174.9
      EL(88)=233.561
      EK(88)=7210.2
      EL(89)=262.827
      EK(89)=7251.0
      EL(90)=292.910
      EK(90)=7293.9
      EL(91)=323.741
      EK(91)=7338.4
      EL(92)=355.293
      EK(92)=7384.7
      EL(93)=387.526
      EK(93)=7433.0
      EL(94)=453.148
      EK(94)=7482.8
      EL(95)=52.42
C     READ K-SHELL ENERGIES FROM FILE. DATA FROM REILMAN & MANSON FILE.
C     H AND HE LIKE EXCLUDED
      OPEN(41,FILE='./ATDAT/kshell.dat',status='old')
      DO K=1,100
        READ(41,*,END=66)KQ,ET
        EK(KQ)=ET
      ENDDO
66    CONTINUE
      DO 300 J=JMIN,JJ
C     O VIII
      IF(E(J).GE.870.) THEN
        SK(6,J)=0.109E-18*(1.287*(870./E1(J))**2.95-.287*(870./E1(J))
     A                                                      **3.95)
      ENDIF

C     C VI
      IF(E(J).GE.490.) THEN
        SIA=1.287*(490./E1(J))**2.95-.287*(490./E1(J))**3.95
        SK(11,J)=0.194E-18*SIA
      ENDIF

C     O III 1D excited state (HENRY 72)
      IF(J.GE.46.AND.J.LE.75) THEN
        SI(95,J)=3.79E-18*(2.777*(52.42/E1(J))**3.-
     &                        1.777*(52.42/E1(J))**3.5)
      ENDIF      

      IF(E(J).GE.667.) THEN
        SK(24,J)=0.142E-18*(1.287*(667./E1(J))**2.95-.287*(667.
     &                                          /E1(J))**3.95)
      ENDIF
C     CA I  FROM SCOTT ET AL J PHYS 16, 3945, REILMAN &.. AND SHINE
      SIGB=0.
      SIRES=0.
      EP=(E1(J)-6.63)/4.2E-2
      IF(E1(J).GT.6.113) SIRES=7.93E-22*(EP+484.)/(1.+EP**2)
      IF(E1(J).GT.7.2.AND.E1(J).LT.9.25) SIGB=(7.38-0.69*E1(J))*1.E-20
      IF(E1(J).GT.9.25.AND.E1(J).LT.15.) SIGB=(.61+.481*(E1(J)-9.)**2
     &-.0413*(E1(J)-9.)**2)*1.E-19
      SI(50,J)=SIRES+SIGB
      SK(50,J)=0.
      IF(E(J).GT.4000.) SK(50,J)=6.85E-20*(4000./E1(J))**2.43
C     APPROX. FIT TO REILMAN AND MANSON
      IF(E1(J).GT.11.87.AND.E1(J).LT.35.) S=0.20*(11.87/E1(J))**0.0001
      IF(E1(J).GE.35..AND.E1(J).LT.45.) S=0.10*(35./E1(J))**1.
      IF(E1(J).GE.45..AND.E1(J).LT.160.) S=1.5
      IF(E1(J).GE.160..AND.E1(J).LT.390.) S=1.15*(160./E1(J))**1.66
      IF(E1(J).GE.390..AND.E1(J).LT.4000.) S=2.5*(390./E1(J))**2.18
      IF(E1(J).GT.4000.) S=6.85E-2*(4000./E1(J))**2.43
      SI(51,J)=S*1.E-18
      SK(51,J)=0.
      IF(E(J).GT.4000.) SK(51,J)=6.85E-20*(4000./E1(J))**2.43
C
C     NA I FROM BUTLER AND MENDOZA J PHYS B 16:L707 (83)
C
      ERY=(E1(J)-5.139)/13.61
      IF(ERY.GT.0..AND.ERY.LT..068) SI(53,J)=12.-177.*ERY
      IF(ERY.GT.0.068.AND.ERY.LT..093) SI(53,J)=0.01
      IF(ERY.GT..093.AND.ERY.LT..562) SI(53,J)=9.06*ERY+10.36-1.044/ERY
      IF(ERY.GT..562.AND.ERY.LT.2.) SI(53,J)=-7.99*ERY+22.57-2.52/ERY
C     ABOVE 100 EV FROM REILMAN AND MANSON
      IF(ERY.GT.2.AND.E1(J).LT.100) SI(53,J)=700.
      IF(E1(J).GT.100) SI(53,J)=700./(E1(J)/100.)**2.439
      SI(53,J)=SI(53,J)*1.E-20
      IF(E1(J).GT.1200.) SK(53,J)=0.2133E-18/(E1(J)/1200.)**2.70
C
C     NA II FROM REILMAN AND MANSON
C
      SI(54,J)=0.
      IF(E1(J).GE.47.29.AND.E1(J).LE.100.) SI(54,J)=7.0E-18
      IF(E1(J).GE.100.) SI(54,J)=7.0E-18*(100./E1(J))**2.34
C
C   Ne I-IV from Reilman and Manson, Ap.J. Suppl. 40, 815 (1979)
C
      XX=DLOG10(E1(J))

  300 CONTINUE
C     SET ALL LIMITS OF THE CROSS SECTIONS TO THE NEAREST ENERGY BIN
      DO K=1,100
        DO J=JMIN,JJ
          IF(EL(K).GT.E(J).AND.EL(K).LT.E(J+1)) THEN
            EMEAN=(E(J)+E(J+1))/2.
            IF(EL(K).GT.EMEAN) THEN 
              JLMIN(K)=J+1
            ELSE
              JLMIN(K)=J
            ENDIF
          ENDIF
          IF(EK(K).GT.E(J).AND.EK(K).LT.E(J+1)) THEN
            EMEAN=(E(J)+E(J+1))/2.
            IF(EK(K).GT.EMEAN) THEN 
              JKMIN(K)=J+1
            ELSE
              JKMIN(K)=J
            ENDIF
          ENDIF
        ENDDO
        DO J=JMIN,JJ
          IF(J.LT.JLMIN(K)) THEN
            SI(K,J)=0.
          ENDIF
          IF(J.LT.JKMIN(K)) THEN
            SK(K,J)=0.
          ENDIF
        ENDDO
        ILQ=0
        IKQ=0
        DO J=JMIN,JJ
          IF(SI(K,J).NE.0.AND.ILQ.EQ.0) THEN
            ILQ=1
          ENDIF
          IF(SK(K,J).NE.0.AND.IKQ.EQ.0) THEN
            IKQ=1
          ENDIF
        ENDDO
      ENDDO
 91   FORMAT(3I5,5F10.2)

      RETURN
      END



      SUBROUTINE CRVERN(iel,imax,ns)
c ns  numer of shells 5 for Si, S, Ar, 7 for Fe       
      IMPLICIT REAL*8(A-H,O-Z)
      integer iel,imax,ns
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/FRE/NINT,JMIN,JJ
c      COMMON/CSI_vern/GSV(14,5,NE1:NE2)
      COMMON/CSI_vern2/GSV2(30,30,7,NE1:NE2)
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/SIK/SK(ncr,NE1:NE2)
      common/ninns/ninn(30)
      common/ntot/ntot(30)
      common/ph1/ph1(6,30,30,7)
      common/ethresh/ethresh(30,30,7)
      common/nshells/nsh(30,30)
      DO  ion=1,imax
c     # of electrons in ion 
         nel=iel+1-ion
         nout=ntot(nel)
         if(iel.eq.nel.and.iel.gt.18)nout=7
         if(iel.eq.(nel+1).and.(iel.eq.20.or.iel.eq.21.or.iel.eq.22.
     &        or.iel.eq.25.or.iel.eq.26))nout=7
         nsh(iel,ion)=nout
         DO  is=1,nout
            ethresh(iel,ion,is)=ph1(1,iel,nel,is)
c            write(6,*)'ethresh(iel,ion,is)',iel,ion,is,
c     &           ethresh(iel,ion,is)
            do j=jmin,jj
               GSV2(iel,ion,is,j)=0.
               call phfit2(iel,nel,is,e1(j),gsv2(iel,ion,is,j))
            enddo
         enddo
      enddo
      DO ion=1,imax
         if(iel==1.and.ion.le.1) then
            ix=0            
         elseif(iel==2.and.ion.le.2) then
            ix=0
         elseif(iel==6.and.ion.le.4) then
            if(ion.le.2) then
               ix=11
            elseif(ion.ge.3) then
               ix=5
            endif
         elseif(iel==7.and.ion.le.7) then
            ix=17
         elseif(iel==8.and.ion<=4) then
            ix=13
         elseif(iel==8.and.ion==5) then
            ix=2
         elseif(iel==8.and.ion.ge.6) then
            ix=-2
         elseif(iel==10.and.ion.le.10) then
            ix=71
         elseif(iel==11.and.ion.le.3) then
            ix=52
         elseif(iel==12.and.ion.le.3) then
            ix=38
         elseif(iel==13.and.ion.le.4) then
            ix=45
         elseif(iel==14.and.ion.le.14) then
            ix=24
         elseif(iel==16.and.ion.le.10) then
            ix=55
         elseif(iel==18.and.ion.le.2) then
            ix=81
         elseif(iel==18.and.ion.ge.3) then
            ix=92
         elseif(iel==20.and.ion.le.3) then
            ix=49
         elseif(iel==26.and.ion.le.4) then
            ix=41
         elseif(iel==26.and.ion.le.15) then
            ix=79
         elseif(iel==26.and.ion.ge.16) then
            ix=85
         endif
         do j=jmin,jj
            SK(ion+ix,J)=GSv2(iel,ion,1,J)
            SI(ion+ix,J)=0.
            DO MI=2,ns
               SI(ion+ix,J)=SI(ion+ix,J)+GSv2(iel,ion,MI,J)
            enddo
         enddo
      enddo
      RETURN
      END
      
      




      subroutine crossint(z)
      implicit real*8(a-h,o-z)
      integer z
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      parameter (nz1=2,nz2=2,nz3=6,nz4=2,nz5=6,nz6=10,nz7=2)
      dimension sg1(26,ne1:ne2),et1(26),sg2(26,ne1:ne2),et2(26),
     $          sg3(26,ne1:ne2),et3(26),sg4(26,ne1:ne2),et4(26),
     $          sg5(26,ne1:ne2),et5(26),sg6(26,ne1:ne2),et6(26),
     $          sg7(26,ne1:ne2),et7(26)
      dimension nl(7)
      data nl/2,2,6,2,6,10,2/
      nk=0
      nls=0
      nlp=0
      nms=0
      nmp=0
      nmd=0
      nns=0
      nnp=0
      if(z.gt.2) then
       nk=2
      else
       nk=z
      endif
      if(z.gt.4) then
       nls=2
      else
       nls=z-2
      endif
      if(z.gt.10) then
       nlp=6
      else
       nlp=z-4
      endif
      if(z.gt.12) then
       nms=2
      else
       nms=z-10
      endif
      if(z.gt.18) then
        nmp=6
        if(z.le.20) then
          nmd=0
          nns=z-18
        elseif(z.le.28) then
          nmd=z-20
          nns=2
          if(z.eq.24) then
c           Cr
            nmd=5
            nns=1
          endif
        endif
      else
       nmp=z-12
      endif
      iz1=z
      iz2=iz1-nns
      iz3=iz2-nmd
      iz4=iz3-nmp
      iz5=iz4-nms
      iz6=iz5-nlp
      iz7=iz6-nls
      nz=z-1
      call manrd(sg1,et1,iz1,sg2,et2,iz2,sg3,et3,iz3,sg4,et4,iz4,
     $                 sg5,et5,iz5,sg6,et6,iz6,sg7,et7,iz7,nz,nsh)
      return
      end

      subroutine manrd(sg1,et1,nz1,sg2,et2,nz2,sg3,et3,nz3,sg4,et4,nz4,
     $                 sg5,et5,nz5,sg6,et6,nz6,sg7,et7,nz7,nz,nsh)
      IMPLICIT REAL*8(A-H,O-Z)
      character*80 head
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/FRE/NINT,JMIN,JJ
c
      common/sigma/ sg(26,ne1:ne2,7),et(26,7)
c
      dimension sg1(26,ne1:ne2),et1(26),sg2(26,ne1:ne2),et2(26),
     $          sg3(26,ne1:ne2),et3(26),sg4(26,ne1:ne2),et4(26),
     $          sg5(26,ne1:ne2),et5(26),sg6(26,ne1:ne2),et6(26),
     $          sg7(26,ne1:ne2),et7(26)
      dimension etmp(200),stmp(200)
      dimension nshell(26),nesh(7)
      dimension nzz(7)
c
      data nshell/2*1,2*2,6*3,2*4,6*5,6*6,2*7/
      data nesh/2,2,6,2,6,6,2/
c
      nzz(1)=nz1
      nzz(2)=nz2
      nzz(3)=nz3
      nzz(4)=nz4
      nzz(5)=nz5
      nzz(6)=nz6
      nzz(7)=nz7
c
c      print 9322,nz,nsh
9322  format (1h , ' in manrd ',2i4)
c
c     step through ion satges
      do 5 mn=1,nz
c
         nelec=nz1+1-mn
         njmx=nshell(nelec)
c         print 9389,mn,nelec,njmx
9389     format (' mn, nel, nj' ,3i4)
c
c        fill the threshold temporary
         et(mn,1)=et1(mn)
         et(mn,2)=et2(mn)
         et(mn,3)=et3(mn)
         et(mn,4)=et4(mn)
         et(mn,5)=et5(mn)
         et(mn,6)=et6(mn)
         et(mn,7)=et7(mn)
c
c        read the ion heading
         read (42,9902)head
9902     format (a)
c
c        step through subshells
         do 6 nj=1,njmx
c
            nelec=nzz(nj)+1-mn
            nelec=min0(nelec,nesh(nj))
            nelec=max0(nelec,0)
            enelec=float(nelec)
c
c           read the threshold.
            read (42,9903)swtch1,swtch2
            et(mn,nj)=swtch1
9903        format (f7.3,3x,f3.1)
            if (swtch1.lt.(-0.9)) go to 5
            if (swtch2.gt.(0.9)) go to 10
c
c           read from manson table
            do 1 jk=1,500
               jkk=jk
               read (42,9901)etmp(jk),stmp(jk)
9901           format (40x,e10.4,18x,e10.4)
               if (etmp(jk).ge.4999.) go to 11
1           continue
11          continue
            numrd=jkk
c
c           interpolate into epi grid
            jkk=1
            do 2 kl=jmin,jj
               epii=e1(kl)
               sg(mn,kl,nj)=0.
               if (epii.lt.et(mn,nj)) go to 2
               if (jkk.ge.(numrd-1)) go to 3
98                if (epii.le.etmp(jkk+1)) go to 99
                     jkk=jkk+1
                     if (jkk.lt.(numrd-1)) go to 98
99                continue
                  sgtmp=stmp(jkk+1)+(stmp(jkk)-stmp(jkk+1))
     $               *(epii-etmp(jkk+1))/(etmp(jkk)-etmp(jkk+1))
                  go to 22
3              continue
               alfa=log(stmp(jkk)/stmp(jkk-10))/
     &                  log(etmp(jkk)/etmp(jkk-10))
               sgtmp=stmp(jkk+1)*(epii/etmp(jkk+1))**alfa
22             continue
                  sg(mn,kl,nj)=sgtmp*(1.e-18)
c                 add Compton
                  csc=1.-(5.11e+5)*et(mn,nj)/(epii*epii)
                  csc=dmax1(csc,-1.d0)   
                  sigmn=(6.65e-25)*(1.+csc*(3.+csc*csc)/4.)/2.
                  sigmn=dmax1(sigmn,0.d0)
                  sg(mn,kl,nj)=dmax1(sg(mn,kl,nj),enelec*sigmn)
2           continue
c
c           use cross sections from neutral values with apprpriate treshold
            go to 6
10          continue
            do 13 kl=jmin,jj
               epii=e1(kl)
               sg(mn,kl,nj)=0.
               if (epii.lt.et(mn,nj)) go to 13
               sg(mn,kl,nj)=sg(1,kl,nj)
13          continue
c
6        continue
c
         read (42,9902)
c
c
5     continue
c
c        fill output arrays
         do 64 ml=jmin,jj
            do 65 ll=1,nz1
65             sg1(ll,ml)=sg(ll,ml,1)
            do 66 ll=1,nz2
66             sg2(ll,ml)=sg(ll,ml,2)
            do 67 ll=1,nz3
67             sg3(ll,ml)=sg(ll,ml,3)
            do 68 ll=1,nz4
68             sg4(ll,ml)=sg(ll,ml,4)
            do 69 ll=1,nz5
69             sg5(ll,ml)=sg(ll,ml,5)
            do 70 ll=1,nz6
70             sg6(ll,ml)=sg(ll,ml,6)
            do 71 ll=1,nz7
71             sg7(ll,ml)=sg(ll,ml,7)
64       continue
c
      return
      end


