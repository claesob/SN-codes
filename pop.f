
      SUBROUTINE POPCA(RS,TSIN,ZS,XQ,PHH,IFPOP)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      COMMON/ABUN/AB(20)
      COMMON/INUT/IUY
      COMMON/ABU/XA(2,nio)
      COMMON/HRAT/RHYD,ZMETAL,XRQ,HYR,HEAT,COOL
      COMMON/NION/IONQ
      COMMON/A25/DHDT,DCODT,DZDT,DPDT
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),PH(NL),DRDT(NL),RECT(5,NL),PHET(NL)
      COMMON/MULTI/IHMUL,IHEMUL,IOMUL,ICAMUL,IFEMUL
      COMMON/INI/INITH,INIT,INITCA,INITHE,INITFE
      COMMON/ITER/ITE
      COMMON/SECX/CSEC(20),DCSDX(20),CISEC(20),DCISDX
      COMMON/NBA/NBACK
      COMMON/IND/IK
      COMMON/ELEC/DEL(MD)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &     BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/ELDEN/IDENS
      COMMON/PHEAT/PHE(NL)
      COMMON/BIN/B1IN,ELD
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/RECAL/RECO(NL)
      COMMON/A12/ITMAX
      DIMENSION XI(NL),XINC(NL),AQ(NLP1,NLP1),AAB(NLP1,NLP1),
     &                                     A(NLP1,NLP1),BS(NL)
      DATA CSAHA/2.0708E-16/
      ifpop=0
      ITE=-1
      ION=1
      IONQ=1
      IIT=0
      NMIN=3
      NMAX=3
      CALL ATDAT
      DO 15 NQ=NMIN,NMAX
      N=NQ
      NP1=N+1
      NP1H=NP1
      NVAR=NP1
      ICO=0
      IW=1
C     IF(IK.EQ.2) INIT=1
      TS=TSIN
      IF(IIT.GT.10) STOP
      RCM=1.E15*RS
      IF(ITE.EQ.-1) DS=DE(RS)
      IF(N.GT.NMIN) GOTO 600
      DO 1963 KK=1,N
      IF(INITCA.ne.1) GOTO 1963
      BB(KK)=BOLCA(KK)*(TS/TOLDCA)**1.5*EXP((E00-E(KK))*
     &1.1609E4*(1./TOLDCA-1./TS))
1963  BS(KK)=BB(KK)
      BB(NP1)=BOLCA(NP1)
      BS(NP1)=BB(NP1)
      IF(IK.NE.2) GOTO 600
 600  NRC=NP1+1
1399  CONTINUE
      EPS1=1.D-300
      EPS2=1.E-3
      IPRINT=2
      J=1
        IQQ=0
      IF(IK.EQ.2.AND.IUY.EQ.1) XELEC=ZS
      IF(IK.EQ.2.AND.IUY.EQ.1) XN1=0.
      IF(IK.NE.2) XELEC=DEL(IK-1)
      XELEC=ZMETAL
C     CALL RAD(TS,XELEC)
C     CALL SPEC
      DO 9 ITERA=1,ITMAX
      IF(ITERA.EQ.ITMAX) WRITE(6,*)' NO CONVERGENCE FOR CA II!!'
      RCM=1.E15*RS
      IF(ITE.EQ.-1) DS=DE(RS)
C
C     CALCULATE COLL EXCITATION RATES
C
        IQW=0
C
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITYXN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C
      C1=CSAHA
      DEN=DE(RS)
      G1=2.
      G2=10.
      G3=6.
      TEV=TS/1.1609E4
      XN1=G1*BB(NP1)*DEL(IK)*C1*DEN*EXP(11.87/TEV)*BB(1)/TS**1.5
      PHN1=XN1*PHE(1)
      XN2=G2*BB(NP1)*DEL(IK)*C1*DEN*EXP(10.174/TEV)*BB(2)/TS**1.5
      PHN2=XN2*PHE(2)
      XN3=G3*BB(NP1)*DEL(IK)*C1*DEN*EXP(8.729/TEV)*BB(3)/TS**1.5
      PHN3=XN3*PHE(3)
      PHH=PHN1+PHN2+PHN3
      CAII=XN1+XN2+XN3
      DEL(IK)=AB(11)*CAII*BB(NP1)+ZMETAL
      IF(INITCA.EQ.1) DEL(IK)=AB(11)*XQ+ZMETAL
      XN1=G1*BB(NP1)*DEL(IK)*C1*DEN*EXP(11.87/TEV)*BB(1)/TS**1.5
      XN2=G2*BB(NP1)*DEL(IK)*C1*DEN*EXP(10.174/TEV)*BB(2)/TS**1.5
      XN3=G3*BB(NP1)*DEL(IK)*C1*DEN*EXP(8.729/TEV)*BB(3)/TS**1.5
      IF(ITERA.GT.1) XELEC=AB(11)*CAII*BB(NP1)+ZMETAL
      CALL SECEX(XELEC)
      CALL COLLEX(ION,TS)
      CALL HPCA(RS,TS)
      IF(N.GT.NMIN) GOTO 6382
      IF(IUY.NE.1) GOTO 6382
      IF(ITERA.NE.1) GOTO 6382
c!!      IF(INITCA.NE.1) GOTO 6382
c!!
cnots      PH(1)=0.
      T4=TS/1.E4
C     FACTOR 0.7 TO AGREE WITH RATES FROM HRATE
      DENEL=XELEC*DEN
      ALO=0.7*6.8E-13/T4**.8
      RI=ALO*EXP((E(1)-E00)/TEV)*TS**1.5/(CSAHA*G(1))
      BB(1)=RI/(PH(1)+DENEL*CI(1)+CISEC(5)/DEN)
      BB(NP1)=1./(1.+ALO*XELEC/(PH(1)/DEN+CI(1)*XELEC+CISEC(5)/DEN**2))
      XIII=BB(NP1)
      A2=RECO(2)*DENEL*XIII + C(1,2)*DENEL*XN1
      A3=RECO(3)*DENEL*XIII + C(1,3)*DENEL*XN1
      BE32=1.
      BE21=1.
      BE31=0.
      AL2=BE32*AS(3,2) + C(3,2)*DENEL
      AL3=C(2,3)*DENEL
      B2= PH(2) + (C(2,3) + C(2,1))*DENEL + BE21*AS(2,1)
      B3= PH(3) + (C(3,2) + C(3,1))*DENEL + BE31*AS(3,1) + 
     &                                              BE32*AS(3,2)
      X2=(B3*A2+AL2*A3)/(B2*B3-AL2*AL3)
      X3=(A3+AL3*X2)/B3
      ZI=XQ*AB(11)
      IDEP=0
      CALL MULTISIMPq(13,2,TS,XI,rl)
      BB(NP1)=XI(NP1)
      DENEL=DEL(IK)*DE(RS)
c use dep coeff?
      IF(IDEP.EQ.1) THEN
         DO J=1,N
            BB(J)=XI(J)
            XN(ION,J)=XI(J)*G(J)*BB(NP1)*DENEL*C1*EXP((E00-E(J))/TEV)
     &           /TS**1.5
         ENDDO
      ELSE
         DO J=1,N
            XN(ION,J)=XI(J)
            BB(J)=XI(J)/(G(J)*BB(NP1)*DENEL*C1*EXP((E00-E(J))/TEV)
     &           /TS**1.5)
         ENDDO
         itcon=1
         goto 999
      ENDIF
      INITCA=0
      DO NK=1,NP1
         BS(NK)=BB(NK)
      enddo
6382  CONTINUE
      IHLO=0
      ITCON=1
      IF(TS.LT.2000..OR.ICAMUL.EQ.0) GOTO 11
      IHLO=1
      IF(N.GT.NMIN.AND.ITERA.EQ.1) BB(NP1)=BB(NP1-1)
      IF(N.GT.NMIN.AND.ITERA.EQ.1) BB(N)=BB(N-1)
      IF(N.GT.NMIN.AND.ITERA.EQ.1) BS(N)=BB(N)
      IF(N.GT.NMIN.AND.ITERA.EQ.1) BS(NP1)=BB(NP1)
      CALL CALCN(ION,ITERA,RS,TS,ZMETAL,XQ,BB,A,NRC)
 7365 NEQ=NP1
      IF(IDENS.LT.0) NEQ=NP1-1
5343  continue
      DETER=SIMUL(NEQ,A,XINC,EPS1,1,NRC)
      IF(DETER.NE.0.) GOTO 3
3     ITCON=1
      DO  I=1,NVAR
         ERRO=XINC(I)/BB(I)
         IF(ABS(ERRO).GT.EPS2) ITCON=0
         write(6,*)'error ',i,xinc(i),bb(i),erro
         BB(I)=BB(I)+XINC(I)
      enddo
      xcaii=0.
      XELEC=XCAII*BB(NP1)*AB(11)+ZMETAL
        IRW=0
      DO 6736 I=1,NVAR
      IF(BB(I).GT.0.) GOTO 6736
 8392 FORMAT(1X,'BB',7E10.4)
 4828 FORMAT(1X,'ITER',I8)
C     BB(I)=BS(I)
      BB(I)=BB(I)-XINC(I)
        IF(0.5*ABS(XINC(I)).LT.BB(I)) BB(I)=BB(I)+XINC(I)/2.
        IF(0.5*ABS(XINC(I)).GT.BB(I)) BB(I)=BB(I)/2.
C     BB(I)=BB(I)/2.
        IF(I.EQ.NP1) GOTO 9473
C     IF(BB(I).LT.1.) BB(I)=1.
9473  IW=2
C     ICO=ICO+1
      IF(ICO.GT.10) IIT=IIT+1
      IF(ICO.GT.20) WRITE(6,*)'NO CONVERGENCE FOR CA II',ts,den,del(ik)
	if(ico.gt.20) ifpop=1
      IF(ICO.GT.20) goto 15
        IRW=1
C     GOTO 1399
 6736 CONTINUE
        ICO=ICO+1
C       IF(IRW.EQ.1) GOTO 1399
      IF(BB(NP1).GT.1.) BB(NP1)=1.
      IF(IDENS.LT.0) BB(NP1)=1.
      IF(ITCON.EQ.0) GOTO 9
      GOTO 11
9     CONTINUE
      IF(ITCON.EQ.0) STOP
 11    CONTINUE
C
C     CALC. THE FRACTIONAL POPULATIONS 
C
      C1=CSAHA
      DEN=DE(RS)
      TEV=TS/1.1609E4
      DO J=1,N
        XN(ION,j)=G(J)*BB(NP1)*DEL(IK)*C1*DEN*EXP((E00-E(J))/TEV)*
     &            BB(J)/TS**1.5
      ENDDO
 999  continue
      PHN1=XN(ION,1)*PHE(1)
      PHN2=XN(ION,2)*PHE(2)
      PHN3=XN(ION,3)*PHE(3)
      PHH=PHN1+PHN2+PHN3
      IF(XN1.LE.0.) WRITE(6,7258)BOLCA(1),TOLD,TS,C1,BB(NP1),RS,DEN,TEV
 7258 FORMAT(1X,'XN1=0',8E12.4)
 7259 FORMAT(1X,'X1,X2,X3',3E12.4)
        TOLDCA=TS
        DO 3767 IS=1,NVAR
3767    BOLCA(IS)=BB(IS)
15    CONTINUE
      IF(ITCON.EQ.0) IFPOP=1
      IF(ITCON.EQ.0) write(6,204)
      IF(IHLO.EQ.1) CALL HLOSS(1,BB,RS,TS,XELEC,XQ,HCA)
      RETURN
201   FORMAT(1X,'MATRIX ILL-CONDITIONED OR SINGULAR',E11.4)
202   FORMAT(1X,I8,E10.4,I8,8E12.5)
203   FORMAT(1X,'SUCCESSFULCONVERGENCE',I5,I5,8E12.5)
204   FORMAT(1X,'NO CONVERGENCE FOR CA II!')
      END

      SUBROUTINE POPO(RS,TSIN,ZS,XQ,PHO,IFPOP)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      COMMON/A19/EM(NL,NL),ESC(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      COMMON/ABUN/AB(20)
      COMMON/INUT/IUY
      COMMON/ABU/XA(2,nio)
      COMMON/HRAT/RHYD,ZMETAL,XRQ,HYR,HEAT,COOL
      COMMON/NION/IONQ
      COMMON/A25/DHDT,DCODT,DZDT,DPDT
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),PH(NL),DRDT(NL),RECT(5,NL),PHET(NL)
      COMMON/MULTI/IHMUL,IHEMUL,IOMUL,ICAMUL,IFEMUL
      COMMON/INI/INITH,INIT,INITCA,INITHE,INITFE
      COMMON/ITER/ITE
      COMMON/SECX/CSEC(20),DCSDX(20),CISEC(20),DCISDX
      COMMON/NBA/NBACK
      COMMON/IND/IK
      COMMON/ELEC/DEL(MD)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &     BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/ELDEN/IDENS
      COMMON/PHEAT/PHE(NL)
      COMMON/BIN/B1IN,ELD
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A12/ITMAX
      COMMON/LYBFLUOR/OIFLUOR,OIFLUORH,OIFLUORO,OILYB,OIREC
      DIMENSION XINC(NL),AQ(NLP1,NLP1),AAB(NLP1,NLP1),XI(NL),
     &                                    A(NLP1,NLP1),BS(NL)
      DATA CSAHA/2.0708E-16/
      ifpop=0
      ITE=-1
      ION=2
      IONQ=2
      N=9
      NMIN=N
      NMAX=N
      CALL ATDATO
      DO 15 NQ=NMIN,NMAX
      NP1=N+1
      NP1H=NP1
      NP2=N+2
      NVAR=NP2
      IIT=0
      ICO=0
      IW=1
c      IF(IK.EQ.2) INIT=1
      TS=TSIN
      IF(IIT.GT.10) STOP
      RCM=1.E15*RS
      IF(ITE.EQ.-1) DS=DE(RS)
      IF(N.GT.NMIN) GOTO 600
      DO 1963 KK=1,N
      IF(INIT.EQ.1) GOTO 1963
      BB(KK)=BOL(KK)*(TS/TOLDO)**1.5*EXP((E00-E(KK))*
     &1.1609E4*(1./TOLDO-1./TS))
1963  BS(KK)=BB(KK)
      BB(NP1)=BOL(NP1)
      BS(NP1)=BB(NP1)
C     BB(NP2)=TSIN-DBLE(IIT)*.5E3
C     BS(NP2)=TSIN-DBLE(IIT)*.5E3
      IF(IK.NE.2) GOTO 600
 600  NRC=NP1+1
1399  CONTINUE
      EPS1=1.D-300
      EPS2=1.E-3
      IPRINT=2
      J=1
      IQQ=0
      IF(IK.EQ.2.AND.IUY.EQ.1) XELEC=ZS
      IF(IK.EQ.2.AND.IUY.EQ.1) XN1=0.
      IF(IK.NE.2) XELEC=DEL(IK-1)
C     CALL RAD(TS,XELEC)
C     CALL SPEC
      DO 9 ITERA=1,ITMAX
      IF(ITERA.EQ.ITMAX) WRITE(6,*)' NO CONVERGENCE FOR O!!'
C     IF(ITE.EQ.1) TS=BB(NP2)
      RCM=1.E15*RS
      IF(ITE.EQ.-1) DS=DE(RS)
C
C     CALCULATE COLL EXCITATION RATES
C
        IQW=0
C
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITY XN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C
      C1=CSAHA
      DEN=DE(RS)
      G1=9.
      G2=5.
      G3=1.
      TEV=TS/1.1609E4
      IF(13.6/TEV.GT.700.) GOTO 777
      XN1=G1*BB(NP1)*DEL(IK)*C1*DEN*EXP(13.60/TEV)*BB(1)/TS**1.5
      PHN1=XN1*PHE(1)
      XN2=G2*BB(NP1)*DEL(IK)*C1*DEN*EXP(11.657/TEV)*BB(2)/TS**1.5
      PHN2=XN2*PHE(2)
C     XN3=G3*BB(NP1)*DEL(IK)*C1*DEN*EXP(9.446/TEV)*BB(3)/TS**1.5
      PHN3=XN3*PHE(3)
777   PHH=PHN1+PHN2+PHN3
      IF(ITERA.GT.1) XELEC=AB(3)*BB(NP1)*XQ+ZMETAL
      CALL COLLEX(ION,TS)
      CALL HORATE(ION,RS,TS)
11    XELOLD=XELEC
      CALL SECEX(XELEC)
      IF(N.GT.NMIN) GOTO 6382
      IF(IUY.NE.1) GOTO 6382
C!!      IF(INIT.NE.1) GOTO 6382
      IF(ITERA.NE.1) GOTO 6382
      do nk=1,n
      BB(nk)=nk*0.1
      enddo
      T4=TS/1.E4
C!!   OTS-APP. DOes not work for external rad field!!
cnots      PH(1)=0.
      ALO=0.
      DO J=1,N
         CALL RECOMB(ION,J,TS,AL,RI)
         ALO=ALO+AL
      ENDDO
      RI=ALO*EXP((E(1)-E00)/TEV)*TS**1.5/(CSAHA*G(1))
      BB(1)=RI/(PH(1)+DEN*XELEC*CI(1)+CISEC(2)/DEN)
      QQ=(PH(1)/DEN+CI(1)*XELEC+CISEC(2)/DEN**2)/(ALO*AB(3)*XQ)
      ZQ=ZS/(AB(3)*XQ)
      BB(NP1)=SQRT(QQ+(QQ+ZQ)**2/4.)-(ZQ+QQ)/2.
      XELEC=ZS+BB(NP1)*AB(3)*XQ
      IF(ABS(XELEC-XELOLD)/XELOLD.GT.0.001) GOTO 11
      INIT=0
      BB(2)=BB(1)/2.
      BB(3)=BB(1)/3.
      ZI=AB(3)*XQ
      IDEP=0
      CALL MULTISIMP(14,5,1,IDEP,ZI,TS,XI)
      BB(NP1)=XI(NP1)
      DENEL=DEL(IK)*DE(RS)
      IF(IDEP.EQ.1) THEN
         DO J=1,N
            BB(J)=XI(J)
            XN(ION,J)=XI(J)*G(J)*BB(NP1)*DENEL*C1*EXP((E00-E(J))/TEV)
     &           /TS**1.5
         ENDDO
      ELSE
         DO J=1,N
            XN(ION,J)=XI(J)
            EXT=DMIN1(700.D0,(E00-E(J))/TEV)
            BB(J)=XI(J)/(G(J)*BB(NP1)*DENEL*C1*EXP(EXT)/TS**1.5)
         ENDDO
      ENDIF
      DO 1827 NK=1,NP1
 1827 BS(NK)=BB(NK)
      IHLO=0
      ITCON=1
      IF(TS.LT.2800..OR.IOMUL.EQ.0) GOTO 111
6382  CONTINUE
      IHLO=1
      IF(N.GT.NMIN.AND.ITERA.EQ.1) BB(NP1)=BB(NP1-1)
      IF(N.GT.NMIN.AND.ITERA.EQ.1) BB(N)=BB(N-1)
      IF(N.GT.NMIN.AND.ITERA.EQ.1) BS(N)=BB(N)
      IF(N.GT.NMIN.AND.ITERA.EQ.1) BS(NP1)=BB(NP1)
c!!!!

      CALL CALCN(ION,ITERA,RS,TS,ZMETAL,XQ,BB,A,NRC)
7365  NEQ=NP2
      IF(ITE.EQ.-1) NEQ=NP1
      IF(ITE.EQ.1) GOTO 406
      DO 407 I=1,NP1
407   A(I,NP2)=A(I,NRC)
406   CONTINUE
      IF(IDENS.LT.0) NEQ=NP1-1
C     PRINT *,' MATRIX A'
      IF(IK.NE.999) GOTO 3922
      DO 3948 I=1,NVAR
 3948  WRITE(6,9376)(A(I,J),J=1,NRC)
9376  FORMAT(5E12.5)
3922  CONTINUE
 5343 continue
      DETER=SIMUL(NEQ,A,XINC,EPS1,1,NRC)
      IF(DETER.NE.0.) GOTO 3
3     ITCON=1
      IF(ITE.EQ.-1) NVAR=NP1
      DO 5 I=1,NVAR
      ERRO=XINC(I)/BB(I)
      IF(ABS(ERRO).GT.EPS2) ITCON=0
 5    BB(I)=BB(I)+XINC(I)
      write(6,922)itera,(xinc(i),i=1,nvar)
      write(6,922)itera,(bb(i),i=1,nvar)
 922  format('xinc,bb ',i5,1pe11.3,20e11.3)
      IRW=0
      DO 6736 I=1,NVAR
      IF(BB(I).GT.0.) GOTO 6736
      BB(I)=BB(I)-XINC(I)
      IF(0.5*ABS(XINC(I)).LT.BB(I)) BB(I)=BB(I)+XINC(I)/2.
      IF(0.5*ABS(XINC(I)).GT.BB(I)) BB(I)=BB(I)/2.
      IF(I.EQ.NP2) BB(NP2)=BB(NP2)*0.9
      IF(I.EQ.NP1) GOTO 9473
      IF(BB(I).LT.1.) BB(I)=1.
9473  IW=2
C     ICO=ICO+1
      IF(ICO.GT.10) IIT=IIT+1
      IF(ICO.GT.20) WRITE(6,*)'NO CONVERGENCE FOR O I  ',ts,den,del(ik)
      If(ico.gt.20) write(6,*)'ifpop',IFPOP
      IF(ICO.GT.20) IFPOP=1
      IF(ICO.GT.20) RETURN
      IRW=1
 6736 CONTINUE
      ICO=ICO+1
      IF(BB(NP1).GT.1.) BB(NP1)=1.
      IF(IDENS.LT.0) BB(NP1)=1.
      IF(ITCON.EQ.0) GOTO 9
      GOTO 111
9     CONTINUE
      IF(ITCON.EQ.0) STOP
111   CONTINUE
C
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITY) XN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C
      C1=CSAHA
      DEN=DE(RS)
      DEL(IK)=AB(3)*XQ*BB(NP1)+ZMETAL
      TEV=TS/1.1609E4
      DO 1928 J=1,N
      EXT=DMIN1(700.D0,(E00-E(J))/TEV)
1928  XN(2,J)=G(J)*BB(NP1)*DEL(IK)*C1*DEN*EXP(EXT)*BB(J)/
     &TS**1.5
      XN1=0.
      XN2=0.
      XN3=0.
      IF(13.6/TEV.GT.700.) GOTO 778
      XN1=G1*BB(NP1)*DEL(IK)*C1*DEN*EXP(13.60/TEV)*BB(1)/TS**1.5
      XN2=G2*BB(NP1)*DEL(IK)*C1*DEN*EXP(11.657/TEV)*BB(2)/TS**1.5
      XN3=G3*BB(NP1)*DEL(IK)*C1*DEN*EXP(9.446/TEV)*BB(3)/TS**1.5
778   PHN1=XN1*PHE(1)
      PHN2=XN2*PHE(2)
      PHN3=XN3*PHE(3)
      PHH=PHN1+PHN2+PHN3
      TOLDO=TS
      DO 3767 IS=1,NVAR
3767  BOL(IS)=BB(IS)
15    CONTINUE
      IF(ITCON.EQ.0) IFPOP=1
      IF(ITCON.EQ.0) write(6,204)
      XELEC=AB(3)*BB(NP1)*XQ+ZMETAL
      IF(IHLO.EQ.1) CALL HLOSS(2,BB,RS,TS,XELEC,XQ,HCA)
      IF(IFPOP.EQ.1.AND.IK.EQ.1) INIT=1
C     LYMAN BETA FLUORESENCE A LA KWAN AND KROLIK
      IF(XN(5,1).GT.0.) THEN
        IF(AB(1).GT.0.01) THEN
          OIFLUOR=1.41E6*6.65*ESC(7,8)/(5.55+1.1*ESC(7,8))
        ELSE
          OIFLUOR=0.
        ENDIF
        OIFLUORH=OIFLUOR*AB(3)*XQ*XN(2,1)/(AB(1)*XN(5,1))
        OIFLUORO=OIFLUOR*XN(5,3)/XN(5,1)
      ENDIF
      RETURN
201   FORMAT(1X,'MATRIX ILL-CONDITIONED OR SINGULAR',E11.4)
202   FORMAT(1X,I8,E10.4,I8,8E12.5)
204   FORMAT(1X,'NO CONVERGENCE FOR O!')
      END

      SUBROUTINE POPSI_I(RS,TSIN,ZS,XQ,PHH,IFPOP)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      common/abl/abn(15)
      common/timecheck/time,itime
      common/ionx/xion(md,14,27)
      parameter (nlp=30000)
      COMMON/EQUIH/FRH(2,500),WEH(500),SR(NFEL),WOBS(NL,NL)
      COMMON/INUT/IUY
      COMMON/HRAT/RHYD,ZMETAL,XRQ,HYR,HEAT,COOL
      COMMON/NION/IONQ
      COMMON/A25/DHDT,DCODT,DZDT,DPDT
      COMMON/IND/IK
      COMMON/ELEC/DEL(MD)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &       BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/ELDEN/IDENS
      COMMON/PHEAT/PHE(NL)
      COMMON/BIN/B1IN,ELD
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/FELEV/NFEII,nfei
      common/multisol/imultih,imultihe,imultio,imultica,imultifei,
     &     imultifeii
      DIMENSION XINC(NL),XI(NL),A(NLP1,NLP1),BS(NL)

      DATA CSAHA/2.0708E-16/
      ION=7
      IONQ=7

      nsii=65

      nsii=56

      n=nsii
      NMIN=nsii
      NMAX=NMIN

      CALL ATDATSI_I
      call si_i_recomb(tsin)
      ifpop=0
      ITE=-1
      IIT=0

      TS=TSIN

C
C     CALCULATE COLL EXCITATION RATES
C
      IQW=0
C
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITYXN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C

      call collsi_i(ts)

      ZI=1.

      CALL MULTISIMPq(10,1,TS,XI,rltot)
      
      return

      END

      subroutine collsi_i(te)
      implicit real*8 (a-h,o-z)
      include "parameters.h"
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      common/abl/abn(15)
      tev=te/1.1609d4
      te5000=te/5000
      t4=te/1.e4
      tesqrt=dsqrt(te)
c *** Si I ***
      ie=iesi
      ion=1
      do i=1,nl
        ci(i)=0.
        do j=1,nl
          c(i,j)=0.
        enddo
      enddo
c
      if (abn(10).gt.0.) then
c     Si I 129.68 68.474 my omegas guessed!! from o i
         o32=9.76e-6*(te-228)+3.46e-11*(te-228.)**2
         o31=3.39e-6*(te-326.)-2.90e-11*(te-326.)**2
         o21=1.89e-6*(te-326.)+8.00e-11*(te-326.)**2
         if(o21.lt.0.) o21=1.e-10
         if(o31.lt.0.) o31=1.e-10
         if(o32.lt.0.) o32=1.e-10
c
         et=(e(3)-e(2))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(2,3)=8.63e-6*o32*etexp/(g(2)*tesqrt)
         c(3,2)=g(2)*c(2,3)/(g(3)*etexp)
c
         et=(e(3)-e(1))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(1,3)=8.63e-6*o31*etexp/(g(1)*tesqrt)
         c(3,1)=g(1)*c(1,3)/(g(3)*etexp)
c
         et=(e(2)-e(1))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(1,2)=8.63e-6*o21*etexp/(g(1)*tesqrt)
         c(2,1)=g(1)*c(1,2)/(g(2)*etexp)
c
c     (Si I) 6591, 10995, 16360  omega for 3p - 1d from
c     m.s. pindzola, a.k. bhatia, a.temkin, phys.rev. 15, 35 (1977)
c     rest scaled from c i with scaling = 3
c      o21=3.97*t4**0.9
         o4_123=3.97*t4**0.9
         o41=o4_123/9.d0
         o42=o4_123/3.d0
         o43=o4_123*5.d0/9.d0
c
         et=(e(4)-e(1))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(1,4)=8.63e-6*o41*etexp/(g(1)*tesqrt)
         c(4,1)=g(1)*c(1,4)/(g(4)*etexp)
c
         et=(e(4)-e(2))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(2,4)=8.63e-6*o42*etexp/(g(2)*tesqrt)
         c(4,2)=g(2)*c(2,4)/(g(4)*etexp)
c
         et=(e(4)-e(3))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(3,4)=8.63e-6*o43*etexp/(g(3)*tesqrt)
         c(4,3)=g(3)*c(3,4)/(g(4)*etexp)
c
c      o31=3.*0.149*(te/5000.)**0.871
         o5_123=3.*0.149*(te5000)**0.871
         o51=o5_123/9.d0
         o52=o5_123/3.d0
         o53=o5_123*5.d0/9.d0
c
         et=(e(5)-e(1))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(1,5)=8.63e-6*o51*etexp/(g(1)*tesqrt)
         c(5,1)=g(1)*c(1,5)/(g(5)*etexp)
c
         et=(e(5)-e(2))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(2,5)=8.63e-6*o52*etexp/(g(2)*tesqrt)
         c(5,2)=g(2)*c(2,5)/(g(5)*etexp)
c
         et=(e(5)-e(3))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(3,5)=8.63e-6*o53*etexp/(g(3)*tesqrt)
         c(5,3)=g(3)*c(3,5)/(g(5)*etexp)
c      o32=3.*0.196*(te/5000.)**0.499
         o54=3.*0.196*(te5000)**0.499
c
         et=(e(5)-e(4))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(4,5)=8.63e-6*o54*etexp/(g(4)*tesqrt)
         c(5,4)=g(4)*c(4,5)/(g(5)*etexp)

         do j=6,n
            do i=1,j-1
               gi=g(i)
               gj=g(j)
               eij=abs(e(j)-e(i))
               et=eij/tev
               if(et.gt.650.) et=650.
               etexp=dexp(-et)
c     omij=omsi1(i,j)*(te5000**alfsi1(i,j))
c     For the remaining transitions omega=0.1. Just a wild guess!
c     Better for the lower transitions from  Vernazza, Avrett, & Loeser 1976, ApJS, 30, 1
c     see Cissis file
               omij=0.1
               c(i,j)=8.63e-6*omij*etexp/(gi*tesqrt)
               c(j,i)=gi*c(i,j)/(gj*etexp)
            enddo
         enddo

c Collisional ionization rate (cm**3 s**-1) from 
c Vernazza, Avrett, & Loeser 1976, ApJS, 30, 1
         et=(e00-e(1))/tev
         if(et.gt.650.) et=650.
         ci(1)=2.41d-8*(te5000**0.649)*dexp(-et)
         ci(2)=ci(1)
         ci(3)=ci(1)
c
         et=(e00-e(4))/tev
         if(et.gt.650.) et=650.
         ci(4)=0.27d-8*(te5000**0.688)*dexp(-et)
c
         et=(e00-e(5))/tev
         if(et.gt.650.) et=650.
         ci(5)=5.39d-8*(te5000**0.353)*dexp(-et)
c
         et=(e00-e(6))/tev
         if(et.gt.650.) et=650.
c change from CK to CF  numbering         
c     ci(6)=0.35d-8*(te5000**0.582)*dexp(-et)
         ci(7)=0.35d-8*(te5000**0.582)*dexp(-et)
         ci(8)=ci(7)
         ci(9)=ci(7)
c
         et=(e00-e(7))/tev
         if(et.gt.650.) et=650.
c     ci(7)=1.23d-8*(te5000**0.547)*dexp(-et)
         ci(10)=1.23d-8*(te5000**0.547)*dexp(-et)
c
         et=(e00-e(8))/tev
         if(et.gt.650.) et=650.
c     ci(8)=6.63d-8*(te5000**0.525)*dexp(-et)
         ci(16)=6.63d-8*(te5000**0.525)*dexp(-et)
         ci(17)=ci(16)
         ci(18)=ci(16)
      endif

      return
      end

      SUBROUTINE ATDATsi_i
      IMPLICIT REAL*8(A-H,O-Z)
      character*80 dum
      character*2 ch2
      character*4 ch4
      character*5 ch5
      SAVE
      include "parameters.h"
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      COMMON/SILEV/NSIII,nsi1
      common/initat/initfeii,initfei,initsi1 
      integer nlsi1
      parameter(nlsi1=66)
      common/csi_idat/gsi1(nlsi1)
      common/si1_reco/tesi1r(81),alsi1r(65,81),totrecsi1(81)
      DIMENSION GS(NL),ES(NL),AS(NL,NL),WLS(NL,NL),OMS(NL,NL)
      ntot=355
      nsi1=65
      N=Nsi1
      E00=8.15168
      IF(Initsi1.EQ.0) THEN
         OPEN(23,FILE='./ATDAT/Si_I_levels_Aij_clean.dat')
         rewind 23
         do i=1,ntot
            read(23,*)nn,wn,ji,gsi
            if(i<=nsi1) then
               gsi1(i)=gsi
               gs(i)=gsi
               es(i)=wn/8065.46d0
            endif
         enddo
         DO 5395 I=1,N
            DO 5394 J=1,N
               WLS(I,J)=0.0
               IF(I.EQ.J) GOTO 5394
               IF(ES(I).EQ.ES(J)) GOTO 5394
               WLS(I,J)=ABS(12398.54/(ES(I)-ES(J)))
 5394       CONTINUE
 5395    CONTINUE
         DO I=1,N
            DO  J=1,N
               OMS(I,J)=0.
               AS(I,J)=0.
            enddo
         enddo

         do i=1,1000
            read(23,*,end=88)iu,il,enl,enu,wlq,aij
            if(iu < n) then
               as(iu,il)=aij+as(iu,il)
            endif
         enddo
 88      continue
         close (23)
         open(23,file='./ATDAT/SiI_recomb_nahar_sorted.dat',status='old')
         read(23,*)(tesi1r(k),k=1,81)
         do i=1,nsi1
            read(23,*)il,(alsi1r(i,k),k=1,81)
         enddo
         read(23,987)dum
         read(23,*)(totrecsi1(k),k=1,81)
 987     format(a)
         Initsi1=1
      ENDIF
      DO I=1,N
         E(I)=ES(I)
         G(I)=GS(I)
         DO J=1,N
            A(I,J)=AS(I,J)
            OM(I,J)=OMS(I,J)
            WL(I,J)=WLS(I,J)
         ENDDO
      ENDDO

      SIG(1)=0.
      SIG(2)=0.
      DO 8457 I=3,N
 8457 SIG(I)=0.
      GA(1)=2.99
      GA(2)=2.5
      DO 1009 I=3,N
1009  GA(I)=3.
      GA(8)=3.56
      GA(9)=3.56
      DO 165 I=1,N
165   C(I,I)=0.
      RETURN
      END


      subroutine si_i_recomb(te)
      implicit real*8 (a-h,o-z)
c
      integer mie,mion

      parameter(mie=15,mion=5)
c
      data ieh/1/,iehe/2/,iec/3/,ien/4/,ieo/5/,iene/6/,iena/7/,
     &     iemg/8/,iesi/9/,ies/10/,iear/11/,ieca/12/,iefe/13/,
     &     ieco/14/,ieni/15/,iemax/15/
cqqqq
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      common/cnlev/nlev(mie,mion)
      COMMON/SILEV/NSIII,nsi1
      integer nlsi1
      parameter(nlsi1=66)
      common/csi_idat/gsi1(nlsi1)
      common/sirecombi/recosi1(nlsi1)
      common/si1_reco/tesi1r(81),alsi1r(65,81),totrecsi1(81)
c     
c     Si I
      do il=1,nsi1
         recosi1(il)=0.d0
      enddo
      if (te.lt.tesi1r(1)) then
         tconst=0.d0
         it1=1
         it2=1
      elseif (te.ge.tesi1r(81)) then
         tconst=0.d0
         it1=81
         it2=81
      else
         if(te.lt.tesi1r(1)) then
            it1=1
            it2=2
         elseif(te.gt.tesi1r(80)) then
            it1=79
            it2=80
         else
            do i=1,80
               if (te.ge.tesi1r(i).and.te.lt.tesi1r(i+1)) then
                  it1=i
                  it2=i+1
                  goto 9
               endif
            enddo
 9          continue
         endif
         tconst=(te-tesi1r(it1))/(tesi1r(it2)-tesi1r(it1))
      endif

      do il=1,nsi1
         recosi1(il)=recosi1(it1)+tconst*(recosi1(it2)-recosi1(it1))
      enddo
c     
      totrec=totrecsi1(it1)+tconst*(totrecsi1(it2)-totrecsi1(it1))
c     
c     Add the remaining recombinations to the 8 highest levels 
c     (114 - 121; y5P, y3F, y3D) weighted
c     with their statistical weights so that the total recombination 
c     coefficient is correct ; gsum=42 (7,5,7,3,5,7,5,3)
c     
      recsum=0.d0
      do il=1,n
         recsum=recsum+recosi1(il)
      enddo
      recadd=totrec-recsum
      if (recadd.lt.0.) then
         write(50,*) 'Error in calc. rec. coeff. for Fe I !!'
      endif
c     add the remaining recomb. contribution to
c     all the levels, not only the highest
      gsum=0.d0
      do il=1,nsi1
         gsum=gsum+gsi1(il)
      enddo
      do il=1,nsi1
         recosi1(il)=recosi1(il)+recadd*gsi1(il)/gsum
      enddo
c     
      return
      end



      SUBROUTINE CALCN(ION,ITER,R,TE,Z,XQ,DXOLD,AA,NRC)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
C
C     CALCULATE THE MATRIX FOR THE NEWTON RAPHSON ROUTINE
C     R = RADIUS
C     TE = TEMPERATURE
C     Z = ELECTRONIC FRACTION FROM IONS OTHER THAN HYDROGEN
C     DXOLD = OLD VALUES OF DEPARTURE COEFF.
C     AA = N-R MATRIX
C     NRC = NUMBER OF COLUMNS ( = NP1+1 = N+2 )
C
      COMMON/CINOUT/INOUT,IPULS
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(20,NE1:NE2),
     &SIGFE(NL,NE1:NE2),SIGH(12,NE1:NE2)
      COMMON/PHEAT/PHE(NL)
      COMMON/DIF/EMQW(MD,NE1:NE2),TAUC(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/CONT/CONR(NL),CONT(NL)
      COMMON/HRAT/RHYD,ZEL,XRQ,HYR,HEAT,COOL
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/A25/DHDT,DCODT,DZDT,DPDT
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/ABUN/AB(20)
      COMMON/LITER/NITER
      COMMON/IND/IKL
      COMMON/NBA/NBACK
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/SPH/ISPH
      COMMON/A14/C(NL,NL),COI(NL),G(NL),E(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/SECX/CSEC(20),DCSDX(20),CISEC(20),DCISDX
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),PH(NL),DRDT(NL),RECT(5,NL),PHET(NL)
      COMMON/A2/SI,T0,TAUL
      COMMON/BIN/B1IN,ELD
      COMMON/COLH/DXI,COLH(5,NL),COLHVT(5,NL),COLHVTT(5,NL),TAUBB
      COMMON/A5/TAU,ALFA,EN,TSN,XL40,TXEV,RMAX
      COMMON/A10/AV
      COMMON/A19/EM(NL,NL),ESC(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      COMMON/CSEX/CSLYA
      COMMON/LYBFLUOR/OIFLUOR,OIFLUORH,OIFLUORO,OILYB,OIREC
      common/test1/qa(5,20,20),rateph(5,20),tsav(5,20,20),esav(5,20,20)
      DIMENSION B(NL),DXOLD(NL),AA(NLP1,NLP1),CSECQ(NL)
      DATA CSAHA/2.0708E-16/
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      EC=ELCH
      SI=AV-1.
      NU=10
      NP1=N+1
C     TOTAL DENSITY
      DEN=DE(R)
      TEV=TE/1.1609E4
      TSEV=TSN/1.1609E4
      C1=TE**1.5/(CSAHA*DEN)
      C11=C1*EXP(-E00/TEV)
      DO I=1,N
      IF(I.LE.20) CSECQ(I)=CSEC(I)
      IF(I.GT.20) CSECQ(I)=0.
      ENDDO
      DO 1 I=1,NP1
      B(I)=DXOLD(I)
      DO 1 J=1,NRC
1     AA(I,J)=0.
C
C     ABXI = AB * FRACTION OF TOTAL IN THE LOWER AND UPPER STAGES
C            = 1 AND NP1
C     FOR NEUTRAL: X(I)+X(II)
C     FOR IONS:    X(II)+X(III)
C     FORMALLY ABXI=AB * (1. - SUM(X(I))) I=1,N, I NOT K, K+1 WHERE K
C           IS THE LOWER IONIZATION STAGE AND K+1 THE UPPER (B+)
C     ABUNDANCE IN UPPER STAGE IS THAN B+ * ABXI
C                  LOWER               SUM(XN) * ABXI
C
      IF(ION.EQ.1) ABXI=AB(11)*XQ
      IF(ION.EQ.2) ABXI=AB(3)*XQ
      IF(ION.EQ.4) ABXI=AB(8)*XQ
      XEL=ABXI*B(NP1)+Z
C     ONCE IONIZED IONS
      IF(ION.EQ.1) XEL=ABXI*(1.+B(NP1))+Z
      IF(ION.EQ.4) XEL=ABXI*(1.+B(NP1))+Z
      IF(ION.EQ.1) SECION=CISEC(5)
      IF(ION.EQ.2) SECION=CISEC(2)
      IF(ION.EQ.4) SECION=CISEC(17)
      C1=1./C1
C     TOTAL ELECTRON DENSITY
2735  CONTINUE
      DEN10=DEN/1.E10
      IF(ISTAT.EQ.1) THEN
C       FOR AGN'S
        VTERM=1.2855E4*SQRT(TE)
        IF(ION.EQ.1) VTERM=VTERM/6.3245
        IF(ION.EQ.2) VTERM=VTERM/4.
        IF(ION.EQ.3) VTERM=VTERM/2.
        IF(ION.EQ.4) VTERM=VTERM/7.483
C       CONST = 2.07E-16 * (1.E-8)**3 / (8 PI**1.5)
        C4=4.646D-42*ABXI*DXI*B(NP1)*XEL*DEN**2/(TE**1.5*VTERM)
      ELSE
C       FOR SN:E
        TIME=8.64E4*TDAY
C       CONST = 2.07E-16 * (1.E-8)**3 / 8 PI
        C4=8.236D-42*ABXI*TIME*B(NP1)*XEL*DEN**2/TE**1.5
      ENDIF
      F6=0.
      F11=0.
      DO 5 I=1,N
      FS=0.
      FRADU=0.
      FCOLLU=0.
      FRADL=0.
      FCOLLL=0.
      FSS=0.
      F5=0.
      F8=0.
      F7=0.
      F14=0.
      F9=0.
      DO 15 J=1,N
      FC=0.
      T0=0.
      IJ=MAX0(I,J)
      EIJ=ABS(E(I)-E(J))
      EDT=EIJ/TEV
      if(b(i).ne.0.) DBIJ=(B(I)-EXP(-EIJ/TEV)*B(J))/b(i)
      DBIJ=B(I)-EXP(-EIJ/TEV)*B(J)
      if(b(j).ne.0.) DBJI=(B(J)-EXP(-EIJ/TEV)*B(I))/b(j)
      DBJI=B(J)-EXP(-EIJ/TEV)*B(I)
      EQ1=ABS(EIJ/TEV)
      EQ2=ABS(E(I)-E00)/TEV
      IF(EIJ.LT.1.E-10) GOTO 16
C     CONTINUUM INTENSITY DIVIDED BY 2*H*V**3/C**2 = OCCUPATION NUMBER
      IF(EIJ.LE.0.) EIJ=1.E-10
      FC=FMEAN(3,EIJ)

      IF(NBACK.EQ.0) FC=0.
C     CALCULATE OPTICAL DEPTHS.
      IF(I.GT.J) GOTO 9481
C     I<J
      C5=C4*WL(J,I)**3*G(J)*EXP((E00-E(I))/TEV)
      T0=C4*WL(J,I)**3.*G(IJ)*A(J,I)*EXP((E00-E(I))/TEV)*B(I)
      IF(ISTAT.EQ.1) THEN
C     2.244E-26 = 1.E-24 / ( 8 PI**1.5 )
        T0=WL(J,I)**3.*G(IJ)*A(J,I)*(C4*EXP((E00-E(I))/TEV)*B(I)
     &      +2.244E-26*COLHVT(ION,I)/G(I))
        T0T=T0*(COLHVTT(ION,J)-COLHVT(ION,J))/COLHVT(ION,J)

        t0t=t0
        IF(I.NE.1) T0T=0.
        A21=A(J,I)
        WL21=WL(J,I)
      ENDIF
      DTADBP=(ABXI*B(NP1)+XEL)*C5*B(I)*A(J,I)/(B(NP1)*XEL)
      GOTO 16
C
C     I>J
C
9481  C5=C4*WL(I,J)**3*G(I)*EXP((E00-E(J))/TEV)
8268  T0=C4*WL(I,J)**3.*G(IJ)*A(I,J)*EXP((E00-E(J))/TEV)*B(J)
      IF(ISTAT.EQ.1) THEN
C     2.244E-26 = 1.E-24 / ( 8 PI**1.5 )
        T0=WL(I,J)**3.*G(IJ)*A(I,J)*(C4*EXP((E00-E(J))/TEV)*B(J)
     &      +2.244E-26*COLHVT(ION,J)/G(J))
        T0T=T0*(COLHVTT(ION,J)-COLHVT(ION,J))/COLHVT(ION,J)
        IF(J.NE.1) T0T=0.
c!!   
        T0T=T0
        A21=A(I,J)
        WL21=WL(I,J)
      ENDIF
      DTADBP=(ABXI*B(NP1)+XEL)*C5*B(J)*A(I,J)/(B(NP1)*XEL)
16    CONTINUE
C     CALCULATE ESCAPE PROBABILITY AND DBE/DTAU
      CALL ESCAPE(T0,T0T,WL21,A21,VTERM,BE,DBEDTA)
      ESC(J,I)=BE
      TOP(J,I)=T0
      IF(ION.NE.4) THEN
        ESAV(ION,J,I)=BE
        TSAV(ION,J,I)=T0
      ENDIF
      IF(I.GT.J) AIJN=A(I,J)/DEN
      IF(I.LE.J) AIJN=A(J,I)/DEN

      IF(I.LT.J) GOTO 6
C
C     I > J
C
C     A(I,J) CONTRIB. I NE J
C
      F1=XEL*C(I,J)
      DB=0.
      DFR=B(J)-B(I)*EXP(-EIJ/TEV)
      IF(ABS(DFR).GT.1.E-10) DB=DBEDTA*C5*A(I,J)
      IF(T0.lT.1.E-20) then
         DB=0.
      endif

      DBEDBP=DBEDTA*DTADBP
      F2=AIJN*(-B(I)*DB+(B(J)*EXP(EIJ/TEV)-B(I))*DB*FC+
     &     EXP(EIJ/TEV)*BE*FC)
      AA(I,J)=F1+F2
C
C     A(I,NRC) CONTRIB.
C
      F3=AIJN*(-B(I)*BE
     &+(B(J)*EXP(EIJ/TEV)-B(I))*BE*FC)
      F33=XEL*C(I,J)*(B(J)-B(I))
      FS=FS+F3+F33
      FCOLLL=FCOLLL+F33
      FRADL=FRADL+F3
      F14=F14+C(I,J)*(B(J)-B(I))*ABXI
     &+AIJN*(-B(I)+(B(J)*EXP(EIJ/TEV)-B(I))*FC)*DBEDBP
C
C     A(I,I) CONTRIB.
C
      DB=0.
      DFR=B(J)-B(I)*EXP(-EIJ/TEV)
      IF(T0.GT.1.) GOTO 7473
      DB=0.
 7473 CONTINUE
      F4=-XEL*C(I,J)+AIJN*(-BE-BE*FC-B(I)*DB+(B(J)*
     &EXP(EIJ/TEV)-B(I))*DB*FC)
      FSS=FSS+F4
      GOTO 15
6     CONTINUE
C     I < J
C
C     A(I,J) CONTRIB. I NE J
C
      DBREL=ABS(B(I)-B(J))/B(J)
      F1=XEL*C(I,J)
      DFR=B(I)-B(J)*EXP(-EIJ/TEV)
      DB=0.
      IF(T0.GT.1.E-20) GOTO 7474
      DB=0.
 7474 CONTINUE
      DBEDBP=DBEDTA*DTADBP
      F2=G(J)*EXP(-EIJ/TEV)*AIJN*(BE+B(J)*DB+(B(J)-B(I)*EXP(EIJ/TEV))
     &*DB*FC+BE*FC)/G(I)
      AA(I,J)=F1+F2
C
C     A(I,NRC) CONTRIB.
C
      F33=XEL*(B(J)-B(I))*C(I,J)
      F3=G(J)*EXP(-EIJ/
     &TEV)*AIJN*(B(J)*BE+(B(J)-B(I)*EXP(EIJ/TEV))*BE*FC)/G(I)
      FS=FS+F3+F33
      FCOLLU=FCOLLU+F33
      FRADU=FRADU+F3
      F7=C(I,J)*(B(J)-B(I))*ABXI+F7
     &+G(J)*EXP(-EIJ/TEV)*AIJN*(B(J)+(B(J)-B(I)*EXP(EIJ/
     &TEV))*FC)*DBEDBP/G(I)
C
C     A(I,I) CONTRIB.
C
      DB=0.
        IF(ABS(DFR).GT.1.E-10) DB=DBEDTA*C5*A(J,I)
      IF(T0.GT.1.E-20) GOTO 7475
      DB=0.
 7475 CONTINUE
      F4=-XEL*C(I,J)+G(J)*EXP(-EIJ/TEV)*AIJN*
     &(B(J)*DB+(B(J)-B(I)*EXP(EIJ/TEV))*DB*FC-EXP(EIJ/TEV)*
     &BE*FC)/G(I)
      FSS=FSS+F4
 15   CONTINUE
      REN=RTE(I)/DEN
      PHN=PH(I)/DEN
ccont
      BEC=1.
      DBECB1=0.
      DBECT=0.
      DBECN=0.
      AA(I,I)=FSS-PHN-XEL*COI(I)
      IF(I.EQ.1) AA(I,I)=AA(I,I)-(CSECQ(1)+SECION)/DEN**2
      AA(I,NRC)=-FS-BEC*REN-XEL*COI(I)+B(I)*(PHN+XEL*COI(I))
      IF(I.EQ.1) AA(I,NRC)=AA(I,NRC)+B(1)*(CSECQ(1)+SECION)/DEN**2
      IF(I.NE.1) AA(I,NRC)=AA(I,NRC)-B(1)*CSECQ(I)
C     O I LY B FLUORESENCE
      IF(ION.EQ.2) THEN
        IF(I.EQ.1) THEN
          AA(1,NRC)=AA(1,NRC)+B(1)*OIFLUORO/DEN
          AA(1,1)=AA(1,1)-OIFLUORO/DEN
        ELSEIF(I.EQ.8) THEN
          OILYB=B(1)*G(1)*EXP(E(8)/TEV)*OIFLUORO/(G(8)*DEN)
          AA(8,NRC)=AA(8,NRC)-OILYB
          OIREC=REN
          AA(8,1)=AA(8,1)+G(1)*EXP(E(8)/TEV)*OIFLUORO/(G(8)*DEN)
        ENDIF
      ENDIF
      AA(I,NP1)=F7+F14+COI(I)*ABXI*(1.-B(I))
      IF(I.EQ.1) AA(I,NP1)=AA(I,NP1)-B(1)*(DCSDX(1)+DCISDX)/DEN**2
      IF(I.NE.1) AA(I,NP1)=AA(I,NP1)+B(1)*DCSDX(1)/DEN**2
	if(ikl.lt.300) goto 5
	if(ikl.gt.30) stop
5     CONTINUE
C
C     A(I,NP1) CONTRIB.
C
      DO  I=1,N
         AA(NP1,I)=B(NP1)*XEL*G(I)*EXP(-E(I)/TEV)
      enddo
      SN=0.
      DO  I=1,N
         SN=SN+B(I)*G(I)*EXP(-E(I)/TEV)
      enddo
      AA(NP1,NRC)=-(SN*B(NP1)*XEL+C11*(B(NP1)-1.))
      AA(NP1,NP1)=(ABXI*B(NP1)+XEL)*SN+C11
      RETURN
      END

      DOUBLE PRECISION FUNCTION FX(X)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      COMMON/IND/I
      COMMON/DIF/EM(MD,NE1:NE2),TAUQ(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/NLEV/e000,NION,N,NP1H,NHMAX,NHMIN
      COMMON/A5/TAU,ALFA,EN,TSN,XL40,TXEV,R15
      COMMON/NION/ION
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/A21/J
      COMMON/A22/R,TE,MQ
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),ipara
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      FX=0.
      IF(TSN.LE.0.) WRITE(6,*)'TS',TSN
      TBEV=TSN/1.1609D4
      TEEV=TE/1.1609D4
      DEN=DE(R)
      E00=E000-E(J)
      ELOG0=LOG10(E00)
      E1=E00*10.**X+1.D-10
      EX1=999.D0
c      write(0,*)'ion,mq,e000,e00,x,e1 ',ion,mq,e000,e00,x,e1
C
C     FJ = MEAN INTENSITY OF RADIATION ( IN CGS UNITS )
C
      IF(QSO.LE.0.) EX1=E1/TBEV

      FJ=FMEAN(2,E1)
c      ENDIF
      IF(E1.LT.E00) GOTO 100
      EX3=E1/TEEV
      IF(ION.NE.2) GOTO 3887
C
C     CROSSECTIONS OF O I FROM DAWIES AND LEWIS
C
      IF(J.EQ.4) ERY=(E1-4.481)/13.61
      IF(J.EQ.5) ERY=(E1-4.110)/13.61
      IF(J.EQ.6) ERY=(E1-2.886)/13.61
      IF(J.EQ.7) ERY=(E1-2.642)/13.61
      IF(J.EQ.8) ERY=(E1-1.545)/13.61
      IF(J.EQ.9) ERY=(E1-1.095)/13.61
      IF(ERY.LE.0.) GOTO 100
      IF(J.EQ.4.AND.ERY.LT..4) SL=(-1.034*ERY+1.365-.0491/ERY)-2.
      IF(J.EQ.4.AND.ERY.GT..4) SL=-1.20-2.71*(ERY/.4)
      IF(J.EQ.5.AND.ERY.LT..5) SL=(-.4038*ERY+1.783-.0449/ERY)-2.
      IF(J.EQ.5.AND.ERY.GT..5) SL=-.483*ERY-.286
      S=10.**SL
      IF(J.EQ.6) S=10.**(-7.204*ERY+1.657-1.)+10.**(-.3748*ERY-.5877)
      IF(J.EQ.7) S=10.**(-5.509*ERY+.2457)+10.**(-.54*ERY-.1715)
      IF(J.EQ.8) S=20.5*(1.545/E1)**3.56
      IF(J.EQ.9) S=20.5*(1.545/E1)**3.56
      SIGM=S
 3887 CONTINUE
      GAUNT=1.
      if(ion.ne.5) GAUNT=1.
      if(ion.eq.5) GAUNT=gbf(j,e1)
      IF(ION.ne.2) SIGM=SIG(J)*GAUNT*(E00/E1)**GA(J)
      IF(ION.eq.2) SIGM=SIGOX(J,E1)
      IF(ION.EQ.1) SIGM=CRPHCA(J,E1)
C
C     IONIZATION RATE (MQ=5)
C     4*PI*LN(10.)*1.E-18/H = 4.36685E9
C
      IF(MQ.EQ.5) FX=4.36685E9*SIGM*FJ
      IF(MQ.EQ.5) GOTO 100
      IF(MQ.NE.2) GOTO 300
C
C     HEATING RATE (MQ=2)
C
C     FOR NET HEATING USE FIRST EXPR.
      FX=ELCH*4.36685E9*SIGM*(E1-E00)*FJ
300   IF(MQ.EQ.2) GOTO 100
C     TOTAL HEATING RATE
      IF(MQ.EQ.8) FX=ELCH*E1*4.36685D9*SIGM*FJ
      IF(MQ.EQ.8) GOTO 100
C
C     RECOMBINATION RATE
C
      IF(EX3.GT.700.) GOTO 100
      FX=4.36685D9*SIGM*(2.0845D-4*E1**3.+FJ)*EXP(-EX3)
      fx=1.637e8 
cCHECK!
      E00=E000-E(J)

      E1=E00*10.**X+1.D-10

      fx=1.637e8*2.306*G(J)*sigm*e1**3*EXP(-E00*(10**x-1.d0)/TEeV)/
     &     (2.*G(N+1)*TE**1.5)

      if(mq==11) then
         fx=sigm*e1**2*EXP(E00*(-10**x+1.d0)/TEeV)
      endif

      
C     DR/DT (MQ=4)
      IF(MQ.EQ.4) FX=FX*EX3/TE
C
C
C     NET RECOMBINATION COOLING (MQ=3)
C
      IF(MQ.EQ.3) FX=ELCH*(E1-E00)*FX
C     TOTAL REC COOLING RATE
      IF(MQ.EQ.7) FX=ELCH*E1*FX
C     DRECCOL/DT (MQ=6)
      IF(MQ.EQ.6) FX=ELCH*(E1-E00)*FX*EX3/TE
100   CONTINUE
      RETURN
      END

      SUBROUTINE HORATE(ION,R,TE)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
C
C     CALCULATE RECOMBINATION AND PHOTOIONIZATION RATES 
C
      COMMON/IND/I
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),PH(NL),DRDT(NL),RECT(5,NL),PHET(NL)
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/PHEAT/PHE(NL)
      COMMON/PHOTOHEAT/PHEAT(5,NL),PHEATT(5,NL)
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(20,NE1:NE2),
     &SIGFE(NL,NE1:NE2),SIGH(12,NE1:NE2)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      real*8 jmean
      COMMON/DTAU/JMEAN(NE1:NE2)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/HEIION/ZHEICO,ZHEILA,ZOTSHI,ZBALM,ZHBALM
      COMMON/RECAL/RECO(NL)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      TS=TE
      DO J=1,N
C     
C     FIRST CALCULATE THE RECOMBINATION RATE FOR ALL LEVELS
C     
         CALL RECOMB(ION,J,TS,AL,RI)
         RECO(J)=AL
C     
C     THEN THE IONIZATION AND HEATING RATE
C     
         E0=E00-EN(J)
         ZA=0.
         HEAT=0.
         HEATT=0.
         DO J1=JMIN,JJ
            IF(E1(J1).ge.E0) then
               IF(ION.EQ.5)  SIGM=SIGH(J,J1)
               IF(ION.EQ.4)  SIGM=SIGFE(J,J1)
               IF(ION.EQ.3)  SIGM=SIGHE(J,J1)
               IF(ION.EQ.2)  SIGM=SIGO(J,J1)
               IF(E1(J1).LE.0.) WRITE(6,*)' HOP ',J1,E1(J1)
c     nots      IF(ION.EQ.5) THEN
C     !    OTS - APPROX. FOR THE LYMAN CONTINUUM
c     nots      JMEAN(J1)=EXP(-TAUTOT(I-1,J1))*FMEAN(1,E1(J1))
c     nots      ELSE 

               ZC=4.*PI*JMEAN(J1)*SIGM*(E(J1+1)-E(J1))/(ELCH*E1(J1))
               ZA=ZA+ZC
               HEAT=ELCH*ZC*(E1(J1)-E0)+HEAT
               HEATT=ELCH*ZC*E1(J1)+HEATT
            endif
         enddo

         PHEAT(ION,J)=HEAT
         PHEATT(ION,J)=HEATT
         PH(J)=ZA
         RTE(J)=RI
      enddo
      IF(ION.EQ.5) ZBALM=PH(2)

      RETURN
      END

      SUBROUTINE HLOSS(ION,B,R,TE,XEL,XI,HC)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
C     **************************************************************
C     *****
C     THIS ROUTINE CALCULATES THE ENERGY LOSS DUE TO H-EMISSION
C     IN ERGS*CM**3 FOR BOUND-BOUND TRANSITIONS (SEE KWAN&KROLIK)
C     AND THE OBSERVED STRENTH OF THE LINE (SEE KLEIN AND CASTOR ).
c     EM12 = AB * N1(cm-3) * E12 * C12 * X(EL) =erg/cm3
c     IN LOW DENSITY LIMIT 
C     *****
C     **************************************************************
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/SPH/ISPH
      COMMON/NBA/NBACK
      COMMON/CICS/RSHOCK,ics
      COMMON/A19/EM(NL,NL),ESC(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A5/TAU,ALFA,EN,TSN,XL40,TXEV,R15
      COMMON/NLEV/e00,NION,N,NP1H,NMA,NMI
      COMMON/ABUN/AB(20)
      COMMON/EQUIH/FRH(2,500),WEQ(500),SR(NFEL),WOBS(NL,NL)
      DIMENSION B(NL)
      DATA CSAHA/2.0708E-16/
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
C
C     ABXI = AB * FRACTION OF TOTAL IN THE LOWER AND UPPER STAGES
C            = 1 AND NP1
C     FOR NEUTRAL: X(I)+X(II)
C     FOR IONS:    X(II)+X(III)
      IF(ION.EQ.1) ABXI=AB(11)*XI
      IF(ION.EQ.2) ABXI=AB(3)*XI
      IF(ION.EQ.3) ABXI=AB(2)*XI
      IF(ION.EQ.4) ABXI=AB(8)*XI
      IF(ION.EQ.5) ABXI=AB(1)
      TSEV=TSN/1.1609E4
      TEV=TE/1.1609E4
C
C     C1= E * SAHA FACTOR
C
      DEN=DE(R)
      NP1=N+1
      C1=CSAHA*ELCH*XEL*ABXI*B(NP1)*DEN/TE**1.5
      HC=0.
      K=0
      DO  I=2,N
         DO J=1,I-1
            EIJ=E(I)-E(J)
            IF(EIJ.LE.0.) THEN
               EIJ=1.E-10
               FIN=0.
            ELSE
               IF(NBACK.EQ.1) THEN
                  FIN=FMEAN(3,EIJ)
                  W=0.
               ELSE
                  FIN=0.
               ENDIF
            ENDIF
            if(ics.eq.1) fin=0.
            IF(QSO.le.0.) then
               EIJX=EIJ/TSEV
C     
C     DILUTION FACTOR
C     
               IF(ISPH.EQ.1) W=0.5*(1.-SQRT(1.-1./R**2.))
C     EVALUATE W AT SURFACE
               IF(ISPH.NE.1) W=.5
            endif
            IF(IAGN.EQ.1) W=1.E-30
C     
C     COLLISIONAL EXCITATION RATE FROM LEVEL J TO I
C     
            FA=G(I)*A(I,J)*EIJ*ESC(I,J)*EXP((E00-E(I))/TEV)
            FQ=B(I)-FIN*(B(J)*EXP((E(I)-E(J))/TEV)-B(I))
            FB=EIJ*G(J)*C(J,I)*EXP((E00-E(J))/TEV)*(B(J)-B(I))
C     
C     COOLING RATE IN ERG CM**3 / S
C     
            EM(I,J)=C1*FB
            HC=HC+EM(I,J)*XEL
            K=K+1
            FD=(1.-W)*B(I)-FIN*(B(J)*EXP((E(I)-E(J))/TEV)-B(I))
C     CONTRIBUTION TO LINE PROFILE
            SOURCE=0.
            IF(EIJ/TEV.LT.700.) SOURCE=2.031E-04*EIJ**3/
     &           (B(J)*EXP(EIJ/TEV)/B(I)-1.)
C     CONVERT TO INTENSITY/ ANGSTROM
            SOURCE=3.E18*SOURCE/WL(J,I)**2
            T0=TOP(J,I)
C     
C     OBSERVED FLUX AT INFINITY
C     
            IF(K.lT.500) then
               SR(K)=(1.-EXP(-T0))*SOURCE
               WEQ(K)=C1*FA*FD/DEN
               WOBS(I,J)=C1*FA*FD/DEN
               DNI=CSAHA*G(I)*B(I)*B(NP1)*XEL*DEN*EXP((E00-E(I))/TEV)/TE**1.5
            endif
         enddo
      enddo
100   CONTINUE
      RETURN
      END

      SUBROUTINE ATDATO
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameters.h'
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      NH=9
      N=NHMAX
      N=NH
      E00=13.618
      E(1)=0.
      E(2)=1.961
      E(3)=4.172
      E(4)=9.137
      E(5)=9.508
      E(6)=10.732
      E(7)=10.976
      E(8)=12.073
      E(9)=12.523
      G(1)=9.
      G(2)=5.
      G(3)=1.
      G(4)=5.
      G(5)=3.
      G(6)=15.
      G(7)=9.
      G(8)=15.
      G(9)=15.
      G(N+1)=4.
      DO 5395 I=1,N
      DO 5394 J=1,N
      WL(I,J)=0.0
      IF(I.EQ.J) GOTO 5394
      IF(E(I).EQ.E(J)) GOTO 5394
      WL(I,J)=ABS(12398.54/(E(I)-E(J)))
5394  CONTINUE
5395  CONTINUE
      DO 7352 I=1,N
      DO 7352 J=1,N
      C(I,J)=0.
 7352 A(I,J)=0.
c     Mendoza -83
      A(2,1)=8.45E-3
c     Mendoza -83
      A(3,1)=7.35E-2
C     Wiese et al SKELTO & SHINE LIST A(4,1)=5.88E3
      A(4,1)=1.70E3
c     Wiese and Martin, NBS -80
      A(5,1)=5.96E8
c     Christenson and Cunningham, J. Geophys. Res. 83, 4393 -78
      A(8,1)=.770E8
c     Christenson and Cunningham, J. Geophys. Res. 83, 4393 -78
      A(9,1)=2.28E8
c     Mendoza -83
      A(3,2)=1.22E0
      A(5,4)=0.
c     Wiese and Martin, NBS -80
      A(6,4)=0.34E8
c     erikson & toft
      A(7,4)=8.52E4
      A(6,5)=0.E0
c     Wiese, Smith  and Glennon, NBS -80
      A(7,5)=2.80E7
c     Christenson and Cunningham, J. Geophys. Res. 83, 4393 -78
      A(8,7)=.303E8
c     Christenson and Cunningham, J. Geophys. Res. 83, 4393 -78
      A(9,7)=9.59E4
7362  CONTINUE
C
      SIG(1)=7.91
      SIG(2)=15.
      DO 8457 I=3,N
8457  SIG(I)=7.91*REAL(I)
      SIG(8)=20.5
      SIG(9)=20.5
      GA(1)=2.99
      GA(2)=2.5
      DO 1009 I=3,N
1009  GA(I)=3.
      GA(8)=3.56
      GA(9)=3.56
      DO 165 I=1,N
165   C(I,I)=0.
      RETURN
      END

      SUBROUTINE CROXY
      IMPLICIT REAL*8(A-H,O-Z)  
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(20,NE1:NE2),
     &SIGFE(NL,NE1:NE2),SIGH(12,NE1:NE2)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      DO  J1=JMIN,JJ
         EN=E1(J1)
C     O I (MCALPINE CORRECTED TO AGREE WITH REILMAN AND MANSON AT
C     HIGH ENERGIES).
c     SIGO(1,J1)=SI(14,J1)
c     use only TopBase cross sections. Note 1-3 are from 1 (LS coupling)      
         DO J=1,13
            S=0.
            S=SIGOX(J,EN)
            SIGO(J,J1)=S*1.E-18
         enddo
      enddo
      END


      SUBROUTINE CRFEII
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(20,NE1:NE2),
     &SIGFE(NL,NE1:NE2),SIGH(12,NE1:NE2)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/SIK/SK(ncr,NE1:NE2)
      DO 100 J1=JMIN,JJ
      EN=E1(J1)
      SIGFE(1,J1)=0.0E0
      SIGFE(1,J1)=SI(43,J1)
      DO 200 J=2,12
      S=0.
 200  SIGFE(J,J1)=S*1.E-18
100   CONTINUE
      RETURN
      END

      SUBROUTINE HPCA(R,TE)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
C
C     CALCULATE RECOMBINATION AND PHOTOIONIZATION RATES FOR CALCIUM
C
      COMMON/IND/I
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),PH(NL),DRDT(NL),RECT(5,NL),PHET(NL)
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/PHEAT/PHE(NL)
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(20,NE1:NE2),
     &SIGFE(NL,NE1:NE2),SIGH(12,NE1:NE2)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      real*8 jmean
      COMMON/DTAU/JMEAN(NE1:NE2)
      COMMON/RECAL/RECO(NL)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      TS=TE
      DO J=1,N
C     
C     FIRST CALCULATE THE RECOMBINATION RATE FOR ALL LEVELS
C     
C     CALL CARECQ(J,TS,RI)
C     
C     THEN THE IONIZATION AND HEATING RATE
C     
         E0=E00-EN(J)
         ZA=0.
         HEAT=0.
         DO J1=JMIN,JJ
            IF(E1(J1).ge.E0) then
               SIGM=SIGCA(J,J1)
               ZC=4.*PI*JMEAN(J1)*SIGM*(E(J1+1)-E(J1))/(ELCH*E1(J1))
               ZA=ZA+ZC
               HEAT=ELCH*ZC*(E1(J1)-E0)+HEAT
            endif
         enddo
         PH(J)=ZA
         
c     recombination from OP corss section (see milne_ca_i_ii_v2 in REC_PHOTO)
c     add rec to levels > 3 ton=2 & 3 according to stat. weights.
         if(j==1) then
            reco(1)=3.6e-15*(1.e4/te)**0.7
         elseif(j==2) then
            reco(2)=(10*1.8e-13/16+2.5e-13)*(1.e4/te)**0.7
         elseif(j==3) then
            reco(3)=(6*1.8e-13/16+4.5e-14)*(1.e4/te)**0.7
         endif
                        
C     RTE(J)=RI
      enddo
      RETURN
      END

      SUBROUTINE CRCA
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(20,NE1:NE2),
     &SIGFE(NL,NE1:NE2),SIGH(12,NE1:NE2)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      DO 100 J1=JMIN,JJ
      DO 200 J=1,3
 200  SIGCA(J,J1)=1.E-18*CRPHCA(J,E1(J1))
      SIGCA(1,J1)=SI(51,J1)
100   CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION CRPHCA(J,E1)
C
C     CROSSECTIONS OF CA II FROM SHINE
C
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameters.h'
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      E0=E00-EN(J)
      S=0.
      IF(E1.LE.E0) GOTO 200
      IF(J.EQ.1) THEN
C     APPROX. FIT TO REILMAN AND MANSON
        IF(E1.LT.35.) S=0.20*(E0/E1)**0.0001
        IF(E1.GE.35..AND.E1.LT.45.) S=0.10*(35./E1)**1.
        IF(E1.GE.45..AND.E1.LT.160.) S=1.5
        IF(E1.GE.160..AND.E1.LT.390.) S=1.15*(160./E1)**1.66
        IF(E1.GE.390..AND.E1.LT.4000.) S=2.5*(390./E1)**2.18
        IF(E1.GT.4000.) S=6.85E-2*(4000./E1)**2.43
      ELSE
        IF(J.EQ.2) S=6.15*(E0/E1)**1.
        IF(J.EQ.3) S=2.38*(E0/E1)**3.65
      ENDIF
 200  CRPHCA=S
      RETURN
      END

      SUBROUTINE ATDAT
c Atomic data for Ca II      
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameters.h'
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      N=NHMAX
      E00=11.871
      E(1)=0.
      E(2)=1.6968
      E(3)=3.1415
      G(1)=2.
      G(2)=10.
      G(3)=6.
c     CA III 1S0
      G(N+1)=1.
      DO 5395 I=1,N
      DO 5394 J=1,N
      WL(I,J)=0.0
      IF(I.EQ.J) GOTO 5394
      IF(E(I).EQ.E(J)) GOTO 5394
      WL(I,J)=ABS(12398.54/(E(I)-E(J)))
5394  CONTINUE
5395  CONTINUE
      DO 7352 I=1,N
      DO 7352 J=1,N
      C(I,J)=0.
 7352 A(I,J)=0.
      A(2,1)=1.3
      A(3,1)=2.35E8
      A(3,2)=7.97E6
7362  CONTINUE
C
      SIG(1)=0.20
      SIG(2)=6.15
      SIG(3)=2.38
      GA(1)=0.
      GA(2)=1.0
      GA(3)=3.65
      DO 165 I=1,N
165   C(I,I)=0.
      RETURN
      END

      SUBROUTINE ATDATFE
      IMPLICIT REAL*8(A-H,O-Z)
      character*80 dum
      character*2 ch2
      character*4 ch4
      character*5 ch5
      SAVE
      include 'parameters.h'
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      COMMON/FELEV/NFEII,nfei
      common/initat/initfeii,initfei,initsi1
      DIMENSION GS(NL),ES(NL),AS(NL,NL),WLS(NL,NL),OMS(NL,NL)
      integer nlf1,nlf2,nlf3,nlf4
      parameter(nlf1=122,nlf2=192,nlf3=113,nlf4=46)
      common/cfedat/gf1(nlf1),gf2(nlf2),gf3(nlf3),gf4(nlf4)
      ntot=191
      nfeii=191
c!! OK      nfeii=130
      N=NFEII      
      E00=16.19791
      IF(INITFEii.EQ.0) THEN
         INITFEii=1
         OPEN(23,FILE='./ATDAT/feIInew2.dat')
         rewind 23
         do i=1,10
            read(23,987)dum
         enddo
 987     format(a)
         do i=1,ntot
            read(23,935)nn,wn,gsi,ch2,ch4,ig
            if(i<=nfeii) then
               gs(i)=gsi
               gf2(i)=gs(i)
 935           format(i3,f13.3,f6.1,1x,a2,1x,a4,i4)
               es(i)=wn/8065.46d0
            endif
         enddo
         DO 5395 I=1,N
            DO 5394 J=1,N
               WLS(I,J)=0.0
               IF(I.EQ.J) GOTO 5394
               IF(ES(I).EQ.ES(J)) GOTO 5394
               WLS(I,J)=ABS(12398.54/(ES(I)-ES(J)))
 5394       CONTINUE
 5395    CONTINUE
         DO I=1,N
            DO J=1,N
               OMS(I,J)=0.
               AS(I,J)=0.
            enddo
         enddo
c     read(23,987)dum
         read(23,987)dum
         read(23,987)dum
c     do i=2,nfeii
         do i=2,ntot
            do j=1,i-1
               read(23,937)i1,i2,ch5,g2,ch5,g1,wlq,asij,OMsq
               if(i<=nfeii) then
                  OMs(J,I)=OMsq
                  as(i,j)=asij
 937              format(2i4,2x,a5,f2.0,2x,a5,f2.0,2x,f12.2,1pe12.4,e12.4)
               endif
            enddo
         enddo
         close (23)
      ENDIF
      DO I=1,N
         E(I)=ES(I)
         G(I)=GS(I)
         DO J=1,N
            A(I,J)=AS(I,J)
            OM(I,J)=OMS(I,J)
            WL(I,J)=WLS(I,J)
         ENDDO
      ENDDO
C
      SIG(1)=0.
      SIG(2)=0.
      DO  I=3,N
         SIG(I)=0.
      enddo
      GA(1)=2.99
      GA(2)=2.5
      DO  I=3,N
         GA(I)=3.
      enddo
      GA(8)=3.56
      GA(9)=3.56
      DO I=1,N
         C(I,I)=0.
      enddo
      RETURN
      END

      SUBROUTINE COLLEX(ION,TS)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameters.h'
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      DIMENSION CQ(NL,NL),CQP(NL,NL)
      NP1=N+1
      DT=0.0001
      TSPDT=(1.+DT)*TS
      DT=DT*TS
      CALL CEX(ION,TS,CQ)
      CALL CEX(ION,TSPDT,CQP)
      DO J=1,N
         IP1=J+1
         DO I=IP1,N
            EIJ=ABS(E(J)-E(I))
            TEV=TSPDT/1.1609E4
            ET=EIJ/TEV
            IF(EIJ/TEV.GT.700.) ET=700.
            CPIJ=G(J)*EXP(ET)*CQP(J,I)/G(I)
            C(J,I)=CQ(J,I)
            TEV=TS/1.1609E4
            ET=EIJ/TEV
            IF(EIJ/TEV.GT.700.) ET=700.
            C(I,J)=G(J)*EXP(ET)*C(J,I)/G(I)
            DCDT(J,I)=(CQP(J,I)-C(J,I))/DT
            DCDT(I,J)=(CPIJ-C(I,J))/DT
         enddo
      enddo
C
C     CALCULATE COLL. IONIZATION RATES
C
      DO I=1,N
         NLQ=I
         NMAX=NP1
         CALL CION(ION,NLQ,NMAX,TS,CI(I))
         CALL CION(ION,NLQ,NMAX,TSPDT,CIP)
         DCIDT(I)=(CIP-CI(I))/DT
      enddo
      RETURN
      END

      SUBROUTINE CION(ION,N,N0,T,CI)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameters.h'
C
C     COLLISIONAL IONIZATION RATES
C
      COMMON/A14/CQ(NL,NL),CIQ(NL),G(NL),E(NL),AQ(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
       COMMON/NLEV/e00,NION,NHY,NP1H,NHMAX,NHMIN
      CI=0.
      IF(ION.NE.1) GOTO 3
C     COLLISIONAL IONIZATION RATES FOR CA II FROM SHINE
      TEV=T/1.1609E4
      EI=E00-E(N)
      EITEV=DMIN1(700.D0,EI/TEV)
      CI=1.14E-10*SQRT(T)*EXP(-EITEV)/EI**2
      GOTO 500
C
C     COLLISIONAL IONIZATION RATES FOR O I FROM SUMMERS
C
3     IF(ION.NE.2) GOTO 4
      T4=T/1.E4
      IF(T.LT.2.E3) GOTO 400
      IF(N.NE.1) GOTO 400
      CI=1.09E-10*SQRT(T)*EXP(-1.58E5/T)/(1.+0.1*T/1.58E5)
      GOTO 500
 400  IF(N.EQ.2) OM=1.68*T4**.956
      IF(N.EQ.3) OM=.496*T4**.952
      IF(N.EQ.4) OM=2.54*T4**.920
      IF(N.EQ.5) OM=1.80*T4**.915
      IF(N.EQ.6) OM=17.4*T4**.892
      IF(N.EQ.7) OM=12.3*T4**.885
      IF(N.EQ.8) OM=54.6*T4**.840
      IF(N.EQ.9) OM=101.*T4**.806
      EITEV=DMIN1(700.D0,(-E(N)+E00)*11609./T)
      CI=8.63E-6*OM*EXP(-EITEV)/(SQRT(T)*G(N))
4     IF(ION.NE.3) GOTO 5
C     HE I
c      CALL CIONHE(N,E(N),T,CI)
5     IF(ION.NE.4) GOTO 6
      CALL CIONFE(N,E(N),T,CI)
6     IF(ION.NE.5) GOTO 500
c      CALL CIONH(N,N0,T,CI)
500   CONTINUE
      RETURN
      END

      SUBROUTINE CEX(ION,T,C)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameters.h'
C
C     COLLISIONAL EXCITATION RATES. Note for Fe I omegas read from file (in atdatfei)
C
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/A14/CQ(NL,NL),CI(NL),G(NL),E(NL),AQ(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)

      COMMON/A27/OM(NL,NL)
      COMMON/FELEV/NFEII,nfei
      common/oicoll/omoi(5,13,13),teoi(5)
      DIMENSION C(NL,NL)
      T4=T/1.E4
      T3=T/1.E3
      TEV=T/1.1609E4
      IF(ION.eq.1) then
C     COLLISIONAL EXCITATION RATES FOR CA II FROM SHINE
         C(1,2)=14.0
         C(1,3)=14.3
         C(2,3)=19.8
      endif
      IF(ION.eq.2) then
c     O I      
         IF(T4.lt.0.5) then
            C(1,2)=.0151*T3**1.31
            C(1,3)=.00184*T3**1.32
            C(2,3)=.031*T3**.534
         elseif(t4.lt.1.) then
            C(1,2)=.266*T4**1.10
            C(1,3)=.034*T4**1.08
            C(2,3)=.105*T4**.52
         else
            C(1,2)=.266*T4**.91
            C(1,3)=.0324*T4**.91
            C(2,3)=.105*T4**.50
 50         C(1,4)=.173*T4**1.25
            C(1,5)=.995*T4**.80
            C(1,6)=.046*T4**.82
            C(1,7)=.041*T4**.50
            C(1,8)=.0017*T4**2.16
            C(1,9)=.0088*T4**2.21
            C(4,6)=18.47*TEV**.64
            C(4,7)=2.5E-2*TEV**.64
            C(5,7)=9.63*TEV**.64
            C(7,8)=20.7*TEV**.64
            C(7,9)=5.3E-2*TEV**.64
C     MEWE FORMULA
            C(1,4)=1.28
            C(1,5)=1.59
            C(1,6)=2.69
            C(1,7)=2.76
            C(1,8)=.51
            C(1,9)=1.34
            C(4,6)=148.*TEV**.64
            C(5,7)=94.*TEV**.64
            C(7,8)=447.*TEV**.64
            C(7,9)=93.*TEV**.64
         endif

c     include coll. exc. of O I from file OI_kb.dat
         do k=1,4
            if(t > teoi(k).and.t<=teoi(k+1)) then
               it=k
            endif
         enddo
         do i=1,12
            do j=i,13
               if(t<teoi(1)) then
                  c(i,j) = omoi(1,i,j)
               elseif(t>teoi(5)) then
                  c(i,j) = omoi(5,i,j)
               else               
                  c(i,j)=omoi(it,i,j) + (t-teoi(it))*
     &                 (omoi(it+1,i,j)-omoi(it,i,j))/(teoi(it+1)-teoi(it))
               endif
            enddo
         enddo
      endif
      
      IF(ION.eq.3) then
C     HE I 
         DO I=1,16
            DO J=I+1,16
               C(I,J)=OM(I,J)
            enddo
         enddo
      endif
      IF(ION.eq.6) then
C     FE I         
         DO  I=1,NFEI
            DO  J=I+1,NFEI
               C(I,J)=OM(I,J)
            enddo
         enddo
      endif
      IF(ION.eq.4) then
C     FE II
         DO  I=1,NFEII
            DO  J=I+1,NFEII
               C(I,J)=OM(I,J)
            enddo
         enddo
      endif
      DO I=1,N
         IP1=I+1
         DO J=IP1,N
            EIJ=ABS(E(J)-E(I))
            TEV=T/1.1609E4
            ET=EIJ/TEV
            IF(EIJ/TEV.GT.700.) ET=700.
            C(I,J)=8.63E-6*C(I,J)*EXP(-ET)/(G(I)*SQRT(T))
         enddo
      enddo
15    CONTINUE
      RETURN
      END


      SUBROUTINE CIONFE(N,EN,T,CI)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     COLLISIONAL IONIZATION RATES FOR FE II
C
      CON=5.465E-11
      Q=0.
      T4=T/1.E4
      YN=1.1609*EION/T4
      CI=CON*SQRT(T)*Q*A1
      RETURN
      END


      SUBROUTINE OXREC(OXR,TE)
C     ************************************************************
C     ******
C     O I RECOMBINATION EMISSION CALCULATED FROM JULIENNE ET AL
C     AND SCALED TO 1.E4 K BY A POWER T4**-0.7
C     WAVELENGTHS IN ORDER 1304, 8446, 11287, 7002, 1356, 7774,
C                          9264, 6157
C     ******
C     *************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION OXR(10)
      T4=TE/1.E4
      OXR(1)=9.44E-25/T4**.7
      OXR(2)=.143
      OXR(3)=.070
      OXR(4)=.0133
      OXR(5)=1.68
      OXR(6)=0.288
      OXR(7)=0.135
      OXR(8)=.033
      DO 1 I=2,8
1     OXR(I)=OXR(I)*OXR(1)
      RETURN
      END

      SUBROUTINE RECOMB(ION,I,TE,A,RI)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameters.h'
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/HEIRECION/AL0(16),BE(16),B0(16),GA(16)
      DATA CSAHA/2.0708E-16/
      A=0.
      TEV=TE/1.1609E4
      TEV0=TEV*10.
      T4=TE/1.E4
      IF(ION.NE.2) GOTO 1
C     O I
c$$$      IF(I.EQ.1) A=1.36E-13/TEV**.7
c$$$      IF(I.EQ.4) A=1.6E-14/TEV0**.7
c$$$      IF(I.EQ.5) A=.14E-13/TEV0**.7
c$$$      IF(I.EQ.6) A=4.77E-13/TEV0**.7
c$$$      IF(I.EQ.7) A=0.47E-13/TEV0**.7
c$$$      IF(I.EQ.8) A=1.68E-13/TEV0**.7
c$$$      IF(I.EQ.9) A=0.52E-13/TEV0**.7

c    rec. rates from Julienne, P.S., Davies, J. & Oran, E. 1974,                                                           
c J. Geophys. Res. 79, 2540                                                                                             
c     3P divided according to statistical weights. Total 4.4e-13
c SHOULD BBE REPLACED BY NAHAR!      
      IF(I.EQ.1) A=2.44
      IF(I.EQ.2) A=1.466
      IF(I.EQ.3) A=0.4888
      IF(I.EQ.4) A=0.
      IF(I.EQ.5) A=0.
      IF(I.EQ.6) A=0.10
      IF(I.EQ.7) A=0.04
      IF(I.EQ.8) A=2.1
      IF(I.EQ.9) A=0.11
      IF(I.EQ.10) A=0.
      IF(I.EQ.11) A=0.
      IF(I.EQ.12) A=2.7
      IF(I.EQ.13) A=1.2

c correction to give same total rate as Chung, Lin & Lee                                                                
c + same temp. dep.                                                                                                     
c     temp in 0.1 eV                                                                                                    
      t01ev=te/1164.
      a=(1.47/1.065)*a*1.e-13/t01ev**0.8433



      
      GOTO 2
1     IF(ION.NE.3) GOTO 2
C     HE I
      T5000=TE/5.E3
      IF(I.EQ.1) A=2.23E-13/T5000**.488
      IF(I.EQ.2) A=1.98E-14/T5000**.42
      IF(I.EQ.3) A=3.6E-14/T5000**.784
      IF(I.EQ.4) A=1.135E-13/T5000**.610
      IF(I.EQ.5) A=7.2E-14/T5000**.784
      IF(I.EQ.6) A=8.45E-15/T5000**.45
      IF(I.EQ.7) A=3.92E-14/T5000**.63
      IF(I.EQ.8) A=4.894E-14/T5000**.85
      IF(I.EQ.9) A=0.343E-14/T5000**.478
      IF(I.EQ.10) A=1.77E-14/T5000**.620
      IF(I.EQ.11) A=3.14E-14/T5000**.851
      IF(I.EQ.12) A=4.87e-14/T5000**1.173
C     FROM ALMOG & NETZER
      A=AL0(I)*T4**BE(I)
2     IF(ION.NE.4) GOTO 3
C     FE II ???
      T5000=TE/5.E3
      IF(I.EQ.1) A=2.23E-13/T5000**.488
      IF(I.EQ.2) A=1.98E-14/T5000**.44
      IF(I.EQ.3) A=0.764E-14/T5000**.461
      IF(I.EQ.4) A=1.34E-13/T5000**.7
      IF(I.EQ.5) A=8.34E-14/T5000**.7
      IF(I.EQ.6) A=3.42E-15/T5000**.7
      IF(I.EQ.7) A=4.13E-14/T5000**.7
      IF(I.EQ.8) A=8.38E-14/T5000**.7
      IF(I.EQ.9) A=0.31E-14/T5000**.7
      IF(I.EQ.10) A=1.77E-14/T5000**.7
      IF(I.EQ.11) A=3.45E-14/T5000**.7
      IF(I.EQ.12) A=0.817e-14/T5000**.7
3     RI=A*EXP((E(I)-E00)/TEV)*TE**1.5/(CSAHA*G(I))
      RETURN
      END


      subroutine readferec
      implicit real*8 (a-h,o-z)
      common/cferec_f1/tef1(81),rectf1(81)
      common/cferecl_f1/rcmpf1(13,81),rcsumf1(10,81)
      common/cferec_f2/tef2(81),rectf2(81),rcmpf2(42,81),rcsumf2(17,81)
      dimension tv(81)
c
c Read total Fe I recombination coefficients
      open(22,file='./ATDAT/rec_tot_fe1.dat',status='old')
      do i=1,4
         read(22,*) 
      enddo
      do i=1,81
         read(22,*) tef1(i),rectf1(i)
      enddo
      close(22)
c Read total Fe II recombination coefficients
      open(22,file='./ATDAT/rec_tot_fe2.dat',status='old')
      do i=1,4
         read(22,*) 
      enddo
      do i=1,81
         read(22,*) tef2(i),rectf2(i)
      enddo
      close(22)
c Read Fe I recombination coefficients to specific multiplets
      open(22,file='./ATDAT/rec_lev_fe1.dat',status='old')
      read(22,*)
      read(22,701) (tv(i),i=1,81)
c
      do imp=1,13
         read(22,*)
         read(22,701) (rcmpf1(imp,it),it=1,81)
      enddo
c Read Fe I sum of recombination coefficients to higher multiplets
      read(22,*)
      do imps=1,10
         read(22,702) (rcsumf1(imps,it),it=1,81)
      enddo
 701  format(19x,85e10.2)
 702  format(11x,85e10.2)
      close(22)
c Read Fe II recombination coefficients to specific multiplets
      open(22,file='./ATDAT/rec_lev_fe2.dat',status='old')
      read(22,*)
      read(22,701) (tv(i),i=1,81)
      do imp=1,42
         read(22,*)
         read(22,701) (rcmpf2(imp,it),it=1,81)
      enddo
c Read Fe II sum of recombination coefficients to higher multiplets
      read(22,*)
      do imps=1,17
         read(22,702) (rcsumf2(imps,it),it=1,81)
      enddo
      close(22)
c
      return
      end

      subroutine ferecomb(ife,te)
      implicit real*8 (a-h,o-z)
c
      integer mie,mion
      integer nlf1,nlf2,nlf3,nlf4
      parameter(nlf1=122,nlf2=192,nlf3=113,nlf4=46)
      parameter(mie=15,mion=5)
c
      data ieh/1/,iehe/2/,iec/3/,ien/4/,ieo/5/,iene/6/,iena/7/,
     &     iemg/8/,iesi/9/,ies/10/,iear/11/,ieca/12/,iefe/13/,
     &     ieco/14/,ieni/15/,iemax/15/
c
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      common/cnlev/nlev(mie,mion)
      common/cferec_f1/tef1(81),rectf1(81)
      common/cferecl_f1/rcmpf1(13,81),rcsumf1(10,81)
      common/cferec_f2/tef2(81),rectf2(81),rcmpf2(42,81),rcsumf2(17,81)
      common/cfedat/gf1(nlf1),gf2(nlf2),gf3(nlf3),gf4(nlf4)
      common/ferecombi/recof1(nlf1),recof2(nlf2)
      if(ife==1) then
c     
c     Fe I
         do il=1,nlf1
            recof1(il)=0.d0
         enddo
         if (te.lt.tef1(1)) then
            tconst=0.d0
            it1=1
            it2=1
         elseif (te.ge.tef1(81)) then
            tconst=0.d0
            it1=81
            it2=81
         else
            do i=1,80
               if (te.ge.tef1(i).and.te.lt.tef1(i+1)) then
                  it1=i
                  it2=i+1
                  goto 9
               endif
            enddo
 9          continue
            tconst=(te-tef1(it1))/(tef1(it2)-tef1(it1))
         endif
c     
         rctfe1=rectf1(it1)+tconst*(rectf1(it2)-rectf1(it1))
c     
c     a 5De ; imp=1, il=1(g=9),2(7),3(5),4(3),5(1); gsum=25
c     Add 5De-sum; imps=2
         rcmp=rcmpf1(1,it1)+tconst*(rcmpf1(1,it2)-rcmpf1(1,it1))
         rcsum=rcsumf1(2,it1)+tconst*(rcsumf1(2,it2)-rcsumf1(2,it1))
         rcmp=rcmp+rcsum
         recof1(1)=rcmp*0.36
         recof1(2)=rcmp*0.28
         recof1(3)=rcmp*0.20
         recof1(4)=rcmp*0.12
         recof1(5)=rcmp*0.04
c     a 5Fe ; imp=2, il=6(11),7(9),8(7),9(5),10(3); gsum=35
c     Add 5Fe-sum; imps=3
         rcmp=rcmpf1(2,it1)+tconst*(rcmpf1(2,it2)-rcmpf1(2,it1))
         rcsum=rcsumf1(3,it1)+tconst*(rcsumf1(3,it2)-rcsumf1(3,it1))
         rcmp=rcmp+rcsum
         recof1(6)=rcmp*0.3143
         recof1(7)=rcmp*0.2571
         recof1(8)=rcmp*0.2
         recof1(9)=rcmp*0.1429
         recof1(10)=rcmp*0.0857
c     a 5Pe ; imp=3, il=14(7),15(5),16(3); gsum=15
c     Add 5Pe-sum; imps=1
         rcmp=rcmpf1(3,it1)+tconst*(rcmpf1(3,it2)-rcmpf1(3,it1))
         rcsum=rcsumf1(1,it1)+tconst*(rcsumf1(1,it2)-rcsumf1(1,it1))
         rcmp=rcmp+rcsum
         recof1(14)=rcmp*0.4667
         recof1(15)=rcmp*0.3333
         recof1(16)=rcmp*0.2
c     z 7Do ; imp=4, il=18(11),21(9),23(7),25(5),26(3); gsum=35
c     Add 7Do-sum; imps=9
         rcmp=rcmpf1(4,it1)+tconst*(rcmpf1(4,it2)-rcmpf1(4,it1))
         rcsum=rcsumf1(9,it1)+tconst*(rcsumf1(9,it2)-rcsumf1(9,it1))
         rcmp=rcmp+rcsum
         recof1(18)=rcmp*0.3143
         recof1(21)=rcmp*0.2571
         recof1(23)=rcmp*0.2
         recof1(25)=rcmp*0.1429
         recof1(26)=rcmp*0.0857
c     z 7Fo ; imp=5, il=34(13),36(11),38(9),40(7),41(5),42(3),43(1); gsum=49
c     Add 7Fo-sum; imps=10
         rcmp=rcmpf1(5,it1)+tconst*(rcmpf1(5,it2)-rcmpf1(5,it1))
         rcsum=rcsumf1(10,it1)+tconst*(rcsumf1(10,it2)-rcsumf1(10,it1))
         rcmp=rcmp+rcsum
         recof1(34)=rcmp*0.2653
         recof1(36)=rcmp*0.2245
         recof1(38)=rcmp*0.1837
         recof1(40)=rcmp*0.1429
         recof1(41)=rcmp*0.1020
         recof1(42)=rcmp*0.0612
         recof1(43)=rcmp*0.0204
c     z 7Po ; imp=6, il=44(9),47(7),50(5); gsum=21
c     Add 7Po-sum; imps=8
         rcmp=rcmpf1(6,it1)+tconst*(rcmpf1(6,it2)-rcmpf1(6,it1))
         rcsum=rcsumf1(8,it1)+tconst*(rcsumf1(8,it2)-rcsumf1(8,it1))
         rcmp=rcmp+rcsum
         recof1(44)=rcmp*0.4286
         recof1(47)=rcmp*0.3333
         recof1(50)=rcmp*0.2381
c     z 5Do ; imp=7, il=54(9),56(7),58(5),61(3),62(1); gsum=25
         rcmp=rcmpf1(7,it1)+tconst*(rcmpf1(7,it2)-rcmpf1(7,it1))
         recof1(54)=rcmp*0.36
         recof1(56)=rcmp*0.28
         recof1(58)=rcmp*0.20
         recof1(61)=rcmp*0.12
         recof1(62)=rcmp*0.04
c     z 5Fo ; imp=8, il=65(11),66(9),67(7),69(5),70(3); gsum=35
         rcmp=rcmpf1(8,it1)+tconst*(rcmpf1(8,it2)-rcmpf1(8,it1))
         recof1(65)=rcmp*0.3143
         recof1(66)=rcmp*0.2571
         recof1(67)=rcmp*0.2
         recof1(69)=rcmp*0.1429
         recof1(70)=rcmp*0.0857
c     z 5Po ; imp=9, il=73(7),78(5),79(3); gsum=15
         rcmp=rcmpf1(9,it1)+tconst*(rcmpf1(9,it2)-rcmpf1(9,it1))
         recof1(73)=rcmp*0.4667
         recof1(78)=rcmp*0.3333
         recof1(79)=rcmp*0.2
c     y 5Do ; imp=10, il=89(9),91(7),94(5),96(3),98(1); gsum=25
c     Add 5Do-sum; imps=5
         rcmp=rcmpf1(10,it1)+tconst*(rcmpf1(10,it2)-rcmpf1(10,it1))
         rcsum=rcsumf1(5,it1)+tconst*(rcsumf1(5,it2)-rcsumf1(5,it1))
         rcmp=rcmp+rcsum
         recof1(89)=rcmp*0.36
         recof1(91)=rcmp*0.28
         recof1(94)=rcmp*0.20
         recof1(96)=rcmp*0.12
         recof1(98)=rcmp*0.04
c     y 5Fo ; imp=11, il=92(11),97(9),99(7),101(5),104(3); gsum=35
c     Add 5Fo-sum; imps=6
         rcmp=rcmpf1(11,it1)+tconst*(rcmpf1(11,it2)-rcmpf1(11,it1))
         rcsum=rcsumf1(6,it1)+tconst*(rcsumf1(6,it2)-rcsumf1(6,it1))
         rcmp=rcmp+rcsum
         recof1(92)=rcmp*0.3143
         recof1(97)=rcmp*0.2571
         recof1(99)=rcmp*0.2
         recof1(101)=rcmp*0.1429
         recof1(104)=rcmp*0.0857
c     z 5Go ; imp=12, il=105(11),106(13),107(9),109(7),111(5); gsum=45
c     Add 5Go-sum; imps=7
         rcmp=rcmpf1(12,it1)+tconst*(rcmpf1(12,it2)-rcmpf1(12,it1))
         rcsum=rcsumf1(7,it1)+tconst*(rcsumf1(7,it2)-rcsumf1(7,it1))
         rcmp=rcmp+rcsum
         recof1(105)=rcmp*0.2444
         recof1(106)=rcmp*0.2889
         recof1(107)=rcmp*0.2
         recof1(109)=rcmp*0.1556
         recof1(111)=rcmp*0.1111
c     y 5Po ; imp=13, il=114(7),115(5),117(3); gsum=15
c     Add 5Po-sum; imps=4
         rcmp=rcmpf1(13,it1)+tconst*(rcmpf1(13,it2)-rcmpf1(13,it1))
         rcsum=rcsumf1(4,it1)+tconst*(rcsumf1(4,it2)-rcsumf1(4,it1))
         rcmp=rcmp+rcsum
         recof1(114)=rcmp*0.4667
         recof1(115)=rcmp*0.3333
         recof1(117)=rcmp*0.2
c     
c     Add the remaining recombinations to the 8 highest levels 
c     (114 - 121; y5P, y3F, y3D) weighted
c     with their statistical weights so that the total recombination 
c     coefficient is correct ; gsum=42 (7,5,7,3,5,7,5,3)
c     
         recsum=0.d0
         do il=1,n
            recsum=recsum+recof1(il)
         enddo
         recadd=rctfe1-recsum
         if (recadd.lt.0.) then
            write(50,*) 'Error in calc. rec. coeff. for Fe I !!'
         endif
c     add the remaining recomb. contribution to
c     all the levels, not only the highest
         gsum=0.d0
         do il=1,n
            gsum=gsum+gf1(il)
         enddo
         do il=1,n
            recof1(il)=recof1(il)+recadd*gf1(il)/gsum
         enddo

      elseif(ife==2) then
c     
c     Fe II
c     
         do il=1,nlf2
            recof2(il)=0.d0
         enddo
         if (te.lt.tef2(1)) then
            tconst=0.d0
            it1=1
            it2=1
         elseif (te.ge.tef2(81)) then
            tconst=0.d0
            it1=81
            it2=81
         else
            do i=1,80
               if (te.ge.tef2(i).and.te.lt.tef2(i+1)) then
                  it1=i
                  it2=i+1
                  goto 10
               endif
            enddo
 10         continue
            tconst=(te-tef2(it1))/(tef2(it2)-tef2(it1))
         endif
c     
         rctfe2=rectf2(it1)+tconst*(rectf2(it2)-rectf2(it1))
c     
c     a 6De ; imp=1, il=1(10),2(8),3(6),4(4),5(2); gsum=30
c     Add 6De-sum; imps=3
         rcmp=rcmpf2(1,it1)+tconst*(rcmpf2(1,it2)-rcmpf2(1,it1))
         rcsum=rcsumf2(3,it1)+tconst*(rcsumf2(3,it2)-rcsumf2(3,it1))
         rcmp=rcmp+rcsum
         recof2(1)=rcmp*0.3333
         recof2(2)=rcmp*0.2667
         recof2(3)=rcmp*0.2
         recof2(4)=rcmp*0.1333
         recof2(5)=rcmp*0.0667
c     a 4Fe ; imp=2, il=6(10),7(8),8(6),9(4); gsum=28
         rcmp=rcmpf2(2,it1)+tconst*(rcmpf2(2,it2)-rcmpf2(2,it1))
         recof2(6)=rcmp*0.3571
         recof2(7)=rcmp*0.2857
         recof2(8)=rcmp*0.2143
         recof2(9)=rcmp*0.1429
c     a 4De ; imp=3, il=10(8),11(6),12(4),13(2); gsum=20
         rcmp=rcmpf2(3,it1)+tconst*(rcmpf2(3,it2)-rcmpf2(3,it1))
         recof2(10)=rcmp*0.4
         recof2(11)=rcmp*0.3
         recof2(12)=rcmp*0.2
         recof2(13)=rcmp*0.1
c     a 4Pe ; imp=4, il=14(6),15(4),16(2); gsum=12
         rcmp=rcmpf2(4,it1)+tconst*(rcmpf2(4,it2)-rcmpf2(4,it1))
         recof2(14)=rcmp*0.5
         recof2(15)=rcmp*0.3333
         recof2(16)=rcmp*0.1667
c     b 4Pe ; imp=5, il=24(6),30(4),31(2); gsum=12
         rcmp=rcmpf2(5,it1)+tconst*(rcmpf2(5,it2)-rcmpf2(5,it1))
         recof2(24)=rcmp*0.5
         recof2(30)=rcmp*0.3333
         recof2(31)=rcmp*0.1667
c     a 4He ; imp=6, il=25(14),27(12),28(10),29(8); gsum=44
c     Add 4He-sum; imps=15
         rcmp=rcmpf2(6,it1)+tconst*(rcmpf2(6,it2)-rcmpf2(6,it1))
         rcsum=rcsumf2(15,it1)+tconst*(rcsumf2(15,it2)-rcsumf2(15,it1))
         rcmp=rcmp+rcsum
         recof2(25)=rcmp*0.3182
         recof2(27)=rcmp*0.2727
         recof2(28)=rcmp*0.2273
         recof2(29)=rcmp*0.1818
c     b 4Fe ; imp=7, il=32(10),33(8),34(6),35(4); gsum=28
         rcmp=rcmpf2(7,it1)+tconst*(rcmpf2(7,it2)-rcmpf2(7,it1))
         recof2(32)=rcmp*0.3571
         recof2(33)=rcmp*0.2857
         recof2(34)=rcmp*0.2143
         recof2(35)=rcmp*0.1429
c     a 6Se ; imp=8, il=36(6) 
c     Add 6Se-sum; imps=1
         rcmp=rcmpf2(8,it1)+tconst*(rcmpf2(8,it2)-rcmpf2(8,it1))
         rcsum=rcsumf2(1,it1)+tconst*(rcsumf2(1,it2)-rcsumf2(1,it1))
         rcmp=rcmp+rcsum
         recof2(36)=rcmp
c     a 4Ge ; imp=9, il=37(12),39(10),40(8),41(6); gsum=36
         rcmp=rcmpf2(9,it1)+tconst*(rcmpf2(9,it2)-rcmpf2(9,it1))
         recof2(37)=rcmp*0.3333
         recof2(39)=rcmp*0.2778
         recof2(40)=rcmp*0.2222
         recof2(41)=rcmp*0.1667
c     b 4De ; imp=10, il=49(4),50(2),51(6),52(8); gsum=20
         rcmp=rcmpf2(10,it1)+tconst*(rcmpf2(10,it2)-rcmpf2(10,it1))
         recof2(49)=rcmp*0.2
         recof2(50)=rcmp*0.1
         recof2(51)=rcmp*0.3
         recof2(52)=rcmp*0.4
c     z 6Do ; imp=11, il=64(10),65(8),66(6),67(4),68(2); gsum=30
c     Add 6Do-sum; imps=4
         rcmp=rcmpf2(11,it1)+tconst*(rcmpf2(11,it2)-rcmpf2(11,it1))
         rcsum=rcsumf2(4,it1)+tconst*(rcsumf2(4,it2)-rcsumf2(4,it1))
         rcmp=rcmp+rcsum
         recof2(64)=rcmp*0.3333
         recof2(65)=rcmp*0.2667
         recof2(66)=rcmp*0.2
         recof2(67)=rcmp*0.1333
         recof2(68)=rcmp*0.0667
c     z 6Fo ; imp=12, il=69(12),70(10),71(8),72(6),73(4),74(2); gsum=42
c     Add 6Fo-sum; imps=5
         rcmp=rcmpf2(12,it1)+tconst*(rcmpf2(12,it2)-rcmpf2(12,it1))
         rcsum=rcsumf2(5,it1)+tconst*(rcsumf2(5,it2)-rcsumf2(5,it1))
         rcmp=rcmp+rcsum
         recof2(69)=rcmp*0.2857
         recof2(70)=rcmp*0.2381
         recof2(71)=rcmp*0.1905
         recof2(72)=rcmp*0.1429
         recof2(73)=rcmp*0.0952
         recof2(74)=rcmp*0.0476
c     z 6Po ; imp=13, il=75(8),76(6),77(4); gsum=18
c     Add 6Po-sum; imps=2
         rcmp=rcmpf2(13,it1)+tconst*(rcmpf2(13,it2)-rcmpf2(13,it1))
         rcsum=rcsumf2(2,it1)+tconst*(rcsumf2(2,it2)-rcsumf2(2,it1))
         rcmp=rcmp+rcsum
         recof2(75)=rcmp*0.4445
         recof2(76)=rcmp*0.3333
         recof2(77)=rcmp*0.2222
c     z 4Fo ; imp=14, il=78(10),80(8),85(6),87(4); gsum=28
         rcmp=rcmpf2(14,it1)+tconst*(rcmpf2(14,it2)-rcmpf2(14,it1))
         recof2(78)=rcmp*0.3571
         recof2(80)=rcmp*0.2857
         recof2(85)=rcmp*0.2143
         recof2(87)=rcmp*0.1429
c     z 4Do ; imp=15, il=79(8),81(6),84(4),86(2); gsum=20
         rcmp=rcmpf2(15,it1)+tconst*(rcmpf2(15,it2)-rcmpf2(15,it1))
         recof2(79)=rcmp*0.4
         recof2(81)=rcmp*0.3
         recof2(84)=rcmp*0.2
         recof2(86)=rcmp*0.1
c     z 4Po ; imp=16, il=88(6),89(4),90(2); gsum=12
         rcmp=rcmpf2(16,it1)+tconst*(rcmpf2(16,it2)-rcmpf2(16,it1))
         recof2(88)=rcmp*0.5
         recof2(89)=rcmp*0.3333
         recof2(90)=rcmp*0.1667
c     c 4Pe ; imp=17, il=93(2),94(4),99(6); gsum=12
         rcmp=rcmpf2(17,it1)+tconst*(rcmpf2(17,it2)-rcmpf2(17,it1))
         recof2(93)=rcmp*0.1667
         recof2(94)=rcmp*0.3333
         recof2(99)=rcmp*0.5
c     c 4Fe ; imp=18, il=95(4),96(6),97(10),98(8); gsum=28
c     Add 4Fe-sum; imps=11
         rcmp=rcmpf2(18,it1)+tconst*(rcmpf2(18,it2)-rcmpf2(18,it1))
         rcsum=rcsumf2(11,it1)+tconst*(rcsumf2(11,it2)-rcsumf2(11,it1))
         rcmp=rcmp+rcsum
         recof2(95)=rcmp*0.1429
         recof2(96)=rcmp*0.2143
         recof2(97)=rcmp*0.3571
         recof2(98)=rcmp*0.2857
c     b 4Ge ; imp=19, il=101(12),102(10),103(6),104(8); gsum=36
c     Add 4Ge-sum; imps=13
         rcmp=rcmpf2(19,it1)+tconst*(rcmpf2(19,it2)-rcmpf2(19,it1))
         rcsum=rcsumf2(13,it1)+tconst*(rcsumf2(13,it2)-rcsumf2(13,it1))
         rcmp=rcmp+rcsum
         recof2(101)=rcmp*0.3333
         recof2(102)=rcmp*0.2778
         recof2(103)=rcmp*0.1667
         recof2(104)=rcmp*0.2222
c     d 4Pe ; imp=20, il=108(6),109(4),110(2); gsum=12
c     Add 4Pe-sum; imps=7
         rcmp=rcmpf2(20,it1)+tconst*(rcmpf2(20,it2)-rcmpf2(20,it1))
         rcsum=rcsumf2(7,it1)+tconst*(rcsumf2(7,it2)-rcsumf2(7,it1))
         rcmp=rcmp+rcsum
         recof2(108)=rcmp*0.5
         recof2(109)=rcmp*0.3333
         recof2(110)=rcmp*0.1667
c     z 4So ; imp=21, il=113(4) 
         rcmp=rcmpf2(21,it1)+tconst*(rcmpf2(21,it2)-rcmpf2(21,it1))
         recof2(113)=rcmp
c     c 4De ; imp=22, il=114(8),115(2),117(4),118(6); gsum=20
c     Add 4De-sum; imps=9
         rcmp=rcmpf2(22,it1)+tconst*(rcmpf2(22,it2)-rcmpf2(22,it1))
         rcsum=rcsumf2(9,it1)+tconst*(rcsumf2(9,it2)-rcsumf2(9,it1))
         rcmp=rcmp+rcsum
         recof2(114)=rcmp*0.4
         recof2(115)=rcmp*0.1
         recof2(117)=rcmp*0.2
         recof2(118)=rcmp*0.3
c     y 4Po ; imp=23, il=116(6),125(2),128(4); gsum=12
         rcmp=rcmpf2(23,it1)+tconst*(rcmpf2(23,it2)-rcmpf2(23,it1))
         recof2(116)=rcmp*0.5
         recof2(125)=rcmp*0.1667
         recof2(128)=rcmp*0.3333
c     z 4Go ; imp=24, il=119(12),120(10),123(8),126(6); gsum=36
         rcmp=rcmpf2(24,it1)+tconst*(rcmpf2(24,it2)-rcmpf2(24,it1))
         recof2(119)=rcmp*0.3333
         recof2(120)=rcmp*0.2778
         recof2(123)=rcmp*0.2222
         recof2(126)=rcmp*0.1667
c     z 4Ho ; imp=25, il=121(14),122(12),124(10),127(8); gsum=44
         rcmp=rcmpf2(25,it1)+tconst*(rcmpf2(25,it2)-rcmpf2(25,it1))
         recof2(121)=rcmp*0.3182
         recof2(122)=rcmp*0.2727
         recof2(124)=rcmp*0.2273
         recof2(127)=rcmp*0.1818
c     z 4Io ; imp=26, il=129(16),130(10),131(14),132(12); gsum=52
c     Add 4Io-sum; imps=17
         rcmp=rcmpf2(26,it1)+tconst*(rcmpf2(26,it2)-rcmpf2(26,it1))
         rcsum=rcsumf2(17,it1)+tconst*(rcsumf2(17,it2)-rcsumf2(17,it1))
         rcmp=rcmp+rcsum
         recof2(129)=rcmp*0.3077
         recof2(130)=rcmp*0.1923
         recof2(131)=rcmp*0.2692
         recof2(132)=rcmp*0.2308
c     y 4Do ; imp=27, il=133(8),138(6),139(2),140(4); gsum=20
         rcmp=rcmpf2(27,it1)+tconst*(rcmpf2(27,it2)-rcmpf2(27,it1))
         recof2(133)=rcmp*0.4
         recof2(138)=rcmp*0.3
         recof2(139)=rcmp*0.1
         recof2(140)=rcmp*0.2
c     y 4Fo ; imp=28, il=134(8),135(6),136(10),137(4); gsum=28
         rcmp=rcmpf2(28,it1)+tconst*(rcmpf2(28,it2)-rcmpf2(28,it1))
         recof2(134)=rcmp*0.2857
         recof2(135)=rcmp*0.2143
         recof2(136)=rcmp*0.3571
         recof2(137)=rcmp*0.1429
c     x 4Do ; imp=29, il=141(8),142(6),143(4),144(2); gsum=20
         rcmp=rcmpf2(29,it1)+tconst*(rcmpf2(29,it2)-rcmpf2(29,it1))
         recof2(141)=rcmp*0.4
         recof2(142)=rcmp*0.3
         recof2(143)=rcmp*0.2
         recof2(144)=rcmp*0.1
c     y 4Go ; imp=30, il=145(12),146(10),147(8),148(6); gsum=36
         rcmp=rcmpf2(30,it1)+tconst*(rcmpf2(30,it2)-rcmpf2(30,it1))
         recof2(145)=rcmp*0.3333
         recof2(146)=rcmp*0.2778
         recof2(147)=rcmp*0.2222
         recof2(148)=rcmp*0.1667
c     x 4Go ; imp=31, il=149(12),150(10),151(8),153(6); gsum=36
         rcmp=rcmpf2(31,it1)+tconst*(rcmpf2(31,it2)-rcmpf2(31,it1))
         recof2(149)=rcmp*0.3333
         recof2(150)=rcmp*0.2778
         recof2(151)=rcmp*0.2222
         recof2(153)=rcmp*0.1667
c     x 4Fo ; imp=32, il=152(10),154(8),157(6),159(4); gsum=28
         rcmp=rcmpf2(32,it1)+tconst*(rcmpf2(32,it2)-rcmpf2(32,it1))
         recof2(152)=rcmp*0.3571
         recof2(154)=rcmp*0.2857
         recof2(157)=rcmp*0.2143
         recof2(159)=rcmp*0.1429
c     y 4Ho ; imp=33, il=155(14),156(12),158(10),160(8); gsum=44
c     Add 4Ho-sum; imps=16
         rcmp=rcmpf2(33,it1)+tconst*(rcmpf2(33,it2)-rcmpf2(33,it1))
         rcsum=rcsumf2(16,it1)+tconst*(rcsumf2(16,it2)-rcsumf2(16,it1))
         rcmp=rcmp+rcsum
         recof2(155)=rcmp*0.3182
         recof2(156)=rcmp*0.2727
         recof2(158)=rcmp*0.2273
         recof2(160)=rcmp*0.1818
c     w 4Po ; imp=34, il=161(6),162(4),164(2); gsum=12
         rcmp=rcmpf2(34,it1)+tconst*(rcmpf2(34,it2)-rcmpf2(34,it1))
         recof2(161)=rcmp*0.5
         recof2(162)=rcmp*0.3333
         recof2(164)=rcmp*0.1667
c     w 4Fo ; imp=35, il=163(4),165(6),166(8),169(10); gsum=28
         rcmp=rcmpf2(35,it1)+tconst*(rcmpf2(35,it2)-rcmpf2(35,it1))
         recof2(163)=rcmp*0.1429
         recof2(165)=rcmp*0.2143
         recof2(166)=rcmp*0.2857
         recof2(169)=rcmp*0.3571
c     w 4Do ; imp=36, il=167(2),168(6),170(8),171(4); gsum=20
         rcmp=rcmpf2(36,it1)+tconst*(rcmpf2(36,it2)-rcmpf2(36,it1))
         recof2(167)=rcmp*0.1
         recof2(168)=rcmp*0.3
         recof2(170)=rcmp*0.4
         recof2(171)=rcmp*0.2
c     v 4Do ; imp=37, il=172(2),173(4),174(6),175(8); gsum=20
         rcmp=rcmpf2(37,it1)+tconst*(rcmpf2(37,it2)-rcmpf2(37,it1))
         recof2(172)=rcmp*0.1
         recof2(173)=rcmp*0.2
         recof2(174)=rcmp*0.3
         recof2(175)=rcmp*0.4
c     w 4Go ; imp=38, il=176(6),177(8),178(10),179(12); gsum=36
c     Add 4Go-sum; imps=14
         rcmp=rcmpf2(38,it1)+tconst*(rcmpf2(38,it2)-rcmpf2(38,it1))
         rcsum=rcsumf2(14,it1)+tconst*(rcsumf2(14,it2)-rcsumf2(14,it1))
         rcmp=rcmp+rcsum
         recof2(176)=rcmp*0.1667
         recof2(177)=rcmp*0.2222
         recof2(178)=rcmp*0.2778
         recof2(179)=rcmp*0.3333
c     y 4So ; imp=39, il=180(4); gsum=4
c     Add 4So-sum; imps=6
         rcmp=rcmpf2(39,it1)+tconst*(rcmpf2(39,it2)-rcmpf2(39,it1))
         rcsum=rcsumf2(6,it1)+tconst*(rcsumf2(6,it2)-rcsumf2(6,it1))
         rcmp=rcmp+rcsum
         recof2(180)=rcmp
c     u 4Po ; imp=40, il=181(2),182(4),183(6); gsum=12
c     Add 4Po-sum; imps=8
         rcmp=rcmpf2(40,it1)+tconst*(rcmpf2(40,it2)-rcmpf2(40,it1))
         rcsum=rcsumf2(8,it1)+tconst*(rcsumf2(8,it2)-rcsumf2(8,it1))
         rcmp=rcmp+rcsum
         recof2(181)=rcmp*0.1667
         recof2(182)=rcmp*0.3333
         recof2(183)=rcmp*0.5
c     t 4Do ; imp=41, il=184(2),185(4),186(6),187(8); gsum=20
c     Add 4Do-sum; imps=10
         rcmp=rcmpf2(41,it1)+tconst*(rcmpf2(41,it2)-rcmpf2(41,it1))
         rcsum=rcsumf2(10,it1)+tconst*(rcsumf2(10,it2)-rcsumf2(10,it1))
         rcmp=rcmp+rcsum
         recof2(184)=rcmp*0.1
         recof2(185)=rcmp*0.2
         recof2(186)=rcmp*0.3
         recof2(187)=rcmp*0.4
c     t 4Fo ; imp=42, il=188(4),189(6),190(10),191(8); gsum=28
c     Add 4Fo-sum; imps=12
         rcmp=rcmpf2(42,it1)+tconst*(rcmpf2(42,it2)-rcmpf2(42,it1))
         rcsum=rcsumf2(12,it1)+tconst*(rcsumf2(12,it2)-rcsumf2(12,it1))
         rcmp=rcmp+rcsum
         recof2(188)=rcmp*0.1429
         recof2(189)=rcmp*0.2143
         recof2(190)=rcmp*0.3571
         recof2(191)=rcmp*0.2857
c     
c     Add the remaining recombinations to the 4 highest levels 
c     (188 - 191; t4F) weighted
c     with their statistical weights so that the total recombination 
c     coefficient is correct ; gsum=28 (4,6,10,8)
c     
         recsum=0.d0
         do il=1,n
            recsum=recsum+recof2(il)
         enddo
         recadd=rctfe2-recsum
         if (recadd.lt.0.) then
            write(50,*) 'Error in calc. rec. coeff. for Fe II !!'
         endif
c     add the remaining recomb. contribution to
c     all the levels, not only the highest
         gsum=0.d0
         do il=1,n
            gsum=gsum+gf2(il)
         enddo
         do il=1,n
            recof2(il)=recof2(il)+recadd*gf2(il)/gsum
         enddo
c
      endif 
      return
      end


      
      SUBROUTINE SIMP(A,B,N,S)
      IMPLICIT REAL*8(A-H,O-Z)
C     CALC. ESCAPE PROB. INTEGRAL
      H=(B-A)/REAL(N-1)
      NA=N-2
      S1=0.
      DO 100 I=2,NA,2
      X1=H*REAL(I-1)
      X2=H*REAL(I)
  100 S1=S1+16.*FU(X1)+14.*FU(X2)
      S1=S1+7.*FU(A)+7.*FU(B)
      S=H*S1/15.+H**2*(FP(A)-FP(B))/15.
      RETURN
      END

      DOUBLE PRECISION FUNCTION FU(X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/A2/SI,TA,TAUL
      F1=ABS(1.+SI*X*X)
      EX=TA/F1
      IF(EX.LT.30.) GOTO 200
      EX=30.
200   CONTINUE
      FU=F1*(1.-EXP(-EX))/TA
      IF(EX.GT.1.E-5) GOTO 300
      FU=F1*EX/TA
300   CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION FP(X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/A2/SI,TA,TAUL
      Q=1.+SI*X*X
      QAB=ABS(Q)
      EX=TA/QAB
      IF(EX.LT.30.) GOTO 200
      EX=30.
200   CONTINUE
      FP=2.*SI*X*(1.-EXP(-EX)*(1.+TA/QAB))/TA
      IF(EX.GT.1.E-5) GOTO 300
      FP=2.*SI*X*(-TA/QAB+EX+EX*TA/QAB)/TA
300   CONTINUE
      IF(Q.GT.0.) GOTO 100
      FP=-FP
  100 CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION ESCCO(TE,TAU)
      IMPLICIT REAL*8(A-H,O-Z)
      TEV=TE/1.1609E4
      A=13.597/TEV
      IF(TAU.GT.A/3) ESCCO=0.5*(A/(3.*TAU))**.25*EXP(-(A/3.)**.75*
     &TAU**.25)
      IF(TAU.LE.A/3.) ESCCO=0.5*EXP(-TAU)
      RETURN
      END

      SUBROUTINE SECEX(X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SECX/CSEC(20),DCSDX(20),CISEC(20),DCISDX
      DIMENSION C(20),CP(20),CI(20)
      XE=X
C     IF(X.GT.1.) XE=0.99
CF 1016
      IF(X.lT.0.) XE=1.e-4
      CALL EXSEC(XE,C,CI)
      DO 10 I=1,20
10    CSEC(I)=C(I)
      DO 20 I=1,20
20    CISEC(I)=CI(I)
      RETURN
      END

      SUBROUTINE SEC(ESEC,XEL,HEFF,SEL)
      IMPLICIT REAL*8(A-H,O-Z)
C      ************************************************************
C      ******
C      THIS ROUTINE CALCULATES THE HEATING EFFICIENCY AND NUMBER OF
C      ADDITIONAL IONIZATIONS DUE TO SECONDARY ELECTRONS, WITH
C      ENERGY ESEC. THE FITS ARE INTERPOLATIONS TO MONTE-CARLO CALC-
C      CULATIONS.
C      ******
C      ************************************************************
      Y=LOG10(XEL)
      E=ESEC
      HEFF=1.
      SEL=1.
      IF(E.LT.10.2) GOTO 10
      IF(XEL.GE.1.) GOTO 10
      IF(Y.LT.-0.5) GOTO 11
      HEFF=81.*Y/E**2.-10.3*Y/E+1.+.2214*Y
      GOTO 10
  11  IF(Y.LT.-1.) GOTO 12
      HEFF=(28.6+138.2*Y)/E**2.+(-4.0-18.3*Y)/E+1.1104+.4422*Y
      GOTO 10
  12  IF(Y.LT.-1.5) GOTO 13
      HEFF=(-169.4-59.8*Y)/E**2.+(17.5+3.2*Y)/E+.9584+.2902*Y
      GOTO 10
  13  IF(Y.LT.-2.0) GOTO 14
      HEFF=(-149.7-46.6*Y)/E**2.+(13.9+0.8*Y)/E+1.1015+.3856*Y
      GOTO 10
  14  IF(Y.LT.-2.5) GOTO 15
      HEFF=(-386.3-164.9*Y)/E**2.+(41.6+14.6*Y)/E+.7027+.1862*Y
      GOTO 10
  15  IF(Y.LT.-3.0) GOTO 16
      HEFF=(138.3+44.9*Y)/E**2.+(-3.4-3.4*Y)/E+.4612+.0896*Y
      GOTO 10
  16  IF(Y.LT.-4.0) GOTO 17
      HEFF=(-3.9-2.5*Y)/E**2.+(13.3+2.2*Y)/E+.3313+.0463*Y
      GOTO 10
  17  HEFF=6.1/E**2.+4.5/E+.1461
  10  IF(HEFF.GT.1.) HEFF=1.
      RETURN
      END

      SUBROUTINE EXSEC(X,CEX,CION)
C
C     1 = HEI, 2 = O I,  3 = O II,   4 = CA I, 5 = CA II, 6 = MG I
C     7 = C I, 8 = NA I, 9 = MG II, 10 = S I, 11 = SI I  12 = NE I
C    13 = NE II 14 = AR I 15 = AR II 16 = FE I 17 = FE II18 = H I
C     ALL RATES SHOULD BE DIVIDED BY DENS**2
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 ION
      COMMON/TPAR/RIN,DRQ,R,TDAY
      COMMON/GAMPAR/GAMLUM,GAELH,GAHE
      COMMON/HYDROS/HSCALE,DEN1,XEL,AMEAN
      COMMON/ABUN/AB(20)
      DIMENSION CEX(20),CION(20),ION(20)
      dimension axsi(20),xa(2,100)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DO K=1,20
      CION(K)=1.e-30
      ENDDO
      IF(X.GT.1.) GOTO 6543
      XL=LOG10(X)
c     HYDROGEN DOMINATED
      ESEC=1.E3
c      CALL SEC(ESEC,X,HEFF,SEL)
c      IF(AB(2).GT.0.3) GOTO 60
c     OXYGEN DOMINATED
      ION(16)=-1.523
      ION(17)=-1.523
      IF(XL.GT.-1.00000) GOTO 50
      HE= 0.85012+0.27247*XL+0.02703*XL**2
      ION( 1)=-2.44824-0.09560*XL-0.01888*XL**2
      ION( 2)=-2.05749-0.24897*XL-0.04278*XL**2
      ION( 3)=-2.46878-0.15229*XL-0.03719*XL**2
      ION( 4)=-1.66678-0.13942*XL-0.01920*XL**2
      ION( 5)=-2.25347-0.07530*XL-0.01100*XL**2
      ION( 6)=-1.83600-0.13622*XL-0.02602*XL**2
      ION( 7)=-1.80086-0.10583*XL-0.01849*XL**2
      ION( 8)=-1.80164-0.09341*XL-0.00618*XL**2
      ION( 9)=-2.36815-0.03717*XL-0.00207*XL**2
      ION(10)=-1.62796-0.06683*XL-0.00966*XL**2
      ION(11)=-1.62177-0.10160*XL-0.01563*XL**2
      ION(12)=-1.93059-0.03954*XL-0.00591*XL**2
      ION(14)=-1.56809-0.01799*XL+0.00091*XL**2
      GOTO  60
50      IF(XL.GT.-0.30103) GOTO 55
      HE= 0.82860+0.21277*XL-0.01115*XL**2
      ION( 1)=-2.57233-0.21134*XL-0.01054*XL**2
      ION( 2)=-2.18307-0.43588*XL-0.10411*XL**2
      ION( 3)=-2.53715-0.26787*XL-0.08441*XL**2
      ION( 4)=-1.74414-0.25944*XL-0.06185*XL**2
      ION( 5)=-2.33686-0.03342*XL+0.11426*XL**2
      ION( 6)=-1.95195-0.34338*XL-0.11724*XL**2
      ION( 7)=-1.99154-0.51203*XL-0.23401*XL**2
      ION( 8)=-1.96099-0.41032*XL-0.16374*XL**2
      ION( 9)=-2.52601-0.28823*XL-0.09526*XL**2
      ION(10)=-1.82242-0.41406*XL-0.16242*XL**2
      ION(11)=-1.75820-0.30483*XL-0.08243*XL**2
      ION(12)=-2.08786-0.17835*XL+0.01254*XL**2
      ION(14)=-1.73749-0.13983*XL+0.04846*XL**2
      GOTO  60
55      IF(XL.GT. 0.00004) GOTO 56
      HE= 0.83097+0.19072*XL-0.11060*XL**2
      ION( 1)=-2.62139-0.55947*XL-0.62560*XL**2
      ION( 2)=-2.23986-0.77837*XL-0.61517*XL**2
      ION( 3)=-2.58118-0.55265*XL-0.54455*XL**2
      ION( 4)=-1.82371-0.62313*XL-0.39196*XL**2
      ION( 5)=-2.45797-0.50263*XL-0.10790*XL**2
      ION( 6)=-1.98587-0.47955*XL-0.19531*XL**2
      ION( 7)=-1.99122-0.71637*XL-0.91636*XL**2
      ION( 8)=-1.98230-0.54215*XL-0.36654*XL**2
      ION( 9)=-2.59339-0.28247*XL+0.66738*XL**2
      ION(10)=-1.83154-0.49657*XL-0.33597*XL**2
      ION(11)=-1.80315-0.28671*XL+0.47372*XL**2
      ION(12)=-2.15078-0.52567*XL-0.44686*XL**2
      ION(14)=-1.82305-0.46815*XL-0.09795*XL**2
      GOTO  60
56    HE= 0.83097+0.35645*XL+0.19409*XL**2
      ION( 1)=-2.62141-0.08512*XL-0.39237*XL**2
      ION( 2)=-2.23989-0.01840*XL-0.47092*XL**2
      ION( 3)=-2.58120-0.09856*XL-1.10769*XL**2
      ION( 4)=-1.82374+0.00649*XL-0.86178*XL**2
      ION( 5)=-2.45799-0.11213*XL-0.27513*XL**2
      ION( 6)=-1.98588-0.16311*XL-0.29836*XL**2
      ION( 7)=-1.99125-0.01754*XL-0.60476*XL**2
      ION( 8)=-1.98231-0.15042*XL-0.18665*XL**2
      ION( 9)=-2.59340-0.10417*XL-0.68904*XL**2
      ION(10)=-1.83156-0.00896*XL-0.82896*XL**2
      ION(11)=-1.80315-0.12404*XL-0.50833*XL**2
      ION(12)=-2.15079-0.19971*XL-0.15258*XL**2
      ION(14)=-1.82306-0.18302*XL-0.30904*XL**2
60    CONTINUE
      IF(AB(2).LT.0.3) GOTO 120
C     HE DOMINATED
      IF(XL.GT.-2.00000) GOTO 31
      HE=+0.76606+0.22527*XL+0.02330*XL**2
      ION( 1)= -1.98021 -0.16043*XL -0.01839*XL**2
      ION( 2)= -1.65366 -0.24070*XL -0.03785*XL**2
      ION( 3)= -1.66609 +0.03123*XL +0.00749*XL**2
      ION( 4)= -1.46773 -0.27453*XL -0.02831*XL**2

      ION( 6)= -1.59638 -0.24246*XL -0.03757*XL**2
      ION( 7)= -1.25043 +0.01179*XL +0.00925*XL**2

c      ION( 9)= -2.46719 -0.41284*XL -0.05766*XL**2
      ION(10)= -0.75351 +0.21784*XL +0.03691*XL**2
      ION(11)= -1.30055 -0.10029*XL -0.00713*XL**2
      ION(12)= -1.50910 -0.03468*XL -0.00268*XL**2
      ION(14)= -1.61882 -0.36051*XL -0.05624*XL**2
      goto 120
31    IF(XL.GT.-1.00000) GOTO 32
      HE=+0.99317+0.44060*XL+0.07419*XL**2
      ION( 1)= -2.20200 -0.37718*XL -0.07132*XL**2
      ION( 2)= -1.31674 +0.09051*XL +0.04353*XL**2
      ION( 3)= -1.81405 +0.01724*XL +0.03748*XL**2
      ION( 4)= -1.29513 -0.15595*XL -0.01217*XL**2

      ION( 6)= -0.56756 +1.01425*XL +0.33358*XL**2
      ION( 7)= -1.52973 -0.24367*XL -0.04865*XL**2

c      ION( 9)= -2.69182 -0.37885*XL +0.01549*XL**2
      ION(10)= -0.39997 +0.96367*XL +0.32144*XL**2
      ION(11)= -2.52918 -1.96409*XL -0.63188*XL**2
      ION(12)= -1.78415 -0.34629*XL -0.08972*XL**2
      ION(14)= -2.13126 -1.14584*XL -0.32079*XL**2
      goto 120
32    IF(XL.GT.-0.30103) GOTO 33
      HE=+1.00812+0.42881*XL+0.04744*XL**2
      ION( 1)= -2.34930 -0.65321*XL -0.20005*XL**2
      ION( 2)= -1.78821 -0.76416*XL -0.33967*XL**2
      ION( 3)= -2.14954 -0.73874*XL -0.38300*XL**2
      ION( 4)= -2.30073 -2.35916*XL -1.20978*XL**2
      ION( 6)= -1.79739 -0.76700*XL -0.21785*XL**2
      ION( 7)= -1.73076 -0.73382*XL -0.33778*XL**2
c      ION( 9)= -0.60010 +4.95778*XL +3.26041*XL**2
      ION(10)= -2.15632 -2.52493*XL -1.41081*XL**2
      ION(11)= -1.23438 +0.14190*XL +0.17932*XL**2
      ION(12)= -1.87810 -0.58727*XL -0.23676*XL**2
      ION(14)= -1.65906 -0.97742*XL -0.62458*XL**2
      goto 120
33    IF(XL.GT.-0.04576) GOTO 120
      HE=+0.99332+0.34986*XL-0.05151*XL**2
      ION( 1)= -2.34542 -0.65243*XL -0.24031*XL**2
      ION( 2)= -1.73438 -0.46256*XL +0.06821*XL**2
      ION( 3)= -2.13966 -0.33522*XL +0.84837*XL**2
      ION( 4)= -1.55621 -0.89530*XL -4.56283*XL**2
      ION( 6)= -1.59713 -1.47803*XL -4.78974*XL**2
      ION( 7)= -1.71968 -0.91049*XL -1.04686*XL**2
c      ION( 9)= -2.60086 -4.70150*XL -6.74825*XL**2
      ION(10)= -1.44071 -1.10230*XL -4.58188*XL**2
      ION(11)= -1.46024 -0.56482*XL +0.32409*XL**2
      ION(12)= -1.91839 -0.97600*XL -1.08341*XL**2
      ION(14)= -1.56207 -1.24450*XL -2.58211*XL**2
c!!!!!
120   IF(AB(1).LT.1.e-30) GOTO 122
C     HYDROGEN DOMINATED PLASMA (SHULL AND VAN STEENBERG 85)
c!!      IF(X.LT.1.D0) HE=0.9971*(1.-(1.-X**0.2663)**1.3163)
      IF(X.LT.1.D0) HIION=0.3908*(1.-X**0.4092)**1.7592
c!!      IF(X.GE.1.) HE=1.
      IF(X.GE.1.) HIION=1.E-10
      ION(18)=-LOG10(13.6/HIION)
C     CORRECT FOR HE ABUNDANCE 0.1 OF H IN THE S&S FIT
      HEION=1.E-10
      IF(X.LT.1.D0) HEION=0.0554*(1.-X**.4614)**1.6660+1.
c!!      ION(1)=-LOG10(24.58/HEION)
122   CONTINUE
      ION(13)=ION(3)
      ION(15)=ION(3)
C     HEATING EFFICIENCY NOT CORRECTED FOR IRON 
      call ionpot(x,axsi)
      sumh=0.
      do i=1,18
      sumh=axsi(i)/10**(-ion(i))+sumh
      enddo
      hyex=.4766*(1.-x**.2735)**1.5221
c      hyex=axsi(18)*hyex/13.6
      sumt=sumh+he
      do i=1,18
      ion(i)=ion(i)-log10(sumt)
      enddo
      he=he/sumt
      sumh=0.
      do i=1,18
      sumh=axsi(i)/10**(-ion(i))+sumh
      enddo
      DO 70 I=1,20
 70   CION(I)=GAHE*10.**ION(I)/ELCH
      GAELH=HE*GAHE
6543  CONTINUE
      RETURN
      END

      SUBROUTINE ionpot(x,axsi)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      COMMON/ABU/XA(2,nio)
      COMMON/ABUN/AB(20)
      dimension axsi(20)
      axsi(1)=ab(2)*xa(2,2)*24.5876
      axsi(2)=ab(3)*xa(2,14)*13.614
      axsi(3)=ab(3)*xa(2,15)*35.108
      axsi(4)=ab(11)*xa(2,50)*6.09
      axsi(5)=ab(11)*xa(2,51)*11.8
      axsi(6)=ab(7)*xa(2,39)*7.6
      axsi(7)=ab(4)*xa(2,12)*11.264
      axsi(8)=ab(10)*xa(2,53)*5.12
      axsi(9)=ab(7)*xa(2,40)*15.0
      axsi(10)=ab(12)*xa(2,56)*10.3
      axsi(11)=ab(6)*xa(2,25)*8.1
      axsi(12)=ab(13)*xa(2,72)*21.559
      axsi(13)=ab(13)*xa(2,73)*41.07
      axsi(14)=ab(14)*xa(2,82)*15.7
      axsi(15)=ab(14)*xa(2,83)*27.8
      axsi(16)=ab(8)*xa(2,42)*7.8
      axsi(17)=ab(8)*xa(2,43)*16.2
      axsi(18)=ab(1)*xa(2,1)*13.595
      axsi(19)=0.
      axsi(20)=0.
      return
      end


      SUBROUTINE ESCAPE(T0,T0T,WL,A21,VTERM,BE,DBEDTA)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
C     EVALUATE ESCAPE PROBABILITY IN BOTH DIRECTIONS
      CALL ESCAP(T0,WL,A21,VTERM,BE1,DBEDTA1)
      BE2=0.
      DBEDTA2=0.
c!!!!      IF(ISTAT.EQ.1) CALL ESCAP(T0T,WL,A21,VTERM,BE2,DBEDTA2)
      BE=BE1+BE2
      DBEDTA=DBEDTA1+DBEDTA2
      RETURN
      END

      SUBROUTINE ESCAP(T0,WL,A21,VTERM,BE,DBEDTA)
      IMPLICIT REAL*8(A-H,O-Z)
C 
C     ESCAPE PROBABILTY
C
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      IF(ISTAT.EQ.1) THEN
C      ESCAPE PROB. FROM FERLAND
c!!
       IF(VTERM.LE.0.) VTERM=1.285E6        
        A=7.96E-10*WL*A21/VTERM
        AT=A*T0
        IF(AT.GT.1.) THEN
          B=1.6 + 3./((2.*A)**.12*(1.+AT))
          DBDT = -3.*A/((2.*A)**.12 * (1. + AT)**2)
        ELSEIF(AT.GT.0.) THEN
          B=1.6 + 3./((2.*A)**.12*(1.+1./SQRT(AT)))
          DBDT = 3./(2.*SQRT(AT)*T0*(2.*A)**.12*(1.+1./SQRT(AT))**2)
        ELSE
          B=1.6
          DBDT=0.
        ENDIF
        IF(B.GT.5) THEN
          B=5.
          DBDT=0.
        ENDIF
        BE = 1./(1. + B * T0)
        DBEDTA = - BE**2*(DBDT*T0 + B)
C!!     ONLY ESCAPE FROM ONE SIDE OF THE SLAB
        BE=BE/2.
        DBEDTA=DBEDTA/2.
        GOTO 123
C     ESCAPE PROB. FROM FERLAND AND NETZER
      IF(T0.LT.1) THEN
        Q=1.11*T0**0.826
        BE=1./(1.+Q)
        DQDT=0.
        IF(T0.GT.0.) DQDT=0.826*Q/T0
        DBEDTA=-DQDT/(1.+Q)**2
      ELSE
        Q=1.11*T0**1.071/(1.+(LOG10(T0)/5.)**5)
        BE=1./(1.+Q)
        DQDT=1.071*Q/T0-6.26E-4*Q**2*(LOG10(T0))**4/T0**2.071
        DBEDTA=-DQDT/(1.+Q)**2
      ENDIF
C     ESCAPE PROB. FROM KWAN AND KROLIK
      GOTO 123
        IF(T0.LT.1.E-5) THEN
          BE=1.-T0
          DBEDTA=-1.
        ELSEIF(T0.LE.1.) THEN
          BE=(1.-EXP(-2.*T0))/(2.*T0)
          DBEDTA=-BE/T0+EXP(-2.*T0)/T0
        ELSE
C     DAMPING CONSTANT
C     A VALID ONLY FOR H AND T=1.E4K. SCALES AS 1./SQRT(T)
c!!   wl in A? see also calling subroutines!
        A=7.96E-10*WL*A21/VTERM
        SQPI=1.77245
        TC=(SQPI*(-LOG(A)+3.102))/A
        BE=1./(SQPI*T0*(1.2+0.5*SQRT(LOG(T0))/(1.+T0/TC)))
        DBEDTA=-BE*(1./T0+0.5*SQPI*SQRT(LOG(T0))*BE*
     &            (1./(2.*LOG(T0))-1./(1.+T0/TC))/(1.+T0/TC))
        ENDIF
C!!   ONLY ESCAPE FROM ONE SIDE OF THE SLAB
      BE=BE/2.
      DBEDTA=DBEDTA/2.
      ELSE
C     SOBOLEV FOR V = R
       IF(T0.LT.1.E-5) THEN
        BE=1.-T0/2.+T0**2/6.
        DBEDTA=-0.5
       ELSE
        IF(T0.LT.100.) THEN
          BE=(1.-EXP(-T0))/T0
          DBEDTA=((T0+1.)*EXP(-T0)-1.)/T0**2
        ELSE
          BE=1./T0
          DBEDTA=-1./T0**2
        ENDIF
       ENDIF
      ENDIF
123   CONTINUE
      RETURN
      END

c$$$      SUBROUTINE RECHC(J,TE,ALI)
c$$$      IMPLICIT REAL*8(A-H,O-Z)
c$$$c     RECOMBINATION COEFFICIENTS ADJUSTED FOR N GE 3 TO GIVE APPROXIMATE
c$$$C     RECOMB. EMISSION SPECTRUM FOR CASE B.SEE PROGRAM atdath.f 
c$$$        T4=TE/1.E4
c$$$        IF(J.EQ.1) ALI=1.58E-13/T4**0.549
c$$$        IF(J.EQ.2) ALI=7.69E-14/T4**0.668
c$$$        IF(J.EQ.3) ALI=7.10E-14/T4**0.769
c$$$        IF(J.EQ.4) ALI=4.63E-14/T4**0.843
c$$$        IF(J.EQ.5) ALI=4.53E-14/T4**1.05
c$$$      RETURN
c$$$      END

      SUBROUTINE MULTISIMP(KI,iel,ion,IDEP,Z,TE,XI)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N
      include "parameters.h"
      COMMON/COLD/OPT(nio),COLTOT(nio),COLTOTT(nio),SRED(NFEL+100)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/A7/C3,C33
      COMMON/PHY/DENS(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/IND/IK
      COMMON/A19/EM(NL,NL),ESC(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      COMMON/EQUIH/FRH(2,500),WEQ(500),SR(NFEL),WOBS(NL,NL)
      COMMON/NLEV/e00,NION,NLEV,NP1H,NMAX,NMIN
      COMMON/A14/CS(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),WLN(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),PH(NL),DRDT(NL),RECT(5,NL),PHET(NL)
      COMMON/RECAL/RECO(NL)
      common/w/wlin(501),wl(NL,NL),taul(501)
      common/kmaxpop/kmaxp(14,27)
      DIMENSION BE(NL,NL),XI(NL),C(NL,NL)
      DIMENSION AA(nlp1,nlp1),X(nl),RL(5,5),XO(NL)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DEN=DENS(IK)
      DENEL=DEN*DEL(IK)
      TEV=TE/1.1609E4

      DO I=1,NLP1
         DO J=1,NLP1
            AA(I,J)=0.
         ENDDO
      ENDDO
C     ADD ONE LEVEL FOR CONTINUUM
      N=NLEV+1
      E(N)=E00

      DO I=1,N
        DO J=1,N
          WL(I,J)=0.0
          IF(I.NE.J.AND.E(I).NE.E(J)) THEN
            WL(I,J)=DABS(12398.54/(E(I)-E(J)))
          ENDIF
        ENDDO
      ENDDO
      TEV=TE/1.1609D4
      TIME=8.64E4*TDAY
      DO I=1,NLEV
        DO J=1,NLEV
           C(J,I)=DENEL*CS(J,I)           
        ENDDO
      ENDDO

      DO L=1,10
         DO I=1,N
            DEI=X(I)*Z*DEN
            DO J=1,N
               AA(I,J)=0.
               IF(I.EQ.N) THEN
                  AA(I,J)=1.D0
               ELSE
                  BE(I,I)=0.
                  IF(I.NE.J) THEN
C     
C     CALCULATE THE ESCAPE PROBABILITY.
C     
                     IF(ISTAT.EQ.1) THEN
                        T0=C33*COLTOT(KI)*X(I)*A(J,I)*G(J)*WL(J,I)**3./G(I)
                        T0T=T0
                     ELSE
                        T0=1.e-24*WL(J,I)**3*A(J,I)*G(J)*DEI*TIME/(8.*PI*G(I))
                     ENDIF
                     IF(L.EQ.1) then
                        T0=.0
                     endif
C!!   ASSUME MEAN ATOMIC WEIGHT = 20. 
                     AU=20.
                     VTERM=1.285E6*SQRT(TE/(AU*1.E4))
                     CALL ESCAPE(T0,T0T,WL(J,I),A(J,I),VTERM,BE(J,I),DBEDTA)
                     ESC(J,I)=BE(J,I)
                  ENDIF
                  IF(I.EQ.J) THEN
                     S=0.
                     DO K=1,Nlev
                        S=S+C(I,K)+BE(I,K)*A(I,K)
                     ENDDO
C     ADD PHOTOIONIZATION
                     AA(I,I)=-S-PH(I)
                  ELSEIF(J.LT.N) THEN
                     AA(I,J)=BE(J,I)*A(J,I)+C(J,I)
                  ELSEIF(J.EQ.N) THEN
C     RECOMBINATION CONTRIBUTION
                     AA(I,N)=RECO(I)*DENEL
                  ENDIF     
               ENDIF     
            ENDDO
         ENDDO

C     NUMBER CONSERVATION, note n=nlev + 1 (cont)
         NRC=N+1
         AA(N,NRC)=1.
         EPS=1.D-30
         NA=N
         DO I=1,N
            XO(I)=X(I)
         ENDDO

         DI=SIMUL(NA,AA,X,EPS,1,NRC)
         ERRMAX=0.
         DO I=1,N
            IF(X(I).NE.0.) THEN
               ERR=ABS((XO(I)-X(I))/X(I))
            ELSE
               ERR=1.
            ENDIF
            IF(ERR.GT.ERRMAX) ERRMAX=ERR
         ENDDO
         IF(ERRMAX.LT.0.01) GOTO 555
      ENDDO
555   CONTINUE
      DO I=1,N
         XI(I)=X(I)
      ENDDO

      K=0
      DO I=2,Nlev
         IM1=I-1
         DO J=1,IM1
            WOBS(I,J)=1.602E-12*Z*(E(I)-E(J))*X(I)*A(I,J)*BE(I,J)/DEN
            EM(I,J)=1.602E-12*Z*(E(I)-E(J))*
     &           (X(J)*C(J,I)-X(I)*C(I,J))
            IF(K.LE.399.and.a(i,j)>0.) THEN
               k=k+1
               SR(K)=0.
               WEQ(K)=WOBS(I,J)
               wlin(k)=wl(i,j)
               kmaxp(iel,ion)=k
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END



      SUBROUTINE POPchianti(iel,ion,te,XQ,rltot)
      IMPLICIT REAL*8(A-H,O-Z)
c     ion = ionel: 1 = Ca II, 2 = O I, 3 = He I, 4 = Fe II, 5 = H, 6 = Fe I,
c     7 = He II ??               

      include "parameters.h"
      parameter (nlp=30000)
      common/abl/abn(15)
      common/ionx/xion(md,14,27)
      common/ichia/ipopch
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/NION/IONQ
      COMMON/INI/INITH,INIT,INITCA,INITHE,INITFE
      COMMON/IND/IK
      COMMON/ELEC/DEL(MD)
      COMMON/PHY/DEN(MD)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &     BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/RECAL/RECO(NL)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),PH(NL),DRDT(NL),RECT(5,NL),PHET(NL)
      COMMON/EQUIH/FRH(2,500),WEH(500),SR(NFEL),WOBS(NL,NL)
      common/timecheck/time,itime
c     DIMENSION XINC(NL),XI(NL),A(NLP1,NLP1),BS(NL)
      real*8 XINC(NL),XI(NL),A(NLP1,NLP1),BS(NL)
      DATA CSAHA/2.0708d-16/


      IONQ=ion

      abxa = abn(iel)*xion(ik,iel,ion)
      ipopch=0
      if(itime.le.3.or.abxa.gt.1.e-8) then
         ipopch=1
         if(iel.eq.3.and.ion.eq.5) then
            N=49
         elseif(iel.eq.3.and.ion.eq.6) then
            N=25
         elseif(iel.eq.4.and.ion.eq.6) then
            N=49
         elseif(iel.eq.4.and.ion.eq.7) then
            N=25
         elseif(iel.eq.5.and.ion.eq.7) then
            N=49
         elseif(iel.eq.5.and.ion.eq.8) then
            N=25
         elseif(iel.eq.6.and.ion.eq.6) then
            N=180
         elseif(iel.eq.6.and.ion.eq.7) then
            N=46
         elseif(iel.eq.6.and.ion.eq.8) then
            N=40
         elseif(iel.eq.6.and.ion.eq.9) then
            N=49
         elseif(iel.eq.6.and.ion.eq.10) then
            N=36
         elseif(iel.eq.8.and.ion.eq.8) then
            N=125
         elseif(iel.eq.8.and.ion.eq.9) then
            N=46
         elseif(iel.eq.8.and.ion.eq.10) then
            N=40
         elseif(iel.eq.8.and.ion.eq.11) then
            N=49
         elseif(iel.eq.8.and.ion.eq.12) then
            N=25
         elseif(iel.eq.10.and.ion.eq.8) then
            N=72
         elseif(iel.eq.10.and.ion.eq.9) then
            N=46
         elseif(iel.eq.10.and.ion.eq.10) then
            N=125
         elseif(iel.eq.10.and.ion.eq.11) then
            N=46
         elseif(iel.eq.10.and.ion.eq.12) then
            N=40
         elseif(iel.eq.10.and.ion.eq.13) then
            N=49
         elseif(iel.eq.10.and.ion.eq.14) then
            N=25
         elseif(iel.eq.11.and.ion.eq.2) then
            N=43
         elseif(iel.eq.11.and.ion.eq.3) then
            N=49
            N=5
         elseif(iel.eq.11.and.ion.eq.13) then
            N=46
         elseif(iel.eq.11.and.ion.eq.14) then
            N=40
         elseif(iel.eq.11.and.ion.eq.15) then
            N=49
         elseif(iel.eq.11.and.ion.eq.16) then
            N=25
         elseif(iel.eq.13.and.ion.eq.9) then
            N=16
         elseif(iel.eq.13.and.ion.eq.10) then
            N=21
         elseif(iel.eq.13.and.ion.eq.11) then
            N=89
         elseif(iel.eq.13.and.ion.eq.12) then
            N=3
         elseif(iel.eq.13.and.ion.eq.13) then
            N=86
         elseif(iel.eq.13.and.ion.eq.14) then
            N=91
         elseif(iel.eq.14.and.ion.eq.10) then
            N=172
         elseif(iel.eq.14.and.ion.eq.11) then
            N=47
         elseif(iel.eq.14.and.ion.eq.12) then
            N=143
         elseif(iel.eq.14.and.ion.eq.13) then
            N=27
         elseif(iel.eq.14.and.ion.eq.14) then
            N=40
         elseif(iel.eq.14.and.ion.eq.15) then
            N=283
c!!!  should be 283!
            n=199
         elseif(iel.eq.14.and.ion.eq.16) then
            N=49
         elseif(iel.eq.14.and.ion.eq.17) then
            N=267
c!!!  should be 267!
            n=199
         elseif(iel.eq.14.and.ion.eq.18) then
            N=337
c!!!  should be 337!
            n=199
         elseif(iel.eq.14.and.ion.eq.19) then
            N=200
c!!!  should be 200!
            n=199            
         elseif(iel.eq.14.and.ion.eq.20) then
            N=302
c!!!  should be 302!
            n=199            
         elseif(iel.eq.14.and.ion.eq.21) then
            N=200
c!!!  should be 283!
            n=199            
         elseif(iel.eq.14.and.ion.eq.22) then
            N=200
c!!!  should be 200!
            n=199            
         elseif(iel.eq.14.and.ion.eq.23) then
            N=58
         elseif(iel.eq.14.and.ion.eq.24) then
            N=40
         elseif(iel.eq.14.and.ion.eq.25) then
            N=49
         elseif(iel.eq.14.and.ion.eq.26) then
            N=25
         endif
         CALL ATDATchianti(iel,ion,te)
         NMIN=n
         NQ=n
         NP1=N+1
         NP1H=NP1
         NVAR=NP1
         NRC=NP1+1
         nmax=nmin
         xelec=del(ik)

c!! no recombination or photoionization to these

         do k=1,nl
            reco(k)=0.
            ph(k)=0.
         enddo


C
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITYXN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C
         C1=CSAHA
         DEN=DEn(ik)
         TEV=Te/1.1609E4
         

         iout=0

         zi=1.
         do i=1,nl
            xi(i)=0.
         enddo
            
         CALL MULTISIMPq(iel,ion,Te,XI,rltot)

         BB(NP1)=XI(NP1)
         DENEL=DEL(IK)*DEn(ik)
      else
         do k=1,500
            weh(k)=1.d-40
         enddo

         rltot = 1.d-40
      endif
      RETURN
      END

      SUBROUTINE POPO_new(RS,TSIN,ZS,XQ,PHO,IFPOP)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      COMMON/A19/EM(NL,NL),ESC(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      common/abl/abn(15) !diff
      COMMON/INUT/IUY
      COMMON/HRAT/RHYD,ZMETAL,XRQ,HYR,HEAT,COOL
      COMMON/NION/IONQ
      COMMON/A25/DHDT,DCODT,DZDT,DPDT
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),PH(NL),DRDT(NL),RECT(5,NL),PHET(NL)
      COMMON/MULTI/IHMUL,IHEMUL,IOMUL,ICAMUL,IFEMUL
      COMMON/INI/INITH,INIT,INITCA,INITHE,INITFE
      COMMON/ITER/ITE
      COMMON/SECX/CSEC(20),DCSDX(20),CISEC(20),DCISDX
      COMMON/NBA/NBACK
      COMMON/IND/IK
      COMMON/ELEC/DEL(MD)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &     BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/ELDEN/IDENS
      COMMON/PHEAT/PHE(NL)
      COMMON/BIN/B1IN,ELD
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A12/ITMAX
      COMMON/LYBFLUOR/OIFLUOR,OIFLUORH,OIFLUORO,OILYB,OIREC
      common/multisol/imultih,imultihe,imultio,imultica,imultifei,
     &     imultifeii
      DIMENSION XINC(NL),XI(NL),bbo(nl),
     &                                    A(NLP1,NLP1),BS(NL)
      DATA CSAHA/2.0708E-16/
      iel=5
      ion=1
      ifpop=0
      ITE=-1

      IONQ=2
      N=13
      NMIN=N
      NMAX=N
      do i=1,nl
         bol(i)=0.
      enddo
      CALL ATDATO_new
      DO NQ=NMIN,NMAX
         NP1=N+1
         NP1H=NP1
         NP2=N+2
         NVAR=NP2
         IIT=0
         ICO=0
         IW=1
c     IF(IK.EQ.2) INIT=1
         TS=TSIN
         IF(IIT.GT.10) STOP
         RCM=1.E15*RS
         IF(ITE.EQ.-1) DS=DE(RS)
         NRC=NP1+1
         EPS1=1.D-300
         EPS2=1.E-3
         IPRINT=2
         J=1
         IQQ=0

         IF(IK.EQ.2.AND.IUY.EQ.1) XELEC=ZS
         IF(IK.EQ.2.AND.IUY.EQ.1) XN1=0.
         IF(IK.ge.2) XELEC=DEL(IK-1)

         RCM=1.E15*RS
         IF(ITE.EQ.-1) DS=DE(RS)
C     
C     CALCULATE COLL EXCITATION RATES
C     
         IQW=0
C     
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITY XN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C     
         C1=CSAHA
         DEN=DE(RS)
         G1=9.
         G2=5.
         G3=1.
         TEV=TS/1.1609E4

         CALL COLLEX(IONq,TS)

         CALL HORATE(IONq,RS,TS)


 11      XELOLD=XELEC



         CALL SECEX(XELEC)

         T4=TS/1.E4
C!!   OTS-APP. DOes not work for external rad field!!
c     nots      PH(1)=0.
         ALO=0.
         DO J=1,N
            CALL RECOMB(IONq,J,TS,AL,RI)
            ALO=ALO+AL
         ENDDO

         RI=ALO*EXP((E(1)-E00)/TEV)*TS**1.5/(CSAHA*G(1))
         BB(1)=RI/(PH(1)+DEN*XELEC*CI(1)+CISEC(2)/DEN)
         QQ=(PH(1)/DEN+CI(1)*XELEC+CISEC(2)/DEN**2)/(ALO*ABn(5)*XQ)
         ZQ=ZS/(ABn(5)*XQ)
         BB(NP1)=SQRT(QQ+(QQ+ZQ)**2/4.)-(ZQ+QQ)/2.
         XELECn=ZS+BB(NP1)*ABn(5)*XQ
         xelec=sqrt(xelecn*xelold)

         IF(ABS(XELEC-XELOLD)/XELOLD.GT.0.01) GOTO 11
         INIT=0

         ZI=ABn(5)*XQ
         IDEP=0

         CALL MULTISIMPq(iel,ion,TS,XI,rltot)

         BB(NP1)=XI(NP1)
         DENEL=DEL(IK)*DE(RS)

         IF(IDEP.EQ.1) THEN
            DO J=1,N
               BB(J)=XI(J)
               XN(ION,J)=XI(J)*G(J)*BB(NP1)*DENEL*C1*EXP((E00-E(J))/TEV)
     &              /TS**1.5
            ENDDO
         ELSE
            DO J=1,N
               XN(ION,J)=XI(J)
               EXT=DMIN1(700.D0,(E00-E(J))/TEV)
               bbo(j)=bb(j)
               BB(J)=XI(J)/(G(J)*BB(NP1)*DENEL*C1*EXP(EXT)/TS**1.5)
               xinc(j)=bb(j)-bbo(j)
            ENDDO
         ENDIF
         DO  NK=1,NP1
            BS(NK)=BB(NK)
         enddo
         IHLO=0
         ITCON=1


C     
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITY) XN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C     
         C1=CSAHA
         DEN=DE(RS)

         TEV=TS/1.1609E4
         DO  J=1,N
            EXT=DMIN1(700.D0,(E00-E(J))/TEV)
            XN(2,J)=G(J)*BB(NP1)*DEL(IK)*C1*DEN*EXP(EXT)*BB(J)/TS**1.5
         enddo
 5555    continue
         XN1=0.
         XN2=0.
         XN3=0.
         IF(13.6/TEV.GT.700.) GOTO 778
         XN1=G1*BB(NP1)*DEL(IK)*C1*DEN*EXP(13.60/TEV)*BB(1)/TS**1.5
         XN2=G2*BB(NP1)*DEL(IK)*C1*DEN*EXP(11.657/TEV)*BB(2)/TS**1.5
         XN3=G3*BB(NP1)*DEL(IK)*C1*DEN*EXP(9.446/TEV)*BB(3)/TS**1.5
 778     PHN1=XN1*PHE(1)
         PHN2=XN2*PHE(2)
         PHN3=XN3*PHE(3)
         PHH=PHN1+PHN2+PHN3
         TOLDO=TS
         DO IS=1,NVAR
            BOL(IS)=BB(IS)
         enddo
      enddo
      IF(ITCON.EQ.0) IFPOP=1
      IF(ITCON.EQ.0) write(6,204)
      XELEC=ABn(5)*BB(NP1)*XQ+ZMETAL
c      IF(IHLO.EQ.1) CALL HLOSS(2,BB,RS,TS,XELEC,XQ,HCA)
      IF(IFPOP.EQ.1.AND.IK.EQ.1) INIT=1
C     LYMAN BETA FLUORESENCE A LA KWAN AND KROLIK
      IF(XN(5,1).GT.0.) THEN
        IF(ABn(1).GT.0.01) THEN
          OIFLUOR=1.41E6*6.65*ESC(7,8)/(5.55+1.1*ESC(7,8))
        ELSE
          OIFLUOR=0.
        ENDIF
        OIFLUORH=OIFLUOR*ABn(5)*XQ*XN(2,1)/(ABn(1)*XN(5,1))
        OIFLUORO=OIFLUOR*XN(5,3)/XN(5,1)
      ENDIF

      RETURN
201   FORMAT(1X,'MATRIX ILL-CONDITIONED OR SINGULAR',E11.4)
202   FORMAT(1X,I8,E10.4,I8,8E12.5)
204   FORMAT(1X,'NO CONVERGENCE FOR O!')
      END

      SUBROUTINE ATDATo_new
      IMPLICIT REAL*8(A-H,O-Z)
      character*8 dum
      include "parameters.h"
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GAP(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      common/initox/initoi
      common/oicoll/omoi(5,13,13),teoi(5)
      dimension eoi(13),goi(14),wloi(13,13),aoi(13,13)
      save eoi,goi,wloi,aoi
      N=NHMAX
      N=13
      E00=13.618
      if(initoi.eq.1) then         
         initoi=0
         
         do i=1,n
            do j=1,n
               aoi(i,j)=0.
            enddo
         enddo

         rewind 25
         OPEN(25,FILE='./ATDAT/OI_kb.dat')
         read(25,987)dum
         read(25,987)dum
         read(25,987)dum
 987     format(a)
         do i=1,n
            read(25,*)nn,wn,goi(i)
c     read(25,*)nn,wn,g(i),al0(i),be(i),b0(i),ga(i)
            eoi(i)=wn/8065.46d0
         enddo


C     O II

         goi(N+1)=4.

         read(25,987)dum
         read(25,987)dum

         do i=1,n-1
            do j=i+1,n
               read(25,*)i1,i2,aoi(j,i),OM5,OM10
            enddo
         enddo
         close(25)
         OPEN(25,FILE='./ATDAT/OI_coll_str_Bhatia_Kastner.dat')
         read(25,*)(teoi(k),k=1,5)
         do i=1,12
            do j=i+1,13
               read(25,*)(omoi(k,i,j),k=1,5)
            enddo
         enddo
         close(25)
      endif

      DO I=1,N
         e(i) = eoi(i)
         g(i) = goi(i)
      enddo


      DO 5395 I=1,N
         DO 5394 J=1,N
            WL(I,J)=0.0
            IF(I.EQ.J) GOTO 5394
            IF(E(I).EQ.E(J)) GOTO 5394
            WL(I,J)=ABS(12398.54/(E(I)-E(J)))
 5394    CONTINUE
 5395 CONTINUE

      DO  I=1,N
         DO J=1,N
            OM(I,J)=0.
            A(I,J)=aoi(i,j)
         enddo
      enddo

C
      SIG(1)=7.91
      SIG(2)=15.
      DO 8457 I=3,N
 8457 SIG(I)=7.91*REAL(I)
      SIG(8)=20.5
      SIG(9)=20.5
      GAP(1)=2.99
      GAP(2)=2.5
      DO 1009 I=3,N
1009  GAP(I)=3.
      GAP(8)=3.56
      GAP(9)=3.56
      DO I=1,N
         C(I,I)=0.
      enddo
      RETURN
      END

      SUBROUTINE POPFEI(RS,TSIN,ZS,XQ,PHH,IFPOP)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      common/abl/abn(15)
      common/timecheck/time,itime
      common/ionx/xion(md,14,27)
      parameter (nlp=15000)
      COMMON/EQUIH/FRH(2,500),WEQ(500),SR(NFEL),WOBS(NL,NL)
      COMMON/INUT/IUY
      COMMON/HRAT/RHYD,ZMETAL,XRQ,HYR,HEAT,COOL
      COMMON/NION/IONQ
      COMMON/A25/DHDT,DCODT,DZDT,DPDT
      COMMON/IND/IK
      COMMON/ELEC/DEL(MD)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &       BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/ELDEN/IDENS
      COMMON/PHEAT/PHE(NL)
      COMMON/BIN/B1IN,ELD
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/FELEV/NFEII,nfei
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      common/multisol/imultih,imultihe,imultio,imultica,imultifei,
     &     imultifeii
      DIMENSION XINC(NL),XI(NL)
      DATA CSAHA/2.0708E-16/
      ION=6
      IONQ=6
      nfei=121
c      nfei=122
      n=nfei
      NMIN=NFEI
      NMAX=NMIN
      CALL ATDATFEI
      call ferecomb(1,tsin)
      ifpop=0
      ITE=-1
      IIT=0

      TS=TSIN

C
C     CALCULATE COLL EXCITATION RATES
C
      IQW=0

      CALL COLLEX(IONq,TS)

      CALL MULTISIMPq(14,1,TS,XI,rltot)

      return

      END
      
      SUBROUTINE ATDATFEI
      IMPLICIT REAL*8(A-H,O-Z)
      character*80 dum
      character*2 ch2
      character*4 ch4
      character*5 ch5
      SAVE
      include "parameters.h"
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/FELEV/NFEII,nfei
      common/initat/initfeii,initfei,initsi1
      integer nlf1,nlf2,nlf3,nlf4
      parameter(nlf1=122,nlf2=192,nlf3=113,nlf4=46)
      common/cfedat/gf1(nlf1),gf2(nlf2),gf3(nlf3),gf4(nlf4)
      DIMENSION GS(NL),ES(NL),AS(NL,NL),WLS(NL,NL),OMS(NL,NL)
      ntot=121
      nfei=121      
      N=NFEI
      E00=7.870
      IF(INITFEi.EQ.0) THEN
      INITFEi=1
      OPEN(23,FILE='./ATDAT/feI_221210.dat')
      rewind 23
      do i=1,2
         read(23,987)dum
      enddo
987   format(a)
      do i=1,ntot
         read(23,*)nn,wn,gs(i)
         gf1(i)=gs(i)
 1       es(i)=wn/8065.46d0
      enddo
      DO I=1,N
         DO 5394 J=1,N
            WLS(I,J)=0.0
            IF(I.EQ.J) GOTO 5394
            IF(ES(I).EQ.ES(J)) GOTO 5394
            WLS(I,J)=ABS(12398.54/(ES(I)-ES(J)))
 5394    CONTINUE
      enddo
      DO I=1,N
         DO J=1,N
            OMS(I,J)=0.
            AS(I,J)=0.
         enddo
      enddo
c     read(23,987)dum
      read(23,987)dum
      read(23,987)dum

      do  i=2,nfei
         do  j=1,i-1
            read(23,*)i1,i2,wlq,as(i,j),OMs(J,I)
         enddo
      enddo
      close (23)
      ENDIF
      DO I=1,N
         E(I)=ES(I)
         G(I)=GS(I)
         DO J=1,N
            A(I,J)=AS(I,J)
            OM(I,J)=OMS(I,J)
            WL(I,J)=WLS(I,J)
         ENDDO
      ENDDO
C
      SIG(1)=0.
      SIG(2)=0.
      DO 8457 I=3,N
 8457 SIG(I)=0.
      GA(1)=2.99
      GA(2)=2.5
      DO 1009 I=3,N
1009  GA(I)=3.
      GA(8)=3.56
      GA(9)=3.56
      DO 165 I=1,N
165   C(I,I)=0.
      RETURN
      END
      

      
