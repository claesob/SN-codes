      SUBROUTINE RLOSS(K,WL,OM,G1,G2,A21,R,XEL,Z,TE,RL,WEQ
     &,CINT,RFL,L,WLI)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MTOT,MNI
c      PARAMETER (NFEL=1000)
      include "parameters.h"
C     **************************************************************
C     *****
C     THIS ROUTINE CALCULATES THE HEAT LOSS AND OBSERVED FLUX FROM
C     A LINE WITH THE ESCAPE PROBABILITY METHOD.
C     WL = WAVELEGTH (A)
C     OM = COLLISION STRENGTH
C     G1,G2 = STAT.WEIGHTS OF THE LOWER AND UPPER LEVELS
C     A21 = TRANSITION PROBABILITY
C     R = RADIUS (IN TERMS OF SHOCK RADIUS)
C     XEL = ELECTRON FRACTION
C     Z = TOTAL ABUNDANCE OF THE ION
C     TE = TEMPERATURE
C     RL = LOSS RATE IN ERG CM**3 /S
C     WEQ = FLUX EMITTED TO INFINITY (ERG/S)
C     CINT = OBSERVED FLUX
c     flux = vol * den**2 * xel * cint
C     RFL = RADIATION FORCE (NOT USED HERE)
C     SRED = (1-EXP(-TAU))*SOURCE FCN FOR LINE PROFILE
C     REVISED 87-05-15
C
C     *****
C     ****************************************************************
      COMMON/NBA/NBACK
      COMMON/A2/SI,T0,TAUL
      COMMON/COLD/OPT(nio),COLTOT(nio),COLTOTT(nio),SRED(NFEL+100)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/A7/C3,C33
      COMMON/A10/AV
      COMMON/SPH/ISPH
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/A11/R11,CV,FLUX
      COMMON/A5/TAU,ALFA,ENE,TSN,XL40,TXEV,RMAX
      COMMON/MPAR/MTOT,MNI,VEXP
      COMMON/HYDROS/HSCALE,DEN1,XQL,AMEAN
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      WLI=WL
      SI=AV-1.
      EN=1.4413E8/WL
      ENCGS=1.989E-8/WL
      EEV=ENCGS/ELCH
      DENEL=XEL*DE(R)
      C21=8.63E-6*OM/(G2*SQRT(TE))
      E=C21*DENEL/A21
C
C     THIS PART APPLIES FOR NO VELOCITY GRADIENT AND NO BACKGROUND
C     CONTINUUM.
C
      IF(ISTAT.EQ.1) THEN
        T0=C33*COLTOT(K)*A21*G2*WL**3./G1
        T0T=(COLTOTT(K)-COLTOT(K))*T0/COLTOTT(K)
        IF(T0T.LE.0.) T0T = 0.
      ENDIF
      OPT(L)=T0
      W=0.E0
      IF(NBACK.EQ.1) THEN 
        FIN=FMEAN(3,EEV)
      ELSEIF(NBACK.EQ.0) THEN
        FIN=0.
      ENDIF
      IF(ISTAT.Eq.0) then
C     
C     THIS PART APPLIES FOR A VELOCITY GRADIENT AND A BACKGROUND
C     CONTINUUM, FC.
C     
         ENX=EN/TSN
         IF(ISPH.EQ.1) W=.5*(1.-SQRT(1.-1./(R*R)))
         TIME=8.64E4*TDAY
C     OPT. DEPTH. COEFF = 1E-24/8 PI
         T00=3.979E-26*TIME*WL**3.*G2*A21*Z*DE(R)/G1
         T0T=0.
      endif
C
C     CALCULATE THE ESCAPE PROBABILITY.
C
      MM=0
      T0=T00
232   BE=1.
      TOP=T0
      ENTE=EN/TE
C!!   ASSUME MEAN ATOMIC WEIGHT = 20. 
      AU=20.
      IF(ISTAT.EQ.1) VTERM=1.285E6*SQRT(TE/(AU*1.E4))
      CALL ESCAPE(TOP,T0T,WL,A21,VTERM,BE,DBEDTA)
      IF(EN/TE.GT.100.) ENTE=100.
C
C     N1*G2/N2*G1
C
      RN1G2N2G1=(BE*(1.+FIN)+E)/(BE*FIN+E*EXP(-ENTE))
      RN1N2=G1*RN1G2N2G1/G2
C
C     CORRECT OPTICAL DEPTH FOR STIMULATED EMISSION
C
      IF(RN1N2.LT.10.) T0NEW=T00*(1.-1./RN1G2N2G1)/(1.+1./RN1N2)
      IF(RN1N2.LT.10.) T0=T0NEW
      IF(RN1N2.LT.10.) MM=MM+1
      IF(RN1N2.LT.10..AND.MM.LT.4) GOTO 232
C
C     COOLING RATE OF ELECTRONS TAKING THERMALIZATION INTO ACCOUNT
C
      RL=ENCGS*G2*C21*Z*EXP(-ENTE)*(1.-EXP(ENTE)/RN1G2N2G1)
     &/(G1*(1.+1./RN1N2))
      XX=G1*(A21*BE*(1.+FIN)+C21*DENEL)/(G2*(C21*DENEL*EXP(-ENTE)+
     &     BE*FIN*A21))
c      write(6,*)'c21,denel,ente,xx',c21,denel,ente,xx
      WEQ=Z*A21*ENCGS*BE*DE(R)*((1.-W)-FIN*
     &     (XX*G2/G1-1.))/(1.+XX)
c      write(6,*)'z,a21,encgs,be,w,fin,xx,g2,g1,weq ',z,a21,encgs,be,w,fin,xx,g2,g1,weq
      IF(EN/TE.LT.100.) BPL=2.031E-04*EEV**3/(EXP(EN/TE)-1.)
      IF(NBACK.EQ.1) THEN 
        FIN=FMEAN(2,EEV)
      ELSEIF(NBACK.EQ.0) THEN
        FIN=0.
      ENDIF
C     SOURCE FUNCTION INCL. STIM. EMISSION
      EP=E
      IF(EN/TE.LT.100.) EP=E*(1.-EXP(-EN/TE))
      IF(EN/TE.LT..001) EP=E*EN/TE
      SOURCE=(EP*BPL+BE*FIN)/(EP+BE)
c      write(6,*)'ep,bpl,be,fin,source ',ep,bpl,be,fin,source
C     OBSERVED INTENSITY
      EXPM1=1.-EXP(-T0)
      IF(T0.LT.1.E-3) EXPM1=T0-0.5*T0**2
      IF(ISTAT.NE.1) THEN
c     CINT=4.*PI*1.E8*EXPM1*SOURCE/(TIME*WL*DENEL*DE(R))
         CINT=4.*PI*1.E8*EXPM1*SOURCE/(TIME*WL*DENEL)
c         write(6,*)'t0,expm1,source,time,wl,cint ',t0,expm1,source,time,wl,cint
      ELSE
c     CINT=WEQ/(DENEL*DE(R))
         CINT=WEQ/DENEL
      ENDIF
C     CONVERT TO INTENSITY/ ANGSTROM
      SOURCE=3.E18*SOURCE/WL**2
c      SRED(L)=EXPM1*SOURCE
      if(abs(cint) > 0.) then
c         write(6,*)' cint > 1 in rloss ',istat,k,r,de(r),denel,wl,weq,cint
      endif
c     simple cint!!!
      c12=g2*exp(-en/te)*c21/g1
      cint2=z*a21*encgs*be*c12*denel/
     &     (be*a21 + c21*denel * (1.+g2*exp(-en/te)/g1))
c      cint2=cint

c!!!  put cont. intensity fc=0!
      fin=0.
      x2_x1=(g2/g1)*(c21*denel*exp(-en/te)+fin*a21)/(c21*denel+a21*(be+fin))
      x2=x2_x1/(1.+x2_x1)
      x1=1.-x2
      cint=a21*encgs*be*z*x2/de(r)
      cint=a21*encgs*be*z*x2/denel
      rl=z*encgs*(x1*c12 - x2*c21)

      therm=1.-c21*denel/(c21*denel+be*a21)
      if(k==26) then
c         write(6,982)wl,te,om,g1,g2,a21,t0,be,en/te,fin,z,denel,de(r),x2,
c     &        therm,rl,cint2,cint
 982     format('Si II ',f10.2,1pe12.3,20e12.3)
      endif
      RETURN
      END
      

      SUBROUTINE DIEL(AR,T0R,BR,T1R,TE,ALD)
      IMPLICIT REAL*8(A-H,O-Z)
c      REAL*4 AR,T0R,BR,T1R
C     ***************************************************************
C     ******
C     DIELECTRIC RECOMBINATION RATE ACCORDING TO FORMULA IN ALD. AND
C     PEQ.
C     ******
C     ***************************************************************
        A=DBLE(AR)
        T0=DBLE(T0R)
        B=DBLE(BR)
        T1=DBLE(T1R)
      EX=T0/TE
      EX2=T1/TE
      ALD=0.
      IF(EX2.GT.100.) EX2=100.
      IF(EX.LT.100)ALD=A*EXP(-EX)*(1.+B*EXP(-EX2))/(TE*SQRT(TE))
      IF(TE.LT.4000.) ALD=0.
c      write(6,*)'diel ',ar,t0r,br,ald
      RETURN
      END

      SUBROUTINE DIELB(AR,BR,CR,DR,FR,T4,ALDB)
      IMPLICIT REAL*8(A-H,O-Z)
c      REAL*4 AR,BR,CR,DR,FR
C     ***************************************************************
C     ******
C     DIELECTRIC RECOMBINATION RATE ACCORDING TO FORMULA IN NUSSBAUMER
C     AND STOREY (A&A 126,75 (1983))
C     ******
C     ***************************************************************
C     DIELECTRIC RECOMBINATION RATE
c      write(6,*)'t4,ar,br,cr,dr,fr ',t4,ar,br,cr,dr,fr
      A=DBLE(AR)
      B=DBLE(BR)
      C=DBLE(CR)
      D=DBLE(DR)
      F=DBLE(FR)
      EX=F/T4
      ALDB=0.
      IF(T4.GT.6.) GOTO 100
      IF(EX.GT.100.) GOTO 100
      ALDB=1.E-12*(A/T4+B+C*T4+D*T4*T4)*EXP(-F/T4)/T4**1.5
 100  CONTINUE
      IF(T4.LT..2) ALDB=0.d0
c      write(6,9)t4,a,b,c,d,f,ald
 9    format('dielb qq ',1pe12.4,10e12.4)
c      write(6,*)'dielb ',ar,br,aldb
      RETURN
      END

      SUBROUTINE FIVELEV_dp(KI,iel,ion,N,E1,E2,E3,E4,E5,G1,G2,
     &     G3,G4,G5,A21,A31,A41,A51,A32,A42,A52,A43,A53,
     &     A54,C21,C31,C41,C51,C32,C42,C52,C43,C53,C54,TE,Z,XI)
      IMPLICIT REAL*8 (A-H,O-Z)
c      REAL*4 A21,A31,A41,A51,A32,A42,A52,A43,A53,A54
c     &,C21,C31,C41,C51,C32,C42,C52,C43,C53,C54
c     &,E1,E2,E3,E4,E5,G1,G2,G3,G4,G5
      INTEGER N
      include "parameters.h"
c      PARAMETER (MD=300,MDP1=MD+1)
c      PARAMETER (NL=130,NLP1=NL+1)
c      PARAMETER (NFEL=1000)
      COMMON/COLD/OPT(nio),COLTOT(nio),COLTOTT(nio),SRED(NFEL+100)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/A7/C3,C33
      COMMON/PHY/DENS(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/IND/IK
      COMMON/A19/EM(NL,NL),ESC(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      COMMON/EQUIH/FRH(2,500),WEQ(500),SR(NFEL),WOBS(NL,NL)
      common/w/wlin(501),wl(NL,NL),taul(501)
      common/kmaxpop/kmaxp(14,27)
      DIMENSION  C(5,5),G(5),E(5),A(5,5),BE(5,5),XI(15)
      DIMENSION AA(nlp1,nlp1),X(nl),RL(5,5),xold(nl),topt(nlp1,nlp1)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DEN=DENS(IK)
      DENEL=DEN*DEL(IK)
c      write(6,9343)KI,N,E1,E2,E3,E4,E5,G1,G2,
c     &G3,G4,G5,A21,A31,A41,A51,A32,A42,A52,A43,A53,
c     &     A54,C21,C31,C41,C51,C32,C42,C52,C43,C53,C54,TE,Z
c      write(6,*)'den i five lev ',iel,ion,ik,den      
 9343 format('fivelev ',2i5,1pe12.3,50e12.3)
      NM1=N-1
      NP1=N+1
      DO I=1,NP1
         DO  J=1,NP1
            AA(I,J)=0.d0
         enddo
      enddo
      DO I=1,5
         DO  J=1,5
            A(I,J)=0.d0
            be(I,J)=0.d0
            c(I,J)=0.d0
         enddo
      enddo
       E(1)=DBLE(E1)
       E(2)=DBLE(E2)
       E(3)=DBLE(E3)
       E(4)=DBLE(E4)
       E(5)=DBLE(E5)
       G(1)=DBLE(G1)
       G(2)=DBLE(G2)
       G(3)=DBLE(G3)
       G(4)=DBLE(G4)
       G(5)=DBLE(G5)
      DO I=1,N
         E(I)=E(I)/8065.46
      enddo
      DO 5395 I=1,N
      DO 5394 J=1,N
      WL(I,J)=0.0
      IF(I.EQ.J) GOTO 5394
      IF(E(I).EQ.E(J)) GOTO 5394
      WL(I,J)=DABS(12398.54/(E(I)-E(J)))
5394  CONTINUE
5395  CONTINUE
      DO 7352 I=1,N
      DO 7352 J=1,N
      C(I,J)=0.
 7352 A(I,J)=0.
      A(2,1)=DBLE(A21)
      A(3,1)=DBLE(A31)
      A(4,1)=DBLE(A41)
      A(5,1)=DBLE(A51)
      A(3,2)=DBLE(A32)
      A(4,2)=DBLE(A42)
      A(5,2)=DBLE(A52)
      A(4,3)=DBLE(A43)
      A(5,3)=DBLE(A53)
      A(5,4)=DBLE(A54)
      C(2,1)=DBLE(C21)
      C(3,1)=DBLE(C31)
      C(4,1)=DBLE(C41)
      C(5,1)=DBLE(C51)
      C(3,2)=DBLE(C32)
      C(4,2)=DBLE(C42)
      C(5,2)=DBLE(C52)
      C(4,3)=DBLE(C43)
      C(5,3)=DBLE(C53)
      C(5,4)=DBLE(C54)
      TEV=TE/1.1609D4
      TIME=8.64E4*TDAY
      DO 5 I=1,NM1
      IP1=I+1
      DO 5 J=IP1,N
      C(J,I)=8.63D-6*DENEL*C(J,I)/(DSQRT(TE)*G(J))
 5    C(I,J)=C(J,I)*DEXP(-(E(J)-E(I))/TEV)*G(J)/G(I)
      do k=1,5
         if(k==1) then
            x(1)=1.
         else
            x(k)=0.
         endif
         xold(k)=x(k)
      enddo
      DO L=1,20
      DO 100 I=1,N
      DO 100 J=1,N
      AA(I,J)=0.
      IF(I.EQ.N) AA(I,J)=1.D0
      IF(I.EQ.N) GOTO 100
      DEI=X(I)*Z*DEN
C
C     CALCULATE THE ESCAPE PROBABILITY.
C
      IF(ISTAT.EQ.1) THEN
        T0=C33*COLTOT(KI)*A(J,I)*G(J)*WL(J,I)**3./G(I)
        T0T=(COLTOTT(KI)-COLTOT(KI))*T0/COLTOTT(KI)
        IF(T0T.LE.0.) T0T = 0.
      ELSE
        T0=1.e-24*WL(J,I)**3*A(J,I)*G(J)*DEI*TIME/(8.*PI*G(I))
        T0T=0.
      ENDIF
      IF(L.EQ.1) T0=.0
C!!   ASSUME MEAN ATOMIC WEIGHT = 20. 
      AU=20.
      VTERM=1.285E6*SQRT(TE/(AU*1.E4))
      CALL ESCAPE(T0,T0T,WL(J,I),A(J,I),VTERM,BE(J,I),DBEDTA)
      topt(j,i)=t0
c      write(6,*)i,j,t0,t0t,be(j,i)
      IF(I.NE.J) GOTO 2
      S=0.
      DO  K=1,N
         S=S+C(I,K)+BE(I,K)*A(I,K)
      enddo
      AA(I,I)=-S
      GOTO 100
2     AA(I,J)=BE(J,I)*A(J,I)+C(J,I)
 100  CONTINUE
      NRC=N+1
      AA(N,NRC)=1.
      EPS=1.D-30
      NA=N

      DI=SIMUL(NA,AA,X,EPS,1,NRC)

      errmax=0.
      do k=1,n
         deltax=abs(x(k)-xold(k))/x(k)
         if(deltax>errmax) then
            errmax=deltax
         endif
      enddo
      if(errmax < 1.e-5) goto 44
      do k=1,5
         xold(k)=x(k)
      enddo
      ENDDO
 44   continue
      DO I=1,N
        XI(I)=X(I)
      ENDDO
      K=0

      DO I=2,N         
         DO J=1,I-1
            WOBS(I,J)=1.602E-12*Z*(E(I)-E(J))*X(I)*A(I,J)*BE(I,J)/DEN
            EM(I,J)=WOBS(I,J)
            IF(K.le.400.and.a(i,j).gt.0.) then
               K=K+1
               SR(K)=0.
               WEQ(K)=WOBS(I,J)
               wlin(k)=wl(i,j)
               taul(k)=topt(i,j)
               kmaxp(iel,ion)=k
            endif
         enddo
         kmaxp(iel,ion)= k
      enddo
      RETURN
      END

      SUBROUTINE FORB3(IJ,K,E21,E32,A21,A31,A32,O21,O31,O32,S1,S2,S3,
     &F21,F31,F32,TE,Z,F)
C     COOLING/VOL= flux/vol = F21*X(EL)*AB(I)*X(I)*DEN**2 etc
c     ij = index for ion 1 - 18 now. only for taufb()
c     k = ion in CF system used for column density
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
c      PARAMETER (MD=300,MDP1=MD+1)
      COMMON/TAUFBO/TAFB(20,3)
      COMMON/PHY/DENS(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/IND/I
      common/wlforb3/wl21,wl31,wl32
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DEN=DENS(I)
      XEL=DEL(I)
      E31=E21+E32
      C21=8.63E-6*DEN*XEL*O21/(S2*SQRT(TE))
      C31=8.63E-6*DEN*XEL*O31/(S3*SQRT(TE))
      C32=8.63E-6*DEN*XEL*O32/(S3*SQRT(TE))
      TEV=TE/1.1609E4
      C12=0.
      C13=0.
      C23=0.
      IF(E21/TEV.LT.100.) C12=S2*C21*EXP(-E21*1.1609E4/TE)/S1
      IF(E31/TEV.LT.100.) C13=S3*C31*EXP(-E31*1.1609E4/TE)/S1
      IF(E32/TEV.LT.100.) C23=S3*C32*EXP(-E32*1.1609E4/TE)/S2
      WL21=12398.*1.E-8/E21
      WL31=12398.*1.E-8/E31
      WL32=12398.*1.E-8/E32
      DENI=Z*DEns(i)
      DN1=1.
      DN2=1.e-10
      DN3=1.e-10
      mw=0
98    DE1=DN1*DENI
      DE2=DN2*DENI
      DE3=DN3*DENI
      mw=mw+1
      DN1O=DN1
      DN2O=DN2
      DN3O=DN3
      TIME=8.64E4*TDAY
      CALL ESCFOR(K,WL21,A21,S1,S2,DE1,TIME,T0,BE21)
      t21=t0
      CALL ESCFOR(K,WL31,A31,S1,S3,DE1,TIME,T0,BE31)
      t31=t0
      CALL ESCFOR(K,WL32,A32,S2,S3,DE2,TIME,T0,BE32)
      t32=t0
      Q1=(C12+C13)*(A32*BE32+C32)+C12*(A31*BE31+C31)
      Q2=(A21*BE21+C21)*(A32*BE32+C32)+(A21*BE21+C21+C23)*(A31*BE31+C31)
      DN21=Q1/Q2
      DN31=(C12+C13-DN21*(A21*BE21+C21))/(A31*BE31+C31)
      DN1=1./(1.+DN21+DN31)
      DN2=DN1*DN21
      DN3=DN1*DN31
      if(dn1o.le.0.and.mw.gt.20)WRITE(6,92)mw,DN1,DN2,DN3
      if(dn2o.le.0..and.mw.gt.20)WRITE(6,92)mw,DN1,DN2,DN3
      if(dn3o.le.0..and.mw.gt.20)WRITE(6,92)mw,wl31,DN1,DN2,DN3
      dr1=0.
      dr2=0.
      dr3=0.
      if(dn1o.gt.0.)DR1=ABS(DN1O-DN1)/DN1O
      if(dn2o.gt.0.)DR2=ABS(DN2O-DN2)/DN2O
      if(dn3o.gt.0.)DR3=ABS(DN3O-DN3)/DN3O
      DNMAX=MAX(DR1,DR2,DR3)
92    FORMAT(' DN ',i5,6E11.4)
      if(mw.gt.20) then
         write(6,*)' no conv in forb',k
         write(6,922)wl21*1.e8,E21,E32,A21,A31,A32,O21,O31,O32,S1,S2,S3,z
 922     format(' i forb ',f12.2,1pe12.3,20e12.3)
      endif
      if(mw.gt.20) goto 231
      IF(DNMAX.GT.0.1) GOTO 98
231   IF(E21/TEV.LT.100.) B21=S1*EXP(E21/TEV)*DN21/S2
      IF(E31/TEV.LT.100.) B31=S1*EXP(E31/TEV)*DN31/S3
      F21=ELCH*DN1*DN21*A21*BE21*E21/(XEL*DEN)
      F31=ELCH*DN1*DN31*A31*BE31*E31/(XEL*DEN)
      F32=ELCH*DN1*DN31*A32*BE32*E32/(XEL*DEN)
      F=F21+F31+F32
      tafb(ij,1)=t21
      tafb(ij,2)=t31
      tafb(ij,3)=t32
      WL21=WL21/1.E-8
      WL31=WL31/1.E-8
      WL32=WL32/1.E-8
      RETURN
      END
      

      SUBROUTINE ESCFOR(K,WL,A21,G1,G2,DE1,TIME,T0,BE)
C
C     CALCULATE THE ESCAPE PROBABILITY.
C
      IMPLICIT REAL*8(A-H,O-Z)
c     PARAMETER (NFEL=1000)
      include "parameters.h"
      COMMON/COLD/OPT(nio),COLTOT(nio),COLTOTT(nio),SRED(NFEL+100)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/A7/C3,C33
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      BE=1.
      IF(ISTAT.EQ.1) THEN
        WL21=WL/1.E-8
        T0=C33*COLTOT(K)*A21*G2*WL21**3./G1
        T0T=(COLTOTT(K)-COLTOT(K))*T0/COLTOTT(K)
        IF(T0T.LE.0.) T0T = 0.
      ELSE
        T0=WL**3*A21*G2*DE1*TIME/(8.*PI*G1)
        T0T=0.
      ENDIF
C!!   ASSUME MEAN ATOMIC WEIGHT = 20. 
      AU=20.
      VTERM=1.285E6*SQRT(TE/(AU*1.E4))
      CALL ESCAPE(T0,T0T,WL,A21,VTERM,BE,DBEDTA)
      RETURN
      END

      SUBROUTINE FORB(E21,E32,A21,A31,A32,O21,O31,O32,S1,S2,S3,
     &F21,F31,F32,TE,F)
C     COOLING/VOL= F*X(EL)*AB(I)*X(I)*DEN**2
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
c     PARAMETER (MD=300,MDP1=MD+1)
      COMMON/PHY/DENS(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/IND/I
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DEN=DENS(I)
      XEL=DEL(I)
      E31=E21+E32
      C21=8.63E-6*DEN*XEL*O21/(S2*SQRT(TE))
      C31=8.63E-6*DEN*XEL*O31/(S3*SQRT(TE))
      C32=8.63E-6*DEN*XEL*O32/(S3*SQRT(TE))
      TEV=TE/1.1609E4
      C12=0.
      C13=0.
      C23=0.
      IF(E21/TEV.LT.100.) C12=S2*C21*EXP(-E21*1.1609E4/TE)/S1
      IF(E31/TEV.LT.100.) C13=S3*C31*EXP(-E31*1.1609E4/TE)/S1
      IF(E32/TEV.LT.100.) C23=S3*C32*EXP(-E32*1.1609E4/TE)/S2
      Q1=(C12+C13)*(A32+C32)+C12*(A31+C31)
      Q2=(A21+C21)*(A32+C32)+(A21+C21+C23)*(A31+C31)
      DN21=Q1/Q2
      DN31=(C12+C13-DN21*(A21+C21))/(A31+C31)
      DN1=1./(1.+DN21+DN31)
      DN2=DN1*DN21
      DN3=DN1*DN31
      IF(E21/TEV.LT.100.) B21=S1*EXP(E21/TEV)*DN21/S2
      IF(E31/TEV.LT.100.) B31=S1*EXP(E31/TEV)*DN31/S3
      F21=ELCH*DN1*DN21*A21*E21/(XEL*DEN)
      F31=ELCH*DN1*DN31*A31*E31/(XEL*DEN)
      F32=ELCH*DN1*DN31*A32*E32/(XEL*DEN)
      F=F21+F31+F32
      RETURN
      END


      SUBROUTINE AUGION(iz,IZMAX,xel,XA)
c     iz = Z (26 for Fe). izmax = max ionization stage.
c      only for aug_fr where numbering is 1-14. translates from z to CF by z2cf(1...26)
      
      IMPLICIT REAL*8(A-H,O-Z)
c     PARAMETER (MD=300,MDP1=MD+1)
      include "parameters.h"
c      COMMON/ABC/AL(7),ALN(8),ALO(9),ALC(7),ALS(15),ALMG(4),ALAL(6)
c     &,COO(8),COC(6),ALCA(4),ALNA(4),COCA(4),CONA(4),ALSUL(11),ALNEO(9)
c     &,CONE(8) 
c      COMMON/COSIL/COS(18)
      parameter (mh=2,mhe=3,mc=7,mn=8,mo=9,mne=11,mna=12,mmg=13,
     &     mal=14,msi=15,ms=17,mar=19,mca=21,mfe=27,mni=29,
     &     mmat=mh+mhe+mc+mn+mo+mne+mna+mmg+mal+msi+ms+mar+mca+mfe+mni,
     &     ispecy=15)
      common/abc1/alc(mc),aln(mn),alo(mo),alne(mne),alna(mna),
     &     almg(mmg),alal(mal),alsi(msi),alsu(ms),alar(mar),alca(mca),
     &     alfe(mfe),alni(mni)
      common/abc2/coh(1),cohe(2),coc(mc-1),con(mn-1),coo(mo-1),
     &     cone(mne-1),cona(mna-1),comg(mmg-1),coal(mal-1),cosi(msi-1),
     &     cosu(ms-1),coar(mar-1),coca(mca-1),cofe(mfe-1),coni(mni-1)      
      COMMON/SECX/CSEC(20),DCSDX(20),CISEC(20),DCISDX
      COMMON/PHY/DEN(MD)
      COMMON/IND/IK
      COMMON/RSI/ZS(18,5)
      common/raug/zion(30,27,7)
      COMMON/AUG/AUG
      common/abl/abn(15)
      common/ionx/xion(md,14,27)
      common/ionx_old/xion_old(14,27)
      common/populations/popul(14,27,400)
      common/ctrates/chtr(14,14)
      parameter(nz=30,nion=27,nshell=10)
      integer iz,ion,shell,epsi,fr(10)      
      integer ns,kmax,n
      real*8 fr_aug,eion,en_aug,en_augi,eioni
      common/augerfrac/eion(nz,nion,nshell),en_aug(nz,nion,nshell),
     &     fr_aug(nz,nion,nshell,10),kmax(nz,nion,nshell),ns(nz,nion),
     &     init_augfrac
      common/rec_coll/alrec(30,30),collion(30,30)
      common/init_ion/inition
c      PARAMETER (NL=130,NLP1=NL+1)
      real*8 x(30),dx(30),aa(31,31),ztot(30)
      DIMENSION XA(30),ZSA(30,18,5),phrate(30,30)
c       dimension nionel(14)
c       data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      integer z2cf(26)      
      data z2cf/1,2,0,0,0,3,4,5,0,6,7,8,9,10,0,11,0,12,0,13,0,0,0,0,0,14/
C     IONIZATION BALANCE EQUATIONS TAKING INTO ACCOUNT AUGER
C      IONIZATION AND COLLISIONAL IONIZATION (FROM WEISHEIT ).
C     XR(I)=N(I+1)/N(I)
c      CALL COSI

c     set up ionization balance
c      iz=14
c     izmax=6
      im1=0
      im2=0
      nionel=iz

      do i=1,30
         x(i)=0.
      enddo
      do i=1,30
         xa(i)=0.
         do j=1,30
            aa(i,j)=0.
         enddo
      enddo

      do ion=1,izmax
         ztot(ion)=0.
         ztot2=0.
         do i=1,ns(z2cf(iz),ion)
            ztot2= zion(iz,ion,i) + ztot2
            ztot(ion)=ztot2
         enddo
      enddo

      do i=1,izmax

         aa(i,izmax+1)=0.

         do k=1,10
            il=i-k
            if(il>0 ) then
               do nsi=1,ns(z2cf(iz),il)
                  if(kmax(z2cf(iz),il,nsi)>0) then                     
                     aa(i,il)=fr_aug(z2cf(iz),il,nsi,k)*zion(iz,il,nsi) +
     &                    aa(i,il)
                     if(il==-1) then
                        aa(i,i-1)=aa(i,i-1) + collion(iz,i-1)*xel*den(ik)
                     endif
                  endif
               enddo
            endif
         enddo

         aa(i,i+1)=alrec(iz,i+1)*xel*den(ik)
         
         aa(i,i)=-alrec(iz,i)*xel*den(ik)-ztot(i) - collion(iz,i)*xel*den(ik)         
c     &              ction(iz,i)*dent - collion(iz,i)*xel*den(ik)         
      enddo

      icon=izmax
c     replace highest stage with number cons
      do i=1,izmax
c         aa(izmax,i)=1.d0
         aa(icon,i)=1.d0
      enddo

      aa(icon,izmax+1)=1.
c      aa(izmax,izmax+1)=1.

c     if(izmax >= 8 .and. inition==0 ) then
      if(izmax >= 8 .and. ik>=3 ) then

         eps=1.e-6
         eps=0.

         im1=0
         do i=1,izmax
            if(xion_old(z2cf(iz),i) < eps) then
               im1=1
            endif
         enddo

         im2=0
         do i=izmax,1,-1
            if(xion_old(z2cf(iz),i) < eps) then
               im2=1
            endif
         enddo
         im1=0
         im2=0

         
         
         izmax1=izmax-im1
         do i=1,izmax	
            i1=i-im1	
            do j=1,izmax+1
               j1=j-im1
               if(i1 > 0.and.j1 > 0) then
                  aa(i1,j1)=aa(i,j)
               endif
            enddo
         enddo

         im2=2
         izmax2=izmax1-im2
         do j=1,izmax2+1
            aa(izmax2,j)=1.
         enddo
      else
         izmax2=izmax
      endif

      di=simul_ion(izmax2,aa,x,eps,1,nrc)

      dxmax =0.

      do i=1,iz+1
         if(i <= im1) then
            xa(i)=0.
         elseif(i > izmax-im2) then
            xa(i)=0.
         else
            xa(i)=x(i-im1)
         endif
         xion(ik,z2cf(iz),i)=xa(i)
      enddo
      RETURN
      END


      subroutine auger(izmax,deel,xa)
      implicit real*8(a-h,o-z)
      include "parameters.h"
      parameter (mh=2,mhe=3,mc=7,mn=8,mo=9,mne=11,mna=12,mmg=13,
     &     mal=14,msi=15,ms=17,mar=19,mca=21,mfe=27,mni=29,
     &     mmat=mh+mhe+mc+mn+mo+mne+mna+mmg+mal+msi+ms+mar+mca+mfe+mni,
     &     ispecy=15)
      common/abc1/alc(mc),aln(mn),alo(mo),alne(mne),alna(mna),
     &     almg(mmg),alal(mal),alsi(msi),alsu(ms),alar(mar),alca(mca),
     &     alfe(mfe),alni(mni)
      common/abc2/coh(1),cohe(2),coc(mc-1),con(mn-1),coo(mo-1),
     &     cone(mne-1),cona(mna-1),comg(mmg-1),coal(mal-1),cosi(msi-1),
     &     cosu(ms-1),coar(mar-1),coca(mca-1),cofe(mfe-1),coni(mni-1)      

      common/phy/den(md)
      common/ind/ik
      common/rar/zar(18,5)
      common/augalu/alaug(27),coaug(26)
      common/aug/aug
      dimension xa(30),zsa(18,5)
c     ionization balance equations taking into account auger
c      ionization and collisional ionization
c     xr(i)=n(i+1)/n(i)
      do  iz=1,izmax
         do  mi=1,5
            zsa(iz,mi)=zar(iz,mi)
         enddo
      enddo
      do iz=1,izmax
         ze=0.
         do mi=1,5
            ze=ze+zar(iz,mi)
         enddo
         if(iz.le.1) then
            zeff=ze
            xa(1)=(coaug(1)+ze/(den(ik)*deel))/alaug(iz+1)
         else
            xs=1.
            ia1=1
            ia2=iz-1
            if(iz.gt.6) then
               ia1=iz-5
            endif
c     ia=i (weisheit)
c     i=z-i (weisheit)
            s=0.
            do ia=ia2,ia1,-1
               i=iz-ia
               if(xa(ia).le.0.) xa(ia)=1.e-20
               xs=xs/xa(ia)
               n=18+1-ia
               call psi(ia,i,n,zsa,ps)
               s=s+xs*ps
            enddo
            zeff=ze+zeff/xa(iz-1)-s
            if(aug.lt.0.) zeff=ze
            xa(iz)=(coaug(iz)+zeff/(deel*den(ik)))/alaug(iz+1)
         endif
      enddo
      return
      end


      SUBROUTINE PSI(IZ,IA,N,ZS,PS)
      IMPLICIT REAL*8(A-H,O-Z)
C     IZ=LOWER INDEX Z
C     IA=UPPER INDEX (1...5)
C     N= C OF EL. LEFT
C     ZS=IONIZATION RATES
C     PS=PSI
      DIMENSION ZS(18,5)
C     ZS(IZ,I), I=1 = K-SHELL, 2= 2S, 3= 2P, 4= 3S. (total M-shell)
      PS=0.
      IF(IA.eq.1) then
         PS=ZS(IZ,5)
         IF(N.le.10) then
c total L skell = 2s+2p            
            PS=PS+ZS(IZ,2)+ZS(IZ,3)
         else
            PS=PS+ZS(IZ,3)
         endif
      elseIF(IA.eq.2) then
         IF(N.ge.12) then
            PS=ZS(IZ,3)
         elseIF(N.eq.12) then
            PS=ZS(IZ,2)
         elseIF(N.eq.11) then
            PS=zs(iz,2)+2.*ZS(IZ,1)/3.
         elseIF(N.le.10) then
            PS=ZS(IZ,1)
         endif
      elseIF(IA.eq.3) then
         IF(N.ge.13) then
            PS=ZS(IZ,2)
         elseIF(N.eq.13) then
            PS=ZS(IZ,2)+0.666667*ZS(IZ,1)
         elseIF(N.eq.12) then
            PS=ZS(IZ,1)
         elseIF(N.eq.11) then
            PS=0.3333333*ZS(IZ,1)
         endif
      elseIF(IA.NE.4) then
         IF(N.eq.13) then
            PS=ZS(IZ,1)/3.
         elseIF(N.eq.14) then
            PS=ZS(IZ,1)
         elseIF(N.ge.15) then
            PS=2.*ZS(IZ,1)/3.
         endif
      elseIF(IA.NE.5) then
         IF(N.ge.15) then
            PS=ZS(IZ,1)/3.
         endif
      endif
      RETURN
      END
      
      SUBROUTINE RATEAUG(iel,izmax)
c Note iel = Z    Fe = 26
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/RSI/ZS(18,5)
      common/raug/zion(30,27,7)
      common/populations/popul(14,27,400)
      COMMON/IND/IK
      real*8 jmean
      COMMON/DTAU/jmean(NE1:NE2)
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/CSI_vern2/GSV2(30,30,7,NE1:NE2)
      common/shells/ns(30,27)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/


c ionization rate. Note that flux is NOT flux but mean intensity J!      
      
      ZA=0.
      DO IZ=1,IZMAX
         DO  IS=1,ns(iel,iz)
            Zion(iel,IZ,IS)=0.
            DO J=JMIN,JJ
               ZA=4.*PI*jmean(J)*GSv2(iel,IZ,IS,J)*(E(J+1)-E(J))/
     &              (ELCH*E1(J))         
               ZION(iel,IZ,IS)=ZION(iel,IZ,IS)+ZA
            enddo
         enddo
         if(iel==14.and.iz==1) then
            DO J=JMIN,JJ
c     photoionization of excited 1D state of Si I (Omand & AJ 2023)
c     add all to the outer shell is=5
               enex=8.15168-0.7810
               if(enex < e1(j) .and. e1(j) < 9.6636) then                     
                  sigex=2.93e-17*(e1(j)/enex)**4.07
               elseif(9.6636 < e1(j) .and. e1(j) < 18.1) then
                  sigex=3.47e-16*(enex/e1(j))**5.48
               elseif(e1(j) > 18.1) then
                  sigex=3.60e-17*(enex/e1(j))**3
               else
                  sigex=0.
               endif
               is=5
               iex=4
               ielcf=10
               xex=popul(ielcf,iz,iex)
               Zex=4.*PI*jmean(J)*xex*sigex*(E(J+1)-E(J))/
     &              (ELCH*E1(J))         
               ZION(iel,IZ,IS)=ZION(iel,IZ,IS)+Zex
            enddo
         endif
      enddo
      RETURN
      END

      SUBROUTINE COLLI(S,Z,E0,CO)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 S,Z,E0
      COMMON/TEQQ/TE
      X=1.578E5*(E0/13.6)/TE
      CO=0.
      F=3.1-1.2/(Z+1.)-.9/(Z+1)**2
      IF(X.LT.100.)CO=2.223E-8*S*F*SQRT(TE)*EXP(-X)*(1.-EXP(-X))
     &/E0**2
      RETURN
      END



