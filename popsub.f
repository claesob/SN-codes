      SUBROUTINE MULTISIMPQ(iel,ion,TE,XI,rltot)
c
c Obs! Keep total number of 'neutral' and 'ionized' atoms constant
c
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N
      include "parameters.h"
      COMMON/NION/IONQ
      COMMON/COLDnew/COLTOTnew(14,27),COLTOTTnew(14,27)
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
      integer nlf1,nlf2,nlf3,nlf4,nlsi1
      parameter(nlf1=122,nlf2=192,nlf3=113,nlf4=46)
      parameter(nlsi1=66)
      common/ferecombi/recof1(nlf1),recof2(nlf2)
      common/sirecombi/recosi1(nlsi1)      
      common/w/wlin(501),wl(NL,NL),taul(501)
      common/line_em/cinx(14,26,401),taulinex(14,26,401),wlix(14,26,401)
     &     ,ilabcfx(14,26,401)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      common/ionx/xion(md,14,27)
      common/abl/abn(15)
      common/populations/popul(14,27,400)
      common/preion/ipre
      common/kmaxpop/kmax(14,27)
      common/timecheck/time,itime
      common/velgrad/dr_dv      
      common/debug/ideb
      common/critdens/dcrit(14,27,401)
      real*8 XI(NL),C(NL,NL)
      real*8 AA(nlp1,nlp1),X(nl),XO(NL),betot(nl,nl)
      dimension kch(nl,nl),tau(nl,nl)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/

      ipr=0

      if(iel.eq.4.and.ion.eq.2) then
         ipr=1
      elseif(iel.eq.5.and.ion.eq.2) then
         ipr=1
      endif
      if(iel==12.and.ion==5) then
         ipr=1
      endif
      if(iel==14.and.ion==7) then
         ipr=1
      endif

      DEN=DENS(IK)
      DENEL=DEN*DEL(IK)
      TEV=TE/1.1609E4
      x0 = xion(ik,iel,ion)
      xplus = xion(ik,iel,ion+1)
      z=abn(iel)
      if(iel==14.and.ion==1) then
         do ii=1,nlf1
            if(ii <= nlev) then
               reco(ii)=recof1(ii)
            endif
         enddo
      elseif(iel==14.and.ion==2) then
         do ii=1,nlf2
            if(ii <= nlev) then
               reco(ii)=recof2(ii)
            endif
         enddo
      elseif(iel==10.and.ion==1) then
         trec=0.
         do ii=1,nlsi1-1
            reco(ii)=recosi1(ii)
            trec=trec+reco(ii)
         enddo
      endif

      if(ik.gt.2.and.(iel.eq.1.or.iel.eq.2)) then
         do i=1,nlev
            x(i)=xn(ionq,i)
         enddo
      endif
      do i=1,nl
         x(i)=1.e-10
      enddo
      x(1)=0.999
      do i=1,500
         weq(i)=0.
         sr(i)=0.
         wlin(i)=0.
         taul(i)=0.
      enddo

      DO I=1,NLP1
         DO J=1,NLP1
            AA(I,J)=0.
         ENDDO
      ENDDO

      DO I=1,NLev
         DO J=1,NLev
            c(I,J)=0.
c            be(I,J)=0.
            betot(I,J)=0.
         ENDDO
      ENDDO
C     ADD ONE LEVEL FOR CONTINUUM

      N=NLEV+1
      E(N)=E00
      DO I=1,N
         ph(i)=0.
         DO J=1,Nlev
            
            WL(I,J)=0.0
            IF(I.NE.J.AND.E(I).NE.E(J)) THEN
               WL(I,J)=DABS(12398.54d0/(E(I)-E(J)))
            ENDIF
         ENDDO
      ENDDO
      TEV=TE/1.1609D4
      TIME=8.64E4*TDAY
      x(1)=x0
      DO I=1,NLEV
        DO J=1,NLEV
           C(J,I)=DENEL*CS(J,I)
       ENDDO
      ENDDO

      dO L=1,20
         DO I=1,Nlev
            DEI=X(I)*Z*DEN
            DO J=1,Nlev
               AA(I,J)=0.
               BEtot(I,I)=0.
               IF(I.NE.J) THEN
C     
C     CALCULATE THE ESCAPE PROBABILITY.
C     
                  if(a(j,i) > 0.) then
                     IF(ISTAT.EQ.1) THEN
                        T0=C33*COLTOTnew(iel,ion)*X(I)*A(J,I)*G(J)*     
     &                       WL(J,I)**3./G(I)
                        T0T=T0
                     ELSE
                        T0=1.e-24*WL(J,I)**3*A(J,I)*G(J)*DEI*TIME/
     &                       (8.*PI*G(I))                        
                     ENDIF
                     tau(j,i)=t0
                  else
                     t0=0.
                  endif
                  AU=4.
                  if(iel.eq.1) then
                     AU=1.
                  elseif(iel.eq.2) then
                     au=4.
                  elseif(iel.eq.3) then
                     au=12.
                  elseif(iel.eq.4) then
                     au=14.
                  elseif(iel.eq.5) then
                     au=16.
                  elseif(iel.eq.10) then
                     au=28.
                  elseif(iel.eq.14) then
                     au=56.
                  elseif(iel.eq.12) then
                     au=40.
                  endif
                  VTERM=1.285E6*SQRT(TE/(AU*1.E4))
                  CALL ESCAPE(T0,T0T,WL(J,I),A(J,I),VTERM,
     &                 BE_ji,DBEDTA)

                  ESC(J,I)=BE_ji
                  
                  tau(j,i) = t0
c$$$                  if((iel==5.or.iel==11).and.t0>0.) then
c$$$                     write(6,*)'tau1 ',iel,ion,i,j,wl(i,j),t0,tau(j,i)
c$$$                  endif
                  
                  if(j.gt.i) then 
                     be_ij=be_ji
                  endif
                  besc=be_ji
                  pd=0.
                  ll=0
                  if(ionq.eq.3.or.ionq.eq.5) then
c     first 584 A line
                     if(ionq.eq.3) then
                        if(i.eq.5.and.j.eq.1) then
                           ll=47
                        endif
                     endif
c     now Lyman alpha
                     if(ionq.eq.5) then
                        if(i.eq.2.and.j.eq.1) then
                           ll=63
                        endif
                     endif
                     if(ll.ne.0) then
c     line center opacity
c     line center opacity
                        opline=3.979E-26*wl(i,j)**3.*g(i)*a(i,j)*
     &                       dei/(g(j)*vterm) 
                        opline2=t0/(vterm*time*den)
c     continuum desruction probability. NOT INCLUDED!!
c     opcon=totopl(ll)
                        opcon=0.
c     H-R expression
                        fhumm=8.5
                        if(opline.gt.0.) then
                           betadest=opcon/opline
                        else
                           betadest=1.d33
                        endif
                        fbeta=min(fhumm*betadest,1.d-3)
                        pd=fbeta/(fbeta+1.)
                     endif
                  endif
                  betot(j,i)=(1.-pd)*besc+pd
               ENDIF
               IF(I.EQ.J) THEN
                  S=0.
                  DO K=1,Nlev
                     S=S+C(I,K)+BEtot(I,K)*A(I,K)
                  ENDDO
C     ADD PHOTOIONIZATION
                  AA(I,I)=-S-PH(I)
               ELSEIF(J.LT.N) THEN
                  AA(I,J)=BEtot(J,I)*A(J,I)+C(J,I)
               ELSEIF(J.EQ.N) THEN
C     RECOMBINATION CONTRIBUTION
                  AA(I,N)=-RECO(I)*DENEL*xplus
               ENDIF     
            ENDDO
         ENDDO

c     replace i=1 eqn with NUMBER CONSERVATION
c     replace first level eqn. by number cons.
         DO J=1,Nlev
            AA(1,J)=1.D0
         enddo

         AA(1,N)=x0

         EPS=1.D-30
c     DO I=1,nlev
         DO I=1,n
            XO(I)=X(I)
         ENDDO
c!!   

         neq=n-1
         nlhs=n

         do i=1,n
            x(i)=0.
         enddo
         DI=SIMUL(Neq,AA,X,EPS,1,Nlhs)
         ERRMAX=0.
         x(n)=xplus
         DO I=1,N
            IF(X(I).NE.0.) THEN
               ERR=ABS((XO(I)-X(I))/X(I))
            ELSE
               ERR=1.
            ENDIF
            IF(ERR.GT.ERRMAX) ERRMAX=ERR
         ENDDO

         IF(ERRMAX.LT.0.01.and.l.ge.2) GOTO 555
      ENDDO
 555  CONTINUE
      DO I=1,N
         XI(I)=X(I)
         popul(iel,ion,i)=x(i)
      ENDDO
      K=0

      rltot=0.
      wtot=0.
      ierr=0

      DO I=2,NLEV
         DO J=1,i-1
c            kch(i,j)=0
            WOBS(I,J)=1.602E-12*Z*(E(I)-E(J))*X(I)*A(I,J)*BEtot(I,J)/DEN
            EM(I,J)=1.602E-12*Z*(E(I)-E(J))*
     &           (X(J)*C(J,I)-X(I)*C(I,J))/DENEL
            rltot=rltot+em(i,j)
            wtot=wtot+wobs(i,j)

            IF(K.LE.400.and.a(i,j).gt.0.) THEN
               if(iel.eq.1.and.(j.eq.2.or.j.eq.3)) then
                  ikk=1
               elseif(iel.eq.1.and.i.le.20) then
                  ikk=1
               elseif(iel.ne.1) then
                  ikk=1
               else
                  ikk=0
               endif
               if(ikk.eq.1) then
                  K=K+1
                  SR(K)=0.
                  WEQ(K)=WOBS(I,J)
                  wlin(k)=wl(i,j)
                  sig=1.e4/wlin(k)
                  refn=1.+1.e-8*(6431.8+2949330./(146.-sig**2)+
     &                 25536./(41.-sig**2))
                  wlair=wl(i,j)/refn
                  if(wl(i,j).gt.2000.) then
                     wlin(k)=wlair
                  endif
                  wlix(iel,ion,k)=wlin(k)
                  taulinex(iel,ion,k)=tau(i,j)
                  dcrit(iel,ion,k)=a(i,j)*denel/c(i,j)
c                  if((iel==5.or.iel==11).and.tau(i,j)>0.) then
c                     write(6,*)'tau2 ',iel,ion,i,j,k,wl(i,j),tau(i,j),tau(j,i)
c                  endif
                  if(weq(k).lt.-1.d-30) then
                     ierr=1
                  endif

                  if(weq(k).lt.-1.d-30.or.weq(k).gt.1.d-10.or.
     &                 (iel.eq.-13.and.ion.eq.2)) then
                     write(6,9229)k,i,j,wl(i,j),x(i),e(i),e(j),te,
     &                    xplus,a(i,j),betot(i,j),c(i,j),em(i,j),del(ik),weq(k)
                  endif
 9228          format('weq ',3i5,f12.2,1pe12.3,10e12.3)
 9229             format('weq<0 ',3i5,f12.2,1pe12.3,10e12.3)
               endif
            endif
         ENDDO
      ENDDO


      xi(n)=xplus
      kmax(iel,ion)= k-1
      
      if(errmax > 0.99) then
c         write(6,*)' No convergence ',iel,ion,te,errmax
         DO I=2,NLEV
            DO J=1,I-1
               em(i,j)=0.
               wobs(i,j)=0.
            enddo
         enddo
      endif
      return
      end

      SUBROUTINE COLLEX_sh(ION,TS)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
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

