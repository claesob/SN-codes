
      SUBROUTINE EMISS(TE)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 ar(100),BD(100),cr(100),dr(100),er(100)
c      PARAMETER (MD=300,MDP1=301)
c      PARAMETER(NFEL=1000)
      PARAMETER (MPAR=310,MPP1=311)
      include "parameters.h"
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
C     *****************************************************************
C     *****
C     THIS ROUTINE CALC. THE EMISSION RATES FOR GIVEN TEMPERATURE
C     AND ABUNDANCES OF THE RELEVANT IONS.FREE-FREE AND RECOMBI-
C     NATION EMISSION ARE CONSIDERED.
C     ALL RATES ARE EMISSION INTEGRATED OVER 4 PI, EXCEPT EM(I,J)
C     WHICH IS EMISSION PER STERADIAN
C     *****
C     *****************************************************************
      COMMON/EMHY/RECEM(5,NL),TWOPHH,TWOPHHEI,COHI,PHI,PHIT,PO,POT
      common/heiires/x2he,he2coll12,he2jla,he2j2g,he2jrec(4),
     &                                                he2jem(7,4)
      COMMON/GSREC/ALGS(nio),ALEX(nio),ALTOT(nio),RECEX(nio),RECGS(nio)
      COMMON/ELEC/DEL(MD)
      COMMON/ION/XB(MD,nio)
      COMMON/ABU/XA(2,nio)
      COMMON/IND/I
      COMMON/TRES/EL(nio),EK(nio),ELord(14,27),EKord(14,27)
      COMMON/EQUIV/RADF,FR(2,100),W(100),CIN(100),WLI(NFEL+100),FB(100)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/PHY/DEN(MD)
      COMMON/REC/AL2(7)
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &     BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTEQ(NL),PHN(NL),DRDT(NL),RECT(5,NL),PHET(NL)
      COMMON/ABUN/AB(20)
      COMMON/SPOT/OTSP(7)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/LINE/XLINE
      COMMON/SIK/SK(ncr,NE1:NE2)
      COMMON/REC1/RE1(7),RHE2B
      COMMON/WFE/WOBFE(NL,NL),IQPOPH
      COMMON/FELEV/NFEII
      COMMON/A14/CQ(NL,NL),CI(NL),G(NL),EI(NL),AQ(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      common/inidiel/initdi
      common/line_em/cinx(14,26,401),taulinex(14,26,401),wlix(14,26,401)
     &     ,ilabcfx(14,26,401)
      common/kmax/kmax(14,27)
      common/gsreco/al_gs(26,26,10)
      integer nionstage
      common/num_ion_stages/nionstage(14)
      common/abl/abn(15)
      common/ethresh/ethresh(30,30,7)
      common/ionx/xion(md,14,27)
      real*8 encont
      common/conten/encont(26,26,10)      
      DIMENSION GS(50),Z(50),KS(50),P(500),ENL(500)
      DIMENSION INDI(50),INDIC(10,50),JMAX(50),JIO(50)
      DIMENSION CHI(nio),X(nio),ALV(nio),EX(nio)
      DIMENSION EMCHECK(NL),WLF(54)
      DIMENSION WLR(8,5),FRR(8,5),XO(8)
      DIMENSION IOION(100),WLDR(100)
      save ar,BD,cr,dr,er,IOION,WLDR
c     recombination emission from O I-V
      DATA WLR/8*0.,
     &        673., 617., 518., 515., 485., 539., 430., 0.,
     &        435., 396., 345., 328., 321., 374., 306., 304., 
     &        279.8, 238.5, 6*0., 
     &        629.7,248.5,172.2,220.3,215.2,192.9,1216.,0./
      DATA FRR/8*0.,
     &            .026,.053,.014,.020,.038,.30,.035,0.,
     &            .014,.074,.013,.023,.032,.27,.069,.040,
     &            .22,.27,6*0.,
     &            .11,.0057,.023,.041,.11,.16,.37,0./
C     GS = (STATISTICAL WEIGHT OF GROUND STATE)/2.
      DATA GS/2.,0.5,2.,2.,.5,2.,.5,.5,2.,.5,2.,.83333,6.,1.25,
     &        .44444,1.5,6.,.44444,1.5,6.,.5,2.,.5,2.,26*0./
C     CHARGE OF THE RECOMBINING ION
      DATA Z/1.,1.,2.,6.,7.,8.,5.,3.,4.,5.,6.,1.,2.,1.,2.,3.,4.,
     &       1.,2.,3.,4.,5.,6.,7.,26*0./
C     IF KS=1 THEN USE K-SHELL CROSSECTION
      DATA KS/1,1,1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,
     &        -1,-1,-1,1,1,26*-1/
      DATA JIO/2,7,14,47*0/
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      dimension nion(14)
      data nion/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      ideb=1
c      write(6,*)'he2jrec(1),he2jrec(2) ',te,he2jrec(1),he2jrec(2)
      emtot=0.
      TEV=TE/1.1609E4
      T4=TE/1.E4
      DO  J=JMIN,JJ
         EM(I,J)=0.
         EMC(I,J)=0.
      ENDDO
      DO K=1,100
         XB(I,K)=XA(2,K)
      ENDDO
      CALL XIONS(X,CHI)
C      ALV(I) = RECOMBINATION COEFFICIENT TO GROUND STATE OF ION I
      DO K=1,100
         EX(K)=0.5
         ALV(K)=ALGS(K)*1.E13*T4**.5
      ENDDO
      ALV(1)=1.54
      EX(1)=0.51
C     HE I CONSISTENT WITH RECOMB()
      ALV(2)=1.59
      EX(2)=0.488
c     do He II separately below
      ALV(3)=0.
      EX(3)=0.503
      ALV(4)=8.98
      EX(4)=0.50
      ALV(7)=3.27
      EX(7)=0.50
      ALV(8)=1.30
      EX(8)=0.50
      ALV(9)=3.61
      EX(9)=0.50
      ALV(10)=28.1
      EX(10)=0.50
      ALV(11)=65.0
      EX(11)=0.50
      ALV(13)=13.35
      EX(13)=0.50
c!!   Same al as in RECOMB()
      ALV(14)=1.22
      EX(14)=0.70
      ALV(42)=2.46E-2
c!!   guesses for 42 ad 43 = Fe I and II
      ALV(42)=0.4*1.42
      ALV(43)=0.4*10.2
      ALV(50)=5.61E-4
      ALV(51)=4.50E-2
      ALV(53)=2.83E-3
      DO K=1,100
         ALV(K)=ALV(K)*1.E-13/T4**EX(K)
      enddo
      if(initdi.eq.1) then
         INITDI=0
         OPEN(39,FILE='./ATDAT/diel.dat')
         rewind (39)
         do k=1,100
            read(39,*,end=88)IOION(k),WLDR(k),ar(k),BD(k),cr(k),dr(k),er(k)
         enddo
 88      ndie=k-1      
      endif

      DO K=2,56
         IF(ALV(K).gt.0.) then
            DO J=JMIN,JJ
c put all the b-f emission into the first bin! ok if kt << chi
               IF(E(J).gt.CHI(K)) then
                  REV=ALV(K)*CHI(K)*ELCH*X(K)*DEL(I)/(4.*PI)
                  EM(I,J-1)=REV/(E(J)-E(J-1))+EM(I,J-1)
                  if(j-1==2.and.ideb==1) then
                     write(6,*)'em rev a0 ',j,k,rev,alv(k),em(i,j-1)
                  endif                   
                  if(em(i,j-1).lt.-1.d-40) write(6,9128)i,j,k,rev,em(i,j-1)
 9128             format(' rev ',3i3,10e12.4)
                  EMC(I,J-1)=REV/(E(J)-E(J-1))+EMC(I,J-1)
 1828             format(' rec ',2i4,9e10.3)
                  GOTO 120
               endif
            enddo
 120        CONTINUE
         endif
      enddo
      do ielcf=3,14
         iel=nion(ielcf)
         do ion=1,26
            if(al_gs(iel,ion,1) > 0.) then
               do n=1,10
                  do j=jmin,jj
                     if(e1(j) > encont(iel,ion,n) .and. encont(iel,ion,n) > 0.) then
                        reclum=al_gs(iel,ion,n)*encont(iel,ion,n)*elch*
     &                       del(i)*abn(ielcf)*xion(i,ielcf,ion)/(4.*pi)
                        em(i,j-1)=em(i,j-1) + reclum/(e(j)-e(j-1))
c                        write(6,9192)iel,ion,n,j,e1(j),encont(iel,ion,n)
c     &                       ,12398./encont(iel,ion,n),reclum,em(i,j-1)
 9192                   format('reclum ',4i5,1pe12.3,10e12.3)
                        goto 77
                     endif
                  enddo
 77               continue
               enddo
            endif
         enddo
      enddo

c$$$
c$$$      
c$$$C     RECOMBINATION EMISSION FROM O II - V
c$$$      DO L=2,5
c$$$        DO K=1,8
c$$$          IF(WLR(K,L).GT.0.) THEN
c$$$            EN=12398.54/WLR(K,L)
c$$$            DO J=JMIN,JJ
c$$$              IF(EN.GT.E(J-1).AND.EN.LT.E(J)) THEN
c$$$                REV=FRR(K,L)*ALEX(12+L)*EN*ELCH*X(12+L)*DEL(I)/(4.*PI)
c$$$                EM(I,J-1)=REV/(E(J)-E(J-1))+EM(I,J-1)
c$$$                if(j-1==0.and.ideb==1) then
c$$$                   write(6,*)'em rev aa ',j,k,rev,em(i,j-1)
c$$$                endif                   
c$$$                GOTO 750
c$$$              ENDIF
c$$$            ENDDO
c$$$          ENDIF
c$$$750       CONTINUE
c$$$        ENDDO
c$$$      ENDDO
c$$$c     excited continua of O III - O IV
c$$$      DO K=1,3
c$$$        IF(K.EQ.1) THEN
c$$$C          O III N=3
c$$$           REV=16.3*ELCH*AB(3)*XA(2,17)*1.11E-12/T4**.769
c$$$           EN=16.3
c$$$        ELSEIF(K.EQ.2) THEN
c$$$C          O IV N=4
c$$$           REV=14.1*ELCH*AB(3)*XA(2,7)*1.23E-12/T4**.843
c$$$           EN=14.10
c$$$        ELSEIF(K.EQ.3) THEN
c$$$C          O IV N=3
c$$$           REV=27.45*ELCH*AB(3)*XA(2,7)*1.53E-12/T4**.769
c$$$           EN=27.45
c$$$        ENDIF
c$$$           REV=REV*DEL(I)/(4.*PI)
c$$$           DO J=JMIN,JJ
c$$$              IF(EN.GT.E(J-1).AND.EN.LT.E(J)) THEN
c$$$                 EM(I,J-1)=REV/(E(J)-E(J-1))+EM(I,J-1)
c$$$                 if(j-1==2.and.ideb==1) then
c$$$                   write(6,*)'em rev a ',j,k,rev,cit,em(i,j-1)
c$$$                endif                   
c$$$                GOTO 752
c$$$              ENDIF
c$$$            ENDDO
c$$$752       CONTINUE
c$$$      ENDDO
c$$$      XO(2)=XA(2,16)
c$$$      XO(3)=XA(2,17)
c$$$      XO(4)=XA(2,7)
c$$$      XO(5)=XA(2,4)
c$$$      do k=1,ndie
c$$$         write(6,9219)k,ar(k),BD(k),cr(k),dr(k),er(k),T4
c$$$         write(0,9219)k,ar(k),BD(k),cr(k),dr(k),er(k),T4
c$$$ 9219    format('k,a,b,c,d,e,t4',i4,1pe11.3,10e11.3)
c$$$        CALL DIELB(ar(k),BD(k),cr(k),dr(k),er(k),T4,ALDB)
c$$$        cit=AB(3)*XO(IOION(k))*ALDB*1.9864E-8/WLDR(k)
c$$$c        write(6,9727)k,ioion(k),AB(3),XO(IOION(k)),ALDB,WLDR(k),cit
c$$$ 9727   format('dielb ',2i5,1pe11.3,10e11.3)
c$$$        EN=12398.54/WLDR(K)
c$$$            DO J=JMIN,JJ
c$$$              IF(EN.GT.E(J-1).AND.EN.LT.E(J)) THEN
c$$$                REV=CIT*DEL(I)/(4.*PI)
c$$$                EM(I,J-1)=REV/(E(J)-E(J-1))+EM(I,J-1)
c$$$                if(j-1==2.and.ideb==1) then
c$$$                   write(6,*)'em rev ',j-1,k,rev,cit,em(i,j-1)
c$$$                endif                   
c$$$                GOTO 751
c$$$              ENDIF
c$$$            ENDDO
c$$$751     CONTINUE
c$$$      enddo
c$$$C     FORBIDDEN LINES 
c$$$      CALL WLFORB(WLF)
c$$$      DO 101 K=1,45
c$$$C     SKIP FORBIDDEN O I LINES
c$$$      IF(K.GE.10.AND.K.LE.12) GOTO 101      
c$$$      IF(WLF(K).LE.0.) GOTO 101
c$$$      IF(WLF(K).GT.0.) EN=12398.54/WLF(K)
c$$$      IF(EN.LT.E(JMIN)) GOTO 101
c$$$            DO 201 J=JMIN,JJ
c$$$            IF(EN.GT.E(J)) GOTO 201
c$$$                  EMJ=DEL(I)*FB(K)/(4.*PI*(E(J)-E(J-1)))
c$$$                  EM(I,J-1)=EM(I,J-1)+EMJ
c$$$                  if(j-1==2.and.ideb==1) then
c$$$                     write(6,*)' em -4 to 6 ',j,k,fb(k),emj,em(i,j-1)
c$$$                  endif
c$$$                  JK=J-1
c$$$                  GOTO 101
c$$$201          CONTINUE
c$$$101   CONTINUE
c$$$C     LINES
c$$$      
      DO iel=3,14
         do ion=1,nionstage(iel)
            do k=1,401
               wl0=wlix(iel,ion,k)
               IF(wl0.gt.0.) then
                  EN=12398.54/WL0
                  IF(EN > E(JMIN)) then
                     DO J=JMIN,JJ
                        IF(en > e(j-1) .and. en < e(j)) then
                           EMJ=DEL(I)*CINx(iel,ion,k)/(4.*PI*(E(J)-E(J-1)))
                           EM(I,J-1)=EM(I,J-1)+EMJ
                           if(abs(emj).gt.1.e-16) then
                              write(6,9127)i,j,iel,ion,k,wl0,DEL(I),CINx(iel,ion,k),emj,em(i,j-1)
 9127                         format(' emjp ',5i4,f7.1,1pe12.3,10e12.3)
                           endif
                        endif
                     enddo
                  endif
               endif
            enddo
         enddo
      enddo
      
      DO KA=1,24
         DO N=1,10
            INDIC(N,KA)=-1
         enddo
         INDI(KA)=-1
      enddo
      EMFFS=0.
      DO 134 K1=1,3
      E2=CHI(K1)+1.5*TE/1.16E4
      DO 135 J=JIO(K1),JJ
      IF(E2.GT.E(J)) GOTO 135
      JMAX(K1)=J-1
      GOTO 136
 135  CONTINUE
 136  IF(JMAX(K1).EQ.JIO(K1)) JMAX(K1)=JIO(K1)+1
 134  CONTINUE
      DO K=1,NL
      EMCHECK(K)=0.
      ENDDO
C     EXCITED O I EMISSION
      CALL ATDATO
      NP1O=10
      DO J=JMIN,JJ
        DO N=4,9
            CALL RECR(2,N,J,CHI(14),AZG,KS(14),TE,RE)
            DRECO=AB(3)*BOL(NP1O)*RE*DEL(I)/(4.*PI)
            EM(I,J)=EM(I,J)+DRECO
c            if(em(i,j).lt.-1.d-40) write(6,9124)i,j,n,dreco,em(i,j)
            EMC(I,J)=EMC(I,J)+DRECO
        ENDDO
      ENDDO
 9124 format(' dreco ',3i3,10e12.4)
      DO 1000 J=JMIN,JJ
C     GAUNT FACTOR
C     FREE-FREE
      F1=E1(J)/TEV
      IF(F1.GT.30.) GOTO 150
C     ONLY H
      zi=1.
      gaff=gaunt(zi,te,e1(J))
      EMFF=1.6543E-23*zi**2*EXP(-F1)*GAFF/SQRT(TE)
      GOTO 175
  150 EMFF=0.0E0
 175  EMFF=EMFF*X(1)*DEL(I)
      TREC=0.
 1000 CONTINUE
      RETURN
      END

      SUBROUTINE RATE(XEL)
C     ****************************************************************
C     *****
C     RATE CALCULATES PHOTOIONIZATION AND PHOTOHEATING RATES
C     OBS! FLUX(J) = MEAN INTENSITY (ERG/CM**2 EV)
c     xel not used!
C
C
C     *****
C     ****************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
c      PARAMETER (MD=300,MDP1=MD+1)
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/IND/IK
      COMMON/PHQ/ZE(nio),GE(nio),ZK(nio)
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/SIK/SK(ncr,NE1:NE2)
      COMMON/CSI_vern2/GSV2(30,30,7,NE1:NE2)
      COMMON/TRES/ EL(nio),EK(nio),ELord(14,27),EKord(14,27)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),ipar
      real*8 jmean
      COMMON/DTAU/JMEAN(NE1:NE2)
      common/chec/gec(nio)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      common/raug/zion(30,27,7)
      COMMON/HEIION/ZHEICO,ZHEILA,ZOTSHI,ZBALM,ZHBALM
      DIMENSION JLOWL(nio),JUPL(nio),JLOWK(nio),JUPK(nio),
     &                                    ZET(nio),ZKT(nio)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      call rateaug(6,6)
      call rateaug(7,7)
      call rateaug(8,8)
      call rateaug(10,10)
      call rateaug(11,11)
      call rateaug(12,12)
      call rateaug(14,14)
      call rateaug(16,16)
      call rateaug(18,7)
      call rateaug(26,26)
c     skip old photrates. did not inclde auger properly
      DO K=1,nio
         ZK(K)=1.D-30
         ZE(K)=1.D-30
         gec(k)=1.d-30
         GE(K)=1.D-30
      enddo

      GAA=0.0E0
      DO K=3,nio
         GAA=0.0E0
         ZA=0.0E0
         ZB=0.E0
c!!!!!
         jlowl(k)=jmin
         jupl(k)=jj
         jlowk(k)=jmin
         jupk(k)=jj

         DO J=JLOWL(K),JUPL(K)
C     ONLY GROUND STATES. 5 EV IS THE LOWEST IONIZATION POTENTIAL = NA
            IF(E(J).GT.5.0) THEN
               ZA=4.*PI*JMEAN(J)*SI(K,J)*(E(J+1)-E(J))/(ELCH*E1(J))
               if(k<=3) then
                  ZE(K)=ZE(K)+ZA
               endif
               if(k<=-3) then
                  write(6,9289)k,j,e1(j),jmean(j),si(k,j),za,ze(k)
c 9289             format('H,He ze ',2i5,1pe11.3,10e11.3)
               endif
               IF(E1(J).GT.EL(K)) THEN
                  GAA=ELCH*ZA*(E1(J)-EL(K))
               ELSE
                  GAA=0.0E0
               ENDIF
               GE(K)=GE(K)+GAA
               if((k>=59.and.k<=60)) then
c                  write(6,9289)k,j,e1(j),jmean(j),si(k,j),za,ge(k)
 9289             format('GE L  ',2i5,1pe11.3,10e11.3)
               endif
            ENDIF
         ENDDO
c      write(6,*)'k,JLOWK(K),JUPK(K) ',k,JLOWK(K),JUPK(K)
         DO J=JLOWK(K),JUPK(K)
C     ONLY GROUND STATES. 5 EV IS THE LOWEST IONIZATION POTENTIAL = NA
            IF(E(J).GT.5.0) THEN
               ZB=4.*PI*JMEAN(J)*SK(K,J)*(E(J+1)-E(J))/(ELCH*E1(J))
c     ZK(K)=ZK(K)+ZB
               IF(E1(J).GT.EK(K)) THEN
                  GAA=ELCH*ZB*(E1(J)-EK(K))
               ELSE
                  GAA=0.0E0
               ENDIF
               GE(K)=GE(K)+GAA
               if((k>=14.and. k<=17).or.(k>=4.and.k<=7)) then
c     write(6,9288)k,j,e1(j),jmean(j),sk(k,j),zb,ge(k)
 9288             format('GE  ',2i5,1pe11.3,10e11.3)
               endif
               if(k==2.or.k==3) then
c     write(6,9282)k,j,e1(j),jmean(j),zb,gaa,ge(k)
 9282             format('ge(k) ',2i5,1pe11.3,10e11.3)
               endif

            ENDIF
         ENDDO
c      enddo
      enddo 
      ZOTSHI=4.*PI*JMEAN(2)*SK(1,2)*(E(3)-E(2))/(ELCH*E1(2)*ZK(1))
      RETURN
      END

      SUBROUTINE RECR(ION,N,J,CHI,AZG,KS,TE,RE)
C
C     RECOMBINATION (USING THE MILNE RELATION)
C     AVERAGE THE EMISSION OVER THE INTERVAL SO THAT THE TOTAL EMISSION
C     IS CORRECT
C
      include 'parameters.h'
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/SIK/SK(ncr,NE1:NE2)
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(20,NE1:NE2),
     &SIGFE(NL,NE1:NE2),SIGH(12,NE1:NE2)
      COMMON/NLEV/e000,NION,NH,NP1H,NHMAX,NHMIN
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      E00=E000-EN(N)
      CHI=E00
      IF(E(J).LT.E00) GOTO 250
      F2=(CHI-E1(J))*1.1609E4/TE
      DELT=E(J+1)-E(J)
      F3=DELT*1.1609E4/(2.*TE)
      F2AB=ABS(F2)
      F3AB=ABS(F3)
      IF(F3AB.GT.75.) GOTO 250
      IF(F2AB.GT.75.) GOTO 250
c$$$            IF(N.LE.6) THEN
c$$$              GAUNT=gbf(N,e1(J))
c$$$            ELSE
c$$$              GAUNT=1.
c$$$            ENDIF
c$$$c            IF(ION.EQ.5) SIGM=1.e-18*SIG(N)*GAUNT*(E00/E1(J))**GA(N)
            IF(ION.EQ.2) SIGM=1.E-18*SIGOX(N,E1(J))
C           AVERAGE OVER THE ENERGY INTERVAL
            RE=G(N)*1.31159E-4*E1(J)**3.*(EXP(F2+F3)-EXP(F2-F3))
     &                              *SIGM/(G(NH+1)*2.*F3*TE*SQRT(TE))
            GOTO 275
  250 RE=0.0E0
  275 CONTINUE
      RETURN
      END

      SUBROUTINE RECEX(Z,J,N,TE,RECO)
      IMPLICIT REAL*8(A-H,O-Z)
C      ***********************************************************
C      *****
C     THIS ROUTINE CALCULATES THE REC. EMISSION DUE TO EXCITED STATES
C     USING A HYDROGENIC APPROXIMATION AND UNIT GAUNT FACTOR.
C     Z = CHARGE
C     J = ENERGY BIN
C     N = LEVEL
C      *****
C      ***********************************************************
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      GAUNT=gbf(n,e1(j))
      RN=DBLE(N)
      RIZ=13.597*Z**2./RN**2.
      TEV=TE/1.1609E4
      EX=ABS((RIZ-E1(J))/TEV)
      RECO=0.
      IF(EX.LT.100.) THEN
            RECO=4.*3.14159*4.1478E-19*Z**4.*GAUNT*
     &                  EXP((RIZ-E1(J))/TEV)/(RN**3.*TE**1.5)
      ENDIF
      RETURN
      END

      SUBROUTINE XIONS(X,CHI)
C
C     CALCULATE IONIZATION POTENTIALS AND TOTAL FRACTIONS OF ALL IONS
C     FOR CALCULATION OF RECOMBINATION RADIATION. 
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      COMMON/TRES/ EL(nio),EK(nio),ELord(14,27),EKord(14,27)
      COMMON/ABU/XA(2,nio)
      COMMON/ABUN/AB(20)
      DIMENSION X(nio),CHI(nio)
      DO K=1,120
         CHI(K)=EL(K)
      enddo
      CHI(1)=13.6
      CHI(2)=24.6
      CHI(3)=54.4
      CHI(5)=EK(5)
      CHI(6)=EK(6)
      CHI(10)=EK(10)
      CHI(11)=EK(11)
      CHI(23)=EK(23)
      CHI(24)=EK(24)
      XO9=1.-XA(2,6)-XA(2,5)-XA(2,4)-XA(2,7)-
     &     XA(2,17)-XA(2,16)-XA(2,15)-XA(2,14)
      XC7=1.-XA(2,11)-XA(2,10)-XA(2,9)-XA(2,8)-XA(2,13)-XA(2,12)
      XN8=1.
      DO  K=18,24
         XN8=XN8-XA(2,K)
      enddo
c H II      
      X(1)=AB(1)*(1.-XA(2,1))
c He II-III      
      X(2)=AB(2)*XA(2,3)
      X(3)=AB(2)*(1.-XA(2,2)-XA(2,3))
c O VI-IX, Strange order 4 = O VII, 5 = O VII, 6 = O IX, 7 = O VI
      X(4)=AB(3)*XA(2,5)
      X(5)=AB(3)*XA(2,6)
      X(6)=AB(3)*XO9
      X(7)=AB(3)*XA(2,4)
c C II - VII. Note strange order. 13 = C II     
      X(8)=AB(4)*XA(2,9)
      X(9)=AB(4)*XA(2,10)
      X(10)=AB(4)*XA(2,11)
      X(11)=AB(4)*XC7
      X(12)=AB(4)*XA(2,13)
      X(13)=AB(4)*XA(2,8)
c O II-V     
      X(14)=AB(3)*XA(2,15)
      X(15)=AB(3)*XA(2,16)
      X(16)=AB(3)*XA(2,17)
      X(17)=AB(3)*XA(2,7)
c N II- N VIII     
      X(18)=AB(5)*XA(2,19)
      X(19)=AB(5)*XA(2,20)
      X(20)=AB(5)*XA(2,21)
      X(21)=AB(5)*XA(2,22)
      X(22)=AB(5)*XA(2,23)
      X(23)=AB(5)*XA(2,24)
      X(24)=AB(5)*XN8
c Si II     
      X(25)=AB(6)*XA(2,26)
c Mg II
      X(39)=AB(7)*XA(2,40)
c Fe II-V
      X(42)=AB(8)*XA(2,43)
      X(43)=AB(8)*XA(2,44)
      X(44)=AB(8)*XA(2,45)
      X(45)=AB(8)*(1.-XA(2,45)-XA(2,44)-XA(2,43)-XA(2,42))
c Ca      
      X(50)=AB(11)*XA(2,51)
      X(51)=AB(11)*XA(2,52)
c Na      
      X(52)=AB(10)*XA(2,53)
      X(53)=AB(10)*XA(2,52)
      X(54)=AB(10)*XA(2,55)
c S II      
      X(56)=AB(12)*XA(2,57)
      DO K=1,100
            IF(X(K).LT.0.) X(K)=0.
      ENDDO
      RETURN
      END

      double precision function gbf(n,e)
C
C     Calculation of bound-free gaunt factor from Mihalas, D., 
C     Ap.J. 149:169 (1967)
C     e in eV
C
      implicit real*8 (a-h,o-z)
      parameter(nl=6)
      dimension a0(nl),a1(nl),a2(nl),a3(nl),am1(nl),am2(nl),am3(nl)
      dimension wl0(nl)
      data wl0/915.3291d0,3704.9034d0,8504.7783d0,15560.594d0,
     &         25260.706d0,38194.187d0/
      data a0/1.2302628d0,1.1595421d0,1.1450949d0,1.1306695d0,
     &        1.1190904d0,1.1168376d0/
      data a1/-2.9094219d-3,-2.0735860d-3,-1.9366592d-3,-1.3482273d-3,
     &        -1.0401085d-3,-8.9466573d-4/
      data a2/7.3993579d-6,2.7033384d-6,2.3572356d-6,-4.6949424d-6,
     &        -6.9943488d-6,-8.8393113e-6/
      data a3/-8.7356966d-9,0.d0,0.d0,2.3548636d-8,2.8496742d-8,
     &        3.4696768d-8/
      data am1/-5.5759888d0,-1.2709045d0,-0.55936432d0,-0.31190730d0,
     &         -0.16051018d0,-0.13075417d0/
      data am2/12.803223d0,2.1325684d0,0.52471924d0,0.19683564d0,
     &         5.5545091d-2,4.1921183d-2/
      data am3/0.d0,-2.0244141d0,-0.23387146d0,-5.4418565d-2,
     &         -8.9182854d-3,-5.5303574d-3/
      wlmy=1.239854/e
      gbf=0.
      if(wlmy.lt.wl0(n)/1.d4) then
            x=1./wlmy
            gbf=a0(n)+a1(n)*x+a2(n)*x*x+a3(n)*x*x*x+am1(n)/x+
     &                                am2(n)/(x*x)+am3(n)/(x*x*x)
      endif
      return
      end

      double precision function gaunt(z,t,e)
c     this routine calculates the free-free gaunt factor according to
c     D.G. Hummer in Ap.J. 327: 477 (1988)
c
      implicit real*8(a-h,o-z)
      dimension cm(11),c(8),d(0:7,0:10)
      data d/
     &  +8.986940175d-0,-4.009515855d-0,+8.808871266D-1,+2.640245111d-2,
     &  -4.580645915d-2,-3.568055702d-3,+2.827798067d-3,+3.365860195d-4,
     &  -8.006936989d-1,+9.466021705d-1,+9.043402532d-2,-9.608451450d-2,
     &  -1.885629865d-2,+1.050313890d-2,+2.800889961d-3,-1.078209202d-3,
     &  -3.781305103d-1,+1.102726332d-1,-1.543619180d-2,+8.310561114d-3,
     &  +2.179620525d-2,+4.259726289d-3,-4.181588794d-3,-1.770208330d-3,
     &  +1.877213132d-2,-1.004885705d-1,-5.483366378d-2,-4.520154409d-3,
     &  +8.366530426d-3,+3.700273930d-3,+6.889320423d-4,+9.460313195d-5,
     &  +7.300158392d-2,+3.576785497d-3,-4.545307025d-3,-1.017965604d-2,
     &  -9.530211924d-3,-3.450186162d-3,+1.040482914d-3,+1.407073544d-3,
     &  -1.744671550d-3,+2.864013856d-2,+1.903394837d-2,+7.091074494d-3,
     &  -9.668371391d-4,-2.999107465d-3,-1.820642230d-3,-3.874082085d-4,
     &  -1.707268366d-2,-4.694254776d-3,+1.311691517d-3,+5.316703136d-3,
     &  +5.178193095d-3,+2.451228935d-3,-2.277321615d-5,-8.182359057d-4,
     &  +2.567331664d-4,-9.155339970d-3,-6.997479192d-3,-3.571518641d-3,
     &  -2.096101038d-4,+1.553822487d-3,+1.509584686d-3,+6.212627837d-4,
     &  +4.098322531d-3,+1.635218463d-3,-5.918883504d-4,-2.333091048d-3,
     &  -2.484138313d-3,-1.359996060d-3,-5.371426147d-5,+5.553549563d-4,
     &  +3.837562402d-5,+2.938325230d-3,+2.393747064d-3,+1.328839809d-3,
     &  +9.135013312d-5,-7.137252303d-4,-7.656848158d-4,-3.504683798d-4,
     &  -8.491991820d-4,-3.615327726d-4,+3.148015257d-4,+8.909207650d-4,
     &  +9.869737522d-4,+6.134671184d-4,+1.068883394d-4,-2.046080100d-4/
      data elch /1.602d-12/      
      data bolz /1.38d-16/
      a1=-1.
      b1=1.
      a2=-10.5/5.5
      b2=0.5/5.5
      m=11
      n=8
      g2=z*15.72e4/t
      u=e*elch/(bolz*t)                
      xg=log10(g2)/3.
      xu=(2.*log10(u)-2.5)/5.5
      gff=0.
      do j=1,n
            do i=1,m
                  cm(i)=d(j-1,i-1)
            enddo
            c(j)=chebev(a1,b1,cm,m,xg)
      enddo            
      gaunt=chebev(a2,b2,c,n,xu)
      return
      end

      DOUBLE PRECISION FUNCTION CHEBEV(A,B,C,M,X)
      implicit real*8(a-h,o-z)
      DIMENSION C(M)
      D=0.
      DD=0.
      Y=(2.*X-A-B)/(B-A)
      Y2=2.*Y
      DO 11 J=M,2,-1
        SV=D
        D=Y2*D-DD+C(J)
        DD=SV
11    CONTINUE
      CHEBEV=Y*D-DD+0.5*C(1)
      RETURN
      END

      SUBROUTINE NONTH(X,FH,FIH,FEH,FIHE,FEHE)
C     NON-THERMAL ELECTRON DEPOSITION FROM SHULL AND VAN STEENBERG (1985)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.LT.0.95) THEN
       FH=0.9971*(1.-(1-X**.2663)**1.3163)
       FIH=0.3908*(1.-X**0.4092)**1.7592
       FIHE=0.0554*(1.-X**0.4614)**1.666
       FEH=0.4766*(1.-X**0.2735)**1.5221
       FEHE=0.0246*(1.-X**0.4049)**1.6594
      ELSE
       FH=1.
       FIH=0.
       FIHE=0.
       FEH=0.
       FEHE=0.
      ENDIF
      RETURN
      END



      subroutine ionization_o_si_s_ar(ik,dt,dent,denel,xel,ze)
      implicit real*8(a-h,o-z)
c     include 'PARAM'
      include 'parameters.h'
      common/ionx/xion(md,14,27)
      common/ionxold/te_old(md),xion_old(md,14,27)
      common/abl/abn(15)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      common/debug/ideb
      common/raug/zion(30,27,7)
      common/rec_coll/alrec(30,30),collion(30,30)
      dimension photok(14,26)
      dimension AA(nlp1,nlp1),X(nl)
      dimension dx(nl)
      real*8 ali(5),zi(5),abi(5),x1(5),x2(5)
      ideb=1
c  dn(i)/dt = (n(i) - nold(i))/dt = 
c             a_(i+1) n(i+1) -(a(i) + c(i) + p(i)) n(i) + (p(i-1)+c(i-1)) n(i-1)

c  n(1) + n(2) +.....+ n(nionel(iel)+1) = 1.

c    n(iel) = 2

c    n(1)/dt - a_(2) n(2) + (c(1) + p(1)) n(1) =  nold(1))/dt

c    n(2)/dt - a_(3) n(3) + (c(2) + p(2)) n(2) - (c(1) + p(1)) n(1) =  nold(2))/dt

c    n(1) + n(2)  + n(3) = 1.
      zi(1)=0.
      zi(2)=0.
      zi(3)=0.
      zi(4)=0.
      do i=1,7
         zi(1)=zion(8,1,i)+zi(1)
         zi(2)=zion(14,1,i)+zi(2)
         zi(3)=zion(16,1,i)+zi(3)
         zi(4)=zion(18,1,i)+zi(4)
         write(6,921)i,zion(8,1,i),zion(14,1,i),zion(16,1,i),zion(18,1,i)
 921     format(' zi O,Si,S,Ar ',i5,1pe12.3,10e12.3)
      enddo
      write(6,922)i,zi(1),zi(2),zi(3),zi(4)
 922  format(' zi O,Si,S,Ar ',i5,1pe12.3,10e12.3)
      ali(1)=alrec(8,2)
      ali(2)=alrec(14,2)
      ali(3)=alrec(16,2)
      ali(4)=alrec(18,2)
      
      abi(1) = abn(5)      
      abi(2) = abn(10)
      abi(3) = abn(11)
      abi(4) = abn(12)
      
      x1(1) = xion(ik,5,1)
      x2(1)=1.-x1(1)
      x1(2) = xion(ik,10,1)
      x2(2)=1.-x1(2)
      x1(3) = xion(ik,11,1)
      x2(3)=1.-x1(3)
      x1(4) = xion(ik,12,1)
      x2(4)=1.-x1(4)
      xei=xel
      xe=0.
      do i=1,4
         xe=xe + abi(i)*x2(i)
         xi1=ali(i)*xei*dent/(zi(i)+ali(i)*xei*dent)
         xi2=1.-xi1
         write(6,98)i,x1(i),x2(i),zi(i),ali(i),abi(i),xe,xei,xi1,xi2
 98      format('i,x1,x2,zi,ali,abi ',i5,1pe12.3,10e12.3)
      enddo

      write(6,93)ik,dent,denel,xe
 93   format('ik,dent,denel,xe ',i5,1pe12.3,10e12.3)

      do it=1,10
         xe_old=xe
         call newtonr(xe,zi,ali,abi,dent,f_xe,df_dxe)      
         dxe=-f_xe/df_dxe
         xe=xe_old+dxe
         write(6,92)it,f_xe,df_dxe,dxe,xe_old,xe
 92      format('it,dxe,xe_old,xe ',i5,1pe12.3,10e12.3)
      enddo

      do i=1,4
         x1(i)=ali(i)*xe*dent/(zi(i)+ali(i)*xe*dent)
         x2(i)=1. - x1(i)
         write(6,96)i,x1(i),x2(i)
 96      format('i,x1, x2 ',i5,1pe12.3,10e12.3)
      enddo

      xel=xe
      
      xion(ik,5,1) = x1(1)      
      xion(ik,5,2) = x2(1)

      xion(ik,10,1) = x1(2)      
      xion(ik,10,2) = x2(2)

      xion(ik,11,1) = x1(3)      
      xion(ik,11,2) = x2(3)

      xion(ik,12,1) = x1(4)      
      xion(ik,12,2) = x2(4) 
      
      return

      end

      subroutine newtonr(xe,zi,ali,abi,den,f_xe,df_dxe)
      implicit real*8(a-h,o-z)
      real*8 ali(5),zi(5),abi(5)
      f_xe=xe
      df_dxe=1.
      do i=1,4
         f_xe= f_xe - abi(i)*zi(i)/(zi(i)+ali(i)*den*xe)
         df_dxe=df_dxe + abi(i)*ali(i)*den*zi(i)/(zi(i)+ali(i)*den*xe)**2
         write(6,9)i,abi(i),ali(i),zi(i),den,xe,f_fxe,df_dxe
 9       format('newt ',i5,1pe12.3,10e12.3)
      enddo
      return
      end

