      DOUBLE PRECISION FUNCTION EXP1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.LT.1.) THEN
            E1=-LOG(X)-.577216+X-.24991*X*X+.0552*X*X*X-.00976
     &            *X*X*X*X
      ELSE
            E1=EXP(-X)*(X*X+2.334733*X+.250621)/((X*X+3.330657*X
     &            +1.681534)*X)
      ENDIF
      EXP1=E1
      RETURN
      END

      DOUBLE PRECISION FUNCTION E2(X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.GT.0.) GOTO 500
      E2=1.
      GOTO 400
500   IF(X.GT.300.) GOTO 100
      IF(X.GT.1.) GOTO 200
      E1=-LOG(X)-.577216+X-.24991*X*X+.0552*X*X*X-.00976
     &*X*X*X*X
      GOTO 300
200   E1=EXP(-X)*(X*X+2.334733*X+.250621)/((X*X+3.330657*X+1.681534
     &)*X)
300   E2=EXP(-X)-X*E1
      GOTO 400
100   E2=0.
400   CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION SIMUL(N,A,X,EPS,II,NRC)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      DIMENSION A(NLP1,NLP1),X(N),B(NLP1)
      DIMENSION LASTN(NLP1),ASAVE(NLP1)                                     
C                                                                       
C -- Find the column index of the last non-zero element in each row --  
C -- Speeds up the solution of loose matrices considerably --           
C
      DO I=1,N
         B(I)=A(I,N+1)
c         write(6,*)'i,b,a',i,b(i),(a(i,j),j=1,n)
      ENDDO

      DO I=1,N                                                      
         DO J=N,I,-1                                                   
            JSAVE=J                                                           
            IF(A(I,J).NE.0.0D0) GOTO 102                                      
         enddo
 102     LASTN(I)=JSAVE
      enddo


      
C     
C -- ************************** --                                      
C -- Forward elimination scheme --                                      
C -- ************************** --                                      
C                 
      DO I=1,N-1    
C     
C     -- Partial pivoting routine --                                        
C     
         AI=0.0D0                                                          
         DO L=I,N                                                      
            ASAVE(L)=0.0D0                                                    
            IF(DABS(A(L,I)).GT.AI) then
               LSAVE=L                                    
            endif
            IF(DABS(A(L,I)).GT.AI) AI=DABS(A(L,I))                            
         enddo
         IF(LSAVE.ne.I) then
C     
C     -- Interchange rows A(I,?) and A(LSAVE,?) if LSAVE.NE.I --            
C     
            LAST=MAX0(LASTN(I),LASTN(LSAVE))                                  
            DO L=I,LAST                                                   
               ASAVE(L)=A(I,L)                                                   
               A(I,L)=A(LSAVE,L)                                                 
               A(LSAVE,L)=ASAVE(L)
            enddo
            BSAVE=B(I)                                                        
            B(I)=B(LSAVE)                                                     
            B(LSAVE)=BSAVE                                                    
            LASTNI=LASTN(I)                                                   
            LASTN(I)=LASTN(LSAVE)                                             
            LASTN(LSAVE)=LASTNI
         endif                                             
C     
C     -- Elimination routine --                                             
C     
            BI=B(I)                                                           
            AII=A(I,I)                                                        
            DO  J=I+1,N                                                    
               IF(A(J,I).EQ.0.0D0.OR.(I+1).GT.LASTN(I)) GOTO 301                
               AJIAII=A(J,I)/A(I,I)                                             
               DO K=I+1,LASTN(I)                                             
                  A(J,K)=A(J,K)-AJIAII*A(I,K)
               enddo
               B(J)=B(J)-AJIAII*BI                                              
               LASTN(J)=MAX0(LASTN(J),LASTN(I))                                
 301           CONTINUE
            enddo
C     
C     -- ************************** --                                      
C     -- End of forward elimination --                                      
C     -- ************************** --                                      
C     
         enddo
C     
C -- *************** --                                                 
C -- Back-substitute --                                                 
C -- *************** --                                                 
C                                                                       
      DO K=N,1,-1                                                   
         BK=B(K)                                                           
         IF((K+1).GT.N) GOTO 401                                           
         DO L=K+1,N                                                    
            BK=BK-A(K,L)*B(L)                                                 
         enddo
 401     B(K)=BK/A(K,K)                                                    
      enddo
C                                                                       
      DO I=1,N
      X(I)=B(I)
      ENDDO
      SIMUL=0.
      RETURN                                                            
      END

      DOUBLE PRECISION FUNCTION SIMUL_ion(N,A,X,EPS,II,NRC)
      IMPLICIT REAL*8(A-H,O-Z)
C     SOLVE A*X = B
C     B CONTAINS AS INPUT THE R.H. SIDE OF THE SYSTEM, AND AS OUTPUT THE
C     SOLUTION
C     EPS, II, NRC NOT USED
      DIMENSION A(31,31),X(30),B(31)
      DIMENSION LASTN(31),ASAVE(31)                                     
C                                                                       
C -- Find the column index of the last non-zero element in each row --  
C -- Speeds up the solution of loose matrices considerably --           
C                                                                       
      DO I=1,N
         B(I)=A(I,N+1)
      ENDDO

      DO I=1,N                                                      
         DO J=N,I,-1                                                   
            JSAVE=J                                                           
            IF(A(I,J).NE.0.0D0) GOTO 102                                      
         enddo
 102     LASTN(I)=JSAVE
      enddo

c      write(6,*)'simul n',n
c      do i=1,n
c         write(6,911)i,(a(i,j),j=1,n),b(i)
c 911     format(i5,1pe12.4,40e12.4)
c      enddo
      
C     
C -- ************************** --                                      
C -- Forward elimination scheme --                                      
C -- ************************** --                                      
C                 
      DO I=1,N-1    
C     
C     -- Partial pivoting routine --                                        
C     
         AI=0.0D0                                                          
         DO L=I,N                                                      
            ASAVE(L)=0.0D0                                                    
            IF(DABS(A(L,I)).GT.AI) then
               LSAVE=L                                    
            endif
            IF(DABS(A(L,I)).GT.AI) AI=DABS(A(L,I))                            
         enddo
         IF(LSAVE.ne.I) then
C     
C     -- Interchange rows A(I,?) and A(LSAVE,?) if LSAVE.NE.I --            
C     
            LAST=MAX0(LASTN(I),LASTN(LSAVE))                                  
            DO L=I,LAST                                                   
               ASAVE(L)=A(I,L)                                                   
               A(I,L)=A(LSAVE,L)                                                 
               A(LSAVE,L)=ASAVE(L)
            enddo
            BSAVE=B(I)                                                        
            B(I)=B(LSAVE)                                                     
            B(LSAVE)=BSAVE                                                    
            LASTNI=LASTN(I)                                                   
            LASTN(I)=LASTN(LSAVE)                                             
            LASTN(LSAVE)=LASTNI
         endif                                             
C     
C     -- Elimination routine --                                             
C     
            BI=B(I)                                                           
            AII=A(I,I)                                                        
            DO  J=I+1,N                                                    
               IF(A(J,I).EQ.0.0D0.OR.(I+1).GT.LASTN(I)) GOTO 301                
               AJIAII=A(J,I)/A(I,I)                                             
               DO K=I+1,LASTN(I)                                             
                  A(J,K)=A(J,K)-AJIAII*A(I,K)
               enddo
               B(J)=B(J)-AJIAII*BI                                              
               LASTN(J)=MAX0(LASTN(J),LASTN(I))                                
 301           CONTINUE
            enddo
C     
C     -- ************************** --                                      
C     -- End of forward elimination --                                      
C     -- ************************** --                                      
C     
         enddo
C     
C -- *************** --                                                 
C -- Back-substitute --                                                 
C -- *************** --                                                 
C                                                                       
      DO K=N,1,-1                                                   
         BK=B(K)                                                           
         IF((K+1).GT.N) GOTO 401                                           
         DO L=K+1,N                                                    
            BK=BK-A(K,L)*B(L)                                                 
         enddo
 401     B(K)=BK/A(K,K)                                                    
      enddo
C                                                                       
      DO I=1,N
      X(I)=B(I)
      ENDDO
      SIMUL=0.
      RETURN                                                            
      END
     

c$$$      SUBROUTINE ROMB(AA,BB,ERR,FI,FAIL,K)
c$$$      IMPLICIT REAL*8(A-H,O-Z)
c$$$      DIMENSION TR(100),W(100)
c$$$      write(6,*)' i romb '
c$$$      A=DBLE(AA)
c$$$      B=DBLE(BB)
c$$$      DX=(B-A)/2.
c$$$      FI1=0.5*(FX(A)+FX(B))
c$$$      FI2=FX(A+DX)
c$$$      TR(1)=DX*(FI1+FI2)
c$$$      N=1
c$$$      K=1
c$$$      W(1)=4.
c$$$1     FI=TR(1)
c$$$      W(K+1)=4.*W(K)
c$$$      N=2*N
c$$$      K=K+1
c$$$      TDX=DX
c$$$      DX=DX/2.
c$$$      X=A+DX
c$$$      DO I=1,N
c$$$         FI2=FI2+FX(X)
c$$$ 2       X=X+TDX
c$$$      enddo
c$$$      TR(K)=DX*(FI1+FI2)
c$$$      KK=K-1
c$$$      DO J=1,KK
c$$$         M=K-J
c$$$         TR(M)=(W(J)*TR(M+1)-TR(M))/(W(J)-1.)
c$$$      enddo
c$$$      FAIL=ABS(TR(1)-FI)/ABS(TR(1))
c$$$      IF(FAIL-ERR)6,6,4
c$$$4     IF(K-100)1,6,6
c$$$6     FI=TR(1)
c$$$      RETURN
c$$$      END

      SUBROUTINE SNPAR(TDAY,RPH,TEFF,RLUM)
      IMPLICIT REAL*8(A-H,O-Z)
C     PARAMETERS FOR SN 1980K D=8MPC EB-V=0.35 T SINCE MAX
c      write(6,*)' i snpar '
      IF(TDAY.LT.90.) BM0=0.0582*TDAY-19.49
      IF(TDAY.GT.90.) BM0=0.01005*TDAY-14.87
      IF(TDAY.LT.60.) BV0=-.05+1.3333E-2*TDAY
      IF(TDAY.GT.60.) BV0=0.75
      IF(TDAY.LT.60.) TEFF=7300./(.55+1.3333E-2*TDAY)
      IF(TDAY.GT.60.) TEFF=5407.
      TEFF=6400.
      BC=-42.54+10.*LOG10(TEFF)+2.9E4/TEFF
      BOL=BM0-BC
      RPH=6.96E10*10.**(8.472-2.*LOG10(TEFF)-0.2*BOL)
      RPH=2.E14
      rph=1.e15
      RLUM=4.*3.1415*RPH**2*5.67E-5*TEFF**4
      RETURN
      END


      SUBROUTINE STRUC_HW(initstruc,filling,pwn_vel,RM,RM1CGS,DENI,RI)
C
C     READS THE STRUCTURE OF THE SUPERNOVA FROM FILE 12
C     DESIGNED FOR THE ENSMAN-WOOSLEY MODELS, BUT MAY BE USED IN GENERAL
C     RM= MASS FROM CENTER
C     IF ISTRU=1 USE DENSITY FROM MODEL, IF NOT USE A UNIFORM
C     DENSITY SHELL
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MR,MSUN,MTOT,MTOTS,MNI,MHE,MOX,MAXM
      real*8 sc,ti,v,cr,mn,co,ni,ni56,ni57,ti44,co60,na22,rhor,rr,
     &     mc(2000)

      include "parameters.h"
      COMMON/IND/II
      COMMON/PHY/DEN(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/TEM/TE(MD)
      COMMON/ABUN/AB(20)
      COMMON/MASSES/RMI(MD)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/AGNPARA/DECON,GAMMION,IPRESS,IDENS,ICDEN
      COMMON/MPAR/MTOTS,MNI,VEXP
      COMMON/HYDROS/HSCALE,DEN1,XEL,AMEAN
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      common/revs/terev,tecs
      COMMON/CHEMCOMP/ABIN(12),ISOLAR,ISNCOMP,INPUTCOMP
      COMMON/CINOUT/INOUT,IPULS
      integer avabund,zone,zams
      common/abzone/rm_abun,avabund,zams,zone
      DIMENSION ABU(15,0:600),MR(0:600),R(0:600),DENS(0:600),
     &            ABI(15),SA(15)
      DIMENSION ABIS(15),ABSUN(12),ABSUNH(12)
      DATA MSUN/1.989E33/
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
c      DATA ABSUN/7.500E-01,2.500E-01,3.981E-03,1.288E-03,1.047E-02,
c     &      1.660E-03,6.457E-04,9.333E-04,5.128614E-04,2.512E-04,
c     &      7.943E-05,2.239E-03/
c     Ferland and Netzer
c      DATA ABSUN/0.909,.0909,3.36E-4,1.09E-4,6.18E-4,1.0E-4,0.30E-4,
c     &0.29E-4,0.145E-4,1.E-30,1.E-30,0.236E-4/
c     NRL test
      DATA ABSUN/0.909,.0909,2.73E-4,0.91E-4,5.45E-4,1.4E-4,0.27E-4,
     &0.27E-4,0.135E-4,1.E-30,1.E-30,0.91E-7/
C     SOLAR COMPOSITION FROM GREVESSE AND ANDERS 1989, COSMIC ABUNDANCES OF
C     MATTER, AIP CONFERENCE PROC 183.
C
      DATA ABSUNH/1.000,.098,3.63E-4,1.12E-4,8.51E-4,1.23E-4,3.80E-5,
     &     3.55E-5,1.62E-5,3.63E-6,2.29E-6,4.68E-5/
      real*8 avab15M(14,3),avab19M(14,3)
c Woosley & Heger 2007 averaged abundancies for Si-S-Ar, Si-O
c 15 M model
c     2      1.920       1.950     
c     3      1.950       2.030
c     4      2.030       2.250      
      data avab15M/1.01205000e-04,1.63040065e-01,2.01455000e-07,2.72144500e-07,
     &     2.28400000e-06,1.55300000e-05,1.55949000e-08,2.05690000e-05,
     &     2.67664500e-05,4.70364465e-02,8.56510595e-02,2.33210845e-02,
     &     2.82279700e-02,6.46845735e-0,
     &     8.74983333e-06,4.74133333e-05,1.05386167e-04,6.19982000e-06,
     &     1.90592173e-01,7.10216667e-05,9.49333333e-10,3.53800000e-04,
     &     2.35483333e-04,3.95733333e-01,3.13716667e-01,4.12950000e-02,
     &     2.23905000e-02,3.37024064e-02,
     &     1.80096111e-07,1.73950000e-05,2.02028889e-03,2.60894444e-05,
     &     8.66511111e-01,2.56902000e-03,1.41692500e-08,3.07298333e-02,
     &     3.14100000e-03,8.19922222e-02,1.17863333e-02,6.80733333e-04,
     &     4.37582778e-05,3.77344617e-04/
c     3      1.792       1.913
c     4      1.913       2.100
c     5      2.100       3.990
      data avab19M/1.1935e-06,5.6786e-05,4.7497e-05,1.8208e-06,
     &     4.1427e-02,3.6125e-05,8.3258e-10,2.2922e-04,3.1072e-04,
     &     4.4441e-01,3.5992e-01,4.6443e-02,2.6946e-02,7.735e-2,
     &     2.0641e-07,2.3027e-05,1.5838e-04,1.9161e-05,5.7963e-01,
     &     1.5167e-04,1.0018e-09,6.4554e-04,1.8519e-04,2.2556e-01,
     &     1.5597e-01,2.7422e-02,8.2943e-03,1.3957e-03,
     &     6.7771e-09,1.0999e-05,4.0353e-02,3.4221e-05,7.4496e-01,
     &     1.6739e-01,1.8483e-06,3.1177e-02,2.8345e-03,6.7909e-03,
     &     1.1701e-03,1.4706e-04,3.5111e-05,4.3356e-04/      

      save imax,mc,abu,dcore,r_pwn
c      write(6,*)'initstruc ',initstruc
      if(initstruc==1) then
         INITstruc=0
         
c     abundaces ffor 19 M solar model.  fractions by number!
         open(12,file='./ATDAT/WH2007_19M_numb_dens.dat',status='old')           
         read(12,*)imax
         read(12,*)tt
         do i=1,imax
            read(12,*)iq,mc(i),velr,(abu(k,i),k=1,13),sc,ti,v,cr,mn,
     &           fe,co,ni,ni56,ni57,ti44,co60,na22,rhor,rr
            abu(14,i)=fe+ni56+ni57
         enddo

c     oxygen core mass = 2.0 M_sun in HW 19 M model
c     assume velocity=2000 km/s for core and uniform mass density
         mcore=2.0
         vcore=2.e3*1.e5
         vol=4.*pi*((vcore*tday*8.64e4)**3-(pwn_vel*tday*8.64e4)**3)
     &        /3.
         dcore=msun*mcore/vol
c     correct for filling factor. THIS IS DONE IN TRANS
c     right after call to struc_hw()
c         dcore=dcore/filling
         r_pwn=pwn_vel*1.e5*tday*8.64e4
         write(6,*)'core density and radius without WITHOUT corr. for filling factor ',dcore,r_pwn
         close(12)
      endif
      if(avabund==0) then
c     select abundance zone. one specific mass
c     rm_abun=1.85
c     rm_abun=2.0
c     rm_abun=1.95
c     rm_abun=1.76
c     write(6,*)' Abundance rm ',rm_abun
         do i=1,imax
            if(rm_abun >= mc(i) .and. rm_abun < mc(i+1)) then
               ik=i
               do k=1,14
                  abi(k)=(abu(k,i)+abu(k,i+1))/2.
               enddo
c               write(6,*)'rm_abund ',i,imax,rm_abun,mc(i),mc(i+1)
            endif
         enddo
c         write(6,*)' abundances at mass ',rm_abun,(abi(k),k=1,14)
      elseif(avabund==1) then
c         write(6,*) 'averaged abund. for zone ',zams,zone
         if(zams==15) then
c     avergaed WH2007 15 M abundances         
            do k=1,14
               abi(k)=avab15M(k,zone)
            enddo         
         elseif(zams==19) then
c     avergaed WH2007 19 M abundances         
            do k=1,14
               abi(k)=avab19M(k,zone)
            enddo
         endif
      endif
      
      rit0=(3.*rm*msun/(4.*pi*dcore) + r_pwn**3)**(0.333333333333)
      rr=(3.*rm*msun/(4.*pi*dcore))**(0.3333333333)
c      write(6,988)ii,rm,pwn_vel,dcore,rr,r_pwn,rit0
 988  format(' i,dcore,r_pwn,rit0 ',i5,1pe15.7,10e15.7)
      RIT0=RIT0/1.E15
c     translate to old abund. numbering. ab() by number
      AB(1)=ABI(1)
      AB(2)=ABI(2)
      AB(3)=ABI(5)
      AB(4)=ABI(3)
      AB(5)=ABI(4)
      AB(6)=ABI(10)
      AB(7)=ABI(8)
      AB(8)=ABI(14)
      AB(9)=ABI(9)
      AB(10)=ABi(7)
      AB(11)=ABI(13)
      AB(12)=ABI(11)
      AB(13)=ABI(6)
      AB(14)=ABI(12)
C     MEAN ATOMIC WEIGHT PER ION
      AMEAN=AB(1)+4.*AB(2)+12.*AB(4)+14.*AB(5)+16.*AB(3)+
     &     28.*AB(6)+56.*AB(8)+26.*AB(9)+24.*AB(7)+23.*AB(10)+
     &     40.*AB(11)+32.*AB(12)+20.*AB(13)+36.*AB(14)
      DENI=dcore/(AMU*AMEAN)
      DEN1=DENI
 444  continue
      DET0=DENI*AMEAN*1.67E-24      
      DEN1=DENI
      ri=rit0
      RETURN
      END

      DOUBLE PRECISION FUNCTION FMEAN(j,M,E_J)
      IMPLICIT REAL*8(A-H,O-Z)
c
c fmean = mean intensity. for point source in sph. geometry J = fmean  = L/(16 pi**2 r**2)
c flux = 4 pi J = 4 pi fmean
c
      SAVE
C
C     IONIZING SPECTRUM
C
      include "parameters.h"
      COMMON/RADIE/R(MD)
      COMMON/IND/I
      COMMON/TPAR/RIN,DR,R1,TDAYS
      COMMON/QSPEC/GAMMA,ALQSO
      COMMON/FNORM/FNORM
      COMMON/HYDROS/HSCALE,DEN1,XEQ,AMEAN
      COMMON/RQW/TEFF,RQ
      COMMON/PULSAR/RLTOT,ALFA,ebreak,EMIN,EMAX,teff_ns,nsterm,ibreak
      COMMON/CINOUT/INOUT,IPULS
      COMMON/CICS/RSHOCK,ics
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/INSPEC/IBLACKB,IPULSSP,IFERNET,INLR,IFERMAT,ICRAB
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      common/ns_spec/tot_lum,FL_ns(NE1:NE2),nstemp
      DIMENSION EI(1000),FSH(1000)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DIMENSION RYCRB(10),FNUCRB(10)
*     Crab Nebula continuum from Davidson and Fesen, 1985 Ann Rev
*     second vector is F_nu as seen from Earth
      DATA RYCRB/-5, -3.284, -2.823, -0.959, -0.137, 0.863, 1.863, 
     & 3.863, 6.176, 6.869 /
     &  FNUCRB/-20.423, -20.860, -20.677, -22.308, -22.721,
     & -23.650, -25.192, -27.395, -30.682, -31.780/
C
C     FMEAN = MEAN INTENSITY IN ERGS / CM**2 S EV ster
C
      FMEAN=0.
      IF(ICS.EQ.1.AND.INITCS.EQ.0.AND.IPULSSP.NE.1) THEN
        INITCS=1  
        READ(19,*)JMAX
        ftot=0.
        DO J=1,JMAX
          READ(19,*)EI(J),FSH(J)
          if(j.gt.1) ftot=ftot+fsh(j)*(ei(j)-ei(j-1))
        ENDDO
      ENDIF
      IF(IBLACKB.EQ.1) THEN
        TEFFEV=TEFF/1.1609E4
        IF(E_J/TEFFEV.LT.100.) FMEAN=5.0403E10*E_J**3/(EXP(E_J/TEFFEV)-1.)
      ENDIF
      IF(IFREEFREE.EQ.1) THEN
C     FREE-FREE SPECTRUM
        RL0=RLTOT*EXP(-E_J/EMAX)/EMAX
        FMEAN=RL0/(4.*PI**2*(R1*1.e15)**2)
        FMEAN=RL0/(4.*PI**2*R(i)**2)
      ENDIF
      IF(IPULSSP.EQ.1.AND.ICRAB.EQ.0) THEN
c     PWN SPECTRUM

c constant power law spectrum
        CON=RLTOT*(1.-ALFA)/(EMIN*((EMAX/EMIN)**(1.-ALFA)-1.))
        RL0=CON*(EMIN/E_J)**ALFA
        con1=con
        rl01=rl0
        
c$$$c broken power law at X-rays
c$$$c normalize to lum between 13.6 and infinity 
c$$$c break for L=1e36 and alpha=1.1
c$$$        ebreak=4.3e3 
c$$$c 
c$$$        ebreak=1.e4
c$$$c break energy for Greco lum and power law with 13.6 - inf lum 1.e36
c$$$        ebreak=4.31e3
c$$$c     d:o for 5e35
c$$$        ebreak=1.34e4
c$$$      ebreak=1.34e14
c$$$c     d:o for 2e36 erg/s
c$$$c        ebreak=1.13e3
c power law for e > ebreak   
        beta=1.8
        e0=13.6
        euv=1.0
        aluv=0.6
        fac=1.-(e0/ebreak)**(alfa-1.)+
     &       (e0/ebreak)**(alfa-1.)*(alfa-1.)/(beta-1.)
        con=rltot*(alfa-1.)/(e0*fac)
        if(e_j<ebreak.and.e_j>euv) then
           fac=(e0/e_j)**alfa
        elseif(ibreak==0 .and. e_j>euv) then
           fac=(e0/e_j)**alfa           
        elseif(e_j<euv) then
           fac=(e0/euv)**alfa*(euv/e_j)**aluv
        elseif(ibreak==1 .and. e_j>ebreak) then
           fac=(e0/ebreak)**alfa * (ebreak/e_j)**beta
        endif
        
c fix the spectrum to the luminosity density at 13.6 iundependent of break
c only depoendent on the total luminosity WITHOUT break
        emax=2.e4
        dlde13p6=rltot*(alfa-1.)/(e0*(1.-(e0/emax)**(alfa-1.)))
        if(e_j<ebreak.and.e_j>euv) then
           spec=(e0/e_j)**alfa
        elseif(ibreak==0 .and. e_j>euv) then
           spec=(e0/e_j)**alfa           
        elseif(e_j<euv) then
           spec=(e0/euv)**alfa*(euv/e_j)**aluv
        elseif(ibreak==1 .and. e_j>ebreak) then
           spec=(e0/ebreak)**alfa * (ebreak/e_j)**beta
        endif        
        rl0=dlde13p6*spec
        
c        FMEAN1=RL01/(16.*PI**2*(R1*1.e15)**2)
        FMEAN=RL0/(16.*PI**2*(r(i)*1.e15)**2)
        if(nsterm==1) then
c     NS Blackbody
           teff_ev=teff_ns/1.1609e4
           dlde=rltot*15.*e_j**3/(pi**4*teff_ev**4*(exp(e_j/teff_ev)-1.))
           if(nstemp==2 .or. nstemp==3) then              
              dlde=fl_ns(j)
           endif           
           fmean=dlde/(16.*PI**2*(R(i)*1.e15)**2)
        endif
        
      ENDIF
      IF(ICS.EQ.1.AND.IPULSSP.NE.1) THEN
C     INTERPOLATE
        DO J=1,JMAX
          IF(E_J.LT.EI(J)) THEN
            FMEAN=FSH(J-1)+(FSH(J)-FSH(J-1))*(E_J-EI(J-1))/
     &                                        (EI(J)-EI(J-1))
            GOTO 33
          ENDIF
        ENDDO
33      CONTINUE
      ENDIF
      IF(INOUT.EQ.1) THEN
        W=0.5*(1.-SQRT(1.-1./RQ**2))
      ELSE
        W=1.
      ENDIF
C     FACTOR OF 2 ALREADY INCLUDED IN CSREM
      IF(IAGN.EQ.1.OR.ICRAB.EQ.1) THEN
        IF(ICRAB.EQ.0) THEN
C         AGN SPECTRUM
          IF(IFERNET.EQ.1) THEN
C           FERLAND & NETZER CALC.
            U=8.7E8
            D0=2.2E9
          ELSEIF(INLR.EQ.1) THEN
c           nlr test
            u=3.e10*0.01
            d0=1.e3
            IF(E_J.GT.1.E5) THEN
              FMEAN=0.
            ELSEIF(E_J.GT.1.) THEN
               FMEAN=1.602E-12*D0*U*(13.6/e_j)/(4.*pi)
            ENDIF
          ELSEIF(IFERMAT.EQ.1) THEN
C           FERLAND AND MATHEWS BLR SPECTRUM
            RYD=E_J/13.5987
*           generic AGN continuum taken from Mathews and Ferland ApJ Dec15 '87
*           except self absorption at 100 microns, nu**-2.5 below that
            IF(RYD.GT.1.E-5.AND.RYD.LT.9.12E-4) THEN
              EC=1.E-5
              FC=0.1115
              SLOPE=-2.5
            ELSEIF(RYD.GT.9.912E-4.AND.RYD.LT.0.206) THEN
              EC=9.91E-4
              FC=5.0115
              SLOPE=1.00
            ELSEIF(RYD.GT.0.206.AND.RYD.LT.1.743) THEN
              EC=0.206
              FC=2.6576
              SLOPE=0.50      
            ELSEIF(RYD.GT.1.743.AND.RYD.LT.4.130) THEN
              EC=1.743
              FC=2.194
              SLOPE=1.00
            ELSEIF(RYD.GT.4.130.AND.RYD.LT.26.84) THEN
              EC=4.13
              FC=1.8196
              SLOPE=3.00      
            ELSEIF(RYD.GT.26.84.AND.RYD.LT.7350.) THEN
              EC=26.84
              FC=-0.6192
              SLOPE=0.70      
            ELSEIF(RYD.GT.7350.) THEN
              EC=7350.
              FC=-7.34
              SLOPE=0.0000000001
            ENDIF      
          ENDIF
        ELSEIF(ICRAB.EQ.1) THEN
C     CRAB SPECTRUM
          RYD=E_J/13.5987
          RYDL=LOG10(RYD)
          DO K=1,10
            IF(RYDL.GT.RYCRB(K).AND.RYDL.LT.RYCRB(K+1)) THEN
              FC=FNUCRB(K)
              SLOPE=(-FNUCRB(K+1) + FNUCRB(K))/(RYCRB(K+1) - RYCRB(K))
              EC=10.**RYCRB(K)
              GOTO 32
            ENDIF
          ENDDO
32        CONTINUE
        ENDIF
        FMEAN=10.**FC*(EC/RYD)**SLOPE
        FMEAN=FNORM*FMEAN
      ENDIF
C!!   ???
      IF(ICS.EQ.1.AND.IPULSSP.NE.1) W=1.
      IF(IPULSSP.eq.1) W=1.
      
      FMEAN=W*FMEAN
C
C     D:O IN ERGS / CM**2 S HZ
C
      IF(M.EQ.2) FMEAN=FMEAN*4.138E-15
C
C     OCCUPATION NUMBER
C
      IF(M.EQ.3) FMEAN=2.00E-11*FMEAN/E_J**3
      

      RETURN
      END

      subroutine ns_atm
      implicit real*8(a-h,o-z)
      COMMON/FRE/NINT,JMIN,JJ
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      real*8 en(1000),flspec(1000)
      integer nstemp
      character*120 fn
      common/ns_spec/tot_lum,FL_ns(NE1:NE2),nstemp
      k=1
      if(nstemp==2) then
         fn='./ATDAT/C_atm_Ho_Heinke_2e6K.dat,'
      elseif(nstemp==3) then
         fn='./ATDAT/C_atm_Ho_Heinke_3e6K.dat'
      endif
      open(11,file=fn,status='old')
      do i=1,1000
         read(11,*,end=22)en(i),flspec(i)
         en(i)=en(i)*1.e3
      enddo
 22   imax=i-1
      do j=jmin,jj
         ei=e1(j)
         if(ei < en(1)) then
            fl_ns(j)=flspec(1)*(ei/en(1))**2
         else
            do i=1,imax-1
               if(ei >= en(i).and.ei <=en(i+1)) then
                  fl_ns(j)=flspec(i) + (ei-en(i))*
     &                 (flspec(i+1)-flspec(i))/(en(i+1)-en(i))
               endif
            enddo
         endif
      enddo
      fltot=0.
      do j=jmin,jj
         fltot=fltot+fl_ns(j)*(e(j+1)-e(j))
      enddo
 
      con=tot_lum/fltot
      do j=jmin,jj
         fl_ns(j)=con*fl_ns(j)
      enddo
      return
      end


      SUBROUTINE AGNION(INIT,GAMMION,DENS,FSUM)
C     NORMALIZE THE AGN SPECTRUM TO GIVE A CERTAIN IONIZATION PARAMETER 
C     GAMMION.
C     IF INIT = 0 CALCULATE NORMALIZATION FACTOR
C     IF INIT = 1 CALCULATE IONIZATION PARAMETER
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/FNORM/FNORM
      COMMON/FRE/NINT,JMIN,JJ
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      FSUM=0.
      DO J=JMIN,JJ
        IF(E1(J).GT.13.5987) THEN
          FSUM=FSUM+FMEAN(j,1,E1(J))*(E(J+1)-E(J))/(ELCH*E1(J))
        ENDIF
      ENDDO
      FSUM=4.*PI*FSUM/(3.E10*DENS)
      IF(INIT.EQ.0) THEN
        FNORM=GAMMION/FSUM
      ENDIF
      RETURN
      END


      DOUBLE PRECISION FUNCTION PLANCK(M,E1,TE)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     PLANCK FUNCTION
C
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
C
C     PLANCK = MEAN INTENSITY IN ERGS / CM**2 S EV
C
      PLANCK=0.
      IF(PLANCK.LE.0.) PLANCK=1.D-30
      TEEV=TE/1.1609E4
      IF(E1/TEEV.LT.100.) PLANCK=5.0403E10*E1**3/(EXP(E1/TEEV)-1.)
C
C     D:O IN ERGS / CM**2 S HZ
C
      IF(M.EQ.2) PLANCK=PLANCK*4.138E-15
C
C     OCCUPATION NUMBER
C
      IF(M.EQ.3) PLANCK=2.00E-11*PLANCK/E1**3
      RETURN
      END

      SUBROUTINE SPEC(TEMP)
C     ***********************************************************
C     *****
C     FOR THE FIRST ITERATION (NITSPH=1) SOLVE THE EQUATION OF TRANSFER
C         IN THE OUTWARD ONLY APPROXIMATION. FOR THE NEXT ITERATIONS
C         (NITSPH > 1) CALCULATE THE CONTINOUS OPACITIES FOR THE SHELL.
C     *****
C     ***********************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      SAVE
      COMMON/INSPEC/IBLACKB,IPULSSP,IFERNET,INLR,IFERMAT,ICRAB
      COMMON/UVBL/IUVBLANK
      COMMON/EQUIV/RADF,FR(2,100),W(100),CIN(100),WLI(NFEL+100),FB(100)
      COMMON/IND/I
      COMMON/SPH/ISPH
      COMMON/CINOUT/INOUT,IPULS
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/PHY/DEN(MD)
      COMMON/DXA/DX(MD)
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/ABU/XA(2,nio)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/HYPOP/XNH(NL)
      COMMON/LITER/N
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(20,NE1:NE2),
     &SIGFE(NL,NE1:NE2),SIGH(12,NE1:NE2)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/TEM/TE(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/SIK/SK(ncr,NE1:NE2)
      COMMON/ABUN/AB(20)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),ipar
      common/sorc/so(NE1:NE2),KM(NE1:NE2)
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/NLEV/e000,NION,NH,NP1H,NHMAX,NHMIN
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),A(NL,NL),WL(NL,NL),
     &     DCDT(NL,NL),DCIDT(NL)
      COMMON/CICS/RSHOCK,ics
      DIMENSION TA(NE1:NE2),XP(nio),tas(50),S(MD,NE1:NE2)
      DIMENSION SIGM(5,15,NE1:NE2),delta_tau(1000)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      CALL IONABUND(XP)
      DO K=1,100
         IF(XP(K).Lt.-1.e-5) WRITE(6,7667)K,XP(K)
         IF(XP(K).LE.0.) XP(K)=1.E-19
 7667    FORMAT(1X,'NEG. ABUND',I5,E12.5)
 7567    FORMAT(1X,'ABUND',I5,E12.5)
      ENDDO
      K=91
      IF(XN1.LE.0.) XN1=1.E-37
      K=92
      IF(XN2.LE.0.) XN2=1.E-37
      K=93
      IF(XN3.LE.0.) XN3=1.E-37
C     
C     ELECTRON SCATTERING DEPTH
C     
C     TEL=1.2*0.665E-24*DEL(I)*DEN(I)*DX(I)+TEL
      RADF=0.
      CALL ATDATO
      DO J=JMIN,JJ
         DO K=4,13
            SIGM(2,K,J)=1.E-18*SIGOX(K,E1(J))
         ENDDO
      ENDDO
      DO J=JMIN,JJ
C     ***********************************************************
C     *****
C     CALCULATE OPTICAL DEPTHS AT EACH ENERGY
C     *****
C     ***********************************************************
c     BF absorption 
         TA(J)=0.0E0
         DTAM=0.

         DO K=1,94
            DTA=XP(K)*(SI(K,J)+SK(K,J))
            IF(DTA.GT.DTAM) THEN
               KM(J)=K
               DTAM=DTA
            ENDIF
            TA(J)=TA(J)+DTA
            delta_tau(k)=dta
            if((k==6.or.k==8.or.k==9).and.j>=-107.and.j<1) then
c               write(83,9252)k,j,xp(k),si(k,j),sk(k,j),dta
 9252          format('k,j,xp,si,sk',2i5,1pe12.3,10e12.3)

            endif
            kk=k
         ENDDO
C     EXCITED HYDROGEN CONTINUA

         DO K=2,5
            DTAB=AB(1)*XNH(K)*SIGM(5,K,J)
            TA(J)=TA(J)+DTAB
            IF(DTAB.GT.DTAM) THEN
               KM(J)=99+K
               DTAM=DTAB
            ENDIF
            kk=kk+1
            delta_tau(k)=dtab
         ENDDO
C     EXCITED O I CONTINUA
         DO K=4,13
            DTAB=AB(3)*XN(2,K)*SIGM(2,K,J)
            TA(J)=TA(J)+DTAB
            IF(DTAB.GT.DTAM) THEN
               KM(J)=101+K
               DTAM=DTAB
            ENDIF
            kk=kk+1
            delta_tau(kk)=dtab
            if(i>150.and.e1(j)>0.5 .and. e1(j) < 10 ) then
c               write(83,*)'ab(3),xn(2,k),sigm(2,k,j)',k,ab(3),xn(2,k),sigm(2,k,j),dtab
            endif
         ENDDO
C     FREE-FREE DEPTH
         ZI=1.
         GFF=gaunt(zi,temp,e1(J))
         EX1=E1(J)*1.1609E4/TEMP
         IF(EX1.LT.300.) DTAB=2.6119E-35*DEN(I)*DEL(I)*(1.-XP(1))*
     &        GFF*(1.-EXP(-EX1))/(E1(J)**3*SQRT(TEMP))
         TA(J)=TA(J)+DTAB
         IF(DTAB.GT.DTAM) THEN
            KM(J)=111
            DTAM=DTAB
         ENDIF
         kk=kk+1
         delta_tau(kk)=dtab

C     
C     TAX=TAU(I)=OPTICAL DEPTH OF SLAB BETWEEN I AND I+1
C     
         TAX=DX(I)*TA(J)*DEN(I)
         TAU(I,J)=TAX
         TAUTOT(I,J)=TAUTOT(I-1,J)+TAU(I,J)

         if(e1(j)>1. .and. e1(j) < 10 .and.i > 150) then
c            write(83,*)'i,j,e1(j) ',i,j,e1(j),ta(j),tau(i,j),tautot(i,j)
            do k=1,kk
               if(delta_tau(k) > 1.e-30) then
c                  write(83,9112)k,delta_tau(k)
 9112             format(i5,1pe12.3)
               endif
            enddo
         endif
         
         IF(TA(J).GT.0.) S(I,J)=EM(I,J)*DEN(I)/TA(J)
         IF(TA(J).LE.0.) S(I,J)=0.
         SO(J)=S(I,J)
         IF(IPULSSP.EQ. 0) THEN
            IF(INOUT.EQ.1.AND.ICS.EQ.0) CALL RADTRANIO(J,TAX,S)      
         ELSE
c            write(6,*)' to radtranpuls',j,tax
            CALL RADTRANPULS(0,J,TAX,S)      
         ENDIF
         if(i >= 30000) then 
            write(6,92)i,j,e1(j),dx(i),den(i),em(i,j),ta(j),s(i,j),tax,
     &           tautot(i,j),tautot(i-1,j),fl(1,j),fl(2,j)
         endif
 92      format('i,j,e1,dx,den,em,ta,s,tax ',2i5,1pe12.3,10e12.3)
      ENDDO
      RETURN
      END

      SUBROUTINE DIFF(IMAX,FL0)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
C      *************************************************************
C     THIS ROUTINE CALCULATES THE DIFFUSE EMISSION FOR GIVEN EMISSION
C     RATE (PER STERADIAN) AND OPTICAL DEPTHS OF THE DIFFERENT SLABS.
C     THE METHOD IS THE SAME AS THAT USED BY WILLIAMS (1967).
C     FL0(MA,J), MA=1,2 IS THE FLUX FROM THE SURFACE AND THE INNER
C     BOUNDARY. THE OPTICAL DEPTH POINTS ARE CALCULATED FROM THE SUR-
C     FACE WITH TAU(1,J)=0.
C      *****
C      *************************************************************
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),ipar
      COMMON/PHY/DEN(MD)
      COMMON/DXA/DX(MD)
      COMMON/ELEC/DEL(MD)
      DIMENSION FL0(2,NE1:NE2)
      write(6,*)' i diff(' 
      DO 500 J=JMIN,JJ
      DO 500 MA=1,2
 500  FL0(MA,J)=0.
9478  FORMAT(' DIFF ',5E12.5)
      DO 400 J=JMIN,JJ
      DO 410 I=2,IMAX
      FD(I,J)=0.0E0
       IMM1=IMAX-1
      DO 430 K=2,IMM1
C
C     CALCULATE THE OPTICAL DEPTH BETWEEN SHELL I AND SHELL K
C
      I1=I-1
      K3=K-1
      K4=K+1
      IF(K.GT.I) GOTO 130
      T2=TAUTOT(I,J)-TAUTOT(K,J)
      T3=TAUTOT(I,J)-TAUTOT(K-1,J)
      GOTO 150
 130  T2=TAUTOT(K-1,J)-TAUTOT(I,J)
      T3=TAUTOT(K,J)-TAUTOT(I,J)
 150  CONTINUE
      IF(DX(K).LE.0.) WRITE(6,*)' DIFX',K,DX(K)
      IF(DEN(K).LE.0.) WRITE(6,*)' DIFX',K,DEN(K)
      TA=(TAUTOT(K,J)-TAUTOT(K-1,J))/(DX(K)*DEN(K))
      TA=ABS(TA)
      ALB=1.2*0.665E-24*DEL(K)/(1.2*0.665E-24*DEL(K)+TA)
      ALB=0.
C
C     ALBEDO IN THE CASE OF A FINITE ELECTRON SCATTERING DEPTH
C
C
C     SOURCE FUNCTION EVALUATED IN THE MIDDLE BETWEEN K AND K+1
C
      IF(TA.LE.0.) WRITE(6,9379)K,J,EM(K,J),TAUTOT(K,J),TAUTOT(K-1,J)
      IF(TA.LE.0.) EMTA=0.
9379  FORMAT(' DIFF',2I5,4E11.4)
      IF(TA.GT.0.) EMTA=EM(K,J)/TA
C!!   SOURC=ALB*FTOT(K,J)+(1.-ALB)*DEN(K)*EMTA
      FDIF=(E2(T2)-E2(T3))*SOURC/2.
      IF(K.EQ.28)WRITE(6,9283)J,EM(K,J),TA,SOURC,FDIF
      IF(FDIF.GT.1.E15)WRITE(6,9284)I,K,J,EM(K,J),TA,SOURC,FDIF
9284  FORMAT(' ERRO ',3I4,6E11.4)
9283  FORMAT(' DI ',I4,6E12.5)
C
C     CALCULATE THE FLUX AT THE SURFACE AND AT THE INNER BOUNDARY.
C
      IF(I.NE.2) GOTO 480
      MA=1
      GOTO 495
 480  IF(I.NE.IMM1) GOTO 490
      MA=2
      E31=0.
      E32=0.
 495  IF(T2.LT.100.) E31=(EXP(-T2)-T2*E2(T2))/2.
      IF(T3.LT.100.) E32=(EXP(-T3)-T3*E2(T3))/2.
      FLUX=(E31-E32)*SOURC/2.
C*****MULTPLY BY AN ARBITARY FACTOR 2
      FLUX=2.*FLUX
      FL0(MA,J)=FL0(MA,J)+FLUX
 490  CONTINUE
      IF(FDIF.LT.1.E15) GOTO 430
      WRITE(6,998)T2,T3,E2(T2),E2(T3),EM(K,J),TAUTOT(K,J)
  998 FORMAT(6(2X,E12.3))
      GOTO 997
  430 FD(I,J)=FD(I,J)+FDIF
  410 CONTINUE
  400 CONTINUE
      WRITE(6,9387)(FD(28,J),J=JMIN,JJ)
9387  FORMAT(' FD ',5E12.5)
      RETURN
  997 WRITE(6,996)
  996 FORMAT(1X,'STOPP')
      STOP
      END


      SUBROUTINE BIS(TE,XEL,IFAIL,IFPOP,ICONV)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/THER/A1,B1,TIN,TOL,E20
      COMMON/HRAT/RHYD,ZEL,XELQ,HYR,HEAT,COOL
      COMMON/IND/Ik
C     IF ICONV = 1 AT EXIT THE FUNCTION AT THE ZERO IS LARGER THAN TOLERATED
C     IF IFAIL = 1 AT EXIT THE FUNCTION HAS THE SAME SIGN AT BOTH LIMITS
C     E20 = TOL. IN FUNCTION VALUE
      a2=0.
      a3=0.
      a4=0.
      b2=0.
      b3=0.
      b4=0.
      eps=1.
      eps4=1.e-2
      ICONV=1
      IFAIL=0
      XE=XEL
      MM=1
      A=A1
      B=B1
      FB=RAD(B1,XE,IFPOP)
      IF(IFPOP.EQ.1) GOTO 4
      mw=0
C
C     LEFT BOUNDARY 
C
 9    FA=RAD(A1,XE,IFPOP)
      if(ifpop.eq.1) a1=(a1+b1)/2.       
      if(ifpop.eq.1) mw=mw+1
      if(mw.ge.10) goto 4
      IF(IFPOP.EQ.1) GOTO 9
c 6    IF(FA*FB.GT.0.) WRITE(6,*)' SAME SIGN',MM,FA,FB
c     SAME SIGN: IFAIL=1
6     IF(FA*FB.GT.0.) IFAIL=1
      IF(MM.GT.100) GOTO 4
C     IF OK CONTINUE
      IF(IFAIL.NE.1) GOTO 1
      a4=a3
      b4=b3
      a3=a2
      b3=b2
      a2=a1
      b2=b1
      if(abs(a2-a3)/a2.lt.eps4.and.abs((fa2-fa)/fa).lt.eps4) then
         a1=0.95*a1
         A=A1
         fa2=fa
         FA=RAD(A1,XE,IFPOP)
         IF(IFPOP.EQ.1) GOTO 4
         IFAIL=0
         MM=MM+1
         GOTO 6
      endif
      if(abs(b2-b3)/b2.lt.eps4.and.abs((fb2-fb)/fb).lt.eps4) then
C     INCREASE UPPER LIMIT
         b1=1.075*b1
         b=b1
         fb2=fb
         Fb=RAD(b1,XE,IFPOP)
         IF(IFPOP.EQ.1) GOTO 4
         IFAIL=0
         MM=MM+1
         GOTO 6
      endif
      if(abs(a2-a4).lt.eps.and.abs(b2-b4).lt.eps) then
         delt=log10(a1/1.e2)/21.
         a10=a1
         do i=1,21
            b1=a1
            a1=a10*10.**(-(i-1)*delt)
            fb=fa
            fa=rad(a1,XE,IFPOP)
            if(fa*fb.lt.0.) goto 55
         enddo
 55      a=a1
         b=b1
         goto 1
      else
      endif
      IF(ABS(FB).LT.ABS(FA)) GOTO 7
      B1=A1
      B=B1
      FB=FA
      IF(A1.LE.3500.) A1=0.95*A1
      IF(A1.GT.3500.) A1=0.95*A1
      A=A1
      fa2=fa
      FA=RAD(A1,XE,IFPOP)
      IF(IFPOP.EQ.1) GOTO 4
      IFAIL=0
      MM=MM+1
      GOTO 6
7     A1=B1
      A=A1
      FA=FB
      IF(A1.LE.3500.) B1=1.05*B1
      IF(A1.GT.3500.) B1=1.051*B1
      B=B1
      fb2=fb
      FB=RAD(B1,XE,IFPOP)
      IF(IFPOP.EQ.1) GOTO 4
      IFAIL=0
      MM=MM+1
      GOTO 6
C     NEW TEMP. BETWEEN A AND B
 1    tme_old=tme
      fm_old=fm
      TME=(A+B)/2.
      mw=0
      a4=a3
      b4=b3
      a3=a2
      b3=b2
      a2=a
      b2=b
 10   FM=RAD(TME,XE,IFPOP)
      if(ifpop.eq.1) tme=(a+tme)/2.       
      if(ifpop.eq.1) mw=mw+1
      if(mw.ge.10) goto 4
      IF(IFPOP.EQ.1) GOTO 10
      IF(FA*FM.GE.0.) THEN
C     IF FM HAS THE SAME SIGN AS FA THEN A=TME
c!            IF(ABS(FM).GT.ABS(FA)) THEN
c!                  A=A-ABS(TME-A)
c!                  FA=RAD(A,XE,IFPOP)
c!            ELSE
                  A=TME
                  FA=FM
c!            ENDIF
            GOTO 3
C     IF FM HAS A DIFFERENT SIGN FROM FA THEN B=TME
      ELSE
c!            IF(ABS(FM).GT.ABS(FB)) THEN
c!                  B=B+ABS(TME-B)
c!                  FB=RAD(B,XE,IFPOP)
c!            ELSE
                  B=TME
                  FB=FM
c!            ENDIF
      ENDIF
      if(abs(a2-a3)/a2.lt.eps4.and.abs((fa2-fa)/fa).lt.eps4) then
            a1=0.95*a1
            A=A1
            fa2=fa
            FA=RAD(A1,XE,IFPOP)
            IF(IFPOP.EQ.1) GOTO 4
            IFAIL=0
            MM=MM+1
            GOTO 6
      endif
      if(abs(b2-b3)/b2.lt.eps4.and.abs((fb2-fb)/fb).lt.eps4) then
C     INCREASE UPPER LIMIT
            b1=1.075*b1
            b=b1
            fb2=fb
            Fb=RAD(b1,XE,IFPOP)
            IF(IFPOP.EQ.1) GOTO 4
            IFAIL=0
            MM=MM+1
            GOTO 6
      endif
      if(abs(a2-a4).lt.1.e-10.and.abs(b2-b4).lt.eps) then
         delt=log10(a1/1.e2)/21.
            a10=a1
            do i=1,21
                  a1=a10*10.**(-(i-1)*delt)
                  fa=rad(a1,XE,IFPOP)
                  if(fa*fb.lt.0.) goto 57
            enddo
57          a=a1
            goto 1
      endif
3     IF(ABS(FM).LE.E20*ABS(HEAT)) ICONV=0
      IF(ABS(FM).LE.E20*ABS(HEAT)) IFAIL=0
      IF(ABS(FM).LE.E20*ABS(HEAT)) GOTO 4
      if(((abs(tme_old-tme)/tme_old.lt.1.e-3)).and.ik > 10) then
         write(6,*)' temperature change very small. stop iter. ',tme,tme_old,fm,fm_old
         iconv=0
         goto 4
      endif
      MM=MM+1
      ITERM=2
      IF(MM.GE.20) then 
         if(abs(a-tme).gt.abs(b-tme)) then
            a1=0.99*a1
            fa=rad(a1,XE,IFPOP)
            a=a1
         else
            b1=1.01*b1
            fb=rad(b1,xe,ifpop)
            b=b1
         endif
         mm=1
         goto 1
      endif
      ITERM=1
C     GO ON WITH A NEW TEMPERATURE
      GOTO 1
4     TE=TME
      IF(ITER.EQ.2) WRITE(6,*)'ITERM'
      RETURN
      END
      

      SUBROUTINE RADTRANIO(J,TAX,S)
C     ***********************************************************
C     *****
C     FOR THE FIRST ITERATION (NITSPH=1) SOLVE THE EQUATION OF TRANSFER
C         IN THE SCHARMER APPROXIMATION. FOR THE NEXT ITERATIONS
C         (NITSPH > 1) CALCULATE THE CONTINOUS OPACITIES FOR THE SHELL.
C     WORKS ONLY INWARD TO OUT
C     *****
C     ***********************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      SAVE
      COMMON/IND/I
      COMMON/CINOUT/INOUT,IPULS
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/ITSPH/NITSPH
      COMMON/LITER/N
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/RADIE/R(MD)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),ipar
      real*8 jmean
      COMMON/DTAU/JMEAN(NE1:NE2)
      DIMENSION S(MD,NE1:NE2)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      write(6,*)' i radtranio '
C     ***********************************************************
C     *****
C     USE THE SCHARMER APRROX. IN THE FIRST ITERATION
C     *****
C     ***********************************************************
      IF(NITSPH.LE.1) THEN
C     BARBIER RELATION
C
C     FOR SECOND AND HIGHER ITERATIONS THE DIFFUSE FLUX IS 
C     ALREADY CALCULATED
C     MEAN INTENSITY = PRIMARY FLUX + DIFFUSE FLUX
C
            IF(IPULS.EQ.1) THEN
                  FPRIM=EXP(-TAUTOT(I,J))*FMEAN(j,1,E1(J))
            ELSE
                  FPRIM=0.
            ENDIF
            FL(2,J)=FPRIM+FD(I,J)
      ENDIF
      JMEAN(J)=FL(2,J)
      RETURN
      END


c$$$      SUBROUTINE RADTRANOI(IH,J,TAX,S)
c$$$C     ***********************************************************
c$$$C     *****
c$$$C     FOR THE FIRST ITERATION (NITSPH=1) SOLVE THE EQUATION OF TRANSFER
c$$$C         IN THE OUTWARD ONLY APPROXIMATION. FOR THE NEXT ITERATIONS
c$$$C         (NITSPH > 1) CALCULATE THE CONTINOUS OPACITIES FOR THE SHELL.
c$$$C     WORKS BOTH IN OUT AND OUT IN
c$$$C     HERE OUT TO IN
c$$$C     USED FOR ALL CS-INT. CALC.
c$$$C     IF IH = 1 CALCULATE FLUX & MEAN INTENSITY
c$$$C     IF IH = 1 CALCULATE MEAN INTENSITY
c$$$C     *****
c$$$C     ***********************************************************
c$$$      IMPLICIT REAL*8(A-H,O-Z)
c$$$      include "parameters.h"
c$$$      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
c$$$      SAVE
c$$$      COMMON/IND/I
c$$$      COMMON/CINOUT/INOUT,IPULS
c$$$      COMMON/FRE/NINT,JMIN,JJ
c$$$      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
c$$$      COMMON/ITSPH/NITSPH
c$$$      COMMON/LITER/N
c$$$      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
c$$$     &           EMC(MD,NE1:NE2)
c$$$      COMMON/RADIE/R(MD)
c$$$      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),ipar
c$$$      real*8 jmean
c$$$      COMMON/DTAU/JMEAN(NE1:NE2)
c$$$      COMMON/DIFFP/FH0(NE1:NE2),FHD(NE1:NE2)
c$$$      DIMENSION S(MD,NE1:NE2),E2Q(2),E3Q(2)
c$$$      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
c$$$      write(6,*)' i radtranoi '
c$$$C
c$$$C     PRIMARY FLUX
c$$$C
c$$$      IF(IPULS.EQ.1) THEN
c$$$            FPRIM=EXP(-TAUTOT(I,J))*FMEAN(j,1,E1(J))
c$$$      ELSE
c$$$            FPRIM=0.
c$$$      ENDIF
c$$$C     ***********************************************************
c$$$C     *****
c$$$C
c$$$C     FOR SECOND AND HIGHER ITERATIONS THE DIFFUSE FLUX IS 
c$$$C     ALREADY CALCULATED
c$$$C     MEAN INTENSITY = PRIMARY FLUX + DIFFUSE FLUX
c$$$c!    IF TAU > 1 USE J = SOURCE FCN.
c$$$C     *****
c$$$C     ***********************************************************
c$$$      IF(N.EQ.1) THEN
c$$$      FD(I,J)=0.
c$$$c      IF(TAUTOT(I,J).GT.0.001.AND.N.EQ.1)  
c$$$c     &                  FL(2,J)=S(I,J)*(1.-exp(-tautot(i,j)))
c$$$C
c$$$c     OUTWARD ONLY APPROX ( TWICE AS MUCH AS IF INTERNAL ABSORPTION)
c$$$C
c$$$            IF(TAX.LT.10.) THEN      
c$$$                  TAUIK=0.
c$$$                  TAUIKM1=0.
c$$$                  DO IK=I,2,-1
c$$$C!!   ASSUME NO DILUTION AND THAT THERE IS NO ABSORPTION INTERNAL TO 
c$$$C           RADIUS I.
c$$$                        TAUIK=TAUIK+TAU(IK,J)
c$$$                        TAUIKM1=TAUIK-TAU(IK,J)
c$$$C
c$$$c     OUTWARD ONLY APPROX ( TWICE AS MUCH AS IF INTERNAL ABSORPTION)
c$$$C
c$$$                  DO KK=1,2
c$$$                    IF(KK.EQ.1) X=TAUIK
c$$$                    IF(KK.EQ.2) X=TAUIKM1
c$$$                    E2Q(KK)=1.
c$$$                    IF(X.GT.0.) THEN
c$$$                      IF(X.GT.0..AND.X.LT.0.1) THEN
c$$$                        E2Q(KK)=1.-X+X*X/2-X*(-LOG(X)-.577216+X)
c$$$                        IF(IH.EQ.1) E3Q(KK)=(1.D0-X+X*X/2.D0
c$$$     &                                     -X**3/6.D0-X*E2Q(KK))/2.
c$$$                      ELSE
c$$$                        E2Q(KK)=E2(X)
c$$$                        IF(IH.EQ.1) E3Q(KK)=(EXP(-X)-X*E2Q(KK))/2.
c$$$                      ENDIF
c$$$                    ENDIF
c$$$                  ENDDO
c$$$                        FD(I,J)=S(IK,J)*(-E2Q(1)+E2Q(2))
c$$$     &                                    +FD(I,J)
c$$$                        IF(IH.EQ.1) FHD(J)=S(IK,J)*(-E3Q(1)+E3Q(2))
c$$$     &                                    +FHD(J)
c$$$                  ENDDO
c$$$            ELSE
c$$$                  FD(I,J)=S(I,J)
c$$$                  IF(IH.EQ.1) FHD(J)=0.
c$$$            ENDIF
c$$$            FL(2,J)=FD(I,J)+FPRIM
c$$$      ELSE
c$$$            FL(2,J)=FPRIM+FD(I,J)
c$$$      ENDIF
c$$$      FH0(J)=FPRIM
c$$$      JMEAN(J)=FL(2,J)
c$$$      RETURN
c$$$      END

      SUBROUTINE IONABUND(XP)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameters.h'
      COMMON/ABU/XA(2,nio)
      COMMON/ABUN/AB(20)
      DIMENSION XP(nio)
      XP(1)=AB(1)*XA(2,1)
      DO K=2,3
         XP(K)=AB(2)*XA(2,K)
      enddo
      DO K=4,7
         XP(K)=AB(3)*XA(2,K)
      enddo
      DO K=8,13
         XP(K)=AB(4)*XA(2,K)
      enddo
      DO K=14,17
         XP(K)=AB(3)*XA(2,K)
      enddo
      DO K=18,24
         XP(K)=AB(5)*XA(2,K)
      enddo
      DO K=25,38
         XP(K)=AB(6)*XA(2,K)
      enddo
      DO K=39,41
         XP(K)=AB(7)*XA(2,K)
      enddo
      DO K=42,45
         XP(K)=AB(8)*XA(2,K)
      enddo
      DO K=46,49
         XP(K)=AB(9)*XA(2,K)
      enddo
      DO K=50,52
         XP(K)=AB(11)*XA(2,K)
      enddo
      DO K=53,54
         XP(K)=AB(10)*XA(2,K)
      enddo
      DO K=56,71
         XP(K)=AB(12)*XA(2,K)
      enddo
      DO  K=72,81
         XP(K)=AB(13)*XA(2,K)
      enddo
      DO K=82,83
         XP(K)=AB(14)*XA(2,K)
      enddo
      DO K=84,94
         XP(K)=AB(8)*XA(2,K)
      enddo
      DO K=95,111
         XP(K)=AB(14)*XA(2,K)
      enddo
      DO K=112,114
         XP(K)=AB(11)*XA(2,K)
      enddo
      RETURN
      END

      SUBROUTINE RADTRANOIq(J,TAX,S)
C     ***********************************************************
C     *****
C     FOR THE FIRST ITERATION (NITSPH=1) SOLVE THE EQUATION OF TRANSFER
C         IN THE OUTWARD ONLY APPROXIMATION. FOR THE NEXT ITERATIONS
C         (NITSPH > 1) CALCULATE THE CONTINOUS OPACITIES FOR THE SHELL.
C     WORKS BOTH IN OUT AND OUT IN
C     HERE OUT TO IN
C     USED FOR ALL CS-INT. CALC.
C     *****
C     ***********************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      SAVE
      COMMON/IND/I
      COMMON/CINOUT/INOUT,IPULS
      COMMON/CICS/RSHOCK,ics
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/ITSPH/NITSPH
      COMMON/LITER/N
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/RADIE/R(MD)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),ipar
      real*8 jmean
      COMMON/DTAU/JMEAN(NE1:NE2)
      DIMENSION S(MD,NE1:NE2)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      write(6,*)' i radtranoiq '
C
C     PRIMARY FLUX
C
      IF(IPULS.EQ.1.AND.ICS.NE.1) THEN
            FPRIM=EXP(-TAUTOT(I,J))*FMEAN(j,1,E1(J))
      ELSEIF(ICS.EQ.1.AND.INOUT.EQ.0) THEN
C           ASSUMES ABSORPTION OF FLUX FROM SHELL BY THE GAS INTERIOR
C             TO R.
C!!           NOT GOOD FOR SECOND AND HIGHER ITERATIONS
            FPRIM=E2(TAUTOT(I,J))*FMEAN(j,1,E1(J))
      ELSE
            FPRIM=0.
      ENDIF
C     ***********************************************************
C     *****
C
C     FOR SECOND AND HIGHER ITERATIONS THE DIFFUSE FLUX IS 
C     ALREADY CALCULATED
C     MEAN INTENSITY = PRIMARY FLUX + DIFFUSE FLUX
c!    IF TAU > 1 USE J = SOURCE FCN.
C     *****
C     ***********************************************************
      IF(N.EQ.1) THEN
      FL(2,J)=0.
      IF(TAUTOT(I,J).GT.0.001.AND.N.EQ.1)  
     &                  FL(2,J)=S(I,J)*(1.-exp(-tautot(i,j)))
c!!!!!!!!!!
c      goto 1234
            IF(TAX.LT.10.) THEN      
                  FL(2,J)=0.
                  TAUIK=0.
                  TAUIKM1=0.
                  DO IK=I,2,-1
C!!   ASSUME NO DILUTION AND THAT THERE IS NO ABSORPTION INTERNAL TO 
C           RADIUS I.
                        TAUIK=TAUIK+TAU(IK,J)
                        TAUIKM1=TAUIK-TAU(IK,J)
C
c     OUTWARD ONLY APPROX ( TWICE AS MUCH AS IF INTERNAL ABSORPTION)
C
                        FL(2,J)=S(IK,J)*(-E2(TAUIK)+E2(TAUIKM1))
     &                                    +FL(2,J)
                  ENDDO
800               CONTINUE
            ELSE
                  FL(2,J)=S(I,J)
            ENDIF
1234  continue
            FL(2,J)=FL(2,J)+FPRIM
      ELSE
            FL(2,J)=FPRIM+FD(I,J)
      ENDIF
      JMEAN(J)=FL(2,J)
      RETURN
      END

      SUBROUTINE RADTRANPULS(IH,J,TAX,S)
C     ***********************************************************
C     *****
C     FOR THE FIRST ITERATION (NITSPH=1) SOLVE THE EQUATION OF TRANSFER
C         IN THE OUTWARD ONLY APPROXIMATION. FOR THE NEXT ITERATIONS
C         (NITSPH > 1) CALCULATE THE CONTINOUS OPACITIES FOR THE SHELL.
C     WORKS BOTH IN OUT AND OUT IN
C     HERE OUT TO IN
C     USED FOR ALL CS-INT. CALC.
C     IF IH = 1 CALCULATE FLUX & MEAN INTENSITY
C     IF IH = 1 CALCULATE MEAN INTENSITY
C     *****
C     ***********************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameters.h"
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      SAVE
      COMMON/IND/I
      COMMON/CINOUT/INOUT,IPULS
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/ITSPH/NITSPH
      COMMON/LITER/N
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/RADIE/R(MD)
      COMMON/DXA/DX(MD)
      COMMON/PHY/DEN(MD)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),ipar
      real*8 jmean
      COMMON/DTAU/JMEAN(NE1:NE2)
      COMMON/DIFFP/FH0(NE1:NE2),FHD(NE1:NE2)
      DIMENSION S(MD,NE1:NE2),E2Q(2),E3Q(2)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
C
C     PRIMARY FLUX
C
      IF(IPULS.EQ.1) THEN
         FPRIM=EXP(-TAUTOT(I,J))*FMEAN(j,1,E1(J))
      ELSE
         FPRIM=0.
      ENDIF
      if(j==jmin) then
         fptot=0.
         fltot=0.
         emtot=0.
      endif
      IF(IH.EQ.1) GOTO 888
C     ***********************************************************
C     *****
C
C     FOR SECOND AND HIGHER ITERATIONS THE DIFFUSE FLUX IS 
C     ALREADY CALCULATED
C     MEAN INTENSITY = PRIMARY FLUX + DIFFUSE FLUX
c!    IF TAU > 1 USE J = SOURCE FCN.
C     *****
C     ***********************************************************
C     ***********************************************************
C     *****
C     USE THE HATCHETT, BUFF & MCCRAY TRANSPORT EQ. FOR SPHERICAL
C     GEOMETRY.
C     *****
C     ***********************************************************

      IF(N.EQ.1) THEN
        IF(TAX.LT.1.E-5) THEN
          EXTAU=1.-TAX+TAX**2/2.
        ELSEIF(TAX.LT.700.) THEN
          EXTAU=EXP(-TAX)
        ELSEIF(TAX.GE.700.) THEN
          EXTAU=0.
       ENDIF
c FL = J       
        FL(2,J)=(FL(1,J)*EXTAU+S(I,J)*(1.-EXTAU))
     &                  /(1.+(R(I)-R(I-1))/R(I-1))**2.
c!! test only!!
c        FL(2,J)=FL(1,J)*EXTAU
        FD(I,J)=FL(2,J)-FPRIM
        dil=1./(1.+(R(I)-R(I-1))/R(I-1))**2.
      ELSE
C
C       FOR SECOND AND HIGHER ITERATIONS THE DIFFUSE FLUX IS 
C       ALREADY CALCULATED
C       MEAN INTENSITY = PRIMARY FLUX + DIFFUSE FLUX
C
        FL(2,J)=FPRIM+FD(I,J)
      ENDIF
 888  continue
      FHD(J)=FD(I,J)
      FH0(J)=FPRIM
      JMEAN(J)=FL(2,J)
      fptot=fptot+16.*pi**2*r(i)**2*fprim*(e(j+1)-e(j))
      emtot=emtot+16.*pi**2*r(i)**2*dx(i)*den(i)*em(i,j)*(e(j+1)-e(j))
      fltot=fltot+16.*pi**2*r(i)**2*fl(2,j)*(e(j+1)-e(j))

      RETURN
      END


