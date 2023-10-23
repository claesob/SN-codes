 99   subroutine enum_cf(cfel,cfion,elcf)
      IMPLICIT REAL*8(A-H,O-Z)
C     ENUMERATION OF THE IONS:
C
C      1 = H I      2 = HE I     3 = HE II    4 = O VI     5 = O VII
C      6 = O VIII   7 = O V      8 = C III    9 = C IV    10 = C V
C     11 = C VI    12 = C I     13 = C II    14 = O I     15 = O II
C     16 = O III   17 = O IV    18 = N I     19 = N II    20 = N III
C     21 = N IV    22 = N V     23 = N VI    24 = N VII   25 = SI I
C     26 = SI II   27 = SI III  28 = SI IV   29 = SI V    30 = SI VI
C     31 = SI VII  32 = SI VIII 33 = SI IX   34 = SI X    35 = SI XI
C     36 = SI XII  37 = SI XIII 38 = SI XIV  39 = MG I    40 = MG II
C     41 = MG III  42 = FE I    43 = FE II   44 = FE III  45 = FE IV
C     46 = AL I    47 = AL II   48 = AL III  49 = AL IV   50 = CA I
C     51 = CA II   52 = CA III  53 = NA I    54 = NA II   55 = NA III
C     56 = S I     57 = S II    58 = S III   59 = S IV    60 = S V
C     61 = S VI    62 = S VII   63 = S VIII  64 = S IX    65 = S X
C     66 = S XI    67 = S XII   68 = S XIII  69 = S XIV   70 = S XV
C     71 = S XVI   72 = NE I    73 = NE II   74 = NE III  75 = NE IV 
C     76 = NE V    77 = NE VI   78 = NE VII  79 = NE VIII 80 = NE IX 
C     81 = NE X    82 = AR I    83 = AR II   84 = FE V    85 = FE VI
C     86 = FE VII  87 = FE VIII 88 = FE IX   89 = FE X    90 = FE XI
C     91 = FE XII  92 = FE XIII 93 = FE XIV  94 = FE XV   95 = O III 1D
c     96 = AR III  97 = AR IV   98 = AR V    99 = AR VI   100 = AR VII
c    101 = FE XVI 102 = FE XVII 103 = FE XVIII 104 = FE XIX 105 = FE XX
c    106 = FE XXI 107 = FE XXII 108 = FE XXIII 109 = FE XXIV 110 = FE XXV
c    111 = FE XXVI 112 = Ca IV  113 Ca V       114 = Ca VI
CK  NEW D

c new numbering
c
c 1 = H, 2 = He, 3 = C, 4 = N, 5 = O, 6 = Ne, 7 = Na, 8 = Mg, 
c 9 = Al, 10 = Si, 11 = S, 12 = Ar, 13 = Ca, 14 = Fe

c      include 'param'
c      common/ionx/xion(md,14,27)
      COMMON/ABUN/AB(20)
      common/abl/abn(15)
      dimension iabcf(14)
      integer cfel(120),cfion(120),elcf(14,27)

c
c conversion from logical to cf  abundance enumeration      
c
c                
      data iabcf/1,2,4,5,3,13,10,7,9,6,12,14,11,8/
c      data icfab/1,2,5,5,3,13,10,7,9,6,12,14,11,8/

c translate from CF radom numbering to ordered      
c H I
      cfel(1)=1
      cfion(1)=1
c He I-II      
      cfel(2)=2
      cfion(2)=1
      cfel(3)=2
      cfion(3)=2
c     O VI-VIII, V
      do i=4,6
         cfel(i)=5
         cfion(i)=i+2
      enddo
      cfel(7)=5
      cfion(7)=5
c C III-VI
      do i=8,11
         cfel(i)=3
         cfion(i)=i-5
      enddo      
c C I-II            
      cfel(12)=3
      cfion(12)=1
      cfel(13)=3
      cfion(13)=2      
c O I-IV
      do i=14,17
         cfel(i)=5
         cfion(i)=i-13
      enddo
c N I-7
      do i=18,24
         cfel(i)=4
         cfion(i)=i-17
      enddo
c Si I-XIV
      do i=25,38
         cfel(i)=10
         cfion(i)=i-24
      enddo
c Mg      
      do i=39,41
         cfel(i)=8
         cfion(i)=i-38
      enddo
c Fe I-V
      do i=42,45
         cfel(i)=14
         cfion(i)=i-41
      enddo
c Al I-IV
      do i=46,49
         cfel(i)=9
         cfion(i)=i-45
      enddo
c     Ca I-III
      do i=50,52
         cfel(i)=13
         cfion(i)=i-49
      enddo
c     Na
      do i=53,55
         cfel(i)=7
         cfion(i)=i-52
      enddo
c S
      do i=56,71
         cfel(i)=11
         cfion(i)=i-55
      enddo
c Ne
      do i=72,81
         cfel(i)=6
         cfion(i)=i-71
      enddo      
c Ar I-II
      do i=82,83
         cfel(i)=12
         cfion(i)=i-81
      enddo
c     Fe V-XV
      do i=84,94
         cfel(i)=14
         cfion(i)=i-79
      enddo

c     O III 1D Obs! ion=33!
      cfel(95)=5
      cfion(95)=33
c Ar III-VII
      do i=96,100
         cfel(i)=12
         cfion(i)=i-93
      enddo
c     Fe XVI-XXVI
      do i=101,111
         cfel(i)=14
         cfion(i)=i-85
      enddo
c     Ca IV-VI
      do i=112,114
         cfel(i)=13
         cfion(i)=i-108
      enddo      
c translate from ordered to old CF
      do iel=1,14
         do ion=1,27
            do k=1,114
               if(cfel(k)==iel.and.cfion(k)==ion) then
                  elcf(iel,ion)=k
               endif
            enddo
         enddo
      enddo


c  CF  normal

c  1 = 1,

c  2 = 2,

c  3 = 5,   O

c  4 = 3,   C

c  5 = 4,   N

c  6 = 10,  Si

c  7 = 8,   Mg

c  8 = 14,  Fe

c  9 = 9,   Al

c 10 = 7,   Na

c 11 = 13,  Ca

c 12 = 11,  S

c 13 = 6,   Ne

c 14 0 12   Ar


      xe = 0.
      do iel=1,14
         iab = iabcf(iel)
         abn(iel)=ab(iab)
      enddo

      return
      end


      subroutine cfenum
c
c translate from normal system to cf enumeration
c note that first index is ionization stage and second is element
c
      integer inum
      common/cfum/inum(27,14)
      data inum/1,26*0,2,3,25*0,12,13,8,9,10,11,21*0,
     &         18,19,20,21,22,23,24,20*0,14,15,16,17,7,4,5,6,19*0,
     &         72,73,74,75,76,77,78,79,80,81,17*0,
     &         53,54,55,24*0,39,40,41,24*0,46,47,48,49,23*0,
     &         25,26,27,28,29,30,31,32,33,34,35,36,37,38,13*0,
     &         56,57,58,59,60,61,62,63,64,65,66,67,68,69,71,12*0,
     &         82,83,96,97,98,99,100,20*0,50,51,52,112,113,114,21*0,
     &     42,43,44,45,84,85,86,87,88,89,90,91,92,93,94,101,102,103,
     &     104,105,106,107,108,109,110,111,0/

      return
      end

      subroutine conv_xa_xion
c     convert from xa() to xion()
      include "parameters.h"
      implicit real*8(a-h,o-z)
      COMMON/ABU/XA(2,nio)
      common/cfum/inum(27,14)
      COMMON/IND/IK
c      include 'PARAM'
      common/ionx/xion(md,14,27)
      COMMON/COLDnew/COLTOTnew(14,27),COLTOTTnew(14,27)
      COMMON/COLD/OPT(nio),COLTOT(nio),COLTOTT(nio),SRED(NFEL+100)
      do iel=1,2
         do ion=1,27
            if(inum(ion,iel) .ne. 0) then
               xion(ik,iel,ion)=xa(2,inum(ion,iel))
            endif
         enddo
      enddo
      do iel=1,14
         do ion=1,27
            if(inum(ion,iel) .ne. 0) then
               icfnum=inum(ion,iel)
               coltotnew(iel,ion)=coltot(icfnum)
               coltottnew(iel,ion)=coltott(icfnum)
               ionmax=ion
            endif
         enddo
      enddo
      return
      end

      subroutine convrtel(iel,nelem)

c convert from elem numbering to atomic number

c eg iel=5 and nelem=8 for O

      implicit real*8(a-h,o-z)
      if(iel.ge.3.and.iel.le.5) then
         nelem=iel+3
      elseif(iel.ge.6.and.iel.le.10) then
         nelem=iel+4
      elseif(iel.eq.11) then
         nelem=iel+5
      elseif(iel.eq.12) then
         nelem=iel+6
      elseif(iel.eq.13) then
         nelem=iel+7
      elseif(iel.eq.14) then
         nelem=iel+12
      endif
      return
      end


      DOUBLE PRECISION FUNCTION RAD(TE,XELEC,IFPOP)
C!!   CHECK CHARGE TRANSFER RECOMB.ALL CT SHOULD BE DIVIDED BY DELA(MI)?
C!!   SEE FE !
      IMPLICIT REAL*8(A-H,O-Z)
      real*4 thpop,th1,th2,tr1,tr2,tr
      real*4 trat1,trat2,trat,tem1,tem2,temi,tsp1,tsp2,tspec
	real*4 tpophe1,tpophe2,tpophe,tpopo1,tpopo2,tpopo,tpopca1,
     & 	tpopca2,tpopca,tpopfe1,tpopfe2,tpopfe
      real*8 jlah,j2gh,jrech,jemh
      character*12 lev
      character*80 lab
      include "parameters.h"
c      PARAMETER (MD=300,MDP1=MD+1)
c      PARAMETER(NFEL=1000)
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
C     ****************************************************************
C     *****
C     THIS SUBROUTINE CALCULATES ALL RADIATIVE AND COLLISIONAL RATES
C     AS WELL AS THE IONIZATION STRUCTURE (EXCEPT FOR SILICON AND
C     MAGNESIUM ), FOR A GIVEN TEMPERATURE TE.
C
C
C     *****
C     ****************************************************************
C     *****
c     old numbering of elements
      
C     AB(1)=ABUNDANCE OF HYDROGEN (ALLEN (73) = 1.0 )
C     AB(2)=    D:O   OF HELIUM (ALLEN (73) = 8.5E-2)
C     AB(3)=    D:O   OF OXYGEN (ALLEN (73) = 6.6E-4)
C     AB(4)=    D:O   OF CARBON (ALLEN (73) = 3.3E-4)
C     AB(5)=    D:O   OF NITROGEN (ALLEN (73) = 9.1E-5)
C     AB(6)=    D:O   OF SILICON (ALLEN (73) = 3.3E-5)
C     AB(7)=    D:O   OF MAGNESIUM (ALLEN (73) = 2.6E-5)
C     AB(8) =    D:O   OF IRON (ALLEN (73) = 4.0E-5)
C     AB(9) =    D:O   OF ALUMINUM (ALLEN (73) =  ???? )
C     AB(10) =    D:O   OF SODIUM (ALLEN (73) =  ???? )
C     AB(11) =    D:O   OF CALCIUM (ALLEN (73) =  ???? )
C     AB(12) =    D:O   OF SULPHUR (ALLEN (73) =  ???? )
C     AB(13) =    D:O   OF NEON (ALLEN (73) =  ???? )
C     AB(14) =    D:O   OF ARGON (ALLEN (73) =  ???? )
C
C     ENUMERATION OF THE IONS:
C
C      1 = H I      2 = HE I     3 = HE II    4 = O VI     5 = O VII
C      6 = O VIII   7 = O V      8 = C III    9 = C IV    10 = C V
C     11 = C VI    12 = C I     13 = C II    14 = O I     15 = O II
C     16 = O III   17 = O IV    18 = N I     19 = N II    20 = N III
C     21 = N IV    22 = N V     23 = N VI    24 = N VII   25 = SI I
C     26 = SI II   27 = SI III  28 = SI IV   29 = SI V    30 = SI VI
C     31 = SI VII  32 = SI VIII 33 = SI IX   34 = SI X    35 = SI XI
C     36 = SI XII  37 = SI XIII 38 = SI XIV  39 = MG I    40 = MG II
C     41 = MG III  42 = FE I    43 = FE II   44 = FE III  45 = FE IV
C     46 = AL I    47 = AL II   48 = AL III  49 = AL IV   50 = CA I
C     51 = CA II   52 = CA III  53 = NA I    54 = NA II   55 = NA III
C     56 = S I     57 = S II    58 = S III   59 = S IV    60 = S V
C     61 = S VI    62 = S VII   63 = S VIII  64 = S IX    65 = S X
C     66 = S XI    67 = S XII   68 = S XIII  69 = S XIV   70 = S XV
C     71 = S XVI   72 = NE I    73 = NE II   74 = NE III  75 = NE IV 
C     76 = NE V    77 = NE VI   78 = NE VII  79 = NE VIII 80 = NE IX 
C     81 = NE X    82 = AR I    83 = AR II   84 = FE V    85 = FE VI
C     86 = FE VII  87 = FE VIII 88 = FE IX   89 = FE X    90 = FE XI
C     91 = FE XII  92 = FE XIII 93 = FE XIV  94 = FE XV   95 = O III 1D
c     96 = AR III  97 = AR IV   98 = AR V    99 = AR VI   100 = AR VII
c    101 = FE XVI 102 = FE XVII 103 = FE XVIII 104 = FE XIX 105 = FE XX
c    106 = FE XXI 107 = FE XXII 108 = FE XXIII 109 = FE XXIV 110 = FE XXV
c    111 = FE XXVI 112 = Ca IV  113 Ca V       114 = Ca VI      
C     
C
C     *****
C     *****
C     ****************************************************************
C     *****
      COMMON/EMHY/RECEM(5,NL),TWOPH,TWOPHHEI,COHI,PHI,PHIT,PO,POT
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      common/heiires/x2he,he2coll12,he2jla,he2j2g,he2jrec(4),
     &                                                he2jem(7,4)
      common/hres/x2h,coll12h,jlah,j2gh,jrech(4),jemh(7,4)
      COMMON/RECMET/recappc,recappo,recappsi,recapps,recappfe
      COMMON/GSREC/ALGS(nio),ALEX(nio),ALTOT(nio),RECEX(nio),RECGS(nio)
      COMMON/CSEX/CSLYA
      COMMON/HYPOP/XNH(NL)
      COMMON/INT/FL(2,NE1:NE2),SI(nio,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/SECX/CSEC(20),DCSDX(20),CISEC(20),DCISDX
      COMMON/SPOT/OTSP(7)
      COMMON/NLEV/e00,NION,NHY,NP1H,NMA,NMI
      COMMON/GAMPAR/GAMLUM,GAELH,GAHE
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/COL/RTE(4),FF,HY,HE,C131,C142,C151,C161,C231,C241,COH1,
     &COHEcollion,C351,C321,C331,C341,C222,C232,C221,C332,C441
      COMMON/EQUIV/RADF,FR(2,100),W(100),CIN(100),WLI(100+NFEL),FB(100)
      COMMON/COLD/OPT(nio),COLTOT(nio),COLTOTT(nio),SRED(NFEL+100)
      COMMON/EQUIH/FRH(2,500),WEH(500),SR(NFEL),WOBS(NL,NL)
      COMMON/PHOTOHEAT/PHEAT(5,NL),PHEATT(5,NL)
      COMMON/HEA/CO,PH,PHEO,PHEI,POX(8),PC(6),PN(7),PMG,PSI,PFE
      COMMON/TEQQ/TEZ
      COMMON/REHEL/REHE21,REHE22,REC31,REN41
      COMMON/PHY/DEN(MD)
      COMMON/IND/IK
      COMMON/PHQ/ZEA(nio),GEA(nio),ZKA(nio)
      COMMON/ABU/XA(2,nio)
      COMMON/COSUL/COSUL(18)
      COMMON/REC/AL2(7)
      COMMON/ELEC/DEL(MD)
      COMMON/ABUN/AB(20)
      common/abl/abn(15)
      COMMON/RADIE/R(MD)
      COMMON/NHY/NH
      COMMON/HPOP/XNQ(6,NL),XN1,XN2,XN3
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &     BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/HRAT/RHYD,ZEL,XEL,HYR,HEAT,COOL
      COMMON/REC1/RE1(7),RHE2B
      COMMON/A14/CQ(NL,NL),CIQ(NL),GQ(NL),EI(NL),AQ(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTEQ(NL),PHN(NL),DRDT(NL),RECT(5,NL),PHET(NL)
      COMMON/RECAL/RECO(NL)
      COMMON/A19/EMCA(NL,NL),ESCCA(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      COMMON/A31/TOPI(10,NL,NL),EMHY(NL,NL)
      COMMON/IOINF/IOSPEC,IOLEVEL,IOEMISS,IOTERM,IOION,IOELITER
      COMMON/IOPAR/IOSHELL,IORAD1,IORAD2,IOFLUX
      COMMON/MULTI/IHMUL,IHEMUL,IOMUL,ICAMUL,IFEMUL
      COMMON/LOWION/ILOWION
      COMMON/CINOUT/INOUT,IPULS
      COMMON/CICS/RSHOCK,ics
      COMMON/FELEV/NFE,nfei
      COMMON/WFE/WOBFE(NL,NL),IQPOPH
      common/w/wlin(501),wlij(NL,NL),taul(501)
      COMMON/FESAVE/INIFE(4),BOLFES(4,NL),TOLDFES(4)
      DIMENSION FI1(7),YI(7),DELA(MD),FD(MD),YEX1(7),CEX(500)
      DIMENSION PF(26),ZE(nio),ZK(nio),GE(nio),COF(5)
      DIMENSION IA(28),rech(10,3),g(100),emfe4(20,20)
      DIMENSION XN(30),ZC(7),XO(30),XR(30),ZA(8),XNA(30),XSU(30),XAR(30)
      dimension OXR(10),FBO(30)
      parameter (mh=2,mhe=3,mc=7,mn=8,mo=9,mne=11,mna=12,mmg=13,
     &mal=14,msi=15,ms=17,mar=19,mca=21,mfe=27,mni=29,
     &mmat=mh+mhe+mc+mn+mo+mne+mna+mmg+mal+msi+ms+mar+mca+mfe+mni,
     &     ispecy=15)
      COMMON/ABC/AL(16)
      common/abc1/alc(mc),aln(mn),alo(mo),alne(mne),alna(mna),
     &almg(mmg),alal(mal),alsi(msi),alsu(ms),alar(mar),alca(mca),
     &     alfe(mfe),alni(mni)
      common/abc2/coh(1),cohe(2),coc(mc-1),con(mn-1),coo(mo-1),
     &cone(mne-1),cona(mna-1),comg(mmg-1),coal(mal-1),cosi(msi-1),
     &     cosu(ms-1),coar(mar-1),coca(mca-1),cofe(mfe-1),coni(mni-1)
      common/line_em/cinx(14,26,401),taulinex(14,26,401),wlix(14,26,401)
     &     ,ilabcfx(14,26,401)
      common/ionx/xion(md,14,27)
      common/ionx_old/xion_old(14,27)
      DIMENSION XC(30),ZB(6),XF(30),XSi(30),XM(30),XAL(30),XCA(30),
     &XNE(30),LHEI(16),WLHEI(16),LFEII(16),WLFEII(16),ZS(16),
     &LHI(10),WLHI(10),X(nio),CHI(nio)
      DIMENSION IAB(nio),ZNEO(10),XI(15),ZSUL(16),A21HNU(5)
      dimension omar(5,5),omca4(3)
      common/rec_coll/ alrec(30,30),collion(30,30)
      integer cfel(120),cfion(120),elcf(14,27)
      common/kmaxpop/kmaxp(14,27)
      common/wlforb3/wl21,wl31,wl32
      common/kmax/kmax(14,27)
      COMMON/datar5/Gar(5),Ear(5),aar(5,5),WLar(5,5)
      common/collar5/omt(12),omfitar(20,20,15),tefitar(15)
      common/ichia/ipopch
      COMMON/SILEV/NSIII,nsi1
      common/gsreco/al_gs(26,26,10)
      common/init_ion/inition
      common/ca4/teca4(50),coll_ca4(3,50)
      common/raug/zion(30,27,7)
      real*8 line_cool
      common/line_c/line_cool(14,27)
      parameter(nchi22=20)
      PARAMETER (nups=65)
      common/chianti_2022_cs/upsil(nchi22,nups,nups)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DATA CSAHA/2.0708E-16/
      DATA A21HNU/1.,7.663E-3,1.301E-4,3.390E-5,1.144E-5/
      DATA LHEI/7,5,17,38,9,14,25,32,49,21,42,28,35,52,44,63/
      DATA WLHEI/584.,10830.,3889.,3188.,20581.,7065.,5876.,4713.,
     &4471.,42944.,12547.,186233.,21120.,17002.,19543.,18648./
      DATA WLHI/1216.,1025.,6563.,972.,4861.,18751.,949.,4340.,
     &12818.,40512./
      DATA LFEII/7,5,17,38,9,14,25,32,49,21,42,28,35,52,44,63/
c     enumeration of elements with respect to ions
      DATA IAB/1,2,2,3,3,3,3,4,4,4,
     &     4,4,4,3,3,3,3,5,5,5,
     &     5,5,5,5,6,6,6,6,6,6,
     &     6,6,6,6,6,6,6,6,7,7,
     &     7,8,8,8,8,9,9,9,9,11,
     &     11,11,10,10,10,12,12,12,12,12,
     &     12,12,12,12,12,12,12,12,12,12,
     &     12,13,13,13,13,13,13,13,13,13,
     &     13,14,14,8,8,8,8,8,8,8,
     &     8,8,8,8,3,14,14,14,14,14,
     &     8,8,8,8,8,8,8,8,8,8,
     &     8,9*20/
      integer z2cf(26)      
      data z2cf/1,2,0,0,0,3,4,5,0,6,7,8,9,10,0,11,0,12,0,13,0,0,0,0,0,14/
      real*8 ar2_te(11),ar2_coll(11)
      data ar2_te/1.58E+03,2.51E+03,3.98E+03,6.31E+03,1.00E+04,1.58E+04,
     &     2.51E+04,3.98E+04,6.31E+04,1.00E+05,1.58E+05/
      data ar2_coll/2.48,2.54,2.63,2.77,2.93,3.09,3.19,3.2,3.13,2.97,2.71/
      real*8 ar6_te(43),ar6_coll(43)
      data ar6_te/2.00E+01,5.00E+01,1.00E+02,2.00E+02,3.00E+02,4.00E+02,
     &     5.00E+02,6.00E+02,7.00E+02,8.00E+02,9.00E+02,1.00E+03,1.50E+03,
     &     2.00E+03,2.50E+03,3.00E+03,3.50E+03,4.00E+03,4.50E+03,5.00E+03,
     &     5.50E+03,6.00E+03,6.50E+03,7.00E+03,7.50E+03,8.00E+03,8.50E+03,
     &     9.00E+03,9.50E+03,1.00E+04,1.10E+04,1.20E+04,1.30E+04,1.40E+04,
     &     1.50E+04,1.60E+04,1.80E+04,2.00E+04,2.40E+04,2.80E+04,3.20E+04,
     &     3.60E+04,4.00E+04/		
      data ar6_coll/3.04,3.02,3.0,2.97,2.95,2.96,2.97,3.0,3.03,3.06,3.08,
     &     3.11,3.23,3.4,3.62,3.86,4.11,4.35,4.57,4.78,4.96,5.13,5.27,5.4,5.51,
     &     5.61,5.7,5.77,5.84,5.9,6.0,6.08,6.14,6.19,6.23,6.26,6.3,6.33,
     &     6.36,6.36,6.36,6.35,6.33/

c     no h or he in this version!
      ab(1)=0.
      ab(2)=0.
      xh=0.
      xhe=0.
      if(te <= 100.) then
         write(6,*)' Te low ',te
c         te=101.
      endif

      call enum_cf(cfel,cfion,elcf)
      ss=rec_colli(ik,te,dens,xelec)
      TLIMO=100.
      TLIMO=10.
c!!   use popca to very low T!
      TLIMCA=100.
      TEZ=TE
      T4=TE/1.E4
      T3=TE/1.E3
      XEL=XELEC
      DEL(IK)=XEL
      MIMAX=50
9287  FORMAT(' RAD ',5E11.4)
      DO I=1,8
         XO(I)=0.0E0
         ZA(I)=0.0E0
      enddo
      DO I=1,9
         XR(I)=0.0E0
         XO(I)=0.0E0
      enddo
      do iel=1,14
         do ion=1,27
            line_cool(iel,ion)=0.
         enddo
      enddo
C     
C     DIVIDE ALL IONIZATION AND HEATING RATES BY THE DENSITY
C
      DO I=1,114
         ZE(I)=ZEA(I)/DEN(IK)
         ZK(I)=ZKA(I)/DEN(IK)
         GE(I)=GEA(I)/DEN(IK)
      enddo

C
C     FIRST CALCULATE THE IONIZATION AND RECOMBINATION RATES
C     ***************************************************************
C     *****
C     RECOMBINATION COEFF. FOR HYDROGEN LIKE IONS AND HELIUM
C     SEATON (1959) AND GOULD AND TAKHUR (1966)
C     RATES ARE ADJUSTED TO AGREE WITH OSTERBROCK (74) FOR T=1E4
C     *****
C     ***************************************************************

      YI(4)=64.*YI(1)
      YI(5)=36.*YI(1)
      YI(6)=49.*YI(1)
      YI(7)=196.*YI(1)

      MQ=0
      MI=1
C
C     FIRST ESTIMATE OF ELECTRON FRACTION.
C
 776  CONTINUE

      IF(IK.EQ.2.AND.IK.NE.IKOLD) THEN
        DELA(1)=DEL(1)
        IKOLD = 2
      ELSE
        IF(IK.EQ.IKOLD.OR.NIT.GT.1) THEN
          DELA(1) = XELEC
        ELSE
          DELA(1) = DEL(IK-1)
          IKOLD = IK
        ENDIF
      ENDIF
      DEEL=DELA(MI)
      CALL SECEX(DELA(MI))
c$$$C     ***************************************************************
c$$$C     *****
c$$$C     CALCULATE GOULDS CORRECTION FOR RADIATIVE REC. TO EXC. STATES
c$$$C     *****
c$$$  C     ***************************************************************
      FIR31=FIR(3,1.d0,TE)
C     ***************************************************************
C     *****
C     IONIZATION EQUIL. FOR OXYGEN INCLUDING AUGER IONIZATIONS
C     *****
C     ***************************************************************
c$$$C
c$$$C     OXYGEN RECOMBINATION COEFF. INCLUDING DIEL. REC.(ALD. AND PEQ.)
c$$$C     AND CHARGE TRANSFER.
c$$$      ALO(9)=AL(4)
c$$$       CALL DIEL(8.6d-2,7.0d6,0.2d0,1.3d6,TE,ALD)
c$$$      SUPR=1.
c$$$c!!      IF(DEN(IK).GT.1.E9) SUPR=0.26
c$$$       ALD=ALD*SUPR
c$$$      ALO(8)=4.1E-11/T4**.742+ALD
c$$$      CALL DIEL(1.1d-1,6.2d6,0.2d0,9.5d5,TE,ALD)
c$$$      ALD=ALD*SUPR
c$$$      ALO(7)=2.3E-11/T4**.802+ALD
c$$$      CALL DIEL(7.1d-3,1.3d5,3.2d0,8.0d5,TE,ALD)
c$$$      CALL DIELB(-2.8425d0,.2283d0,4.0407d1,-3.4956d0,1.7558d0,T4,ALDB)
c$$$      ALD=(ALD+ALDB)*SUPR
c$$$      ALO(6)=1.2E-11/T4**.779+ALD
c$$$      CALL DIEL(1.7d-2,2.2d5,2.0d0,5.9d5,TE,ALD)
c$$$      CALL DIELB(6.1d-3,.2269d0,3.2142d1,1.9939d0,-6.46d-2,T4,ALDB)
c$$$      ALD=(ALD+ALDB)*SUPR
c$$$C
c$$$C     C-T BD 80
c$$$C
c$$$C*********
c$$$      XCT=XH
c$$$      CT=((.16+.096*T4**1.14)*1.E-9*XH+6.5E-10*XHE)/DELA(MI)
c$$$C     ALRAD FROM GOULD (1978)
c$$$      ALTRAD=1.03E-11/T4**.5
c$$$      ALGS(17)=0.314*ALTRAD
c$$$      ALEX(17)=ALTRAD*(1.-0.314)*FIR31
c$$$      ALRAD=ALGS(17)+ALEX(17)     
c$$$      ALO(5)=ALRAD+ALD+CT
c$$$      CALL DIEL(2.8d-3,1.8d5,6.d0,9.1d4,TE,ALD)
c$$$      CALL DIELB(0.0d0,21.88d0,1.6273d1,-.7020d0,1.1899d0,T4,ALDB)
c$$$       ALD=(ALD+ALDB)*SUPR
c$$$C
c$$$C     C-T (BHD 80)
c$$$C
c$$$      CT=(.45+8.18*T4**.47)*1.E-9*XH/DELA(MI)
c$$$      CTHE=1.E-9*XHE/DELA(MI)
c$$$C     ALRAD FROM GOULD (1978)
c$$$      ALTRAD=5.43E-12/T4**.5
c$$$      ALGS(16)=0.350*ALTRAD
c$$$      ALEX(16)=ALTRAD*(1.-0.350)*FIR31
c$$$      ALRAD=ALGS(16)+ALEX(16)     
c$$$      ALO(4)=ALRAD+ALD+CT+CTHE
c$$$      CALL DIEL(1.4d-3,1.7d5,3.3d0,5.8d4,TE,ALD)
c$$$      CALL DIELB(-3.6d-3,.7519d0,1.5252d0,-8.38d-2,.2769d0,T4,ALDB)
c$$$c!!      IF(DEN(IK).GT.1.E9) SUPR=0.20
c$$$c!!      IF(DEN(IK).LT.1.E9) SUPR=1.000
c$$$       ALD=(ALD+ALDB)*SUPR
c$$$C
c$$$C     C-T (BHD 80)
c$$$C
c$$$      XH=XA(2,1)*AB(1)
      CT=(.28+.49*T4**.61)*1.E-9*XH/DELA(MI)
c$$$      CTHE=T4*2.E-10*XHE/DELA(MI)
c$$$      ALTRAD=2.05E-12/T4**.5
c$$$      ALGS(15)=0.346*ALTRAD
c$$$      ALEX(15)=ALTRAD*(1.-0.346)*FIR31
c$$$      ALRAD=ALGS(15)+ALEX(15)     
c$$$      ALO(3)=ALRAD+ALD+CT+CTHE
      ALCTO3=CT
c$$$      CALL DIEL(1.4d-3,1.7d5,2.5d0,1.3d5,TE,ALD)
c$$$      CALL DIELB(0.0d0,2.28d-2,6.59d-2,3.49d-2,.5334d0,T4,ALDB)
c$$$       ALD=(ALD+ALDB)*SUPR
c$$$C
c$$$C     C-T FS 71
c$$$C
c$$$      CTH=(.66+.25*T4**.5)*1.E-9*XH/DELA(MI)
       CTCA=7.6E-10*SQRT(TE/300.)
       CTCA=CTCA*XCA(1)*AB(11)/DELA(MI)
       ALCTO2=CTCA
c$$$      ALTRAD=3.31E-13/T4**.5
c$$$      ALGS(14)=0.260*ALTRAD
c$$$      ALEX(14)=ALTRAD*(1.-0.260)*FIR31
c$$$      ALRAD=ALGS(14)+ALEX(14)     
c$$$      ALO(2)=ALRAD+ALD+ALCTO2
c$$$      ALTOT(14)=ALO(2)-ALCTO2


      IF(ILOWION.EQ.1) KMAXi=2
      IF(ILOWION.EQ.0) KMAXi=8

      ibugion=0
      if(ibugion.eq.1) then
         do l=1,8
            if(l<=4) lq=13+l
            if(l==5) lq=7
            if(l>=6) lq=l-2
        enddo
      endif
c     *************
c     Oxygen ion eqiul. incl. Auger. Include ion. of 1D level (see above!)
c     *************

      call augion(8,9,dela(mi),xo)

      XA(2,6)=XO(8)
      XA(2,5)=XO(7)
      XA(2,4)=XO(6)
      XA(2,7)=XO(5)
      XA(2,14)=XO(1)
      XA(2,15)=XO(2)
      XA(2,16)=XO(3)
      XA(2,17)=XO(4)
      X8=1.0E0
      DO KK=1,KMAXi
         X8=X8-XO(KK)
      enddo
      XO2=AB(3)*XO(2)
C     ***************************************************************
C     *****
C     IONIZATION EQUIL. FOR CARBON INCLUDING AUGER IONIZATION
C     *****
C     ***************************************************************
c      ALC(1)=0.0E0
C      ALRAD FROM GOULD (78)
c$$$      ALTRAD=4.66E-13/T4**.5
c$$$      ALGS(12)=0.500*ALTRAD
c$$$      ALRAD=ALTRAD*(0.500+(1.-0.500)*FIR31)
c$$$      CALL DIEL(6.9d-4,1.1d5,3.d0,4.9d4,TE,ALD)
c$$$      CALL DIELB(1.08d-2,-.1075d0,.281d0,-1.93d-2,-.1127d0,T4,ALDB)
c$$$c!!      IF(DEN(IK).GT.1.E9) SUPR=0.17
c$$$c!!      IF(DEN(IK).LT.1.E9) SUPR=1.
c$$$       ALD=(ALD+ALDB)*SUPR
c$$$      IF(TE.LT.1000.) ALD=0.
c$$$      ALC2=ALRAD+ALD
c$$$      ALTOT(12)=ALC2
c$$$C     ALRAD FROM GOULD (1978)
c$$$      ALRAD=2.45E-12*(.474+.526/T4**.27)/SQRT(T4)
c$$$      CALL DIEL(7.0d-3,1.5d5,0.5d0,2.3d5,TE,ALD)
c$$$      CALL DIELB(1.8267d0,4.1012d0,4.8443d0,.2261d0,.596d0,T4,ALDB)
c$$$       ALD=(ALD+ALDB)*SUPR
c$$$C
c$$$C     C-T BHD 80
c$$$C
c$$$      CT=1.18E-12*XH/DELA(MI)
c$$$      ALC3=ALRAD+ALD+CT
c$$$      CALL DIEL(3.8d-3,9.1d4,2.d0,3.7d5,TE,ALD)
c$$$      CALL DIELB(2.3196d0,1.0733d1,6.883d0,-.1824d0,.4101d0,T4,ALDB)
c$$$       ALD=(ALD+ALDB)*SUPR
c$$$C
c$$$C     C-T BHD 80
c$$$C
c$$$      CT=(1.49+2.09*T4**.39)*1.E-9*XH/DELA(MI)
c$$$C
c$$$C     REC (GOULD 78)
c$$$C
c$$$      ALC4=5.05E-12/T4**.770+ALD+CT
c$$$      CALL DIEL(4.8d-2,3.4d6,0.2d0,5.1d5,TE,ALD)
c$$$       ALD=ALD*SUPR
c$$$C
c$$$C     C-T (B D 80)
c$$$C
c$$$      CT=(-.05+.81*T4**1.353)*1.E-9*XH/DELA(MI)
c$$$      ALC5=8.45E-12/T4**.817+ALD+CT
c$$$      CALL DIEL(4.8d-2,4.1d6,0.2d0,7.6d5,TE,ALD)
c$$$       ALD=ALD*SUPR
c$$$      ALC6=1.7E-11/T4**.721+ALD
c$$$      ALC7=AL5

C
C      COLL. IONIZATION (LOTZ)
C
c$$$      COC(1)=1.4E-10*SQRT(TE)*EXP(-13.1/T4)
c$$$      COC(2)=4.0E-10*SQRT(TE)*EXP(-35.8/T4)
c$$$      COC(3)=2.3E-11*SQRT(TE)*EXP(-55.6/T4)
c$$$      COC(4)=6.4E-12*SQRT(TE)*EXP(-74.8/T4)
c$$$      COC(5)=0.
c$$$      IF(T4.LT.10.) GOTO 891
c$$$      COC(5)=3.5E-13*SQRT(TE)*EXP(-454./T4)
c$$$  891 COC(6)=0.
c$$$      IF(T4.LT.10.) GOTO 892
c$$$      COC(6)=1.1E-13*SQRT(TE)*EXP(-568./T4)
c$$$  892 CONTINUE
c$$$C
c$$$C     AUGER IONIZATION EQUATIONS
c$$$C
c$$$      AUG=0.
c$$$      ZB(1)=ZE(12)+ZK(12)+CISEC(7)/DEN(IK)**2
c$$$      ZB(2)=ZE(13)+ZK(13)+(ZK(12)*AUG*ALC(2)*DELA(MI))/ZB(1)
c$$$      ZB(3)=ZE(8)+ZK(8)+(ZK(13)*AUG*ALC(3)*DELA(MI))/ZB(2)
c$$$      ZB(4)=ZE(9)+ZK(9)+(ZK(8)*AUG*ALC(4)*DELA(MI))/ZB(3)
c$$$c      ZB(2)=ZE(13)+ZK(13)
c$$$c      ZB(3)=ZE(8)+ZK(8)
c$$$c      ZB(4)=ZE(9)+ZK(9)
c$$$      ZB(5)=ZE(10)+ZK(10)
c$$$      ZB(6)=ZE(11)+ZK(11)
c$$$      DO  L=1,6
c$$$         XR(L)=(COC(L)+ZB(L)/DELA(MI))/ALC(L+1)
c$$$      enddo
c$$$c      CALL SOLVE(1,6,7,XR,XC)

      call augion(6,6,dela(mi),xc)            
      DO L=3,6
         XA(2,L+5)=XC(L)
      enddo
      XA(2,12)=XC(1)
      XA(2,13)=XC(2)
      X6=1.0E0
      DO 334 KK=1,6
  334 X6=X6-XC(KK)
C     ***************************************************************
C     *****
C     IONIZATION EQUIL. FOR NITOGEN INCLUDING AUGER IONIZATIONS.
C     *****
C     ***************************************************************
C
C     NITROGEN COLL. ION. RATES (LOTZ)
C
c$$$      CON(1)=6.52E-11*SQRT(TE)*EXP(-16.8/T4)
c$$$      CON(2)=4.33E-11*SQRT(TE)*EXP(-34.3/T4)
c$$$      CON(3)=1.25E-11*SQRT(TE)*EXP(-55./T4)
c$$$      CON(4)=0.
c$$$      IF(T4.LT.2.) GOTO 885
c$$$      CON(4)=8.96E-12*SQRT(TE)*EXP(-89.9/T4)
c$$$  885 CON(5)=0.
c$$$      IF(T4.LT.2.) GOTO 886
c$$$      CON(5)=2.83E-12*SQRT(TE)*EXP(-113.6/T4)
c$$$  886 CON(6)=0.
c$$$      CON(7)=0.
c$$$      IF(T4.LT.10.) GOTO 823
c$$$      CON(6)=1.44E-13*SQRT(TE)*EXP(-641./T4)/(1.+0.1*T4/641.)
c$$$      CON(7)=4.93E-14*SQRT(TE)*EXP(-774./T4)/(1.+0.1*T4/774.)
c$$$ 823  CONTINUE
c$$$C
c$$$C     NITROGEN RECOMB. COEFF (ALD.&PEQ)
c$$$C
c$$$      ALN1=0.
c$$$      CALL DIEL(5.2d-4,1.3d5,3.8d0,4.8d4,TE,ALD)
c$$$      CALL DIELB(0.d0,.631d0,.199d0,-1.97d-2,.4398d0,T4,ALDB)
c$$$c!!      IF(DEN(IK).GT.1.E9) SUPR=0.20
c$$$       ALD=(ALD+ALDB)*SUPR
c$$$C
c$$$C     C-T BD 79
c$$$C
c$$$      CT=1.0E-12*XH/DELA(MI)
c$$$      ALN2=4.1E-13/T4**.608+ALD+CT
c$$$      CALL DIEL(1.7d-3,1.4d5,4.1d0,6.8d4,TE,ALD)
c$$$      CALL DIELB(3.2d-2,-.6624d0,4.3191d0,3.d-4,.5946d0,T4,ALDB)
c$$$       ALD=(ALD+ALDB)*SUPR
c$$$C
c$$$C     C-T (BHD 1980 )
c$$$C
c$$$      CT=(.57+.29*T4**.46)*1.E-9*XH/DELA(MI)
c$$$      ALN3=2.2E-12/T4**.639+ALD+CT
c$$$      CALL DIEL(1.2d-2,1.8d5,1.4d0,3.8d5,TE,ALD)
c$$$      CALL DIELB(-.8806d0,1.124d1,3.071d1,-1.1721d0,.6127d0,T4,ALDB)
c$$$       ALD=(ALD+ALDB)*SUPR
c$$$C
c$$$C     C-T BHD 80
c$$$C
c$$$      CT=(-.82+3.75*T4**.67)*1.E-9*XH/DELA(MI)
c$$$      ALN4=5.0E-12/T4**.676+ALD+CT
c$$$      CALL DIEL(5.5d-3,1.1d5,3.0d0,5.9d5,TE,ALD)
c$$$      CALL DIELB(.4134d0,-4.6319d0,2.5917d1,-2.229d0,.236d0,T4,ALDB)
c$$$       ALD=(ALD+ALDB)*SUPR
c$$$C
c$$$C     CT (BUTLER AND DALG. 1980 )
c$$$C
c$$$      CT=(-.0036+.164*T4**1.455)*1.E-9*XH/DELA(MI)
c$$$      ALN5=9.55E-12/T4**.743+ALD+CT
c$$$      CALL DIEL(7.6d-2,4.7d6,0.2d0,7.2d5,TE,ALD)
c$$$       ALD=ALD*SUPR
c$$$      ALN6=1.5E-11/T4**.850+ALD
c$$$      CALL DIEL(6.6d-2,5.4d6,0.2d0,9.8d5,TE,ALD)
c$$$       ALD=ALD*SUPR
c$$$      ALN7=2.9E-11/T4**.750+ALD
c$$$      ALN8=AL(6)
c$$$C
c$$$C     AUGER IONIZATION EQUATIONS
c$$$C
c$$$      ZC(1)=ZE(18)+ZK(18)
c$$$      IF(ZC(1).LE.1.E-33) GOTO 1276
c$$$      ZC(2)=ZE(19)+ZK(19)+ALN(2)*DELA(MI)*ZK(18)/ZC(1)
c$$$      ZC(3)=ZE(20)+ZK(20)+ALN(3)*DELA(MI)*ZK(19)/ZC(2)
c$$$      ZC(4)=ZE(21)+ZK(21)+ALN(4)*DELA(MI)*ZK(20)/ZC(3)
c$$$      ZC(5)=ZE(22)+ZK(22)+ALN(5)*DELA(MI)*ZK(21)/ZC(4)
c$$$      ZC(6)=ZE(23)+ZK(23)
c$$$      ZC(7)=ZE(24)+ZK(24)
c$$$ 1276 CONTINUE
c$$$      IF(ILOWION.EQ.1) KMAXi=2
c$$$      IF(ILOWION.EQ.0) KMAXi=7
c$$$      DO  K=1,KMAXi
c$$$         XR(K)=(CON(K)+ZC(K)/DELA(MI))/ALN(K+1)
c$$$      enddo
c$$$c     CALL SOLVE(1,KMAX,8,XR,XN)

      call augion(7,7,dela(mi),xn)
      DO K=1,7
         XA(2,17+K)=XN(K)
      enddo
      X7=1.
      DO K=1,KMAXi
         X7=X7-XN(K)
      enddo
C     ***************************************************************
C     *****
C     IONIZATION EQUIL. FOR SILICON INCLUDING AUGER IONIZATION.
C     *****

C     ***************************************************************
C
C     SILICON RECOMB. COEFF. (ALD.&PEQ.)
C
C      ALRAD FROM GOULD (78)
      ALTRAD=6.48E-13/T4**.5
      ALGS(25)=0.677*ALTRAD
      ALRAD=ALTRAD*(0.677+(1.-0.677)*FIR31)
      CALL DIEL(6.2d-3,1.1d5,0.0d0,0.0d0,TE,ALD)
c!!       IF(DEN(IK).GT.1.E9) SUPR=0.25
c!!       IF(DEN(IK).LT.1.E9) SUPR=1.
       ALD=ALD*SUPR
      ALS2=ALRAD+ALD
      ALTOT(25)=ALS2
C     ALRAD FROM GOULD (1978)
      ALRAD=1.44E-12*(.205+.795/T4**.31)/SQRT(T4)
      CALL DIEL(1.4d-2,1.2d5,0.d0,0.d0,TE,ALD)
       ALD=ALD*SUPR
C
C     CT (MCCARROL AND VALANIAN 76 )
C
      CT=5.16E-9*XH/DELA(MI)
C
C     REC. FROM GOULD 78
C

      CALL DIEL(1.1d-2,1.0d5,0.d0,0.d0,TE,ALD)
C     DIELECTRIC REC. AT LOW T SCALED FROMCIV-CIII
      IF(T4.LT.2.) ALDB=1.45E-11*EXP(-0.41/T4)/T4**1.5
       ALD=(ALD+ALDB)*SUPR
C
C     CT (BUTLER AND DALGARNO 1980 )
C
      CT=(4.1E-10*XH+XHE*9.6E-10*T4**0.64)/DELA(MI)

      CALL DIEL(1.4d-2,1.2d6,0.d0,0.d0,TE,ALD)
       ALD=ALD*SUPR
C
C     C-T BD 80
C
      CT=(2.3E-9*XH+1.2E-9*XHE)/DELA(MI)

      DEEL=DELA(MI)
C
C     CALCULATE IONIZATION RATIOS INCLUDING AUGER IONIZATION.
C
      IF(ILOWION.EQ.1) KMAXi=6
      IF(ILOWION.EQ.0) KMAXi=14
c      CALL AUGSI(KMAX,DEEL,XR)
c     CALL SOLVE(1,KMAX,15,XR,XSi)

c     *************
c     Silicon eqiul. incl. Auger
c     *************
      call augion(14,14,deel,xsi)
      DO  K=1,14
         XA(2,24+K)=XSi(K)
      enddo
C     ***************************************************************
C     *****
C     IONIZATION EQUIL. FOR MAGNESIUM INCLUDING AUGER IONIZATION.
C     *****
C     ***************************************************************

      call augion(12,12,dela(mi),xm)
      DO K=1,3
         XA(2,38+K)=XM(K)
      enddo
      XMGI=XM(1)
C     ***************************************************************
C     *****
C     IONIZATION EQUIL. FOR CALCIUM.
C     *****
C     ***************************************************************

      call augion(20,20,dela(mi),xca)

      DO  K=1,6
         if(k.le.3) then
            XA(2,49+K)=XCA(K)
         elseif(k.ge.4) then
            XA(2,108+K)=XCA(K)
         endif
      enddo
      XCAII=XCA(2)
      
c     *************
c     Neon ion eqiul. incl. Auger
c     *************
      call augion(10,10,dela(mi),xne)
      DO K=1,10
         XA(2,71+K)=XNE(K)
      enddo

C     ***************************************************************
C     *****
C     IONIZATION EQUIL. FOR ARGON
C     *****
C     ***************************************************************

c     *************
c     Argon ion eqiul. incl. Auger
c     *************
      call augion(18,7,deel,xar)
      DO K=1,7
         if(k.le.2) then
            XA(2,81+K)=XAR(K)
         else
            XA(2,93+K)=XAR(K)
         endif
      enddo
C     ***************************************************************
C     *****
C     IONIZATION EQUIL. FOR SODIUM.
C     *****
C     ***************************************************************

      call augion(11,11,dela(mi),xna)
      DO K=1,2
         XA(2,52+K)=XNA(K)
      enddo
CK     ***************************************************************   
CK     *****                                                             
CK     IONIZATION EQUIL. FOR SULPHUR INCLUDING AUGER IONIZATIONS.       
CK     *****                                                             
CK     ***************************************************************   

      call augion(16,16,dela(mi),xsu)
      DO  K=1,16
         XA(2,55+K)=XSU(K)
      enddo

C     ***************************************************************
C     *****
C     IONIZATION EQUIL. FOR ALUMINUM INCLUDING AUGER IONIZATION.
C     *****
C     ***************************************************************

      call augion(13,13,dela(mi),xal)
      DO K=1,4
         XA(2,45+K)=XAL(K)
      enddo
CK     ***************************************************************   
CK     *****                                                             
CK     IONIZATION EQUIL. FOR IRON (NO AUGER IONIZATIONS).
CK     *****                                                             
CK     ***************************************************************   
 7361 CONTINUE

      call augion(26,12,dela(mi),xf)
      DO  K=1,4
         XA(2,K+41)=XF(K)
      enddo
      DO K=1,11
         XA(2,K+83)=XF(K+4)
      enddo
      DO K=16,26
         XA(2,K+85)=XF(K)
      enddo

      call conv_xa_xion
C     ***************************************************************
C     *****
C     CALCULATE THE FRACTION OF FREE ELECTRONS.
C     ITERATE UNTIL THE RELATIVE CHANGE IN X(EL) IS LESS THAN 0.01
C     *****
C     ***************************************************************
      
      DC=AB(4)*(6.*X6+5.*XC(6)+4.*XC(5)+3.*XC(4)+2.*XC(3)+1.*XC(2))
      
      DNI=AB(5)*(7.*X7+6.*XN(7)+5.*XN(6)+4.*XN(5)+3.*XN(4)
     &                                          +2.*XN(3)+XN(2))
      DOX=AB(3)*(8.*X8+7.*
     &      XO(8)+6.*XO(7)+5.*XO(6)+4.*XO(5)+3.*XO(4)+
     &      2.*XO(3)+XO(2))

      DMG=AB(7)*(XM(2)+2.*XM(3)+3.*(1.-XM(1)-XM(2)-XM(3)))

c Ne+10 = Ne XI
      X11=1.
      DO IO=1,10
         X11=X11-XNe(IO)
      enddo
      DNe=0.
      DO IO=2,10
         DNE=DNE+DBLE(IO-1)*XNe(IO)
      enddo
      DNE=DNE+X11*10.
      DNE=DNE*abn(6)

c Si+14 = Si XV
      X15=1.
      DO IO=1,14
         X15=X15-XSi(IO)
      enddo
      DSI=0.
      DO IO=2,14
         DSI=DSI+DBLE(IO-1)*XSi(IO)
      enddo
      DSI=DSI+X15*14.
      DSI=DSI*ABN(10)      

c S+16 = Si XVII
      X17=1.
      DO IO=1,16
         X17=X17-XSi(IO)
      enddo
      DSU=0.
      DO IO=2,16
         DSU=DSU+DBLE(IO-1)*XSU(IO)
      enddo
      DSU=DSU+X17*16.
      DSU=DSU*ABn(11)      
      X8=1.
      DO IO=1,7
         X17=X8-XAr(IO)
      enddo
      DAR=0.
      DO IO=2,7
         DAR=DAR+DBLE(IO-1)*XAR(IO)
      enddo
      DAR=DAR+X8*7.
      DAR=DAR*ABn(12)            
      
      DCA=AB(11)*(XCA(2)+2.*(1.-XCA(1)-XCA(2)))
      
      DFE=0.
      DO IO=2,26
         DFE=DFE+DBLE(IO-1)*XF(IO)
      ENDDO
      DFE=AB(8)*DFE
      
      DNEW=DC+DNI+dox+dne+DMG+DSI+DSU+DAR+DCA+DFE
       
      DELTEL=ABS(DEL(IK-1)-DNEW)/DNEW
      IF(DELTEL.LT.0.003) THEN
        DELA(MI)=DNEW
        XEL=DNEW
        GOTO 7782
      ENDIF
      FD(MI)=DNEW-DELA(MI)
C     IF(IK.EQ.1) GOTO 7782
      IF(FD(MI).EQ.0.) then
         DELA(MI+1)=DNEW
         GOTO 7328
      endif
      IF(MI.GE.2) then
         GOTO 7778
      endif
C     IF(IK.GT.3) GOTO 7790
      DELA(2)=0.95*DELA(1)
      GOTO 7779
C7790 DELA(2)=10.**(2.*LOG10(DEL(IK-1))-LOG10(DEL(IK-2)))
 7779 MI=2
      GOTO 776
 7778 DELA(MI+1)=(DELA(MI-1)-FD(MI-1)*DELA(MI)/FD(MI))/
     &(1.-FD(MI-1)/FD(MI))
7328  XEL=DELA(MI+1)
9090  FORMAT(1X,'XEL',2I5,8E11.4)
      IF(MI.LT.MIMAX) then
         RATIO=ABS((DELA(MI+1)-DELA(MI))/DELA(MI))
         IF(RATIO.GT.0.5.OR.DELA(MI+1).LT.0.) DELA(MI+1)=DNEW
         IF(RATIO.LT.0.003) GOTO 7782
         MI=MI+1
         GOTO 776
      endif
      WRITE(6,7781)RATIO,(DELA(MIK),MIK=1,49)
 7781 FORMAT(1X,'NO CONVERGENCE',50E11.3)
 7782 DEL(IK)=DELA(MI)
      ZEL=DEL(IK)-AB(11)*XCA(2)
C     ***************************************************************
C     *****
C     RECOMBINATION COOLING (CASE A) FOR HYDROGEN AND HELIUM
C     RATES ARE ADJUSTED TO AGREE WITH OSTERBROCKS VALUES FOR
C     T<2.E4*Z**2
C     *****
C     ***************************************************************
c$$$      YI(1)=1.57E5/TE
c$$$      YI(2)=YI(1)
      YI(3)=4.*YI(1)
      YI(4)=64.*YI(1)
      YI(5)=36.*YI(1)
      YI(6)=49.*YI(1)
      YI(7)=196.*YI(1)

C
C     CARBON
C
      RECCA=0.
      DO K=1,6
      RECCA=RECCA+AB(4)*XC(K+1)*1.5*1.38E-16*TE*ALC(K+1)
      ENDDO
C
C     OXYGEN
C
      RECOX=0.
      DO K=1,8
         ALFO=ALO(K+1)
         IF(K.EQ.1) ALFO=ALO(2)-ALCTO2
         IF(K.EQ.2) ALFO=ALO(3)-ALCTO3
         RECOX=RECOX+AB(3)*XO(K+1)*1.5*1.38E-16*TE*ALFO
      ENDDO
C
C     SILICON
C
      RECSI=0.
      DO K=1,14
         CSI=ABN(10)*XSi(K+1)*1.5*1.38E-16*TE*ALSi(K+1)
         RECSI=RECSI+CSI
      enddo
C     SULPHUR
C
C
      RECS=0.
      DO  K=1,10
         CSU=AB(12)*XSU(K+1)*1.5*1.38E-16*TE*ALSU(K+1)
         RECS=RECS+CSU
      enddo
C     IRON
C
C
      RECFE=0.
      DO K=1,15
         CFE=AB(8)*XF(K+1)*1.5*1.38E-16*TE*ALFE(K+1)
         RECFE=RECFE+CFE
      enddo
C     ***************************************************************
C     *****
C     FREE-FREE COOLING,(COX&TUCKER AP.J. 157:1157) FOR T>1+7
C     AND SPITZER TABLE 3.3 FOR GAUNT FACTOR
C     *****
C     ***************************************************************
C      IF(TE.GT.1.E7) GOTO 8436
      FF=0.
      DO 8437 KK=4,7
      ZZ=1.
      IF(KK.GE.4) ZZ=(KK-3.)**2
      GAUNT=-1.08+0.925*LOG10(TE/ZZ)-0.085*(LOG10(TE/ZZ))**2.
      IF(KK.EQ.7) ABX=AB(8)*(1.-XF(1)-XF(2)-XF(3)-XF(4))
      IF(KK.GE.4.AND.KK.LE.6) ABX=AB(8)*XF(KK-2)
8437  FF=FF+1.426E-27*ZZ*ABX*GAUNT*SQRT(TE)
      GOTO 8438
8436  FF=1.426E-27*(1.23*AB(1)+4.*1.3587*AB(2))*SQRT(TE)
8438  CONTINUE


C     **************************************************************
C     *****
C     FORBIDDEN LINE COOLING
C     *****
C     **************************************************************
C     (C I) 4619,8729,9812 A. NUSSBAUMER & RUSCA AA 72, 129.
      O21=0.603*(TE/5000.)**0.96
      O31=0.149*(TE/5000.)**0.871
      O32=0.196*(TE/5000.)**0.499
      Z=AB(4)*XC(1)
       CALL FORB3(1,12,1.262D0,1.420D0,3.26D-4,2.73D-3,0.528D0,O21
     &,O31,O32,9.D0,5.D0,1.D0,FB(20),FB(21),FB(22),TE,Z,F)
      CALL DIELB(-2.02d-2,.3799d0,0.0890d0,-0.57d-2,.9237d0,T4,ALDB)
      DIFB=AB(4)*XC(2)*ALDB*1.2633*1.602E-12
      FB(20)=FB(20)+AB(4)*XC(2)*ALDB*1.2633*1.602E-12
      C211=AB(4)*XC(1)*F

      do k=1,3
         cinx(3,1,k)=fb(k)
         if(k==1) then
            wlix(3,1,k)=wl21
         elseif(k==2) then
            wlix(3,1,k)=wl31
         elseif(k==3) then
            wlix(3,1,k)=wl32
         endif
         ilabcfx(3,1,k)=301
      enddo
      kmaxp(3,1)=3
      
C     
C     COLL. EXIT. OF CARBON
C
C     C I 2966-68 A  OMEGA FROM HAYES & NUSSBAUMER AA 134:193
C
      Z=AB(4)*XA(2,12)
      OM=0.475*(TE/5000.)**0.5
      CALL RLOSS(12,2967.D0,OM,9.D0,5.D0,6.D8,RS,XEL,Z,TE,C212
     &     ,W(27),CIN(27),FR(2,27),27,WLI(27))
      cinx(3,1,4)=cin(27)
      wlix(3,1,4)=wli(27)
      ilabcfx(3,1,4)=301
      kmaxp(3,1)=kmaxp(3,1)+1
C     C I 370.4 609.1 MY MAZDA COLLISION STRENGTHS FROM O I!!
      O32=9.76E-6*(TE-228)+3.46E-11*(TE-228.)**2
      O31=3.39E-6*(TE-326.)-2.90E-11*(TE-326.)**2
      O21=1.89E-6*(TE-326.)+8.00E-11*(TE-326.)**2
      IF(O21.LT.0.) O21=1.e-10
      IF(O31.LT.0.) O31=1.e-10
      IF(O32.LT.0.) O32=1.e-10
      Z=AB(4)*XC(1)
      CALL FORB3(13,12,3.348D-3,2.03D-3,7.93D-8,1.71D-14,2.65D-7,O21,
     &O31,O32,1.D0,3.D0,5.D0,FB(32),FBQQ,FB(33),TE,Z,F)
      DO NM=32,33
         FB(NM)=AB(4)*XC(1)*FB(NM)
      enddo
      do k=1,3
         if(k==1) then
            wlix(3,1,k+4)=wl21
            cinx(3,1,k+4)=fb(32)
         elseif(k==2) then
            wlix(3,1,k+4)=wl31
            cinx(3,1,k+4)=AB(4)*Xc(1)*fbqq
         elseif(k==3) then
            wlix(3,1,k+4)=wl32
            cinx(3,1,k+4)=fb(33)
         endif
         ilabcfx(3,1,k+4)=301
      enddo
      kmaxp(3,1)=kmaxp(3,1)+3
      C213=AB(4)*XC(1)*F
      line_cool(3,1)=c211+c212+c213
      
C     C II 157.73 MY HAYES & NUSSBAUMER  A A 134, 193.
      Z=AB(4)*XC(2)
      IF(TE.LT.1.E3) OM=1.82
      IF(TE.GT.1.E3) OM=1.87*T3**0.21
      CALL RLOSS(12,1.58D6,OM,2.D0,4.D0,2.29D-6,RS,XEL,Z,TE,C223
     &     ,W(26),CIN(26),FRQQ,26,WLI(26))
      cinx(3,2,1)=cin(26)
      wlix(3,2,1)=wli(26)
      ilabcfx(3,2,1)=302
C
C     C II 1334 A  OMEGA FROM HAYES & NUSSBAUMER AA 134:193
C
      Z=AB(4)*XA(2,13)
      OM=5.79*T4**0.05
      CALL RLOSS(13,1334.D0,OM,6.D0,10.D0,6.D8,RS,XEL,Z,TE,C221
     &     ,W(1),CIN(1),FR(2,1),1,WLI(1))
      wlix(3,2,2)=wli(1)
      cinx(3,2,2)=cin(1)
      ilabcfx(3,2,2)=302
C
C     C II] 2326 OMEGA FROM HAYES & NUSSBAUMER AA 134:193 A=44.1
C
      CALL RLOSS(13,2326.D0,2.8D0,6.D0,12.D0,4.41D1,RS,XEL,Z,TE,C222
     &     ,W(2),CIN(2),FR(2,2),2,WLI(2))
      wlix(3,2,3)=wli(2)
      cinx(3,2,3)=cin(2)
      ilabcfx(3,2,3)=302
      kmaxp(3,2)=3
      line_cool(3,2)=c221+c222+c223
C
C     C III 977 A (FLOWER & LAUNAY -73)
C
      Z=AB(4)*XA(2,8)
      OM=4.34*T4**0.112
      CALL RLOSS(8,977.D0,OM,1.D0,3.D0,1.79D9,RS,XEL,Z,TE,C232
     &     ,W(3),CIN(3),FR(2,3),3,WLI(3))
      wlix(3,3,1)=wli(3)
      cinx(3,3,1)=cin(3)
      ilabcfx(3,3,1)=303
C
C     CIII) 1909 (MENDOZA)
C
      Z=AB(4)*XA(2,8)
      CALL RLOSS(8,1909.D0,1.00D0,1.D0,9.D0,3.2D1,RS,XEL,Z,TE,C231
     &     ,W(4),CIN(4),FR(2,4),4,WLI(4))
      wlix(3,3,2)=wli(4)
      cinx(3,3,2)=cin(4)
      ilabcfx(3,3,2)=303
      kmaxp(3,3)=2
      line_cool(3,3)=c231+c232
C
C     DIEL. REC. TO C III 1909 (STOREY 81),WITH SURPRESSION 0.17
C
c!!      IF(DEN(IK).GT.1.E9) SUPR=0.17
c!!      IF(DEN(IK).LT.1.E9) SUPR=1.
C     IF(XA(2,9).GT.1.E-20) REC31=AB(4)*XA(2,9)*SUPR*7.73E-23/T4**0.89
C
C     C IV 1546 A,(OSTERBROCK 77) OMEGA=8.80
C
      Z=AB(4)*XA(2,9)
      CALL RLOSS(9,1550.D0,8.66D0,2.D0,6.D0,2.65D8,RS,XEL,Z,TE,C241
     &     ,W(5),CIN(5),FR(2,5),5,WLI(5))
      cinx(3,4,1)=cin(5)
      taulinex(3,4,1)=t0
      wlix(3,4,1)=wli(5)
      ilabcfx(3,4,1)=304
      kmaxp(3,4)=1
      line_cool(3,4)=c241
      
C
C     NITROGEN
C
C     (N I ) 3468, 5201, 10406 (MENDOZA 83)
      O21=(0.48*T4)
      O31=(0.17*T4)
      O32=(0.62*T4)
      Z=AB(5)*XN(1)
      CALL FORB3(7,18,2.401D0,1.200D0,1.24D-5,5.28D-3,8.85D-2,O21,
     &O31,O32,4.D0,10.D0,6.D0,FB(13),FB(14),FB(15),TE,Z,F)
      C311=AB(5)*XN(1)*F
      DO NM=13,15
         FB(NM)=AB(5)*XN(1)*FB(NM)
      enddo
      do k=1,3         
         cinx(4,1,k)=fb(k+12)
         if(k==1) then
            wlix(4,1,k)=wl21
         elseif(k==2) then
            wlix(4,1,k)=wl31
         elseif(k==3) then
            wlix(4,1,k)=wl32
         endif
         ilabcfx(4,1,k)=401
      enddo
      kmaxp(4,1)=3
      line_cool(4,1)=c331
C     (N II)  3063,5755,6548-83 OSTERBROCK
      Z=AB(5)*XN(2)
      CALL FORB3(8,19,1.89D0,2.15D0,3.0D-3,3.4D-2,1.1D0,2.99D0
     &,0.36D0,0.39D0,9.D0,5.D0,1.D0,FB(7),FB(8),FB(9),TE,Z,F)
      C322=AB(5)*XN(2)*F
      DO NM=7,9
         FB(NM)=AB(5)*XN(2)*FB(NM)
      enddo
      do k=1,3         
         cinx(4,2,k)=fb(k+6)
         if(k==1) then
            wlix(4,2,k)=wl21
         elseif(k==2) then
            wlix(4,2,k)=wl31
         elseif(k==3) then
            wlix(4,2,k)=wl32
         endif
      enddo
      ilabcfx(4,2,k)=402
C
C     N II 2143 .OMEGA=1.29 (JACKSON 73)  ONLY GUESS OF A21 FROM OIII 16
C
      Z=AB(5)*XN(2)
      CALL RLOSS(19,2143.D0,1.29D0,9.D0,5.D0,2.D2,RS,XEL,Z,TE,C321
     &     ,W(10),CIN(10),FR(2,10),10,WLI(10))
      wlix(4,2,4)=wli(10)
      cinx(4,2,4)=cin(10)
      ilabcfx(4,2,4)=402
      kmaxp(4,2)=4
      line_cool(4,2)=c321+c322
C
C     N III 1750.OMEGA=1.89 (JACKSON 73)
C
      Z=AB(5)*XN(3)
      CALL RLOSS(20,1750.D0,1.9D0,6.D0,12.D0,2.2D2,RS,XEL,Z,TE,C331
     &     ,W(11),CIN(11),FR(2,11),11,WLI(11))
      wlix(4,3,1)=wli(11)
      cinx(4,3,1)=cin(11)
      ilabcfx(4,2,1)=403
C     N III 57.33 MAZDA
      Z=AB(5)*XN(3)
      CALL RLOSS(20,5.733D5,0.701d0,2.D0,4.D0,4.77D-5,RS,XEL,Z,TE,C333
     &     ,W(82),CIN(82),FRQQ,82,WLI(82))
      wlix(4,3,2)=wli(82)
      cinx(4,3,2)=cin(82)
      ilabcfx(4,3,2)=403
C
C     DIEL. REC. TO N IV 1486 (STOREY 81),WITH SURPRESSION 0.20
C                  
C
C     NIII 990 A (MENDOZA)C ALPINE) OMEGA=5.22 G1=6
C
      Z=AB(5)*XN(3)
      CALL RLOSS(20,990.D0,5.22D0,6.D0,10.D0,7.3D8,RS,XEL,Z,TE,C332
     &     ,W(12),CIN(12),FR(2,12),12,WLI(12))
      wlix(4,3,3)=wli(12)
      cinx(4,3,3)=cin(12)
      ilabcfx(4,3,3)=403
      kmaxp(4,3)=3
      line_cool(4,3)=c331+c332+c333            
            
C
C     N IV 1486.OMEGA=0.82 (OSTERBROCK 70)
C
      Z=AB(5)*XN(4)
      CALL RLOSS(21,1486.D0,0.82D0,1.D0,9.D0,1.88D2,RS,XEL,Z,TE,C341
     &     ,W(13),CIN(13),FR(2,13),13,WLI(13))
      wlix(4,4,1)=wli(13)
      cinx(4,4,1)=cin(13)
      ilabcfx(4,4,1)=404
      kmaxp(4,4)=1
C
C     N IV 765.15 A (MENDOZA)
C
      C342=AB(5)*XN(4)*0.776E-15*EXP(-1.88E5/TE)/SQRT(TE)
      wlix(4,4,2)=765.15
      cinx(4,4,2)=c342
      ilabcfx(4,4,2)=404
      kmaxp(4,4)=2
      line_cool(4,4)=c341+c342
C
C     N V 1238 A .OMEGA=6.61 (OSTERBROCK&WALLACE)
C
      Z=AB(5)*XN(5)
      CALL RLOSS(22,1238.D0,6.65D0,2.D0,6.D0,3.38D8,RS,XEL,Z,TE,C351
     &     ,W(14),CIN(14),FR(2,14),14,WLI(14))
      cinx(4,5,1)=cin(14)
      taulinex(4,5,1)=t0
      wlix(4,5,1)=wli(14)
      ilabcfx(4,5,1)=405
      kmaxp(4,5)=1
      line_cool(4,5)=c351


C     COLL. EXIT. OF OXYGEN
C
C     (O I ) 2964, 5581, 6302-64 (MENDOZA 83)
      IN=N
      INP=NP
      TEV=TE/1.1609E4
      IF(T4.GT.0.5) GOTO 10
      O21=.0151*T3**1.31
      O31=.00184*T3**1.32
      O32=.031*T3**.534
      GOTO 50
 10   IF(T4.GT.1.) GOTO 20
      O21=.266*T4**1.10
      O31=.034*T4**1.08
      O32=.105*T4**.52
      GOTO 50
 20   O21=.266*T4**.91
      O31=.0324*T4**.91
      O32=.105*T4**.50
 50   CONTINUE
      Z=AB(3)*XO(1)
      CALL FORB3(4,14,1.957D0,2.223D0,8.45D-3,7.35D-2,1.22D0,O21,
     &O31,O32,9.D0,5.D0,1.D0,FB(10),FB(11),FB(12),TE,Z,F)
      IF(AB(3).LT.1.E-5.OR.TE.LT.TLIMO) then
         C111=AB(3)*XO(1)*F
      endif
      DO NM=10,12
         FB(NM)=AB(3)*XO(1)*FB(NM)
      enddo
      do k=1,3         
         cinx(5,1,k)=fb(k+9)
         if(k==1) then
            wlix(5,1,k)=wl21
         elseif(k==2) then
            wlix(5,1,k)=wl31
         elseif(k==3) then
            wlix(5,1,k)=wl32
         endif
         ilabcfx(5,1,k)=501
      enddo
C     O I 63.18, 145.5 MY MAZDA
      O21=(.0018*(TE/1000.)**1.169)
      IF(TE.LE.1.E3) O31=(.0022*(TE/1000.)**1.874)
      IF(TE.GT.1.E3) O31=(.0292*(TE/10000.)**1.123)
      IF(TE.GT.1.E3) O32=(.0987*(TE/10000.)**1.113)
      IF(TE.LE.1.E3) O32=(.0076*(TE/1000.)**1.493)
      Z=AB(3)*XO(1)
      CALL FORB3(14,14,1.963D-2,8.521D-3,8.92D-5,0.D0,1.74D-5,O21,O31,
     &O32,5.D0,3.D0,1.D0,FB(17),FB(18),FB(19),TE,Z,F)
      DO NM=17,19
         FB(NM)=AB(3)*XO(1)*FB(NM)
      enddo
      do k=1,3         
         cinx(5,1,k+3)=fb(k+16)
         if(k==1) then
            wlix(5,1,k+3)=wl21
         elseif(k==2) then
            wlix(5,1,k+3)=wl31
         elseif(k==3) then
            wlix(5,1,k+3)=wl32
         endif
         ilabcfx(5,1,k+3)=501
      enddo      
      C113=AB(3)*XO(1)*F
      kmaxp(5,1)=6
      line_cool(5,1)=c111+c113
c     O I 7949. DIEL. REC. 
      CALL DIELB(0.d0,0.0400d0,0.d0,0.d0,0.5587d0,T4,ALDB)
      CIN(87)=AB(3)*XA(2,15)*ALDB*1.9864E-8/7949.
      WLI(87)=7949.
c     O I 6319. REC. 
      CALL DIELB(0.d0,0.0280d0,0.d0,0.d0,1.0257d0,T4,ALDB)
      CIN(88)=AB(3)*XA(2,15)*ALDB*1.9864E-8/6319.
      WLI(88)=6319.
C     (O II) 2471,3727,7327 (MENDOZA 83)
      Z=AB(3)*XO(2)
      CALL FORB3(3,15,3.32D0,1.69D0,8.89D-5,4.53D-2,1.73D-1,1.335D0,
     &0.405D0,1.30D0,4.D0,10.D0,6.D0,FB(4),FB(5),FB(6),TE,Z,F)
      C121=AB(3)*XO(2)*F
      DO NM=4,6
         FB(NM)=AB(3)*XO(2)*FB(NM)
      enddo
      do k=1,3         
         cinx(5,2,k)=fb(k+3)
         if(k==1) then
            wlix(5,2,k)=wl21
         elseif(k==2) then
            wlix(5,2,k)=wl31
         elseif(k==3) then
            wlix(5,2,k)=wl32
         endif
         ilabcfx(5,2,k)=502
      enddo      
c     O II 4651 REC. (BORKOWSKI AND SHULL 1990 AP J 348,169)
      CIN(84)=AB(3)*XA(2,16)*ALEX(15)*1.9864E-8/4651.*0.16
      WLI(84)=4651.
      wlix(5,2,4)=wli(84)
      cinx(5,2,4)=cin(84)
      ilabcfx(5,2,4)=502
      kmaxp(5,2)=4
      line_cool(5,2)=c121
C     (O III) 2321,4363,4959 & 5007 MENDOZA
      Z=AB(3)*XO(3)
c      CALL FORB3(2,16,2.49D0,2.84D0,2.63D-2,2.23D-1,1.78D0,2.30D0,
c     &0.299D0,0.64D0,9.D0,5.D0,1.D0,FB(1),FB(2),FB(3),TE,Z,F)
C
C     O III 3P,1D,1S FROM MAZDA OMEGAS AT 15000 K
C

      CALL FIVELEV_dp(16,5,3,5,0.d0,113.2d0,306.2d0,20273.3d0,43185.7d0,
     &     1.d0,3.d0,5.d0,5.d0,1.d0,2.62d-5,3.02d-11,2.74d-6,0.d0,
     &     9.76d-5,6.74d-3,0.223d0,1.96d-2,7.85d-4,1.78d0,0.553d0,
     &     0.281d0,0.256d0,0.0332d0,1.32d0,0.767d0,0.0997d0,1.278d0,
     &     0.166d0,0.638d0,TE,Z,XI)
c     POPULATION OF EXCITED STATES
      XO31D=XI(4)
C     TOTAL O III 1D POPULATION
      XA(2,95)=XI(4)*XA(2,16)
      DO I=1,5
         DO J=1,5
            EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      COOIII=0.
      DO  I=1,5
        DO  J=I+1,5
          COOIII=COOIII+EMCA(J,I)
        ENDDO
      ENDDO
      C131=COOIII
      FB(1)=EMCA(4,2)+EMCA(4,3)
      FB(2)=EMCA(5,2)
      FB(3)=EMCA(5,4)
      FB(49)=EMCA(2,1)
      FB(50)=EMCA(3,2)

      do k=1,kmaxp(5,3)
         wlix(5,3,k)=wlin(k)
         cinx(5,3,k)=weh(k)/(del(ik))
         ilabcfx(5,3,k)=503
      enddo
      
C
C     O III 1664 A=204 (DOSHEK ET AL 78 AP.J.) OM=1.36 JACKSO 73
C
      Z=AB(3)*XO(3)
      CALL RLOSS(16,1664.D0,1.36D0,9.D0,1.D0,2.04D2,RS,XEL,Z,TE,C132
     &     ,W(9),CIN(9),FR(2,9),9,WLI(9))
      wlix(5,3,kmaxp(5,3)+1)=wli(9)
      cinx(5,3,kmaxp(5,3)+1)=cin(9)
      ilabcfx(5,3,kmaxp(5,3)+1)=503
c     O III 3762. REC. (BORKOWSKI AND SHULL 1990 AP J 348,169)
      CALL DIELB(5.30d-3,-9.830d-2,.4693d0,5.40d-3,0.5684d0,T4,ALDB)
      CIN(85)=AB(3)*XA(2,17)*(ALEX(16)*0.17+ALDB)*1.9864E-8/3762.
      WLI(85)=3762.
      wlix(5,3,kmaxp(5,3)+2)=wli(85)
      cinx(5,3,kmaxp(5,3)+2)=cin(85)
      ilabcfx(5,3,kmaxp(5,3)+2)=503
c     O III 3266. DIEL. REC. 
      CALL DIELB(1.d-4,.9493d0,.1777d0,.0623d0,1.2049d0,T4,ALDB)
      CIN(89)=AB(3)*XA(2,17)*ALDB*1.9864E-8/3266.
      WLI(89)=3266.
      wlix(5,3,kmaxp(5,3)+3)=wli(89)
      cinx(5,3,kmaxp(5,3)+3)=cin(89)
      ilabcfx(5,3,kmaxp(5,3)+3)=503
      kmaxp(5,3)=kmaxp(5,3)+3
      line_cool(5,3)=cooiii+c132
      
C     O IV 788 A (MC ALPINE) OMEGA=3.2 G1=6
C     C141=1.18E-16*AB(3)*XO(4)*EXP(-1.83E5/TE)/SQRT(TE)
C
C     O IV 1400 OMEGA=1.37 (O&W -77)
C
      Z=AB(3)*XO(4)
      CALL RLOSS(17,1400.D0,1.37D0,6.D0,12.D0,1.6D3,RS,XEL,Z,TE,C142
     &     ,W(8),CIN(8),FR(2,8),8,WLI(8))
      cinx(5,4,1)=cin(8)
      taulinex(5,4,1)=t0
      wlix(5,4,1)=wli(8)
      ilabcfx(5,4,1)=504
C     O IV 25.88 MU MAZDA
      Z=AB(3)*XO(4)
      CALL RLOSS(17,2.588D5,2.33D0,2.D0,4.D0,5.20D-4,RS,XEL,Z,TE,C143
     &     ,W(83),CIN(83),FRQQ,83,WLI(83))
      cinx(5,4,2)=cin(83)
      taulinex(5,4,2)=t0
      wlix(5,4,2)=wli(83)
      ilabcfx(5,4,2)=504
c     O IV 3066. REC. (BORKOWSKI AND SHULL 1990 AP J 348,169)
      CIN(86)=AB(3)*XA(2,7)*ALEX(17)*1.9864E-8/3066.*0.16
      WLI(86)=3066.
      cinx(5,4,3)=cin(86)
      taulinex(5,4,3)=t0
      wlix(5,4,3)=wli(86)
      ilabcfx(5,4,3)=504
c     O IV 3027. DIEL. REC. 
      CALL DIELB(0.0d0,1.5217d0,0.d0,0.d0,0.7085d0,T4,ALDB)
      CIN(90)=AB(3)*XA(2,7)*ALDB*1.9864E-8/3027.
      WLI(90)=3027.
      cinx(5,4,4)=cin(90)
      taulinex(5,4,4)=t0
      wlix(5,4,4)=wli(90)
      ilabcfx(5,4,4)=504
      kmaxp(5,4)=4
      line_cool(5,4)=c142+c143
C
C     EXIT. OF OV (MENDOZA)
C
      Z=AB(3)*XO(5)
      CALL RLOSS(7,1216.D0,0.721D0,1.D0,9.D0,2.25d3,RS,XEL,Z,TE,C151
     &     ,W(7),CIN(7),FR(2,7),7,WLI(7))
      cinx(5,5,1)=cin(7)
      wlix(5,5,1)=wli(7)
      ilabcfx(5,5,1)=505
C     O V 629.7 A Pradhan compil 95
      om=2.76*t4**0.046
      CALL RLOSS(7,629.7D0,om,1.D0,3.D0,2.80D9,RS,XEL,Z,TE,C152
     &     ,WX,CINq,FRX3,7,WLIq)
      
      cinx(5,5,2)=cinq
      wlix(5,5,2)=wliq
      ilabcfx(5,5,2)=505
c     O V 6488. DIEL. REC.
      if(te.gt.1.e3.and.te.lt.6.e4) then
         CALL DIELB(-.0765d0,4.168d0,1.2966d0,-.1261d0,2.8234d0,T4,ALDB)
      else
         aldb=0.
      endif
      CIN(91)=AB(3)*XO(6)*ALDB*1.9864E-8/6488.
      WLI(91)=6488.
      cinx(5,5,3)=cin(91)
      wlix(5,5,3)=wli(91)
      ilabcfx(5,5,3)=505
      kmaxp(5,5)=3
      line_cool(5,5)=c151+c152
                  
C     EXIT. OF O VI, A & Omega from Pradhan
C
      om=5.00*t4**0.014
      CALL RLOSS(4,1034.D0,om,2.D0,6.D0,4.12D8,RS,XEL,Z,TE,c161
     &     ,W(6),CIN(6),FR(2,6),6,WLI(6))
      cinx(5,6,1)=cin(6)
      wlix(5,6,1)=wli(6)
      ilabcfx(5,6,1)=506
      kmaxp(5,6)=1
      line_cool(5,6)=c161
C
C     NEON
C
C     NE II 12.814 MY Pradhan
      Z=ABN(6)*XNE(2)
      OM=0.303*t4**0.065
      CALL RLOSS(6,12.814D4,OM,4.D0,2.D0,8.55D-3,RS,XEL,Z,TE,c1322
     &     ,W(46),CIN(46),FRQQ,46,WLI(46))
      cinx(6,2,1)=cin(46)
      taulinex(6,2,1)=t0
      wlix(6,2,1)=wli(46)
      ilabcfx(6,2,1)=602
      kmaxp(6,2)=1
      line_cool(6,2)=c1322

c     New Ne III

c     Coll. strengths from Butler & Zeippen 1994, A:s from Pradhan comp.

      call FIVELEV_dp(74,6,3,5,
     &     0.0d0,642.9d0,920.4d0,25840.8d0,55750.6d0,
     &     5.d0,3.d0,1.d0,1.d0,1.d0,
     &     5.97d-3,2.18d-8,1.71d-1,3.94d-3,1.15d-3,5.42d-2,
     &     2.00d0,8.51d-6,0.0d0,0.271d0,
     &     0.774d0,0.208d0,0.754d0,0.084d0,0.244d0,0.452d0,
     &     0.050d0,0.151d0,0.017d0,0.269d0,
     &     TE,Z,XI)

      DO I=1,5
        DO J=1,5
         EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      CONEIII=0.
      DO  I=1,5
        DO  J=I+1,5
          CONEIII=CONEIII+EMCA(J,I)
        ENDDO
      ENDDO
      do k=1,kmaxp(6,3)
         cinx(6,3,k)=weh(k)/del(ik)
         wlix(6,3,k)=wlin(k)
         ilabcfx(6,3,k)=603
      enddo

      FB(37)=EMCA(4,2)+EMCA(4,1)
      FB(38)=EMCA(5,2)
      FB(39)=EMCA(5,4)
      FB(46)=EMCA(2,1)
      FB(47)=EMCA(3,1)
      line_cool(6,3)=coneiii
C                                                                       
C       Ne IV    (1602, 2423-25, 4714-24)                                 
C                                                                       
      Z=ABn(6)*XNE(4)                                                     
      CALL FORB3(11,75,5.12D0,2.63D0,7.28D-1,1.02D0,2.59D-3,1.36D0                
     &,0.46D0,1.77D0,4.D0,10.D0,6.D0,FB(40),FB(41),FB(42),TE,Z,F)      
      DO NM=40,42                                                     
         FB(NM)=Z*FB(NM)
      enddo
      do k=1,3         
         cinx(6,4,k)=fb(k+39)
         if(k==1) then
            wlix(6,4,k)=wl21
         elseif(k==2) then
            wlix(6,4,k)=wl31
         elseif(k==3) then
            wlix(6,4,k)=wl32
         endif
         ilabcfx(6,4,k)=604
      enddo
      kmaxp(6,4)=3
      C1341=Z*F
      line_cool(6,4)=c1341
C                                                                       
C       Ne V    (1575, 2975, 3346-3426)                                
C                                                                       
      Z=ABn(6)*XNE(5)                                                     
C      CALL FORB3(10,76,3.71D0,4.17D0,5.2D-1,4.2D0,2.6D0,1.92D0,           
C     &4.2D0,2.6D0,9.D0,5.D0,1.D0,FB(43),FB(44),FB(45),TE,Z,F)      
      DO  NM=43,45                                                     
         FB(NM)=Z*FB(NM)
      enddo
C      C1351=Z*F
C     A:s FROM MENDOZA COLL. STRENGTHS FROM LENNON & BURKE MN 251, 628 (91)
C           AT 10,000 K.
c      CALL FIVELEV_dp(76,5,0.d0,412.4d0,1110.1d0,30291.5d0,63913.6d0,
c     &1.d0,3.d0,5.d0,5.d0,1.d0,1.28e-3,5.08E-9,2.37E-5,0.,4.59E-3,
c     &1.31E-1,4.21d0,0.365d0,6.69E-3,2.85d0,1.401d0,1.766d0,0.231d0,
c     &     0.0660d0,5.725d0,0.692d0,0.198d0,1.154d0,0.330d0,0.5939d0,TE,Z,XI)
c     New Ne V
      CALL FIVELEV_dp(76,6,5,5,0.d0,411.227d0,1109.5d0,30290.7d0,
     &     63915.4d0,1.d0,3.d0,5.d0,5.d0,1.d0,1.28d-3,5.08d-9,2.37d-5,
     &     0.d0,4.59d-3,1.31d-1,4.21d0,0.365d0,6.69d-3,2.85d0,
     &     1.41d0,1.81d0,0.232d0,0.027d0,5.83d0,0.695d0,0.082d0,1.159d0,
     &     0.137d0,0.58d0,TE,Z,XI)
      DO I=1,5
        DO J=1,5
         EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      CONEV=0.
      DO  I=1,5
        DO  J=I+1,5
          CONEV=CONEV+EMCA(J,I)
        ENDDO
      ENDDO
      FB(45)=EMCA(4,2)+EMCA(4,3)
      FB(43)=EMCA(5,2)
      FB(44)=EMCA(5,4)
      C1351=CONEV
      do k=1,kmaxp(6,5)
         cinx(6,5,k)=weh(k)
         wlix(6,5,k)=wlin(k)
         ilabcfx(6,5,k)=605
      enddo
      line_cool(6,5)=conev
C
C     NE VII 895.1
C
      Z=AB(13)*XNE(7)
      OM=0.172*T4**0.41
      CALL RLOSS(78,895.1D0,OM,1.D0,9.D0,6.6D3,RS,XEL,Z,TE,C1371
     &     ,W(75),CIN(75),FR(2,75),75,WLI(75))    
      
      cinx(6,7,1)=cin(75)
      wlix(6,7,1)=wli(75)
      ilabcfx(6,7,1)=607
      kmaxp(6,7)=1
      line_cool(6,7)=c1371
C     
C     SODIUM
C
C
C     NA I 5889 OM FROM KENNEDY ET AL J. PHYS. B 10, 3759
C
      Z=AB(10)*XNA(1)
      OM=10.45*T4**0.837
      CALL RLOSS(50,5889.D0,OM,2.D0,6.D0,1.07D8,RS,XEL,Z,TE,C1011
     &     ,W(22),CIN(22),FR(2,22),22,WLI(22))
      cinx(7,1,1)=cin(22)
      wlix(7,1,1)=wli(22)
      ilabcfx(7,1,1)=701
      kmaxp(7,1)=1
      line_cool(7,1)=c1011
C
C     NA I 5889 RECOMBINATION (ALL REC. COMES OUT IN 5889)
C
      CIN(41)=AB(10)*XNA(2)*2.52E-13*3.37E-12/T4**.682
      cinx(7,1,2)=cin(41)+cin(22)
      wlix(7,1,2)=5889
      ilabcfx(7,1,2)=701
      kmaxp(7,1)=1

C
C     MAGNESIUM
C
C     MG I 5177 DIELECTRONIC REC. (NUSSBAUMERY&.STOREY 1985)
C
      CALL DIELB(0.0d0,.1658d0,0.5114d0,-5.69d-2,2.1817d0,T4,ALDB)
      CIN(28)=AB(7)*XM(2)*ALDB*3.87E-12
      WLI(28)=5177.
      cinx(8,1,1)=cin(28)
      wlix(8,1,1)=wli(28)
      ilabcfx(8,1,1)=801
C
C     MG I 8806 DIELECTRONIC REC. (NUSSBAUMERY&.STOREY 1985)
C
      CALL DIELB(.1259d0,-.6903d0,1.1785d0,-9.60d-2,2.0651d0,T4,ALDB)
      CIN(29)=AB(7)*XM(2)*ALDB*2.27E-12
      WLI(29)=8806.
      cinx(8,1,2)=cin(29)
      wlix(8,1,2)=wli(29)
      ilabcfx(8,1,2)=801

C
C     MG I 3834 DIELECTRONIC REC. (NUSSBAUMERY&.STOREY 1985)
C
      CALL DIELB(.1283d0,-1.2733d0,4.6334d0,-4.73d-1,2.7691d0,T4,ALDB)
      CIN(30)=AB(7)*XM(2)*ALDB*5.22E-12
      WLI(30)=3834.
      cinx(8,1,3)=cin(30)
      wlix(8,1,3)=wli(30)
      ilabcfx(8,1,3)=801
C
C     MG I 4571 OM FROM FABRIKANT (J. PHYS. B 7:91 -74), AS GIVEN
C          BY OSTERBROCK (-74) P 245.
C               A FROM
C
      Z=AB(7)*XM(1)
      OM=1.60*T4**0.56
      CALL RLOSS(39,4572.D0,OM,1.D0,9.D0,7.20D1,RS,XEL,Z,TE,C511
     &     ,W(20),CIN(20),FR(2,20),20,WLI(20))
      cinx(8,1,4)=cin(20)
      wlix(8,1,4)=wli(20)
      ilabcfx(8,1,4)=801
C
C     MG I 4572 RECOMBINATION (SEE NOTES). DIEL. FROM NUSSBAUMER & STORE
C
      CALL DIELB(.5116d0,-2.8906d0,7.445d0,-.72340d0,2.414d0,T4,ALDB)
      CIN(40)=AB(7)*XM(2)*(1.96E-13/T4**0.855+ALDB)*4.34E-12
      cinx(8,1,5)=cin(40)
      wlix(8,1,5)=4572.
      ilabcfx(8,1,5)=801
C
C     MG I 2852 OM FROM FABRIKANT (J. PHYS. B 7:91)
C
      Z=AB(7)*XM(1)
      OM=2.11*T4**1.187
      CALL RLOSS(39,2852.D0,OM,1.D0,3.D0,4.93D8,RS,XEL,Z,TE,C512
     &     ,W(21),CIN(21),FR(2,21),21,WLI(21))
      cinx(8,1,6)=cin(21)
      wlix(8,1,6)=wli(21)
      ilabcfx(8,1,6)=801
      kmaxp(8,1)=6
      line_cool(8,1)=c511+c512
C
C     MG II 2800 OM (MENDOZA 83)
C
      Z=AB(7)*XM(2)
      om=16.9*t4**0.1
      CALL RLOSS(40,2795.D0,om,2.D0,6.D0,2.55D8,RS,XEL,Z,TE,C521
     &     ,W(15),CIN(15),FR(2,15),15,WLI(15))
      cinx(8,2,6)=cin(21)
      wlix(8,2,6)=wli(21)
      ilabcfx(8,2,6)=801
      kmaxp(8,2)=1
      line_cool(8,2)=c521

C
C     Al II 1671 Pradhan compli.
C
      Z=ABn(9)*XAL(2)
      OM=3.251*T4**0.112
      CALL RLOSS(47,1671.D0,OM,1.D0,3.D0,1.46D9,RS,XEL,Z,TE,c921
     &     ,W(42),CIN(42),FR(2,42),42,WLI(42))
      cinx(9,2,1)=cin(42)
      taulinex(9,2,1)=t0
      wlix(9,2,1)=wli(42)
      ilabcfx(9,2,1)=902
C
C     AL II 2670. Pradhan comp.
C
      Z=ABn(9)*XAL(2)
      om=3.251*t4**0.25
      CALL RLOSS(47,2670.D0,om,1.D0,9.D0,1.11d3,RS,XEL,Z,TE,c922
     &     ,W(81),CIN(81),FR(2,81),81,WLI(81))
      cinx(9,2,2)=cin(81)
      taulinex(9,2,2)=t0
      wlix(9,2,2)=wli(81)
      ilabcfx(9,2,2)=902
      kmaxp(9,2)=2
      line_cool(9,2)=c921+c922

C
C     AL III 1854.7 1862.8 OM and A Kingdon 
C
      Z=ABn(9)*XAL(3)
      om=16.0*t4**0.1
      CALL RLOSS(48,1858.D0,om,2.D0,6.D0,5.55D8,RS,XEL,Z,TE,c931
     &     ,W(96),CIN(96),FR(2,96),96,WLI(96))
      cinx(9,3,1)=cin(96)
      taulinex(9,3,1)=t0
      wlix(9,3,1)=wli(96)
      ilabcfx(9,3,1)=903
      kmaxp(9,3)=1
      line_cool(9,3)=c931
C
C     SILICON
C
C     (SI I) 6591, 10995, 16360 A. MENDOZA. OMEGA FOR 3P - 1D FROM
C     M.S. PINDZOLA, A.K. BHATIA, A.TEMKIN, PHYS.REV. 15, 35 (1977)
C     REST SCALED FROM C I WITH SCALING = 3
      O21=3.97*T4**0.9
      O31=3.*0.149*(TE/5000.)**0.871
      O32=3.*0.196*(TE/5000.)**0.499
      Z=ABN(10)*XSi(1)
      CALL FORB3(5,25,0.7580D0,1.1277D0,3.043D-3,3.22D-2,1.14D0,O21
     &,O31,O32,9.D0,5.D0,1.D0,FB(30),FBQQ,FB(31),TE,Z,F)
      C412=ABn(10)*XSi(1)*F
      DO NM=30,31
         FB(NM)=ABn(10)*XSi(1)*FB(NM)
      enddo
      do k=1,3
         if(k==1) then
            wlix(10,1,k)=wl21
            cinx(10,1,k)=fb(30)
         elseif(k==2) then
            wlix(10,1,k)=wl31
            cinx(10,1,k)=ABn(10)*XSi(1)*fbqq
         elseif(k==3) then
            wlix(10,1,k)=wl32
            cinx(10,1,k)=fb(31)
         endif
         ilabcfx(10,1,k)=1001
      enddo



C     SI I 129.68 68.474 MY A:s FROM MAZDA OMEGAS GUESSED!! FROM O I
      O32=9.76E-6*(TE-228)+3.46E-11*(TE-228.)**2
      O31=3.39E-6*(TE-326.)-2.90E-11*(TE-326.)**2
      O21=1.89E-6*(TE-326.)+8.00E-11*(TE-326.)**2
      IF(O21.LT.0.) O21=1.e-10
      IF(O31.LT.0.) O31=1.e-10
      IF(O32.LT.0.) O32=1.e-10
      Z=ABn(10)*XSi(1)
      CALL FORB3(15,25,9.559D-3,1.810D-2,8.25D-6,0.D0,4.21D-5,O21,O31,
     &O32,1.D0,3.D0,5.D0,FB(28),FBQQ,FB(29),TE,Z,F)
      DO  NM=28,29
         FB(NM)=ABN(10)*XSi(1)*FB(NM)
      enddo
      C411=ABN(10)*XSi(1)*F
      do k=4,6
         if(k==4) then
            wlix(10,1,k)=wl21
            cinx(10,1,k)=fb(28)
         elseif(k==5) then
            wlix(10,1,k)=wl31
            cinx(10,1,k)=ABN(10)*XSi(1)*fbqq
         elseif(k==6) then
            wlix(10,1,k)=wl32
            cinx(10,1,k)=fb(29)
         endif
         ilabcfx(10,1,k)=1001
      enddo
      kmaxp(10,1)=6
      line_cool(10,1)=c411+c412
      
C
C     SI II) 2335 (MENDOZA 83)
C
      Z=ABn(10)*XA(2,26)
      CALL RLOSS(26,2335.D0,5.14D0,6.D0,12.D0,4.14D1,RS,XEL,Z,TE,C421
     &     ,W(17),CIN(17),FR(2,17),17,WLI(17))
      cinx(10,2,1)=cin(17)
      wlix(10,2,1)=wli(17)
      ilabcfx(10,2,1)=1002
C     SI II 34.81 MY A FROM(MENDOZA)
C     omega from Keenan et al MNRAS 214, 37p, 1985.
      Z=ABn(10)*XSi(2)
      OM=5.6
      CALL RLOSS(26,3.481D5,OM,2.D0,4.D0,2.17D-4,RS,XEL,Z,TE,C422
     &     ,W(44),CIN(44),FRQQ,44,WLI(44))
      cinx(10,2,2)=cin(44)
      wlix(10,2,2)=wli(44)
      ilabcfx(10,2,2)=1002
      kmaxp(10,2)=2
      line_cool(10,2)=c421+c422
C
C     SI III 1892 OMEGA FROM MENDOZA 83
C
      Z=ABn(10)*XA(2,27)
      OM=5.43/T4**0.34
      CALL RLOSS(27,1892.D0,OM,2.D0,9.D0,4.20D3,RS,XEL,Z,TE,C431
     &     ,W(73),CIN(73),FR(2,73),72,WLI(73))
      cinx(10,3,1)=cin(73)
      wlix(10,3,1)=wli(73)
      ilabcfx(10,3,1)=1003
      kmaxp(10,3)=1
      line_cool(10,3)=c431

C
C     SI IV 1400 OMEGA=21.?
C
      Z=ABN(10)*xsi(4)
      CALL RLOSS(28,1400.D0,17.D0,2.D0,6.D0,9.22D8,RS,XEL,Z,TE,C441
     &     ,W(16),CIN(16),FR(2,16),16,WLI(16))
      cinx(10,4,1)=cin(16)
      wlix(10,4,1)=wli(16)
      ilabcfx(10,4,1)=1004
      kmaxp(10,4)=1
      line_cool(10,4)=c441
C
C     SULPHUR
C
C     S I 25.25 , 56.31 MY A:s FROM MAZDA OMEGAS GUESSED!! FROM O I
      O32=9.76E-6*(TE-228)+3.46E-11*(TE-228.)**2
      O31=3.39E-6*(TE-326.)-2.90E-11*(TE-326.)**2
      O21=1.89E-6*(TE-326.)+8.00E-11*(TE-326.)**2
      IF(O21.LT.0.) O21=1.e-10
      IF(O31.LT.0.) O31=1.e-10
      IF(O32.LT.0.) O32=1.e-10
      Z=ABn(11)*XSU(1)
      CALL FORB3(16,56,4.912D-2,2.202D-2,1.39D-3,6.71D-8,3.02D-4,O21,
     &O31,O32,5.D0,3.D0,1.D0,FB(26),FBQQ,FB(27),TE,Z,F)
      DO  NM=26,27
         FB(NM)=ABn(11)*XSU(1)*FB(NM)
      enddo
      C1211=ABn(11)*XSU(1)*F
      do k=1,3
         if(k==1) then
            wlix(11,1,k)=wl21
            cinx(11,1,k)=fb(26)
         elseif(k==2) then
            wlix(11,1,k)=wl31
            cinx(11,1,k)=AB(12)*XSu(1)*fbqq
         elseif(k==3) then
            wlix(11,1,k)=wl32
            cinx(11,1,k)=fb(27)
         endif
         ilabcfx(11,1,k)=1101
      enddo
      kmaxp(11,1)=3
      line_cool(11,1)=c1211


C     *************************************************************
C     *****
c   S II 43 level atom      
C     *****
C     *************************************************************
      ilcf=1102
      call popchianti(11,2,te,XQ,cos2)

      do k=1,kmaxp(11,2)
         ll=ll+1
         if(ipopch.eq.1) then
            cinx(11,2,k)=weh(k)/del(ik)
            taulinex(11,2,k)=taul(k)
            wlix(11,2,k)=wlin(k)
         endif
         ilabcfx(11,2,k)=1102
      enddo
      line_cool(11,2)=cos2

C     *************************************************************
C     *****
c   S III 49 level atom      
C     *****
C     *************************************************************
      ilcf=1103
      z=abn(11)*xsu(3)
      
      call popchianti(11,3,te,XQ,cos3)
      kmaxp(11,3)=kmaxp(11,3)+1
      do k=1,kmaxp(11,3)
         cinx(11,3,k)=weh(k)/del(ik)
         taulinex(11,3,k)=taul(k)
         wlix(11,3,k)=wlin(k)
         ilabcfx(11,3,k)=1103
      enddo
      line_cool(11,3)=cos3

c S IV 
      ilcf=1104
      z=abn(11)*xsu(4)
      call popchianti_new(11,4,te,cos4)
      do k=1,kmaxp(11,4)
         if(ipopch.eq.1) then
            cinx(11,4,k)=weh(k)/del(ik)
            taulinex(11,4,k)=taul(k)
            wlix(11,4,k)=wlin(k)
         endif
         ilabcfx(11,4,k)=1104
      enddo
      line_cool(11,4)=cos4

c S V 
      ilcf=1105
      z=abn(11)*xsu(5)
      call popchianti_new(11,5,te,cos5)
      do k=1,kmaxp(11,5)
         if(ipopch.eq.1) then
            cinx(11,5,k)=weh(k)/del(ik)
            taulinex(11,5,k)=taul(k)
            wlix(11,5,k)=wlin(k)
         endif
         ilabcfx(11,5,k)=1105
      enddo
      line_cool(11,5)=cos5

C     S V 786.5
C
C      Z=AB(12)*XA(2,60)
C      CALL RLOSS(60,786.5D0,OM,G1,G2,A21,RS,XEL,Z,TE,C1251
C     &,W(76),CIN(76),FR(2,76),76,WLI(76)) 
C                                                                       
C     S VI 937.1     (G & S)
C
      Z=AB(12)*XA(2,61)
      OM=9.65D0*T4**0.156
      CALL RLOSS(61,937.1D0,OM,2.D0,6.D0,1.61D9,RS,XEL,Z,TE,C1261
     &     ,W(77),CIN(77),FR(2,77),77,WLI(77))
      wlix(11,6,2)=wli(77)
      cinx(11,6,2)=cin(77)
      ilabcfx(11,6,2)=1106
      kmaxp(11,6)=1
      line_cool(11,6)=c1261

      
C     AR II 6.985 MY from Pelan & Berrington 1995
      Z=AB(14)*XAR(2)
      OM=0.635
      telog=log10(te)
      if(te < ar2_te(1)) then
         al_ar2=log10(ar2_coll(2)/ar2_coll(1))/log10(ar2_te(2)/ar2_te(1))
         om=ar2_coll(1)*(te/ar2_te(1))**al_ar2
         ii=1
      elseif(te > ar2_te(11)) then
         al_ar2=log10(ar2_coll(11)/ar2_coll(10))/log10(ar2_te(11)/ar2_te(10))
         om=ar2_coll(11)*(te/ar2_te(11))**al_ar2
      else
         do i=1,10
            if(te > ar2_te(i) .and. te <= ar2_te(i+1)) then
               al_ar2=log10(ar2_coll(i+1)/ar2_coll(i))/log10(ar2_te(i+1)/ar2_te(i))
               om=ar2_coll(i)*(te/ar2_te(i))**al_ar2
               ii=i
            endif
         enddo
      endif
      CALL RLOSS(60,6.985D4,OM,4.D0,2.D0,5.27D-2,RS,XEL,Z,TE,C1422
     &     ,W(45),CIN(45),FRQQ,45,WLI(45))
      wlix(12,2,1)=wli(45)
      cinx(12,2,1)=cin(45)
      ilabcfx(12,2,1)=1202
      kmaxp(12,2)=1
      line_cool(12,2)=c1422
      
c!!! Ar III-IV New
C
C     Ar III 3P,1D,1S OMEGAS and A:s from Pradhan
C
      om21=2.24e0
      om31=0.531e0
      om32=1.18e0
      om1d3p=4.74e0
      om43=1.*om1d3p/9.
      om42=3.*om1d3p/9.
      om41=5.*om1d3p/9.
      om1s3p=0.680e0
      om53=1.*om1s3p/9.
      om52=3.*om1s3p/9.
      om51=5.*om1s3p/9.
      om54=0.823e0
      z=ab(14)*xar(3)
      CALL FIVELEV_dp(96,12,3,5,0.d0,1112.175d0,1570.229d0,14010.004d0,
     &     33265.724d0,5.d0,3.d0,1.d0,5.d0,1.d0,
     &     3.08d-2,2.37d-6,3.14d-1,4.17d-2,5.17d-3,8.23d-2,3.91d0,
     &     2.21d-5,0.d0,2.59d0,
     &     2.d0,0.53d0,0.5255d0,3.77d-1,1.18d0,1.58d0,2.266d-1,
     &     0.5266d0,7.55d-2,0.823d0,TE,Z,XI)
      DO I=1,5
        DO J=1,5
         EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      coarIII=0.
      DO  I=1,5
        DO  J=I+1,5
          COariii=COariii+EMCA(J,I)
        ENDDO
      ENDDO


      fb(55)=emca(2,1)
      fb(56)=emca(3,2)
      fb(57)=emca(4,1)
      fb(58)=emca(4,2)
      fb(59)=emca(5,1)
      fb(60)=emca(5,2)
      fb(61)=emca(5,4)

      do k=1,kmaxp(12,3)
         wlix(12,3,k)=wlin(k)
         cinx(12,3,k)=weh(k)/(del(ik))
         ilabcfx(12,3,k)=1203
      enddo
      line_cool(12,3)=coariii
            
C     Ar IV 3P,1D,1S OMEGAS from Ramsbottom et al '97 and A:s from Pradhan
C
      om21=1.144*t4**0.027
      om31=0.762*t4**0.028
      om32=7.055*t4**0.086
      om1d3p=2.29*t4**0.132
      om41=0.785*t4**0.140
      om42=3.939*t4**0.025
      om43=2.139*t4**0.027
      om1s3p=0.293*t4**0.167
      om51=0.393*t4**0.140
      om52=1.533*t4**0.028
      om53=1.507*t4**0.018
      om54=2.065*t4**0.383
      z=abn(12)*xar(4)

      CALL FIVELEV_dp(97,12,4,5,0.d0,21090.4d0,21219.3d0,34855.5d0,
     &     35032.6d0,4.d0,4.d0,6.d0,2.d0,4.d0,
     &     1.77d-3,2.23d-2,2.11d0,0.862d0,2.30d-5,0.598d0,
     &     0.119d0,0.789d0,0.603d0,4.94d-5,
     &     om21,om31,om41,om51,om32,om42,om52,
     &     om43,om53,om54,TE,Z,XI)
      DO I=1,5
        DO J=1,5
         EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      coarIv=0.
      DO  I=1,5
        DO  J=I+1,5
          COariv=COariv+EMCA(J,I)
        ENDDO
      ENDDO

      do k=1,kmaxp(12,4)
         wlix(12,4,k)=wlin(k)
         cinx(12,4,k)=weh(k)/(del(ik))
         ilabcfx(12,4,k)=1204
      enddo
      line_cool(12,4)=coariv
      
C     ************************************************************
C     *****
c   Ar V 5 level atom      
C     *****
C     *************************************************************
      ilcf=1205
c     Omegas from  Galavis, M.~E., Mendoza, C., \& Zeippen, C.~J.\ 1995, \aaps, 111, 347
      telog=log10(te)
      n=5
      nte=11
      do i=1,n
         DO J=i+1,N
            do k=1,nte-1 

               if(telog.gt.tefitar(nte)) then
                  omint = omfitar(j,i,nte)
               elseif(telog.lt.tefitar(1)) then
                  omint = omfitar(j,i,1)
               elseif(telog.gt.tefitar(k).and.
     &                 telog.le.tefitar(k+1)) 
     &                 then
                  omint = omfitar(j,i,k) + (telog-tefitar(k))*
     &                 (omfitar(j,i,k+1)-omfitar(j,i,k))/
     &                 (tefitar(k+1)-tefitar(k))

               endif
            enddo
            omar(j,i)=omint
         enddo
      enddo
      z=abn(12)*xar(5)
      CALL FIVELEV_dp(98,12,5,5,0.d0,765.23d0,2028.80d0,16298.9d0,
     &     37912.0d0,1.d0,3.d0,5.d0,5.d0,1.d0,
     &     8.05d-3,1.33d-6,4.9d-5,0.d0,
     &     2.73d-2,2.23d-1,6.8d0,
     &     5.2d-1,8.1d-2,
     &     3.8d0,
     &     omar(2,1),omar(3,1),omar(4,1),omar(5,1),omar(3,2),omar(4,2),
     &     omar(5,2),omar(4,3),omar(5,3),omar(5,4),TE,Z,XI)
      DO I=1,5
        DO J=1,5
         EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      do k=1,kmaxp(12,5)
         wlix(12,5,k)=wlin(k)
         cinx(12,5,k)=weh(k)/(del(ik))
         ilabcfx(12,5,k)=1205
      enddo
      
      coarv=0.
      DO  I=1,5
        DO  J=I+1,5
          coarv=coarv+EMCA(J,I)
        ENDDO
      ENDDO
      line_cool(12,5)=coarv

C     ************************************************************
C     *****
c   Ar VI 2 level atom      
C     *****
C     *************************************************************
      Z=AB(14)*xar(6)
      OM=0.303*t4**0.065
      telog=log10(te)
c Omegas from Saraph, H. E., & Storey, P. J. 1996, A&AS, 115, 151      
      if(te < ar6_te(1)) then
         al_ar6=log10(ar6_coll(2)/ar6_coll(1))/log10(ar6_te(2)/ar6_te(1))
         om=ar6_coll(1)*(te/ar6_te(1))**al_ar6
         ii=1
      elseif(te > ar6_te(43)) then
         al_ar6=log10(ar6_coll(43)/ar6_coll(42))/log10(ar6_te(43)/ar6_te(42))
         om=ar6_coll(43)*(te/ar6_te(43))**al_ar6
      else
         do i=1,42
            if(te > ar6_te(i) .and. te <= ar6_te(i+1)) then
               al_ar6=log10(ar6_coll(i+1)/ar6_coll(i))/log10(ar6_te(i+1)/ar6_te(i))
               om=ar6_coll(i)*(te/ar6_te(i))**al_ar6
               ii=i
            endif
         enddo
      endif
      CALL RLOSS(99,4.52922D4,OM,2.D0,4.D0,9.7D-2,RS,XEL,Z,TE,coarvi
     &     ,War6,CINar6,FRQQ,46,WLIar6)
      cinx(12,6,1)=cinar6
      taulinex(12,6,1)=t0
      wlix(12,6,1)=wliar6
      ilabcfx(12,6,1)=1206
      kmaxp(12,6)=1
      line_cool(12,6)=coarvi

C     ************************************************************
C     *****
c   Ca IV 3 level atom up to 3p6 2S     
C     *****
C     *************************************************************

c     interpolate coll strength
      if(te.lt.teca4(1)) then
         o21=coll_ca4(1,1)
         o31=coll_ca4(2,1)
         o32=coll_ca4(3,1)
      elseif(te.ge.teca4(50)) then
         o21=coll_ca4(1,50)
         o31=coll_ca4(2,50)
         o32=coll_ca4(3,50)
      else
         do kc=1,49
            if(te>=teca4(kc) .and. te<=teca4(kc+1)) then
               do k=1,3
                  omca4(k)=coll_ca4(k,kc) + (te-teca4(kc))*
     &                 (coll_ca4(k,kc+1)-coll_ca4(k,kc))/(teca4(kc+1)-teca4(kc))
               enddo
               om21=omca4(1)
               om31=omca4(2)
               om32=omca4(3)
            endif
         enddo
      endif
      Z=ABn(13)*xca(4)
      e21=0.38661d0
      e32=18.90010d0-e21
      CALL FORB3(18,114,e21,e32,0.545d0,5.11d8,2.466d8,om21,om31,om32, 
     &     4.D0,2.D0,2.D0,FB21,FB31,FB32,TE,Z,F)      
      DO NM=64,66                                                     
         FB(NM)=Z*FB(NM)
      enddo
      do k4=1,3         
         if(k4==1) then
            cinx(13,4,k4)=z*fb21
            wlix(13,4,k4)=wl21
         elseif(k4==2) then
            cinx(13,4,k4)=z*fb31                     
            wlix(13,4,k4)=wl31
         elseif(k4==3) then
            cinx(13,4,k4)=z*fb32            
            wlix(13,4,k4)=wl32
         endif
         ilabcfx(13,4,k4)=1304
      enddo
      kmaxp(13,4)=3
      Cocaiv=Z*F
      cocal=cocal+cocaiv
      line_cool(13,4)=cocaiv


C
C     Ca V
C
      om21=2.8e0
      om31=0.66e0
      om32=0.996e0
      om1d3p=3.1e0
      om43=1.*om1d3p/9.
      om42=3.*om1d3p/9.
      om41=5.*om1d3p/9.
      om1s3p=0.147e0
      om53=1.*om1s3p/9.
      om52=3.*om1s3p/9.
      om51=5.*om1s3p/9.
      om54=1.09e0

      call upsilon_chianti(0,17,13,5,5,te)
      om21=upsil(17,1,2)
      om31=upsil(17,1,3)
      om41=upsil(17,1,4)
      om51=upsil(17,1,5)
      om32=upsil(17,2,3)
      om42=upsil(17,2,4)
      om52=upsil(17,2,5)
      om43=upsil(17,3,4)
      om53=upsil(17,3,5)
      om54=upsil(17,4,5)
      z=abn(13)*xca(5)

      CALL FIVELEV_dp(113,13,5,5,0.d0,2404.7d0,3275.6d0,18830.3d0,
     &     43836.5d0,5.d0,3.d0,1.d0,5.d0,1.d0,
     &     3.079d-01,3.960d-05,1.955d+00,1.224d-01,3.681d-02,4.356d-01,
     &     2.345d+01,7.986d-05,0.d0, 3.697d+00,
     &     om21,om31,om41,om51,om32,om42,om52,
     &     om43,om53,om54,TE,Z,XI)
      DO I=1,5
        DO J=1,5
         EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      coca5=0.
      DO  I=1,5
        DO  J=I+1,5
          COca5=COca5+EMCA(J,I)
        ENDDO
      ENDDO
      do k=1,kmaxp(13,5)
         wlix(13,5,k)=wlin(k)
         cinx(13,5,k)=weh(k)/(del(ik))
         ilabcfx(13,5,k)=1305
      enddo
      line_cool(13,5)=coca5
      cocal=cocal+coca5
      
C
C     FE I FIT TO AXELROD  AT N=1.E7
C
      W(18)=ABN(14)*XF(1)*5.1E-19*EXP(-2.83/T4)/SQRT(T4)
C     FINE STRUCTURE
      FEIFS=1.15E-23*EXP(-0.6499/T3)/(T3**0.358*(DEN(IK)*XEL/1.E7)**
     &                                                         0.96)
      FEIFS=FEIFS*ABN(14)*XF(1)
      WLI(18)=-1.
      CIN(18)=W(18)+FEIFS
      
C
C     FE II FIT TO AXELROD AT N=1.E7
C
      W(19)=ABN(14)*XF(2)*3.7E-19*EXP(-2.91/T4)/SQRT(T4)
C     FINE STRUCTURE
      FEIIFS=1.03E-23*EXP(-0.5077/T3)/(T3**0.162*(DEN(IK)*XEL/1.E7)**
     &                                                         0.96)
      FEIIFS=FEIIFS*ABN(14)*XF(2)
      WLI(19)=-1.
      CIN(19)=W(19)+FEIIFS
C     FE III

      IF(ABN(14)*XF(3).GT.1.E-6) THEN
         denel=den(ik)*xel
         z=abn(14)*xf(3)
         if(ilowion.eq.0) then
            call feiiiatom(te,denel,z,rtot)
         elseif(ilowion.eq.1) then

C     
C     FE III FIT TO AXELROD AT N=1.E7
C     
            CIN42=ABN(14)*XF(3)*3.7E-19*EXP(-31.82/T3)/T3**0.634
            FEIIIFS=2.92E-21*EXP(-0.6313/T3)/T3**0.0313
            FEIIIFS=FEIIIFS/(1.+(DEN(IK)*XEL)/4.04E4)
            FEIIIFS=FEIIIFS*ABN(14)*XF(3)
            CIN42=CIN42+FEIIIFS      
         endif
         feiii=rtot*abn(14)*xf(3)
         cin42=feiii
      endif

      WLI42=-1.
C     *************************************************************
C     *****
c   Fe III 55 level atom (should be larger!)
C     *****
C     *************************************************************
      ilcf=1402
      call popchianti_new(14,3,te,cofe3)
      do k=1,kmaxp(14,3)
         if(ipopch.eq.1) then
            cinx(14,3,k)=weh(k)/del(ik)
            taulinex(14,3,k)=taul(k)
            wlix(14,3,k)=wlin(k)
         endif
         ilabcfx(14,3,k)=1403
      enddo
      line_cool(14,3)=cofe3

      do io=3,9
         if(io==3) then
            ilcf=1403
         elseif(io==4) then
            ilcf=1404
         elseif(io==5) then
            ilcf=1405
         elseif(io==6) then
            ilcf=1406
         elseif(io==7) then
            ilcf=1407
         elseif(io==8) then
            ilcf=1408
         elseif(io==9) then
            ilcf=1409
         endif

         if(te>500. .or. io<=4) then
            call popchianti_new(14,io,te,cofeio)
            do k=1,kmaxp(14,io)
               if(ipopch.eq.1) then
                  cinx(14,io,k)=weh(k)/del(ik)
                  taulinex(14,io,k)=taul(k)
                  wlix(14,io,k)=wlin(k)
               endif
               ilabcfx(14,io,k)=ilcf
            enddo
         else
            do k=1,kmaxp(14,io)
               if(ipopch.eq.1) then
                  cinx(14,io,k)=0.
                  taulinex(14,io,k)=0.
                  wlix(14,io,k)=wlin(k)
               endif
               ilabcfx(14,io,k)=ilcf
               cofeio=0.
            enddo
         endif
         if(io==3) then
            cofe3=cofeio
         elseif(io==4) then
            cofe4=cofeio
         elseif(io==5) then
            cofe5=cofeio
         elseif(io==6) then
            cofe6=cofeio
         elseif(io==7) then
            cofe7=cofeio
         elseif(io==8) then
            cofe8=cofeio
         elseif(io==9) then
            cofe9=cofeio
         endif
         line_cool(14,io)=cofeio
      enddo

C     *************************************************************
C     *****
C     CALCULATE POPULATIONS AND LINE EMISSION FROM O I
C     *****
C     *************************************************************
5005  ZEL=DEL(IK)-AB(3)*XO(2)
      IFPOP=0
      XQ=XO(1)+XO(2)
      DO I=1,65
         SR(I)=0.
 3976    WEH(I)=0.
      enddo
      DO I=1,12
         DO J=I+1,12
            EMCA(J,I)=0.
         enddo
      enddo
      WLI(31)=6300.
      WLI(32)=2972.
      WLI(33)=5572.
      WLI(34)=1356.
      WLI(35)=1302.
      WLI(36)=7774.
      WLI(37)=8446.
      WLI(38)=1027.
      WLI(39)=990.
      TOI=TE
      IF(TE.LT.TLIMO) TOI=TLIMO
C      IF(AB(3).GT.1.E-5) THEN

6901  format(' O I-III ',10e12.4)
      IF(AB(3)*XO(1).GT.1.E-10) then
         call POPO_new(RS,Toi,zel,XQ,PHO,IFPOP)
      endif
      if(xo(2).gt.0.1) then
         recto=0.
         do k=1,9
            recto=recto+reco(k)
         enddo
      endif
      IF(IFPOP.EQ.1) GOTO 8278
      DO IJ=1,4
         SRED(IJ+30)=SR(IJ)
         CIN(IJ+30)=WEH(IJ)/DEL(IK)
      enddo
      CIN(35)=WEH(7)/DEL(IK)
      CIN(36)=WEH(14)/DEL(IK)
      CIN(37)=WEH(20)/DEL(IK)
      CIN(38)=WEH(22)/DEL(IK)
      CIN(39)=WEH(29)/DEL(IK)

      SRED(35)=SR(7)
      SRED(36)=SR(14)
      SRED(37)=SR(20)
      SRED(38)=SR(22)
      SRED(39)=SR(29)

      do k=1,kmaxp(5,1)
         wlix(5,1,k)=wlin(k)
         cinx(5,1,k)=weh(k)/(den(ik)*del(ik))
         ilabcfx(5,1,k)=501
      enddo
      COOI=0.
      DO I=1,8
         DO J=I+1,9
            COOI=COOI+EMCA(J,I)
         ENDDO
      ENDDO
      line_cool(5,1)=cooi
      IF(TE.LT.TLIMO) THEN
C     LOW TEMPERATURE APPROX.
         T3OI=TOI/1.E3
         O21=.0151*T3OI**1.31
         O31=.00184*T3OI**1.32
         O32=.031*T3OI**.534
         Z=AB(3)*XO(1)
         CALL FORB3(4,14,1.957D0,2.223D0,8.45D-3,7.35D-2,1.22D0,O21,
     &        O31,O32,9.D0,5.D0,1.D0,FBO(10),FBO(11),FBO(12),TOI,Z,F)
         DO NM=10,12
            FBO(NM)=AB(3)*XO(1)*FBO(NM)
         ENDDO
         CALL OXREC(OXR,TE)
         DO IJ=1,8
            OXR(IJ)=AB(3)*XA(2,15)*OXR(IJ)
         enddo
         CIN(31)=FB(10)*CIN(31)/FBO(10)
         CIN(32)=FB(11)*CIN(32)/FBO(11)
         CIN(33)=FB(12)*CIN(33)/FBO(12)
         CIN(34)=OXR(5)
         CIN(35)=OXR(1)
         CIN(36)=OXR(6)
         CIN(37)=OXR(2)
         CIN(38)=0.
         CIN(39)=0.
         DO IJ=1,8
            SRED(IJ+30)=0.
         enddo
         COOI=COOI*(FB(10)+FB(11)+FB(12))/(FBO(10)+FBO(11)+FBO(12))
         line_cool(5,1)=cooi
      ENDIF
CF1016
C     RECOMB. EMISSION
      RECEMO=0.
      NP1H=10
      DO J=1,9
         CALL RECOMB(2,J,TE,ALOXI,RI)      
         RECEM(2,J)=AB(3)*BOL(NP1H)*ALOXI*(E00-EI(J))*1.602e-12
         RECEMO=RECEMO+RECEM(2,J)
      ENDDO
      C111=COOI
c!!!! 333200!!
7327  IF(AB(3).LT.1.E-5.OR.TE.LT.-333200.) goto 7127
      C111=0.      
      DO IJ=31,39
         C111=C111+CIN(IJ)
      enddo
C     *************************************************************
C     *****
C     CALCULATE POPULATIONS AND LINE EMISSION FROM CA II
C     *****
C     *************************************************************
 7127 ZEL=DEL(IK)-AB(11)*XCA(2)
      IFPOP=0
      WEH(1)=0.
      WEH(2)=0.
      WEH(3)=0.
      EMCA(2,1)=0.
      EMCA(3,1)=0.
      EMCA(3,2)=0.
      SR(1)=0.
      SR(2)=0.
      SR(3)=0.
      WLI(23)=7291.
      WLI(24)=3950.
      WLI(25)=8662.
      XQ=XCA(2)+XCA(3)
      TCA=TE
      IF(TE.LT.TLIMCA) TCA=TLIMCA
      IF(AB(11)*XCA(2).GT.1.E-7) THEN
          CALL POPCA(RS,TCA,ZEL,XQ,PHCA,IFPOP)

        IF(IFPOP.EQ.1) THEN
          DTQ=(TOLDCA-TCA)/11.
          TQ=TOLDCA
          DO KQ=1,11
            TQ=TQ-DTQ
            IF(AB(11)*XCA(2).GT.1.E-7) then

               CALL POPCA(RS,TQ,ZEL,XQ,PHCA,IFPOP)
            endif
          ENDDO
       ENDIF

        IF(IFPOP.EQ.1) GOTO 8278
        CIN(23)=WEH(1)/DEL(IK)
        CIN(24)=WEH(2)/DEL(IK)
        CIN(25)=WEH(3)/DEL(IK)
        SRED(23)=SR(1)
        SRED(24)=SR(2)
        SRED(25)=SR(3)

        do k=1,kmaxp(13,2)
           wlix(13,2,k)=wlin(k)
           cinx(13,2,k)=weh(k)/(den(ik)*del(ik))
           ilabcfx(13,2,k)=1302
        enddo
        
        COCAL=EMca(2,1)+EMCA(3,1)+EMCA(3,2)
        line_cool(13,2)=cocal
      ELSE
        IF(TE.LT.TLIMCA) THEN
        Z=AB(11)*XCA(2)
        CALL RLOSS(51,7291.D0,14.D0,2.D0,10.D0,1.3D0,RS,
     &         XEL,Z,TE,CXXX,W(23),CIN23,FR(2,23),23,WLI(23))
        CALL RLOSS(51,7291.D0,14.D0,2.D0,10.D0,1.3D0,RS,
     &         XEL,Z,TCA,CXCA,W(23),CIN23CA,FR(2,23),23,WLI(23))
        cin(23)=cin23*CIN(23)/CIN23CA
  	  EMCA(2,1)=CXXX*COCAL/CXCA
          COCAL=EMCA(2,1)
          line_cool(13,2)=cocal
        ENDIF
      ENDIF
      CALL RLOSS(51,7291.D0,14.D0,2.D0,10.D0,1.3D0,RS,
     &     XEL,Z,TE,CXXX,W23,CIN23,FR23,23,WLI23)
      
      IF(AB(11)*XCA(2).LT.1.E-7) THEN
        CIN(23)=0.
        CIN(24)=0.
        CIN(25)=0.
        SRED(23)=0.
        SRED(24)=0.
        SRED(25)=0.
        COCAL=0.
      ENDIF
C     *************************************************************
C     *****
C     CALCULATE POPULATIONS AND LINE EMISSION FROM Si I
C     *****
C     *************************************************************
      ZEL=DEL(IK)-ABn(10)*(Xsi(1)+2.*Xsi(2))
      IFPOP=0

      DO I=2,Nsi1
        DO J=1,I-1
          WOBS(I,J)=0.
          EMCA(I,J)=0.
        ENDDO
      ENDDO
      XFEII=1.-XF(1)-XF(3)-XF(4)

      XQ=Xsi(2)+Xsi(3)

      call POPsi_i(RS,TE,ZEL,XQ,PHHE,IFPOP)

      COSI1=0.
      SI1TOT=0.
      DO I=1,NSI1
         DO J=1,NSI1
c            WOBSI1(J,I)=WOBS(J,I)
            SI1TOT=SI1TOT+WOBS(J,I)/DEL(IK)
            COSI1=COSI1+EMCA(J,I)
         enddo
      enddo

      do k=1,kmaxp(10,1)
         wlix(10,1,k)=wlin(k)
         cinx(10,1,k)=weh(k)/(del(ik))
         ilabcfx(10,1,k)=1001
      enddo
      line_cool(10,1)=cosi1

C     *************************************************************
C     *****
C     CALCULATE POPULATIONS AND LINE EMISSION FROM FE I
C     *****
C     *************************************************************
      ZEL=DEL(IK)-ABN(14)*(XF(2)+2.*xf(3))

      IFPOP=0

      DO I=2,NFE
         DO J=1,I-1
            WOBS(I,J)=0.
            EMCA(I,J)=0.
         ENDDO
      ENDDO
      XFEI=1.-XF(2)-XF(3)-XF(4)
      TEFEi=TE
c!!
c     IF(TE.LT.500.) TEFE=500.
      IF(TE.LT.100.) TEFE=100.

c!!
      IF(IFEMUL.LT.0) THEN
        EPSFE=1.
      ELSE
        EPSFE=1.E-10
      ENDIF
      XQ=XF(1)+XF(2)

c      call POPFEi(RS,TE,ZEL,XQ,PHHE,IFPOP)

      COFEI=0.
      FEITOT=0.
      DO I=1,NFE
         DO J=1,NFE
            WOBFE(J,I)=WOBS(J,I)
            FEITOT=FEITOT+WOBS(J,I)/DEL(IK)
            COFEI=COFEI+EMCA(J,I)
         enddo
      enddo

      call popchianti_new(14,1,te,cofe1)      
      do k=1,kmaxp(14,1)
         wlix(14,1,k)=wlin(k)
         cinx(14,1,k)=weh(k)/(del(ik))
         ilabcfx(14,1,k)=1401
      enddo
      line_cool(14,1)=cofe1
      
C     *************************************************************
C     *****
C     CALCULATE POPULATIONS AND LINE EMISSION FROM FE II
C     *****
C     *************************************************************


c Fe II
      ilcf=1402
      z=abn(14)*xf(2)
      call popchianti_new(14,2,te,cofe2)

      do k=1,kmaxp(14,2)
         wlix(14,2,k)=wlin(k)
         cinx(14,2,k)=weh(k)/(del(ik))
         ilabcfx(14,2,k)=1402
      enddo
      line_cool(14,2)=cofe2

C     **************************************************************
C     *****
C     PHOTOIONIZATION HEATING
C     *****
C     **************************************************************
 8372 HYOLD=HY
      CALL EMISS(TE)
      CALL SPEC(TE)
      ZOLD=ZQ
      ZQ=DNEW
c      call rateaug(1,1)
c      call rateaug(2,2)
      call rateaug(6,6)
      call rateaug(7,7)
      call rateaug(8,8)
      call rateaug(10,10)
      call rateaug(11,11)
      call rateaug(12,12)
      call rateaug(14,14)
      call rateaug(16,16)
      call rateaug(18,7)
      call rateaug(20,3)
      call rateaug(26,26)

C     CALCULATE NON-THERMAL IONIZATION
      CALL NONTH(XEL,FHEAT,FIH,FEH,FIHE,FEHE)
c     TOTAL HEATING
      GATOT=0.
      DO K=1,100
         GATOT=GATOT+GEA(K)*AB(IAB(K))*XA(2,K)
      ENDDO
C     ADD NON-THERMAL IONIZATION
      ZKA(1)=ZKA(1)+FIH*GATOT/(13.6*1.602E-12)
      ZKA(2)=ZKA(2)+FIHE*GATOT/(25.4*1.602E-12)
      CSLYA=FEH*GATOT/(10.2*1.602E-12)

      DO I=1,100
         ZE(I)=ZEA(I)/DEN(IK)
         ZK(I)=ZKA(I)/DEN(IK)
         GE(I)=GEA(I)/DEN(IK)
      enddo
      
      MQ=MQ+1
C      IF(ZQ.NE.0.) DZ=ABS(ZQ-ZOLD)/ZQ
      MI=1
c!! mq.le.1
      IF(MQ.LE.1) GOTO 776
C     IF(DZ.GT.0.1) GOTO 776
C!!!  OTS FOR ALL IONS!!
cnots      GE=0.
cnots      PH=0.
c$$$C     HYDROGEN
c$$$      PHI=0.
c$$$      PHIT=0.
c$$$      DO J=1,5
c$$$         PHIJ=AB(1)*XNQ(5,J)*PHEAT(5,J)/DEN(IK)
c$$$C     CORRECT FOR NONTHERMAL IONIZATIONS AND EXCITATIONS
c$$$         IF(J.EQ.1) PHIJ=FHEAT*PHIJ
c$$$         PHITJ=AB(1)*XNQ(5,J)*PHEATT(5,J)/DEN(IK)
c$$$         PHI=PHI+PHIJ
c$$$         PHIT=PHIT+PHITJ
c$$$      ENDDO
c$$$      PHEI=FHEAT*AB(2)*XA(2,2)*GE(2)
c$$$  PHEII=FHEAT*AB(2)*XA(2,3)*GE(3)
      pheii=0.
      PO=0.
      POT=0.
      DO J=1,9
         POJ=FHEAT*AB(3)*XNQ(2,J)*PHEAT(2,J)/DEN(IK)
         POTJ=FHEAT*AB(3)*XNQ(2,J)*PHEATT(2,J)/DEN(IK)
         PO=PO+POJ
         POT=POT+POTJ
      ENDDO
      poxt=0.
      DO L=1,4
         POX(L)=FHEAT*AB(3)*XA(2,L+13)*GE(L+13)
         poxt=poxt + pox(l)
      enddo
      POX(5)=FHEAT*AB(3)*XA(2,7)*GE(7)
      poxt=poxt + pox(5)
      DO L=6,8
         POX(L)=FHEAT*AB(3)*XA(2,L-2)*GE(L-2)
         poxt=poxt + pox(l)
      enddo
      
C     OTS
      IF(OTSP(4).LE.1) POX(5)=0.
      DO L=1,4
         PC(L)=FHEAT*AB(4)*XA(2,L+7)*GE(L+7)
      enddo
      DO L=1,2
         PC(L+4)=FHEAT*AB(4)*XA(2,L+11)*GE(L+11)
      enddo
C     OTS
      IF(OTSP(5).LE.1) PC(5)=0.
      DO L=18,24
         PN(L-17)=FHEAT*AB(5)*XA(2,L)*GE(L)
      enddo
C     OTS
      IF(OTSP(6).LE.1) PN(1)=0.      
      pne=0.
      do l=1,10
         pne=pne + fheat*abn(6)*xa(2,71+l)*ge(71+l)
      enddo
      PMG=0.
      DO L=1,3
         PMG=PMG+FHEAT*AB(7)*XA(2,L+38)*GE(L+38)
      enddo
      PSI=0.
      DO L=1,14
         PSI=PSI+FHEAT*abn(10)*XA(2,24+L)*GE(24+L)
      enddo
      PSU=0.
      DO L=1,16
         PSU=PSU+FHEAT*abn(11)*XA(2,55+L)*GE(55+L)
      enddo
      par=0.
      do l=1,2
         par=par + fheat*abn(12)*xa(2,81+l)*ge(81+l)
      enddo
      do l=3,7
         par=par + fheat*abn(12)*xa(2,93+l)*ge(93+l)
      enddo
      pca=0.
      do l=1,2
         pca=pca + fheat*abn(13)*xa(2,49+l)*ge(49+l)
      enddo
      PFE=0.
      DO L=1,4
         PF(L)=FHEAT*ABN(14)*XF(L)*GE(L+41)
         PFE=PF(L)+PFE
      enddo
      DO L=5,15
         PF(L)=FHEAT*ABN(14)*XF(L)*GE(L+79)         
         PFE=PF(L)+PFE
      enddo
      DO L=16,26
         PF(L)=FHEAT*ABN(14)*XF(L)*GE(L+85)
         PFE=PF(L)+PFE
      enddo

C     **************************************************************
C     *****
C     TOTAL COOLING AND HEATING RATES
C     *****
C     **************************************************************
      CEXTOT=0.
      DO K=1,9
         CEXTOT=CEXTOT+CEX(K)
      enddo
      cextot=0.
c      HE=HEI+CHEII
c      line_cool(2,1)=hei
c      line_cool(2,2)=cheii
      COCAR=C211+C212+c213+C221+C222+C231+C232+C241
      CONIT=C321+C331+C332+C333+C341+C342+C351+C311+C322
      COOX=C161+C152+C151+C131+C142+C143+C111+C113+C121+C132
      COONE=C1322+C1371+C1351+C1341+CONEIII
      CONAT=C1011
      COMAG=C521+C511+C512
      COSIL=C411+C412+C421+C422+C431+C441
      COOSU=C1211+cos2+cos3+cos4+cos5+C1261
      COArg=C1422 + coariii + coariv
      COOFE=CIN(18)+COfe2+cofe3+cofe4+cofe5+cofe6+cofe7+cofe8+cofe9
      coolt=0.
c no H, He      
      do iel=3,14
         do ion=1,27
            coolt=coolt+line_cool(iel,ion)
            if(abs(line_cool(iel,ion)) > 0.) then
            endif
         enddo
      enddo

      COOL=FF+RECCA+RECOX+RECSI+RECS+RECFE + coolt
     
c     Photelectric heating
      PCAR=0.0E0
      DO L=1,6
         PCAR=PCAR+PC(L)
      enddo
      PNI=0.0E0
      DO L=1,7
         PNI=PNI+PN(L)
      enddo
      POXY=PO
c     obs! O I SEPARATELY
      DO L=2,8
         POXY=POXY+POX(L)
      enddo
c!!! Not Na ! pn() = N
      PNA=0.
c      DO L=1,7
c         PNA=PNA+PN(L)
c      enddo
      HEAT=CO+PCAR+PNI+POXY+pne+PNA+PMG+PSI+psu+par+pca+pfe
C     ADD GAMMA RAY HEATING
C
      HEAT=HEAT+GAELH/DEN(IK)**2
      HY=0.
      RAD=(HEAT-(HY+COOL)*DEL(IK))
9323  format(5e12.4)
C     APPROX. REC. EMISSION FROM C, O, S, SI, FE
      RECAPPC=DEL(IK)*RECCA*1.307e5/(1.5*TE)
      RECAPPO=DEL(IK)*RECOX*1.581e5/(1.5*TE)
      RECAPPSI=DEL(IK)*RECSI*0.946e5/(1.5*TE)
      RECAPPS=DEL(IK)*RECS*1.194e5/(1.5*TE)
      RECAPPFE=DEL(IK)*RECFE*0.9136e5/(1.5*TE)
C     RECOMBINATION COOLING
      CALL XIONS(X,CHI)
      DO K=1,100
            IF(algs(k).gt.0.) THEN
                  RECEX(K)=X(K)*CHI(K)*ELCH*(ALTOT(K)-ALGS(K))*DEL(IK)
                  RECGS(K)=X(K)*CHI(K)*ELCH*ALGS(K)*DEL(IK)
            ENDIF
      ENDDO

      write(6,9311)ff,cocar,conit,coox,coone,conat,comag,
     &     cosil,coosu,coarg,cocal,coofe,cool
c      write(6,*)'Fractional cooling '
c      write(6,9311)ff/cool,cocar/cool,
c     &     conit/cool,coox/cool,coone/cool,conat/cool,comag/cool,
c     &     cosil/cool,coosu/cool,coarg/cool,cocal/cool,coofe/cool,cool
 9311 format('F-F ',1Pe11.3,' C',e11.3,
     &     ' N',e11.3,' O',e11.3,' Ne',e11.3,' Na',e11.3,' Mg',e11.3,
     &     ' Si',e11.3,' S',e11.3,' Ar',e11.3,' Ca',e11.3,' Fe',e11.3,' Total',e11.3)

      XELEC=DEL(IK)
c      write(6,*)'Fractional heating '
c      write(6,9313)pcar/heat,pni/heat,
c     &     poxy/heat,pne/heat,pna/heat,pmg/heat,psi/heat,psu/heat,
c     &     par/heat,pca/heat,pfe/heat,heat/heat      
c      write(6,*)' '
 9313 format(' C',1pe11.3,
     &     ' N',e11.3,' O',e11.3,' Ne',e11.3,' Na',e11.3,' Mg',e11.3,
     &     ' Si',e11.3,' S',e11.3,' Ar',e11.3,' Ca',e11.3,
     &     ' Fe',e11.3,' Total',e11.3)
      WRITE(6,7362)ik,TE,DEL(IK),r(ik),den(ik),HEAT,COOL,del(ik)*cool,
     &     gahe/(DEL(IK)*den(ik)**2),RAD

 8278 IF(IFPOP.EQ.1) WRITE(6,*)' NO CONVERGENCE IN POP-ROUTINE'
 7362 FORMAT(1X,'I,T,X,R,DEN,H,C,X*C,H-XC ',i5,1pe12.3,e12.3,e15.7,8E11.4)
      inition=0
c no H, He
      do iel=3,14
         do ion=1,27
            xion_old(iel,ion)=xion(ik,iel,ion)
         enddo
      enddo

      RETURN
       END



      double precision function fir(n,z,te)
C     CALCULATE GOULDS FUNCTION PSI IN GOULD EQ (31) 
      IMPLICIT REAL*8(A-H,O-Z)      
            y=z**2*15.79e4/te
            y4=z**2*15.79
            fir=fin(n,y)/fin(n,y4)
      return
      end

      double precision function fin(N,Y)
      IMPLICIT REAL*8(A-H,O-Z)      
      IF(Y.LT.0.10) THEN
            FIN=Y*(-1.202*LOG(Y)-.298)
      ELSE
C            FIN=0.4288+0.5*LOG(Y)+0.469/(Y**.33333)
            fiN=0.5*(log(y)+1.735+1./(6.*y))
      ENDIF
      DO M=1,N-1
            X=Y/M**2
      IF(X.LT.1.) THEN
            E1=-LOG(X)-.577216+X-.24991*X*X+.0552*X*X*X-.00976
     &            *X*X*X*X
            E1=E1*EXP(X)
      ELSE
            E1=(X*X+2.334733*X+.250621)/((X*X+3.330657*X
     &            +1.681534)*X)
      ENDIF
            FIM=X*E1/M
            FIN=FIN-FIM
      ENDDO
      RETURN
      END


      subroutine feiiiatom(te,den,z,rtot)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NMAX=30)
      SAVE
      character*80 lab
      character*11 lev
      common/wfeiv/wl(45,45),fbfeiii(41,41)
      COMMON/FE3/WLFE3(30,30),AFE3(30,30),EMFE3(30,30)
      DIMENSION FB(30),RL(41,41),X(NMAX)
      dimension e(45),g(45),a(45,45),c(45,45)
      if(initfe3.eq.0) then
         initfe3=1
c         REWIND 28
         read(28,1921)lab
         read(28,1921)lab
 1921    format(a)
         n=nmax
         do  i=1,n
            read(28,912)ix,e(i),g(i),lev
         enddo
 912     format(i3,f12.2,f7.1,a11)
         fe4=0.
         read(28,1921)lab
         read(28,1921)lab
         do i=2,n
            im1=i-1
            do  j=1,im1
               read(28,*)iq,jq,dlair,a(i,j),c(i,j)
            enddo
         enddo
         close(28)
      endif
      t4=te/1.e4
      CALL FORBN(NMAX,E,G,A,C,TE,Z,DEN,RL,X)
      rtot=0.
      do  i=2,nmax
         do  j=1,i-1
            wlfe3(i,j)=wl(i,j)
            emfe3(i,j)=z*rl(i,j)
            rtot=rtot+rl(i,j)
         enddo
      enddo
      return
      END

      subroutine feivatom(te,den,z,rtot)
c!!cfq hela rutinen och naesta!
C
C       THIS PROGRAM CALCULATES THE EMISSIVITY OF FORBIDDEN LINES
C       USING A 10-LEVEL SCHEME FOR CIII
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NMAX=22)
      SAVE
      character*80 lab
      character*11 lev
      common/wfeiv/wl(45,45),fbfeiv(41,41)
      COMMON/FE4/WLFE4(30,30),AFE4(30,30),EMFE4(30,30)
      DIMENSION FB(30),RL(41,41),X(NMAX)
      dimension e(45),g(45),a(45,45),c(45,45)
      if(initfe.eq.0) then
         initfe=1
c         REWIND 27
         read(27,1921)lab
         read(27,1921)lab
 1921    format(a)
         n=nmax
         do i=1,n
            read(27,912)ix,e(i),g(i),lev
         enddo
 912     format(i3,f12.2,f7.1,a11)
         fe4=0.
         read(27,1921)lab
         read(27,1921)lab
         do i=2,n
            im1=i-1
            do j=1,im1
               read(27,*)iq,jq,dlair,a(i,j),c(i,j)
            enddo
         enddo
         close(27)
      endif
      t4=te/1.e4
c!!
      CALL FORBN(NMAX,E,G,A,C,TE,Z,DEN,RL,X)
      rtot=0.
c!!
      do i=2,nmax
         do j=1,i-1
            wlfe4(i,j)=wl(i,j)
            emfe4(i,j)=z*rl(i,j)
            rtot=rtot+rl(i,j)
         enddo
      enddo
      
      return
      END

      SUBROUTINE FORBN(N,E,G,A,CQ,TE,Z,DEN,RL,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N
c      PARAMETER (MD=300,MDP1=MD+1)
      include "parameters.h"
      COMMON/PHY/DENS(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/IND/IK
      common/wfeiv/wl(45,45),fbfeiv(41,41)
      DIMENSION  C(45,45),CQ(45,45),G(45),E(45),EEV(45),A(45,45)
      DIMENSION AA(NLP1,NLP1),X(N),RL(41,41),BE(45,45)
      DATA PI/3.1415926E0/
      DEN=DENS(IK)
      DENEL=DEN*DEL(IK)
      NM1=N-1
      NP1=N+1
      DO  I=1,NP1
         DO  J=1,NP1
            AA(I,J)=0.
         enddo
      enddo
      DO I=1,N
         EEV(I)=E(I)/8065.46
      enddo
      DO 5395 I=1,N
         DO 5394 J=1,N
            WL(I,J)=0.0
            be(I,J)=0.0
            IF(I.EQ.J) GOTO 5394
            IF(E(I).EQ.E(J)) GOTO 5394
            WL(I,J)=DABS(12398.54/(EEV(I)-EEV(J)))
5394  CONTINUE
5395  CONTINUE
      TEV=TE/1.1609D4
      TIME=8.64E4*TDAY
      DO I=1,NM1
         IP1=I+1
         DO J=IP1,N
            C(J,I)=8.63D-6*DEN*CQ(J,I)/(DSQRT(TE)*G(J))
            C(I,J)=C(J,I)*DEXP(-(EEV(J)-EEV(I))/TEV)*G(J)/G(I)
         enddo
      enddo
      DO L=1,4
         DO I=1,N
            DEI=X(I)*Z*DEN
            DO J=1,N
               IF(I.NE.J) then
C
C     CALCULATE THE ESCAPE PROBABILITY.
C
                  T0=1.e-24*WL(J,I)**3*A(J,I)*G(J)*DEI*TIME/(8.*PI*G(I))
                  IF(L.EQ.1) T0=0.0
                  IF(T0.LT.0.01) BE(J,I)=1.-0.5*T0+T0**2/6.
                  IF(T0.GE.0.01) BE(J,I)=(1.-EXP(-T0))/T0
                  AA(I,J)=BE(J,I)*A(J,I)+C(J,I)
               endif
            enddo
         enddo
         
         do i=1,n-1
            S=0.
            DO K=1,N
               S=S+C(I,K)+BE(I,K)*A(I,K)
            enddo       
            AA(I,I)=-S
         enddo

         DO J=1,N
            AA(n,J)=1.D0
         enddo

         NRC=N+1
         AA(N,NRC)=1.
         EPS=1.D-30
         NA=N

         DI=SIMUL(NA,AA,X,EPS,1,NRC)

      ENDDO
      RTOT=0.
      DO I=1,NM1
         IP1=I+1
         DO  J=IP1,N
            IF(A(J,I).GT.0.) THEN
               RL(J,I)=1.602E-12*(EEV(J)-EEV(I))*X(J)*BE(J,I)*A(J,I)/DEN
               RTOT=RTOT+RL(J,I)
            ENDIF
         enddo
      enddo

25    CONTINUE
      RETURN
      END

      subroutine  heiiatom(z,t,den,abhe,xe,x0,xp,t0)
c     He II atom including Ly alpha and 2-gamma emission
c     also applies to H I with same levels
c     910420
c
c     t = temp
c     den = total density
c     abhe = He abundance
c     xe = el. frac.
c     x0 = n = 1 frac
c     xp = he III abundance
c     tday = time in days
c
c     results are the emissivities (incl. the 4 pi division)
c     jrec(k),k=1,4  = recomb. cont. corrected for thermal contr.
c     jla = Lyman alpha 
c     j2g = two-gamma
c     jem(l,k) = line emission 
c
      implicit real*8(a-h,o-z)
      real*8 n2,n1,ne,nhe3,jla,j2g,jrec,jem,jlahe,j2ghe,jreche,jemhe
      real*8 jlah,j2gh,jrech,jemh
c     PARAMETER(NFEL=1000)
      include "parameters.h"
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/COLD/OPT(nio),COLTOT(nio),COLTOTT(nio),SRED(NFEL+100)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/A7/C3,C33
      common/heiires/x2he,coll12he,jlahe,j2ghe,jreche(4),jemhe(7,4)
      common/hres/x2h,coll12h,jlah,j2gh,jrech(4),jemh(7,4)
      dimension aln(4),jem(7,4),jrec(4)
      data pi/3.14159/
      tev=t/1.1609e4
      t4=t/1.e4
      iz=nint(z)
      n1=x0*abhe*den
      ne=xe*den
      nhe3=xp*abhe*den
      wl=z**2*1216.*1.e-8
      g2=6.
      g1=2.
      e21=z**2*0.75*13.397
      hnu=1.602e-12*e21
      theii=t/z**2
      theii4=theii/1.e4
c     effective recomb. to n=2 for Case B from Pengelly (MN 127, 145)
      al2eff=z*(0.837+1.76)*1.e-13/theii4**0.806
      a21=z**4*6.25e8
      a2g=z**4*8.23/4.
c     from daresbury comp.
      if(iz.eq.1) omega=0.259*t4**.073+0.419*t4**.215
      if(iz.eq.2) omega=0.142/t4**.059+0.319*t4**.166
      c21=8.63e-6*omega/(sqrt(t)*g2)
      c12=g2*c21*exp(-e21/tev)/g1
      IF(ISTAT.EQ.1) then
       IF(IZ.EQ.1) IJ=1
       IF(IZ.EQ.2) IJ=3
       T0=C33*COLTOT(IJ)*A21*G2*WL**3./G1
       T0T=(COLTOTT(IJ)-COLTOT(IJ))*T0/COLTOTT(IJ)
       IF(T0T.LE.0.) T0T = 0.
C
C      CALCULATE THE ESCAPE PROBABILITY.
C
       BE=1.
       AU=4.
       VTERM=1.285E6*SQRT(T4/(AU))
       CALL ESCAPE(T0,T0T,WL,A21,VTERM,BE,DBEDTA)
       a21beta=a21*be
      else
       time=8.64e4*tday
       a21beta=8.*pi*g1/(3.*wl**3*n1*g2*time)
       t0=(3.*wl**3*a21*g2*n1*time)/(8.*pi*g1)
      endif
      n2=ne*(al2eff*nhe3+g2*exp(-e21/tev)*c21*n1)/
     &      (a2g+a21beta+c21*ne)
      x2=n2/(n1+n2+nhe3)
      jla=0.75*n2*hnu*a21beta/(4.*pi*ne*den)
      j2g=n2*hnu*a2g/(4.*pi*ne*den)
      coll12=n1*ne*hnu*c12/(4.*pi*ne*den)
c     rec cont emission
      call hrec(theii,aln,ala,alb)
      do k=1,4
         enk=z**2*13.597/k**2
         alheii=z*aln(k)
         jrec(k)=nhe3*ne*alheii*enk*1.602e-12/(4.*pi*ne*den)
c     correct for the thermal energy of the electrons (from Osterbrock)
         if(k.eq.1) thermcorr=1.+0.949*t/(enk*1.1609e4)
         if(k.gt.1) thermcorr=1.+0.67*t/(enk*1.1609e4)
         jrec(k)=thermcorr*jrec(k)

      enddo
c rec line emission
      call hrecem(theii,jem)
      do k=2,4
         do l=k+1,5
            jem(l,k)=nhe3*ne*z**3*jem(l,k)/(4.*pi*ne*den)
         enddo
      enddo
      totrec=z**2*13.597*1.602e-12*z*ala*nhe3*ne/(4.*pi*ne*den)
      if(iz.eq.1) then
        x2h=x2
        coll12h=coll12
        jlah=jla
        j2gh=j2g
	do k=1,4
           jrech(k)=jrec(k)
           do l=1,7
              jemh(l,k)=jem(l,k)
           enddo
	enddo
      elseif(iz.eq.2) then
        x2he=x2
        coll12he=coll12
        jlahe=jla
        j2ghe=j2g
	do k=1,4
           jreche(k)=jrec(k)
           do l=1,7
              jemhe(l,k)=jem(l,k)
           enddo
	enddo
      endif
      return
      end

      subroutine hrec(t,aln,ala,alb)
c
c     Hydrogen recombination rates taken from Osterbrock at 1E4 K
c     powerlaw fit between 5000 and 10000 K
c
      implicit real*8(a-h,o-z)
      dimension all(4,4),aln(4)
      t4=t/1.e4
      do i=1,4
         do j=1,4
            all(i,j)=0.
         enddo
      enddo
      all(1,1)=1.58e-13/t4**.529
      all(2,1)=2.34e-14/t4**.365
      all(2,2)=5.35e-14/t4**.639
      all(3,1)=7.81e-15/t4**.532
      all(3,2)=2.04e-14/t4**.643
      all(3,3)=1.73e-14/t4**.809
      all(4,1)=3.59e-15/t4**.543
      all(4,2)=9.66e-15/t4**.644
      all(4,3)=1.08e-14/t4**.815
      all(4,4)=5.54e-15/t4**.976
      do k=1,4
         aln(k)=0.
         do l=1,4
            aln(k)=all(k,l)+aln(k)
         enddo
      enddo
      ala=4.18e-13/t4**.706
      alb=2.59e-13/t4**.810
      return
      end

      subroutine hrecem(t,jem)
c
c     Hydrogen recombination emission rates taken from Osterbrock at 1E4 K
c     Table 4.4 T=1.e4 n=1.e6
c     case B
c     powerlaw fit between 5000 and 10000 K.
c
      implicit real*8(a-h,o-z)
      real*8 jem,jhb
      dimension jem(7,4)
      t4=t/1.e4
      jem(3,2)=2.81
      jem(4,2)=1.
      jem(5,2)=0.471
      jem(6,2)=0.262
      jem(7,2)=0.163
      jem(4,3)=0.317
      jem(5,3)=0.335
      jem(6,3)=0.339
      jem(7,3)=0.333
      jem(5,4)=0.154
      jem(6,4)=0.163
      jem(7,4)=0.163
      jhb=1.25e-25/t4**0.828
      do k=2,4
      do l=k+1,5
      jem(l,k)=jhb*jem(l,k)
      enddo
      enddo
      return
      end

      SUBROUTINE ATDATar5
      IMPLICIT REAL*8(A-H,O-Z)
      character*1 dum
      COMMON/datar5/Gar(5),Ear(5),aar(5,5),WLar(5,5)
      common/collar5/omt(12),omfitar(20,20,15),tefitar(15)
      n=5
      rewind (45)

      do i=1,n
         gar(i)=0.
         ear(i)=0.
         do j=1,n
            aar(j,i)=0.
         enddo
      enddo
      open(45,file='./ATDAT/ar5_lev_coll_rad.dat')
      read(45,987)  dum
      read(45,987)  dum
987   format(a)
      do i=1,n
         read(45,*)iq,jlev,wn
         gar(i)=2.*jlev+1.
         ear(i)=wn/8065.46d0
      enddo

      do i=1,n
         DO J=i+1,N
            wlar(j,i)=12398.54/(ear(j)-ear(i))
            wlar(i,j)=wlar(j,i)
         enddo
      enddo
      do i=1,1000
         read(45,*,err=11,end=11)il,iu,wlq,el,eu,gl,gu,a21
         aar(iu,il)=a21
      enddo

 11   continue
      read(45,987)  dum
      nte=11
      read(45,*)(tefitar(k),k=1,nte)
      do i=1,1000
         read(45,*,err=12,end=12)iu,il,(omt(k),k=1,nte)
         do k=1,nte
            if(iu.eq.4.or.iu.eq.5.and.il.ne.4) then
c  fine structure splitting for 4th and 5th levels
               omfitar(iu,1,k)=1.*omt(k)/9.
               omfitar(iu,2,k)=3.*omt(k)/9.
               omfitar(iu,3,k)=5.*omt(k)/9.
            else
               omfitar(iu,il,k)=omt(k)
            endif
         enddo
      enddo
 12   continue
      close(45)

      return
      end      

      subroutine read_coll
      implicit real*8(a-h,o-z)
      common/ca4/teca4(50),coll_ca4(3,50)
      open(11,file='./ATDAT/Ca_IV_coll_strengths_Nahar_v2.dat',status='old')
      do i=1,50
         read(11,*)ii,teca4(i),(coll_ca4(k,i),k=1,3)
      enddo
      close(11)
      return
      end
      
