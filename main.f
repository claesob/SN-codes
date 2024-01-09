C     CIRCUMSTELLAR:
C     SET ICS = 1
C         INOUT = 0
C     SPECIFY RSHOCK, MTOT
C     FOR REVERSE SHOCK 
C     ICS = 1
C     INOUT=1
C     IREV=1
C     
C***********************************************************************
C**********
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MTOT,MNI56,MSUN
      character*72 text1,text2,text3,TEXT4,INMOD,MODNRSTR,intext(1000)
      character*6 FILE1
      character*6 FILE2
      character*7 FILE3
      character*6 FILE4      
      character*9 FILE5
      character*8 FILE6
      character*9 FILE7
      character*20 FILE8
      character*15 infil
      CHARACTER*25 LABEL
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/MOD/INMOD
      COMMON/TEXT/TEXT1,TEXT2,TEXT3,TEXT4
      COMMON/ITER/ITE
      COMMON/A5/TAUE,ALFA2,EN,TSN,XL40,TEXA,RMAX
      COMMON/PULSAR/RLTOT,ALFA,ebreak,EMIN,EMAX,teff_ns,nsterm,ibreak
      COMMON/DENS/DEN0,R0,RN
      COMMON/MPAR/MTOT,MNI56,VEXP
      COMMON/RQW/TEFF,RQ
      COMMON/INUT/IUY
      COMMON/ELDEN/IDENS
      COMMON/NBA/NBACK
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/FRE/NINQ,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/AUG/AUG
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/INSPEC/IBLACKB,IPULSSP,IFERNET,INLR,IFERMAT,ICRAB
      COMMON/AGNPARA/DECON,GAMMION,IPRESS,IDEN,ICDEN
      COMMON/MULTI/IHMUL,IHEMUL,IOMUL,ICAMUL,IFEMUL
      COMMON/UVBL/IUVBLANK
      COMMON/SPH/ISPH
      COMMON/GAMPAR/GAMLUM,GAELH,GAHE
      COMMON/HYDROS/HSCALE,DEN1,XEL,AMEAN
      COMMON/CHG/CHGH,CHGHE,CHTE
      COMMON/THER/A1,B1,TIN,E10,E20
      COMMON/ITSPH/NITSPH
      COMMON/ABUN/AB(20)
      COMMON/A12/ITER
      COMMON/A10/AV
      COMMON/IOPAR/IOSHELL,IORAD1,IORAD2,IOFLUX
      COMMON/CICS/RSHOCK,ICS
      common/revs/terev,tecs
      COMMON/CINOUT/INOUT,IPULS
      COMMON/LOWION/ILOWION
      COMMON/FELEV/NFEII,nfei
      COMMON/IOINF/IOSPEC,IOLEVEL,IOEMISS,IOTERM,IOION,IOELITER
      COMMON/CHEMCOMP/ABIN(12),ISOLAR,ISNCOMP,INPUTCOMP
      COMMON/THBAL/TMIN,TMAX,IBAL,IKBAL(1000)
      integer avabund,zams,zone
      common/abzone/rm_abun,avabund,zams,zone
      integer nz,nion,nshell,i,k
      parameter(nz=30,nion=27,nshell=10)
      integer ns,kmax
      real*8 fr_aug,eion,en_aug,en_augi,eioni,min_vel
      integer nstemp
      common/ns_spec/tot_lum,FL_ns(NE1:NE2),nstemp
      common/augerfrac/eion(nz,nion,nshell),en_aug(nz,nion,nshell),
     &     fr_aug(nz,nion,nshell,10),kmax(nz,nion,nshell),ns(nz,nion),
     &     init_augfrac
      common/initchianti/initch
      common/inidiel/initdi
      common/initox/initoi
      common/initat/initfeii,initfei,initsi1
      common/init_ion/inition      
      character*100 dum(1000)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DATA SOLMA/1.989E33/

      INOUT=2
      IPULS=1
      call flush()
      
C
C     INITIALIZE ALL PARAMETERS
C
C
C	OPEN FILES
C
c      OPEN(11,FILE='./ATDAT/cordat.dat',status='old')
      open(27,file='./ATDAT/feIV_22.dat',status='old')
      open(28,file='./ATDAT/feIII_30.dat',status='old')
      
C     INPUT PARAMETERS IN THIS FILE
      WRITE(0,*)' INFIL?'
      READ(5,9899)INFIL(8:15)
      INFIL(1:7)='pulin3.'
      write(0,*)infil
      OPEN(99,FILE=infil)
      do i=1,1000
         read(99,9281,end=1234)intext(i)
         write(6,9281)intext(i)
 9281    format(a72)
      enddo
 1234 continue
      inmax=i-1
      rewind(99)
      READ(99,9899)TEXT1
      READ(99,9899)TEXT2
      READ(99,9899)TEXT4
      READ(99,9899)TEXT3      
      READ(99,*)MODNR
      MN1=IFIX(real(modnr)/100)
      MN2=MODNR-100*MN1
      MN2=IFIX(real(MN2)/10)
      MN3=MODNR-MN1*100-MN2*10
      MN3=IFIX(real(MN3))

      MODNRSTR(1:1)=CHAR(MN1+48)
      MODNRSTR(2:2)=CHAR(MN2+48)
      MODNRSTR(3:3)=CHAR(MN3+48)

C     FILE1= GENERAL OUTPUT FILE = UNIT 6
      FILE1(1:3)='slt'
      FILE1(4:6)=MODNRSTR(1:3)
      
      OPEN(6,FILE=FILE1)
      rewind(99)
      do i=1,inmax
         write(6,9181)intext(i)
         write(116,9181)intext(i)
 9181    format(a72)
      enddo
      write(6,*)'Input file for model ',MODNRSTR(1:3)
      write(6,*)' '
c      write(6,*)'rm_abun,avabund,zams,zone ',rm_abun,avabund,zams,zone
      write(6,9378)rm_abun
 9378 format('Mass for abundances, rm (1.76, 1.85, 1.95, 2.0): ',1pe12.3)
      do i=1,1000
         read(99,932,end=323)dum(i)
 932     format(a100)
      enddo
 323  lmax=i-1
      rewind(99)

      READ(99,9899)TEXT1
      READ(99,9899)TEXT2
      READ(99,9899)TEXT4
      READ(99,9899)TEXT3      
      READ(99,*)MODNR
c      write(6,*)' Model ',modnr
      FILE8(1:13)='strong_lines_'
      FILE8(14:16)=MODNRSTR(1:3)
      FILE8(17:20)='.txt'
      OPEN(76,FILE=FILE8)
9899  FORMAT(A)
      FILE4(1:3)='slf'
      FILE4(4:6)=MODNRSTR(1:3)
      FILE6(1:5)='struc'
      FILE6(6:8)=MODNRSTR(1:3)
      OPEN(59,FILE=FILE6)
      JMIN=-199
      JJ=129
      NINQ=2
      CALL ENINT
      init_augfrac=1
      initch=1
      initoi=1
      initdi=1
      inition=1
      call auger_fr
      call atdatar5
c Auger data      
      call readferec
c data for recomb. of Si, S and Ar      
      call recomb_adas
c call recomb data from Badnell et al for Na, Mg, Al, Si, P sequencies      
      call badnell_et_al
c      call badnell_dr_fits      
C     DEFAULT PARAMETERS 
      MAX=299
      ICEN=35
      MNI56=1.E-20
      A1IN=3000.
      B1IN=10000.
      NFEII=116
      tday=12927.
      TOTCOLUMN=1.E33
      RIN=1.E14/1.E15
      FILLING=1.
      CGGH=0.1
      CHTE=0.2
      ISTEPION=1
      VEXP=1.00E9
      IPH=0
      IOSHELL=0
      IORAD1=0
      IORAD2=0
      IOFLUX=0
      IOSPEC=0
      IOLEVEL=0
      IOTERM=0
      IOION=0
      IOELITER=0
      IDEP=0
      IHMUL=1
      IHEMUL=1
      IOMUL=1
      ICAMUL=1
      IFEMUL=1
      iagn=0
      nmax=1
      TOTCOLUMN=1.d33
      MNI56=0.07
cINITIAL FRACTION OF STROEMGREN DEPTH FOR FIRST SHELL      
      DMINIT = 0.0000015                                                                
C     IF IDENS < 0 USE A CONSTANT ELECTRON FRACTION IN POPULATION CALC.
      IDENS=1
C     IF NBACK = 0 PUT BACKGROUND FLUX IN MULTI LEVEL CALC. TO ZERO.
      NBACK=1
C     IF AUG>0 INCLUDE AUGER IONIZATIONS
      AUG=999.
C     QSO=INDEX SPECIFYING IF THE CALCULATION APPLIES TO A QUASAR
C         OR NOT. QSO>0 : NO BACKGROUND RADIATION FIELD, QSO<0 :
       QSO=-999.
      DO K=1,100
        READ(99,9)LABEL
 9      FORMAT(A20)
        IF(LABEL(1:7).EQ.'INMODEL') THEN
          READ(99,9899)INMOD
        ELSEIF(LABEL(1:5).EQ.'INOUT') THEN
C         IF INOUT = 1 START FROM CENTER, 0 from outer boundary
          read(99,*)INOUT
        ELSEIF(LABEL(1:7).EQ.'TEMPLIM') THEN
          READ(99,*)A1IN,B1IN
        ELSEIF(LABEL(1:7).EQ.'BALANCE') THEN
C         SCAN THE HEATING AND COOLING CURVES FROM TMIN TO TMAX FOR
C         ALL I = IBAL(K)
          READ(99,*)TMIN,TMAX
          IBAL=1
          DO KQ=1,1000
            READ(99,*)IKBAL(KQ)
            IF(IKBAL(KQ).EQ.0) GOTO 569
          ENDDO
569       CONTINUE
        ELSEIF(LABEL(1:4).EQ.'MTOT') THEN
C         MTOT = TOTAL MASS IN G
          READ(99,*)MTOT
        ELSEIF(LABEL(1:5).EQ.'RIN15') THEN
C         RIN= RADIUS OF INNER BOUNDARY IN 1E15 CM
          READ(99,*)RIN
          IRINNER=1
        ELSEIF(LABEL(1:7).EQ.'FILLING') THEN
C         FILLING FACTOR
           READ(99,*)FILLING
        ELSEIF(LABEL(1:7).EQ.'MIN_VEL') THEN
C         veolcity (in km/s) at inner boundar to set radius at time t
           READ(99,*)min_vel
        ELSEIF(LABEL(1:9).EQ.'DIFF ITER') THEN
C         NMAX=NUMBER OF ITERATIONS FOR RADIATION FIELD. EG. NMAX=1 MEANS
C          NO ITERATION , NMAX=2 , ONE ITERATION ETC.
          READ(99,*)NMAX
        ELSEIF(LABEL(1:6).EQ.'SHELLS') THEN
C         MAX=TOTAL NUMBER OF RADIAL SHELLS AND NUMBER OF SHELL IN CENTRE
          READ(99,*)MAX
        ELSEIF(LABEL(1:7).EQ.'CSHELLS') THEN
C         ICEN=NUMBER OF SHELLS IN CENTRE
          READ(99,*)ICEN
        ELSEIF(LABEL(1:5).EQ.'RMCEN') THEN
C         MAX=TOTAL NUMBER OF RADIAL SHELLS AND NUMBER OF SHELL IN CENTRE
          READ(99,*)RMCEN
        ELSEIF(LABEL(1:7).EQ.'RMINNER') THEN
          READ(99,*)rminner
        ELSEIF(LABEL(1:9).EQ.'UNIF STEP') THEN
          ISTEPION=0
        ELSEIF(LABEL(1:9).EQ.'CH LOGOPT') THEN
C         CHGH=MAXIMUM LOGARITMIC CHANGE IN OPTICAL DEPTH AT J=2
          READ(99,*)CHGH
        ELSEIF(LABEL(1:7).EQ.'CH LOGT') THEN
C         CHTE=RATIO BETWEEN ALLOWED CHANGE IN THE LOG(TEMPERATURE) TO THE
          READ(99,*)CHTE
        ELSEIF(LABEL(1:11).EQ.'ENERGY BINS') THEN
C         JMIN=INDEX OF MINIMUM ENERGY BIN
C         JJ= D:O OF MAXIMUM BIN
          READ(99,*)ENMIN,ENMAX
        ELSEIF(LABEL(1:5).EQ.'FE II') THEN
C         NUMBER OF FE II LEVELS
          READ(99,*)NFEII
        ELSEIF(LABEL(1:19).EQ.'ITERATIVE MULTI H I') THEN
C         SIMPLE MULTILEVEL CALCN. FOR H I
          IHMUL=0
        ELSEIF(LABEL(1:20).EQ.'ITERATIVE MULTI HE I') THEN
C         SIMPLE MULTILEVEL CALCN. FOR HE I
          IHEMUL=0
        ELSEIF(LABEL(1:19).EQ.'ITERATIVE MULTI O I') THEN
C         SIMPLE MULTILEVEL CALCN. FOR O I
          IOMUL=0
        ELSEIF(LABEL(1:21).EQ.'ITERATIVE MULTI CA II') THEN
C         SIMPLE MULTILEVEL CALCN. FOR CA II
          ICAMUL=0
        ELSEIF(LABEL(1:21).EQ.'ITERATIVE MULTI FE II') THEN
C         SIMPLE MULTILEVEL CALCN. FOR FE II
          IFEMUL=0
        ELSEIF(LABEL(1:18).EQ.'NO MULTI FOR FE II') THEN
C         SIMPLE MULTILEVEL CALCN. FOR FE II
          IFEMUL=-1
        ELSEIF(LABEL(1:5).EQ.'DECON') THEN
C         AGN DENSITY AT INNER BOUNDARY
          READ(99,*)DECON
          ICDEN=1
        ELSEIF(LABEL(1:6).EQ.'IONPAR') THEN
C         AGN IONIZATION PARAMETER
          READ(99,*)GAMMION
        ELSEIF(LABEL(1:9).EQ.'BLACKBODY') THEN
C         BLACK-BODY SPECTRUM
          READ(99,*)TEFFBB
          IBLACKB=1
        ELSEIF(LABEL(1:19).EQ.'PHOTOSPHERIC RADIUS') THEN
C         PHOTOSPHERIC RADIUS IN CM
          READ(99,*)RPH
          IPH=1
        ELSEIF(LABEL(1:11).EQ.'PULSAR SPEC') THEN
C         PULSAR SPECTRUM
          READ(99,*)IPULSSP
        ELSEIF(LABEL(1:4).EQ.'CRAB') THEN
C         CRAB SPECTRUM
          ICRAB = 1
        ELSEIF(LABEL(1:14).EQ.'FERLAND NETZER') THEN
C         FERLAND & NETZER AGN SPECTRUM
          READ(99,*)IFERNET
        ELSEIF(LABEL(1:15).EQ.'FERLAND MATHEWS') THEN
C         FERLAND & MATHEWS AGN SPECTRUM
          READ(99,*)IFERMAT
        ELSEIF(LABEL(1:3).EQ.'NLR') THEN
C         NLR AGN SPECTRUM
          READ(99,*)INLR
        ELSEIF(LABEL(1:3).EQ.'NIC') THEN
C         56 NI MASS
          READ(99,*)
        ELSEIF(LABEL(1:11).EQ.'COMPOSITION') THEN
C         SPECIFY CHEMICAL COMPOSITION
          READ(99,9)LABEL
          IF(LABEL(1:5).EQ.'SOLAR') THEN
C           SOLAR COMPOSITION
            ISOLAR=1
          ELSEIF(LABEL(1:6).EQ.'SNCOMP') THEN
C           SUPERNOVA COMPOSITION
             ISNCOMP=1
          ELSEIF(LABEL(1:6).EQ.'AVCOMP') THEN
C     Averaged SHELL COMP from Woosley & Heger 2007 15 of 19 M models
             ISNCOMP=1
c      write(0,*)' Abundance zone? Averaged abundaces=1, Zone 1=Si-S-Ar,'
c     &     ,' 2 = O-Si-S, 3 = O-Ne-Mg'
             read(99,*)rm_abun,avabund,zams,zone
          ELSEIF(LABEL(1:5).EQ.'INPUT') THEN
C           SPECIFY CHEMICAL COMPOSITION
            INPUTCOMP=1
            READ(99,*)(ABIN(L),L=1,12)      
          ENDIF
        ELSEIF(LABEL(1:3).EQ.'VEL') THEN
C         EXPANSION VELOCITY AT THE OUTER BOUNDARY IN CM/S
          READ(99,*)VEXP
        ELSEIF(LABEL(1:7).EQ.'NO BACK') THEN
C         NO BACKGROUND CONT. IN LINE EXC.
          NBACK=0
        ELSEIF(LABEL(1:6).EQ.'OUTPUT') THEN
C         READ PARAMETERS CONTROLLING THE OUTPUT
          READ(99,*)IOSHELL,IORAD1,IORAD2,IOFLUX
        ELSEIF(LABEL(1:10).EQ.'PRINT SPEC') THEN
C         PRINT SPECTRUM FOR EACH ZONE
          READ(99,*)IOSPEC
        ELSEIF(LABEL(1:9).EQ.'DEPARTURE') THEN
C         PRINT DEPARTURE COEFFICIENTS RATHER THEN NUMBER FRACTIONS
          IDEP=1
        ELSEIF(LABEL(1:9).EQ.'PRINT LEV') THEN
C         PRINT LEVEL INFORMATION FOR H I AND He I
          READ(99,*)IOLEVEL
        ELSEIF(LABEL(1:11).EQ.'PRINT EMISS') THEN
C         PRINT EMISSIVITIES FOR EACH ZONE
          READ(99,*)IOEMISS
        ELSEIF(LABEL(1:12).EQ.'PRINT TERMAL') THEN
C         PRINT COOLING AND HEATING FOR EACH ITERATION
          READ(99,*)IOTERM
        ELSEIF(LABEL(1:16).EQ.'PRINT IONIZATION') THEN
C         PRINT IONIZATION RATES FOR EACH ZONE
          READ(99,*)IOION
        ELSEIF(LABEL(1:19).EQ.'PRINT ELECTRON ITER') THEN
C         PRINT IONIZATION RATES FOR EACH ZONE
          READ(99,*)IOELITER
        ENDIF
C       CONSTANT PRESSURE IPRESS = 1
        IF(LABEL(1:3).EQ.'PRE') IPRESS=1
C       CONSTANT DENSITY IDEN = 1
        IF(LABEL(1:3).EQ.'DEN') IDEN=1
C       AGN SPECTRUM 
        IF(LABEL(1:3).EQ.'AGN') IAGN=1
C       PULSAR 
        IF(LABEL(1:3).EQ.'PUL') THEN
C         SPECTRUM PARAMETERS FOR PULSAR AND FREE FREE SPECTRA
           READ(99,*)RLTOT,ALFA,ibreak,ebreak,EMIN,EMAX,teff_ns,nsterm,nstemp
           if(nstemp==2 .or. nstemp==3) then
              tot_lum=rltot
              call ns_atm
           endif
        ENDIF
C
C       IPULS=0 ; only low ionized ions   
C       IPULS=1 ; higher ionized ions included
C
        IF(LABEL(1:3).EQ.'PUL') IPULS=1
C       CIRCUMSTELLAR:
C         SET ICS = 1
C         INOUT = 0
C         SPECIFY RSHOCK, MTOT
C       FOR REVERSE SHOCK 
C         ICS = 1
C         INOUT=1
        IREV=0
        IF(LABEL(1:3).EQ.'CS ') ICS=1
        IF(LABEL(1:3).EQ.'REV') THEN
          IREV=1
          READ(99,*)TECOOL
        ENDIF
C       ONLY LOW IONIZATION IONS INCLUDED: ILOWION = 1
        IF(LABEL(1:3).EQ.'LOW') ILOWION=1
C       STATIC MEDIUM
        IF(LABEL(1:3).EQ.'STA') ISTAT=1
C       SOBOLEV APP. 
        IF(LABEL(1:3).EQ.'SOB') ISOBOL=1
C       ISPH=INDEX SPECIFYING WHETER A SPHERICAL (ISPH=1)
C         OR PLANE (ISPH=0) GEOMETRY SHOULD BE USED
        IF(LABEL(1:3).EQ.'SPH') ISPH=1
C       IPAR=INDEX SPECIFYING WHETHER A PARALLELL (IPAR=1) OR ISOTROPIC
C       BEAM OF IONIZING RAD. SHOULD BE USED.
        IF(LABEL(1:3).EQ.'PAR') IPAR=1
C       IF IUVBLANK = 1 NO UV-CONTINUUM BELOW 3000 A
        IF(LABEL(1:3).EQ.'UVB') IUVBLANK=1
        IF(LABEL(1:4).EQ.'STOP') GOTO 33
      ENDDO
33    CONTINUE
      IPAR=1
      ISPH=1
      IF(ISTAT.EQ.0.AND.ISOBOL.EQ.0) WRITE(0,*)' Sobolev or static?'
      IF(ISTAT.EQ.0.AND.ISOBOL.EQ.0) STOP
      IF(INOUT.EQ.2) WRITE(0,*)' You must specify in or out!'
      IF(INOUT.EQ.2) STOP
C     ************************************************************
C     *****
C     ENERGY INTERVALS
C     *****
C     ************************************************************
      JMIN0=JMIN
      JJ0=JJ
      DO J=JMIN0,JJ0
        IF(E(J).LT.ENMIN) JMIN=J+1
        IF(E(J).LT.ENMAX) JJ=J+1
      ENDDO
C     NINQ= NUMBER OF INTERVALS -1 BETWEEN 13.6 AND 15.3 EV
      NINQ=2
      JJ=jj+NINQ
      IF(INOUT.EQ.0.AND.ICS.EQ.0) CHGH=0.25
C     IF ITE=1 SOLV THERMAL EQ. IF ITE=+1 SKIP THIS.
      ITE=+1
      DO 100 IQ=1,1
      IF(ICS.EQ.1.AND.IPULSSP.NE.1) THEN
        open(19,file='./ATDAT/shockspec.dat',status='old')
        read(19,*)tday
      ENDIF
      DO 100 IU=1,1
      IUY=1
c!!
C     ITER=NUMBER OF ITERATIONS OF STATISTICAL EQUATIONS
      ITER=100
      IGAMM=1
C     AV=POWERLAW OF VELOCITY LAW, V(R)=R**AV
      AV=1.
C     R1= PHOTOSPHERIC RADIUS IN 1E15 CM
C     TEFF = EFFECTIVE TEMPERATURE OF  BLACKBODY CONTINUUM
      CALL SNPAR(TDAY,R1,TEFF,RLUM)
      IF(IBLACKB.EQ.1) THEN
        TEFF=TEFFBB
      ENDIF
      IF(IPH.EQ.1) THEN
         R1=RPH
      ENDIF
C     CIRCUMSTELLAR INTERACTION
      IF(ICS.EQ.1.AND.IPULSSP.NE.1) THEN
C       SHOCK RADIUS IN 1E15 CM
        read(19,*)rshock
        IF(IREV.EQ.1) THEN
          R1=RSHOCK
        ENDIF
      ENDIF
      R1=1.E-15*R1
      TIME=8.64E4*TDAY
C     RADIUS OF THE OUTER BOUNDARY IN CM
      RMAX=VEXP*TIME
      IF(RIN.LT.R1) RIN=R1+1.E-5
C
      IF(ICS.EQ.1.AND.IPULSSP.NE.1) THEN
        read(19,*)mtot
c       temperature of reverse and c-s shocks
        read(19,*)terev,tecs
        read(19,*)revmass
        read(19,*)EJDENS
      ENDIF
      IF(IAGN.EQ.1) MTOT=1.D66
      IF(IAGN.EQ.1) REVMASS=1.D66
      MSUN=MTOT/SOLMA
      RMASS=MTOT/SOLMA
C     DELM = INITIAL MASS INTERVAL
      DELM=(MSUN-RMCEN)/FLOAT(MAX-1-ICEN)
      RMEN=MTOT*2./3.
      RMCO=RMEN/3.
      R0=RMAX
C     GAMMA RAY LUMINOSITY DUE TO THE DECAY OF 56-COBOLT
C     56   CO DECAYS
      GAMLUM=1.36D43*MNI56*EXP(-TDAY/111.26)
C     57 CO DECAYS
C     ASSUME SOLAR RATIO OF 56CO/57CO = 0.0243
      GAM57=3.47D38*(MNI56/0.075)*EXP(-TDAY/391.)
      GAMLUM=GAMLU+GAM57
      RMAX=1.E-15*RMAX
      WRITE(6,9023)R1,RMAX
 9023 FORMAT(' R1= ',E11.4,' RMAX= ',E11.4)
C     NPRINT=INTERVALS BETWEEN PRINTOUTS
      NPRINT=0
C     A1=LOWER LIMIT OF START TEMPERATURE
C     B1=UPPER        D:O
      A1=A1IN
      B1=B1IN

C     TIN=INITIAL GUESS OF TEMPERATURE
      TIN=0.55E4
C     E10=TOLERANCE IN TEMPERATURE
      E10=5.E1
C     E20 = RELATIVE ERROR IN FUNCTION RAD=HEAT-XEL*COOL
      E20=0.01
      MQMAX=MAX
      R1Q=R1
      CALL TRANS(MODNR,NMAX,MAX,ICEN,IREV,mqmax,n,IDEP,IRINNER,ISTEPION,ipar,
     &     RMCEN,NPRINT,DELM,R1,T00,A1IN,B1IN,dminit,rminner,revmass,EJDENS,
     &     TECOOL,TOTCOLUMN,FILLING,min_vel)
 100  CONTINUE
      END


      DOUBLE PRECISION FUNCTION DE(R)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/DENS/D0,R0,RN
      COMMON/TPAR/RIN,DR,R1,TDAY
      COMMON/HYDROS/HSCALE,DEN1,XEL,AMEAN
      RCGS=1.E15*R1*R
      RINCGS=1.E15*RIN
      IF(RCGS.LT.RINCGS) DE=0.
C     DENSITY FROM STRUC
      DE=DEN1
      RETURN
      END

      SUBROUTINE TRANS(MODNR,NMAX,MAX,ICEN,IREV,MQMAX,NIND,IDEP,IRINNER,ISTEPION,IPAR,      
     &     RMCEN,NPRINT,DELM,R1,T00,A1IN,B1IN,dminit,rminner,revmass,EJDENS,
     &     TECOOL,totcolumn,FILLING,min_vel)
      IMPLICIT REAL*8(A-H,O-Z)
c      real*4 etime,tb1,tb2,tb,ts1,ts2,tarr(2)
      character*72 text1,text2,text3,TEXT4,INMOD
      CHARACTER*8 LAB(200)
      REAL*8 MTOT,MNI56
      include "parameters.h"
      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      COMMON/TEXT/TEXT1,TEXT2,TEXT3,TEXT4
      COMMON/MOD/INMOD
      COMMON/IOINF/IOSPEC,IOLEVEL,IOEMISS,IOTERM,IOION,IOELITER
      COMMON/THBAL/TMIN,TMAX,IBAL,IKBAL(1000)
      common/sorc/sor(NE1:NE2),KM(NE1:NE2)
      COMMON/CONT/CONR(NL),CONT(NL)
      COMMON/SECX/CSEC(20),DCSDR(20),CISEC(20),DCISDR
      COMMON/DENS/DEN0,R0,RN
      COMMON/TAUFBO/TAFB(20,3)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/INSPEC/IBLACKB,IPULSSP,IFERNET,INLR,IFERMAT,ICRAB
      COMMON/AGNPARA/DECON,GAMMION,IPRESS,IDEN,ICDEN
      COMMON/PULSAR/RLTOT,ALFA,ebreak,EMIN,EMAX,teff_ns,nsterm,ibreak
      COMMON/FNORM/FNORM
      COMMON/LOWION/ILOWION
      COMMON/SPH/ISPH
      COMMON/INUT/IUY
      COMMON/RQW/TEFF,RQ
      COMMON/TPAR/RIN,DRQ,R1Q,TDAYS
      COMMON/INI/INITH,INIT,INITCA,INITHE,INITFE
      common/initat/initfeii,initfei,initsi1
      common/inttopbase/inittb
      COMMON/ITER/ITE
      COMMON/MPAR/MTOT,MNI56,VEXP
      COMMON/TERMA/TERMAL
      COMMON/QSPEC/GAMMA,ALQSO
      COMMON/LITER/N
      COMMON/IND/I
      COMMON/FRE/NINQ,JMIN,JJ
      COMMON/ITSPH/NITSPH
      COMMON/THER/A1,B1,TIN,E10,E20
      COMMON/SPOT/OTSP(7)
      COMMON/COL/RTE(4),FF,HY,HE,C131,C142,C151,C161,C231,C241,COH,
     &     COHE,C351,C321,C331,C341,C222,C232,C221,C332,C441
      COMMON/EQUIV/RADF,FR(2,100),W(100),CIN(100),WLI(100+NFEL),FB(100)
      COMMON/EQUIH/FRH(2,500),WEH(500),SR(NFEL),WOBS(NL,NL)
      COMMON/GAMPAR/GAMLUM,GAELH,GAHE
      COMMON/T/TES,SS
      COMMON/REHEL/REHE21,REHE22,REC31,REN41
      COMMON/CHG/CHGH,CHGHE,CHTE
      COMMON/GSREC/ALGS(nio),ALEX(nio),ALTOT(nio),RECEX(nio),RECGS(nio)
      COMMON/HEA/CO,PH,PHEO,PHEI,PO(8),PC(6),PN(7),PMG,PSI,PFE
      COMMON/INT/FL(2,NE1:NE2),SI(ncr,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/DIFFP/FH0(NE1:NE2),FHD(NE1:NE2)
      COMMON/TRES/ EL(nio),EK(nio),ELord(14,27),EKord(14,27)
      COMMON/SIK/SK(ncr,NE1:NE2)
      COMMON/PHQ/ZE(nio),GE(nio),ZK(nio)
      common/chec/gec(nio)
      COMMON/OTS/NOTS
      COMMON/EDDFLUX/EDDFLUX(NE1:NE2)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/REC/AL2(7)
      parameter (mh=2,mhe=3,mc=7,mn=8,mo=9,mne=11,mna=12,mmg=13,
     &     mal=14,msi=15,ms=17,mar=19,mca=21,mfe=27,mni=29,
     &     mmat=mh+mhe+mc+mn+mo+mne+mna+mmg+mal+msi+ms+mar+mca+mfe+mni,
     &     ispecy=15)
      common/abc1/alc(mc),aln(mn),alo(mo),alne(mne),alna(mna),
     &     almg(mmg),alal(mal),alsi(msi),alsu(ms),alar(mar),alca(mca),
     &     alfe(mfe),alni(mni)
      COMMON/HYDROS/HSCALE,DEN1,XQL,AMEAN
      COMMON/RADIE/R(MD)
      COMMON/DXA/DR(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/PHY/DEN(MD)
      COMMON/MASSES/RM(MD)
      COMMON/ABUN/AB(20)
      COMMON/TEM/TE(MD)
      COMMON/ION/XB(MD,nio)
      COMMON/ABU/XA(2,nio)
      COMMON/NHY/NH
      COMMON/NLEV/e00,NION,NHY,NP1H,NMA,NMI
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &     BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/HYPOP/XNH(NL)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/A5/TAUE,ALFAQ,EN,TSN,XL40,TEXA,RMAX
      COMMON/A7/C3,C33
      COMMON/A10/AV
      COMMON/A11/R11,CV,FLUX
      COMMON/A12/ITER
      COMMON/PHEAT/PHE(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),REC(NL),PHOT(NL),DRDT(NL),RECT(5,NL),PHET(NL)
      COMMON/HRAT/RHYD,ZEL,XEL,HYR,HEAT,COOL
      COMMON/A14/CQ(NL,NL),CI(NL),G(NL),EI(NL),AQ(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/COLH/DXI,COLH(5,NL),COLHVT(5,NL),COLHVTT(5,NL),TAUB
      COMMON/A19/EMIS(NL,NL),ESC(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      COMMON/A31/TOPI(10,NL,NL),EMHI(NL,NL)
      COMMON/COLD/OPT(nio),COLTOT(nio),COLTOTT(nio),SRED(NFEL+100)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),ipara
      real*8 jmean
      COMMON/DTAU/JMEAN(NE1:NE2)
      COMMON/IOPAR/IOSHELL,IORAD1,IORAD2,IOFLUX
      COMMON/CINOUT/INOUT,IPULS
      COMMON/CICS/RSHOCK,ics
      common/revs/terev,tecs
      COMMON/FE3/WLFE3(30,30),AFE3(30,30),EMFE3(30,30)
      COMMON/FE4/WLFE4(30,30),AFE4(30,30),EMFE4(30,30)
      COMMON/PROFFE/WLFEIII(300),WLFEIV(300),FEIII(300,MD),FEIV(300,MD)
      COMMON/FELEV/NFEII,nfei
      COMMON/WFE/WOBFE(NL,NL),IQPOPH
      COMMON/EMHY/RECEM(5,NL),TWOPH,TWOPHHEI,COHI,PHQ,PHT,POQ,POT
      common/heiires/x2he,he2coll12,he2jla,he2j2g,he2jrec(4),
     &                                                he2jem(7,4)
      common/test1/qa(5,20,20),rateph(5,20),tsav(5,20,20),esav(5,20,20)
      COMMON/LYBFLUOR/OIFLUOR,OIFLUORH,OIFLUORO,OILYB,OIREC
      DIMENSION DRA(MD),COLD(nio),COLDN(nio),FBQ(100),FBI(100,MD)
      DIMENSION DLO(100),WEQ(100),WEQH(100),FTOTAL(101),FREK(101),
     &     ILI(100)
      DIMENSION FL0(2,NE1:NE2),WREC(30),WRHE(10),HER(10),OXR(10),
     &                                                       WROX(10)
      DIMENSION EMTOT(NL,NL),CNRTOT(NL),CONTOT(NL),CITOT(100),
     &            WFEII(NL,NL),DWFEII(NL,NL),CITOTFE(500),WLFET(500)
      DIMENSION DRECHI(NL),RECHITOT(NL),BX(NL)
      DIMENSION DWFEIII(30,30),WFEIII(30,30),DWFEIV(30,30),
     &            WFEIV(30,30)
      DIMENSION DCIT(100),DFBQ(100),RMO(MD),TEINT(MD),
     &          SIRED(100+NFEL,MD)
      DIMENSION BN(50,5),THA(50),THB(50),SO(50),FU(251),ABUN(MD,20)
      DIMENSION BQH(NL),BQ(NL),BQCA(NL),BQHE(NL),BQFE(NL),BQMG(NL)
      DIMENSION DLW(100),OLDFLU(NE1:NE2)
      DIMENSION BBOL(10,MD,NL),axsi(20),gaion(20)
      DIMENSION WLF1(71),IABU(20),ABOLD(20),KLI(300),KKMAX(5),SPL(100)
      DIMENSION TREC(100)
      DIMENSION EJMIN(7),EJMAX(7),SION(7),SIONO(7)
      DATA EJMIN/3.4,13.5987,24.5,54.4,100.,200.,1.E3/
      DATA EJMAX/13.5987,24.5,54.4,200.,1.E3,1.E4,1.E6/
      DATA IABU/0,0,1,1,0,1,1,1,0,0,1,0,0,0,6*0/
      DATA KKMAX/3,9,12,0,5/
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DATA SOLMA/1.989E33/
      real*8 tot_line_lum,min_vel
      common/linelum/tot_line_lum(14,27,401)      
      integer nionstage
      common/num_ion_stages/nionstage(14)      
      data nionstage/1,2,6,7,8,10,3,3,4,14,16,7,6,26/
      common/line_em/cinx(14,26,401),taulinex(14,26,401),wlix(14,26,401)
     &     ,ilabcfx(14,26,401)
      common/kmaxpop/kmaxp(14,27)
      common/ionx/xion(md,14,27)
      common/populations/popul(14,27,400)
      real*8 line_cool
      common/line_c/line_cool(14,27)
      integer z2cf(26)      
      data z2cf/1,2,0,0,0,3,4,5,0,6,7,8,9,10,0,11,0,12,0,13,0,0,0,0,0,14/

      KKMAX(4)=NFEII
C     ***********************************************************
C     *****
C     INITIALIZE PARAMETERS
C     *****
C     ***********************************************************
      initstruc=1
      inittb=1
      npri=0

c     nummering of ions from old CF system

      call cfenum
      do iel=1,14
         do ion=1,nionstage(iel)
c     do k=1,kmax(iel,ion)
            do k=1,401
               tot_line_lum(iel,ion,k)=0.
               wlix(iel,ion,k)=0.
            enddo
         enddo
      enddo

C     ***********************************************************
C     *****
C     INITIALIZE PARAMETERS
C     *****
C     ***********************************************************
      write(0,*)' max number of shells?'
      read(5,*)iqmax
      INI=0

C
C     ELECTRON DENSITY EXCEPT FOR H
C
      ZEL=AB(2)+AB(3)+AB(4)+AB(6)+AB(7)+AB(8)+AB(9)
      TSN=1.E-10
      DO J=JMIN,JJ
         DO KK=1,2
            FL0(KK,J)=0.
            FL(KK,J)=0.
         enddo
         DO KK=1,100
            SI(KK,J)=0.
            SK(KK,J)=0.
         enddo
         F0(J)=0.
         DO I=1,MAX
            JMEAN(J)=0.
            EM(I,J)=0.
            TAU(I,J)=0.
            TAUTOT(I,J)=0.
            FD(I,J)=0.
         enddo
      enddo
      IF(IUY.eq.1) then
         DO IKA=1,NL
            BB(IKA)=1.E-10
            BOLCA(IKA)=1.
            BOL(IKA)=1.+0.1*FLOAT(IKA-1)
         enddo
      else
         BB(1)=1.E0
      endif
      NPR=NPRINT
      IPARA=IPAR
      IOIND=1
C
C     R15=PHOTOSPHERIC RADIUS IN 1E15 CM
C
      R15=R1
      R11=1.E4*R15
C
C     CV = VELOCITY OF LIGHT / VELOCITY AT THE POSITION OF THE  BOUNDARY
C
      CV=3.E10/VEXP
      RMAXCGS=RMAX*1.E15	
      RMCV=RMAX*1.E15*CV
C
C     OPTICAL DEPTH CONSTANT
C
      C33=1.3263E-21*CV*RMAX
      VTERM=1.2856E4
      IF(ISTAT.EQ.1) C33=2.2448E-26/VTERM
C
C     COUNTER OF ITERATIONS IN RAD. TRANSFER CALC.
C
      N=1
      NITSPH=N
      IF(INOUT.EQ.0) NITSPH=2
      KQ=1
C
C     NUMBER OF LEVELS IN H-MODEL ATOM
C
      NH=3
      NP1=NH+1
      C4=1.44E8
      NUM=28
      NOTS=1
      DO K=1,100
         ZK(K)=0.
         ZE(K)=0.0E0
         GE(K)=0.0E0
         DO I=1,MAX
            XB(I,K)=0.
         enddo
         XA(1,K)=0.0E0
         XA(2,K)=0.0E0
      enddo
C
C     INITIAL GUESSES OF IONIC ABUNDANCIES
C
      XA(2,2)=1.
      XA(2,12)=1.
      XA(2,14)=1.
      XA(2,18)=1.
      XA(2,26)=1.
      XA(2,40)=1.
      XA(2,43)=1.
      XA(2,51)=1.
      XA(2,54)=1.
      XA(2,56)=1.
C
C     ASSUME CASE B INITIALLY
C
      OTSP(1)=0.
      OTSP(2)=0.
      OTSP(3)=0.
      DEL(1)=AB(1)+AB(2)+AB(3)+AB(4)+AB(5)+AB(6)+AB(8)
      DEL(1)=0.1
      I=1
C
C     DENSITY AT INNER BOUNDARY
C
C     **************************************************************
C     *****
C     R(1)=  BOUNDARY RADIUS IN TERMS OF THE PHOTOSPHERIC RADIUS
C     R(I) = R(CGS)/R(PHOT)
C     DR(I) = SHELL THICKNESS IN CM
C     DRA(I) = SHELL THICKNESS RELATIVE TO THE  BOUNDARY RADIUS.
C     DEN(I) = NUMBER DENSITY IN SHELL I.
C     *****
C     **************************************************************
C
C     ISTRU=1 IF DENSITY FROM WOOSLEY WR-MODEL, 0 IF UNIFORM OX & HE
C
      ISTRU=1
      RSUN=MTOT/SOLMA
      RM(1)=RSUN/1.0001
      IF(INOUT.EQ.1) RM(1)=RMINNER
      IF(IREV.EQ.1) THEN
        RM(1)=RSUN/1.0001
        RI15=RSHOCK
      ENDIF
      IF(IRINNER.EQ.1) THEN
        RI15=RIN
      ENDIF
      RM1CGS=RI15*1.E15
      CALL STRUC_HW(initstruc,filling,min_vel,RM(1),RM1CGS,DENR,RI15)
      write(6,*)'AB(6),AB(12),AB(14),AB(11),AB(8) c ',AB(6),AB(12),
     &        AB(14),AB(11),AB(8)
      WRITE(6,9725)(AB(IK),IK=1,2),AB(4),AB(5),AB(3),AB(13),AB(10),
     &     AB(7),AB(9),AB(6),AB(12),AB(14),AB(11),AB(8)
      write(6,*)' '
 9725 format('Abundances: H',1pe12.3,' He',e12.3,' C',e12.3,' N',e12.3,
     &     ' O ',e12.3,' Ne ',e12.3,' Na',e12.3,' Mg ',e12.3,
     &     ' Al',e12.3,' Si',e12.3,' S',e12.3,' Ar',e12.3,
     &     ' Ca',e12.3,' Fe',e12.3)
      DENR=DENR/FILLING
      IF(IRINNER.EQ.1) THEN
         RI15=RIN
         write(6,*)'rin,r1,ri15,rin15 ',rin,r1,ri15,rin15
      ENDIF
      
      IF(ICDEN.EQ.1) THEN
        DENR=DECON
        DENR=DENR/FILLING
        DEN1=DENR
      ENDIF
      if(irev.eq.1) RM(1)=rminner
      R(1)=RI15
      r(1)=min_vel*1.e5*tdays*8.64e4
c     all r() in units of r15*1.e15. Here r15=1.0, so r(1)=r/1e15
      r(1)=r(1)/(r15*1.e15)
      RQ=R(1)
      RCGS=1.E15*R(1)*R15
      write(6,*)'tdays,  inner vel (km/s), R(1)/1.e15 ',tdays,min_vel,R(1)*R15
      RMO(1)=RM(1)
      RMASS=0.
      DEN(1)=DENR
      TE(1)=TIN
      IF(ICS.NE.1) THEN
c     only for CSM ICS=1
C     TOTAL COLUMN DENSITY
        COTOT=0.
        COTOTN=0.
        RMI=RSUN/1.0000001
        RM1CGS=RIN*1.E15
        CALL STRUC_HW(initstruc,filling,min_vel,RMi,RM1CGS,DENR,RI15)
        DENR=DENR/FILLING
        DO 18 Iq=1,1000
          ROLD=RI15
          RMI=-Iq*RSUN/(1000)+RSUN+1.E-5
          RM1CGS=RI15*1.E15
          CALL STRUC_HW(initstruc,filling,min_vel,RMi,RM1CGS,DENR,RI15)
          DENR=DENR/FILLING
          DEN1=DENR
          IF(ICDEN.EQ.1) THEN
            DENR=DECON
            DENR=DENR/FILLING
            DEN1=DENR
         ENDIF
          COTOTN=1.E15*DENR*(-RI15+ROLD)+COTOT
18      COTOT=1.E15*DENR*(-RI15+ROLD)*AMEAN*AMU+COTOT
      ENDIF
      INI=0
C     *************************************************************
C     EVALUATTION OF CROSSECTIONS FOR EACH ENERGY
C     *************************************************************
      E00=11.871
      NHY=NH
      NMA=NH
      NION=1
      CALL CROXY
      CALL CROSSSECT
      CALL CRFEII
      CALL ATDATO
      CALL ATDATFE
      CALL IONLABEL(LAB)
      NI=1
      DO  Iq=2,NFEII
         DO  J=1,Iq-1
            IF (AQ(iq,J).GT.0..and. ni.le.1015)THEN
               WLI(81+NI)=WL(Iq,J)
               NI=NI+1
            ENDIF
         ENDDO
      ENDDO
      CALL ATDAT
      CALL CRCA
      call read_coll
      do 1534 k=1,1
      iqsx=1
      do j=jmin,jj
         if(iqsx.eq.1.and.sk(k,j).gt.0.) iqsx=0
      enddo
      iqsxc=1
      do j=jmin,jj
      if(iqsxc.eq.1.and.si(k,j).gt.0.) iqsxc=0
      enddo
 1534 continue
C     ***********************************************************
C     *****
C     INFALLING SPECTRUM
C     *****
C     ***********************************************************
C     FL=MEAN INTENSITY PER UNIT OF ENERGY (ERG/(CM**2 STER EV))
      IF(IPAR.EQ.1) GEOM=4.
      IF(IPAR.NE.1) GEOM=2.
c      write(6,*)' IPAR,GEOM,IPULSSP,ICRAB',ipar,geom,ipulssp,icrab
      totl=0.
      DO J=JMIN,JJ
C
C       MEAN INTENSITY. Note fl = J = fmean
C
        FL(1,J)=FMEAN(j,1,E1(J))
        FL(2,J)=FL(1,J)
        JMEAN(J)=FL(1,J)        
        OLDFLU(J)=JMEAN(J)
        if(e1(j)>13.6.and.e1(j)<2.e4) then
           totl=totl+16.*pi**2*(r15*r(1)*1.e15)**2*fl(1,j)*(e(j+1)-e(j))
        endif
        WRITE(116,9245)J,r15*r(1),E1(J),FL(2,J),totl
9245    FORMAT('J,R,E,FMEAN=J,tot_ion_lum ',I4,1E12.3,10E12.3)
        F0(J)=FL(1,J)
      ENDDO
      write(6,9362)totl
 9362 format('Total ionizing luminosity 13.6 - 20keV ',1pe12.3)

c      stop
C
C     CALC. NUMBER OF IONIZING PHOTONS
C
      DO K=1,7
        SIONO(K)=SION(K)
        SION(K)=0.
        DO J=JMIN,JJ
          IF(E1(J).GT.EJMIN(K).AND.E1(J).LT.EJMAX(K)) THEN
            SION(K)=SION(K)+FL(1,J)*(E(J+1)-E(J))/E1(J)
          ENDIF
        ENDDO
      ENDDO
      SIONXO=SIONX
      SIONX=SION(5)+SION(6)
C
C     IONIZATION PARAMETER
C     6.242E41 = 1.E15**2/ELCH
C
      V=SION(2)+SION(3)
      VV=SION(4)+SION(5)+SION(6)+SION(7)
      VP=6.242D41*4.*GEOM*PI**2.*(R15*R(1))**2.*V
      VVP=6.242D41*4.*GEOM*PI*PI*(R15*R(1))**2.*VV
      U1=4.*PI*V/(den(1)*ELCH*3.E10)
      U2=4.*PI*VV/(DEN(1)*ELCH*3.E10)
      WRITE(6,7778) VP,VVP
C     WRITE(6,525)
      WRITE(6,7779) U1,U2

 7778 FORMAT(1X,'NUMBER OF IONIZING PHOTONS: S(13.6-54 EV)=',1PE13.6
     &,' S( >54 EV)=',E13.6)
 7779 FORMAT(1X,'IONIZATION PARAMETER: U1(13.6-54 EV)=',1PE13.6
     &,' U2( >54 EV)=',E13.6)

C
C     ESTIMATE THE INITIAL STEP IN MASS (IN SOLAR UNITS)
C
C     H II REC AT 2.E4 K
      ALREC=1.4E-13
      DMINITH=VP*AMEAN*AMU/(ALREC*DEN(1))
C     HE III REC AT 2.E4 K
      ALREC=9.08E-13
      DMINITHE=VVP*AMEAN*AMU/(ALREC*DEN(1))
      DMINIT=DMINIT*DMIN1(DMINITH,DMINITHE)/SOLMA
c      WRITE(6,9278)DMINIT,DMINITH,DMINITHE
c!!!  new
      dminit=0.00001*dminit
      dminit=1.e-11

C     CALCULATE FLUX
      DO J=JMIN,JJ
        CALL RADTRANPULS(1,J,TAU(I,J),S)
      ENDDO
      V=0.
      VV=0.
c mean intensity above 54 eV            
      DO J=47,JJ
         VV=VV+(FH0(J)+FHD(J))*(E(J+1)-E(J))
      enddo
c mean intensity 13.6-54 eV      
      DO J=2,46
         V=V+(FH0(J)+FHD(J))*(E(J+1)-E(J))
      enddo
c ionizing       
      vion=0.
      v13p6_1e4=0.
      DO J=2,jj
         Vion=Vion+(FH0(J)+FHD(J))*(E(J+1)-E(J))
         if(e(j)< 1.e4) then
            v13p6_1e4=vion
         endif
      enddo
c Luminosity between these energies      
      V=1.E30*4.*GEOM*PI**2.*(R15*R(1))**2.*V
      VV=1.E30*4.*GEOM*PI*PI*(R15*R(1))**2.*VV
      Vion=1.E30*4.*GEOM*PI*PI*(R15*R(1))**2.*Vion
      V13p6_1e4=1.E30*4.*GEOM*PI*PI*(R15*R(1))**2.*v13p6_1e4
c Total luminosity including diffuse emission      
      VVV=0.
c Total luminosity from primary source only
      VVV1=0.
      DO J=JMIN,JJ
        VVV=VVV+(FH0(J)+FHD(J))*(E(J+1)-E(J))
        VVV1=VVV1+FH0(J)*(E(J+1)-E(J))
      ENDDO
      ERADTOT=1.E30*4.*GEOM*PI*PI*(R15*R(1))**2.*VVV
      ERADTOT1=1.E30*4.*GEOM*PI*PI*(R15*R(1))**2.*VVV1
      WRITE(6,903) ERADTOT,ERADTOT1,V,VV,vion,v13p6_1e4,DERAD,DERADDV,
     &     (CIT+FBT),R(I)*R15*1.e15
 903  FORMAT(1X,'ENERGY (TOT) ',1PE12.3,' (PRIM)',E12.3,' 13.6-54 ',E12.3,' > 54 EV '
     &     ,E12.3,' > 13.6 eV ',e12.3,' 13.6 - 1e4 eV ',e12.3' DE ',E12.5,' DEDV ',E12.3,
     &     ' DEDV(LINES) ',E12.3,' R=',E13.6)
      WRITE(6,525)
 525  FORMAT(' ')
      DO Iq=1,NL
        DO K=1,5
          COLHVTT(K,Iq)=2.
          IF(Iq.EQ.1) COLHVTT(K,Iq)=1.D33
        ENDDO
      ENDDO
C
C     START OF ITERATION LOOP FOR THE DIFFUSE FIELD
      ICONTI=0
555   CONTINUE
      INIT=1
      INITH=1
      INITHE=1
      INITCA=1
      INITFE=1
      INITMG=1
      initsi1=0
      initfei=0
      initfeii=0

C
C     IF ICONTI = 1 READ OLD VALUES OF N , IMAXP, RM AND FD FROM FILE 14
C
      IF(ICONTI.NE.1) GOTO 7378
      REWIND 14
      READ(14) N
      NITSPH=N
      N=N+1
      write(6,*)'read imaxp',imax
      READ(14)IMAXP
      MAX=IMAXP
      READ(14)RM
      READ(14)FD
      ICONTI=0
      WRITE(6,*)' FLUXES READ FROM UTIN'
7378  CONTINUE
      DXI=DR(2)
      REHE1=0.
      REN1=0.
      REC1=0.
      REHE2=0.
      rehe3=0.
      rehe4=0.
      REC31=0.
      REN41=0.
      DO J=JMIN,JJ
         TAU(1,J)=0.
         FPRIM=FMEAN(j,1,E1(J))
         FL(2,J)=FPRIM+FD(1,J)
         FL(1,J)=FL(2,J)
      enddo
      TEL=0.
      COLHYD=0.
      DO Iq=1,100
C       ASSUME INFINITE OPTICAL DEPTHS IN THE DIRECTION OUTWARDS
C            FROM THE SOURCE
        COLTOTT(iq)=1.D100
        COLTOT(Iq)=0.
        COLD(Iq)=0.
        COLDN(Iq)=0.
        CITOT(Iq)=0.
      ENDDO
      DO IU=2,NFEII
        DO IL=1,IU-1
          WFEII(IU,IL)=0.
        ENDDO
      ENDDO
      NFEIII=30
      DO IU=2,NFEIII
        DO IL=1,IU-1
          WFEIII(IU,IL)=0.
        ENDDO
      ENDDO
      NFEIV=22
      DO IU=2,NFEIV
        DO IL=1,IU-1
          WFEIV(IU,IL)=0.
        ENDDO
      ENDDO
      RECHETOT=0.
      DO K=1,5
        RECHITOT(K)=0.
      ENDDO
      TWOPHTOT=0.
      TWOPHHEIT=0.
      FREEFREE=0.
      HEATTOT=0.
      COOLTOT=0.
      GAMTOT=0.
      DO I=1,NL
        CONTOT(I)=0.
        CNRTOT(I)=0.
        PHOT(I)=0.
        REC(I)=0.
        DO K=1,5
          COLH(k,I)=0.
          COLHVT(k,I)=1.
        ENDDO
        DO J=1,NL
          ESC(I,J)=0.
          TTOT(I,J)=1.E20
          TOP(I,J)=0.
          EMTOT(I,J)=0.
          EMIS(I,J)=0.
        ENDDO
      ENDDO
      DO K=1,30
        WEQ(K)=0.
        WEQH(K)=0.
        WREC(K)=0.
        IF(K.LE.10) WRHE(K)=0.
        IF(K.LE.10) WROX(K)=0.
        W(K)=0.
      ENDDO
      DO K=1,100
        FBQ(K)=0.
        CIN(K)=0.
        FB(K)=0.
      ENDDO
      DO K=1,40
        WEH(K)=0.
      ENDDO
      TAUB=0.
      citt1=0.
      citt2=0.
      citt3=0.
      citt4=0.
      citt5=0.
      citt6=0.
      citt7=0.
      citt8=0.
      citt9=0.
      citt10=0.
C     ***********************************************************
C     *****
C     CALCULATE PHOTOIONIZATION RATES FOR SLAB 1
C     *****
C     ***********************************************************
      CALL RATE(1.d0)
      DO J=JMIN,JJ
        TAU(1,J)=0.0E0
      ENDDO
      RMASS=1.E-6
c      RMASS=RM(1)
      COLUMN=COTOT
      IF(INOUT.EQ.1) COLUMN=0.
      IF(ICS.EQ.1) COLUMN=0.
      COLUMNN=COTOTN
      IF(INOUT.EQ.1) COLUMNN=0.
      IF(ICS.EQ.1) COLUMNN=0.
      TAUELE=0.
      M=2
      MQW=0
      MOXY=0
      IFE=0
      ISI=0
      IOMG=0
      IOC=0
      IF(IUY.EQ.1) TOLD=TIN
      IF(IUY.EQ.1) TOLDH=TIN
      IF(IUY.EQ.1) TOLDHE=TIN
      IF(IUY.EQ.1) TOLDO=TIN
      IF(IUY.EQ.1) TOLDCA=TIN
      IF(IUY.EQ.1) TOLDFE=TIN
      IF(IUY.EQ.1) TOLDMG=TIN
      IF(IUY.EQ.1) TOLDQ=TIN
      IF(IUY.EQ.1) TOLDHQ=TIN
      IF(IUY.EQ.1) TOLDHEQ=TIN
      IF(IUY.EQ.1) TOLDOQ=TIN
      IF(IUY.EQ.1) TOLDCAQ=TIN
      IF(IUY.EQ.1) TOLDFEQ=TIN
      IF(IUY.EQ.1) TOLDMGQ=TIN
C     ***********************************************************
C     *****
C     MAIN LOOP IN DEPTH FOR EACH ITERATION
C     *****
C     ***********************************************************
      DO 700 I=2,MAX
         A1OLD=A1
         B1OLD=B1
         NIPOP=0
         COLOLD=COLUMN
C     H DOMINATED MODELS
C     IF(INOUT.EQ.0)DELRM=-DELM
C     MIXED MODELS
         IF(INOUT.EQ.0.AND.I.EQ.2) DELRM=-DMINIT
c     IF(AB(1).GT.0.5) DELRM=-0.2
c     IF(I.LE.(ICEN+1))DELRM=-.5/FLOAT(ICEN)
c     IF(I.LE.(ICEN+1).AND.INOUT.EQ.1)DELRM=RMCEN/FLOAT(ICEN)
c     IF(INOUT.EQ.1.AND.TAUTOT(I-1,2).LT.0.001) 
c     &                              delrm=10.**((i-1)*0.025)*dminit
         IF(INOUT.EQ.1.AND.I.EQ.2) DELRM=DMINIT
         RM0=RM(1)
         IDEL=0
C     ***********************************************************
C     *****
C     CALC. NEW STEP SIZE IN MASS
C     *****
C     ***********************************************************
C     
C     CALC. NUMBER OF IONIZING PHOTONS
C     
         DO K=1,7
            SIONO(K)=SION(K)
            SION(K)=0.
            DO J=JMIN,JJ
               IF(E1(J).GT.EJMIN(K).AND.E1(J).LT.EJMAX(K)) THEN
                  SION(K)=SION(K)+FL(2,J)*(E(J+1)-E(J))/E1(J)
               ENDIF
            ENDDO
         ENDDO
         SIONXO=SIONX
         SIONX=SION(5)+SION(6)
C     FLUX WEIGHTED CHANGE IN FLUX
         if(i >=4) then
            DMFLOP=0.
            DO J=JMIN,JJ
               IF(TAUTOT(I-1,J).LT.700.) THEN
                  OPAC=EXP(-TAUTOT(I-1,J)/2.)
               ELSE
                  OPAC=0.
               ENDIF
               DFLOP=ABS((1.-EXP(-TAU(I-1,J)))*OPAC)
               IF(DFLOP.GT.DMFLOP) THEN
                  DMFLOP=DFLOP
                  STFAC=1.
                  IF(TAUTOT(I-1,J).LT.5.) THEN
                     STFAC=5.
                  ENDIF
                  EMJ=E1(J)
                  JMOP=J
               ENDIF
            ENDDO
         endif
         DMFLOP=DMFLOP*STFAC
         IF(N.GT.1) GOTO 805
         IF(I.LE.3) GOTO 805
c!!   IF(TAUTOT(I-1,3).LT..0001) GOTO 805
         write(6,*)i,te(i-1),te(i-2)
         write(6,*)del(i-1),del(i-2)
         if(te(i-2).ne.0.) then
            DLTE=CHTE*ABS(TE(I-1)-TE(I-2))/TE(I-2)
            DLXE=0.5*CHTE*ABS(DEL(I-1)-DEL(I-2))/DEL(I-2)
            DLA=ABS(LOG10(TAUTOT(I-1,2))-LOG10(TAUTOT(I-2,2)))
         endif
         IF (tautot(i-1,3).gt.1.e-4) THEN
            DLA=0.
            DLB=0.
            DLC=0.
            IF(TAUTOT(I-1,3).LT.1.E3) 
     &           DLA=ABS(LOG10(SION(5))-LOG10(SIONO(5)))
            IF(TAUTOT(I-1,23).LT.1.E3) 
     &           DLB=ABS(LOG10(SION(3))-LOG10(SIONO(3)))
            IF(TAUTOT(I-1,47).LT.1.E3) 
     &           DLC=ABS(LOG10(SION(4))-LOG10(SIONO(4)))
            DLD=ABS(LOG10(SIONX)-LOG10(SIONXO))
         ENDIF
         CHG=CHGH
         IF (INOUT.EQ.1)DL=DMAX1(DLTE,DLA,DLB,DLC,DLD,DLXE)
         IF (INOUT.EQ.0)DL=DMAX1(DLTE,DLA)
C!!   
         DL=DMAX1(DMFLOP,DLTE,DLXE)
         
         WRITE(6,9628)i,te(i-1),EMJ,DMFLOP,DLTE,DLXE,DL,
     &        TAUTOT(I-1,JMOP),TAU(I-1,JMOP)
         if(dl<0.5e-3) then
            dl=0.5e-3
            write(6,*)' DL adjusted ',dl
         elseif(dl > 0.05) then
            dl=0.05
            write(6,*)' DL adjusted down ',i,dl
         endif
 9628    format(' STEP: ',i4,1pe12.3,20e11.3)
         IF(DL.EQ.0.) DL=CHG
         DELOLD=DELRM
         DELRM=DELRM*DMIN1(1.2D0,CHG/DL)
C!!   FOR MIXED MODEL
c     UNIFORM INCREASE IN STEP SIZE
         IUNIF=1
         IF(IPULS.EQ.1) IUNIF=0
         IF(IUNIF.EQ.1) DELRM=1.25*DELOLD
         RAT=DELRM/DELOLD

         IF(ABS(DELRM).GT.(RSUN/200.).AND.INOUT.EQ.0) DELRM=-RSUN/200.
         IF(ABS(DELRM).GT.(RSUN/200.).AND.INOUT.EQ.1) DELRM=RSUN/200.
         IF(TOTCOLUMN.LE.1.E30) THEN
            DCOL=DELRM*SOLMA/(4.*PI*RCGS**2*AMEAN*1.67E-24)
            IF(DCOL.GT.TOTCOLUMN/25.) DELRM=DELRM*TOTCOLUMN/(25.*DCOL)
         ENDIF

         if(i>2 .and. taumax > 1.) then
            write(6,*)'RMASS test ',i,taumax,rmass,delrm/rmass
            delrmp=delrm
            if(delrm > 0.5d-2*rmass) then
               delrm=0.5e-2*rmass
c test
               delrm=1.5e-2*rmass
               write(6,*)'RMASS taumax ',i,taumax,rmass,delrmp,delrm
            endif
         endif
 805     CONTINUE


 7832    CONTINUE

c!!! test only. below works.,
c     delrm=dmax1(0.025*rm(i-1),delrm)
         if(del(i-1) > 1.1.or.i==2) then
c     933 etc
c         if(del(i-1) > 1.05.or.i==2) then
c         if(del(i-1) > 1.0.or.i==2) then
            delrm=dmax1(0.03*rm(i-1),delrm)
c     test
            delrm_new=dmax1(0.04*rm(i-1),delrm)          
         else
            delrm_old=delrm
            if(del(i-1) < 1.0) then
               delrm=dmax1(0.005*rm(i-1),delrm)
            endif               
            delta_xe=abs(del(i-1)/del(i-2)-1.)
            delrm_new=1.e-2*delrm/delta_xe
c     test
            delrm_new=2.e-2*delrm/delta_xe            
            if(te(i-1) < 175.) then
c small change in T_e for T_e < 150K               
               dl_te=ABS(TE(I-1)-TE(I-2))/TE(I-2)
               delrm=0.015*delrm_old/dl_te
            else
               dl_te=0.
               delrm=dmin1(delrm_new,2.*delrm_old)
            endif
         endif

         IF(N.EQ.1) RM(I)=RM(I-1)+DELRM
         IF(RM(I).LT.0.) GOTO 9365
         DELRM=ABS(RM(I)-RM(I-1))
         AMOLD=AMEAN
         IF(I.EQ.2) AMOLD=0.
         DO K=1,20
            ABOLD(K)=AB(K)
            IF(I.EQ.2) ABOLD(K)=0.
         ENDDO
C     RADIUS AND DENSITY AT MASS = RM(I)
         RM1CGS=R15*1.E15*R(I-1)
         CALL STRUC_HW(initstruc,filling,min_vel,RM(i),RM1CGS,DENR,RI15)
         DENR=DENR/FILLING
         IF(ICDEN.EQ.1) THEN
            DENR=DECON
            DENR=DENR/FILLING
            DEN1=DENR
         ENDIF
         IF(IREV.EQ.1.AND.I.EQ.2) THEN
            DENREV=DENR
         ENDIF

         RMO(I)=RM(I)
         R(I)=RI15/R15      

         IF(ICS.EQ.1) DELRM=1.1*DELRM
         IF(IPULS.EQ.1.AND.I.EQ.2.OR.IDEN.EQ.1) THEN
            DRCGS=ABS(DELRM)*SOLMA/(4.*PI*DENR*AMEAN*AMU*RCGS**2)
            R(I)=DRCGS/(R15*1.E15)+R(I-1)
         ENDIF
         IF(R(I).LE.1.AND.INOUT.EQ.0) GOTO 7832
         DRA(I)=ABS(R(I)-R(I-1))
         DRCGS=DRA(I)*R15*1.E15
         DR(I)=DRA(I)*R15*1.E15
         DXI=DR(I)
         RCGS=1.E15*R(I)*R15
         IF(ISTRU.EQ.0) GOTO 3828
c!!   don't do this
         IF(ICDEN.EQ.1.AND.I.LT.3) DENR=DECON 
         IF(IREV.EQ.1.AND.I.EQ.2) THEN
            DENR=DENREV
         ENDIF
C     CONSTANT PRESSURE
         IF(IPRESS.EQ.1.AND.I.GE.3) DENR=DEN(2)*(1.+DEL(2))*TE(2)/
     &        ((1.+DEL(I-1))*TE(I-1))
c     const dens
         IF(IDEN.EQ.1.AND.I.GE.3) DENR=DEN(2)
 3828    DECGS=DENR*AMU*AMEAN
         DENQ=DENR
         DEN1=DENR
         RI=RCGS
C     
C     CALCULATE THE COLUMN DENSITY, THE GAMMA RAY OPTICAL DEPTH
C     AND THE GAMMA RAY HEATING PER MAS,S GAHE.
C     
         IF(INOUT.EQ.1) COLUMN=COLOLD+DENQ*DR(I)*AMU*AMEAN
         IF(INOUT.EQ.0) COLUMN=COLOLD-DENQ*DR(I)*AMU*AMEAN
         IF(INOUT.EQ.1) COLUMNN=COLOLDN+DENQ*DR(I)
         IF(INOUT.EQ.0) COLUMNN=COLOLDN-DENQ*DR(I)
C     
C!!   THIS MAY NOT ALWAYS BE APPLICABLE, ONLY WHEN YOU WANT TO 
C     IGNORE THE ABSORPTION OF THE GAMMA RAYS BY THE INTERIOR MASS
C     
         IF(I.EQ.2.AND.INOUT.EQ.1) COLUMN=0.
         IF(I.EQ.2.AND.INOUT.EQ.1) COLUMNN=0.
         DMASS=4.*PI*RCGS**2*DR(I)*AMU*AMEAN*DENQ/SOLMA
         DMASS=FILLING*DMASS
         IF(INOUT.EQ.1) RMASS=RMASS+DMASS
         IF(INOUT.EQ.0) RMASS=RMASS-DMASS
         IDEL=IDEL+1
 7536    TAUR=0.03*COLUMN
         GAHE=0.03*AMEAN*DENQ*GAMLUM*EXP(-TAUR)/(4.*3.14*RCGS**2)
         GAHE=GAHE*AMU
c         WRITE(6,9286)R(I),DRA(I),DRcgs,RCGS,RMASS,rm(i),TAUR,DENR
 9286    FORMAT(' RI, DRA, DRCGS, RCGS, M(R), RM(R), T(R)',4E11.4,/,
     &        4E11.4)
         RQ=R(I)
         DEN(I)=DE(R(I))
         XELEC=DEL(I-1)
C     
C     UPDATE THE MEAN INTENSITY
C     
         IF(I.EQ.2) CALL SPEC(TE(I-1))
         IF(I.EQ.2.AND.N.GT.1) CALL SPEC(TE(I))
         DO J=JMIN,JJ
            JMEAN(J)=OLDFLU(J)
         ENDDO
         CALL RATE(1.d0)
         DO 6372 J=JMIN,JJ
            IF(I.NE.2)  FL(1,J)=FL(2,J)
 6372    CONTINUE
C     ***********************************************************
C     *****
C     NEW TEMPERATURES AND IONIZATION DEGREES FOR SLAB I
C     *****
C     ***********************************************************
         IF(N.GT.1.AND.I.EQ.2) THEN
            TOLD=TE(I)
            TOLDH=TE(I)
            TOLDCA=TE(I)
            TOLDO=TE(I)
            TOLDHE=TE(I)
            TOLDFE=TE(I)
            TOLDMG=TE(I)
         ENDIF
         IF(IBAL.EQ.1) THEN
            IF(I.EQ.IKBAL(KQ)) THEN
               WRITE(6,*)' SCAN OF HEATING AND COOLING FOR I= ',I
               KQ=KQ+1
               DTLOG=LOG10(TMAX/TMIN)/21.
               DO K=1,21
                  T2=TMIN*10.**(DTLOG*(K-1))
                  write(6,*)'rad 1',te
                  FQ=RAD(T2,XELEC,IFPOP)
               ENDDO
            ENDIF
         ENDIF
         IOIND=1
         ABSAA=0.
         DO K=3,14
            IF(IABU(K).EQ.1) THEN
               ABSDAB=ABS((AB(K)-ABOLD(K))/AB(K))
               ABSAA=DMAX1(ABSAA,ABSDAB)
            ENDIF
         ENDDO
         IF(N.GT.1) THEN
            A1=TEINT(I)*0.975
            B1=TEINT(I)*1.03
            XELEC=DEL(I)
         ENDIF
         IF(ABSAA.GT.0.5.AND.I.GE.3) THEN
            WRITE(6,*)' ABRUPT CHANGE IN COMPOSITION'
            WRITE(6,9721)(AB(IK),IK=1,2),AB(4),AB(5),AB(3),AB(10),AB(7),
     &           AB(9),AB(6),AB(12),AB(11),AB(8)
 9721       FORMAT(1X,'ABUND ',7E10.3)
            DO K=1,5
               write(6,*)'rad 2',te
               FA=RAD(B1IN,XELEC,IFPOP)
            ENDDO
            TP=B1IN
            if(te(i-1) < 500.) then
               tp=3.*te(i-1)
            endif
            DO K=1,100
               TP=0.95*TP
               FBB=RAD(TP,XELEC,IFPOP)
               IF(FA*FBB.LT.0.) GOTO 2211
            ENDDO
 2211       A1=0.95*TP
            B1=TP/(0.95*0.95)
         ENDIF
         IF(I.EQ.2) A1=A1IN
         IF(I.EQ.2) B1=B1IN
         IQWD=0
         if(i.lt.ii.and.i.gt.2) goto 700
         IF(I.NE.2) then
            do k=1,5
               FA=RAD(TE(I-1),XELEC,IFPOP)
            enddo
         endif
         do k=1,5
            FA=RAD(A1,XELEC,IFPOP)
         enddo
         T2=TE(I-1)
c No convergence in pop         
         IF(IFPOP.EQ.1) THEN 
 1923       T2=0.975*T2
            FA=RAD(T2,XELEC,IFPOP)
            IF(T2.LT.A1) GOTO 1924
            GOTO 1923
         ENDIF
 1924    CALL BIS(TE(I),XELEC,IFAIL,IFPOP,ICONV)
         IF(IFPOP.EQ.1.OR.ICONV.EQ.1) DELRM=DELRM/2.
         IF(INOUT.EQ.0) DELRM=-DELRM
         IF(IFPOP.EQ.1.OR.ICONV.EQ.1) then
            do j=jmin,jj
               write(6,9111)i,j,e1(j),FL(2,J),FPRIM,FD(1,J)
 9111          format('spec ',2i5,1pe12.3,10e12.3)
            enddo
            WRITE(6,*)' NEW DELRM = ',DELRM         
         endif
         IF(IFPOP.EQ.1.OR.ICONV.EQ.1)WRITE(6,*)' IFPOP,ICONV=',IFPOP,
     &        ICONV
         IF(IFPOP.EQ.1.OR.ICONV.EQ.1) NIPOP=NIPOP+1
         IF(IFPOP.NE.1.AND.ICONV.NE.1) GOTO 3646
         A1=A1OLD
         B1=B1OLD
         DO K=1,NL
            BOLH(K)=BQH(K)
            BOL(K)=BQ(K)
            BOLCA(K)=BQCA(K)
            BOLHE(K)=BQHE(K)
            BOLFE(K)=BQFE(K)
         ENDDO
         TOLD=TOLDQ
         TOLDH=TOLDHQ
         TOLDCA=TOLDCAQ
         TOLDO=TOLDOQ
         TOLDHE=TOLDHEQ
         TOLDFE=TOLDFEQ
         TOLDMG=TOLDMGQ
 3646    CONTINUE
         IF(NIPOP.GE.8) then
            write(6,*)'nipop > 8 ',nipop            
            goto 9365
         endif
         IF(IFPOP.EQ.1.OR.ICONV.EQ.1) GOTO 7832
         IF(IFAIL.EQ.1) IMAX=I-1
         IF(IFAIL.EQ.1) IQWD=IQWD+1
         IF(IFAIL.EQ.1) B1=2.*B1OLD
c     IF(IFAIL.EQ.1) A1=A1/2.
         IF(IFAIL.EQ.1) B1OLD=B1
         IF(IFAIL.EQ.1) WRITE(6,*)' NO CONVERGECE IN TEMPERATURE ',
     &        'FOR SHELL I=',I
         IF(IFAIL.EQ.1.AND.IQWD.LT.4) GOTO 1924
         IF(IFAIL.EQ.1) then
            write(6,*)'ifail=1'
            goto 9365
         endif
         IMAX=I
C     NEW LIMITS  
         B1=1.05*TE(I)
         A1=0.95*TE(I)
         B1=2.5*TE(I)
         A1=0.7*TE(I)
         IF(IFAIL.EQ.1) NIFA=NIFA+1
         IF(NIFA.GE.5) then
            WRITE(6,*)'stop at nifa ',nifa
            STOP
         endif
         TOLD1=TE(I)
         TOLD=TOLD1
         DO K=1,100
            XB(I,K)=XA(2,K)
         enddo
         DO K=1,NL
            BBOL(1,I,K)=BOLCA(K)
            BBOL(2,I,K)=BOL(K)
            BBOL(3,I,K)=BOLHE(K)
            BBOL(4,I,K)=BOLFE(K)
            BBOL(5,I,K)=BOLH(K)
            BQH(K)=BOLH(K)
            BQ(K)=BOL(K)
            BQCA(K)=BOLCA(K)
            BQHE(K)=BOLHE(K)
            BQFE(K)=BOLFE(K)
         ENDDO
         TOLDQ=TOLD
         TOLDHQ=TOLDH
         TOLDCAQ=TOLDCA
         TOLDOQ=TOLDO
         TOLDHEQ=TOLDHE
         TOLDFEQ=TOLDFE
         TOLDMGQ=TOLDMG
         DO J=JMIN,JJ
            OLDFLU(J)=JMEAN(J)
         ENDDO
C     
C     CALCULATE THE TOTAL FLUX IN THE LINES (ERG/SEC)
C
         K=0
         DO  IJU=2,30
            DO  IJL=1,IJU-1
               IF(AFE3(IJU,IJL).LE.0.) GOTO 8373
               K=K+1
               WLFEIII(K)=WLFE3(IJU,IJL)
               FEIII(K,I)=EMFE3(IJU,IJL)
 8373          CONTINUE
            enddo
         enddo

         DO  IJU=2,22
            DO  IJL=1,IJU-1
               IF(AFE4(IJU,IJL).LE.0.) GOTO 8374
               K=K+1
               WLFEIV(K)=WLFE4(IJU,IJL)
               FEIV(K,I)=EMFE4(IJU,IJL)
 8374          CONTINUE
            enddo
         enddo
         DO IJ=1,71
            DFBQ(IJ)=4.*PI*RCGS**2*DEN(I)**2.*DEL(I)*FB(IJ)
            FBQ(IJ)=DFBQ(IJ)*DR(I)+FBQ(IJ)
            DFBQ(IJ)=DFBQ(IJ)*DR(I)/ABS(DELRM)
            FBI(IJ,I)=FB(IJ)
            WEQH(IJ)=4.*PI*RCGS**2*DEN(I)**2.*WEH(IJ)*DR(I)+WEQH(IJ)
            WEQH(IJ)=WEQH(IJ)+WEH(IJ)
            WEQ(IJ)=4.*PI*RCGS**2*W(IJ)*DEN(I)**2*DR(I)+WEQ(IJ)
         enddo
         DO IJ=1,100
            DCIT(IJ)=4.*PI*RCGS**2*DEN(I)**2*DEL(I)*CIN(IJ)
            CITOT(IJ)=CITOT(IJ)+DCIT(IJ)*DR(I)
            DCIT(IJ)=DCIT(IJ)*DR(I)/ABS(DELRM)
         enddo
C     FE II LINES

         DO IU=2,NFEII
            DO IL=1,IU-1
               DWFEII(IU,IL)=4.*PI*RCGS**2*DEN(I)**2*WOBFE(IU,IL)
               WFEII(IU,IL)=WFEII(IU,IL)+DWFEII(IU,IL)*DR(I)
            enddo
         enddo
C     FE III LINES
         DO IU=2,NFEIII
            DO IL=1,IU-1
               DWFEIII(IU,IL)=4.*PI*RCGS**2*DEN(I)**2*DEL(I)*
     &              EMFE3(IU,IL)
               EMFE3(IU,IL)=0.
               WFEIII(IU,IL)=WFEIII(IU,IL)+DWFEIII(IU,IL)*DR(I)
            enddo
         enddo
C     FE IV LINES
         DO IU=2,NFEIV
            DO IL=1,IU-1
               DWFEIV(IU,IL)=4.*PI*RCGS**2*DEN(I)**2*DEL(I)*EMFE4(IU,IL)
               EMFE4(IU,IL)=0.
               WFEIV(IU,IL)=WFEIV(IU,IL)+DWFEIV(IU,IL)*DR(I)
            enddo
         enddo

C     TOTAL HEATING AND COOLING
         DCCOOL=4.*PI*RCGS**2*DEN(I)**2*COOL*DEL(I)*DR(I)
         DCHEAT=4.*PI*RCGS**2*DEN(I)**2*HEAT*DR(I)
         DGAMM=4.*PI*RCGS**2*GAHE*DR(I)
         HEATTOT=HEATTOT+DCHEAT
         COOLTOT=COOLTOT+DCCOOL
         GAMTOT=GAMTOT+DGAMM
         DCHEATQ=DCHEAT/ABS(DELRM)
c     O I rec.
         CALL OXREC(OXR,TE(I))
         DO IJ=1,8
            OXR(IJ)=AB(3)*XA(2,15)*DEL(I)*OXR(IJ)
            WROX(IJ)=4.*PI*RCGS**2*DEN(I)**2*DR(I)*OXR(IJ)+WROX(IJ)
         enddo
         FBI(51,I)=OXR(3)/DEL(I)
         FBI(52,I)=OXR(4)/DEL(I)
         FBI(53,I)=OXR(7)/DEL(I)
         FBI(54,I)=OXR(8)/DEL(I)
c         REHE1=REHE1+4.*PI*RCGS**2.*DR(I)*DEN(I)**2.*REHE21
         REN1=REN1+4.*PI*RCGS**2.*DR(I)*DEN(I)**2.*REN41
         REC1=REC1+4.*PI*RCGS**2.*DR(I)*DEN(I)**2.*REC31
         RECOX1=4.*PI*RCGS**2*DEN(I)**2*DR(I)*ALO(2)*13.6*ELCH*
     &        AB(3)*XB(I,15)*DEL(I)
         RECOX2=4.*PI*RCGS**2*DEN(I)**2*DR(I)*ALO(3)*35.1*ELCH*
     &        AB(3)*XB(I,16)*DEL(I)
         HEATT=4.*PI*RCGS**2*DR(I)*GAELH
         TOTI1=4.*PI*RCGS**2*DR(I)*AB(3)*XB(I,14)*CISEC(2)*13.6*ELCH
         TOTI2=4.*PI*RCGS**2*DR(I)*AB(3)*XB(I,15)*CISEC(3)*35.1*ELCH
c calculate total line lum.. 
         call emiss_calc(npri,rcgs,drcgs,den(i),del(i))
c         call sort(1)
C     TOTAL ENERGY
C     CALCULATE FLUX
         DO J=JMIN,JJ
            CALL RADTRANPULS(1,J,TAU(I,J),S)
         ENDDO
         V=0.
         VV=0.
         DO J=47,JJ
            VV=VV+(FH0(J)+FHD(J))*(E(J+1)-E(J))
         ENDDO
         DO  J=2,46
            V=V+(FH0(J)+FHD(J))*(E(J+1)-E(J))
         ENDDO
         V=1.E30*4.*GEOM*PI**2.*(R15*R(1))**2.*V
         VV=1.E30*4.*GEOM*PI*PI*(R15*R(1))**2.*VV
         VVV=0.
         VVV1=0.
         DO J=JMIN,JJ
            VVV=VVV+(FH0(J)+FHD(J))*(E(J+1)-E(J))
            VVV1=VVV1+FH0(J)*(E(J+1)-E(J))
         ENDDO
         ERADTOT=1.E30*4.*GEOM*PI*PI*(R15*R(1))**2.*VVV
         ERADTOT1=1.E30*4.*GEOM*PI*PI*(R15*R(1))**2.*VVV1
         DERAD=ERADOLD-(V+VV)
         ERADOLD=V+VV
         WRITE(6,525)
         DO IK=1,14
            ABUN(I,IK)=9.999
            IF(AB(IK).GT.0)  ABUN(I,IK)=-LOG10(AB(IK))
         enddo
c$$$C     ***********************************************************
c$$$C     *****
c$$$C     CHECK RECOMBINATION TIME SCALE
c$$$C     *****
c$$$C     ***********************************************************
         IRECWARN=0
         TIME=TDAYS*8.64E4
         TREC(1)=1./AL2(1)
         TREC(2)=1./AL2(2)
         TREC(3)=1./AL2(3)
         DO K=2,9
            TREC(K+3)=1./ALO(K)
         ENDDO
         DO K=2,9
            TREC(K+11)=1./ALsi(K)
         ENDDO
         DO K=1,20
            TREC(K)=TREC(K)/(DEL(I)*DEN(I))
            IF(TREC(K).GT.TIME/2.) IRECWARN=1
            IF(TREC(K).GT.TIME/2.) IRECWARNQ=1
         ENDDO

 9323    format(1pe12.4,5e12.4)
         IF(N.LE.NMAX.AND.IOSHELL.EQ.0) GOTO 543
         RCGS=1.E15*R(I)*R15
         VELI=RCGS*VEXP/RMAXCGS
         IF(I.EQ.2) VELIMAX=VELI
         DRSCGS=1.E15*R15*(R(I)-R(1))
         IF(ISTAT.NE.1) THEN
            WRITE(6,900)I,VELI/1.E5,DRSCGS,DEN(I),DEL(I),TE(I)
         ENDIF
         IF(ISTAT.EQ.1.AND.IREV.NE.1) THEN
            WRITE(6,9910)I,RCGS,DEN(I),DEL(I),TE(I)
         ENDIF
         IF(IREV.EQ.1) WRITE(6,9910)I,DRSCGS,DEN(I),DEL(I),TE(I)
         WRITE(6,9721)(AB(IK),IK=1,2),AB(4),AB(5),AB(3),AB(13),AB(10),
     &        AB(7),AB(9),AB(6),AB(12),AB(14),AB(11),AB(8)
         TAUELE=TAUELE+DENQ*DR(I)*DEL(I)*0.665E-24
         WRITE(6,9283)RMASS,COLUMN,TAUR,DENQ,AMEAN
c     WRITE(6,9283)RM(I),COLUMN,TAUR,DENQ,AMEAN
         DERADDV=DERAD/(4.*PI*RCGS**2*DR(I)*DEN(I)**2)
         WRITE(6,903) ERADTOT,ERADTOT1,V,VV,vion,v13p6_1e4,DERAD,DERADDV,
     &     (CIT+FBT),R(I)*R15*1.e15
c         WRITE(6,903) ERADTOT,ERADTOT1,V,VV,DERAD,DERADDV,
c     &        (CIT+FBT),R(I)*R15
         RSUN=MTOT/SOLMA
c     
         etau1=0.
         DO J=JJ,JMIN,-1
            IF(TAUTOT(I,J).GT.1.) GOTO 2876
         ENDDO
 2876    ETAU1=E1(J+1)
         taumax=0.
         DO J=jmin,JJ
            if(tautot(i,j)>taumax) then
               taumax=tautot(i,j)
               jtaum=j
            endif
         ENDDO
         WRITE(6,3966)TAUTOT(I,2),TAUTOT(I,23),
     &        TAUTOT(I,47),tautot(i,63),tautot(i,87),tautot(i,99),
     &        ETAU1,taumax,e1(jtaum)
         ZPR=BOL(1)*PHOT(1)-REC(1)
         ZHE=BOL(1)*PHE(1)-RECNET(1)
         IF(IDEP.EQ.1) THEN
            WRITE(6,2758)(BOL(KK),KK=1,10)
            WRITE(6,2759)(BOLCA(KK),KK=1,4)
         ELSE
            DO L=1,5
               SUM=0.
               DO KK=1,KKMAX(L)
                  SUM=SUM+XN(L,KK)
               ENDDO
               XN(L,KKMAX(L)+1)=1.-SUM
            ENDDO
            WRITE(6,2758)(XN(2,KK),KK=1,10)
            WRITE(6,2759)(XN(1,KK),KK=1,4)
         ENDIF
         CMH=HEAT-XEL*COOL
         WRITE(6,904)XEL,ZEL,HEAT,COOL,CMH
         DO IH=2,NH
            IPM=IH-1
 9922       FORMAT(' EMN ',6E11.4)
C     LINE COOLING RATE
            DO IP=1,IPM
               EMIS(IH,IP)=DEN(I)**2*EMIS(IH,IP)
 4382          EMTOT(IH,IP)=EMTOT(IH,IP)+DR(I)*EMIS(IH,IP)
            enddo
         enddo
 9374    FORMAT(1PE15.7,6E15.7)
         NPR=NPRINT
         IF(NUM.LT.NPR) GOTO 543
         WRITE(6,901)RCGS
         write(6,9573)(xion(i,3,l),l=1,7)
         write(6,9574)(xion(i,4,l),l=1,8)
         write(6,9575)(xion(i,5,l),l=1,9)
         write(6,9584)(xion(i,6,l),l=1,11)
         write(6,9582)(xion(i,7,l),l=1,12)
         write(6,9578)(xion(i,8,l),l=1,13)
         write(6,9580)(xion(i,9,l),l=1,14)
         write(6,9576)(xion(i,10,l),l=1,15)
         write(6,9583)(xion(i,11,l),l=1,17)
         write(6,9585)(xion(i,12,l),l=1,19)
         write(6,9581)(xion(i,13,l),l=1,21)
         WRITE(6,9579)(xion(i,14,l),l=1,27)


         write(59,9291)I,VELI/1.E5,DRSCGS,DEN(I),DEL(I),TE(I),
     &        (xion(i,5,ion),ion=1,8),(xion(i,10,ion),ion=1,14),
     &        (xion(i,11,ion),ion=1,16),(xion(i,12,ion),ion=1,18),
     &        (xion(i,13,ion),ion=1,5),(xion(i,14,ion),ion=1,14),
     &        (xion(i,6,ion),ion=1,10)
c O 6-13, Si 14-27, S 28-43, Ar 44-61, Ca 62-66, Fe 67-80, Ne 81-90         
 9291    format(i4,f12.4,1pe13.5,3e11.3,111e11.3)
         IF(IRECWARN.EQ.1) THEN
            WRITE(6,*)' Recombination time scale is long!!'
            WRITE(6,9028)(TREC(K),K=1,12)
 9028       FORMAT(5E12.5)            
         ENDIF

         NUM=0
 543     NUM=NUM+1
         M=M+1
         IF(IAGN.EQ.1) THEN
            DIST=DRSCGS
         ELSE
            DIST=RCGS
         ENDIF
c         WRITE(17,9374)DIST,RMASS,(DCIT(IJ),IJ=1,100),DCHEATQ
c         WRITE(17,7467)(DFBQ(JS),JS=1,45)
c         WRITE(21,9374)DIST,RMASS,(DCIT(IJ),IJ=1,100),DCHEATQ
c         WRITE(21,7467)(DFBQ(JS),JS=1,45)
c Total line cooling
         coolt=0.
         do iel=1,14
            do ion=1,27
               coolt=coolt+line_cool(iel,ion)
               if(abs(line_cool(iel,ion)) > 0.) then
c                  write(6,9241)iel,ion,line_cool(iel,ion),coolt
c 9241             format('line cooling ',2i5,1pe12.3,e12.3)
               endif
            enddo
         enddo
C     
         RSUN=MTOT/SOLMA
         IF(rmass.ge.revmass.and.n.eq.1.and.irev.eq.1) then
            write(6,*)'rmass>revmass'
            goto 9365
         endif
         IF(veli.ge.2.e8.and.n.eq.1) then
            write(6,*)'veli>2000km/s'
            goto 9365
         endif
         IF(i.ge.iqmax.and.n.eq.1) then
            write(6,*)'i>iqmax'
            goto 9365
         endif
         IF(RMASS.GT.RSUN) then
            write(6,*)'rmass>rsun'
            goto 9365
         endif
         IF(TE(I).LT.100.) GOTO 3288
c         IF(TE(I).LT.50.) GOTO 3288
C     IF(NMAX.NE.0) GOTO 700
C     STOP
         IF(TE(I).GT.9.E3) GOTO 700
         MQW=MQW+1
C     IF(MQW.GE.MQMAX) GOTO 9365
 700  CONTINUE 
C     ***********************************************************
C     *****
C     SOLVE THE EQUATION OF TRANSFER FOR THE DIFFUSE EMISSION
C     *****
C     ***********************************************************
 9365 IF(ISPH.EQ.0) CALL DIFF(IMAX,FL0)
      DO I=1,IMAX
         TEINT(I)=TE(I)
      enddo
2727  FORMAT(' TI ',5E10.4)
      IF(N.EQ.NMAX.OR.IMAX.EQ.2) GOTO 7327
      DO J=JMIN,JJ
         CALL SPHTRP(IMAX,J,R15,TEINT,RM,IMAXP)
      enddo
      write(6,*)'sph',(ts2-ts1)
      WRITE(6,2727)TEINT
7327  MAX=IMAXP
C
C     SAVE IMAX, RM, AND DIFFUSE FLUXES FD  FOR THIS ITERATION ON UNIT
C        14 = FILE "UTIN"
C
      IF(IMAX.eq.2) THEN
c      IF(IMAX.NE.2) THEN
      WRITE(6,*)' FLUXES WRITTEN ON FILE UTIN'
      REWIND 14
      WRITE(14)N
      WRITE(14)IMAXP
      WRITE(14)RM
      WRITE(14)FD
      WRITE(6,9282)(RM(IK),IK=1,IMAX)
      ENDIF
9282  FORMAT(' RM ',5E11.4)
C     PRINT OUTGOING FLUX ON FILE 31
c      DO J=JMIN,JJ
c        WRITE(31,9128)J,E1(J),E(J),EDDFLUX(J)
c      ENDDO
9128  FORMAT(I5,3E13.5)
C     SAVE TOTAL COLUMN DENSITY (DIVIDED BY SQRT(T) AND A) FOR ESCAPE PROB.
      DO K=1,100
        COLTOTT(K)=COLTOT(K)
      ENDDO
      DO K=1,5
        DO I=1,NL
            COLHVTT(K,I)=COLHVT(K,I)
        ENDDO
      ENDDO
C     ***********************************************************
C     *****
C     CALCULATE LINE PROFILES
C     *****
C     ***********************************************************
      R1CGS=R15*1.E15
      CALL ATDATFE
      NI=1
      DO 221 I=2,NFEII
      DO 222 J=1,I-1
      IF (AQ(I,J).GT.0.)THEN
      WLI(81+NI)=WL(I,J)
      NI=NI+1
      ENDIF
 222  CONTINUE
 221  CONTINUE      

C     ***********************************************************
C     *****
C     CALCULATE THE TOTAL ENERGY OF THE DIFFUSE EMISSION ABOVE 13.6
C     EV AND BELOW.
C     *****
C     ***********************************************************
      V=0.
      VV=0.
      DO 8187 J=JMIN,1
 8187 V=V+FL0(1,J)*(E(J+1)-E(J))
      DO 8178 J=2,JJ
 8178 VV=VV+FL0(1,J)*(E(J+1)-E(J))
      V=(4.*PI)**2.*1.E30*R15**2.*V
      VV=(4.*PI)**2.*1.E30*R15**2.*VV
      WX=0.
      WW=0.
      DO 8185 J=JMIN,1
 8185 WX=WX+FL0(2,J)*(E(J+1)-E(J))
      DO 8186 J=2,JJ
 8186 WW=WW+FL0(2,J)*(E(J+1)-E(J))
      WX=(4.*PI)**2.*1.E30*R15**2.*WX
      WW=(4.*PI)**2.*1.E30*R15**2.*WW
      EDIL=0.
      EDIH=0.
      DO 8182 I=2,IMAX
      DO 8183 J=JMIN,1
 8183 EDIL=16.*PI**2.*R15**2.*DR(I)*DEN(I)**2.*EM(I,J)*
     &(E(J+1)-E(J))+EDIL
      DO 8184 J=2,JJ
 8184 EDIH=16.*PI**2.*R15**2.*DR(I)*DEN(I)**2.*EM(I,J)*
     &(E(J+1)-E(J))+EDIH
 8182  CONTINUE
 3288 WRITE(6,699)
      WRITE(6,*)'to final '
 9433 CONTINUE
C     SORT LINES IN LUMINOSITY
      write(6,*)' to sort '
      call sort(1)

      WRITE(6,525)
      WRITE(6,9117)SUMCO,SUMRE
 9117 FORMAT(' TOTAL COOLING =',1PE11.4,' TOTAL REC. EM. ',E11.4)
      WRITE(6,525)
      WRITE(6,9567)COOLTOT,HEATTOT,GAMTOT
 9567 FORMAT(' TOTAL COOLING =',1PE11.4,' TOTAL HEATING ',E11.4,
     &      ' TOTAL INPUT =',E11.4)
C     ***********************************************************
C     *****
C     SHOULD THE RADIATION FIELD BE ITERATED ONCE MORE?
C     *****
C     ***********************************************************

C     SPECTRUM
C     DO 1928 J=JMIN,JJ
C1928  F0(J)=F0(J)/R(IMAX)**2
      WRITE(6,*)'INPUT SPECTRUM'
      WRITE(6,9420)(F0(J),J=JMIN,JJ)
      WRITE(6,*)' '
      WRITE(6,*)'TRANSMITTED SPECTRUM'
      WRITE(6,9420)(FL(2,J),J=JMIN,JJ)
      WRITE(6,*)'OPTICAL DEPTHS'
      WRITE(6,9420)(TAUTOT(IMAX,J),J=JMIN,JJ)
9100  FORMAT(1X,'*******************************************************
     &******************')
9899  FORMAT(A)
2758  FORMAT(1X,'N(O)=',5E10.3)
2759  FORMAT(1X,'N(CA)=',5E10.3)
9573  FORMAT(' C     ',1PE10.2,9E10.2)
9574  FORMAT(' N     ',1PE10.2,9E10.2)
9575  FORMAT(' O     ',1PE10.2,9E10.2)
9576  FORMAT(' SI    ',1PE10.2,9E10.2)
9577  FORMAT(' SI    ',1PE10.2,9E10.2)
9578  FORMAT(' MG    ',1PE10.2,9E10.2)
9579  FORMAT(' FE    ',1PE10.2,9E10.2)
9580  FORMAT(' AL    ',1PE10.2,9E10.2)
9581  FORMAT(' CA    ',1PE10.2,9E10.2)
9582  FORMAT(' NA    ',1PE10.2,9E10.2)
9583  FORMAT(' S     ',1PE10.2,9E10.2)
 9584 FORMAT(' NE    ',1PE10.2,9E10.2)
 9585 FORMAT(' AR    ',1PE10.2,9E10.2)

 
  900 FORMAT(1X,'I=',I3,1X,'V=',F9.3,1X,'R(i)-R(1)=',1PE9.3,1x,
     &'N=',E9.3,1X,'X=',E9.3,1X,'T=',0PF9.1)
 9910 FORMAT(1X,'I=',I3,1X,'DR=',1PE11.5,' CM',1X,'DENS=',E10.4,1X,
     &      'EL=',E10.4,1X,'TEMP=',E10.4)
 9283 FORMAT(1X,'M(R)=',1PE10.3,1X,'N=',E10.3,1X,'TAU=',E10.3
     &,' DEN=',E10.3,' <A>=',0PF6.2)

 904  FORMAT(1X,'COOLING=',8E12.5)
912   FORMAT(1X,'FLUX AT INNER BOUNDARY BELOW 13.6 EV=',1PE12.5,
     &' ABOVE 13.6 EV',E12.5)
914   FORMAT(1X,'FLUX AT SURFACE BELOW 13.6 EV=',1PE12.5,
     &' ABOVE 13.6 EV',E12.5)
913   FORMAT(1X,'DIFFUSE EMISSION BELOW 13.6 EV=',1PE12.5,
     &' ABOVE 13.6 EV',E12.5)
 3966 FORMAT(1X,'T(13.6)=',1PE12.4,1X,'T(24.5)=',E12.4,1X,
     &     'T(54.0)=',E12.4,1X,'T(100)=',E12.4,1X,'T(500)=',E12.4,1X,
     &     'T(1keV)=',E12.4,1X,'E(TAU=1)',E12.4,'TAUMAX=',e12.4,
     &     ' E(TAUMAX)=',e12.3)
 1008 FORMAT(1X,'HE II (1640) =',1pE10.3,'  HE II (4686)  =',E10.3)
 1098 FORMAT(1X,'He II (304)  =',1pe10.3,'  HE II (Two-g) =',E10.3)
c     &'  C III (1909) REC.=',E10.3,'  N IV (1486) REC.=',E10.3)
7468  FORMAT(' REC=',8E11.4)
7469  FORMAT(' RES=',8E11.4)
 7467 FORMAT(1X,7E12.4)
  901 FORMAT(1X,'RADIUS= ',1PE14.7)
 9420 FORMAT(1X,5E12.5)
  699 FORMAT(1X,'FINAL RESULTS')

      RETURN
      END




      SUBROUTINE IONLABEL(LAB1)
      CHARACTER*8 LAB(200),LAB1(200)
c               1234567800012345678000123456780001234567800012345678
      DATA LAB/'H I     ','HE I    ','HE II   ','O VI    ','O VII   '
     &        ,'O VIII  ','O V     ','C III   ','C IV    ','C V     '
     &        ,'C VI    ','C I     ','C II    ','O I     ','O II    '
     &        ,'O III   ','O IV    ','N I     ','N II    ','N III   '
     &        ,'N IV    ','N V     ','N VI    ','N VII   ','SI I    '
     &        ,'SI II   ','SI III  ','SI IV   ','SI V    ','SI VI   '
     &        ,'SI VII  ','SI VIII ','SI IX   ','SI X    ','SI XI   '
     &        ,'SI XII  ','SI XIII ','SI XIV  ','MG I    ','MG II   '
     &        ,'MG III  ','FE I    ','FE II   ','FE III  ','FE IV   '
     &        ,'AL I    ','AL II   ','AL III  ','AL IV   ','CA I    '
     &        ,'CA II   ','CA III  ','NA I    ','NA II   ','NA III  '
     &        ,'S I     ','S II    ','S III   ','NE I    ','NE II   '
     &        ,'AR I    ','AR II   ',38*'        ','BALMER  '
     &        ,'PASCHEN ','BRACKETT','PFUND   ','O I 4   ','O I 5   '
     &        ,'O I 6   ','O I 7   ','O I 8   ','O I 9   ','FF      '
     &        ,89*' '/
      DO K=1,111
            LAB1(K)=LAB(K)
      ENDDO
      RETURN
      END

      SUBROUTINE WLFORB(WLF)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION WLF(71),WLF1(71)
C     WAVELENGTHS OF ALL FORBIDDEN LINES
      DATA WLF1/5007.,2321.,4363.,3726.,2470.,7320.,6548.,3063.,
     &         5755.,6300.,2964.,5581.,3468.,5201.,10406.,0.,630000.,
     &         0.,1455600.,9818.,4619.,8729.,6718.,4071.,10320.,
     &         252460.,563060.,1297016.,684932.,16360.,10995.,
     &         3704000.,6091000.,3722.,6312.,9069.,1815.,3342.,
     &         3869.,1602.,2423.,4714.,1575.,2975.,3346.,     
     &     5*0.,11287.,7002.,9264.,6157.,89913.84,218314.71,7137.75,
     &     7753.24,3006.10,3110.07,3155.02,
     &     4741.49,4712.69,2868.99,2854.48,775794.63,564652.33,7264.74,
     &     7333.42,7172.46,7239.39/
      DO K=1,71
            WLF(K)=WLF1(K)
      ENDDO
      RETURN
      END


      subroutine emiss_calc(npri,rcgs,drcgs,deni,xel)
      implicit real*8(a-h,o-z)
      common/line_em/cinx(14,26,401),taulinex(14,26,401),wlix(14,26,401)
     &     ,ilabcfx(14,26,401)
      common/kmax/kmax(14,27)
c     number of ionization stages on ordered scale 1-14      
      integer nionstage
      common/num_ion_stages/nionstage(14)
      COMMON/IND/I
      real*8 tot_line_lum
      common/linelum/tot_line_lum(14,27,401)      
      data pi/3.14159/
      totn=4.*pi*rcgs**2*drcgs*deni**2*xel
c      write(74,*)'I, rcgs, drcgs ',i,rcgs,drcgs
c      write(75,*)i,rcgs,drcgs
      n=1
      npri=npri+1
      if(npri==1) then
         do iel=1,14
            do ion=1,nionstage(iel)
c     do k=1,kmax(iel,ion)
               do k=1,401
                  tot_line_lum(iel,ion,k)=totn*cinx(iel,ion,k) +
     &                 tot_line_lum(iel,ion,k)
c     if(iel==14.and.ion<=3.and.k<40) then
                  if(tot_line_lum(iel,ion,k).ne.0.d0) then
                     n=n+1
c                     write(74,91)iel,ion,k,wlix(iel,ion,k),totn*cinx(iel,ion,k),
c     &                    tot_line_lum(iel,ion,k)
c                     write(75,*)iel,ion,k,wlix(iel,ion,k),totn*cinx(iel,ion,k),
c     &                    tot_line_lum(iel,ion,k)
 91                  format(3i5,f12.1,1pe12.3,10e12.3)
                  endif
               enddo
            enddo
         enddo
         npri=0         

c         call sort(0)
      endif
      return
      end


      subroutine sort(iend)
      implicit real*8 (a-h,o-z)
      common/line_em/cinx(14,26,401),taulinex(14,26,401),wlix(14,26,401)
     &     ,ilabcfx(14,26,401)
c      common/kmax/kmax(14,27)
      common/kmaxpop/kmaxp(14,27)
c     number of ionization stages on ordered scale 1-14      
      integer done(14,27,401)
      integer nionstage
      common/num_ion_stages/nionstage(14)
      real*8 tot_line_lum
      common/linelum/tot_line_lum(14,27,401)      
      data pi/3.14159/
      character*2 elid(14)
      data elid/'H ','He','C ','N ','O ','Ne','Na','Mg','Al','Si','S ','Ar','Ca','Fe'/
      character*5 ionid(14)
      data ionid/'I   ','II   ','III  ','IV   ','V    ','VI   ','VII  ',
     &     'VIII ','IX   ','X    ','XI   ','XII  ','XIII ','XIV  '/

      write(6,*)' 500 Strongest lines'
      tot_lum=0.
      do iel=1,14
         do ion=1,27
            do k=1,401
               done(iel,ion,k)=0
               tot_lum=tot_lum + tot_line_lum(iel,ion,k)
            enddo
         enddo
      enddo
      write(6,9276)tot_lum
      if(iend==1) then
         write(76,9276)tot_lum
      endif
 9276 format('Total line luminosity ',1pe12.3)
      CIMAX=1.D100
      ar2lum=tot_line_lum(12,2,1)
      nstr=500
      DO KL=1,nstr
         CMAX=0.
         do iel=1,14
            do ion=1,nionstage(iel)
               do k=1,kmaxp(iel,ion)
c                  if(iel.eq.-14.and.ion==2)write(6,*)'Fe II ', 
c     &                 iel,ion,k,done(iel,ion,k),wlix(iel,ion,k),
c     &                 tot_line_lum(iel,ion,k),cmax
                  if(tot_line_lum(iel,ion,k) > cmax .and. 
     &                 done(iel,ion,k).eq.0) then
                     cmax=tot_line_lum(iel,ion,k)
                     ielx=iel
                     ionx=ion
                     kx=k
                     if(iel.eq.-14.and.ion==2)write(6,*)'Fe II max ', 
     &                    ielx,ionx,kx,wlix(iel,ion,k)/1.e4,cmax
                  endif
               enddo
            enddo
         enddo
         done(ielx,ionx,kx)=1
         write(6,9)elid(ielx),ionid(ionx),wlix(ielx,ionx,kx),cmax,cmax/ar2lum
c         write(74,9)elid(ielx),ionid(ionx),wlix(ielx,ionx,kx),cmax,cmax/ar2lum
         if(iend==1) then
            write(76,9)elid(ielx),ionid(ionx),wlix(ielx,ionx,kx),cmax,cmax/ar2lum
         endif
 9       format(2a6,f12.2,1pe12.3,e12.3)
      enddo
      return
      end
      
