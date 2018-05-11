      SUBROUTINE HATCH
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
cale: egs4comm.for per pathname
      include 'egs4comm.for'
      character*80 Matfilename

      character*80 totMatfilename
cale

c
      CHARACTER*4 MBUF(72),MDLABL(8)
      real*4 ACD ,  ADEV ,  ASD ,  COST ,  CTHET ,  DEL ,  DFACT ,  DFAC
     *TI,  DUNITO,  DUNITR,  FNSSS ,  P ,  PZNORM,  RDEV ,  S2C2 ,  S2C2
     *MN,  S2C2MX,  SINT ,  SX ,  SXX ,  SXY ,   SY ,   WID ,  XS ,  XS0
     * ,  XS1 ,  XSI ,  WSS ,  YS ,  ZEROS(3)
      integer*4 I ,  I1ST ,  IB ,  ID ,  IE ,  IL ,  IM ,  IRAYL ,  IRN
     *,  ISTEST,  ISUB ,  ISS ,  IZ ,   IZZ ,  J ,  JR ,  LCTHET,  LMDL
     *,  LMDN ,  LTHETA,  MD ,  MXSINC,  NCMFP ,   NEKE ,   NGE ,   NGRI
     *M ,  NISUB ,  NLEKE ,    NM ,  NRANGE,    NRNA ,  NSEKE ,   NSGE ,
     *   NSINSS,  LOK(21)
      DATA MDLABL/' ','M','E','D','I','U','M','='/,LMDL/8/,LMDN/24/,DUNI        
     *TO/1./
      DATA I1ST/1/,NSINSS/37/,MXSINC/1002/,ISTEST/0/,NRNA/1000/
3240  FORMAT(1X,14I5)
3250  FORMAT(1X,1PE14.5,4E14.5)
3260  FORMAT(72A1)
      IF ((I1ST.NE.0)) THEN
        I1ST=0
          DO 3271 J=1,numgeom
          IF ((SMAXIR(J).LE.0.0)) THEN
            SMAXIR(J)=1E10
          END IF
          IF ((ESTEPR(J).LE.0.0)) THEN
            ESTEPR(J)=1.0
          END IF
3271    CONTINUE
3272    CONTINUE
        PRM=RM
        PRMT2=2.D0*PRM
        PZERO=0.0D0
        NISUB=MXSINC-2
        FNSSS=NSINSS
        WID=PI5D2/FLOAT(NISUB)
        WSS=WID/(FNSSS-1.0)
        ZEROS(1)=0.
        ZEROS(2)=PI
        ZEROS(3)=TWOPI
          DO 3281 ISUB=1,MXSINC
          SX=0.
          SY=0.
          SXX=0.
          SXY=0.
          XS0=WID*FLOAT(ISUB-2)
          XS1=XS0+WID
          IZ=0
            DO 3291 IZZ=1,3
            IF (((XS0.LE.ZEROS(IZZ)).AND.(ZEROS(IZZ).LE.XS1))) THEN
              IZ=IZZ
              GO TO3292
            END IF
3291      CONTINUE
3292      CONTINUE
          IF ((IZ.EQ.0)) THEN
            XSI=XS0
          ELSE
            XSI=ZEROS(IZ)
          END IF
            DO 3301 ISS=1,NSINSS
            XS=WID*FLOAT(ISUB-2)+WSS*FLOAT(ISS-1)-XSI
            YS=SIN(XS+XSI)
            SX=SX+XS
            SY=SY+YS
            SXX=SXX+XS*XS
            SXY=SXY+XS*YS
3301      CONTINUE
3302      CONTINUE
          IF ((IZ.NE.0)) THEN
            SIN1(ISUB)=SXY/SXX
            SIN0(ISUB)=-SIN1(ISUB)*XSI
          ELSE
            DEL=FNSSS*SXX-SX*SX
            SIN1(ISUB)=(FNSSS*SXY-SY*SX)/DEL
            SIN0(ISUB)=(SY*SXX-SX*SXY)/DEL - SIN1(ISUB)*XSI
          END IF
3281    CONTINUE
3282    CONTINUE
        SINC0=2.0
        SINC1=1.0/WID
        IF ((ISTEST.NE.0)) THEN
          ADEV=0.
          RDEV=0.
          S2C2MN=10.
          S2C2MX=0.
            DO 3311 ISUB=1,NISUB
              DO 3321 ISS=1,NSINSS
              THETA=WID*FLOAT(ISUB-1)+WSS*FLOAT(ISS-1)
              CTHET=PI5D2-THETA
              SINTHE=sin(THETA)
              COSTHE=sin(CTHET)
              SINT=SIN(THETA)
              COST=COS(THETA)
              ASD=ABS(SINTHE-SINT)
              ACD=ABS(COSTHE-COST)
              ADEV=max(ADEV,ASD,ACD)
              IF((SINT.NE.0.0))RDEV=max(RDEV,ASD/ABS(SINT))
              IF((COST.NE.0.0))RDEV=max(RDEV,ACD/ABS(COST))
              S2C2=SINTHE**2+COSTHE**2
              S2C2MN=min(S2C2MN,S2C2)
              S2C2MX=max(S2C2MX,S2C2)
              IF ((ISUB.LT.11)) THEN
                WRITE(6,3330)THETA,SINTHE,SINT,COSTHE,COST
3330            FORMAT(1PE20.7,4E20.7)
              END IF
3321        CONTINUE
3322        CONTINUE
3311      CONTINUE
3312      CONTINUE
          WRITE(6,3340)MXSINC,NSINSS
3340      FORMAT(' SINE TESTS,MXSINC,NSINSS=',2I5)                              
          WRITE(6,3350)ADEV,RDEV,S2C2MN,S2C2MX
3350      FORMAT(' ADEV,RDEV,S2C2(MN,MX) =',1PE16.8,3E16.8)                     
          ADEV=0.
          RDEV=0.
          S2C2MN=10.
          S2C2MX=0.
            DO 3361 IRN=1,NRNA
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            THETA = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            THETA=THETA*PI5D2
            CTHET=PI5D2-THETA
            SINTHE=sin(THETA)
            COSTHE=sin(CTHET)
            SINT=SIN(THETA)
            COST=COS(THETA)
            ASD=ABS(SINTHE-SINT)
            ACD=ABS(COSTHE-COST)
            ADEV=max(ADEV,ASD,ACD)
            IF((SINT.NE.0.0))RDEV=max(RDEV,ASD/ABS(SINT))
            IF((COST.NE.0.0))RDEV=max(RDEV,ACD/ABS(COST))
            S2C2=SINTHE**2+COSTHE**2
            S2C2MN=min(S2C2MN,S2C2)
            S2C2MX=max(S2C2MX,S2C2)
3361      CONTINUE
3362      CONTINUE
          WRITE(6,3370)NRNA
3370      FORMAT(' TEST AT ',I7,' RANDOM ANGLES IN (0,5*PI/2)')                 
          WRITE(6,3380)ADEV,RDEV,S2C2MN,S2C2MX
3380      FORMAT(' ADEV,RDEV,S2C2(MN,MX) =',1PE16.8,3E16.8)                     
        END IF
        P=1.
          DO 3391 I=1,50
          PWR2I(I)=P
          P=P/2.
3391    CONTINUE
3392    CONTINUE
      END IF
        DO 3401 J=1,NMED
3410    CONTINUE
          DO 3411 I=1,numgeom
          IF ((IRAYLR(I).EQ.1.AND.MED(I).EQ.J)) THEN
            IRAYLM(J)=1
            GO TO 3412
          END IF
3411    CONTINUE
3412    CONTINUE
3401  CONTINUE
3402  CONTINUE
      REWIND KMPI
cale NB: fort.12 e' il file dei materiali by PEGS4

      Matfilename='Mat900Kev.dat'

      totMatfilename=pathname(1:pathlength) // Matfilename
      write(*,*) 'cale',pathname, totMatfilename
      OPEN (UNIT=KMPI,FILE=totMatfilename, STATUS='OLD')                             
      NM=0
        DO 3421 IM=1,NMED
        LOK(IM)=0
        IF ((IRAYLM(IM).EQ.1)) THEN
          WRITE(6,3430)IM
3430      FORMAT(' RAYLEIGH OPTION REQUESTED FOR MEDIUM NUMBER',I3,/)           
        END IF
3421  CONTINUE
3422  CONTINUE
3440  CONTINUE
3441    CONTINUE
3450    CONTINUE
3451      CONTINUE
          READ(KMPI,3260,END=3460)MBUF
            DO 3471 IB=1,LMDL
            IF((MBUF(IB).NE.MDLABL(IB)))GO TO 3451
3471      CONTINUE
3472      CONTINUE
3480      CONTINUE
            DO 3481 IM=1,NMED
              DO 3491 IB=1,LMDN
              IL=LMDL+IB
              IF((MBUF(IL).NE.MEDIA(IB,IM)))GO TO 3481
              IF((IB.EQ.LMDN))GO TO 3452
3491        CONTINUE
3492        CONTINUE
3481      CONTINUE
3482      CONTINUE
        GO TO 3451
3452    CONTINUE
        IF((LOK(IM).NE.0))GO TO 3450
        LOK(IM)=1
        NM=NM+1
        READ(KMPI,1,ERR=3500) (MBUF(I),I=1,5),RHO(IM),NNE(IM),IUNRST(IM)
     *  ,EPSTFL(IM),IAPRIM(IM)
1       FORMAT(5A1,5X,F11.0,4X,I2,9X,I1,9X,I1,9X,I1)
        GO TO 3510
3500    BACKSPACE(KMPI)
        READ(KMPI,2)(MBUF(I),I=1,5),RHO(IM),NNE(IM),IUNRST(IM),EPSTFL(IM
     *  ), IAPRIM(IM)
2       FORMAT(5A1,5X,F11.0,4X,I2,26X,I1,9X,I1,9X,I1)
3510    CONTINUE
          DO 3511 IE=1,NNE(IM)
          READ(KMPI,3520)(MBUF(I),I=1,6),(ASYM(IM,IE,I),I=1,2), ZELEM(IM
     *    ,IE),WA(IM,IE),PZ(IM,IE),RHOZ(IM,IE)
3520      FORMAT (6A1,2A1,3X,F3.0,3X,F9.0,4X,F12.0,6X,F12.0)
3511    CONTINUE
3512    CONTINUE
        READ(KMPI,3250) RLC(IM),AE(IM),AP(IM),UE(IM),UP(IM)
        TE(IM)=AE(IM)-RM
        THMOLL(IM)=TE(IM)*2. + RM
        READ(KMPI,3240) MSGE(IM),MGE(IM),MSEKE(IM),MEKE(IM),MLEKE(IM),MC
     *  MFP(IM),MRANGE(IM),IRAYL
        NSGE=MSGE(IM)
        NGE=MGE(IM)
        NSEKE=MSEKE(IM)
        NEKE=MEKE(IM)
        NLEKE=MLEKE(IM)
        NCMFP=MCMFP(IM)
        NRANGE=MRANGE(IM)
        READ(KMPI,3250)(DL1(I,IM),DL2(I,IM),DL3(I,IM),DL4(I,IM),DL5(I,IM
     *  ),DL6(I,IM),I=1,6)
        READ(KMPI,3250)DELCM(IM),(ALPHI(I,IM),BPAR(I,IM),DELPOS(I,IM),I=
     *  1,2)
        READ(KMPI,3250)XR0(IM),TEFF0(IM),BLCC(IM),XCC(IM)
        READ(KMPI,3250)EKE0(IM),EKE1(IM)
        READ(KMPI,3250) (ESIG0(I,IM),ESIG1(I,IM),PSIG0(I,IM),PSIG1(I,IM)
     *  ,EDEDX0(I,IM),EDEDX1(I,IM),PDEDX0(I,IM),PDEDX1(I,IM),EBR10(I,IM)
     *  ,EBR11(I,IM),PBR10(I,IM),PBR11(I,IM),PBR20(I,IM),PBR21(I,IM),TMX
     *  S0(I,IM),TMXS1(I,IM),I=1,NEKE)
        READ(KMPI,3250)EBINDA(IM),GE0(IM),GE1(IM)
        READ(KMPI,3250)(GMFP0(I,IM),GMFP1(I,IM),GBR10(I,IM),GBR11(I,IM),
     *  GBR20(I,IM),GBR21(I,IM),I=1,NGE)
        IF ((IRAYLM(IM).EQ.1.AND.IRAYL.NE.1)) THEN
          WRITE(6,3530)IM
3530      FORMAT(' STOPPED IN HATCH: REQUESTED RAYLEIGH OPTION FOR MEDIU        
     *M',I3, /,' BUT RAYLEIGH DATA NOT INCLUDED IN DATA CREATED BY PEGS.        
     *')                                                                        
          STOP
        END IF
        IF ((IRAYL.EQ.1)) THEN
          READ(KMPI,3240) NGR(IM)
          NGRIM=NGR(IM)
          READ(KMPI,3250)RCO0(IM),RCO1(IM)
          READ(KMPI,3250)(RSCT0(I,IM),RSCT1(I,IM),I=1,NGRIM)
          READ(KMPI,3250)(COHE0(I,IM),COHE1(I,IM),I=1,NGE)
          IF ((IRAYLM(IM).NE.1)) THEN
            WRITE(6,3540)IM
3540        FORMAT(' RAYLEIGH DATA AVAILABLE FOR MEDIUM',I3, ' BUT OPTIO        
     *N NOT REQUESTED.',/)                                                      
          END IF
        END IF
        IF((NM.GE.NMED))GO TO3442
      GO TO 3441
3442  CONTINUE
      DUNITR=DUNIT
      IF ((DUNIT.LT.0.0)) THEN
        ID=MAX0(1,MIN0(5,IFIX(-DUNIT)))
        DUNIT=RLC(ID)
      END IF
      IF ((DUNIT.NE.1.0)) THEN
        WRITE(6,3550)DUNITR,DUNIT
3550    FORMAT(' DUNIT REQUESTED&USED ARE:',1PE14.5,E14.5,'(CM.)')              
      END IF
        DO 3561 IM=1,NMED
        DFACT=RLC(IM)/DUNIT
        DFACTI=1.0/DFACT
        I=1
          GO TO 3573
3571      I=I+1
3573      IF(I-(MEKE(IM)).GT.0)GO TO 3572
          ESIG0(I,IM)=ESIG0(I,IM)*DFACTI
          ESIG1(I,IM)=ESIG1(I,IM)*DFACTI
          PSIG0(I,IM)=PSIG0(I,IM)*DFACTI
          PSIG1(I,IM)=PSIG1(I,IM)*DFACTI
          EDEDX0(I,IM)=EDEDX0(I,IM)*DFACTI
          EDEDX1(I,IM)=EDEDX1(I,IM)*DFACTI
          PDEDX0(I,IM)=PDEDX0(I,IM)*DFACTI
          PDEDX1(I,IM)=PDEDX1(I,IM)*DFACTI
          TMXS0(I,IM)=TMXS0(I,IM)*DFACT
          TMXS1(I,IM)=TMXS1(I,IM)*DFACT
        GO TO 3571
3572    CONTINUE
        TEFF0(IM)=TEFF0(IM)*DFACT
        BLCC(IM)=BLCC(IM)*DFACTI
        XCC(IM)=XCC(IM)*SQRT(DFACTI)
        RLDU(IM)=RLC(IM)/DUNIT
        I=1
          GO TO 3583
3581      I=I+1
3583      IF(I-(MGE(IM)).GT.0)GO TO 3582
          GMFP0(I,IM)=GMFP0(I,IM)*DFACT
          GMFP1(I,IM)=GMFP1(I,IM)*DFACT
        GO TO 3581
3582    CONTINUE
3561  CONTINUE
3562  CONTINUE
      VACDST=VACDST*DUNITO/DUNIT
      DUNITO=DUNIT
        DO 3591 JR=1,numgeom
        MD=MED(JR)
        IF (((MD.GE.1).AND.(MD.LE.NMED))) THEN
          ECUT(JR)=max(ECUT(JR),AE(MD))
          PCUT(JR)=max(PCUT(JR),AP(MD))
          IF ((RHOR(JR).EQ.0.0)) THEN
            RHOR(JR)=RHO(MD)
          END IF
        END IF
3591  CONTINUE
3592  CONTINUE
      IF ((IBRDST.EQ.1)) THEN
          DO 3601 IM=1,NMED
          ZBRANG(IM)=0.0
          PZNORM=0.0
            DO 3611 IE=1,NNE(IM)
            ZBRANG(IM)= ZBRANG(IM)+PZ(IM,IE)*ZELEM(IM,IE)*(ZELEM(IM,IE)+
     *      1.0)
            PZNORM=PZNORM+PZ(IM,IE)
3611      CONTINUE
3612      CONTINUE
          ZBRANG(IM)=(8.116224E-05)*(ZBRANG(IM)/PZNORM)**(1./3.)
3601    CONTINUE
3602    CONTINUE
      END IF
      IF ((IPRDST.GT.0)) THEN
          DO 3621 IM=1,NMED
          ZBRANG(IM)=0.0
          PZNORM=0.0
            DO 3631 IE=1,NNE(IM)
            ZBRANG(IM)= ZBRANG(IM)+PZ(IM,IE)*ZELEM(IM,IE)*(ZELEM(IM,IE)+
     *      1.0)
            PZNORM=PZNORM+PZ(IM,IE)
3631      CONTINUE
3632      CONTINUE
          ZBRANG(IM)=(8.116224E-05)*(ZBRANG(IM)/PZNORM)**(1./3.)
3621    CONTINUE
3622    CONTINUE
      END IF
      call mscati
      call EDGSET(1,1)
      call init_compton
      call fix_brems
      IF (( ibr_nist .EQ. 1 )) THEN
        call init_nist_brems
      END IF
      IF ((NMED.EQ.1)) THEN
        WRITE(6,3640)
3640    FORMAT(' EGSnrc SUCCESSFULLY ''HATCHED'' FOR ONE MEDIUM.')              
      ELSE
        WRITE(6,3650)NMED
3650    FORMAT(' EGSnrc SUCCESSFULLY ''HATCHED'' FOR ',I5,' MEDIA.')            
      END IF
      CLOSE(UNIT=KMPI,STATUS='KEEP')                                            
      RETURN
3460  WRITE(6,3660)KMPI
3660  FORMAT(' END OF FILE ON UNIT ',I2,//, ' PROGRAM STOPPED IN HATCH B        
     *ECAUSE THE',/, ' FOLLOWING NAMES WERE NOT RECOGNIZED:',/)                 
        DO 3671 IM=1,NMED
        IF ((LOK(IM).NE.1)) THEN
          WRITE(6,3680)(MEDIA(I,IM),I=1,LMDN)
3680      FORMAT(40X,'''',24A1,'''')                                            
        END IF
3671  CONTINUE
3672  CONTINUE
      STOP
      END
