      SUBROUTINE PHOTON(IRCODE)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      integer*4 IRCODE
      DOUBLE PRECISION PEIG
      real*4 EIG, RNNO35, DPMFP, GMFPR0, GMFP, COHFAC, RNNO37, XXX, X2
     *, Q2, CSQTHE, REJF, RNNORJ, RNNO36, GBR1, GBR2, T
      integer*4 IARG, IDR, IRL, LGLE, LXXX
      IRCODE=1
      PEIG=E(NP)
      EIG=PEIG
      IRL=IR(NP)
      MEDIUM=MED(IRL)
      IF ((EIG.LE.PCUT(IRL))) THEN
        GO TO 4860
      END IF
4870  CONTINUE
4871    CONTINUE
        GLE=LOG(EIG)
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        RNNO35 = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        IF ((RNNO35.EQ.0.0)) THEN
          RNNO35=1.E-30
        END IF
        DPMFP=-LOG(RNNO35)
        IROLD=IR(NP)
4880    CONTINUE
4881      CONTINUE
          IF ((MEDIUM.NE.0)) THEN
            LGLE=GE1(MEDIUM)*GLE+GE0(MEDIUM)
            GMFPR0=GMFP1(LGLE,MEDIUM)*GLE+GMFP0(LGLE,MEDIUM)
          END IF
4890      CONTINUE
4891        CONTINUE
            IF ((MEDIUM.EQ.0)) THEN
              TSTEP=VACDST
            ELSE
              RHOF=RHOR(IRL)/RHO(MEDIUM)
              GMFP=GMFPR0/RHOF
              IF ((IRAYLR(IRL).EQ.1)) THEN
                COHFAC=COHE1(LGLE,MEDIUM)*GLE+COHE0(LGLE,MEDIUM)
                GMFP=GMFP*COHFAC
              END IF
              TSTEP=GMFP*DPMFP
            END IF
            IRNEW=IR(NP)
            IDISC=0
            USTEP=TSTEP
            TUSTEP=USTEP
            IF (( ustep .GT. dnear(np) .OR. wt(np) .LE. 0 )) THEN
              call howfar
            END IF
            IF ((IDISC.GT.0)) THEN
              GO TO 4900
            END IF
            VSTEP=USTEP
            TVSTEP=VSTEP
            EDEP=PZERO
            IARG=0
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
            X(NP)=X(NP)+U(NP)*USTEP
            Y(NP)=Y(NP)+V(NP)*USTEP
            Z(NP)=Z(NP)+W(NP)*USTEP
            DNEAR(NP)=DNEAR(NP)-USTEP
            IF ((MEDIUM.NE.0)) THEN
              DPMFP=MAX(0.,DPMFP-USTEP/GMFP)
            END IF
            IROLD=IR(NP)
            MEDOLD=MEDIUM
            IF ((IRNEW.NE.IROLD)) THEN
              IR(NP)=IRNEW
              IRL=IRNEW
              MEDIUM=MED(IRL)
            END IF
            IARG=5
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
            IF ((EIG.LE.PCUT(IRL))) THEN
              GO TO 4860
            END IF
            IF((IDISC.LT.0))GO TO 4900
            IF((MEDIUM.NE.MEDOLD))GO TO 4892
            IF ((MEDIUM.NE.0.AND.DPMFP.LE.1.E-5)) THEN
              GO TO 4882
            END IF
          GO TO 4891
4892      CONTINUE
        GO TO 4881
4882    CONTINUE
        IF ((IRAYLR(IRL).EQ.1)) THEN
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          RNNO37 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF ((RNNO37.LE.(1.0-COHFAC))) THEN
            IARG=23
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
            NPold = NP
4910        CONTINUE
4911          CONTINUE
              IF (( rng_seed .GT. 24 )) THEN
                call ranlux(rng_array)
                rng_seed = 1
              END IF
              XXX = rng_array(rng_seed)
              rng_seed = rng_seed + 1
              LXXX=RCO1(MEDIUM)*XXX+RCO0(MEDIUM)
              X2=RSCT1(LXXX,MEDIUM)*XXX+RSCT0(LXXX,MEDIUM)
              Q2=X2*RMSQ/(20.60744*20.60744)
              COSTHE=1.-Q2/(2.*E(NP)*E(NP))
              IF((ABS(COSTHE).GT.1.0))GO TO 4910
              CSQTHE=COSTHE*COSTHE
              REJF=(1.0+CSQTHE)/2.0
              IF (( rng_seed .GT. 24 )) THEN
                call ranlux(rng_array)
                rng_seed = 1
              END IF
              RNNORJ = rng_array(rng_seed)
              rng_seed = rng_seed + 1
              IF(((RNNORJ .LE. REJF)))GO TO4912
            GO TO 4911
4912        CONTINUE
            SINTHE=SQRT(1.0-CSQTHE)
            CALL UPHI(2,1)
            IARG=24
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
            GOTO 4870
          END IF
        END IF
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        RNNO36 = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        GBR1=GBR11(LGLE,MEDIUM)*GLE+GBR10(LGLE,MEDIUM)
        IF (((RNNO36.LE.GBR1).AND.(E(NP).GT.RMT2) )) THEN
          IARG=15
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          CALL PAIR
          IARG=16
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          IF (( iq(np) .NE. 0 )) THEN
            GO TO 4872
          ELSE
            goto 4920
          END IF
        END IF
        GBR2=GBR21(LGLE,MEDIUM)*GLE+GBR20(LGLE,MEDIUM)
        IF ((RNNO36.LT.GBR2)) THEN
          IARG=17
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          CALL COMPT
          IARG=18
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          IF((IQ(NP).NE.0))GO TO 4872
        ELSE
          IARG=19
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          CALL PHOTO
          IF ((NP.EQ.0)) THEN
            IRCODE=2
            RETURN
          END IF
          IARG=20
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          IF((IQ(NP) .NE. 0))GO TO 4872
        END IF
4920    PEIG=E(NP)
        EIG=PEIG
        IF((EIG.LT.PCUT(IRL)))GO TO 4860
      GO TO 4871
4872  CONTINUE
      RETURN
4860  IF ((EIG.GT.AP(MEDIUM))) THEN
        IDR=1
      ELSE
        IDR=2
      END IF
      EDEP=PEIG
      IARG=IDR
      IF ((IAUSFL(IARG+1).NE.0)) THEN
        CALL AUSGAB(IARG)
      END IF
      IRCODE=2
      NP=NP-1
      RETURN
4900  EDEP=PEIG
      IARG=3
      IF ((IAUSFL(IARG+1).NE.0)) THEN
        CALL AUSGAB(IARG)
      END IF
      IRCODE=2
      NP=NP-1
      RETURN
      END
