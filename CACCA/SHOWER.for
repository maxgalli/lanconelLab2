      SUBROUTINE SHOWER(IQI,EI,XI,YI,ZI,UI,VI,WI,IRI,WTI)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      real*4 EI,  XI,YI,ZI, UI,VI,WI, WTI
      integer*4 IQI,  IRI
      DOUBLE PRECISION DEG,  DPGL,  DEI,  DPI,  DCSTH,  DCOSTH,  PI0MSQ
      real*4 DNEARI,  CSTH
      integer*4 IRCODE
      DATA PI0MSQ/1.8215416D4/
      NP=1
      NPold = NP
      DNEARI=0.0
      IQ(1)=IQI
      E(1)=EI
      U(1)=UI
      V(1)=VI
      W(1)=WI
      X(1)=XI
      Y(1)=YI
      Z(1)=ZI
      IR(1)=IRI
      WT(1)=WTI
      DNEAR(1)=DNEARI
      LATCH(1)=LATCHI
      IF ((IQI.EQ.2)) THEN
        IF ((EI**2.LE.PI0MSQ)) THEN
          WRITE(6,4930)EI
4930      FORMAT(//,' STOPPED IN SUBROUTINE SHOWER---PI-ZERO OPTION INVO        
     *KED', /,' BUT THE TOTAL ENERGY WAS TOO SMALL (EI=',G15.5,' MEV)')         
          STOP
        END IF
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        CSTH = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        DCSTH=CSTH
        DEI=EI
        DPI=DSQRT(DEI*DEI-PI0MSQ)
        DEG=DEI+DPI*DCSTH
        DPGL=DPI+DEI*DCSTH
        DCOSTH=DPGL/DEG
        COSTHE=DCOSTH
        SINTHE=DSQRT(1.D0-DCOSTH*DCOSTH)
        IQ(1)=0
        E(1)=DEG/2.
        CALL UPHI(2,1)
        NP=2
        DEG=DEI-DPI*DCSTH
        DPGL=DPI-DEI*DCSTH
        DCOSTH=DPGL/DEG
        COSTHE=DCOSTH
        SINTHE=-DSQRT(1.D0-DCOSTH*DCOSTH)
        IQ(2)=0
        E(2)=DEG/2.
        CALL UPHI(3,2)
      END IF
4940  CONTINUE
4941    CONTINUE
        IF((IQ(NP).EQ.0))GO TO 4950
4961    CONTINUE
4970      CALL ELECTR(IRCODE)
          IF((IRCODE.EQ.2))GO TO4962
4950      CALL PHOTON(IRCODE)
          IF((IRCODE.EQ.2))GO TO4962
        GO TO 4961
4962    CONTINUE
        IF((NP.LE.0))GO TO4942
      GO TO 4941
4942  CONTINUE
      RETURN
      END
