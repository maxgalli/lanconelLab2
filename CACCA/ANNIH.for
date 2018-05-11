      SUBROUTINE ANNIH
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      DOUBLE PRECISION PAVIP,  PESG1,  PESG2
      real*4 AVIP, A, G, T, P, POT, EP0, WSAMP, RNNO01, RNNO02, EP, REJF
     *, ESG1, ESG2, aa, bb, cc, sinpsi, sindel, cosdel, us, vs,cphi,sphi
      integer*4 ibr
      real*4 xphi, xphi2, yphi, yphi2, rhophi2
      NPold = NP
      PAVIP=E(NP)+PRM
      AVIP=PAVIP
      A=AVIP/RM
      G=A-1.0
      T=G-1.0
      P=SQRT(A*T)
      POT=P/T
      EP0=1.0/(A+P)
      WSAMP=LOG((1.0-EP0)/EP0)
      aa = u(np)
      bb = v(np)
      cc = w(np)
      sinpsi = aa*aa + bb*bb
      IF (( sinpsi .GT. 1e-20 )) THEN
        sinpsi = sqrt(sinpsi)
        sindel = bb/sinpsi
        cosdel = aa/sinpsi
      END IF
      IF (( nbr_split .GT. 1 )) THEN
        wt(np) = wt(np)/nbr_split
      END IF
        DO 2811 ibr=1,nbr_split
        IF (( np+1 .GT. 50 )) THEN
          WRITE(6,2820)np+1
2820      FORMAT(//' Stack overflow in ANNIH! np = ',i6, ' Increase $MXS        
     *TACK and try again'//)                                                    
          stop
        END IF
2831    CONTINUE
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          RNNO01 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          EP=EP0*EXP(RNNO01*WSAMP)
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          RNNO02 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          REJF = 1 - (EP*A-1)**2/(EP*(A*A-2))
          IF(((RNNO02 .LE. REJF)))GO TO2832
        GO TO 2831
2832    CONTINUE
        ESG1=AVIP*EP
        PESG1=ESG1
        E(NP)=PESG1
        IQ(NP)=0
        X(np)=X(NPold)
        Y(np)=Y(NPold)
        Z(np)=Z(NPold)
        IR(np)=IR(NPold)
        WT(np)=WT(NPold)
        DNEAR(np)=DNEAR(NPold)
        LATCH(np)=LATCH(NPold)
        COSTHE=MIN(1.0,(ESG1-RM)*POT/ESG1)
        SINTHE=SQRT(1.0-COSTHE*COSTHE)
2841    CONTINUE
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          xphi = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          xphi = 2*xphi - 1
          xphi2 = xphi*xphi
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          yphi = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          yphi2 = yphi*yphi
          rhophi2 = xphi2 + yphi2
          IF(rhophi2.LE.1)GO TO2842
        GO TO 2841
2842    CONTINUE
        rhophi2 = 1/rhophi2
        cphi = (xphi2 - yphi2)*rhophi2
        sphi = 2*xphi*yphi*rhophi2
        IF (( sinpsi .GE. 1e-10 )) THEN
          us = sinthe*cphi
          vs = sinthe*sphi
          u(np) = cc*cosdel*us - sindel*vs + aa*costhe
          v(np) = cc*sindel*us + cosdel*vs + bb*costhe
          w(np) = cc*costhe - sinpsi*us
        ELSE
          u(np) = sinthe*cphi
          v(np) = sinthe*sphi
          w(np) = costhe
        END IF
        np = np + 1
        PESG2=PAVIP-PESG1
        esg2 = pesg2
        e(np) = pesg2
        iq(np) = 0
        X(np)=X(NPold)
        Y(np)=Y(NPold)
        Z(np)=Z(NPold)
        IR(np)=IR(NPold)
        WT(np)=WT(NPold)
        DNEAR(np)=DNEAR(NPold)
        LATCH(np)=LATCH(NPold)
        COSTHE=MIN(1.0,(ESG2-RM)*POT/ESG2)
        SINTHE=-SQRT(1.0-COSTHE*COSTHE)
        IF (( sinpsi .GE. 1e-10 )) THEN
          us = sinthe*cphi
          vs = sinthe*sphi
          u(np) = cc*cosdel*us - sindel*vs + aa*costhe
          v(np) = cc*sindel*us + cosdel*vs + bb*costhe
          w(np) = cc*costhe - sinpsi*us
        ELSE
          u(np) = sinthe*cphi
          v(np) = sinthe*sphi
          w(np) = costhe
        END IF
        np = np + 1
2811  CONTINUE
2812  CONTINUE
      np = np-1
      RETURN
      END
