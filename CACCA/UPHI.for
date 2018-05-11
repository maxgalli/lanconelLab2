      SUBROUTINE UPHI(IENTRY,LVL)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      integer*4 IENTRY,LVL
      real*4 CTHET,  RNNO38,  PHI,  CPHI,  A,B,C,  SINPS2,  SINPSI,  US,
     *VS,  SINDEL,COSDEL
      integer*4 IARG,  LPHI,LTHETA,LCTHET,LCPHI
      real*4 xphi,xphi2,yphi,yphi2,rhophi2
      double precision norm
      IARG=21
      IF ((IAUSFL(IARG+1).NE.0)) THEN
        CALL AUSGAB(IARG)
      END IF
      GO TO (4980,4990,5000),IENTRY
      GO TO 5010
4980  CONTINUE
      SINTHE=sin(THETA)
      CTHET=PI5D2-THETA
      COSTHE=sin(CTHET)
4990  CONTINUE
5021  CONTINUE
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
        IF(rhophi2.LE.1)GO TO5022
      GO TO 5021
5022  CONTINUE
      rhophi2 = 1/rhophi2
      cosphi = (xphi2 - yphi2)*rhophi2
      sinphi = 2*xphi*yphi*rhophi2
5000  GO TO (5030,5040,5050),LVL
      GO TO 5010
5030  A=U(NP)
      B=V(NP)
      C=W(NP)
      GO TO 5060
5050  A=U(NP-1)
      B=V(NP-1)
      C=W(NP-1)
5040  X(NP)=X(NP-1)
      Y(NP)=Y(NP-1)
      Z(NP)=Z(NP-1)
      IR(NP)=IR(NP-1)
      WT(NP)=WT(NP-1)
      DNEAR(NP)=DNEAR(NP-1)
      LATCH(NP)=LATCH(NP-1)
5060  SINPS2=A*A+B*B
      IF ((SINPS2.LT.1.0E-20)) THEN
        U(NP)=SINTHE*COSPHI
        V(NP)=SINTHE*SINPHI
c
c known bug corrected - mirko
c
c        W(NP)=COSTHE
c
        W(NP)=C*COSTHE
cnico normalize u,v,w
	norm = u(np)*u(np)+v(np)*v(np)+w(np)*w(np)
	norm = 1./sqrt(norm)
	u(np) = u(np) * norm
	v(np) = v(np) * norm
	w(np) = w(np) * norm
cnico
c
c-------------------------
c
      ELSE
        SINPSI=SQRT(SINPS2)
        US=SINTHE*COSPHI
        VS=SINTHE*SINPHI
        SINDEL=B/SINPSI
        COSDEL=A/SINPSI
        U(NP)=C*COSDEL*US-SINDEL*VS+A*COSTHE
        V(NP)=C*SINDEL*US+COSDEL*VS+B*COSTHE
        W(NP)=-SINPSI*US+C*COSTHE
cnico normalize u,v,w
	norm = u(np)*u(np)+v(np)*v(np)+w(np)*w(np)
	norm = 1./sqrt(norm)
	u(np) = u(np) * norm
	v(np) = v(np) * norm
	w(np) = w(np) * norm
cnico
      END IF
      IARG=22
      IF ((IAUSFL(IARG+1).NE.0)) THEN
        CALL AUSGAB(IARG)
      END IF
      RETURN
5010  WRITE(6,5070)IENTRY,LVL
5070  FORMAT(' STOPPED IN UPHI WITH IENTRY,LVL=',2I6)                           
      STOP
      END
