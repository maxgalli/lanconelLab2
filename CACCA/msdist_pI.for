      subroutine msdist_pI(e0,eloss,tustep_orig,rhof_orig,medium_orig,
     *qel,spin_effects_orig,u0,v0,w0,x0,y0,z0,us,vs,ws,xf,yf,zf,
     *ustep_orig)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none      
c
      include 'egs4fcomm.for'
c

      real*4 e0, eloss, rhof_orig, tustep_orig, u0, v0, w0, x0, y0, z0
      integer*4 medium_orig, qel
      logical spin_effects_orig
      real*4 us,  vs,  ws,  xf,  yf,  zf,  ustep_orig
      real*4 blccc,xcccc,z_orig,r,z2,r2,r2max,chia2,chilog,cphi0,
     * cphi, sphi, e_orig, elke_orig, beta2, etap, xi_corr, ms_corr, 
     *epsilon,  temp,  factor,lambda,p2,p2i,q1,  rhophi2,  sint,  
     *sint0,  sint02,  sint0i,sphi0,u2p,ut,vt,wt_orig,xi,xphi,  
     *xphi2,  yphi,  yphi2
      logical find_index,  spin_index
      integer*4 lelke
      blccc = blcc(medium)
      xcccc = xcc(medium)
      e_orig = e0 - 0.5*eloss
      p2 = e_orig*(e_orig + rmt2)
      p2i = 1/p2
      chia2 = xcccc*p2i/(4*blccc)
      beta2 = p2/(p2 + rmsq)
      lambda = tustep_orig*rhof_orig*blccc/beta2
      factor = 1/(1 + 0.9784671*e_orig)
      epsilon= eloss/e0
      epsilon= epsilon/(1-0.5*epsilon)
      temp = 0.25*(1 - factor*(1 - 0.333333*factor))*epsilon**2
      lambda = lambda*(1 + temp)
      IF (( spin_effects_orig )) THEN
        elke_orig = Log(e_orig)
        Lelke=eke1(MEDIUM)*elke_orig+eke0(MEDIUM)
        IF (( lelke .LT. 1 )) THEN
          lelke = 1
          elke_orig = (1 - eke0(medium))/eke1(medium)
        END IF
        IF (( qel .EQ. 0 )) THEN
          etap=etae_ms1(Lelke,MEDIUM)*elke_orig+etae_ms0(Lelke,MEDIUM)
          xi_corr=q1ce_ms1(Lelke,MEDIUM)*elke_orig+
     *q1ce_ms0(Lelke,MEDIUM)
        ELSE
          etap=etap_ms1(Lelke,MEDIUM)*elke_orig+etap_ms0(Lelke,MEDIUM)
          xi_corr=q1cp_ms1(Lelke,MEDIUM)*elke_orig+
     *q1cp_ms0(Lelke,MEDIUM)
        END IF
        ms_corr=blcce1(Lelke,MEDIUM)*elke_orig+blcce0(Lelke,MEDIUM)
      ELSE
        etap = 1
        xi_corr = 1
        ms_corr = 1
      END IF
      chia2 = xcccc*p2i/(4*blccc)*etap
      lambda = lambda/etap/(1+chia2)*ms_corr
      chilog = Log(1 + 1/chia2)
      q1 = 2*chia2*(chilog*(1 + chia2) - 1)
      xi = q1*lambda
      find_index = .true.
      spin_index = .true.
      call mscat(lambda,chia2,xi,elke_orig,beta2,qel,medium,
     * spin_effects_orig,find_index,spin_index, ws,sint)
4481  CONTINUE
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
        IF(rhophi2.LE.1)GO TO4482
      GO TO 4481
4482  CONTINUE
      rhophi2 = 1/rhophi2
      cphi = (xphi2 - yphi2)*rhophi2
      sphi = 2*xphi*yphi*rhophi2
      us = sint*cphi
      vs = sint*sphi
      xi = xi*xi_corr
      IF (( xi .LT. 0.1 )) THEN
        z_orig = 1 - xi*(0.5 - xi*(0.166666667 - 0.041666667*xi))
      ELSE
        z_orig = (1 - Exp(-xi))/xi
      END IF
      r = 0.5*sint
      r2 = r*r
      z2 = z_orig*z_orig
      r2max = 1 - z2
      IF (( r2max .LT. r2 )) THEN
        r2 = r2max
        r = Sqrt(r2)
      END IF
      ut = r*cphi
      vt = r*sphi
      wt_orig = z_orig
      ustep_orig = Sqrt(z2 + r2)*tustep_orig
      sint02 = u0**2 + v0**2
      IF ((sint02 .GT. 1e-20)) THEN
        sint0 = sqrt(sint02)
        sint0i = 1/sint0
        cphi0 = sint0i*u0
        sphi0 = sint0i*v0
        u2p = w0*us + sint0*ws
        ws = w0*ws - sint0*us
        us = u2p*cphi0 - vs*sphi0
        vs = u2p*sphi0 + vs*cphi0
        u2p = w0*ut + sint0*wt_orig
        wt_orig = w0*wt_orig - sint0*ut
        ut = u2p*cphi0 - vt*sphi0
        vt = u2p*sphi0 + vt*cphi0
      END IF
      xf = x0 + tustep_orig*ut
      yf = y0 + tustep_orig*vt
      zf = z0 + tustep_orig*wt_orig
      return
      end
