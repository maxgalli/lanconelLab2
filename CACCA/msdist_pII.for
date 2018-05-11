      subroutine msdist_pII(e0,eloss,tustep_orig,rhof_orig,medium_orig,
     *qel,spin_effects_orig,u0,v0,w0,x0,y0,z0,  us,vs,ws,xf,yf,zf,
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
      real*4 b,  blccc,  xcccc,  c,  eta,eta1,  chia2,  chilog,  cphi0,
     * cphi1, cphi2, w1, w2, w1v2,  delta,e_orig,elke_orig,beta2,etap
     *,  xi_corr,  ms_corr, tau,  epsilon,  epsilonp,  temp,temp1, temp2
     *,  factor,  gamma,  lambda,   p2,  p2i,  q1,  rhophi2,  sint0,  si
     *nt02,  sint0i,  sint1,  sint2,  sphi0,   sphi1,  sphi2,  u2p,  u2,
     *  v2,  ut,  vt,  wt_orig,  xi,  xphi,  xphi2,  yphi,  yphi2
      logical find_index,  spin_index
      integer*4 lelke
      count_pII_steps = count_pII_steps + 1
      blccc = blcc(medium_orig)
      xcccc = xcc(medium_orig)
      e_orig = e0 - 0.5*eloss
      tau = e_orig/0.5110034
      epsilon = eloss/e0
      epsilonp= eloss/e_orig
      e_orig = e_orig * (1 - epsilonp*epsilonp*((6+tau*(10+5*tau))/
     *(tau+1)/(tau+2))/24)
      p2 = e_orig*(e_orig + rmt2)
      p2i = 1/p2
      beta2 = p2/(p2 + rmsq)
      chia2 = xcccc*p2i/(4*blccc)
      lambda = 0.5*tustep_orig*rhof_orig*blccc/beta2
      temp2 = 0.166666*(4+tau*(6+tau*(7+tau*(4+tau))))* (epsilonp/(tau+1
     *)/(tau+2))**2
      lambda = lambda*(1 - temp2)
      IF (( spin_effects_orig )) THEN
        elke_orig = Log(e_orig)
        Lelke=eke1(MEDIUM_orig)*elke_orig+eke0(MEDIUM_orig)
        IF (( lelke .LT. 1 )) THEN
          lelke = 1
          elke_orig = (1 - eke0(medium_orig))/eke1(medium_orig)
        END IF
        IF (( qel .EQ. 0 )) THEN
         etap=etae_ms1(Lelke,MEDIUM_orig)*elke_orig+etae_ms0(Lelke,
     aMEDIUM_orig)
         xi_corr=q1ce_ms1(Lelke,MEDIUM_orig)*elke_orig+q1ce_ms0(Lelke,
     b MEDIUM_orig)
         gamma=q2ce_ms1(Lelke,MEDIUM_orig)*elke_orig+q2ce_ms0(Lelke,
     c MEDIUM_orig)
        ELSE
         etap=etap_ms1(Lelke,MEDIUM_orig)*elke_orig+etap_ms0(Lelke,
     d MEDIUM_orig)
         xi_corr=q1cp_ms1(Lelke,MEDIUM_orig)*elke_orig+q1cp_ms0(Lelke,
     e MEDIUM_orig)
         gamma=q2cp_ms1(Lelke,MEDIUM_orig)*elke_orig+q2cp_ms0(Lelke,
     f MEDIUM_orig)
        END IF
        ms_corr=blcce1(Lelke,MEDIUM_orig)*elke_orig+blcce0(Lelke,
     g MEDIUM_orig)
      ELSE
        etap = 1
        xi_corr = 1
        gamma = 1
        ms_corr = 1
      END IF
      chia2 = chia2*etap
      lambda = lambda/etap/(1+chia2)*ms_corr
      chilog = Log(1 + 1/chia2)
      q1 = 2*chia2*(chilog*(1 + chia2) - 1)
      gamma = 6*chia2*(1 + chia2)*(chilog*(1 + 2*chia2) - 2)/q1*gamma
      xi = q1*lambda
      find_index = .true.
      spin_index = .true.
      call mscat(lambda,chia2,xi,elke_orig,beta2,qel,medium_orig, 
     *spin_effects_orig,find_index,spin_index, w1,sint1)
4461  CONTINUE
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
        IF(rhophi2.LE.1)GO TO4462
      GO TO 4461
4462  CONTINUE
      rhophi2 = 1/rhophi2
      cphi1 = (xphi2 - yphi2)*rhophi2
      sphi1 = 2*xphi*yphi*rhophi2
      call mscat(lambda,chia2,xi,elke_orig,beta2,qel,medium_orig, 
     *spin_effects_orig,find_index,spin_index, w2,sint2)
4471  CONTINUE
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
        IF(rhophi2.LE.1)GO TO4472
      GO TO 4471
4472  CONTINUE
      rhophi2 = 1/rhophi2
      cphi2 = (xphi2 - yphi2)*rhophi2
      sphi2 = 2*xphi*yphi*rhophi2
      u2 = sint2*cphi2
      v2 = sint2*sphi2
      u2p = w1*u2 + sint1*w2
      us = u2p*cphi1 - v2*sphi1
      vs = u2p*sphi1 + v2*cphi1
      ws = w1*w2 - sint1*u2
      xi = 2*xi*xi_corr
      IF (( rng_seed .GT. 24 )) THEN
        call ranlux(rng_array)
        rng_seed = 1
      END IF
      eta = rng_array(rng_seed)
      rng_seed = rng_seed + 1
      eta = Sqrt(eta)
      eta1 = 0.5*(1 - eta)
      delta = 0.9082483-(0.1020621-0.0263747*gamma)*xi
      temp1 = 2 + tau
      temp = (2+tau*temp1)/(tau+1)/temp1
      temp = temp - (tau+1)/(tau+2)/(chilog*(1+chia2)-1)
      temp = temp * epsilonp
      temp1 = 1 - temp
      delta = delta + 0.40824829*(epsilon*(tau+1)/(tau+2)/ (chilog*(1+ch
     *ia2)-1)/(chilog*(1+2*chia2)-2) - 0.25*temp*temp)
      b = eta*delta
      c = eta*(1-delta)
      w1v2 = w1*v2
      ut = b*sint1*cphi1 + c*(cphi1*u2 - sphi1*w1v2) + eta1*us*temp1
      vt = b*sint1*sphi1 + c*(sphi1*u2 + cphi1*w1v2) + eta1*vs*temp1
      wt_orig = eta1*(1+temp) + b*w1 + c*w2 + eta1*ws*temp1
      ustep_orig = tustep_orig*sqrt(ut*ut + vt*vt + wt_orig*wt_orig)
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
