      SUBROUTINE BREMS
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      DOUBLE PRECISION PEIE,  PESG,  PESE
      real*4 EIE,  EKIN,  brmin,  waux,  aux,  r1,  ajj,  alias_sample1,
     * RNNO06,  RNNO07,  BR,  ESG,  ESE,  DELTA,  phi1,  phi2,  REJF
      real*4 a, b, c, sinpsi, sindel, cosdel, us, vs, ztarg, tteie, beta
     *, y2max, y2maxi, ttese, rjarg1, rjarg2, rjarg3, rejmin, rejmid
     *, rejmax, rejtop, rejtst, esedei, y2tst, y2tst1, rtest, xphi, yphi
     *, xphi2, yphi2, rhophi2, cphi, sphi
      integer*4 L, L1, ibr, jj, j
      IF((nbr_split .LT. 1))return
      NPold = NP
      PEIE=E(NP)
      EIE=PEIE
      IF ((EIE.LT.50.0)) THEN
        L=1
      ELSE
        L=3
      END IF
      L1 = L+1
      ekin = peie-prm
      brmin = ap(medium)/ekin
      waux = elke - log_ap(medium)
      IF (( ibrdst .GE. 0 )) THEN
        a = u(np)
        b = v(np)
        c = w(np)
        sinpsi = a*a + b*b
        IF (( sinpsi .GT. 1e-20 )) THEN
          sinpsi = sqrt(sinpsi)
          sindel = b/sinpsi
          cosdel = a/sinpsi
        END IF
        ztarg = zbrang(medium)
        tteie = eie/rm
        beta = sqrt((tteie-1)*(tteie+1))/tteie
        y2max = 2*beta*(1+beta)*tteie*tteie
        y2maxi = 1/y2max
      END IF
      IF (( ibr_nist .EQ. 1 )) THEN
        ajj = 1 + (waux + log_ap(medium) - nb_lemin(medium))*nb_dlei(med
     *  ium)
        jj = ajj
        ajj = ajj - jj
        IF (( jj .GT. 100 )) THEN
          jj = 100
          ajj = -1
        END IF
      END IF
        DO 2861 ibr=1,nbr_split
        IF (( ibr_nist .EQ. 1 )) THEN
          IF (( ekin .GT. nb_emin(medium) )) THEN
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            r1 = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            IF (( r1 .LT. ajj )) THEN
              j = jj+1
            ELSE
              j = jj
            END IF
            br = alias_sample1(50,nb_xdata(0,j,medium), nb_fdata(0,j,med
     *      ium), nb_wdata(1,j,medium),nb_idata(1,j,medium))
          ELSE
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            br = rng_array(rng_seed)
            rng_seed = rng_seed + 1
          END IF
          esg = ap(medium)*exp(br*waux)
          pesg = esg
          pese = peie - pesg
          ese = pese
        ELSE
2871      CONTINUE
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            rnno06 = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            rnno07 = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            br = brmin*exp(rnno06*waux)
            esg = ekin*br
            pesg = esg
            pese = peie - pesg
            ese = pese
            delta = esg/eie/ese*delcm(medium)
            aux = ese/eie
            IF (( delta .LT. 1 )) THEN
              phi1 = dl1(l,medium)+delta*(dl2(l,medium)+delta*dl3(l,medi
     *        um))
              phi2 = dl1(l1,medium)+delta*(dl2(l1,medium)+ delta*dl3(l1,
     *        medium))
            ELSE
              phi1 = dl4(l,medium)+dl5(l,medium)*log(delta+dl6(l,medium)
     *        )
              phi2 = phi1
            END IF
            rejf = (1+aux*aux)*phi1 - 2*aux*phi2/3
            IF(((rnno07 .LT. rejf)))GO TO2872
          GO TO 2871
2872      CONTINUE
        END IF
        np=np+1
        IF (( np .GT. 50 )) THEN
          WRITE(6,2880)np
2880      FORMAT(//' Stack overflow in BREMS! np = ',i6, ' Increase $MXS        
     *TACK and try again'//)                                                    
          stop
        END IF
        e(np) = pesg
        iq(np) = 0
        X(np)=X(NPold)
        Y(np)=Y(NPold)
        Z(np)=Z(NPold)
        IR(np)=IR(NPold)
        WT(np)=WT(NPold)
        DNEAR(np)=DNEAR(NPold)
        LATCH(np)=LATCH(NPold)
        wt(np) = wt(np)/nbr_split
        IF (( ibrdst .LT. 0 )) THEN
          u(np) = u(npold)
          v(np) = v(npold)
          w(np) = w(npold)
        ELSE
          IF (( ibrdst .EQ. 1 )) THEN
            ttese = ese/rm
            esedei = ttese/tteie
            rjarg1 = 1+esedei*esedei
            rjarg2 = 3*rjarg1 - 2*esedei
            rjarg3 = ((1-esedei)/(2*tteie*esedei))**2
            Y2TST1=(1.+0.0)**2
            REJMIN= (4.+LOG(RJARG3+ZTARG/Y2TST1))*(4.*ESEDEI*0.0/Y2TST1-
     *      RJARG1)+RJARG2
            Y2TST1=(1.+1.0)**2
            REJMID= (4.+LOG(RJARG3+ZTARG/Y2TST1))*(4.*ESEDEI*1.0/Y2TST1-
     *      RJARG1)+RJARG2
            Y2TST1=(1.+y2max)**2
            REJMAX= (4.+LOG(RJARG3+ZTARG/Y2TST1))*(4.*ESEDEI*y2max/Y2TST
     *      1-RJARG1)+RJARG2
            rejtop = max(rejmin,rejmid,rejmax)
2891        CONTINUE
              IF (( rng_seed .GT. 24 )) THEN
                call ranlux(rng_array)
                rng_seed = 1
              END IF
              y2tst = rng_array(rng_seed)
              rng_seed = rng_seed + 1
              y2tst = y2tst/(1-y2tst+y2maxi)
              Y2TST1=(1.+Y2TST)**2
              REJTST= (4.+LOG(RJARG3+ZTARG/Y2TST1))*(4.*ESEDEI*Y2TST/Y2T
     *        ST1-RJARG1)+RJARG2
              IF (( rng_seed .GT. 24 )) THEN
                call ranlux(rng_array)
                rng_seed = 1
              END IF
              rtest = rng_array(rng_seed)
              rng_seed = rng_seed + 1
              IF(((rtest*rejtop .LE. REJTST)))GO TO2892
            GO TO 2891
2892        CONTINUE
          ELSE
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            y2tst = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            y2tst = y2tst/(1-y2tst+y2maxi)
          END IF
          costhe = 1 - 2*y2tst*y2maxi
          sinthe = sqrt(max((1-costhe)*(1+costhe),0.0))
2901      CONTINUE
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
            IF(rhophi2.LE.1)GO TO2902
          GO TO 2901
2902      CONTINUE
          rhophi2 = 1/rhophi2
          cphi = (xphi2 - yphi2)*rhophi2
          sphi = 2*xphi*yphi*rhophi2
          IF (( sinpsi .GE. 1e-10 )) THEN
            us = sinthe*cphi
            vs = sinthe*sphi
            u(np) = c*cosdel*us - sindel*vs + a*costhe
            v(np) = c*sindel*us + cosdel*vs + b*costhe
            w(np) = c*costhe - sinpsi*us
          ELSE
            u(np) = sinthe*cphi
            v(np) = sinthe*sphi
            w(np) = costhe
          END IF
        END IF
2861  CONTINUE
2862  CONTINUE
      e(npold) = pese
      RETURN
      END
