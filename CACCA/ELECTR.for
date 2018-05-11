      SUBROUTINE ELECTR(IRCODE)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      integer*4 IRCODE
      DOUBLE PRECISION  demfp,  peie,  total_tstep,  total_de
      real*4 ekems,  elkems,  chia2,  etap,  lambda,  blccl,  xi,  xi_co
     *rr,  ms_corr, p2,  beta2,  de,  save_de,  dedx,  dedx0,  dedxmid,
     * ekei,  elkei,  aux,  ebr1,  eie,  ekef,  elkef,  ekeold,  eketmp,
     *  elktmp,  fedep,  fedep1,  tuss,  flip,  pbr1,  pbr2,  range,  rf
     *ict,  rnne1,  rnno24,  rnno25,  rnnotu,  rnnoss,  sig,  sig0,  sig
     *f,  skindepth,  ssmfp,  tmxs,  tperp,  ustep0,  uscat,  vscat,  ws
     *cat,  xtrans,  ytrans,  ztrans,  cphi,sphi
      real*4 xphi,xphi2,yphi,yphi2,rhophi2
      integer*4 iarg,  idr,  ierust,  irl,  lelec,  qel,  lelke,  lelkem
     *s,  lelkef,  lelktmp,  ibr
      logical  callhowfar,   domultiple,  dosingle,   callmsdist,
     *                findindex,
     *              spin_index,                                   comput
     *e_tstep
     *
      real*4 tau,xxx,bbb,aux1,lambda_max, sigratio
      data ierust/0/
      save ierust
      ircode = 1
      irold = ir(np)
      irl = irold
      medium = med(irl)
cmirko
c              write (*,*) 'electr1  ', x(np), y(np), z(np), vstep, u(np)
c     a, v(np), w(np),medium,sig,compute_tstep,ekef,E_array(1,medium)
cmirko
3090  CONTINUE
cmirko
c              write (*,*) 'electr2  ', x(np), y(np), z(np), vstep, u(np)
c     a, v(np), w(np),medium,sig,compute_tstep,ekef,E_array(1,medium)
cmirko
3091    CONTINUE
cmirko
c              write (*,*) 'electr3  ', x(np), y(np), z(np), vstep, u(np)
c     a, v(np), w(np),medium,sig,compute_tstep,ekef,E_array(1,medium)
cmirko
        lelec = iq(np)
        qel = (1+lelec)/2
        peie = e(np)
        eie = peie
c
        IF ((eie .LE. ecut(irl))) THEN
          go to 3100
        END IF
c
        medium = med(irl)
3110    CONTINUE
3111      CONTINUE
          compute_tstep = .true.
          eke = eie - rm
          IF ((medium .NE. 0)) THEN
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            RNNE1 = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            IF ((RNNE1.EQ.0.0)) THEN
              RNNE1=1.E-30
            END IF
            DEMFP=MAX(-LOG(RNNE1),1.E-5)
            elke = log(eke)
            Lelke=eke1(MEDIUM)*elke+eke0(MEDIUM)
            IF (( lelec .LT. 0 )) THEN
              sig0 = esig_e(medium)
            ELSE
              sig0 = psig_e(medium)
            END IF
          END IF
3120      CONTINUE
3121        CONTINUE
            IF ((medium .EQ. 0)) THEN
              tstep = vacdst
              ustep = tstep
              tustep = ustep
              callhowfar = .true.
            ELSE
              RHOF=RHOR(IRL)/RHO(MEDIUM)
              sig = sig0
              IF ((sig .LE. 0)) THEN
                tstep = vacdst
                sig0 = 1.E-15
              ELSE
                IF (( compute_tstep )) THEN
                  total_de = demfp/sig
                  fedep = total_de
                  ekef = eke - fedep
                  IF (( ekef .LE. E_array(1,medium) )) THEN
                    tstep = vacdst
                  ELSE
                    elkef = Log(ekef)
                    Lelkef=eke1(MEDIUM)*elkef+eke0(MEDIUM)
                    IF (( lelkef .EQ. lelke )) THEN
                      fedep = 1 - ekef/eke
                      elktmp = 0.5*(elke+elkef+0.25*fedep*fedep*(1+fedep
     *                *(1+0.875*fedep)))
                      lelktmp = lelke
                      IF ((lelec .LT. 0)) THEN
                        dedxmid=ededx1(Lelktmp,MEDIUM)*elktmp+ededx0(Lel
     *                  ktmp,MEDIUM)
                        aux = ededx1(lelktmp,medium)/dedxmid
                      ELSE
                        dedxmid=pdedx1(Lelktmp,MEDIUM)*elktmp+pdedx0(Lel
     *                  ktmp,MEDIUM)
                        aux = pdedx1(lelktmp,medium)/dedxmid
                      END IF
                      aux = aux*(1+2*aux)*(fedep/(2-fedep))**2/6
                      tstep = fedep*eke/dedxmid*(1+aux)
                    ELSE
                      ekei = E_array(lelke,medium)
                      elkei = (lelke - eke0(medium))/eke1(medium)
                      fedep = 1 - ekei/eke
                      elktmp = 0.5*(elke+elkei+0.25*fedep*fedep*(1+fedep
     *                *(1+0.875*fedep)))
                      lelktmp = lelke
                      IF ((lelec .LT. 0)) THEN
                        dedxmid=ededx1(Lelktmp,MEDIUM)*elktmp+ededx0(Lel
     *                  ktmp,MEDIUM)
                        aux = ededx1(lelktmp,medium)/dedxmid
                      ELSE
                        dedxmid=pdedx1(Lelktmp,MEDIUM)*elktmp+pdedx0(Lel
     *                  ktmp,MEDIUM)
                        aux = pdedx1(lelktmp,medium)/dedxmid
                      END IF
                      aux = aux*(1+2*aux)*(fedep/(2-fedep))**2/6
                      tuss = fedep*eke/dedxmid*(1+aux)
                      ekei = E_array(lelkef+1,medium)
                      elkei = (lelkef + 1 - eke0(medium))/eke1(medium)
                      fedep = 1 - ekef/ekei
                      elktmp = 0.5*(elkei+elkef+0.25*fedep*fedep*(1+fede
     *                p*(1+0.875*fedep)))
                      lelktmp = lelkef
                      IF ((lelec .LT. 0)) THEN
                        dedxmid=ededx1(Lelktmp,MEDIUM)*elktmp+ededx0(Lel
     *                  ktmp,MEDIUM)
                        aux = ededx1(lelktmp,medium)/dedxmid
                      ELSE
                        dedxmid=pdedx1(Lelktmp,MEDIUM)*elktmp+pdedx0(Lel
     *                  ktmp,MEDIUM)
                        aux = pdedx1(lelktmp,medium)/dedxmid
                      END IF
                      aux = aux*(1+2*aux)*(fedep/(2-fedep))**2/6
                      tstep = fedep*ekei/dedxmid*(1+aux)
                      tstep=tstep+tuss+ range_ep(qel,lelke,medium)-range
     *                _ep(qel,lelkef+1,medium)
                    END IF
                  END IF
                  total_tstep = tstep
                  compute_tstep = .false.
                END IF
                tstep = total_tstep/rhof
              END IF
              IF ((lelec .LT. 0)) THEN
                dedx0=ededx1(Lelke,MEDIUM)*elke+ededx0(Lelke,MEDIUM)
              ELSE
                dedx0=pdedx1(Lelke,MEDIUM)*elke+pdedx0(Lelke,MEDIUM)
              END IF
              dedx = rhof*dedx0
              tmxs=tmxs1(Lelke,MEDIUM)*elke+tmxs0(Lelke,MEDIUM)
              tmxs = tmxs/rhof
              ekei = E_array(lelke,medium)
              elkei = (lelke - eke0(medium))/eke1(medium)
              fedep = 1 - ekei/eke
              elktmp = 0.5*(elke+elkei+0.25*fedep*fedep*(1+fedep*(1+0.87
     *        5*fedep)))
              lelktmp = lelke
              IF ((lelec .LT. 0)) THEN
                dedxmid=ededx1(Lelktmp,MEDIUM)*elktmp+ededx0(Lelktmp,MED
     *          IUM)
                aux = ededx1(lelktmp,medium)/dedxmid
              ELSE
                dedxmid=pdedx1(Lelktmp,MEDIUM)*elktmp+pdedx0(Lelktmp,MED
     *          IUM)
                aux = pdedx1(lelktmp,medium)/dedxmid
              END IF
              aux = aux*(1+2*aux)*(fedep/(2-fedep))**2/6
              range = fedep*eke/dedxmid*(1+aux)
              range = (range + range_ep(qel,lelke,medium))/rhof
              IF ((.false.)) THEN
                IF (( rng_seed .GT. 24 )) THEN
                  call ranlux(rng_array)
                  rng_seed = 1
                END IF
                rnnotu = rng_array(rng_seed)
                rng_seed = rng_seed + 1
                tmxs = rnnotu*min(tmxs,smaxir(irl))
              ELSE
                tmxs = min(tmxs,smaxir(irl))
              END IF
              tustep = min(tstep,tmxs,range)
c
c mirko - modified hownear - call
c
c              call hownear(tperp,z(np),ir(np))
              call hownear(tperp)
c              
              dnear(np) = tperp
              IF (( i_do_rr(irl) .EQ. 1 .AND. e(np) .LT. e_max_rr(irl) )
     *        ) THEN
                IF ((tperp .GE. range)) THEN
                  idisc = 50 + 49*iq(np)
                  go to 3130
                END IF
              END IF
              blccl = rhof*blcc(medium)
              p2 = eke*(eke+rmt2)
              beta2 = p2/(p2 + rmsq)
              IF (( spin_effects )) THEN
                IF ((lelec .LT. 0)) THEN
                  etap=etae_ms1(Lelke,MEDIUM)*elke+etae_ms0(Lelke,MEDIUM
     *            )
                ELSE
                  etap=etap_ms1(Lelke,MEDIUM)*elke+etap_ms0(Lelke,MEDIUM
     *            )
                END IF
                ms_corr=blcce1(Lelke,MEDIUM)*elke+blcce0(Lelke,MEDIUM)
                blccl = blccl/etap/(1+0.25*etap*xcc(medium)/blcc(medium)
     *          /p2)*ms_corr
              END IF
              ssmfp=beta2/blccl
              skindepth = skindepth_for_bca*ssmfp
              tustep = min(tustep,max(tperp,skindepth))
              count_all_steps = count_all_steps + 1
              IF (((tustep .LE. tperp) .AND. ((.NOT.exact_bca) .OR. (tus
     *        tep .GT. skindepth)))) THEN
                callhowfar = .false.
                domultiple = .false.
                dosingle = .false.
                callmsdist = .true.
                tuss = range - range_ep(qel,lelke,medium)/rhof
                IF (( tuss .GE. tustep )) THEN
                  IF (( lelec .LT. 0 )) THEN
                    dedxmid=ededx1(Lelke,MEDIUM)*elke+ededx0(Lelke,MEDIU
     *              M)
                    aux = ededx1(lelke,medium)/dedxmid
                  ELSE
                    dedxmid=pdedx1(Lelke,MEDIUM)*elke+pdedx0(Lelke,MEDIU
     *              M)
                    aux = pdedx1(lelke,medium)/dedxmid
                  END IF
                  de = dedxmid*tustep
                  fedep = de/eke
                  de = de*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1- 0.2
     *            5*fedep*(2-aux*(4-aux)))))
                ELSE
                  lelktmp = lelke
                  tuss = (range - tustep)*rhof
                  IF (( tuss .LE. 0 )) THEN
                    de = eke
                  ELSE
3141                IF(tuss.GE.range_ep(qel,lelktmp,medium))GO TO 3142
                      lelktmp = lelktmp - 1
                    GO TO 3141
3142                CONTINUE
                    elktmp = (lelktmp+1-eke0(medium))/eke1(medium)
                    eketmp = E_array(lelktmp+1,medium)
                    tuss = range_ep(qel,lelktmp+1,medium) - tuss
                    IF (( lelec .LT. 0 )) THEN
                      dedxmid=ededx1(Lelktmp,MEDIUM)*elktmp+ededx0(Lelkt
     *                mp,MEDIUM)
                      aux = ededx1(lelktmp,medium)/dedxmid
                    ELSE
                      dedxmid=pdedx1(Lelktmp,MEDIUM)*elktmp+pdedx0(Lelkt
     *                mp,MEDIUM)
                      aux = pdedx1(lelktmp,medium)/dedxmid
                    END IF
                    de = dedxmid*tuss
                    fedep = de/eketmp
                    de = de*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1- 0
     *              .25*fedep*(2-aux*(4-aux)))))
                    de = de + eke - eketmp
                  END IF
                END IF
                tvstep = tustep
                IF ((transport_algorithm .EQ. 0)) THEN
                  call msdist_pII (  eke,de,tustep,rhof,medium,qel,spin_
     *            effects, u(np),v(np),w(np),x(np),y(np),z(np),  uscat,v
     *            scat,wscat,xtrans,ytrans,ztrans,ustep )
                ELSE
                  call msdist_pI (  eke,de,tustep,rhof,medium,qel,spin_e
     *            ffects, u(np),v(np),w(np),x(np),y(np),z(np),  uscat,vs
     *            cat,wscat,xtrans,ytrans,ztrans,ustep )
                END IF
              ELSE
                callmsdist = .false.
                IF ((exact_bca)) THEN
                  domultiple = .false.
                  IF (( rng_seed .GT. 24 )) THEN
                    call ranlux(rng_array)
                    rng_seed = 1
                  END IF
                  rnnoss = rng_array(rng_seed)
                  rng_seed = rng_seed + 1
                  lambda = - Log(1 - rnnoss)
                  lambda_max = 0.5*blccl*rm/dedx*(eke/rm+1)**3
                  IF (( lambda .LT. lambda_max )) THEN
                    tuss = lambda * ssmfp * (1 - 0.5*lambda/lambda_max)
                    IF ((tuss .LT. tustep)) THEN
                      tustep = tuss
                      dosingle = .true.
                    ELSE
                      dosingle = .false.
                    END IF
                  ELSE
cale
                    write(6,*) ' lambda > lambda_max: ',lambda,lambda_ma        
     *              x
                    write(6,*) ' eke dedx: ',eke,dedx                           
                    write(6,*) ' medold,medium, blcc: ',medold,medium,
     *                blcc(medium) 
                    write(6,*) ' np, charge,energy: ',np,iq(np),e(np)                           
                    write(6,*) ' u,v,w: ',u(np),v(np),w(np),
     *               'x,y,z:',x(np),y(np),z(np)                          
                    write(6,*) ' ir(np): ',ir(np),'rhof:',rhof                          
                    write(6,*)'dedx0:',dedx0
                    write(6,*)'Lelke:',Lelke,'eke1(medium)',eke1(medium)
                    write(6,*)'eke0(medium)',eke0(medium)
                    write(6,*)'eke',eke,'elke',elke

cale


             
                    dosingle = .false.
                  END IF
                  ustep = tustep
                ELSE
                  dosingle = .false.
                  domultiple = .true.
                  ekems = eke - 0.5*tustep*dedx
                  p2 = ekems*(ekems+rmt2)
                  beta2 = p2/(p2 + rmsq)
                  blccl = blcc(medium)
                  chia2 = xcc(medium)/(4*blccl*p2)
                  xi = 0.5*xcc(medium)/p2/beta2*tustep
                  IF (( spin_effects )) THEN
                    elkems = Log(ekems)
                    Lelkems=eke1(MEDIUM)*elkems+eke0(MEDIUM)
                    IF ((lelec .LT. 0)) THEN
                      etap=etae_ms1(Lelkems,MEDIUM)*elkems+etae_ms0(Lelk
     *                ems,MEDIUM)
                      xi_corr=q1ce_ms1(Lelkems,MEDIUM)*elkems+q1ce_ms0(L
     *                elkems,MEDIUM)
                    ELSE
                      etap=etap_ms1(Lelkems,MEDIUM)*elkems+etap_ms0(Lelk
     *                ems,MEDIUM)
                      xi_corr=q1cp_ms1(Lelkems,MEDIUM)*elkems+q1cp_ms0(L
     *                elkems,MEDIUM)
                    END IF
                    chia2 = chia2*etap
                    xi = xi*xi_corr
                    ms_corr=blcce1(Lelkems,MEDIUM)*elkems+blcce0(Lelkems
     *              ,MEDIUM)
                    blccl = blccl*ms_corr
                  ELSE
                    xi_corr = 1
                    etap = 1
                  END IF
                  xi = xi*(Log(1+1./chia2)-1/(1+chia2))
                  IF (( xi .LT. 0.1 )) THEN
                    ustep = tustep*(1 - xi*(0.5 - xi*0.166667))
                  ELSE
                    ustep = tustep*(1 - Exp(-xi))/xi
                  END IF
                END IF
                IF ((ustep .LT. tperp)) THEN
                  callhowfar = .false.
                ELSE
                  callhowfar = .true.
                END IF
              END IF
            END IF
            irnew = ir(np)
            idisc = 0
            ustep0 = ustep
            IF ((callhowfar .OR. wt(np) .LE. 0)) THEN
              call howfar
            END IF
            IF ((idisc .GT. 0)) THEN
              go to 3130
            END IF
            IF ((ustep .LE. 0)) THEN
              IF ((ustep .LT. -1e-4)) THEN
                ierust = ierust + 1
                WRITE(6,3150)ierust,ustep,dedx, ir(np),irnew,irold,x(np)
     *          ,y(np),z(np)
3150            FORMAT(i6,' Negative ustep = ',2e14.6, ' ir,irnew,irold         
     *=',3i4,'x,y,z =',4e10.3)                                                  
                IF ((ierust .GT. 1000)) THEN
                  WRITE(6,3160)
3160              FORMAT(////' Called exit---too many ustep errors'///)         
                  stop
                END IF
              END IF
              ustep = 0
            END IF
            IF ((ustep .EQ. 0 .OR. medium .EQ. 0)) THEN
              IF ((ustep .NE. 0)) THEN
                vstep = ustep
                tvstep = vstep
                edep = pzero
                e_range = vacdst
                IARG=0
                IF ((IAUSFL(IARG+1).NE.0)) THEN
                  CALL AUSGAB(IARG)
                END IF
                x(np) = x(np) + u(np)*vstep
                y(np) = y(np) + v(np)*vstep
                z(np) = z(np) + w(np)*vstep
                dnear(np) = dnear(np) - vstep
                irold = ir(np)
              END IF
              ir(np) = irnew
              irl = irnew
              medium = med(irl)
              IF ((ustep .NE. 0)) THEN
                IARG=5
                IF ((IAUSFL(IARG+1).NE.0)) THEN
                  CALL AUSGAB(IARG)
                END IF
              END IF
              IF ((eie .LE. ecut(irl))) THEN
                go to 3100
              END IF
              IF ((ustep .NE. 0 .AND. idisc .LT. 0)) THEN
                go to 3130
              END IF
              GO TO 3111
            END IF
            vstep = ustep
            IF ((callhowfar)) THEN
              IF ((exact_bca)) THEN
                tvstep = vstep
                IF ((tvstep .NE. tustep)) THEN
                  dosingle = .false.
                END IF
              ELSE
                IF (( vstep .LT. ustep0 )) THEN
                  ekems = eke - 0.5*tustep*vstep/ustep0*dedx
                  p2 = ekems*(ekems+rmt2)
                  beta2 = p2/(p2 + rmsq)
                  blccl = blcc(medium)
                  chia2 = xcc(medium)/(4*blccl*p2)
                  xi = 0.5*xcc(medium)/p2/beta2*vstep
                  IF (( spin_effects )) THEN
                    elkems = Log(ekems)
                    Lelkems=eke1(MEDIUM)*elkems+eke0(MEDIUM)
                    IF ((lelec .LT. 0)) THEN
                      etap=etae_ms1(Lelkems,MEDIUM)*elkems+etae_ms0(Lelk
     *                ems,MEDIUM)
                      xi_corr=q1ce_ms1(Lelkems,MEDIUM)*elkems+q1ce_ms0(L
     *                elkems,MEDIUM)
                    ELSE
                      etap=etap_ms1(Lelkems,MEDIUM)*elkems+etap_ms0(Lelk
     *                ems,MEDIUM)
                      xi_corr=q1cp_ms1(Lelkems,MEDIUM)*elkems+q1cp_ms0(L
     *                elkems,MEDIUM)
                    END IF
                    chia2 = chia2*etap
                    xi = xi*xi_corr
                    ms_corr=blcce1(Lelkems,MEDIUM)*elkems+blcce0(Lelkems
     *              ,MEDIUM)
                    blccl = blccl*ms_corr
                  ELSE
                    xi_corr = 1
                    etap = 1
                  END IF
                  xi = xi*(Log(1+1./chia2)-1/(1+chia2))
                  IF (( xi .LT. 0.1 )) THEN
                    tvstep = vstep*(1 + xi*(0.5 + xi*0.333333))
                  ELSE
                    IF (( xi .LT. 0.999999 )) THEN
                      tvstep = -vstep*Log(1 - xi)/xi
                    ELSE
                      write(6,*) ' Stoped in SET-TVSTEP because xi > 1!         
     *'                                                                         
                      write(6,*) ' Medium: ',medium                             
                      write(6,*) ' Initial energy: ',eke                        
                      write(6,*) ' Average step energy: ',ekems                 
                      write(6,*) ' tustep: ',tustep                             
                      write(6,*) ' ustep0: ',ustep0                             
                      write(6,*) ' vstep:  ',vstep                              
                      write(6,*) ' ==> xi = ',xi                                
                      stop
                    END IF
                  END IF
                ELSE
                  tvstep = tustep
                END IF
              END IF
              tuss = range - range_ep(qel,lelke,medium)/rhof
              IF (( tuss .GE. tvstep )) THEN
                IF (( lelec .LT. 0 )) THEN
                  dedxmid=ededx1(Lelke,MEDIUM)*elke+ededx0(Lelke,MEDIUM)
                  aux = ededx1(lelke,medium)/dedxmid
                ELSE
                  dedxmid=pdedx1(Lelke,MEDIUM)*elke+pdedx0(Lelke,MEDIUM)
                  aux = pdedx1(lelke,medium)/dedxmid
                END IF
                de = dedxmid*tvstep
                fedep = de/eke
                de = de*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1- 0.25*
     *          fedep*(2-aux*(4-aux)))))
              ELSE
                lelktmp = lelke
                tuss = (range - tvstep)*rhof
                IF (( tuss .LE. 0 )) THEN
                  de = eke
                ELSE
3171              IF(tuss.GE.range_ep(qel,lelktmp,medium))GO TO 3172
                    lelktmp = lelktmp - 1
                  GO TO 3171
3172              CONTINUE
                  elktmp = (lelktmp+1-eke0(medium))/eke1(medium)
                  eketmp = E_array(lelktmp+1,medium)
                  tuss = range_ep(qel,lelktmp+1,medium) - tuss
                  IF (( lelec .LT. 0 )) THEN
                    dedxmid=ededx1(Lelktmp,MEDIUM)*elktmp+ededx0(Lelktmp
     *              ,MEDIUM)
                    aux = ededx1(lelktmp,medium)/dedxmid
                  ELSE
                    dedxmid=pdedx1(Lelktmp,MEDIUM)*elktmp+pdedx0(Lelktmp
     *              ,MEDIUM)
                    aux = pdedx1(lelktmp,medium)/dedxmid
                  END IF
                  de = dedxmid*tuss
                  fedep = de/eketmp
                  de = de*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1- 0.2
     *            5*fedep*(2-aux*(4-aux)))))
                  de = de + eke - eketmp
                END IF
              END IF
            ELSE
              tvstep = tustep
              IF (( .NOT.callmsdist )) THEN
                tuss = range - range_ep(qel,lelke,medium)/rhof
                IF (( tuss .GE. tvstep )) THEN
                  IF (( lelec .LT. 0 )) THEN
                    dedxmid=ededx1(Lelke,MEDIUM)*elke+ededx0(Lelke,MEDIU
     *              M)
                    aux = ededx1(lelke,medium)/dedxmid
                  ELSE
                    dedxmid=pdedx1(Lelke,MEDIUM)*elke+pdedx0(Lelke,MEDIU
     *              M)
                    aux = pdedx1(lelke,medium)/dedxmid
                  END IF
                  de = dedxmid*tvstep
                  fedep = de/eke
                  de = de*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1- 0.2
     *            5*fedep*(2-aux*(4-aux)))))
                ELSE
                  lelktmp = lelke
                  tuss = (range - tvstep)*rhof
                  IF (( tuss .LE. 0 )) THEN
                    de = eke
                  ELSE
3181                IF(tuss.GE.range_ep(qel,lelktmp,medium))GO TO 3182
                      lelktmp = lelktmp - 1
                    GO TO 3181
3182                CONTINUE
                    elktmp = (lelktmp+1-eke0(medium))/eke1(medium)
                    eketmp = E_array(lelktmp+1,medium)
                    tuss = range_ep(qel,lelktmp+1,medium) - tuss
                    IF (( lelec .LT. 0 )) THEN
                      dedxmid=ededx1(Lelktmp,MEDIUM)*elktmp+ededx0(Lelkt
     *                mp,MEDIUM)
                      aux = ededx1(lelktmp,medium)/dedxmid
                    ELSE
                      dedxmid=pdedx1(Lelktmp,MEDIUM)*elktmp+pdedx0(Lelkt
     *                mp,MEDIUM)
                      aux = pdedx1(lelktmp,medium)/dedxmid
                    END IF
                    de = dedxmid*tuss
                    fedep = de/eketmp
                    de = de*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1- 0
     *              .25*fedep*(2-aux*(4-aux)))))
                    de = de + eke - eketmp
                  END IF
                END IF
              END IF
            END IF
            save_de = de
            edep = de
            ekef = eke - de
            eold = eie
            enew = eold - de
            IF (( .NOT.callmsdist )) THEN
              IF (( domultiple )) THEN
                lambda = blccl*tvstep/beta2/etap/(1+chia2)
                xi = xi/xi_corr
                findindex = .true.
                spin_index = .true.
                call mscat(lambda,chia2,xi,elkems,beta2,qel,medium, spin
     *          _effects,findindex,spin_index, costhe,sinthe)
              ELSE
                IF ((dosingle)) THEN
                  ekems = Max(ekef,ecut(irl)-rm)
                  p2 = ekems*(ekems + rmt2)
                  beta2 = p2/(p2 + rmsq)
                  chia2 = xcc(medium)/(4*blcc(medium)*p2)
                  IF (( spin_effects )) THEN
                    elkems = Log(ekems)
                    Lelkems=eke1(MEDIUM)*elkems+eke0(MEDIUM)
                    IF ((lelec .LT. 0)) THEN
                      etap=etae_ms1(Lelkems,MEDIUM)*elkems+etae_ms0(Lelk
     *                ems,MEDIUM)
                    ELSE
                      etap=etap_ms1(Lelkems,MEDIUM)*elkems+etap_ms0(Lelk
     *                ems,MEDIUM)
                    END IF
                    chia2 = chia2*etap
                  END IF
                  call sscat(chia2,elkems,beta2,qel,medium, spin_effects
     *            ,costhe,sinthe)
                ELSE
                  theta = 0
                  sinthe = 0
                  costhe = 1
                END IF
              END IF
            END IF
            e_range = range
            IARG=0
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
            IF (( callmsdist )) THEN
              u(np) = uscat
              v(np) = vscat
              w(np) = wscat
              x(np) = xtrans
              y(np) = ytrans
              z(np) = ztrans
            ELSE
              x(np) = x(np) + u(np)*vstep
              y(np) = y(np) + v(np)*vstep
              z(np) = z(np) + w(np)*vstep
              IF (( domultiple .OR. dosingle )) THEN
                call uphi(2,1)
              END IF
            END IF
            dnear(np) = dnear(np) - vstep
            irold = ir(np)
            peie = peie - edep
            eie = peie
            e(np) = peie
            IF ((eie .LE. ecut(irl))) THEN
cmirko
c          if ((abs(x(np)).gt.100.).or.(abs(y(np)).gt.100.).or.
c     a (abs(z(np)).gt.100.)) write (*,*) 'pippo c overflow', callmsdist
cmirko
              go to 3100
            END IF
            medold = medium
            IF ((medium .NE. 0)) THEN
              ekeold = eke
              eke = eie - rm
              elke = log(eke)
              Lelke=eke1(MEDIUM)*elke+eke0(MEDIUM)
            END IF
            IF ((irnew .NE. irold)) THEN
              ir(np) = irnew
              irl = irnew
              medium = med(irl)
            END IF
            IARG=5
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
            IF ((eie .LE. ecut(irl))) THEN
              go to 3100
            END IF
            IF ((idisc .LT. 0)) THEN
              go to 3130
            END IF
            IF((medium .NE. medold))GO TO 3111
            demfp = demfp - save_de*sig
            total_de = total_de - save_de
            total_tstep = total_tstep - tvstep*rhof
            IF (( total_tstep .LT. 1e-9 )) THEN
              demfp = 0
            END IF
            IF(((demfp .LT. 1.E-5)))GO TO3122
          GO TO 3121
3122      CONTINUE
          IF ((lelec .LT. 0)) THEN
            sigf=esig1(Lelke,MEDIUM)*elke+esig0(Lelke,MEDIUM)
            dedx0=ededx1(Lelke,MEDIUM)*elke+ededx0(Lelke,MEDIUM)
            sigf = sigf/dedx0
          ELSE
            sigf=psig1(Lelke,MEDIUM)*elke+psig0(Lelke,MEDIUM)
            dedx0=pdedx1(Lelke,MEDIUM)*elke+pdedx0(Lelke,MEDIUM)
            sigf = sigf/dedx0
          END IF
          sigratio = sigf/sig0
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          rfict = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF(((rfict .LE. sigratio)))GO TO3112
        GO TO 3111
3112    CONTINUE
        IF ((lelec .LT. 0)) THEN
          ebr1=ebr11(Lelke,MEDIUM)*elke+ebr10(Lelke,MEDIUM)
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          rnno24 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
c mirko
c                write (*,*) 'ciao 2a', u(np),v(np),w(np),rnno24,ebr1
c     a,e(np),thmoll(medium),iausfl(9),iausfl(10)
c
          IF ((rnno24 .LE. ebr1)) THEN
            go to 3190
          ELSE
            IF ((e(np) .LE. thmoll(medium))) THEN
              IF ((ebr1 .LE. 0)) THEN
c mirko
c                write (*,*) 'ciao 1'
c
                go to 3090
              END IF
              go to 3190
            END IF
            IARG=8
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
            call moller
            IARG=9
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
          END IF
c mirko
c                write (*,*) 'ciao 2b', u(np),v(np),w(np)
c
          go to 3090
        END IF
        pbr1=pbr11(Lelke,MEDIUM)*elke+pbr10(Lelke,MEDIUM)
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        rnno25 = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        IF ((rnno25 .LT. pbr1)) THEN
          go to 3190
        END IF
        pbr2=pbr21(Lelke,MEDIUM)*elke+pbr20(Lelke,MEDIUM)
        IF ((rnno25 .LT. pbr2)) THEN
          IARG=10
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          call bhabha
          IARG=11
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
        ELSE
          IARG=12
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          call annih
          IARG=13
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          GO TO 3092
        END IF
      GO TO 3091
3092  CONTINUE
      return
3190  IARG=6
      IF ((IAUSFL(IARG+1).NE.0)) THEN
        CALL AUSGAB(IARG)
      END IF
      call brems
      IARG=7
      IF ((IAUSFL(IARG+1).NE.0)) THEN
        CALL AUSGAB(IARG)
      END IF
      IF ((iq(np) .EQ. 0)) THEN
        return
      ELSE
c mirko
c                write (*,*) 'ciao 3'
c
        go to 3090
      END IF
3100  IF ((eie .GT. ae(medium))) THEN
        idr = 1
        IF ((lelec .LT. 0)) THEN
          edep = e(np) - prm
        ELSE
          EDEP=PEIE-PRM
        END IF
      ELSE
        idr = 2
        edep = e(np) - prm
      END IF
      IARG=idr
      IF ((IAUSFL(IARG+1).NE.0)) THEN
        CALL AUSGAB(IARG)
      END IF
3200  CONTINUE
      IF ((lelec .GT. 0)) THEN
        IF ((edep .LT. peie)) THEN
          NPold = NP
          IF (( nbr_split .GT. 1 )) THEN
            wt(np) = wt(np)/nbr_split
          END IF
            DO 3211 ibr=1,nbr_split
            IF (( np+1 .GT. 50 )) THEN
              WRITE(6,3220)np+1
3220          FORMAT(//' Stack overflow in ANNIH at rest! np = ',i6, ' I        
     *ncrease $MXSTACK and try again'//)                                        
              stop
            END IF
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            costhe = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            costhe = 2*costhe-1
            sinthe = sqrt(max(0.0,(1-costhe)*(1+costhe)))
3231        CONTINUE
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
              IF(rhophi2.LE.1)GO TO3232
            GO TO 3231
3232        CONTINUE
            rhophi2 = 1/rhophi2
            cphi = (xphi2 - yphi2)*rhophi2
            sphi = 2*xphi*yphi*rhophi2
            e(np) = prm
            iq(np) = 0
            X(np)=X(npold)
            Y(np)=Y(npold)
            Z(np)=Z(npold)
            IR(np)=IR(npold)
            WT(np)=WT(npold)
            DNEAR(np)=DNEAR(npold)
            LATCH(np)=LATCH(npold)
            u(np) = sinthe*cphi
            v(np) = sinthe*sphi
            w(np) = costhe
            np = np+1
            e(np) = prm
            iq(np) = 0
            X(np)=X(npold)
            Y(np)=Y(npold)
            Z(np)=Z(npold)
            IR(np)=IR(npold)
            WT(np)=WT(npold)
            DNEAR(np)=DNEAR(npold)
            LATCH(np)=LATCH(npold)
            u(np) = -u(np-1)
            v(np) = -v(np-1)
            w(np) = -w(np-1)
            np = np+1
3211      CONTINUE
3212      CONTINUE
          np = np-1
          IARG=14
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          return
        END IF
      END IF
      np = np - 1
      ircode = 2
      return
3130  idisc = abs(idisc)
      IF (((lelec .LT. 0) .OR. (idisc .EQ. 99))) THEN
        edep = e(np) - prm
      ELSE
        edep = e(np) + prm
      END IF
      IARG=3
      IF ((IAUSFL(IARG+1).NE.0)) THEN
        CALL AUSGAB(IARG)
      END IF
      IF((idisc .EQ. 99))goto 3200
      ircode = 2
      np = np - 1
      return
      end
