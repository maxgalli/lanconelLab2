      SUBROUTINE PHOTO
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      real*4 EELEC,  BETA,  GAMMA,  ALPHA,  RATIO,  RNPHT,  FKAPPA, XI,
     * SINTH2, RNPHT2
      DOUBLE PRECISION PEIG
      real*4 BR,  sigma,  aux,aux1,  probs(50),  sigtot,  e_vac,  rnno_R
     *R
      integer*4 IARG,  iZ,   irl,  ints(50),  j,ip,  n_warning,  k
      logical do_relax
      save n_warning
      data n_warning / 0 /
      NPold = NP
      PEIG=E(NP)
      irl = ir(np)
      IF (( peig .LT. edge_energies(2,1) )) THEN
        IF (( n_warning .LT. 100 )) THEN
          n_warning = n_warning + 1
          write(6,*) ' Subroutine PHOTO called with E = ',peig, ' which         
     *is below the current min. energy of 1 keV! '                              
          write(6,*) ' Converting now this photon to an electron, '             
          write(6,*) ' but you should check your code! '                        
        END IF
        iq(np) = -1
        e(np) = peig + prm
        return
      END IF
      iZ = iedgfl(irl)
      do_relax = .false.
      edep = pzero
      IF (( iedgfl(irl) .NE. 0 )) THEN
        IF (( nne(medium) .EQ. 1 )) THEN
          iZ = int( zelem(medium,1) + 0.5 )
            DO 4541 j=1,edge_number(iZ)
            IF((peig .GE. edge_energies(j,iZ)))GO TO4542
4541      CONTINUE
4542      CONTINUE
        ELSE
          aux = peig*peig
          aux1 = aux*peig
          aux = aux*Sqrt(peig)
          sigtot = 0
            DO 4551 k=1,nne(medium)
            iZ = int( zelem(medium,k) + 0.5 )
            IF (( iZ .LT. 1 .OR. iZ .GT. 100 )) THEN
              write(6,*) ' Error in PHOTO: '                                    
              write(6,*) '   Atomic number of element ',k, ' in medium '        
     *        ,medium,' is not between 1 and ', 100                             
              stop
            END IF
            IF (( peig .GT. edge_energies(1,iZ) )) THEN
              j = 1
              sigma = (edge_a(1,iZ) + edge_b(1,iZ)/peig + edge_c(1,iZ)/a
     *        ux + edge_d(1,iZ)/aux1)/peig
            ELSE
                DO 4561 j=2,edge_number(iZ)
                IF((peig .GE. edge_energies(j,iZ)))GO TO4562
4561          CONTINUE
4562          CONTINUE
              sigma = edge_a(j,iZ) + gle*(edge_b(j,iZ) + gle*(edge_c(j,i
     *        Z) + gle*edge_d(j,iZ) ))
              sigma = Exp(sigma)
            END IF
            sigma = sigma * pz(medium,k)
            sigtot = sigtot + sigma
            probs(k) = sigma
            ints(k) = j
4551      CONTINUE
4552      CONTINUE
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          br = rng_array(rng_seed)
          rng_seed = rng_seed + 1
            DO 4571 k=1,nne(medium)
            sigma = probs(k)/sigtot
            br = br - sigma
            IF((br .LE. 0))GO TO4572
4571      CONTINUE
4572      CONTINUE
          iZ = int( zelem(medium,k) + 0.5 )
          j = ints(k)
        END IF
        IF (( peig .LE. binding_energies(6,iZ) )) THEN
          edep = peig
          e(np) = pzero
          wt(np) = 0
        ELSE
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          br = rng_array(rng_seed)
          rng_seed = rng_seed + 1
            DO 4581 k=1,5
            IF (( peig .GT. binding_energies(k,iZ) )) THEN
              IF((br .LT. interaction_prob(k,iZ)))GO TO4582
              br = (br - interaction_prob(k,iZ))/(1-interaction_prob(k,i
     *        Z))
            END IF
4581      CONTINUE
4582      CONTINUE
          e_vac = binding_energies(k,iZ)
          e(np) = peig - e_vac + prm
          do_relax = .true.
          iq(np) = -1
        END IF
      ELSE
        e(np) = peig + prm
        iq(np) = -1
      END IF
      IF (( iq(np) .EQ. -1 )) THEN
        IF ((IPHTER(IR(NP)).EQ.1)) THEN
          EELEC=E(NP)
          IF ((EELEC.GT.ECUT(IR(NP)))) THEN
            BETA=SQRT((EELEC-RM)*(EELEC+RM))/EELEC
            GAMMA=EELEC/RM
            ALPHA=0.5*GAMMA-0.5+1./GAMMA
            RATIO=BETA/ALPHA
4591        CONTINUE
              IF (( rng_seed .GT. 24 )) THEN
                call ranlux(rng_array)
                rng_seed = 1
              END IF
              RNPHT = rng_array(rng_seed)
              rng_seed = rng_seed + 1
              RNPHT=2.*RNPHT-1.
              IF ((RATIO.LE.0.2)) THEN
                FKAPPA=RNPHT+0.5*RATIO*(1.-RNPHT)*(1.+RNPHT)
                COSTHE=(BETA+FKAPPA)/(1.+BETA*FKAPPA)
                XI=1./(1.-BETA*COSTHE)
              ELSE
                XI=GAMMA*GAMMA*(1.+ALPHA*(SQRT(1.+RATIO*(2.*RNPHT+RATIO)
     *          )-1.))
                COSTHE=(1.-1./XI)/BETA
              END IF
              SINTH2=MAX(0.,(1.-COSTHE)*(1.+COSTHE))
              IF (( rng_seed .GT. 24 )) THEN
                call ranlux(rng_array)
                rng_seed = 1
              END IF
              RNPHT2 = rng_array(rng_seed)
              rng_seed = rng_seed + 1
              IF(RNPHT2.LE.0.5*(1.+GAMMA)*SINTH2*XI/GAMMA)GO TO4592
            GO TO 4591
4592        CONTINUE
            SINTHE=SQRT(SINTH2)
            CALL UPHI(2,1)
          END IF
        END IF
      END IF
      IF (( do_relax )) THEN
        call relax(e_vac,k,iZ)
      END IF
      IF (( EDEP .GT. 0 )) THEN
        IARG=4
        IF ((IAUSFL(IARG+1).NE.0)) THEN
          CALL AUSGAB(IARG)
        END IF
      END IF
      i_survived_RR = 0
      IF (( i_play_RR .EQ. 1 )) THEN
        IF (( prob_RR .LE. 0 )) THEN
          IF (( n_RR_warning .LT. 50 )) THEN
            n_RR_warning = n_RR_warning + 1
            WRITE(6,4600)prob_RR
4600        FORMAT('**** Warning, attempt to play Roussian Roulette with        
     * prob_RR<=0! ',g14.6)                                                     
          END IF
        ELSE
          ip = NPold
4611      CONTINUE
            IF (( iq(ip) .NE. 0 )) THEN
              IF (( rng_seed .GT. 24 )) THEN
                call ranlux(rng_array)
                rng_seed = 1
              END IF
              rnno_RR = rng_array(rng_seed)
              rng_seed = rng_seed + 1
              IF (( rnno_RR .LT. prob_RR )) THEN
                wt(ip) = wt(ip)/prob_RR
                ip = ip + 1
              ELSE
                i_survived_RR = i_survived_RR + 1
                IF ((ip .LT. np)) THEN
                  e(ip) = e(np)
                  iq(ip) = iq(np)
                  wt(ip) = wt(np)
                  u(ip) = u(np)
                  v(ip) = v(np)
                  w(ip) = w(np)
                END IF
                np = np-1
              END IF
            ELSE
              ip = ip+1
            END IF
            IF(((ip .GT. np)))GO TO4612
          GO TO 4611
4612      CONTINUE
          IF (( np .EQ. 0 )) THEN
            np = 1
            e(np) = 0
            iq(np) = 0
            wt(np) = 0
          END IF
        END IF
      END IF
      return
      end
