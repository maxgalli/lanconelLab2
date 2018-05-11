      SUBROUTINE RELAX(energy,n,iZ)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      integer*4 n,iZ
      real*4 energy
      integer*4 vac_array(50),  n_vac,  shell
      integer*4 final,finala,  final1,final2,   iql,  irl
      integer*4 first_transition(5), last_transition(5)
      integer*4 final_state(39)
      integer*4 k, np_old, ip, iarg
      real*4 e_array_orig(50),Ei,Ef,Ex,eta,e_check,e_cut,ekcut,pkcut,
     *elcut
      real*4 xphi,yphi,xphi2,yphi2,rhophi2, cphi,sphi
      data first_transition/1,20,27,33,38/
      data last_transition/19,26,32,37,39/
      data final_state/ 3,4,5,6,  202,302,402,404,403,303,  502,503,504,
     *602,603,604,  505,605,606,  13,14,  5,6,  505,605,606,  14,  5,6,
     * 505,605,606,  5,6,  505,605,606,  6,  606/
      save first_transition,last_transition,final_state
      IF (( n .LT. 1 .OR. n .GT. 6 )) THEN
        return
      END IF
      irl = ir(np)
      ekcut = ecut(irl)-rm
      pkcut = pcut(irl)
      e_cut = Min(ekcut,pkcut)
      e_cut = Max(0.001,e_cut)
      IF (( energy .LE. e_cut )) THEN
        edep = edep + energy
        return
      END IF
      n_vac = 1
      vac_array(n_vac) = n
      np_old = np
      e_check = 0
      e_array_orig(n_vac) = energy
4620  CONTINUE
4621    CONTINUE
        shell = vac_array(n_vac)
        Ei = e_array_orig(n_vac)
        n_vac = n_vac - 1
        IF (( Ei .LE. e_cut )) THEN
          edep = edep + Ei
          IF((n_vac .GT. 0))goto 4620
          GO TO4622
        END IF
        IF (( shell .EQ. 6 )) THEN
          IF (( Ei .GT. e_cut )) THEN
            np = np + 1
            IF (( np .GT. 50 )) THEN
              write(6,'(a)') ' *****************************************        
     ***********'                                                               
              write(6,'(a,a)') ' In subroutine ','RELAX',' stack size ex        
     *ceeded! '                                                                 
              write(6,'(a,i6,a,i6)') ' $MXSTACK = ',50,' np = ',np              
              write(6,'(a)') ' Increase $MXSTACK and try again '                
              write(6,'(a)') ' Terminating execution '                          
              write(6,'(a)') ' *****************************************        
     ***********'                                                               
              stop
            END IF
            e(np) = Ei + prm
            iq(np) = -1
            X(np)=X(np_old)
            Y(np)=Y(np_old)
            Z(np)=Z(np_old)
            IR(np)=IR(np_old)
            WT(np)=WT(np_old)
            DNEAR(np)=DNEAR(np_old)
            LATCH(np)=LATCH(np_old)
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            eta = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            eta = 2*eta - 1
            w(np) = eta
            eta = (1-eta)*(1+eta)
            IF (( eta .GT. 1e-20 )) THEN
              eta = Sqrt(eta)
4631          CONTINUE
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
                IF(rhophi2.LE.1)GO TO4632
              GO TO 4631
4632          CONTINUE
              rhophi2 = 1/rhophi2
              cphi = (xphi2 - yphi2)*rhophi2
              sphi = 2*xphi*yphi*rhophi2
              u(np) = eta*cphi
              v(np) = eta*sphi
            ELSE
              u(np) = 0
              v(np) = 0
              w(np) = 1
            END IF
            IARG=27
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
          ELSE
            edep = edep + Ei
          END IF
          IF((n_vac .GT. 0))goto 4620
          GO TO4622
        END IF
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        eta = rng_array(rng_seed)
        rng_seed = rng_seed + 1
          DO 4641 k=first_transition(shell),last_transition(shell)-1
          eta = eta - relaxation_prob(k,iZ)
          IF((eta .LE. 0))GO TO4642
4641    CONTINUE
4642    CONTINUE
        final = final_state(k)
        finala = final
        IF (( final .LT. 100 )) THEN
          IF (( final .LT. 10 )) THEN
            iql = 0
            elcut = pkcut
          ELSE
            final = final - 10
            iql = -1
            elcut = ekcut
          END IF
          Ef = binding_energies(final,iZ)
          Ex = Ei - Ef
          n_vac = n_vac + 1
          vac_array(n_vac) = final
          e_array_orig(n_vac) = Ef
        ELSE
          final1 = final/100
          final2 = final - final1*100
          n_vac = n_vac + 1
          vac_array(n_vac) = final1
          e_array_orig(n_vac) = binding_energies(final1,iZ)
          n_vac = n_vac + 1
          vac_array(n_vac) = final2
          e_array_orig(n_vac) = binding_energies(final2,iZ)
          iql = -1
          Ex = Ei - e_array_orig(n_vac) - e_array_orig(n_vac-1)
          elcut = ekcut
        END IF
        IF (( Ex .LE. elcut )) THEN
          edep = edep + Ex
        ELSE
          np = np + 1
          IF (( np .GT. 50 )) THEN
            write(6,'(a)') ' *******************************************        
     *********'                                                                 
            write(6,'(a,a)') ' In subroutine ','RELAX',' stack size exce        
     *eded! '                                                                   
            write(6,'(a,i6,a,i6)') ' $MXSTACK = ',50,' np = ',np                
            write(6,'(a)') ' Increase $MXSTACK and try again '                  
            write(6,'(a)') ' Terminating execution '                            
            write(6,'(a)') ' *******************************************        
     *********'                                                                 
            stop
          END IF
          iq(np) = iql
          IF (( iql .EQ. 0 )) THEN
            e(np) = Ex
          ELSE
            e(np) = Ex + rm
          END IF
          X(np)=X(np_old)
          Y(np)=Y(np_old)
          Z(np)=Z(np_old)
          IR(np)=IR(np_old)
          WT(np)=WT(np_old)
          DNEAR(np)=DNEAR(np_old)
          LATCH(np)=LATCH(np_old)
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          eta = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          eta = 2*eta - 1
          w(np) = eta
          eta = (1-eta)*(1+eta)
          IF (( eta .GT. 1e-20 )) THEN
            eta = Sqrt(eta)
4651        CONTINUE
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
              IF(rhophi2.LE.1)GO TO4652
            GO TO 4651
4652        CONTINUE
            rhophi2 = 1/rhophi2
            cphi = (xphi2 - yphi2)*rhophi2
            sphi = 2*xphi*yphi*rhophi2
            u(np) = eta*cphi
            v(np) = eta*sphi
          ELSE
            u(np) = 0
            v(np) = 0
            w(np) = 1
          END IF
          IF (( finala .LT. 10 )) THEN
            IARG=25
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
          ELSE IF(( finala .LT. 100 )) THEN
            IARG=26
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
          ELSE
            IARG=27
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
          END IF
        END IF
      GO TO 4621
4622  CONTINUE
      return
      end
