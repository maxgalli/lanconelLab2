      SUBROUTINE COMPT
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      DOUBLE PRECISION PEIG,  PESG,  PESE
      real*4 ko,  broi,  broi2,  bro,  bro1,  alph1,  alph2,  alpha,  rn
     *no15,rnno16,rnno17,rnno18,rnno19,  br,  temp,  rejf3,  rejmax,  Uj
     *,  Jo,  br2,  fpz,fpz1, qc,  qc2,  af,  Fmax,  frej,  eta_incoh, e
     *ta,  aux,aux1,aux2,aux3,aux4,  pzmax,  pz_orig,  pz2,  rnno_RR
      integer*4 irl,  i,  j,  iarg,  ip
      NPold = NP
      peig=E(NP)
      ko = peig/rm
      broi = 1 + 2*ko
      irl = ir(np)
      IF (( ibcmp(irl) .EQ. 1 )) THEN
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        rnno17 = rng_array(rng_seed)
        rng_seed = rng_seed + 1
          DO 2911 i=1,n_shell(medium)
          rnno17 = rnno17 - eno_array(i,medium)
          IF((rnno17 .LE. 0))GO TO2912
2911    CONTINUE
2912    CONTINUE
        j = shell_array(i,medium)
        Uj = be_array(j)
        IF (( ko .LE. Uj )) THEN
          goto 2920
        END IF
        Jo = Jo_array(j)
      END IF
2930  CONTINUE
      IF (( ko .GT. 2 )) THEN
        broi2 = broi*broi
        alph1 = Log(broi)
        alph2 = ko*(broi+1)/broi2
        alpha = alph1/(alph1+alph2)
2941    CONTINUE
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          rnno15 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          rnno16 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF (( rnno15 .LT. alpha )) THEN
            br = Exp(alph1*rnno16)/broi
          ELSE
            br = Sqrt(rnno16 + (1-rnno16)/broi2)
          END IF
          temp = (1-br)/ko/br
          sinthe = Max(0.,temp*(2-temp))
          rejf3 = 1 - br*sinthe/(1+br*br)
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          rnno19 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF((rnno19.le.rejf3))GO TO2942
        GO TO 2941
2942    CONTINUE
      ELSE
        bro = 1./broi
        bro1 = 1 - bro
        rejmax = broi + bro
2951    CONTINUE
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          rnno15 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          rnno16 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          br = bro + bro1*rnno15
          temp = (1-br)/ko/br
          sinthe = Max(0.,temp*(2-temp))
          rejf3 = (br + 1./br - sinthe)/rejmax
          IF((rnno16.le.rejf3))GO TO2952
        GO TO 2951
2952    CONTINUE
      END IF
      IF ((br .LT. 1./broi .OR. br .GT. 1)) THEN
        IF (( br .LT. 0.99999/broi .OR. br .GT. 1.00001 )) THEN
          write(6,*) ' sampled br outside of allowed range! ',ko,1./broi        
     *    ,br
        END IF
        goto 2930
      END IF
      costhe = 1 - temp
      IF (( ibcmp(irl) .EQ. 0 )) THEN
        Uj = 0
        goto 2960
      END IF
      br2 = br*br
      aux = ko*(ko-Uj)*temp
      pzmax = (aux - Uj)/sqrt(2*aux + Uj*Uj)
      IF((pzmax .LE. -1))goto 2920
      qc2 = 1 + br*br - 2*br*costhe
      qc = sqrt(qc2)
      IF (( pzmax .GT. 1 )) THEN
        pzmax = 1
        af = 0
        Fmax = 1
        fpz = 1
        goto 2970
      END IF
      aux3 = 1 + 2*Jo*abs(pzmax)
      aux4 = 0.5*(1-aux3*aux3)
      fpz = 0.5*exp(aux4)
      af = qc*(1+br*(br-costhe)/qc2)
      IF (( af .LT. 0 )) THEN
        IF((pzmax .GT. 0))fpz = 1 - fpz
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        eta_incoh = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        IF((eta_incoh .GT. fpz))goto 2920
        af = 0
        Fmax = 1
        goto 2970
      END IF
      IF (( pzmax .LT. -0.15 )) THEN
        Fmax = 1-af*0.15
        fpz1 = fpz*Fmax
      ELSE IF(( pzmax .LT. 0.15 )) THEN
        Fmax = 1 + af*pzmax
        aux3 = 1/(1+0.33267252734*aux3)
        aux4 = fpz*aux3*(0.3480242+aux3*(-0.0958798+aux3*0.7478556)) + e
     *  rfJo_array(j)
        IF (( pzmax .GT. 0 )) THEN
          fpz1 = 1 - Fmax*fpz - 0.31332853433*af/Jo_array(j)*aux4
          fpz = 1 - fpz
        ELSE
          fpz1 = Fmax*fpz - 0.31332853433*af/Jo_array(j)*aux4
        END IF
      ELSE
        Fmax = 1 + af*0.15
        fpz1 = 1 - Fmax*fpz
        fpz = 1 - fpz
      END IF
      IF (( rng_seed .GT. 24 )) THEN
        call ranlux(rng_array)
        rng_seed = 1
      END IF
      eta_incoh = rng_array(rng_seed)
      rng_seed = rng_seed + 1
      IF((eta_incoh .GT. fpz1))goto 2920
2970  CONTINUE
      IF (( rng_seed .GT. 24 )) THEN
        call ranlux(rng_array)
        rng_seed = 1
      END IF
      rnno18 = rng_array(rng_seed)
      rng_seed = rng_seed + 1
      rnno18 = rnno18*fpz
      IF (( rnno18 .LT. 0.5 )) THEN
        rnno18 = Max(1e-30,2*rnno18)
        pz_orig = 0.5*(1-Sqrt(1-2*Log(rnno18)))/Jo
      ELSE
        rnno18 = 2*(1-rnno18)
        pz_orig = 0.5*(Sqrt(1-2*Log(rnno18))-1)/Jo
      END IF
      IF((abs(pz_orig) .GT. 1))goto 2970
      IF (( pz_orig .LT. 0.15 )) THEN
        IF (( pz_orig .LT. -0.15 )) THEN
          frej = (1 - af*0.15)/Fmax
        ELSE
          frej = (1 + af*pz_orig)/Fmax
        END IF
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        eta = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        IF((eta .GT. frej))goto 2970
      END IF
      pz2 = pz_orig*pz_orig
      IF (( abs(pz_orig) .LT. 0.01 )) THEN
        br = br*(1 + pz_orig*(qc + (br2-costhe)*pz_orig))
      ELSE
        aux = 1 - pz2*br*costhe
        aux1 = 1 - pz2*br2
        aux2 = qc2 - br2*pz2*sinthe
        IF (( aux2 .GT. 1e-10 )) THEN
          br = br/aux1*(aux+pz_orig*Sqrt(aux2))
        END IF
      END IF
      Uj = Uj*prm
2960  pesg = br*peig
      pese = peig - pesg - Uj + prm
      sinthe = Sqrt(sinthe)
      call uphi(2,1)
      e(np) = pesg
      aux = 1 + br*br - 2*br*costhe
      IF (( aux .GT. 1e-8 )) THEN
        costhe = (1-br*costhe)/Sqrt(aux)
        sinthe = (1-costhe)*(1+costhe)
        IF (( sinthe .GT. 0 )) THEN
          sinthe = -Sqrt(sinthe)
        ELSE
          sinthe = 0
        END IF
      ELSE
        costhe = 0
        sinthe = -1
      END IF
      np = np + 1
      IF (( np .GT. 50 )) THEN
        write(6,'(a)') ' ***********************************************        
     *****'                                                                     
        write(6,'(a,a)') ' In subroutine ','COMPT',' stack size exceeded        
     *! '                                                                       
        write(6,'(a,i6,a,i6)') ' $MXSTACK = ',50,' np = ',np                    
        write(6,'(a)') ' Increase $MXSTACK and try again '                      
        write(6,'(a)') ' Terminating execution '                                
        write(6,'(a)') ' ***********************************************        
     *****'                                                                     
        stop
      END IF
      call uphi(3,2)
      e(np) = pese
      iq(np) = -1
      IF (( ibcmp(irl) .EQ. 1 )) THEN
        IF (( Uj .GT. 1e-3 )) THEN
          edep = 0
          call relax(Uj,shn_array(j),iz_array(j))
        ELSE
          edep = Uj
        END IF
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
            WRITE(6,2980)prob_RR
2980        FORMAT('**** Warning, attempt to play Roussian Roulette with        
     * prob_RR<=0! ',g14.6)                                                     
          END IF
        ELSE
          ip = NPold+1
2991      CONTINUE
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
            IF(((ip .GT. np)))GO TO2992
          GO TO 2991
2992      CONTINUE
          IF (( np .EQ. 0 )) THEN
            np = 1
            e(np) = 0
            iq(np) = 0
            wt(np) = 0
          END IF
        END IF
      END IF
      return
2920  return
      end
