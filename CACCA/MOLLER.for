      SUBROUTINE MOLLER
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      DOUBLE PRECISION PEIE, PEKSE2, PESE1, PESE2, PEKIN, H1, DCOSTH
      real*4 EIE,  EKIN,  T0,  E0,  EXTRAE,  E02,  EP0,  G2,G3,  GMAX,
     *BR,  R,  REJF4,  RNNO27,  RNNO28,  ESE1,  ESE2
      NPold = NP
      PEIE=E(NP)
      EIE=PEIE
      PEKIN=PEIE-PRM
      EKIN=PEKIN
      T0=EKIN/RM
      E0=T0+1.0
      EXTRAE = EIE - THMOLL(MEDIUM)
      E02=E0*E0
      EP0=TE(MEDIUM)/EKIN
      G2=T0*T0/E02
      G3=(2.*T0+1.)/E02
      GMAX=(1.+1.25*G2)
3921  CONTINUE
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        RNNO27 = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        BR = TE(MEDIUM)/(EKIN-EXTRAE*RNNO27)
        R=BR/(1.-BR)
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        RNNO28 = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        REJF4=(1.+G2*BR*BR+R*(R-G3))
        RNNO28=GMAX*RNNO28
        IF((RNNO28.LE.REJF4))GO TO3922
      GO TO 3921
3922  CONTINUE
      PEKSE2=BR*EKIN
      PESE1=PEIE-PEKSE2
      PESE2=PEKSE2+PRM
      ESE1=PESE1
      ESE2=PESE2
      E(NP)=PESE1
      IF (( np+1 .GT. 50 )) THEN
        write(6,'(a)') ' ***********************************************        
     *****'                                                                     
        write(6,'(a,a)') ' In subroutine ','MOLLER',' stack size exceede        
     *d! '                                                                      
        write(6,'(a,i6,a,i6)') ' $MXSTACK = ',50,' np = ',np+1                  
        write(6,'(a)') ' Increase $MXSTACK and try again '                      
        write(6,'(a)') ' Terminating execution '                                
        write(6,'(a)') ' ***********************************************        
     *****'                                                                     
        stop
      END IF
      E(NP+1)=PESE2
      H1=(PEIE+PRM)/PEKIN
      DCOSTH=H1*(PESE1-PRM)/(PESE1+PRM)
      SINTHE=DSQRT(1.D0-DCOSTH)
      COSTHE=DSQRT(DCOSTH)
c mirko
c       write (*,*) ' ' 
c      write (22,*) 'moller1 ', u(np),v(np),w(np),np
c
      CALL UPHI(2,1)
c mirko
c      write (*,*) 'moller2 ', u(np),v(np),w(np),np
c
c 987  format (a8, 1x,3(1x,f15.12))
c
      NP=NP+1
      IQ(NP)=-1
      DCOSTH=H1*(PESE2-PRM)/(PESE2+PRM)
      SINTHE=-DSQRT(1.D0-DCOSTH)
      COSTHE=DSQRT(DCOSTH)
c mirko
c      write (*,*) 'moller3 ', u(np),v(np),w(np)
c
      CALL UPHI(3,2)
c mirko
c      write (*,*) 'moller4 ', u(np),v(np),w(np)
c
      RETURN
      END
