      SUBROUTINE BHABHA
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      DOUBLE PRECISION PEIP, PEKIN, PEKSE2, PESE1, PESE2, H1, DCOSTH
      real*4 EIP, EKIN, T0, E0, E02, YY, Y2, YP, YP2, BETA2, EP0
     *, EP0C, B1, B2, B3, B4, RNNO03, RNNO04, BR, REJF2, ESE1, ESE2
      NPold = NP
      PEIP=E(NP)
      EIP=PEIP
      PEKIN=PEIP-PRM
      EKIN=PEKIN
      T0=EKIN/RM
      E0=T0+1.
      YY=1./(T0+2.)
      E02=E0*E0
      BETA2=(E02-1.)/E02
      EP0=TE(MEDIUM)/EKIN
      EP0C=1.-EP0
      Y2=YY*YY
      YP=1.-2.*YY
      YP2=YP*YP
      B4=YP2*YP
      B3=B4+YP2
      B2=YP*(3.+Y2)
      B1=2.-Y2
2851  CONTINUE
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        RNNO03 = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        BR=EP0/(1.-EP0C*RNNO03)
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        RNNO04 = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        REJF2=(1.0-BETA2*BR*(B1-BR*(B2-BR*(B3-BR*B4))))
        IF((RNNO04.LE.REJF2))GO TO2852
      GO TO 2851
2852  CONTINUE
      IF (( np+1 .GT. 50 )) THEN
        write(6,'(a)') ' ***********************************************        
     *****'                                                                     
        write(6,'(a,a)') ' In subroutine ','BHABHA',' stack size exceede        
     *d! '                                                                      
        write(6,'(a,i6,a,i6)') ' $MXSTACK = ',50,' np = ',np+1                  
        write(6,'(a)') ' Increase $MXSTACK and try again '                      
        write(6,'(a)') ' Terminating execution '                                
        write(6,'(a)') ' ***********************************************        
     *****'                                                                     
        stop
      END IF
      IF ((BR.LT.0.5)) THEN
        IQ(NP+1)=-1
      ELSE
        IQ(NP)=-1
        IQ(NP+1)=1
        BR=1.-BR
      END IF
      BR=max(BR,0.0)
      PEKSE2=BR*EKIN
      PESE1=PEIP-PEKSE2
      PESE2=PEKSE2+PRM
      ESE1=PESE1
      ESE2=PESE2
      E(NP)=PESE1
      E(NP+1)=PESE2
      H1=(PEIP+PRM)/PEKIN
      DCOSTH=MIN(1.0D0,H1*(PESE1-PRM)/(PESE1+PRM))
      SINTHE=DSQRT(1.D0-DCOSTH)
      COSTHE=DSQRT(DCOSTH)
      CALL UPHI(2,1)
      NP=NP+1
      DCOSTH=H1*(PESE2-PRM)/(PESE2+PRM)
      SINTHE=-DSQRT(1.D0-DCOSTH)
      COSTHE=DSQRT(DCOSTH)
      CALL UPHI(3,2)
      RETURN
      END
