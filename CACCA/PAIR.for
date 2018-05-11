      SUBROUTINE PAIR
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      DOUBLE PRECISION PEIG,  PESE1,  PESE2
      real*4 EIG,  ESE2,  RNNO30, RNNO31, rnno32, rnno33, rnno34, DELTA
     *, REJF, rejmax, aux1, aux2, Amax, Bmax, del0, br, Eminus, Eplus
     *, Eavail, rnno_RR
      integer*4 L, L1
      real*4 ESE,  PSE, ZTARG, TTEIG, TTESE, TTPSE, ESEDEI, ESEDER
     *, XIMIN, XIMID, REJMIN, REJMID, REJTOP, YA, XITRY, GALPHA, GBETA
     *, XITST, REJTST_on_REJTOP, REJTST, RTEST
      integer*4 ICHRG
      NPold = NP
      IF (( i_play_RR .EQ. 1 )) THEN
        i_survived_RR = 0
        IF (( prob_RR .LE. 0 )) THEN
          IF (( n_RR_warning .LT. 50 )) THEN
            n_RR_warning = n_RR_warning + 1
            WRITE(6,4490)prob_RR
4490        FORMAT('**** Warning, attempt to play Russian Roulette with         
     *prob_RR<0! ',g14.6)                                                       
          END IF
        ELSE
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          rnno_RR = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF (( rnno_RR .GT. prob_RR )) THEN
            i_survived_RR =2
            IF (( np .GT. 1 )) THEN
              np = np-1
            ELSE
              wt(np) = 0
              e(np) = 0
            END IF
            return
          ELSE
            wt(np) = wt(np)/prob_RR
          END IF
        END IF
      END IF
      PEIG=E(NP)
      EIG=PEIG
      IF ((EIG.LE.2.1)) THEN
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        RNNO30 = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        ESE2 = PRM + 0.5*RNNO30*(PEIG-2*PRM)
      ELSE
        IF ((EIG.LT.50.)) THEN
          L = 5
          L1 = L + 1
          delta = 4*delcm(medium)/eig
          IF (( delta .LT. 1 )) THEN
            Amax = dl1(l,medium)+delta*(dl2(l,medium)+delta*dl3(l,medium
     *      ))
            Bmax = dl1(l1,medium)+delta*(dl2(l1,medium)+delta*dl3(l1,med
     *      ium))
          ELSE
            aux2 = log(delta+dl6(l,medium))
            Amax = dl4(l,medium)+dl5(l,medium)*aux2
            Bmax = dl4(l1,medium)+dl5(l1,medium)*aux2
          END IF
          aux1 = 1 - rmt2/eig
          aux1 = aux1*aux1
          aux1 = aux1*Amax/3
          aux1 = aux1/(Bmax+aux1)
        ELSE
          L = 7
          Amax = dl1(l,medium)
          Bmax = dl1(l+1,medium)
          aux1 = bpar(2,medium)*(1-bpar(1,medium)*rm/eig)
        END IF
        del0 = eig*delcm(medium)
        Eavail = eig - rmt2
4501    CONTINUE
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          RNNO30 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          RNNO31 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          RNNO34 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF (( rnno30 .GT. aux1 )) THEN
            br = 0.5*rnno31
            rejmax = Bmax
            l1 = l+1
          ELSE
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            rnno32 = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            rnno33 = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            br = 0.5*(1-max(rnno31,rnno32,rnno33))
            rejmax = Amax
            l1 = l
          END IF
          Eminus = br*Eavail + rm
          Eplus = eig - Eminus
          delta = del0/(Eminus*Eplus)
          IF (( delta .LT. 1 )) THEN
            rejf = dl1(l1,medium)+delta*(dl2(l1,medium)+delta*dl3(l1,med
     *      ium))
          ELSE
            rejf = dl4(l1,medium)+dl5(l1,medium)*log(delta+dl6(l1,medium
     *      ))
          END IF
          IF((( rnno34*rejmax .LE. rejf )))GO TO4502
        GO TO 4501
4502    CONTINUE
        ese2 = Eminus
      END IF
      PESE2=ESE2
      PESE1=PEIG-PESE2
      IF (( np+1 .GT. 50 )) THEN
        write(6,'(a)') ' ***********************************************        
     *****'                                                                     
        write(6,'(a,a)') ' In subroutine ','PAIR',' stack size exceeded!        
     * '                                                                        
        write(6,'(a,i6,a,i6)') ' $MXSTACK = ',50,' np = ',np+1                  
        write(6,'(a)') ' Increase $MXSTACK and try again '                      
        write(6,'(a)') ' Terminating execution '                                
        write(6,'(a)') ' ***********************************************        
     *****'                                                                     
        stop
      END IF
      E(NP)=PESE1
      E(NP+1)=PESE2
      IF (((IPRDST.EQ.1).OR.((IPRDST.EQ.2).AND.(EIG.LT.4.14)))) THEN
          DO 4511 ICHRG=1,2
          IF ((ICHRG.EQ.1)) THEN
            ESE=PESE1
          ELSE
            ESE=ESE2
          END IF
          PSE=SQRT(MAX(0.0,(ESE-RM)*(ESE+RM)))
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          COSTHE = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          COSTHE=1.0-2.0*COSTHE
          SINTHE=RM*SQRT((1.0-COSTHE)*(1.0+COSTHE))/(PSE*COSTHE+ESE)
          COSTHE=(ESE*COSTHE+PSE)/(PSE*COSTHE+ESE)
          IF ((ICHRG.EQ.1)) THEN
            CALL UPHI(2,1)
          ELSE
            NP=NP+1
            SINTHE=-SINTHE
            CALL UPHI(3,2)
          END IF
4511    CONTINUE
4512    CONTINUE
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        RNNO34 = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        IF ((RNNO34.LE.0.5)) THEN
          IQ(NP)=1
          IQ(NP-1)=-1
        ELSE
          IQ(NP)=-1
          IQ(NP-1)=1
        END IF
        RETURN
      ELSE IF(((IPRDST.EQ.2).AND.(EIG.GE.4.14))) THEN
        ZTARG=ZBRANG(MEDIUM)
        TTEIG=EIG/RM
          DO 4521 ICHRG=1,2
          IF ((ICHRG.EQ.1)) THEN
            ESE=PESE1
          ELSE
            ESE=ESE2
          END IF
          TTESE=ESE/RM
          TTPSE=SQRT((TTESE-1.0)*(TTESE+1.0))
          ESEDEI=TTESE/(TTEIG-TTESE)
          ESEDER=1.0/ESEDEI
          XIMIN=1.0/(1.0+(3.141593*TTESE)**2)
          REJMIN = 2.0+3.0*(ESEDEI+ESEDER) - 4.00*(ESEDEI+ESEDER+1.0-4.0
     *    *(XIMIN-0.5)**2)*( 1.0+0.25*LOG( ((1.0+ESEDER)*(1.0+ESEDEI)/(2
     *    .*TTEIG))**2+ZTARG*XIMIN**2 ) )
          YA=(2.0/TTEIG)**2
          XITRY=MAX(0.01,MAX(XIMIN,MIN(0.5,SQRT(YA/ZTARG))))
          GALPHA=1.0+0.25*LOG(YA+ZTARG*XITRY**2)
          GBETA=0.5*ZTARG*XITRY/(YA+ZTARG*XITRY**2)
          GALPHA=GALPHA-GBETA*(XITRY-0.5)
          XIMID=GALPHA/(3.0*GBETA)
          IF ((GALPHA.GE.0.0)) THEN
            XIMID=0.5-XIMID+SQRT(XIMID**2+0.25)
          ELSE
            XIMID=0.5-XIMID-SQRT(XIMID**2+0.25)
          END IF
          XIMID=MAX(0.01,MAX(XIMIN,MIN(0.5,XIMID)))
          REJMID = 2.0+3.0*(ESEDEI+ESEDER) - 4.00*(ESEDEI+ESEDER+1.0-4.0
     *    *(XIMID-0.5)**2)*( 1.0+0.25*LOG( ((1.0+ESEDER)*(1.0+ESEDEI)/(2
     *    .*TTEIG))**2+ZTARG*XIMID**2 ) )
          REJTOP=1.02*MAX(REJMIN,REJMID)
4531      CONTINUE
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            XITST = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            REJTST = 2.0+3.0*(ESEDEI+ESEDER) - 4.00*(ESEDEI+ESEDER+1.0-4
     *      .0*(XITST-0.5)**2)*( 1.0+0.25*LOG( ((1.0+ESEDER)*(1.0+ESEDEI
     *      )/(2.*TTEIG))**2+ZTARG*XITST**2 ) )
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            RTEST = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            THETA=SQRT(1.0/XITST-1.0)/TTESE
            REJTST_on_REJTOP = REJTST/REJTOP
            IF((((RTEST .LE. REJTST_on_REJTOP) .AND. (THETA .LT. PI) )))
     *      GO TO4532
          GO TO 4531
4532      CONTINUE
          SINTHE=SIN(THETA)
          COSTHE=COS(THETA)
          IF ((ICHRG.EQ.1)) THEN
            CALL UPHI(2,1)
          ELSE
            NP=NP+1
            SINTHE=-SINTHE
            CALL UPHI(3,2)
          END IF
4521    CONTINUE
4522    CONTINUE
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        RNNO34 = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        IF ((RNNO34.LE.0.5)) THEN
          IQ(NP)=1
          IQ(NP-1)=-1
        ELSE
          IQ(NP)=-1
          IQ(NP-1)=1
        END IF
        RETURN
      ELSE
        THETA=RM/EIG
      END IF
      CALL UPHI(1,1)
      NP=NP+1
      SINTHE=-SINTHE
      CALL UPHI(3,2)
      IF (( rng_seed .GT. 24 )) THEN
        call ranlux(rng_array)
        rng_seed = 1
      END IF
      RNNO34 = rng_array(rng_seed)
      rng_seed = rng_seed + 1
      IF ((RNNO34.LE.0.5)) THEN
        IQ(NP)=1
        IQ(NP-1)=-1
      ELSE
        IQ(NP)=-1
        IQ(NP-1)=1
      END IF
      RETURN
      END
