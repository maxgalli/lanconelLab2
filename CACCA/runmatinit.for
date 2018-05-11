      subroutine runmatinit
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
c
      integer ispin
      integer i,j,k
c
c-------------------------------------------------------c
c
      read(ltydat,80,err=40,end=45) nmed,pcut(1),ecut(1),IRAYLR(1)
     * ,IEDGFL(1),IPHTER(1),ibcmp(1),iprdst,ibrdst,ibr_nist,ispin
      read(ltydat,81,err=40,end=45) estepe,ximax,transport_algorithm
     * ,bca_algorithm,skindepth_for_bca,i_do_rr(1),e_max_rr(1)
     * ,nbr_split,i_play_RR
c
      IF((iraylr(1).NE.1).AND.(iraylr(1).NE.0)) goto 40
      IF((iedgfl(1).NE.1).AND.(iedgfl(1).NE.0)) goto 40
      IF((iphter(1).NE.1).AND.(iphter(1).NE.0)) goto 40
      IF((ibcmp(1) .NE.1).AND.(ibcmp(1) .NE.0)) goto 40
      IF((iprdst   .LT.0).OR. (iprdst   .GT.2)) goto 40
      IF((ibrdst   .LT.0).OR. (ibrdst   .GT.1)) goto 40
      IF((ibr_nist .LT.0).OR. (ibr_nist .GT.1)) goto 40
      IF((ispin    .NE.1).AND.(ispin    .NE.0)) goto 40
c
      IF((estepe   .LE.0).OR. (estepe   .GE.1)) THEN
        estepe = 0.25
        WRITE(*,1540)estepe
1540    FORMAT(' using default value, estepe = : ',f10.3)
      ENDIF
c
      IF((ximax    .LE.0).OR. (ximax    .GE.1)) THEN
        ximax = 0.5
        WRITE(*,1570)ximax
1570    FORMAT(' using default, ximax = ',f10.3)
      END IF
c
      IF ((transport_algorithm.LT.0).OR.(transport_algorithm.GT.1))
     c goto 40
      IF ((bca_algorithm .LT. 0).OR.(bca_algorithm .GT. 1)) goto 40 
c
      IF (skindepth_for_bca .LE. 0) goto 40
c
c-------------------------------------------------------c
c
c  SKINDEPTH_FOR_BCA comment
c
c  In case of exact boundary crossing ( bca_algorithm.eq.0 ) 
c
c  This is the distance from a boundary (measured in elastic mean-free-paths)
c  at which the simulation switches to single scattering mode.
c  Best choice for efficiency is 3. If you set this parameter to a very large
c  number (e.g. 1e10), you can force single scattering simulation in the 
c  entire geometry (this is very slow).
c
c-------------------------------------------------------c
c 
c  In case of PRESTA algorithm ( bca_algorithm.eq.1 )
c
c  This is the distance from a boundary (measured in elastic mean-free-paths)
c  at which lateral deflections will be turned off. If you select a very 
c  large number (e.g. 1e10), standard EGS4 behaviour (no PRESTA) will result.
c  If you input a number < 1, this parameter will be determined in the way 
c  it was with PRESTA (i.e. depending on ECUT).                      
c
c-------------------------------------------------------c
c
      IF ((nbr_split .LE. 0)) THEN
         WRITE(6,1830)
1830  FORMAT(' Negative or zero value of nbr_split made 1=> no splitting
     *')
         nbr_split = 1
      END IF
c
      IF ((i_play_RR .NE. 0)) THEN
        i_play_RR = 1
        prob_RR = 1./float(nbr_split)
      END IF
c
      do 10 j=1,nmed
         read(ltydat,70,err=40,end=45) (media(i,j),i=1,24)
   10 continue
c
c----------- add to the  e-cut the mass energy ---------c
c
      ecut(1) = ecut(1) + rm
c
c same "pcut" and "ecut" for all regions ( and materials )
c
      do 15 j=2,maxreg
      	 pcut(j)     = pcut(1)
         ecut(j)     = ecut(1)
         IRAYLR(j)   = IRAYLR(1)
         IEDGFL(j)   = IEDGFL(1)
         IPHTER(j)   = IPHTER(1)
         ibcmp(j)    = ibcmp(1)
         i_do_rr(j)  = i_do_rr(1)
         e_max_rr(j) = e_max_rr(1)
   15 continue
c
      do 20 j=1,nmed
         IRAYLM(j) = IRAYLR(1)
   20 continue   
c
      IF (( ispin .EQ. 0 )) THEN
        spin_effects = .false.
      ELSE
        spin_effects = .true.
      END IF
c
      write (*,*) ' '
      do 30 j=1,nmed
        write(*,50) j,(media(i,j),i=1,24),pcut(j),ecut(j),IRAYLM(j)
     *,IEDGFL(j),IPHTER(j),ibcmp(j),iprdst,ibrdst,ibr_nist,ispin,estepe
     *,ximax,transport_algorithm,bca_algorithm,skindepth_for_bca
     *,i_do_rr(j),e_max_rr(j),nbr_split,i_play_RR
   30 continue  
c
      return
c
c-------------------------------------------------------c
c
   40 write(*,*) ' ERROR IN GEOMETRY FILE - RUNMATINIT '
   45 write(*,*) ' END IN GEOMETRY FILE - RUNMATINIT '
      close(ltydat)
c
      return
c
c-------------------------------------------------------c
c
   50 FORMAT(2X,'MED=',I3,2X,24A1,2X,'PCUT',E14.3,2X,'ECUT',E14.3,/,2X,
     * 'RAYLEIGH',I3,2X,'ATOMIC RELAXATION',I3,2X,
     * 'PHOTO-ELECTR.ANG.DISTR.',I3,2X,'BINDING',I3,/,2X,
     * 'PAIR ANG. DISTR.: 0=FIXED, 1=LEADING TERMS, 2=KOCH & MOTZ ', I6,
     * /,2X, 'BREMSST. ANG. DISTR.: 0=LEADING TERMS, 1=KOCH & MOTZ 2BS  
     *',I6,/,2X, 'BREMMST. DIFFER. PHOT. XSECT.: 0=BETHE-HEITLER, 1=NIST
     *-ICRU37',I3,/,2X,'SPIN EFFECTS',I3,4X,'MAXIMUM FRACTIONAL ENERGY L
     *OSS PER STEP',F14.7,/,2X,'MAXIMUM 1ST ELASTIC SCATTERING MOMENT PE
     *R STEP',F26.7,/,2X,'ELECTRON-STEP ALGORITHM: 0=DEFAULT, 1=PRESTA'
     *,I10,/,2X,'BOUNDARY CROSSING ALGORITHM: 0=EXACT, 1=PRESTA',I8,/
     *,2X,'SKIN-DEPTH FOR BCA ',F15.0,/
     *,2X,'ELECTRON RANGE REJECTION',I3,6X,'MAX ENERGY for REJECTION'
     *,F10.0,/,2X,'NUMBER OF BREM PHOTONS FOR EVENT',I20,/,2X
     *,'RUSSIAN ROULETTE ALL SECONDARY CHARGED PARTICLES',I4,/)
   70 FORMAT(24A1)
   80 FORMAT(I10,2F10.0,8I5)
   81 FORMAT(2F10.0,2I5,F15.0,I5,F10.0,2I5)
      END
c
