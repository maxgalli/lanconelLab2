Copyright Stanford  Mortran3.1   (FORTRAN 77 11JUN85)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
cc
cc.....................................................................c
cc                                                                     c
cc          Modified version for CT - SPECT - PET studies              c
cc                                                                     c
cc  Dante Bollini, Mirko Gombia, Nico Lanconelli, Alessandro Riccardi  c
cc                                                                     c
cc     Physics Department - University of Bologna                      c
cc     Lab. of Montecuccolino - DIENCA  - University of Bologna        c
cc                                                                     c
cc                          June 2001                                  c
cc.....................................................................c
cc
      program main
c
      implicit none
c

      include 'egs4comm.for'
      include 'egs4fcomm.for'

      include 'BLOCK_DATA.for'

c----------------------------------------------------------------------c
c----- timing of execution --------------------------------------------c
c
c  step 1.  pre-hatchscan histin.for-call-initialization
c
      call timest(totim)
      call nrc_init
      call histin
c
c----------------------------------------------------------------------c
c
c  step 2.  hatch-call
c
c     initialize the maximum number of regions
      numgeom=120001
      write (*,*)  ' ENTERING HATCH SUBROUTINE ...'
      write (*,*) ' totim = ' ,totim,' tevent = ',
     a              tevent,' tlim =  ' , tlim
c
c      write (*,*)  ' prima di HATCH SUBROUTINE ...'
      call hatch
c      write (*,*)  ' dopo di HATCH SUBROUTINE ...'
c
      WRITE (*,1020) AE(1)-0.511, AP(1)
1020  FORMAT(/' knock-on electrons can be created and any electron follo
     *wed down to' /T40,F8.3,' MeV kinetic energy'/ ' brem photons can b
     *e created and any photon followed down to      ', /T40,F8.3,' MeV
     *')
c
      close ( 08 )
      close ( 12 )
c
c      kmpi      =  material   for hatch  unit = 12
c      kmpo      =  echo unit for hatch   unit =  8
c
c----------------------------------------------------------------------c
c
c  step 3.  shower-call
c
c 
      WRITE(*,*) ' STARTING SHOWER SIMULATION ...'
      tlim = 0.
c
      call timed(tevent)
      call timex(totim)
      write(*,*) ' totim = ' ,totim,' tevent = ',
     a                                         tevent,' tlim =  ' , tlim
c
c-----  loop over sources
c
      do 30 nactsx = 1, isxnum
c
         ncases = isxev (nactsx )
         write(*,*) ' ISX ', nactsx ,' , EVENTI ' , ncases
         if (ncases.le.1) ncases=1
c
c -----  loop over events
c
         nevent = 0
         iloop  = 0
         ilsx   = 0
c
   10    continue
c
c--------   new event  --  increment of event in source  ----c
c
         tlim = 0.0
c
c-----   timing every itcxt events   ------------------------c
c
         if  (  mod(nevent,itctx).eq.0 .or. nevent.eq.1 ) then
cnico   TIMED  returns the execution time interval since the last call to TIMED
             call timed ( tevent )
             tlim = tlim + tevent
             tmed = tlim / float (itctx + 1)
cnico    TIMEX returns the CPU  seconds used so far (execution time)
cnico    The subroutine TIMEST may be needed to initialise the clock
             call timex ( totim )
             write (*,60) nevent,totim,tlim,tmed,ilsx,iloop
         endif
c
         call evinit
c
         call source
c
         eireal = ein
         call shower(iqin,eireal,xin, yin, zin, uin, vin, win,irin,wtin)
c
         call evsum
c
c
c---- end of event loop
c
         if( nevent - ncases ) 10,20,20
   20        write(*,*) ' END OF SX ',nactsx,'   NEVENT:',
     a                    nevent,'(',iloop,ilsx,')'
c
c---- end of sources loop
c
   30 continue
c
c-----------------------------------------------------------c
c
      continue
c
      call histout
      tlim = 0.
      call timed(tevent)
      call timex(totim)
      write(*,*) ' '
      write(*,*) ' ************** END OF RUN ', totim ,tevent
      write(*,*) '-----------------------------------------------------'
      write(*,70) totein, toteoob, 100*toteoob/totein , totedeodx
     a, 100*totedeodx/totein, totemod, 100*totemod/totein, totebio
     b, 100*totebio/totein, totedet, 100*totedet/totein, toteew
     c, 100*toteew/totein, 100*(toteoob+totedeodx+totemod)/totein
     d, 100*(toteoob+totebio+totedet+toteew)/totein
      write(*,*) '-----------------------------------------------------'
c
   60 FORMAT(1X,' EVENT RUN ',1X,I10,2(2X,F10.3),2X,E12.3,1X,I6,1X,I10)
   70 FORMAT(/,1X, ' INPUT TOTAL ENERGY      =', F14.4,1X,'MeV',//
     a,1X, ' ENERGY OUT OF BOUNDS    =',F14.4,1X,'MeV',4X,F8.4,1X,'%',/
     b,1X, ' ENERGY IN NON-MODULAR   =',F14.4,1X,'MeV',4X,F8.4,1X,'%',/
     b,1X, ' ENERGY IN MODULAR       =',F14.4,1X,'MeV',4X,F8.4,1X,'%',//
     c,1X, ' ENERGY IN BIOMASS       =',F14.4,1X,'MeV',4X,F8.4,1X,'%',/
     d,1X, ' ENERGY IN DETECTOR      =',F14.4,1X,'MeV',4X,F8.4,1X,'%',/
     e,1X, ' ENERGY ELSEWHERE        =',F14.4,1X,'MeV',4X,F8.4,1X,'%',//
     f        ,1X, ' TOTAL 1 should be 100 % = ',F8.4, 1X,'%',/
     g        ,1X, ' TOTAL 2 should be 100 % = ',F8.4, 1X,'%',/)
      END
c
c-----------------------------------------------------------c
c
      BLOCK DATA
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      CHARACTER*4 MEDIA1(24)
      EQUIVALENCE(MEDIA1(1),MEDIA(1,1))
c
      include 'BLOCK_DATA2.for'
c
      DATA IBRDST/1/
      data ibr_nist/0/
      DATA IPRDST/1/
      DATA binding_energies/ 600*0.0/
      DATA IAUSFL/5*1,23*0/,RHOF/1.0/
      data ximax /0.5/,  estepe/0.25/, skindepth_for_bca/3/
      data transport_algorithm/0/, bca_algorithm/0/, exact_bca/.true./,
     *spin_effects/.true./
      data count_pII_steps/0./, count_all_steps/0./
      DATA NMED /1/, MEDIA1/'N','A','I',' ',' ',' ',' ',' ',' ',' ',' ',
     *' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '/
      DATA IRAYLM/21*0/
      DATA KMPI/12/,KMPO/8/,DUNIT/1./
      DATA rng_seed /100/
      DATA LATCHI /0/
      DATA RMT2/1.0220068/,RMSQ/.26112447/
      DATA PI/3.141593/,TWOPI/6.283185/,PI5D2/7.853982/
      DATA RM/.5110034/
      data nbr_split/1/
      data i_play_RR/0/
      data i_survived_RR/0/
      data prob_RR/-1.0/
      data n_RR_warning/0/
      END
c
      subroutine nrc_init
c
c-------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
c
      character*80 geomfilename
cale

      integer*4 lnblnk1
cale
c-----  open run description file  ---------------------c
c
      write(*,*) 'Enter file with geometry'
c
      read(*,90) geomfilename
      write(*,100) geomfilename
c
      open(unit=ltydat,file=geomfilename,form='formatted',status='old')
cale
      read(ltydat,90) pathname
      pathlength = lnblnk1(pathname)

cale      write(*,*)'cale:', pathlength
cale
      read(ltydat,90) histofilename

c
c-------------------------------------------------------c
c
      write(*,*)' ltydat : ENTERING EGS4INIT/randomnminit'
      call randomnminit
c
c-------------------------------------------------------c
c---  read running data  -------------------------------c
c
      write(*,*)' ltydat : ENTERING runparinit'
      call runparinit 
c
c-------------------------------------------------------c
c---  media  to be processed in hatch  -----------------c
c
      write(*,*)' ltydat : ENTERING runmatinit'
      call runmatinit
c
c-------------------------------------------------------c
c
      write(*,*)' ltydat : ENTERING rungeominit'
      call rungeominit
c
c-------------------------------------------------------c
c
      write(*,*)' ltydat : ENTERING runsxinit'
      call runsxinit
c
c-------------------------------------------------------c
c
      iqin=ichar
      irin=0
      wtin=1.0
c
c-------------------------------------------------------c
c
      write(*,*)' ltydat : ENTERING EGS4INIT/close(ltydat)'
      close(ltydat)
      write(*,*) ' END OF  EGS4 INIT '
      write(*,*) '-------------------------------------'
      return
c
c-------------------------------------------------------c
c
   10 write(*,*) ' ERROR IN  GEOMETRY DATA FILE - INIT '
      close(ltydat)
      return
c
   90 FORMAT(A80)
  100 FORMAT(/,'  GEOMETRY FILE IS ', a20,/)
      END
c


      subroutine randomnminit
c
c-------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
c
      integer idm1, idm2, iseed
      real float_temp
      integer*4  i
c
      read (ltydat,7800,err=700,end=700) iseed
c
c----- INIT TO TIME AND DATE ( HOPEFULLY RANDOM )
c      START A NEW RANDOM SEQUENCE

      call datime(idm1,idm2)
      float_temp = float (iseed) / float(isl(1)*isl(6)+isl(2)*isl(3)
     a*isl(4)*isl(5))
      write (*,6700) isl(1),isl(2),isl(3),isl(4),isl(5),isl(6),iseed,
     afloat_temp
      iseed = int( ( 10 ** 9 ) * (float_temp - int(float_temp)))
      write (*,6750) iseed
      call init_ranlux(4,iseed)
      call ranlux(rng_array)
      write (*,*) ' RNG_array = '
      write (*,6770) (rng_array(i),i=1,24)
      rng_seed = 1
c
      return
c
c-------------------------------------------------------c
c
  700 write(*,*) ' ERROR IN RANDOM SEED  '
      return
c
 6700 FORMAT (/,1X, ' ISL1   = ', i10 , 5X , 'ISL2       = ', i10 , / ,
     c          1X, ' ISL3   = ', i10 , 5X , 'ISL4       = ', i10 , / ,
     c          1X, ' ISL5   = ', i10 , 5X , 'ISL6       = ', i10 , / ,
     c          1X, ' ISEED  = ', i10 , 5X , 'FLOAT_TEMP = ', f20.10)
 6750 FORMAT (1X, ' ISEED2 = ', I10 )
 6770 FORMAT ((4(2X,F12.10)),/,(4(2X,F12.10)),/,(4(2X,F12.10)),/,
     c        (4(2X,F12.10)),/,(4(2X,F12.10)),/,(4(2X,F12.10)),/)
 7800 FORMAT(I10)
      END
c
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      subroutine ranlux(rng_array)
      implicit none
      real*4 rng_array(24)
      integer*4 seedin,luxury_level
      integer*4 state(25)
      integer*4 ounit
      character*(*) fmt_flags
      integer*4 seeds(24),carry
      integer*4 i24,j24
      integer*4 next(24)
      integer*4 jseed_dflt,nskip,icon,j,k,status,jseed,nskipll(0:4),icar
     *ry
      logical not_initialized
      real*4 twom24,twop24
      integer*4 uni
      save seeds,carry,i24,j24,next,twom24,not_initialized, nskip,twop24
     *,nskipll
      data nskipll/0,24,73,199,365/
      data jseed_dflt/314159265/, icon/2147483563/
      data not_initialized/.true./
      IF (( not_initialized )) THEN
        not_initialized = .false.
        nskip = nskipll(1)
        twom24 = 1
        twop24 = 1
        jseed = jseed_dflt
          DO 1981 j=1,24
          twom24 = twom24 * 0.5
          twop24 = twop24 * 2
          k = jseed/53668
          jseed = 40014*(jseed-k*53668)-k*12211
          IF (( jseed .LT. 0 )) THEN
            jseed = jseed + icon
          END IF
          seeds(j) = mod(jseed,16777216)
          next(j) = j-1
1981    CONTINUE
1982    CONTINUE
        next(1) = 24
        i24 = 24
        j24 = 10
        carry = 0
        IF (( seeds(24) .EQ. 0 )) THEN
          carry = 1
        END IF
      END IF
        DO 1991 j=1,24
        uni = seeds(j24) - seeds(i24) - carry
        IF (( uni .LT. 0 )) THEN
          uni = uni + 16777216
          carry = 1
        ELSE
          carry = 0
        END IF
        seeds(i24) = uni
        i24 = next(i24)
        j24 = next(j24)
        IF (( uni .GE. 4096 )) THEN
          rng_array(j) = uni*twom24
        ELSE
          rng_array(j) = uni*twom24 + seeds(j24)*twom24*twom24
        END IF
1991  CONTINUE
1992  CONTINUE
      IF (( nskip .GT. 0 )) THEN
          DO 2001 j=1,nskip
          uni = seeds(j24) - seeds(i24) - carry
          IF (( uni .LT. 0 )) THEN
            uni = uni + 16777216
            carry = 1
          ELSE
            carry = 0
          END IF
          seeds(i24) = uni
          i24 = next(i24)
          j24 = next(j24)
2001    CONTINUE
2002    CONTINUE
      END IF
      return
      entry init_ranlux(luxury_level,seedin)
      jseed = seedin
      IF((jseed .LE. 0))jseed = jseed_dflt
      IF (( luxury_level .LT. 0 .OR. luxury_level .GT. 4 )) THEN
        luxury_level = 1
      END IF
      nskip = nskipll(luxury_level)
      WRITE(6,2010)luxury_level,jseed
2010  FORMAT(//' ***************** RANLUX initialization ***************        
     ****'/, ' luxury level: ',i2,/, ' initial seed: ',i12,/, '*********        
     ***************************************************'//)                    
      not_initialized = .false.
      twom24 = 1
      twop24 = 1
        DO 2021 j=1,24
        twom24 = twom24 * 0.5
        twop24 = twop24 * 2
        k = jseed/53668
        jseed = 40014*(jseed-k*53668)-k*12211
        IF (( jseed .LT. 0 )) THEN
          jseed = jseed + icon
        END IF
        seeds(j) = mod(jseed,16777216)
        next(j) = j-1
2021  CONTINUE
2022  CONTINUE
      next(1) = 24
      i24 = 24
      j24 = 10
      carry = 0.
      IF (( seeds(24) .EQ. 0 )) THEN
        carry = 1
      END IF
      return
      entry get_ranlux_state(state)
        DO 2031 j=1,24
        state(j) = seeds(j)
2031  CONTINUE
2032  CONTINUE
      state(25) = i24 + 100*(j24 + 100*nskip)
      IF((carry .GT. 0))state(25) = -state(25)
      return
      entry set_ranlux_state(state)
      twom24 = 1
      twop24 = 1
        DO 2041 j=1,24
        twom24 = twom24 * 0.5
        twop24 = twop24 * 2
        next(j) = j-1
2041  CONTINUE
2042  CONTINUE
      next(1) = 24
        DO 2051 j=1,24
        seeds(j) = state(j)
2051  CONTINUE
2052  CONTINUE
      IF (( state(25) .LE. 0 )) THEN
        status = -state(25)
        carry = 1
      ELSE
        status = state(25)
        carry = 0
      END IF
      nskip = status/10000
      status = status - nskip*10000
      j24 = status/100
      i24 = status - 100*j24
      IF (( j24 .LT. 1 .OR. j24 .GT. 24 .OR. i24 .LT. 1 .OR. i24 .GT. 24
     * )) THEN
        WRITE(6,2060)state(25),nskip,i24,j24
2060    FORMAT('// *********** Error in set_ranlux_state: seeds outside         
     *of allowed range!'/, '   status = ',i8/, '   nskip  = ',i8/, '   i        
     *24    = ',i8/, '   j24    = ',i8/, '******************************        
     ******************************************'//)                             
        stop
      END IF
      not_initialized = .false.
      return
      entry show_ranlux_seeds(ounit)
      IF (( carry .GT. 0 )) THEN
        icarry = 1
      ELSE
        icarry = 0
      END IF
      write(ounit,'(a,i4,a,2i3,a,i2,$)') ' skip = ',nskip,' ix jx = ',i2        
     *4,j24,' carry = ',icarry                                                  
      return
      entry print_ranlux_seeds(ounit,fmt_flags)
      IF (( carry .GT. 0 )) THEN
        icarry = 1
      ELSE
        icarry = 0
      END IF
      write(ounit,fmt_flags) nskip,i24,j24,icarry
      return
      end


       subroutine runparinit
c
c+------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
c
c
c-------------------------------------------------------c
c
      read (ltydat,60,err=10,end=15) ichar,labeldet,labelbiomass,enkin
      read (ltydat,60,err=10,end=15) ncases,itctx
c
      evkin = enkin + abs(float(ichar)) * rm
c
      write (*,20) enkin,evkin,ichar
      write (*,30) labeldet,labelbiomass
c
c-------------------------------------------------------c
c
c SPREN = Absolute energy resolution (in MeV ):
c a characteristic value for Germanium is 0.00112  MeV (@ 0.1405 MeV)
c a characteristic value for  NaI(Tl) is 0.0091325 MeV (@ 0.1405 MeV)
c
      read(ltydat,70,err=10,end=15) maxhit , trigger, treshold, esnr
     a, spren, ispres
c
      write(*,50) maxhit, trigger, treshold, esnr, spren, ispres
c
      return
c
c-------------------------------------------------------c
c
   10 write (*,*) ' ERROR IN GEOMETRY DATA FILE - RUNPARINIT '
   15 write (*,*) ' END IN GEOMETRY DATA FILE - RUNPARINIT '
      close (ltydat)
c
      return
c
c-------------------------------------------------------c
c
   20 FORMAT(/,2X,'ENKIN     ',E15.4,1X,'MeV',
     +       /,2X,'ETOT      ',E15.4,1X,'MeV',
     +       /,2X,'CHARGE    ',I23)
   30 FORMAT(  2X,'KOUNTLABEL',I23,
     +       /,2X,'BIOMLABEL ',I23)
   50  FORMAT('  MAXIMUM NUMBER OF HITS       ', I5  ,/,
     +        '  TRIGGER                      ', F10. 4,1X,'MeV',/,
     +        '  THRESHOLD ENERGY DETECTOR    ', F10. 4,1X,'MeV',/,
     +        '  SIGNAL/NOISE VALUE           ', F10. 4,/,
     +        '  ENERGY RESOLUTION (absolute) ', F10. 4,1X,'MeV',/,
     +        '  SPATIAL RESOLUTION           ', I5 ,/)
   60  FORMAT(3I10,E20.5)
   70  FORMAT(I10,4F10.0,I5)
c
       END
c
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
      subroutine rungeominit
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
c
c
      integer  k, j
      integer itmpidm,numslab,cont_cell,mod_tempx,mod_tempy
      real tempinf, tempsup
      doubleprecision xdet1,ydet1,xdet2,ydet2
      integer numcelly_det1mod,numcelly_det2mod
c
      read(ltydat,80,err=41,end=46) nplan
c
      if(nplan.ge.maxreg) then
         write(*,*) ' error in region number ' , nplan,' > ',maxreg
      endif
c
      do 10 k=1,nplan
         read(ltydat,100,err=41,end=46) indmat(k), zinit(k),zend(k), 
     + xinf(k),xsup(k),yinf(k),ysup(k),radius(k),index(k),iregtype(k)
c
      write(*,100) indmat(k),
     + zinit(k),zend(k),xinf(k),xsup(k),yinf(k),ysup(k),radius(k)
     + ,index(k),iregtype(k)
c
c ale iregtype,radius 09 04 2001 
c
   10 continue
c
      nreg = nplan
      indmat(maxreg+1) = 0
      xmax  = -10000.0
      xmin  = +10000.0
      ymax  = -10000.0
      ymin  = +10000.0
      zmin  = +10000.0
      zmax  = -10000.0
c
      numslab = 0
c
c------------------------------------------------------c
c
      write (*,*) ' '
      write (*,*) ' DETECTORS '
      write (*,*) ' '
c
      do 30 k=1,nplan
c
c---  sort the geometrical limits of the detector
c
      tempinf  = amin1(zinit(k),zend(k))
      tempsup  = amax1(zinit(k),zend(k))
      zinit(k) = tempinf
      zend(k)  = tempsup
      tempinf  = amin1(xinf(k),xsup(k))
      tempsup  = amax1(xinf(k),xsup(k))
      xinf(k)  = tempinf
      xsup(k)  = tempsup
      tempinf  = amin1(yinf(k),ysup(k))
      tempsup  = amax1(yinf(k),ysup(k))
      yinf(k)  = tempinf
      ysup(k)  = tempsup
c
c-------
c
      if ( (index(k)/10000).eq.labeldet) then 
c
         write(*,100) indmat(k), zinit(k), zend(k), xinf(k), xsup(k)
     +, yinf(k), ysup(k), index(k), iregtype(k)
c
         zdetend = zend(k)
         ydetinf = yinf(k)
         xdetinf = xinf(k)
         ydetsup = ysup(k)
         xdetsup = xsup(k)
c
      endif
c
c ale iregtype, raysqrt, radius 11 04 2001
c
      raysqrt(k)= radius(k)**2
   30 continue
c
      write(*,60) nplan
c
c-------------------------------------------------------c
c
      do 20 k=1,nplan
c
c---  set geometrical limits of apparatus
c---
c
      if(zinit(k).le.zmin)     zmin  = zinit(k)
      if(zend(k).ge.zmax)      zmax  = zend (k)
c
      if(xinf(k).le.xmin)      xmin  = xinf(k)
      if(xsup(k).ge.xmax)      xmax  = xsup(k)
      if(yinf(k).le.ymin)      ymin  = yinf(k)
      if(ysup(k).ge.ymax)      ymax  = ysup(k)
c
      med(k)=indmat(k)
      itmpidm = index(k)/10000
c
      write(*,70) k,indmat(k),(media(j,indmat(k)),j=1,24), med(k),
     +pcut(indmat(k)),ecut(indmat(k)), zinit(k),zend(k),xinf(k),xsup(k),
     +yinf(k),ysup(k), index(k),itmpidm
 
c
   20 continue
c
c----
c
      tracklim  = sqrt ((zmax-zmin)*(zmax-zmin)+
     a                   (xmax-xmin)*(xmax-xmin)+
     b                   (ymax-ymin)*(ymax-ymin))
      steplim = 1.e-04
 
      write(*,*) ' GEOM_LIMIT ',zmin,zmax,xmin,xmax,ymin,ymax
      write(*,*) ' TRACKLIMIT ',tracklim
      write(*,*) ' STEPLIM    ',steplim
      write(*,*) ' '
cc
cale
ccc   NICO aggiunge x geometria modulare detector
ccc 
       read(ltydat,95,err=41,end=46) modlabel,zinit_det1mod,
     a zend_det1mod,xinf_det1mod,xsup_det1mod,yinf_det1mod,ysup_det1mod
     b ,xcell_det1mod,ycell_det1mod,xair_det1mod,yair_det1mod,
     c indlead_det1mod,indair_det1mod
c
      xdet1 = xsup_det1mod - xinf_det1mod
      ydet1 = ysup_det1mod - yinf_det1mod
      xcellinv_det1mod = 1./xcell_det1mod
      ycellinv_det1mod = 1./ycell_det1mod
      numcellx_det1mod = dint (xdet1 / xcell_det1mod)
      numcelly_det1mod = dint (ydet1 / ycell_det1mod)
      numcell_det1mod = numcellx_det1mod * numcelly_det1mod
      regoffset_det1mod = nreg
cnico cambia x dead zone      nreg = nreg + numcell_det1mod
      nreg = nreg + numcell_det1mod * 3
      write(*,*) 'xdet1mod ',xcell_det1mod,xdet1/xcell_det1mod
      write(*,*) 'xdet1 ydet1 ncellxdet1 ncellydet1 regoffsdet1 nreg'
     a,xdet1,ydet1,numcellx_det1mod,numcelly_det1mod,regoffset_det1mod
     b,nreg
      cont_cell = 0
      do 11 k = regoffset_det1mod+1,regoffset_det1mod+numcell_det1mod
     	    med(k) = indair_det1mod
	    indmat(k) = indair_det1mod
	    index(k) = labeldet * 10000
cnico aggiunge x dead zone
     	    med(k+numcell_det1mod) = indlead_det1mod
	    indmat(k+numcell_det1mod) = indlead_det1mod
	    index(k+numcell_det1mod) = 0
     	    med(k+numcell_det1mod+numcell_det1mod) = indlead_det1mod
	    indmat(k+numcell_det1mod+numcell_det1mod) = indlead_det1mod
	    index(k+numcell_det1mod+numcell_det1mod) = 0
cnico
	    mod_tempx = mod(cont_cell,numcellx_det1mod)
	    mod_tempy = int(cont_cell/numcellx_det1mod)
	    xinf(k) = xinf_det1mod + mod_tempx * xcell_det1mod
	    xsup(k) = xinf(k) + xcell_det1mod
	    yinf(k) = yinf_det1mod + mod_tempy * ycell_det1mod
	    ysup(k) = yinf(k) + ycell_det1mod
	    zinit(k) = zinit_det1mod
	    zend(k) = zend_det1mod
	    cont_cell = cont_cell + 1
            if(zend(k).gt.zdetend) zdetend = zend(k)
            if(yinf(k).lt.ydetinf) ydetinf = yinf(k)
            if(xinf(k).lt.xdetinf) xdetinf = xinf(k)
            if(ysup(k).gt.ydetsup) ydetsup = ysup(k)
            if(xsup(k).gt.xdetsup) xdetsup = xsup(k)
   11 continue	    
ccc
ccc   NICO aggiunge x geometria modulare collimatore
ccc 
       read(ltydat,95,err=41,end=46) modlabel,zinitcoll,zendcoll,
     a xinfcoll,xsupcoll,yinfcoll,ysupcoll,xcell,ycell,xair,yair,
     b indlead,indair
c
       if(modlabel.eq.1) then 
c
c     moduli esagonali
c       
c
       septa = xcell
       xcoll = xsupcoll - xinfcoll
       ycoll = ysupcoll - yinfcoll
c       xair = sqrt(3.) * xair / 2
       xcell = 2. * (xair + septa / 2.)
       xcell_2 = xcell / 2.
       xlead = xcell_2 - xair
       xcellinv = 1./xcell_2
       ycell = xcell * 2. / sqrt(3.)
       numcellx = dint (xcoll / xcell_2) +1
       numcelly = dint (ycoll / (ycell/4) + 2)
       numcell = numcellx * numcelly
       regoffset = nreg 
       u3offset = (int((numcellx + 1)/2) ) * xcell_2
       u3celloffset = int((numcellx + 1)/2) + 1
       nreg = nreg + 2 * numcell
       write(*,*) 'zinitcoll zendcoll, xinfcoll xsupcoll, yinfcoll
     a ysupcoll xcell ycell xair yair indlead indair',zinitcoll,zendcoll
     b ,xinfcoll,xsupcoll,yinfcoll,ysupcoll,xcell,ycell,xair,yair,
     c indlead,indair
       write(*,*) 'xcoll ycoll numcellex numcelley regoffset nreg'
     a ,xcoll,ycoll,numcellx,numcelly,regoffset,nreg
       write(*,*) 'u3off u3celloff xcell_2', u3offset,
     a u3celloffset,xcell_2 
	if(mod(u3celloffset,2).eq.0) then
		label_hex = 0
	else
		label_hex = 1
	endif
	label_hex1 = mod(int((2*u3celloffset-1)/2),3)
       write(*,*) 'label_0 label_1',label_hex,label_hex1
       do 12 k = 1,numcellx
          do 12 j = 1,numcelly
		iu1 = k
		iu2 = int((j+k)/2)
		iu3 = int((j-k+2*u3celloffset-1)/2)
		if(label_hex.eq.0) then
		icell_temp = regoffset + (iu1-1)*numcelly + 
     a		2*(iu2-int((iu1+1)/2)) + mod(abs(iu3-iu2+1),2) + 
     b		2*mod(iu1,2)*mod(abs(iu3-iu2),2)
     		else
		icell_temp = regoffset + (iu1-1)*numcelly + 
     a		2*(iu2-int((iu1+1)/2)) + mod(abs(iu3-iu2),2) + 
     b		2*mod(iu1,2)*mod(abs(iu3-iu2+1),2)
		endif
		med(icell_temp) = indair
		indmat(icell_temp) = indair
		med(icell_temp+numcell) = indlead
		indmat(icell_temp+numcell) = indlead
		index(icell_temp) = 0
		index(icell_temp+numcell) = 0
   12  continue
c
c    fine moduli esagonali       
c
       elseif(modlabel.eq.2) then 
c
c     moduli rettangolari
c       
      xcoll = xsupcoll - xinfcoll
      ycoll = ysupcoll - yinfcoll
      xcellinv = 1./xcell
      ycellinv = 1./ycell
      numcellx = dint (xcoll / xcell)
      numcelly = dint (ycoll / ycell)
      numcell = numcellx * numcelly
      regoffset = nreg
      nreg = nreg + 3 * numcell
      do 13 k = regoffset+1,regoffset+numcell
     	    med(k) = indair
	    indmat(k) = indair
	    med(k+numcell) = indlead
	    indmat(k+numcell) = indlead
	    med(k+numcell+numcell) = indlead
	    indmat(k+numcell+numcell) = indlead
   13 continue	    
       write(*,*) 'zinitcoll zendcoll, xinfcoll xsupcoll, yinfcoll
     a ysupcoll xcell ycell xair yair indlead indair',zinitcoll,zendcoll
     b ,xinfcoll,xsupcoll,yinfcoll,ysupcoll,xcell,ycell,xair,yair,
     c indlead,indair
       write(*,*) 'xcoll ycoll numcellex numcelley regoffset nreg'
     a ,xcoll,ycoll,numcellx,numcelly,regoffset,nreg
       write(*,*) 'u3off u3celloff xcell_2', u3offset,
     a u3celloffset,xcell_2 
c
c    fine moduli rettangolari       
c
       elseif(modlabel.eq.0) then 
         regoffset = nreg
       write(*,*) 'zinitcoll zendcoll, xinfcoll xsupcoll, yinfcoll
     a ysupcoll xcell ycell xair yair indlead indair',zinitcoll,zendcoll
     b ,xinfcoll,xsupcoll,yinfcoll,ysupcoll,xcell,ycell,xair,yair,
     c indlead,indair
       write(*,*) 'xcoll ycoll numcellex numcelley regoffset nreg'
     a ,xcoll,ycoll,numcellx,numcelly,regoffset,nreg
       write(*,*) 'u3off u3celloff xcell_2', u3offset,
     a u3celloffset,xcell_2 
       endif
ccc
ccc  fine aggiunta x geometria modulare collimatore
ccc	
cale
c------------------------------------------------------c
c
      return
c
c------------------------------------------------------c
c
   41 write(*,*) ' ERROR IN GEOMETRY FILE - RUNGEOMINIT '
   46 write(*,*) ' END IN GEOMETRY FILE - RUNGEOMINIT '
      close(ltydat)
c
      return
c
c-----------------------------------------------------c
c
   50 FORMAT(1X,I9,1X,6F10.3)
   60 FORMAT(/,2X,' NUMERO PIANI ',I10,/)
   70 FORMAT( I5,' INDMED =',I3,', MED =',2X,24A1,2X,I3,/,5X
     +,' CUTS _ P & E ',2X,E16.4,2X,E16.4,/, 5X ,' COUNT SIZE  '
     +, 6(1X,F 10.3),/,5X,' COUNT INDEX  ', I10 , ' INDX/10000 ', I10,/)
   80 FORMAT(I10)
   90 FORMAT(I10,6F10.0,I10)
   95 FORMAT(I2,10F10.4,2I2)
c
c ale iregtype radius format 09 04 2001
c
  100 FORMAT(I10,7F10.4,I10,I10)
c
c ale format iregtype radius09 04 2001
c
      END
c
      subroutine runsxinit
c
c--------------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
      integer k
c
c--------------------------------------------------------------c
c
      write (*,*) ' '
      read(ltydat,50,err=20,end=20) isxnum
c
      if (isxnum.gt.20) then
        write(*,*)' ISXNUM is TOO large', isxnum
        isxnum = 20
        write(*,*)' ISXNUM HAS LIMITED TO', isxnum
c
      else
        write(*,*)' ISXNUM = ', isxnum
      endif
c
c-------
c
      do 10 k=1,isxnum
c
      read(ltydat,60,err=20,end=25) isxtype(k),isxev(k), sx_cent(k),
     + sy_cent(k),sz_cent(k),rsx(k),rsy(k),rsz(k),esource(k),
     + thetaangsx(k), phiangsx(k)
c
      rsxsq(k) = rsx(k)*rsx(k)
      rsysq(k) = rsy(k)*rsy(k)
      rszsq(k) = rsz(k)*rsz(k)
      if(rsx(k).eq.0.) then
         rsxsqinv(k) = 0.
      else
         rsxsqinv(k) = 1./rsxsq(k)
      endif
      if(rsy(k).eq.0.) then
         rsysqinv(k) = 0.
      else
         rsysqinv(k) = 1./rsysq(k)
      endif
      if(rsz(k).eq.0.) then
         rszsqinv(k) = 0.
      else
         rszsqinv(k) = 1./rszsq(k)
      endif
c
c----    
c
       isxreg(k) = kgeom(sx_cent(k),sy_cent(k),sz_cent(k))
c
c---- 
c
      write(*,30)k,isxtype(k),isxreg(k),isxev(k), sx_cent(k),
     + sy_cent(k),sz_cent(k),rsx(k),rsxsq(k),rsy(k),rsysq(k),
     + rsz(k),rszsq(k),esource(k), thetaangsx(k), phiangsx(k) 
c
c---- find out the angular and geometric source type 
c
c------------- isxangtype:
c
c------------- 1:cone-beam source (isotropic)
c------------- 2:fan-beam source (isotropic on a plane)
c
c
c------------- isxgeomtype:
c
c------------- 1:elliptic source
c------------- 2:cilindrical source //z axis
c ale-190401-- 3:cilindrical source //x axis
c ale-190401-- 4:box source
c
c
      isxangtype(k)  = int(isxtype(k)/10)
      isxgeomtype(k) = mod(isxtype(k),10)
      write(*,*) 'Source ang type:',isxangtype(k),
     +	 '  Source geom type:',isxgeomtype(k) 
      if (isxangtype(k).gt.2) then
         write(*,*)' Source ang type is not correct', isxangtype(k)
      endif 
      if (isxgeomtype(k).gt.4) then
         write(*,*)' Source geom type is not correct', isxgeomtype(k)
      endif 
c
      sx_inf(k) = sx_cent(k) - rsx(k)
      sx_sup(k) = 2. * rsx(k)
c
      sy_inf(k) = sy_cent(k) - rsy(k)
      sy_sup(k) = 2. * rsy(k)
c
      sz_inf(k) = sz_cent(k) - rsz(k)
      sz_sup(k) = 2. * rsz(k)
c
      write(*,40) k ,isxev(k),sx_inf(k),sx_sup(k)
     +, sy_inf(k),sy_sup(k), sz_inf(k),sz_sup(k),esource(k)
ccc      read(ltydat,70,err=20,end=25) thetaangsx(k), phiangsx(k)
c
c------  convert phi and theta in radiants
c
      phiangsx(k) = phiangsx(k) * pi / 180.
      thetaangsx(k) = thetaangsx(k) * pi / 180.
      write(*,*) ' '
      write(*,*)' PHI ANG ',phiangsx(k) ,' THETA ANG ',thetaangsx(k)
c
      cthsup(k) = 1 - cos (thetaangsx(k))
      cthinf(k) = 1.    
      write(*,*)' CTHSUP ', cthsup(k),' CTHINF ', cthinf(k)
      write(*,*) ' '
c
   10 continue
c
c
c--------------------------------------------------------c
c      
      return
c
c--------------------------------------------------------c
c
   20 write(*,*) ' ERROR IN GEOMETRY FILE - RUNSXINIT '
   25 write(*,*) ' END IN GEOMETRY FILE - RUNSXINIT '
      close(ltydat)
c
      return
c
   30 FORMAT( /,1X,'SX PARM,INPUT ' , 3I5,I11,/,1X,10(1x,F9.3))
   40 FORMAT( 1X,'SX PARM,TRASF ' ,I5,2X,I11,/,8(1x,F10.3))
   50 FORMAT(I10)
   60 FORMAT(I10,I10,11F10.0)
   70 FORMAT(2F10.0)
c
       END
c
       function kgeom(xpos,ypos,zpos)
c-------
c-----------------------------------------------------------c
c
      implicit none

       include 'egs4comm.for'
c
c
      integer k,iu1cell,iu2cell,iu3cell
      integer iu1temp,iu2temp,iu3temp,ind1temp
      integer ind2temp,ind3temp,ind4temp,indcell_temp
      real xpos,ypos,zpos,u2,u3
      real xpostemp,ypostemp,zpostemp
      doubleprecision xrel,yrel,xtemp0,ytemp0
      integer rawcell,colcell,cell
c
c-----------------------------------------------------------c
c
c          test for system  boundaries
c
c          cycle testing slabs
c
c           write (*,*) 'entro kgeom'
	   kgeom  = maxreg + 1 
c
c-------------- test if the position is inside the apparatus
c
           if ( (zmin.gt.zpos ) .or. (zmax.lt.zpos ) .or.
     a          (xmin.gt.xpos ) .or. (xmax.lt.xpos ) .or.
     b          (ymin.gt.ypos ) .or. (ymax.lt.ypos ) )then
c
                 kgeom  = maxreg +1
c           write (*,*) 'out of limits'

c
                 return
c
           endif
c	   
c--------------
c

cnico
	   dnear_temp = 0.
cnico	   
             if (zinitcoll.lt.zpos ) then
	        if (zendcoll.ge.zpos) then
	           if (xinfcoll.lt.xpos) then
	              if (xsupcoll.gt.xpos) then
	                 if (yinfcoll.lt.ypos) then
	                    if (ysupcoll.gt.ypos) then
ccc
ccc	particle inside the collimator
ccc
              if(modlabel.eq.1) then 		
c
c       moduli esagonali
c
		xrel = xpos - xinfcoll
		yrel = ypos - yinfcoll
		xtemp = 0.5 * xrel 
		ytemp = 0.8660254 * yrel
		u2 = xtemp + ytemp
		u3 = -xtemp + ytemp + u3offset
c	remind: xcellinv = 1/xcell_2
		iu1cell = dint(xrel*xcellinv)
		iu2cell = dint(u2*xcellinv)
		iu3cell = dint(u3*xcellinv) 
		iu1temp = mod(iu1cell,3)
		iu2temp = mod(iu2cell,3)
		iu3temp = mod(iu3cell,3)
		ind1temp = iu1cell + 1
		ind2temp = iu3cell - iu2cell
		ind3temp = abs(ind2temp)
		ind4temp = abs(ind2temp+1)
		if(label_hex.eq.0) then
		indcell_temp = regoffset + iu1cell * numcelly 
     a	        + 2*((iu2cell+1)-int((ind1temp+1)*.5)) + mod(ind4temp,2)
     c		+ 2*mod(ind1temp,2)*mod(ind3temp,2)
     		else
		indcell_temp = regoffset + iu1cell * numcelly 
     a	        + 2*((iu2cell+1)-int((ind1temp+1)*.5)) + mod(ind3temp,2)
     c		+ 2*mod(ind1temp,2)*mod(ind4temp,2)
		endif
c
	if(label_hex1.eq.0) then
		if(iu1temp.eq.0) then
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.0) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.0) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.1) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		elseif(iu1temp.eq.1) then     
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.1) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.0) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
			  return
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.0) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
			  return
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		else
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.0) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.1) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.0) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		endif
c     fine ciclo label_hex1 = 0 
	else if(label_hex1.eq.1) then 
		if(iu1temp.eq.0) then
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.1) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.1) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.2) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		elseif(iu1temp.eq.1) then     
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.2) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.1) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
			  return
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.1) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
			  return
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		else
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.1) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.2) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.1) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		endif
c     fine ciclo label_hex1 = 1 
	else
		if(iu1temp.eq.0) then
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.2) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.2) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.0) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		elseif(iu1temp.eq.1) then     
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.0) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.2) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
			  return
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.2) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
			  return
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		else
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.2) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.0) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.2) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		endif
	endif 
c     fine ciclo label_hex1 = 2 

c
c    fine moduli esagonali       
c
              elseif(modlabel.eq.2) then 
c
c     moduli rettangolari
c       
		xrel = xpos - xinfcoll
		yrel = ypos - yinfcoll
		xtemp0 = xrel * xcellinv
		ytemp0 = yrel * ycellinv
		rawcell = dint (xtemp0)
		colcell = dint (ytemp0)
		xtemp = xrel - rawcell * xcell
		ytemp = yrel - colcell * ycell
		cell = rawcell + 1 + colcell * numcellx
		if   (xtemp.gt.xair) then
			kgeom = regoffset + cell + numcell
			return
		elseif (ytemp.gt.yair) then
			kgeom = regoffset + cell + numcell + numcell
			return
		endif
		kgeom = regoffset + cell
		return
c
c    fine moduli rettangolari       
c
	      endif
                            endif
                         endif
                      endif
                   endif
                endif
             endif
ccc
ccc	particle outside the collimator
ccc
             if (zinit_det1mod.lt.zpos ) then
	        if (zend_det1mod.ge.zpos) then
	           if (xinf_det1mod.lt.xpos) then
	              if (xsup_det1mod.gt.xpos) then
	                 if (yinf_det1mod.lt.ypos) then
	                    if (ysup_det1mod.gt.ypos) then
ccc
ccc	particle inside the detector 
ccc
			xrel = xpos - xinf_det1mod
			yrel = ypos - yinf_det1mod
			xtemp0 = xrel * xcellinv_det1mod
			ytemp0 = yrel * ycellinv_det1mod
			rawcell = dint (xtemp0)
			colcell = dint (ytemp0)
			xtemp = xrel - rawcell * xcell_det1mod
			ytemp = yrel - colcell * ycell_det1mod
			cell = rawcell + 1 + colcell * numcellx_det1mod
cnico aggiunge x dead zone
		if   (xtemp.gt.xair_det1mod) then
			kgeom = regoffset_det1mod + cell + numcell_det1mod
			return
		elseif (ytemp.gt.yair_det1mod) then
			kgeom = regoffset_det1mod + cell + numcell_det1mod 
     a                          + numcell_det1mod
			return
		endif
cnico
			kgeom = regoffset_det1mod + cell
			return
                            endif
                         endif
                      endif
                   endif
                endif
             endif
ccc
ccc	particle outside the detector
ccc

             k = 0
c
   10        k=k+1
c
c-----     test for zpos
c
cale modifica per regioni cilindriche 09 04 2001
c              write (*,*) 'iregtype', iregtype(k)
      if(iregtype(k).eq.1)then
c                              write (*,*) 'inside parall1'
c      
cale: regioni parallelepipedi  
             if (zinit(k).lt.zpos ) then
	        if (zend(k).ge.zpos) then
	           if (xinf(k).lt.xpos) then
	              if (xsup(k).ge.xpos) then
	                 if (yinf(k).lt.ypos) then
	                    if (ysup(k).ge.ypos) then
c-----
c
c   this is the actual volume of detector space
c   and exit the slabs loop
c
                              kgeom = k
c                   
c                              write (*,*) 'inside parall2'
                              return
c
                            endif
                         endif
                      endif
                   endif
                endif
             endif
c
      elseif(iregtype(k).eq.2)then
c b
cale: regioni cilindriche parallele asse z
c123456 
       if(zinit(k).lt.zpos)then
       if(zend(k).ge.zpos)then
       if((((xsup(k)-xpos)*(xsup(k)-xpos))+((ypos-yinf(k))*       
     a (ypos-yinf(k)))).le.raysqrt(k))then
c        write(*,*) ((xpos-xsup(k))**2)+(ypos-yinf(k))**2
c         write(*,*) raysqrt(k)
                                                  kgeom = k
c
                                                  return
                                            endif
                                        endif
                                  endif
      else
c c
cale: regioni cilindriche parallele asse x

      if(xinf(k).lt.xpos)then
      if(xsup(k).ge.xpos)then
      if((((zpos-zend(k))*(zpos-zend(k)))+(ypos-yinf(k))*
     a (ypos-yinf(k))).le.raysqrt(k))then
                                                  kgeom = k
c
                                                  return
                                            endif
                                        endif
                                  endif
         endif
cale: end if iniziale
c
c       try next slab
c-------
c
          if( nplan - k ) 20,20,10
c
   20     continue
ccc	else
c
c     no region found : reject immediate
c
ccc          kgeom   =  maxreg +1
c
ccc          return
c
       end
      subroutine histin
c
c-----------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
c
      character*80 title, indstring, tuple_title
      integer k_loop, id_seq, lpage, nchan_x, nchan_y, kt
      real x_min, x_max, y_min, y_max, valmax
c
      data tuple_title / ' SLAB DETECTOR TUPLE 20  '/
      data k_loop, nchan_x, nchan_y / 0 , 0 , 0 /
      data x_min, x_max, y_min, y_max, valmax / 5 * 0.0 /
      data id_seq, k_loop, lpage / 0 , 0 , 60 /
c
cnico 
cc      open(unit=outtxt,file='out_txt',form='formatted',status='new')
cc      write(outtxt,*) 'x	y' 
cnico 

c-------------------------------------------------------c
c         set parameter for hbook package --------------c
c
      call hlimit( maxmem )
c
      open (22, file='output_histin')
c
      call houtpu( 22 )
      call hermes( 22 )
c
      call hpagsz( lpage  )
c
c----------------------------------------------------------
c
c
      write(*,50)
c
c----------------------------------------------------------------------c
c
      open( unit=ltyfil,file='xtest.hist'
     a,         form='formatted',status='old'     )
c
c----------------------------------------------------------------------c
c
      write(*,*) 'opening histograms file',
     a                           hislun,histofilename
c
c-----------------------------------------------------------
      open(unit=hislun,file=histofilename,
     a   form='unformatted',recl=1024,access='direct',status='new')
c
      call hropen(hislun,'x_ray',histofilename,'N',1024,istat)
      call hcdir (xray_memdir,'r')
      write ( * , 90) istat,xray_memdir
c
      k_loop=0
c
   10 read(ltyfil,110,err=30,end=40) id_seq,nchan_x,x_min,x_max,
     +                                  nchan_y,y_min,y_max,valmax
c
      if(id_seq.gt.0.and.k_loop.lt.20) then
         k_loop=k_loop+1
         read(ltyfil,150,err=30,end=40) title
         read(ltyfil,150,err=30,end=40) indstring
         call hbook1(id_seq,title,nchan_x,x_min,x_max,valmax)
         id_hist(k_loop)=id_seq
         write(*,*) id_seq,nchan_x,x_min,x_max
     *,                           nchan_y,y_min,y_max,valmax
         write(*,*) title
         write(*,*) indstring
         goto 10
      endif
c
      histo_count = k_loop
      k_loop = 0
c
   20 read(ltyfil,110,err=30,end=40) id_seq,nchan_x,x_min,x_max,
     +                                 nchan_y,y_min,y_max,valmax
c
      if(id_seq.gt.0.and.k_loop.lt.20) then
c
         k_loop=k_loop+1
         read(ltyfil,150,err=30,end=40) title
         read(ltyfil,150,err=30,end=40) indstring
         call hbook2(id_seq,title,nchan_x,x_min,x_max
     a,                             nchan_y,y_min,y_max,valmax)
         write(*,*) id_seq,nchan_x,x_min,x_max
     a,                     nchan_y,y_min,y_max,valmax
         write(*,*) title
         write(*,*) indstring
         bd_hist(k_loop)=id_seq
         goto 20
      endif
      plot_count = k_loop
c
c---------
c
      write(*,130) nt_hist,size_of_ntup,num_of_ntup, (kt,
     +   nt_tags (kt),kt=1,size_of_ntup)
c
      call hbookn(nt_hist,tuple_title,size_of_ntup
     a,  xray_memdir,num_of_ntup,nt_tags )
c
c---------------------------------------------------c
 
c
      write(*,60)
c
      close(ltyfil)
c
c-------------------------------------------------------c
c
      return
c
c-----  error section  ---------------------------------c
c
   30 write(*,70)
      close(ltyfil)
      return
c
   40 write(*,80)
      close(ltyfil)
      return
c
c------------------------------------------------------------------c
c 
   50 FORMAT(1X,'   ENTERING HISTIN  ')
   60 FORMAT(1X,'   END OF HISTIN  ')
   70 FORMAT(1H ,' ERROR IN READING FILE LTYFIL ',//)
   80 FORMAT(1H ,' EOF IN READING FILE LTYFIL ',// )
   90 FORMAT(1X,'STATUS',I10,2X,'X RAY DIRECTORY ' , A64)
  100 FORMAT(I5,I5,2F10.0,10X,20X,F10.0)
  110 FORMAT(I10,I10,2F10.0,I10,2F10.0,F10.0)
  130 FORMAT(1X , ' NTUPLE ID ',I8,' SIZE OF NTUP ' ,I8, ' NUM OF NTUP '
     +, I8 , /,' TUPLE TAG ' , (1H ,I4,2X,A10) )
  150 FORMAT(A80)
c
c
  180 FORMAT(A80)
  120 FORMAT(3I10)
  140 FORMAT(80A1)
  160 FORMAT(25I2)
  170 FORMAT( 1H , ' IAUSFL( ' , I3 , ' ) ' , 2X,I3)
c
      END
c
c------------------------------------------------------------------c
c
      subroutine histout
c
      implicit none
      include 'egs4comm.for'
      integer k_loop
c
      write(*,*) ' entering histout '
      call histdo
c
      do 10 k_loop = 1, histo_count
         call hrout(id_hist(k_loop),id_cycle(k_loop), ' ')
cale         write(*,*)' mono',id_hist(k_loop)
   10 continue
c
      do 20 k_loop=1,plot_count
         call hrout(bd_hist(k_loop),bd_cycle(k_loop), ' ')
cale         write(*,*)' bidi',bd_hist(k_loop)
   20 continue
c
      call hcdir (xray_memdir,'r')
      write ( * , 30) istat,xray_memdir
c
      call hrout ( 0,nt_cycle,'T')
      write( * ,*)' TUPLE ',nt_hist,nt_ident,nt_cycle
c
      call hldir(xray_memdir,'t')
      call hrend(xray_memdir)
c
      close( hislun )
      write(*,*) ' HISTO FILE CLOSED '
cc      close (outtxt)
      return
c
   30 FORMAT(1X,'STATUS',I10,2X,'X RAY DIRECTORY ' , A64)
      end
c
      SUBROUTINE HATCH
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
cale: egs4comm.for per pathname
      include 'egs4comm.for'
      character*80 Matfilename

      character*80 totMatfilename
cale

c
      CHARACTER*4 MBUF(72),MDLABL(8)
      real*4 ACD ,  ADEV ,  ASD ,  COST ,  CTHET ,  DEL ,  DFACT ,  DFAC
     *TI,  DUNITO,  DUNITR,  FNSSS ,  P ,  PZNORM,  RDEV ,  S2C2 ,  S2C2
     *MN,  S2C2MX,  SINT ,  SX ,  SXX ,  SXY ,   SY ,   WID ,  XS ,  XS0
     * ,  XS1 ,  XSI ,  WSS ,  YS ,  ZEROS(3)
      integer*4 I ,  I1ST ,  IB ,  ID ,  IE ,  IL ,  IM ,  IRAYL ,  IRN
     *,  ISTEST,  ISUB ,  ISS ,  IZ ,   IZZ ,  J ,  JR ,  LCTHET,  LMDL
     *,  LMDN ,  LTHETA,  MD ,  MXSINC,  NCMFP ,   NEKE ,   NGE ,   NGRI
     *M ,  NISUB ,  NLEKE ,    NM ,  NRANGE,    NRNA ,  NSEKE ,   NSGE ,
     *   NSINSS,  LOK(21)
      DATA MDLABL/' ','M','E','D','I','U','M','='/,LMDL/8/,LMDN/24/,DUNI        
     *TO/1./
      DATA I1ST/1/,NSINSS/37/,MXSINC/1002/,ISTEST/0/,NRNA/1000/
3240  FORMAT(1X,14I5)
3250  FORMAT(1X,1PE14.5,4E14.5)
3260  FORMAT(72A1)
      IF ((I1ST.NE.0)) THEN
        I1ST=0
          DO 3271 J=1,numgeom
          IF ((SMAXIR(J).LE.0.0)) THEN
            SMAXIR(J)=1E10
          END IF
          IF ((ESTEPR(J).LE.0.0)) THEN
            ESTEPR(J)=1.0
          END IF
3271    CONTINUE
3272    CONTINUE
        PRM=RM
        PRMT2=2.D0*PRM
        PZERO=0.0D0
        NISUB=MXSINC-2
        FNSSS=NSINSS
        WID=PI5D2/FLOAT(NISUB)
        WSS=WID/(FNSSS-1.0)
        ZEROS(1)=0.
        ZEROS(2)=PI
        ZEROS(3)=TWOPI
          DO 3281 ISUB=1,MXSINC
          SX=0.
          SY=0.
          SXX=0.
          SXY=0.
          XS0=WID*FLOAT(ISUB-2)
          XS1=XS0+WID
          IZ=0
            DO 3291 IZZ=1,3
            IF (((XS0.LE.ZEROS(IZZ)).AND.(ZEROS(IZZ).LE.XS1))) THEN
              IZ=IZZ
              GO TO3292
            END IF
3291      CONTINUE
3292      CONTINUE
          IF ((IZ.EQ.0)) THEN
            XSI=XS0
          ELSE
            XSI=ZEROS(IZ)
          END IF
            DO 3301 ISS=1,NSINSS
            XS=WID*FLOAT(ISUB-2)+WSS*FLOAT(ISS-1)-XSI
            YS=SIN(XS+XSI)
            SX=SX+XS
            SY=SY+YS
            SXX=SXX+XS*XS
            SXY=SXY+XS*YS
3301      CONTINUE
3302      CONTINUE
          IF ((IZ.NE.0)) THEN
            SIN1(ISUB)=SXY/SXX
            SIN0(ISUB)=-SIN1(ISUB)*XSI
          ELSE
            DEL=FNSSS*SXX-SX*SX
            SIN1(ISUB)=(FNSSS*SXY-SY*SX)/DEL
            SIN0(ISUB)=(SY*SXX-SX*SXY)/DEL - SIN1(ISUB)*XSI
          END IF
3281    CONTINUE
3282    CONTINUE
        SINC0=2.0
        SINC1=1.0/WID
        IF ((ISTEST.NE.0)) THEN
          ADEV=0.
          RDEV=0.
          S2C2MN=10.
          S2C2MX=0.
            DO 3311 ISUB=1,NISUB
              DO 3321 ISS=1,NSINSS
              THETA=WID*FLOAT(ISUB-1)+WSS*FLOAT(ISS-1)
              CTHET=PI5D2-THETA
              SINTHE=sin(THETA)
              COSTHE=sin(CTHET)
              SINT=SIN(THETA)
              COST=COS(THETA)
              ASD=ABS(SINTHE-SINT)
              ACD=ABS(COSTHE-COST)
              ADEV=max(ADEV,ASD,ACD)
              IF((SINT.NE.0.0))RDEV=max(RDEV,ASD/ABS(SINT))
              IF((COST.NE.0.0))RDEV=max(RDEV,ACD/ABS(COST))
              S2C2=SINTHE**2+COSTHE**2
              S2C2MN=min(S2C2MN,S2C2)
              S2C2MX=max(S2C2MX,S2C2)
              IF ((ISUB.LT.11)) THEN
                WRITE(6,3330)THETA,SINTHE,SINT,COSTHE,COST
3330            FORMAT(1PE20.7,4E20.7)
              END IF
3321        CONTINUE
3322        CONTINUE
3311      CONTINUE
3312      CONTINUE
          WRITE(6,3340)MXSINC,NSINSS
3340      FORMAT(' SINE TESTS,MXSINC,NSINSS=',2I5)                              
          WRITE(6,3350)ADEV,RDEV,S2C2MN,S2C2MX
3350      FORMAT(' ADEV,RDEV,S2C2(MN,MX) =',1PE16.8,3E16.8)                     
          ADEV=0.
          RDEV=0.
          S2C2MN=10.
          S2C2MX=0.
            DO 3361 IRN=1,NRNA
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            THETA = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            THETA=THETA*PI5D2
            CTHET=PI5D2-THETA
            SINTHE=sin(THETA)
            COSTHE=sin(CTHET)
            SINT=SIN(THETA)
            COST=COS(THETA)
            ASD=ABS(SINTHE-SINT)
            ACD=ABS(COSTHE-COST)
            ADEV=max(ADEV,ASD,ACD)
            IF((SINT.NE.0.0))RDEV=max(RDEV,ASD/ABS(SINT))
            IF((COST.NE.0.0))RDEV=max(RDEV,ACD/ABS(COST))
            S2C2=SINTHE**2+COSTHE**2
            S2C2MN=min(S2C2MN,S2C2)
            S2C2MX=max(S2C2MX,S2C2)
3361      CONTINUE
3362      CONTINUE
          WRITE(6,3370)NRNA
3370      FORMAT(' TEST AT ',I7,' RANDOM ANGLES IN (0,5*PI/2)')                 
          WRITE(6,3380)ADEV,RDEV,S2C2MN,S2C2MX
3380      FORMAT(' ADEV,RDEV,S2C2(MN,MX) =',1PE16.8,3E16.8)                     
        END IF
        P=1.
          DO 3391 I=1,50
          PWR2I(I)=P
          P=P/2.
3391    CONTINUE
3392    CONTINUE
      END IF
        DO 3401 J=1,NMED
3410    CONTINUE
          DO 3411 I=1,numgeom
          IF ((IRAYLR(I).EQ.1.AND.MED(I).EQ.J)) THEN
            IRAYLM(J)=1
            GO TO 3412
          END IF
3411    CONTINUE
3412    CONTINUE
3401  CONTINUE
3402  CONTINUE
      REWIND KMPI
cale NB: fort.12 e' il file dei materiali by PEGS4

      Matfilename='Mat900Kev.dat'

      totMatfilename=pathname(1:pathlength) // Matfilename
      write(*,*) 'cale',pathname, totMatfilename
      OPEN (UNIT=KMPI,FILE=totMatfilename, STATUS='OLD')                             
      NM=0
        DO 3421 IM=1,NMED
        LOK(IM)=0
        IF ((IRAYLM(IM).EQ.1)) THEN
          WRITE(6,3430)IM
3430      FORMAT(' RAYLEIGH OPTION REQUESTED FOR MEDIUM NUMBER',I3,/)           
        END IF
3421  CONTINUE
3422  CONTINUE
3440  CONTINUE
3441    CONTINUE
3450    CONTINUE
3451      CONTINUE
          READ(KMPI,3260,END=3460)MBUF
            DO 3471 IB=1,LMDL
            IF((MBUF(IB).NE.MDLABL(IB)))GO TO 3451
3471      CONTINUE
3472      CONTINUE
3480      CONTINUE
            DO 3481 IM=1,NMED
              DO 3491 IB=1,LMDN
              IL=LMDL+IB
              IF((MBUF(IL).NE.MEDIA(IB,IM)))GO TO 3481
              IF((IB.EQ.LMDN))GO TO 3452
3491        CONTINUE
3492        CONTINUE
3481      CONTINUE
3482      CONTINUE
        GO TO 3451
3452    CONTINUE
        IF((LOK(IM).NE.0))GO TO 3450
        LOK(IM)=1
        NM=NM+1
        READ(KMPI,1,ERR=3500) (MBUF(I),I=1,5),RHO(IM),NNE(IM),IUNRST(IM)
     *  ,EPSTFL(IM),IAPRIM(IM)
1       FORMAT(5A1,5X,F11.0,4X,I2,9X,I1,9X,I1,9X,I1)
        GO TO 3510
3500    BACKSPACE(KMPI)
        READ(KMPI,2)(MBUF(I),I=1,5),RHO(IM),NNE(IM),IUNRST(IM),EPSTFL(IM
     *  ), IAPRIM(IM)
2       FORMAT(5A1,5X,F11.0,4X,I2,26X,I1,9X,I1,9X,I1)
3510    CONTINUE
          DO 3511 IE=1,NNE(IM)
          READ(KMPI,3520)(MBUF(I),I=1,6),(ASYM(IM,IE,I),I=1,2), ZELEM(IM
     *    ,IE),WA(IM,IE),PZ(IM,IE),RHOZ(IM,IE)
3520      FORMAT (6A1,2A1,3X,F3.0,3X,F9.0,4X,F12.0,6X,F12.0)
3511    CONTINUE
3512    CONTINUE
        READ(KMPI,3250) RLC(IM),AE(IM),AP(IM),UE(IM),UP(IM)
        TE(IM)=AE(IM)-RM
        THMOLL(IM)=TE(IM)*2. + RM
        READ(KMPI,3240) MSGE(IM),MGE(IM),MSEKE(IM),MEKE(IM),MLEKE(IM),MC
     *  MFP(IM),MRANGE(IM),IRAYL
        NSGE=MSGE(IM)
        NGE=MGE(IM)
        NSEKE=MSEKE(IM)
        NEKE=MEKE(IM)
        NLEKE=MLEKE(IM)
        NCMFP=MCMFP(IM)
        NRANGE=MRANGE(IM)
        READ(KMPI,3250)(DL1(I,IM),DL2(I,IM),DL3(I,IM),DL4(I,IM),DL5(I,IM
     *  ),DL6(I,IM),I=1,6)
        READ(KMPI,3250)DELCM(IM),(ALPHI(I,IM),BPAR(I,IM),DELPOS(I,IM),I=
     *  1,2)
        READ(KMPI,3250)XR0(IM),TEFF0(IM),BLCC(IM),XCC(IM)
        READ(KMPI,3250)EKE0(IM),EKE1(IM)
        READ(KMPI,3250) (ESIG0(I,IM),ESIG1(I,IM),PSIG0(I,IM),PSIG1(I,IM)
     *  ,EDEDX0(I,IM),EDEDX1(I,IM),PDEDX0(I,IM),PDEDX1(I,IM),EBR10(I,IM)
     *  ,EBR11(I,IM),PBR10(I,IM),PBR11(I,IM),PBR20(I,IM),PBR21(I,IM),TMX
     *  S0(I,IM),TMXS1(I,IM),I=1,NEKE)
        READ(KMPI,3250)EBINDA(IM),GE0(IM),GE1(IM)
        READ(KMPI,3250)(GMFP0(I,IM),GMFP1(I,IM),GBR10(I,IM),GBR11(I,IM),
     *  GBR20(I,IM),GBR21(I,IM),I=1,NGE)
        IF ((IRAYLM(IM).EQ.1.AND.IRAYL.NE.1)) THEN
          WRITE(6,3530)IM
3530      FORMAT(' STOPPED IN HATCH: REQUESTED RAYLEIGH OPTION FOR MEDIU        
     *M',I3, /,' BUT RAYLEIGH DATA NOT INCLUDED IN DATA CREATED BY PEGS.        
     *')                                                                        
          STOP
        END IF
        IF ((IRAYL.EQ.1)) THEN
          READ(KMPI,3240) NGR(IM)
          NGRIM=NGR(IM)
          READ(KMPI,3250)RCO0(IM),RCO1(IM)
          READ(KMPI,3250)(RSCT0(I,IM),RSCT1(I,IM),I=1,NGRIM)
          READ(KMPI,3250)(COHE0(I,IM),COHE1(I,IM),I=1,NGE)
          IF ((IRAYLM(IM).NE.1)) THEN
            WRITE(6,3540)IM
3540        FORMAT(' RAYLEIGH DATA AVAILABLE FOR MEDIUM',I3, ' BUT OPTIO        
     *N NOT REQUESTED.',/)                                                      
          END IF
        END IF
        IF((NM.GE.NMED))GO TO3442
      GO TO 3441
3442  CONTINUE
      DUNITR=DUNIT
      IF ((DUNIT.LT.0.0)) THEN
        ID=MAX0(1,MIN0(5,IFIX(-DUNIT)))
        DUNIT=RLC(ID)
      END IF
      IF ((DUNIT.NE.1.0)) THEN
        WRITE(6,3550)DUNITR,DUNIT
3550    FORMAT(' DUNIT REQUESTED&USED ARE:',1PE14.5,E14.5,'(CM.)')              
      END IF
        DO 3561 IM=1,NMED
        DFACT=RLC(IM)/DUNIT
        DFACTI=1.0/DFACT
        I=1
          GO TO 3573
3571      I=I+1
3573      IF(I-(MEKE(IM)).GT.0)GO TO 3572
          ESIG0(I,IM)=ESIG0(I,IM)*DFACTI
          ESIG1(I,IM)=ESIG1(I,IM)*DFACTI
          PSIG0(I,IM)=PSIG0(I,IM)*DFACTI
          PSIG1(I,IM)=PSIG1(I,IM)*DFACTI
          EDEDX0(I,IM)=EDEDX0(I,IM)*DFACTI
          EDEDX1(I,IM)=EDEDX1(I,IM)*DFACTI
          PDEDX0(I,IM)=PDEDX0(I,IM)*DFACTI
          PDEDX1(I,IM)=PDEDX1(I,IM)*DFACTI
          TMXS0(I,IM)=TMXS0(I,IM)*DFACT
          TMXS1(I,IM)=TMXS1(I,IM)*DFACT
        GO TO 3571
3572    CONTINUE
        TEFF0(IM)=TEFF0(IM)*DFACT
        BLCC(IM)=BLCC(IM)*DFACTI
        XCC(IM)=XCC(IM)*SQRT(DFACTI)
        RLDU(IM)=RLC(IM)/DUNIT
        I=1
          GO TO 3583
3581      I=I+1
3583      IF(I-(MGE(IM)).GT.0)GO TO 3582
          GMFP0(I,IM)=GMFP0(I,IM)*DFACT
          GMFP1(I,IM)=GMFP1(I,IM)*DFACT
        GO TO 3581
3582    CONTINUE
3561  CONTINUE
3562  CONTINUE
      VACDST=VACDST*DUNITO/DUNIT
      DUNITO=DUNIT
        DO 3591 JR=1,numgeom
        MD=MED(JR)
        IF (((MD.GE.1).AND.(MD.LE.NMED))) THEN
          ECUT(JR)=max(ECUT(JR),AE(MD))
          PCUT(JR)=max(PCUT(JR),AP(MD))
          IF ((RHOR(JR).EQ.0.0)) THEN
            RHOR(JR)=RHO(MD)
          END IF
        END IF
3591  CONTINUE
3592  CONTINUE
      IF ((IBRDST.EQ.1)) THEN
          DO 3601 IM=1,NMED
          ZBRANG(IM)=0.0
          PZNORM=0.0
            DO 3611 IE=1,NNE(IM)
            ZBRANG(IM)= ZBRANG(IM)+PZ(IM,IE)*ZELEM(IM,IE)*(ZELEM(IM,IE)+
     *      1.0)
            PZNORM=PZNORM+PZ(IM,IE)
3611      CONTINUE
3612      CONTINUE
          ZBRANG(IM)=(8.116224E-05)*(ZBRANG(IM)/PZNORM)**(1./3.)
3601    CONTINUE
3602    CONTINUE
      END IF
      IF ((IPRDST.GT.0)) THEN
          DO 3621 IM=1,NMED
          ZBRANG(IM)=0.0
          PZNORM=0.0
            DO 3631 IE=1,NNE(IM)
            ZBRANG(IM)= ZBRANG(IM)+PZ(IM,IE)*ZELEM(IM,IE)*(ZELEM(IM,IE)+
     *      1.0)
            PZNORM=PZNORM+PZ(IM,IE)
3631      CONTINUE
3632      CONTINUE
          ZBRANG(IM)=(8.116224E-05)*(ZBRANG(IM)/PZNORM)**(1./3.)
3621    CONTINUE
3622    CONTINUE
      END IF
      call mscati
      call EDGSET(1,1)
      call init_compton
      call fix_brems
      IF (( ibr_nist .EQ. 1 )) THEN
        call init_nist_brems
      END IF
      IF ((NMED.EQ.1)) THEN
        WRITE(6,3640)
3640    FORMAT(' EGSnrc SUCCESSFULLY ''HATCHED'' FOR ONE MEDIUM.')              
      ELSE
        WRITE(6,3650)NMED
3650    FORMAT(' EGSnrc SUCCESSFULLY ''HATCHED'' FOR ',I5,' MEDIA.')            
      END IF
      CLOSE(UNIT=KMPI,STATUS='KEEP')                                            
      RETURN
3460  WRITE(6,3660)KMPI
3660  FORMAT(' END OF FILE ON UNIT ',I2,//, ' PROGRAM STOPPED IN HATCH B        
     *ECAUSE THE',/, ' FOLLOWING NAMES WERE NOT RECOGNIZED:',/)                 
        DO 3671 IM=1,NMED
        IF ((LOK(IM).NE.1)) THEN
          WRITE(6,3680)(MEDIA(I,IM),I=1,LMDN)
3680      FORMAT(40X,'''',24A1,'''')                                            
        END IF
3671  CONTINUE
3672  CONTINUE
      STOP
      END
      subroutine mscati
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      real*4 ededx,ei,eil,eip1,eip1l,si,sip1,eke_orig,elke_orig,aux,
     *ecutmn,tstbm,tstbmn
      real*4 p2,beta2,dedx0,ekef,elkef,estepx,ektmp,elktmp,chi_a2
      integer*4 i,leil,leip1l,neke,lelke,lelkef,lelktmp
      real*4 sigee,sigep,sig
      IF (( bca_algorithm .EQ. 0 )) THEN
        exact_bca = .true.
      ELSE
        exact_bca = .false.
      END IF
      IF (( estepe .LE. 0 .OR. estepe .GE. 1)) THEN
        estepe = 0.25
      END IF
      IF (( ximax .LE. 0 .OR. ximax .GE. 1 )) THEN
        IF (( exact_bca )) THEN
          ximax = 0.5
        ELSE
          ximax = 0.5
        END IF
      END IF
      IF ((transport_algorithm .NE. 0 .AND. transport_algorithm .NE. 1 .
     *AND. transport_algorithm .NE. 2 )) THEN
        transport_algorithm = 0
      END IF
      IF (( skindepth_for_bca .LE. 1 )) THEN
        IF (( transport_algorithm .EQ. 1 .AND. .NOT.exact_bca )) THEN
          write(6,*) ' PRESTA calculates default min. step-size for BCA:        
     * '                                                                        
          ecutmn = 1e30
            DO 3931 i=1,numgeom
            IF (( med(i) .GT. 0 .AND. med(i) .LE. nmed )) THEN
              ecutmn = Min(ecutmn,ecut(i))
            END IF
3931      CONTINUE
3932      CONTINUE
          write(6,*) '     minimum ECUT found: ',ecutmn                         
          tstbmn = 1e30
            DO 3941 medium=1,nmed
            tstbm = (ecutmn-0.5110034)*(ecutmn+0.5110034)/ecutmn**2
            tstbm = blcc(medium)*tstbm*(ecutmn/xcc(medium))**2
            tstbm = Log(tstbm/Log(tstbm))
            tstbmn = Min(tstbmn,tstbm)
3941      CONTINUE
3942      CONTINUE
          write(6,*) '     default BLCMIN is: ',tstbmn                          
          skindepth_for_bca = Exp(tstbmn)
          write(6,*) '     this corresponds to ',skindepth_for_bca, ' el        
     *astic MFPs '                                                              
        ELSE
          skindepth_for_bca = 3
        END IF
      END IF
      call init_ms_SR
        DO 3951 medium=1,nmed
        blcc(medium) = 1.16699413758864573*blcc(medium)
        xcc(medium) = xcc(medium)**2
3951  CONTINUE
3952  CONTINUE
      IF (( spin_effects )) THEN
        call init_spin
      END IF
      write(6,*)
      esige_max = 0
      psige_max = 0
        DO 3961 medium=1,nmed
        sigee = 1E-15
        sigep = 1E-15
        neke = meke(medium)
          DO 3971 i=1,neke
          ei = exp((float(i) - eke0(medium))/eke1(medium))
          eil = log(ei)
          leil = i
          ededx=ededx1(Leil,MEDIUM)*eil+ededx0(Leil,MEDIUM)
          sig=esig1(Leil,MEDIUM)*eil+esig0(Leil,MEDIUM)
          sig = sig/ededx
          IF((sig .GT. sigee))sigee = sig
          ededx=pdedx1(Leil,MEDIUM)*eil+pdedx0(Leil,MEDIUM)
          sig=psig1(Leil,MEDIUM)*eil+psig0(Leil,MEDIUM)
          sig = sig/ededx
          IF((sig .GT. sigep))sigep = sig
3971    CONTINUE
3972    CONTINUE
cale        write(6,*) ' Medium ',medium,' sige = ',sigee,sigep                     
        esig_e(medium) = sigee
        psig_e(medium) = sigep
        IF((sigee .GT. esige_max))esige_max = sigee
        IF((sigep .GT. psige_max))psige_max = sigep
3961  CONTINUE
3962  CONTINUE
      write(6,*)
      write(6,*) ' Initializing tmxs for estepe = ',estepe,' and ximax =        
     * ',ximax                                                                  
      write(6,*)
      rm = 0.5110034
        DO 3981 medium=1,nmed
        ei = exp((1 - eke0(medium))/eke1(medium))
        eil = log(ei)
        leil = 1
        E_array(1,medium) = ei
        expeke1(medium) = Exp(1./eke1(medium))-1
        range_ep(0,1,medium) = 0
        range_ep(1,1,medium) = 0
        neke = meke(medium)
          DO 3991 i=1,neke - 1
          eip1 = exp((float(i + 1) - eke0(medium))/eke1(medium))
          E_array(i+1,medium) = eip1
          eke_orig = 0.5*(eip1+ei)
          elke_orig = Log(eke_orig)
          Lelke=eke1(MEDIUM)*elke_orig+eke0(MEDIUM)
          ededx=pdedx1(Lelke,MEDIUM)*elke_orig+pdedx0(Lelke,MEDIUM)
          aux = pdedx1(i,medium)/ededx
          range_ep(1,i+1,medium) = range_ep(1,i,medium) + (eip1-ei)/eded
     *    x*(1+aux*(1+2*aux)*((eip1-ei)/eke_orig)**2/24)
          ededx=ededx1(Lelke,MEDIUM)*elke_orig+ededx0(Lelke,MEDIUM)
          aux = ededx1(i,medium)/ededx
          range_ep(0,i+1,medium) = range_ep(0,i,medium) + (eip1-ei)/eded
     *    x*(1+aux*(1+2*aux)*((eip1-ei)/eke_orig)**2/24)
          ei = eip1
3991    CONTINUE
3992    CONTINUE
        eil = (1 - eke0(medium))/eke1(medium)
        ei = Exp(eil)
        leil = 1
        p2 = ei*(ei+2*rm)
        beta2 = p2/(p2+rm*rm)
        chi_a2 = Xcc(medium)/(4*p2*blcc(medium))
        dedx0=ededx1(Leil,MEDIUM)*eil+ededx0(Leil,MEDIUM)
        estepx = 2*p2*beta2*dedx0/ei/Xcc(medium)/(Log(1+1./chi_a2)*(1+ch
     *  i_a2)-1)
        estepx = estepx*ximax
        IF (( estepx .GT. estepe )) THEN
          estepx = estepe
        END IF
        si = estepx*ei/dedx0
          DO 4001 i=1,neke - 1
          elke_orig = (i + 1 - eke0(medium))/eke1(medium)
          eke_orig = Exp(elke_orig)
          lelke = i+1
          p2 = eke_orig*(eke_orig+2*rm)
          beta2 = p2/(p2+rm*rm)
          chi_a2 = Xcc(medium)/(4*p2*blcc(medium))
          ededx=ededx1(Lelke,MEDIUM)*elke_orig+ededx0(Lelke,MEDIUM)
          estepx = 2*p2*beta2*ededx/eke_orig/ Xcc(medium)/
     *    (Log(1+1./chi_a2)*(
     *    1+chi_a2)-1)
          estepx = estepx*ximax
          IF (( estepx .GT. estepe )) THEN
            estepx = estepe
          END IF
          ekef = (1-estepx)*eke_orig
          IF (( ekef .LE. E_array(1,medium) )) THEN
            sip1 = (E_array(1,medium) - ekef)/dedx0
            ekef = E_array(1,medium)
            elkef = (1 - eke0(medium))/eke1(medium)
            lelkef = 1
          ELSE
            elkef = Log(ekef)
            Lelkef=eke1(MEDIUM)*elkef+eke0(MEDIUM)
            leip1l = lelkef + 1
            eip1l = (leip1l - eke0(medium))/eke1(medium)
            eip1 = E_array(leip1l,medium)
            aux = (eip1 - ekef)/eip1
            elktmp = 0.5*(elkef+eip1l+0.25*aux*aux*(1+aux*(1+0.875*aux))
     *      )
            ektmp = 0.5*(ekef+eip1)
            lelktmp = lelkef
            ededx=ededx1(Lelktmp,MEDIUM)*elktmp+ededx0(Lelktmp,MEDIUM)
            aux = ededx1(lelktmp,medium)/ededx
            sip1 = (eip1 - ekef)/ededx*( 1+aux*(1+2*aux)*((eip1-ekef)/ek
     *      tmp)**2/24)
          END IF
          sip1 = sip1 + range_ep(0,lelke,medium) - range_ep(0,lelkef+1,m
     *    edium)
          tmxs1(i,medium) = (sip1 - si)*eke1(medium)
          tmxs0(i,medium) = sip1 - tmxs1(i,medium)*elke_orig
          si = sip1
4001    CONTINUE
4002    CONTINUE
        tmxs0(neke,medium) = tmxs0(neke - 1,medium)
        tmxs1(neke,medium) = tmxs1(neke - 1,medium)
3981  CONTINUE
3982  CONTINUE
      return
      end
      subroutine init_ms_SR
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      integer*4 i,j,k
cale
      include 'egs4comm.for'
      include 'egs4fcomm.for'

      character*80 msfilename

      character*80 totmsfilename

      msfilename='msnew.dat'

      totmsfilename=pathname(1:pathlength) // msfilename
cale     write(*,*)'cale:',totmsfilename
      open(unit=11,file=totmsfilename,form='formatted',status='old')
cale

      write(6,'(a,$)') '  Reading screened Rutherford MS data ..........        
     *..... '                                                                   
        DO 4061 i=0,63
          DO 4071 j=0,7
          read(11,*) (ums_array(i,j,k),k=0,31)
          read(11,*) (fms_array(i,j,k),k=0,31)
          read(11,*) (wms_array(i,j,k),k=0,31-1)
          read(11,*) (ims_array(i,j,k),k=0,31-1)
            DO 4081 k=0,31-1
            fms_array(i,j,k) = fms_array(i,j,k+1)/fms_array(i,j,k)-1
            ims_array(i,j,k) = ims_array(i,j,k)-1
4081      CONTINUE
4082      CONTINUE
          fms_array(i,j,31)=fms_array(i,j,31-1)
4071    CONTINUE
4072    CONTINUE
4061  CONTINUE
4062  CONTINUE
      write(6,'(a)') ' done '                                                   
      llammin = Log(1.)
      llammax = Log(1e5)
      dllamb = (llammax-llammin)/63
      dllambi = 1./dllamb
      dqms = 0.5/7
      dqmsi = 1./dqms
      return
      end
      subroutine init_spin
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c

      real*4 eta_array(0:1,0: 31), c_array(0:1,0: 31),g_array(0:1,0: 31)
     *, earray(0: 31),tmp_array(0: 31), sum_Z2,sum_Z,sum_A,sum_pz,Z_orig
     *,tmp,Z23,g_m,g_r,sig,dedx, dum1,dum2,dum3,aux_o,tau,tauc,beta2,eta
     *,gamma,fmax, eil,e_orig,si1e,si2e,si1p,si2p,aae,etap, elarray
     *(0: 31),farray(0: 31), af(0: 31),bf(0: 31),cf(0: 31), df(0: 31),
     * spline
      integer*4 medium_orig,iq_orig,i,j,k,i_ele,iii,iZ,iiZ,n_ener,n_q,
     *n_point,je,neke, ndata,leil,length

      character*80 spin_file
cale
      include 'egs4comm.for'

      character*80 prespin_file
cale
      character*6 string
      integer*4 lnblnk1
      real*4 fine,TF_constant
      parameter (fine=137.03604,TF_constant=0.88534138)
cale      call getenv('HEN_HOUSE',spin_file)
      prespin_file='spinms/z000'
      spin_file=pathname(1:pathlength)//prespin_file
cale      write(*,*)'cale:',spin_file
cale
      length = lnblnk1(spin_file)
c
c mirko - modified spin_files location
c
cale      IF (( spin_file(length:length) .NE. 's' )) THEN                           
cale        length = length + 1
cale        spin_file(length:length) = 's'                                          
cale      END IF
cale      spin_file(length+1:length+11) = 'pinms/z000'                             
cale      length = lnblnk1(spin_file)
c
c mirko - end of modified spin_files location
c
cale


        write(6,*) '  Initializing spin data for media '

        DO 4091 medium_orig=1,NMED
cale        write(6,'(a,i4,a,$)') '  Initializing spin data for
cale     * medium_orig ', medium_orig, ' ..................... '
          DO 4101 iq_orig=0,1
            DO 4111 i=0, 31
            eta_array(iq_orig,i)=0
            c_array(iq_orig,i)=0
            g_array(iq_orig,i)=0
              DO 4121 j=0,15
                DO 4131 k=0,31
                spin_rej(medium_orig,iq_orig,i,j,k) = 0
4131          CONTINUE
4132          CONTINUE
4121        CONTINUE
4122        CONTINUE
4111      CONTINUE
4112      CONTINUE
4101    CONTINUE
4102    CONTINUE
        sum_Z2=0
        sum_A=0
        sum_pz=0
        sum_Z=0
          DO 4141 i_ele=1,NNE(medium_orig)
          Z_orig = ZELEM(medium_orig,i_ele)
          iZ = int(Z_orig+0.5)
          tmp = PZ(medium_orig,i_ele)*Z_orig*(Z_orig+1)
          iii = iZ/100
          spin_file(length-2:length-2) = char(iii+48)
          iiZ = iZ - iii*100
          iii = iiZ/10
          spin_file(length-1:length-1) = char(iii+48)
          iiZ = iiZ - 10*iii
          spin_file(length:length) = char(iiZ+48)
          open(61,file=spin_file,status='old',err=4150)   
          read(61,*) espin_min,espin_max,b2spin_min,b2spin_max
          read(61,*) n_ener,n_q,n_point
          IF (( n_ener .NE. 15 .OR. n_q .NE. 15 .OR. n_point .NE. 31)) 
     *    THEN
            write(6,*) ' Wrong spin file for Z = ',iZ                           
            stop
          END IF
          sum_Z2 = sum_Z2 + tmp
          sum_Z = sum_Z + PZ(medium_orig,i_ele)*Z_orig
          sum_A = sum_A + PZ(medium_orig,i_ele)*WA(medium_orig,i_ele)
          sum_pz = sum_pz + PZ(medium_orig,i_ele)
          Z23 = Z_orig**0.6666667
            DO 4161 iq_orig=0,1
            read(61,*)
            read(61,*)
              DO 4171 i=0, 31
              read(61,'(a,g14.6)') string,earray(i)                             
              read(61,*) dum1,dum2,dum3,aux_o
              eta_array(iq_orig,i)=eta_array(iq_orig,i)+tmp*Log
     * (Z23*aux_o)
              tau = earray(i)/511.0034
              beta2 = tau*(tau+2)/(tau+1)**2
              eta = Z23/(fine*TF_constant)**2*aux_o/4/tau/(tau+2)
              c_array(iq_orig,i)=c_array(iq_orig,i)+ tmp*(Log(1+1/eta)
     *        -1/(1+eta))*dum1*dum3
              g_array(iq_orig,i)=g_array(iq_orig,i)+tmp*dum2
                DO 4181 j=0,15
                read(61,*) tmp_array
                  DO 4191 k=0,31
                  spin_rej(medium_orig,iq_orig,i,j,k) = spin_rej(
     *            medium_orig,iq_orig,i,j,k)+ tmp*tmp_array(k)
4191            CONTINUE
4192            CONTINUE
4181          CONTINUE
4182          CONTINUE
4171        CONTINUE
4172        CONTINUE
4161      CONTINUE
4162      CONTINUE
          close(61)
4141    CONTINUE
4142    CONTINUE
          DO 4201 iq_orig=0,1
            DO 4211 i=0, 31
              DO 4221 j=0,15
              fmax = 0
                DO 4231 k=0,31
                IF (( spin_rej(medium_orig,iq_orig,i,j,k) .GT. fmax ))
     * THEN
                  fmax = spin_rej(medium_orig,iq_orig,i,j,k)
                END IF
4231          CONTINUE
4232          CONTINUE
                DO 4241 k=0,31
                spin_rej(medium_orig,iq_orig,i,j,k) = spin_rej(medium
     *          _orig,iq_orig,i,j,k)/fmax
4241          CONTINUE
4242          CONTINUE
4221        CONTINUE
4222        CONTINUE
4211      CONTINUE
4212      CONTINUE
4201    CONTINUE
4202    CONTINUE
          DO 4251 i=0, 31
          tau = earray(i)/511.0034
          beta2 = tau*(tau+2)/(tau+1)**2
            DO 4261 iq_orig=0,1
            aux_o = Exp(eta_array(iq_orig,i)/sum_Z2)/(fine*TF_constant)
     * **2
            eta_array(iq_orig,i) = 0.26112447*aux_o*blcc(medium_orig)
     * /xcc(medium_orig)
            eta = aux_o/4/tau/(tau+2)
            gamma = 3*(1+eta)*(Log(1+1/eta)*(1+2*eta)-2)/ (Log(1+1/eta)
     *      *(1+eta)-1)
            g_array(iq_orig,i) = g_array(iq_orig,i)/sum_Z2/gamma
            c_array(iq_orig,i) = c_array(iq_orig,i)/sum_Z2/(Log(1+1/eta)
     *      -1/(1+eta))
4261      CONTINUE
4262      CONTINUE
4251    CONTINUE
4252    CONTINUE
        espin_min = espin_min/1000
        espin_max = espin_max/1000
        dlener = Log(espin_max/espin_min)/15
        dleneri = 1/dlener
        espml = Log(espin_min)
        dbeta2 = (b2spin_max-b2spin_min)/15
        dbeta2i = 1/dbeta2
        eil = (1 - eke0(medium_orig))/eke1(medium_orig)
        e_orig = Exp(eil)
        IF (( e_orig .LE. espin_min )) THEN
          si1e = eta_array(0,0)
          si1p = eta_array(1,0)
        ELSE
          IF (( e_orig .LE. espin_max )) THEN
            aae = (eil-espml)*dleneri
            je = aae
            aae = aae - je
          ELSE
            tau = e_orig/0.5110034
            beta2 = tau*(tau+2)/(tau+1)**2
            aae = (beta2 - b2spin_min)*dbeta2i
            je = aae
            aae = aae - je
            je = je + 15 + 1
          END IF
          si1e = (1-aae)*eta_array(0,je) + aae*eta_array(0,je+1)
          si1p = (1-aae)*eta_array(1,je) + aae*eta_array(1,je+1)
        END IF
        neke = meke(medium_orig)
          DO 4271 i=1,neke - 1
          eil = (i+1 - eke0(medium_orig))/eke1(medium_orig)
          e_orig = Exp(eil)
          IF (( e_orig .LE. espin_min )) THEN
            si2e = eta_array(0,0)
            si2p = eta_array(1,0)
          ELSE
            IF (( e_orig .LE. espin_max )) THEN
              aae = (eil-espml)*dleneri
              je = aae
              aae = aae - je
            ELSE
              tau = e_orig/0.5110034
              beta2 = tau*(tau+2)/(tau+1)**2
              aae = (beta2 - b2spin_min)*dbeta2i
              je = aae
              aae = aae - je
              je = je + 15 + 1
            END IF
            si2e = (1-aae)*eta_array(0,je) + aae*eta_array(0,je+1)
            si2p = (1-aae)*eta_array(1,je) + aae*eta_array(1,je+1)
          END IF
          etae_ms1(i,medium_orig) = (si2e - si1e)*eke1(medium_orig)
          etae_ms0(i,medium_orig) = si2e - etae_ms1(i,medium_orig)*eil
          etap_ms1(i,medium_orig) = (si2p - si1p)*eke1(medium_orig)
          etap_ms0(i,medium_orig) = si2p - etap_ms1(i,medium_orig)*eil
          si1e = si2e
          si1p = si2p
4271    CONTINUE
4272    CONTINUE
        etae_ms1(neke,medium_orig) = etae_ms1(neke-1,medium_orig)
        etae_ms0(neke,medium_orig) = etae_ms0(neke-1,medium_orig)
        etap_ms1(neke,medium_orig) = etap_ms1(neke-1,medium_orig)
        etap_ms0(neke,medium_orig) = etap_ms0(neke-1,medium_orig)
          DO 4281 i=0,15
          elarray(i) = Log(earray(i)/1000)
          farray(i) = c_array(0,i)
4281    CONTINUE
4282    CONTINUE
          DO 4291 i=15+1, 31-1
          elarray(i) = Log(earray(i+1)/1000)
          farray(i) = c_array(0,i+1)
4291    CONTINUE
4292    CONTINUE
        ndata =  31+1
        IF (( ue(medium_orig) .GT. 1e5 )) THEN
          elarray(ndata-1) = Log(ue(medium_orig))
        ELSE
          elarray(ndata-1) = Log(1e5)
        END IF
        farray(ndata-1) = 1
        call set_spline(elarray,farray,af,bf,cf,df,ndata)
        eil = (1 - eke0(medium_orig))/eke1(medium_orig)
        si1e = spline(eil,elarray,af,bf,cf,df,ndata)
          DO 4301 i=1,neke-1
          eil = (i+1 - eke0(medium_orig))/eke1(medium_orig)
          si2e = spline(eil,elarray,af,bf,cf,df,ndata)
          q1ce_ms1(i,medium_orig) = (si2e - si1e)*eke1(medium_orig)
          q1ce_ms0(i,medium_orig) = si2e - q1ce_ms1(i,medium_orig)*eil
          si1e = si2e
4301    CONTINUE
4302    CONTINUE
        q1ce_ms1(neke,medium_orig) = q1ce_ms1(neke-1,medium_orig)
        q1ce_ms0(neke,medium_orig) = q1ce_ms0(neke-1,medium_orig)
          DO 4311 i=0,15
          farray(i) = c_array(1,i)
4311    CONTINUE
4312    CONTINUE
          DO 4321 i=15+1, 31-1
          farray(i) = c_array(1,i+1)
4321    CONTINUE
4322    CONTINUE
        call set_spline(elarray,farray,af,bf,cf,df,ndata)
        eil = (1 - eke0(medium_orig))/eke1(medium_orig)
        si1e = spline(eil,elarray,af,bf,cf,df,ndata)
          DO 4331 i=1,neke-1
          eil = (i+1 - eke0(medium_orig))/eke1(medium_orig)
          si2e = spline(eil,elarray,af,bf,cf,df,ndata)
          q1cp_ms1(i,medium_orig) = (si2e - si1e)*eke1(medium_orig)
          q1cp_ms0(i,medium_orig) = si2e - q1cp_ms1(i,medium_orig)*eil
          si1e = si2e
4331    CONTINUE
4332    CONTINUE
        q1cp_ms1(neke,medium_orig) = q1cp_ms1(neke-1,medium_orig)
        q1cp_ms0(neke,medium_orig) = q1cp_ms0(neke-1,medium_orig)
          DO 4341 i=0,15
          farray(i) = g_array(0,i)
4341    CONTINUE
4342    CONTINUE
          DO 4351 i=15+1, 31-1
          farray(i) = g_array(0,i+1)
4351    CONTINUE
4352    CONTINUE
        call set_spline(elarray,farray,af,bf,cf,df,ndata)
        eil = (1 - eke0(medium_orig))/eke1(medium_orig)
        si1e = spline(eil,elarray,af,bf,cf,df,ndata)
          DO 4361 i=1,neke-1
          eil = (i+1 - eke0(medium_orig))/eke1(medium_orig)
          si2e = spline(eil,elarray,af,bf,cf,df,ndata)
          q2ce_ms1(i,medium_orig) = (si2e - si1e)*eke1(medium_orig)
          q2ce_ms0(i,medium_orig) = si2e - q2ce_ms1(i,medium_orig)*eil
          si1e = si2e
4361    CONTINUE
4362    CONTINUE
        q2ce_ms1(neke,medium_orig) = q2ce_ms1(neke-1,medium_orig)
        q2ce_ms0(neke,medium_orig) = q2ce_ms0(neke-1,medium_orig)
          DO 4371 i=0,15
          farray(i) = g_array(1,i)
4371    CONTINUE
4372    CONTINUE
          DO 4381 i=15+1, 31-1
          farray(i) = g_array(1,i+1)
4381    CONTINUE
4382    CONTINUE
        call set_spline(elarray,farray,af,bf,cf,df,ndata)
        eil = (1 - eke0(medium_orig))/eke1(medium_orig)
        si1e = spline(eil,elarray,af,bf,cf,df,ndata)
          DO 4391 i=1,neke-1
          eil = (i+1 - eke0(medium_orig))/eke1(medium_orig)
          si2e = spline(eil,elarray,af,bf,cf,df,ndata)
          q2cp_ms1(i,medium_orig) = (si2e - si1e)*eke1(medium_orig)
          q2cp_ms0(i,medium_orig) = si2e - q2cp_ms1(i,medium_orig)*eil
          si1e = si2e
4391    CONTINUE
4392    CONTINUE
        q2cp_ms1(neke,medium_orig) = q2cp_ms1(neke-1,medium_orig)
        q2cp_ms0(neke,medium_orig) = q2cp_ms0(neke-1,medium_orig)
        tauc = te(medium_orig)/0.5110034
        si1e = 1
          DO 4401 i=1,neke-1
          eil = (i+1 - eke0(medium_orig))/eke1(medium_orig)
          e_orig = Exp(eil)
          leil=i+1
          tau=e_orig/0.5110034
          IF (( tau .GT. 2*tauc )) THEN
            sig=esig1(Leil,medium_orig)*eil+esig0(Leil,medium_orig)
            dedx=ededx1(Leil,medium_orig)*eil+ededx0(Leil,medium_orig)
            sig = sig/dedx
            IF (( sig .GT. 1e-6 )) THEN
              etap=etae_ms1(Leil,medium_orig)*eil+etae_ms0(Leil,
     *  medium_orig)
              eta = 0.25*etap*xcc(medium_orig)/blcc(medium_orig)/
     *  tau/(tau+2)
              g_r = (1+2*eta)*Log(1+1/eta)-2
              g_m = Log(0.5*tau/tauc)+ (1+((tau+2)/(tau+1))**2)*Log(2*(t
     *        au-tauc+2)/(tau+4))- 0.25*(tau+2)*(tau+2+2*(2*tau+1)/(tau+
     *        1)**2)* Log((tau+4)*(tau-tauc)/tau/(tau-tauc+2))+ 0.5*(tau
     *        -2*tauc)*(tau+2)*(1/(tau-tauc)-1/(tau+1)**2)
              IF (( g_m .LT. g_r )) THEN
                g_m = g_m/g_r
              ELSE
                g_m = 1
              END IF
              si2e = 1 - g_m*sum_Z/sum_Z2
            ELSE
              si2e = 1
            END IF
          ELSE
            si2e = 1
          END IF
          blcce1(i,medium_orig) = (si2e - si1e)*eke1(medium_orig)
          blcce0(i,medium_orig) = si2e - blcce1(i,medium_orig)*eil
          si1e = si2e
4401    CONTINUE
4402    CONTINUE
        blcce1(neke,medium_orig) = blcce1(neke-1,medium_orig)
        blcce0(neke,medium_orig) = blcce0(neke-1,medium_orig)
cale        write(6,'(a)') ' done'                                                  
4091  CONTINUE
4092  CONTINUE
      return
4150  write(6,*) ' ******************** Error in init_spin *************        
     ******* '                                                                  
      write(6,'(a,a)') '  could not open file ',spin_file                       
      write(6,*) ' terminating execution '                                      
      write(6,*) ' *****************************************************        
     ********'                                                                  
      stop
      end
      SUBROUTINE EDGSET(NREGLO,NREGHI)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
cale
      include 'egs4comm.for'

      character*80 phabsfilename

      character*80 totphabsfilename

      character*80 photofilename

      character*80 totphotofilename

cale

c
      integer*4 NREGLO,NREGHI
      integer*4 i,j,k,jj,iz
      logical do_relax
      logical got_data
      save got_data
      data got_data/.false./
      IF((got_data))return
      WRITE(6,4660)
4660  FORMAT(' Output from subroutine EDGSET:'/ ' ======================        
     *========')                                                                
      do_relax = .false.
        DO 4671 jj=1,numgeom
        IZ=IEDGFL(JJ)
        IF (( iz .GT. 0 .AND. iz .LE. 100 )) THEN
          do_relax = .true.
          GO TO4672
        END IF
4671  CONTINUE
4672  CONTINUE
      IF (( .NOT.do_relax )) THEN
        WRITE(6,4680)
4680    FORMAT(' Atomic relaxations not requested! '/)                          
        return
      END IF
      WRITE(6,4690)
4690  FORMAT(/'  Atomic relaxations requested! ')                               
      WRITE(6,4700)
4700  FORMAT(/' Reading photo-absorption data .....',$)                         
      got_data = .true.
cale
      phabsfilename='photo_relax.dat'

      totphabsfilename=pathname(1:pathlength) // phabsfilename
      write(*,*)'cale:',totphabsfilename
      open(unit=77,file=totphabsfilename,form='formatted',status='old')
cale
cale      write(*,*)'cale: dopo open'
        DO 4711 i=1,100
        read(77,*) j,(binding_energies(k,i),k=1,6)
          DO 4721 k=1,6
          binding_energies(k,i) = binding_energies(k,i)*1e-6
4721    CONTINUE
4722    CONTINUE
4711  CONTINUE
4712  CONTINUE
      read(77,*)
        DO 4731 i=1,100
        read(77,*) j,(interaction_prob(k,i),k=1,5-1)
4731  CONTINUE
4732  CONTINUE
      WRITE(6,4740)
4740  FORMAT(' Done')                                                           
      WRITE(6,4750)
4750  FORMAT(' Reading relaxation data .... ',$)                                
      read(77,*)
        DO 4761 i=1,100
        read(77,*) j,(relaxation_prob(k,i),k=1,19)
4761  CONTINUE
4762  CONTINUE
      read(77,*)
        DO 4771 i=1,100
        read(77,*) j,(relaxation_prob(k,i),k=20,26)
4771  CONTINUE
4772  CONTINUE
      read(77,*)
        DO 4781 i=1,100
        read(77,*) j,(relaxation_prob(k,i),k=27,32)
4781  CONTINUE
4782  CONTINUE
      read(77,*)
        DO 4791 i=1,100
        read(77,*) j,(relaxation_prob(k,i),k=33,37)
4791  CONTINUE
4792  CONTINUE
      read(77,*)
        DO 4801 i=1,100
        read(77,*) j,relaxation_prob(38,i)
4801  CONTINUE
4802  CONTINUE
      WRITE(6,4810)
4810  FORMAT(' Done')                                                           
      WRITE(6,4820)
4820  FORMAT(' Reading photo cross section data .... ',$)
cale
      photofilename='photo_cs.dat'

      totphotofilename=pathname(1:pathlength) // photofilename
cale     write(*,*)'cale:',totphotofilename
      open(unit=79,file=totphotofilename,form='formatted',status='old')
cale                       
        DO 4831 i=1,100
        read(79,*) j,edge_number(i)
          DO 4841 j=1,edge_number(i)
          read(79,*) edge_a(j,i),edge_b(j,i),edge_c(j,i), edge_d(j,i),ed
     *    ge_energies(j,i)
4841    CONTINUE
4842    CONTINUE
4831  CONTINUE
4832  CONTINUE
      WRITE(6,4850)
4850  FORMAT(' Done')                                                           
      RETURN
      END
      subroutine init_compton
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
cale
      include 'egs4comm.for'

      character*80 incohfilename

      character*80 totincohfilename

cale

c
      integer*4 i,j,iz,nsh
      real*4 aux,pztot
      real*4 aux_erf,erf
        DO 3711 j=1,numgeom
        medium = med(j)
        IF (( medium .GT. 0 .AND. medium .LE. nmed)) THEN
          IF (( ibcmp(j) .GT. 0 )) THEN
            goto 3720
          END IF
        END IF
3711  CONTINUE
3712  CONTINUE
      WRITE(6,3730)
3730  FORMAT(' Bound Compton scattering not requested! '/)                      
      return
3720  WRITE(6,3740)
3740  FORMAT(' Bound Compton scattering requested, reading data ......',        
     *$)
cale
      incohfilename='incoh.dat'

      totincohfilename=pathname(1:pathlength) // incohfilename
cale     write(*,*)'cale:',totincohfilename
      open(unit=78,file=totincohfilename,form='formatted',status='old')
cale
        DO 3751 j=1,18
        read(78,*)
3751  CONTINUE
3752  CONTINUE
        DO 3761 j=1,1538
        read(78,*) iz_array(j),shn_array(j),ne_array(j), Jo_array(j),be_
     *  array(j)
        Jo_array(j) = Jo_array(j)*137.
        be_array(j) = be_array(j)*1e-6/rm
        aux_erf = 0.70710678119*(1+0.3*Jo_array(j))
        erfJo_array(j) = 0.82436063535*(Erf(aux_erf)-1)
3761  CONTINUE
3762  CONTINUE
      WRITE(6,3770)
3770  FORMAT(' Done')                                                           
      WRITE(6,3780)
3780  FORMAT(/' Initializing Bound Compton scattering ......')                  
        DO 3791 medium=1,nmed
        pztot = 0
        nsh = 0
          DO 3801 i=1,nne(medium)
          iz = int(zelem(medium,i))
            DO 3811 j=1,1538
            IF (( iz .EQ. iz_array(j) )) THEN
              nsh = nsh + 1
              IF (( nsh .GT. 50 )) THEN
                WRITE(6,3820)medium,50
3820            FORMAT(/' For medium ',i3,' the number of shells is > ',        
     *          i4,'!')                                                         
                WRITE(6,3830)
3830            FORMAT(' Increase the parameter $MXMDSH! ')                     
                stop
              END IF
              shell_array(nsh,medium) = j
              aux = pz(medium,i)*ne_array(j)
              eno_array(nsh,medium) = aux
              pztot = pztot + aux
            END IF
3811      CONTINUE
3812      CONTINUE
3801    CONTINUE
3802    CONTINUE
        IF (( nsh .EQ. 0 )) THEN
          WRITE(6,3840)medium
3840      FORMAT(' Medium ',i3,' has zero shells! ')                            
          stop
        END IF
        n_shell(medium) = nsh
        WRITE(6,3850)medium,nsh
3850    FORMAT(' Medium ',i3,' has ',i3,' shells: ')                            
          DO 3861 i=1,nsh
          j = shell_array(i,medium)
          eno_array(i,medium) = eno_array(i,medium)/pztot
cale          WRITE(6,3870)i,j,shn_array(j),eno_array(i,medium), Jo_array(j)
cale     *    ,be_array(j)*rm*1000.
3870      FORMAT(i3,i4,i3,f9.5,e10.3,f10.3)
3861    CONTINUE
3862    CONTINUE
3791  CONTINUE
3792  CONTINUE
      WRITE(6,3880)
3880  FORMAT('...... Done.'/)                                                   
        DO 3891 j=1,numgeom
        medium = med(j)
        IF (( medium .GT. 0 .AND. medium .LE. nmed)) THEN
          IF (( iedgfl(j) .GT. 0 .AND. iedgfl(j) .LE. 100 )) THEN
            goto 3900
          END IF
        END IF
3891  CONTINUE
3892  CONTINUE
      WRITE(6,3910)
3910  FORMAT(' In subroutine init_compton: ',/ '   fluorescence not set         
     *but relaxation data are required for ',/ '   bound Compton scatter        
     *ing. ',/ '   calling EDGSET. '//)                                         
      iedgfl(1) = 1
      call edgset(1,1)
      iedgfl(1) = 0
3900  CONTINUE
      return
      end
      subroutine fix_brems
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c

      integer*4 medium_orig,i
      real*4 Zt,Zb,Zf,Zg,Zv,fmax1,fmax2,Zi,pi_orig,fc,xi,aux,XSIF,FCOULC
        DO 3691 medium_orig=1,nmed
        log_ap(medium_orig) = log(ap(medium_orig))
        Zt = 0
        Zb = 0
        Zf = 0
          DO 3701 i=1,NNE(medium_orig)
          Zi = ZELEM(medium_orig,i)
          pi_orig = PZ(medium_orig,i)
          fc = FCOULC(Zi)
          xi = XSIF(Zi)
          aux = pi_orig*Zi*(Zi + xi)
          Zt = Zt + aux
          Zb = Zb - aux*Log(Zi)/3
          Zf = Zf + aux*fc
3701    CONTINUE
3702    CONTINUE
        Zv = (Zb - Zf)/Zt
        Zg = Zb/Zt
        fmax1 = 2*(20.863 + 4*Zg) - 2*(20.029 + 4*Zg)/3
        fmax2 = 2*(20.863 + 4*Zv) - 2*(20.029 + 4*Zv)/3
        dl1(1,medium_orig) = (20.863 + 4*Zg)/fmax1
        dl2(1,medium_orig) = -3.242/fmax1
        dl3(1,medium_orig) = 0.625/fmax1
        dl4(1,medium_orig) = (21.12+4*Zg)/fmax1
        dl5(1,medium_orig) = -4.184/fmax1
        dl6(1,medium_orig) = 0.952
        dl1(2,medium_orig) = (20.029+4*Zg)/fmax1
        dl2(2,medium_orig) = -1.93/fmax1
        dl3(2,medium_orig) = -0.086/fmax1
        dl4(2,medium_orig) = (21.12+4*Zg)/fmax1
        dl5(2,medium_orig) = -4.184/fmax1
        dl6(2,medium_orig) = 0.952
        dl1(3,medium_orig) = (20.863 + 4*Zv)/fmax2
        dl2(3,medium_orig) = -3.242/fmax2
        dl3(3,medium_orig) = 0.625/fmax2
        dl4(3,medium_orig) = (21.12+4*Zv)/fmax2
        dl5(3,medium_orig) = -4.184/fmax2
        dl6(3,medium_orig) = 0.952
        dl1(4,medium_orig) = (20.029+4*Zv)/fmax2
        dl2(4,medium_orig) = -1.93/fmax2
        dl3(4,medium_orig) = -0.086/fmax2
        dl4(4,medium_orig) = (21.12+4*Zv)/fmax2
        dl5(4,medium_orig) = -4.184/fmax2
        dl6(4,medium_orig) = 0.952
        dl1(5,medium_orig) = (3*(20.863 + 4*Zg) - (20.029 + 4*Zg))
        dl2(5,medium_orig) = (3*(-3.242) - (-1.930))
        dl3(5,medium_orig) = (3*(0.625)-(-0.086))
        dl4(5,medium_orig) = (2*21.12+8*Zg)
        dl5(5,medium_orig) = (2*(-4.184))
        dl6(5,medium_orig) = 0.952
        dl1(6,medium_orig) = (3*(20.863 + 4*Zg) + (20.029 + 4*Zg))
        dl2(6,medium_orig) = (3*(-3.242) + (-1.930))
        dl3(6,medium_orig) = (3*0.625+(-0.086))
        dl4(6,medium_orig) = (4*21.12+16*Zg)
        dl5(6,medium_orig) = (4*(-4.184))
        dl6(6,medium_orig) = 0.952
        dl1(7,medium_orig) = (3*(20.863 + 4*Zv) - (20.029 + 4*Zv))
        dl2(7,medium_orig) = (3*(-3.242) - (-1.930))
        dl3(7,medium_orig) = (3*(0.625)-(-0.086))
        dl4(7,medium_orig) = (2*21.12+8*Zv)
        dl5(7,medium_orig) = (2*(-4.184))
        dl6(7,medium_orig) = 0.952
        dl1(8,medium_orig) = (3*(20.863 + 4*Zv) + (20.029 + 4*Zv))
        dl2(8,medium_orig) = (3*(-3.242) + (-1.930))
        dl3(8,medium_orig) = (3*0.625+(-0.086))
        dl4(8,medium_orig) = (4*21.12+16*Zv)
        dl5(8,medium_orig) = (4*(-4.184))
        dl6(8,medium_orig) = 0.952
        bpar(2,medium_orig) = dl1(7,medium_orig)/(3*dl1(8,medium_orig)
     *  + dl1(7,medium_orig))
        bpar(1,medium_orig) = 12*dl1(8,medium_orig)/(3*dl1(8, 
     *  medium_orig)+ dl1(7,medium_orig))
3691  CONTINUE
3692  CONTINUE
      return
      end
      real*4 function FCOULC(Z)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      real*4 Z
      real*4 fine,asq
      data fine/137.03604/
      asq = Z/fine
      asq = asq*asq
      FCOULC = asq*(1.0/(1.0+ASQ)+0.20206+ASQ*(-0.0369+ASQ*(0.0083+ASQ*(
     *-0.002))))
      return
      end
      real*4 function XSIF(Z)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      real*4 Z
      integer*4 iZ
      real*4 alrad(4),alradp(4),a1440,a183,FCOULC
      data alrad/5.31,4.79,4.74,4.71/
      data alradp/6.144,5.621,5.805,5.924/
      data a1440/1194.0/,A183/184.15/
      IF (( Z .LE. 4 )) THEN
        iZ = Z
        xsif = alradp(iZ)/(alrad(iZ) - FCOULC(Z))
      ELSE
        xsif = Log(A1440*Z**(-0.666667))/(Log(A183*Z**(-0.33333))-FCOULC
     *  (Z))
      END IF
      return
      end
      subroutine init_nist_brems
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none


      real*4 energy_array(57),x_array(30), cs_array(57,30,100)
      real*4 xi_array(30)
      real*8 x_gauss(64),w_gauss(64)
      integer*4 nmix,kmix,i,n,k,j,ii
      integer*4 ngauss,i_gauss
      integer*4 ifirst,ilast,nener,neke,leil
      real*4 cs(57,30),ee(57),ele(57)
      real*4 csx(30),afx(30),bfx(30),cfx(30),dfx(30)
      real*4 cse(57),afe(57),bfe(57),cfe(57),dfe(57)
      real*4 Z_orig,sumA
      real*4 emin,xi,res,spline,eil,ei,beta2,aux,sigb,sigt,ebr1,ebr2
      real*4 sigee,sigep,sige,si_esig,si1_esig,si_ebr1,si1_ebr1,ededx, s
     *ig_bhabha,si_psig,si1_psig,si_pbr1,si1_pbr1,si_pbr2,si1_pbr2
      integer*4 iz
      real*4 ple,qle,x_orig,f,error,max_error,x_max_error,f_max_error
      integer*4 ndat,k_max_error
      real*4 amu
      parameter (amu = 1660.5655)
cale
      include 'egs4comm.for'
      include 'egs4fcomm.for'

      character*80 bremfilename

      character*80 totbremfilename

      bremfilename='nist_brems.dat'

      totbremfilename=pathname(1:pathlength) // bremfilename
cale     write(*,*)'cale:',totbremfilename
      open(unit=76,file=totbremfilename,form='formatted',status='old')
cale
      read(76,*)
      read(76,*) nmix,kmix
      IF (( kmix .NE. 30 .OR. nmix .NE. 57)) THEN
        write(6,*) ' init_nist_brems: wrong data file!'                         
        stop
      END IF
      read(76,*) energy_array
      read(76,*) x_array
      read(76,*)
        DO 5081 i=1,100
        read(76,*) ((cs_array(n,k,i),n=1,nmix),k=1,kmix)
5081  CONTINUE
5082  CONTINUE
      close(76)
        DO 5091 k=1,kmix
        xi_array(k)=Log(1-x_array(k)+1e-6)
5091  CONTINUE
5092  CONTINUE
      ngauss = 64
      call gauss_legendre(0d0,1d0,x_gauss,w_gauss,ngauss)
      write(6,*)
      write(6,*) ' Using NIST brems cross sections! ' 
      write(6,*) ' Initializing brems data for media '                          
      write(6,*)
        DO 5101 medium=1,nmed
        log_ap(medium) = log(ap(medium))
cale        write(6,*) ' Initializing brems data for medium ',medium,'...'          
        emin = max(ae(medium) - rm, ap(medium))
          DO 5111 i=1,nmix
          IF((energy_array(i) .GE. emin))GO TO5112
5111    CONTINUE
5112    CONTINUE
        ifirst = i
          DO 5121 i=nmix,1,-1
          IF((energy_array(i) .LT. ue(medium) - rm))GO TO5122
5121    CONTINUE
5122    CONTINUE
        ilast = i+1
        IF (( ifirst .LT. 1 .OR. ilast .GT. nmix )) THEN
          write(6,*) ' init_nist_brems: data available only for '               
          write(6,*) energy_array(i),' <= E <= ',energy_array(nmix)             
          write(6,*) ' will use spline interpolations to get cross '            
          write(6,*) ' sections beyond the available data but this may'         
          write(6,*) ' produce nonsense!'                                       
          IF((ifirst .LT. 1))ifirst=1
          IF((ilast .GT. nmix))ilast = nmix
        END IF
          DO 5131 i=ifirst,ilast
          ii = i+1 - ifirst
          ee(ii) = energy_array(i)
          ele(ii) = log(ee(ii))
          sumA = 0
            DO 5141 j=1,NNE(medium)
            sumA = sumA + pz(medium,j)*wa(medium,j)
5141      CONTINUE
5142      CONTINUE
          sumA = sumA*amu
            DO 5151 k=1,kmix
            cs(ii,k) = 0
              DO 5161 j=1,NNE(medium)
              Z_orig = zelem(medium,j)
              iz = int(Z_orig+0.1)
              Z_orig = Z_orig*Z_orig/sumA
              cs(ii,k) = cs(ii,k) + pz(medium,j)*Z_orig*cs_array(i,k,iz)
5161        CONTINUE
5162        CONTINUE
            csx(k) = Log(cs(ii,k))
5151      CONTINUE
5152      CONTINUE
          call set_spline(xi_array,csx,afx,bfx,cfx,dfx,kmix)
          cse(ii) = 0
          aux = Log(ee(ii)/ap(medium))
            DO 5171 i_gauss=1,ngauss
            xi = log(1 - ap(medium)/ee(ii)*exp(x_gauss(i_gauss)*aux)+1e-
     *      6)
            res = spline(xi,xi_array,afx,bfx,cfx,dfx,kmix)
            cse(ii) = cse(ii) + w_gauss(i_gauss)*exp(res)
5171      CONTINUE
5172      CONTINUE
5131    CONTINUE
5132    CONTINUE
        nener = ilast - ifirst + 1
        call set_spline(ele,cse,afe,bfe,cfe,dfe,nener)
        neke = meke(medium)
        sigee = 1E-15
        sigep = 1E-15
          DO 5181 i=1,neke
          eil = (float(i) - eke0(medium))/eke1(medium)
          ei = exp(eil)
          leil = i
          beta2 = ei*(ei+2*rm)/(ei+rm)**2
          IF (( ei .LE. ap(medium) )) THEN
            sigb = 1e-30
          ELSE
            sigb = spline(eil,ele,afe,bfe,cfe,dfe,nener)
            sigb = sigb*log(ei/ap(medium))/beta2*rho(medium)
          END IF
          sigt=esig1(Leil,MEDIUM)*eil+esig0(Leil,MEDIUM)
          ebr1=ebr11(Leil,MEDIUM)*eil+ebr10(Leil,MEDIUM)
          IF((sigt .LT. 0))sigt = 0
          IF((ebr1 .GT. 1))ebr1 = 1
          IF((ebr1 .LT. 0))ebr1 = 0
          IF (( i .GT. 1 )) THEN
            si_esig = si1_esig
            si_ebr1 = si1_ebr1
            si1_esig = sigt*(1 - ebr1) + sigb
            si1_ebr1 = sigb/si1_esig
            esig1(i-1,medium) = (si1_esig - si_esig)*eke1(medium)
            esig0(i-1,medium) = si1_esig - esig1(i-1,medium)*eil
            ebr11(i-1,medium) = (si1_ebr1 - si_ebr1)*eke1(medium)
            ebr10(i-1,medium) = si1_ebr1 - ebr11(i-1,medium)*eil
          ELSE
            si1_esig = sigt*(1 - ebr1) + sigb
            si1_ebr1 = sigb/si1_esig
          END IF
          sigt=psig1(Leil,MEDIUM)*eil+psig0(Leil,MEDIUM)
          ebr1=pbr11(Leil,MEDIUM)*eil+pbr10(Leil,MEDIUM)
          ebr2=pbr21(Leil,MEDIUM)*eil+pbr20(Leil,MEDIUM)
          IF((sigt .LT. 0))sigt = 0
          IF((ebr1 .GT. 1))ebr1 = 1
          IF((ebr1 .LT. 0))ebr1 = 0
          IF((ebr2 .GT. 1))ebr2 = 1
          IF((ebr2 .LT. 0))ebr2 = 0
          sig_bhabha = sigt*(ebr2 - ebr1)
          IF((sig_bhabha .LT. 0))sig_bhabha = 0
          IF (( i .GT. 1 )) THEN
            si_psig = si1_psig
            si_pbr1 = si1_pbr1
            si_pbr2 = si1_pbr2
            si1_psig = sigt*(1 - ebr1) + sigb
            si1_pbr1 = sigb/si1_psig
            si1_pbr2 = (sigb + sig_bhabha)/si1_psig
            psig1(i-1,medium) = (si1_psig - si_psig)*eke1(medium)
            psig0(i-1,medium) = si1_psig - psig1(i-1,medium)*eil
            pbr11(i-1,medium) = (si1_pbr1 - si_pbr1)*eke1(medium)
            pbr10(i-1,medium) = si1_pbr1 - pbr11(i-1,medium)*eil
            pbr21(i-1,medium) = (si1_pbr2 - si_pbr2)*eke1(medium)
            pbr20(i-1,medium) = si1_pbr2 - pbr20(i-1,medium)*eil
          ELSE
            si1_psig = sigt*(1 - ebr1) + sigb
            si1_pbr1 = sigb/si1_psig
            si1_pbr2 = (sigb + sig_bhabha)/si1_psig
          END IF
          ededx=ededx1(Leil,MEDIUM)*eil+ededx0(Leil,MEDIUM)
          sige = si1_esig/ededx
          IF((sige .GT. sigee))sigee = sige
          ededx=pdedx1(Leil,MEDIUM)*eil+pdedx0(Leil,MEDIUM)
          sige = si1_psig/ededx
          IF((sige .GT. sigep))sigep = sige
5181    CONTINUE
5182    CONTINUE
        esig1(neke,medium) = esig1(neke-1,medium)
        esig0(neke,medium) = esig0(neke-1,medium)
        ebr11(neke,medium) = ebr11(neke-1,medium)
        ebr10(neke,medium) = ebr10(neke-1,medium)
        psig1(neke,medium) = psig1(neke-1,medium)
        psig0(neke,medium) = psig0(neke-1,medium)
        pbr11(neke,medium) = pbr11(neke-1,medium)
        pbr10(neke,medium) = pbr10(neke-1,medium)
        pbr21(neke,medium) = pbr21(neke-1,medium)
        pbr20(neke,medium) = pbr20(neke-1,medium)
cale        write(6,*) ' Max. new cross sections per energy loss: ',sigee,si        
cale     *  gep
        esig_e(medium) = sigee
        psig_e(medium) = sigep
        IF((sigee .GT. esige_max))esige_max = sigee
        IF((sigep .GT. psige_max))psige_max = sigep
        nb_emin(medium) = energy_array(ifirst)
        IF (( nb_emin(medium) .LE. ap(medium) )) THEN
          nb_emin(medium) = energy_array(ifirst+1)
        END IF
        nb_emax(medium) = energy_array(ilast)
        nb_lemin(medium) = log(nb_emin(medium))
        nb_lemax(medium) = log(nb_emax(medium))
        nb_dle(medium) = (nb_lemax(medium) - nb_lemin(medium))/(100-1)
        nb_dlei(medium) = 1/nb_dle(medium)
        eil = nb_lemin(medium) - nb_dle(medium)
          DO 5191 i=1,100
          eil = eil + nb_dle(medium)
          ei = exp(eil)
            DO 5201 ii=1,nener
            IF((ei .LT. ee(ii)))GO TO5202
5201      CONTINUE
5202      CONTINUE
          ii = ii-1
          IF((ii .GT. nener-1))ii = nener-1
          ple = (eil - ele(ii))/(ele(ii+1)-ele(ii))
          qle = 1 - ple
            DO 5211 k=1,30
            csx(k) = log(qle*cs(ii,k) + ple*cs(ii+1,k))
5211      CONTINUE
5212      CONTINUE
          call set_spline(xi_array,csx,afx,bfx,cfx,dfx,kmix)
          x_orig = ap(medium)/ei
          aux = -log(x_orig)
          xi = log(1 - x_orig+1e-6)
          res = spline(xi,xi_array,afx,bfx,cfx,dfx,kmix)
          nb_xdata(0,i,medium) = 0
          nb_fdata(0,i,medium) = exp(res)
            DO 5221 k=1,30
            IF((x_array(k) .GT. x_orig))GO TO5222
5221      CONTINUE
5222      CONTINUE
          IF((k .GT. 30))k = 30
          ndat = 0
            DO 5231 j=k+1,30-1
            ndat = ndat+1
            nb_xdata(ndat,i,medium) = log(x_array(j)/x_orig)/aux
            nb_fdata(ndat,i,medium) = exp(csx(j))
5231      CONTINUE
5232      CONTINUE
          ndat = ndat+1
          nb_xdata(ndat,i,medium) = 1
          nb_fdata(ndat,i,medium) = exp(csx(30))
          IF((ndat .EQ. 50))goto 5240
5251      CONTINUE
            x_max_error = 0
            f_max_error = 0
            k_max_error = 0
            max_error = 0
              DO 5261 k=0,ndat-1
              x_orig = 0.5*(nb_xdata(k,i,medium) +
     a         nb_xdata(k+1,i,medium))
              f = 0.5*(nb_fdata(k,i,medium) + nb_fdata(k+1,i,medium))
              xi = log(1 - ap(medium)/ei*exp(x_orig*aux)+1e-6)
              res = spline(xi,xi_array,afx,bfx,cfx,dfx,kmix)
              res = exp(res)
              error = abs(1-f/res)
              IF (( error .GT. max_error )) THEN
                x_max_error = x_orig
                f_max_error = res
                max_error = error
                k_max_error = k
              END IF
5261        CONTINUE
5262        CONTINUE
            ndat = ndat+1
              DO 5271 k=ndat,k_max_error+2,-1
              nb_xdata(k,i,medium) = nb_xdata(k-1,i,medium)
              nb_fdata(k,i,medium) = nb_fdata(k-1,i,medium)
5271        CONTINUE
5272        CONTINUE
            nb_xdata(k_max_error+1,i,medium) = x_max_error
            nb_fdata(k_max_error+1,i,medium) = f_max_error
            IF(((ndat .EQ. 50)))GO TO5252
          GO TO 5251
5252      CONTINUE
5240      call prepare_alias_table(50,nb_xdata(0,i,medium), nb_fdata(0,i
     *    ,medium),nb_wdata(1,i,medium),nb_idata(1,i,medium))
5191    CONTINUE
5192    CONTINUE
5101  CONTINUE
5102  CONTINUE
      write(6,*)
      write(6,*)
      return
      end
      subroutine gauss_legendre(x1,x2,x,w,n)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      integer*4 n
      real*8 x1,x2,x(n),w(n)
      real*8 eps,Pi
      parameter (eps = 3.D-14,Pi=3.141592654D0)
      integer*4 i,m,j
      real*8 xm,xl,z,z1,p1,p2,p3,pp
      m = (n + 1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
        DO 5331 i=1,m
        z=cos(Pi*(i-.25d0)/(n+.5d0))
5341    CONTINUE
          p1=1.d0
          p2=0.d0
            DO 5351 j=1,n
            p3 = p2
            p2 = p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
5351      CONTINUE
5352      CONTINUE
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
          IF(((abs(z-z1) .LT. eps)))GO TO5342
        GO TO 5341
5342    CONTINUE
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
5331  CONTINUE
5332  CONTINUE
      return
      end
      subroutine set_spline(x,f,a,b,c,d,n)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      integer*4 n
      real*4 x(n),f(n),a(n),b(n),c(n),d(n)
      integer*4 m1,m2,m,mr
      real*4 s,r
      m1 = 2
      m2 = n-1
      s = 0
        DO 4411 m=1,m2
        d(m) = x(m+1) - x(m)
        r = (f(m+1) - f(m))/d(m)
        c(m) = r - s
        s = r
4411  CONTINUE
4412  CONTINUE
      s=0
      r=0
      c(1)=0
      c(n)=0
        DO 4421 m=m1,m2
        c(m) = c(m) + r*c(m-1)
        b(m) = 2*(x(m-1) - x(m+1)) - r*s
        s = d(m)
        r = s/b(m)
4421  CONTINUE
4422  CONTINUE
      mr = m2
        DO 4431 m=m1,m2
        c(mr) = (d(mr)*c(mr+1) - c(mr))/b(mr)
        mr = mr - 1
4431  CONTINUE
4432  CONTINUE
        DO 4441 m=1,m2
        s = d(m)
        r = c(m+1) - c(m)
        d(m) = r/s
        c(m) = 3*c(m)
        b(m) = (f(m+1)-f(m))/s - (c(m)+r)*s
        a(m) = f(m)
4441  CONTINUE
4442  CONTINUE
      return
      end
      real*4 function spline(s,x,a,b,c,d,n)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      integer*4 n
      real*4 s,x(n),a(n),b(n),c(n),d(n)
      integer*4 m_lower,m_upper,direction,m,ml,mu,mav
      real*4 q
      IF (( x(1) .GT. x(n) )) THEN
        direction = 1
        m_lower = n
        m_upper = 0
      ELSE
        direction = 0
        m_lower = 0
        m_upper = n
      END IF
      IF (( s .GE. x(m_upper + direction) )) THEN
        m = m_upper + 2*direction - 1
      ELSE IF(( s .LE. x(m_lower+1-direction) )) THEN
        m = m_lower - 2*direction + 1
      ELSE
        ml = m_lower
        mu = m_upper
4451    IF(iabs(mu-ml).LE.1)GO TO 4452
          mav = (ml+mu)/2
          IF (( s .LT. x(mav) )) THEN
            mu = mav
          ELSE
            ml = mav
          END IF
        GO TO 4451
4452    CONTINUE
        m = mu + direction - 1
      END IF
      q = s - x(m)
      spline = a(m) + q*(b(m) + q*(c(m) + q*d(m)))
      return
      end
      subroutine prepare_alias_table(nsbin,xs_array,fs_array,ws_array,ib
     *in_array)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c

      integer*4 nsbin,ibin_array(nsbin)
      real*4 xs_array(0:nsbin),fs_array(0:nsbin),ws_array(nsbin)
      integer*4 i,j_l,j_h
      real*4 sum,aux
      sum = 0
        DO 5281 i=1,nsbin
        aux = 0.5*(fs_array(i)+fs_array(i-1))*(xs_array(i)-xs_array(i-1)
     *  )
        IF((aux .LT. 1e-30))aux = 1e-30
        ws_array(i) = -aux
        ibin_array(i) = 1
        sum = sum + aux
5281  CONTINUE
5282  CONTINUE
      sum = sum/nsbin
        DO 5291 i=1,nsbin-1
          DO 5301 j_h=1,nsbin
          IF (( ws_array(j_h) .LT. 0 )) THEN
            IF((abs(ws_array(j_h)) .GT. sum))GOTO 5310
          END IF
5301    CONTINUE
5302    CONTINUE
        j_h = nsbin
5310    CONTINUE
          DO 5311 j_l=1,nsbin
          IF (( ws_array(j_l) .LT. 0 )) THEN
            IF((abs(ws_array(j_l)) .LT. sum))GOTO 5320
          END IF
5311    CONTINUE
5312    CONTINUE
        j_l = nsbin
5320    aux = sum - abs(ws_array(j_l))
        ws_array(j_h) = ws_array(j_h) + aux
        ws_array(j_l) = -ws_array(j_l)/sum
        ibin_array(j_l) = j_h
        IF((i .EQ. nsbin-1))ws_array(j_h) = 1
5291  CONTINUE
5292  CONTINUE
      return
      end
      real*4 function alias_sample1(nsbin,xs_array,fs_array,ws_array,ibi
     *n_array)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c

      integer*4 nsbin,ibin_array(nsbin)
      real*4 xs_array(0:nsbin),fs_array(0:nsbin),ws_array(nsbin)
      integer*4 j
      real*4 r1,r2,aj,x_orig,dx,a,rnno1
      IF (( rng_seed .GT. 24 )) THEN
        call ranlux(rng_array)
        rng_seed = 1
      END IF
      r1 = rng_array(rng_seed)
      rng_seed = rng_seed + 1
      IF (( rng_seed .GT. 24 )) THEN
        call ranlux(rng_array)
        rng_seed = 1
      END IF
      r2 = rng_array(rng_seed)
      rng_seed = rng_seed + 1
      aj = 1 + r1*nsbin
      j = aj
      aj = aj - j
      IF((aj .GT. ws_array(j)))j = ibin_array(j)
      x_orig = xs_array(j-1)
      dx = xs_array(j)-x_orig
      IF (( fs_array(j-1) .GT. 0 )) THEN
        a = fs_array(j)/fs_array(j-1)-1
        IF (( abs(a) .LT. 0.2 )) THEN
          rnno1 = 0.5*(1-r2)*a
          alias_sample1 = x_orig + r2*dx*(1+rnno1*(1-r2*a))
        ELSE
          alias_sample1 = x_orig - dx/a*(1-sqrt(1+r2*a*(2+a)))
        END IF
      ELSE
        alias_sample1 = x_orig + dx*sqrt(r2)
      END IF
      return
      end
      subroutine vmc_electron(ircode)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      integer*4 ircode
      WRITE(6,5360)
      WRITE(1,5360)
5360  FORMAT(//' ********* VMC Transport option not in this distribution        
     * ****** '//)                                                              
      stop
      end
      FUNCTION lnblnk1(c)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      CHARACTER c*(*)
      INTEGER j
        DO 5371 j=LEN(c),1,-1
        IF ((c(j:j) .NE. ' ')) THEN                                             
          lnblnk1=j
          RETURN
        END IF
5371  CONTINUE
5372  CONTINUE
      lnblnk1=0
      RETURN
      end
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      subroutine evinit
c
      implicit none
      include 'egs4comm.for'
c
      integer mh, nt
c
      first_int_label = .true.
c
      nhit         = 0
      elost        = 0.0
      biosource    = 0.0
      multiplicity = 0       
c
c---  set o zero energy loss in all regions
c
cnico nico cambia nplan con regoffset
c      do 10 mh=1,nplan
      do 10 mh=1,regoffset
            deodx(mh)=0.d0
   10 continue
c
      do 20 nt = 1,size_of_ntup
            xtup3(nt) =  0.0
   20 continue
c
      eavg=0.0
      xavg=0.0
      yavg=0.0
      zavg=0.0
c
      ecountmax=0.0
      xcountmax=0.0
      ycountmax=0.0
      zcountmax=0.0
c
      efirst=0.0
      xfirst=0.0
      yfirst=0.0
      zfirst=0.0
c
      epixel=0.0
      xpixel=0.0
      ypixel=0.0
      zpixel=0.0
c
      epixmax=0.0
      xpixmax=0.0
      ypixmax=0.0
      zpixmax=0.0
c
      return
c
      end
c
      subroutine source
c
c-----------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
c
c-----------------------------------------------------------c
c      sets the kinetic and total energy for monoenergetic
c      gamma  decays
c
      call ranlux(rng_array)
      rng_seed = 1
c
cnico      evkin = enkin
      evkin = esource(nactsx) + abs(float(ichar)) * rm
c
cnico      ein    = enkin
      ein    = evkin
c
      eireal = ein
c
      call hfill (100,eireal,0.,1.)
c
      totein  =  totein + ein
c
c------------- if different geometrical sources
c
c------------- isxgeomtype:
c------------- 1:elliptic source
c------------- 2:cilindrical source
c------------- 3:box source
c
c
      if(isxgeomtype(nactsx).eq.1) then
c
c------------- elliptic source
c	     
   11     continue
          iloop = iloop + 1
c
          xi0 = sx_sup(nactsx)*rng_array(rng_seed)+sx_inf(nactsx)
          rng_seed = rng_seed + 1
          yi0 = sy_sup(nactsx)*rng_array(rng_seed)+sy_inf(nactsx)
          rng_seed = rng_seed + 1
          zi0 = sz_sup(nactsx)*rng_array(rng_seed)+sz_inf(nactsx)
          rng_seed = rng_seed + 1
c
c------------- check if the point(xi0,yi0,zi0) is inside the sphere
c
          if ( (((xi0-sx_cent(nactsx))*(xi0-sx_cent(nactsx)))
     +	                            * rsxsqinv(nactsx)
     +         + ((yi0-sy_cent(nactsx))*(yi0-sy_cent(nactsx)))
     +	                            * rsysqinv(nactsx)
     +         + ((zi0-sz_cent(nactsx))*(zi0-sz_cent(nactsx)))
     +	                           * rszsqinv(nactsx)) - 1 )    21,21,99
c
 99       IF (( rng_seed .GT. 20 )) THEN
               call ranlux(rng_array)
               rng_seed = 1
          END IF
          goto 11
c
c------------------------------------------------
c
   21     nevent = nevent + 1
c
      elseif(isxgeomtype(nactsx).eq.2) then
c
c------------- cilindrical source
c
   12     continue
          iloop = iloop + 1
c
          xi0 = sx_sup(nactsx)*rng_array(rng_seed)+sx_inf(nactsx)
          rng_seed = rng_seed + 1
          yi0 = sy_sup(nactsx)*rng_array(rng_seed)+sy_inf(nactsx)
          rng_seed = rng_seed + 1
          zi0 = sz_sup(nactsx)*rng_array(rng_seed)+sz_inf(nactsx)
          rng_seed = rng_seed + 1
c
c------------- check if the point(xi0,yi0,zi0) is inside the disk
c
          if ( (((xi0-sx_cent(nactsx))*(xi0-sx_cent(nactsx)))
     +	                            * rsxsqinv(nactsx)
     +         + ((yi0-sy_cent(nactsx))*(yi0-sy_cent(nactsx))) 
     +	                            * rsysqinv(nactsx)) - 1 )    22,22,98
c
 98       IF (( rng_seed .GT. 20 )) THEN
               call ranlux(rng_array)
               rng_seed = 1
          END IF
          goto 12
c
c------------------------------------------------
c
   22     nevent = nevent + 1
c
c ale aggiunge per avere anche sorgenti //asse x. 19 04 2001
c 
      elseif(isxgeomtype(nactsx).eq.3) then
c
c------------- cilindrical source // x axis
c
  212     continue
          iloop = iloop + 1
c
          xi0 = sx_sup(nactsx)*rng_array(rng_seed)+sx_inf(nactsx)
          rng_seed = rng_seed + 1
          yi0 = sy_sup(nactsx)*rng_array(rng_seed)+sy_inf(nactsx)
          rng_seed = rng_seed + 1
          zi0 = sz_sup(nactsx)*rng_array(rng_seed)+sz_inf(nactsx)
          rng_seed = rng_seed + 1
c
c------------- check if the projected point(0,yi0,zi0) is inside the disk
c
          if ( (((zi0-sz_cent(nactsx))*(zi0-sz_cent(nactsx)))
     +	                         * rszsqinv(nactsx)
     +         + ((yi0-sy_cent(nactsx))*(yi0-sy_cent(nactsx)))
     +	                         * rsysqinv(nactsx)) - 1 )    222,222,97
c
 97       IF (( rng_seed .GT. 20 )) THEN
               call ranlux(rng_array)
               rng_seed = 1
          END IF
          goto 212	 
c
c------------------------------------------------
c
  222     nevent = nevent + 1
c
c ale cambia box source da 3 a 4. 19 04 2001
c
      elseif(isxgeomtype(nactsx).eq.4) then	
c
c------------- box source
c
          iloop = iloop + 1
c
          xi0 = sx_sup(nactsx)*rng_array(rng_seed)+sx_inf(nactsx)
          rng_seed = rng_seed + 1
          yi0 = sy_sup(nactsx)*rng_array(rng_seed)+sy_inf(nactsx)
          rng_seed = rng_seed + 1
          zi0 = sz_sup(nactsx)*rng_array(rng_seed)+sz_inf(nactsx)
          rng_seed = rng_seed + 1
c
c------------------------------------------------
c
          nevent = nevent + 1
c
      endif
c
c----------
c
      xin      =   xi0
      yin      =   yi0
      zin      =   zi0
c
      call hfill (101,xin,0.,1.)
      call hfill (102,yin,0.,1.)
      call hfill (103,zin,0.,1.)
      call hfill (200,xin,yin,1.)
      call hfill (201,xin,zin,1.)
      call hfill (202,yin,zin,1.)
cc
c      call hfill (300,xi,yi,zi)
cc
c------------- parallel-beam source
c
      if(thetaangsx(nactsx).eq.0.) then
c
          cth = 1.
          sth = 0.
          phibis = 0.
          goto 31
c
      endif
c
c------------- if different angular sources
c
c------------- isxangtype:
c------------- 1:cone-beam source (isotropic)
c------------- 2:fan-beam source (isotropic on a plane)
c
c
      IF (( rng_seed .GT. 18 )) THEN
              call ranlux(rng_array)
              rng_seed = 5
      END IF
c
      if(isxangtype(nactsx).eq.1) then
c
c------------- cone-beam source (isotropic)
c
          cth = cthinf(nactsx) - cthsup(nactsx) * rng_array(rng_seed)
          rng_seed = rng_seed + 1
          sth = sqrt ( 1.- cth*cth )
c
          phibis = - pi +  twopi * rng_array(rng_seed)
          rng_seed = rng_seed + 1
c
      elseif(isxangtype(nactsx).eq.2) then
c
c------------- fan-beam source (isotropic on a plane)
c
          thetarandom=rng_array(rng_seed)*thetaangsx(nactsx)
          rng_seed = rng_seed + 1
          cth = cos(thetarandom)
          sth = sin(thetarandom)
c
          if (rng_array(rng_seed).ge.0.5) then
                phibis = phiangsx(nactsx)
          else
                phibis = phiangsx(nactsx) + pi
          endif
          rng_seed = rng_seed + 1
c
      endif
c
c-----------------------------------------------------------
c
 31   win = cth
      uin = sth * cos (phibis)
      vin = sth * sin (phibis)             
c
      dnorm = sqrt(uin*uin+vin*vin+win*win)

      win = win / dnorm
      uin = uin / dnorm
      vin = vin / dnorm
c
      call hfill (110,uin,0.,1.)
      call hfill (111,vin,0.,1.)
      call hfill (112,win,0.,1.)
c
c-----------------------------------------------------------
c
      irin =  kgeom(xin,yin,zin)
      wtin =  1.00000
c
      return
      end
c
c-----------------------------------------------------
c
      SUBROUTINE SHOWER(IQI,EI,XI,YI,ZI,UI,VI,WI,IRI,WTI)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      real*4 EI,  XI,YI,ZI, UI,VI,WI, WTI
      integer*4 IQI,  IRI
      DOUBLE PRECISION DEG,  DPGL,  DEI,  DPI,  DCSTH,  DCOSTH,  PI0MSQ
      real*4 DNEARI,  CSTH
      integer*4 IRCODE
      DATA PI0MSQ/1.8215416D4/
      NP=1
      NPold = NP
      DNEARI=0.0
      IQ(1)=IQI
      E(1)=EI
      U(1)=UI
      V(1)=VI
      W(1)=WI
      X(1)=XI
      Y(1)=YI
      Z(1)=ZI
      IR(1)=IRI
      WT(1)=WTI
      DNEAR(1)=DNEARI
      LATCH(1)=LATCHI
      IF ((IQI.EQ.2)) THEN
        IF ((EI**2.LE.PI0MSQ)) THEN
          WRITE(6,4930)EI
4930      FORMAT(//,' STOPPED IN SUBROUTINE SHOWER---PI-ZERO OPTION INVO        
     *KED', /,' BUT THE TOTAL ENERGY WAS TOO SMALL (EI=',G15.5,' MEV)')         
          STOP
        END IF
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        CSTH = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        DCSTH=CSTH
        DEI=EI
        DPI=DSQRT(DEI*DEI-PI0MSQ)
        DEG=DEI+DPI*DCSTH
        DPGL=DPI+DEI*DCSTH
        DCOSTH=DPGL/DEG
        COSTHE=DCOSTH
        SINTHE=DSQRT(1.D0-DCOSTH*DCOSTH)
        IQ(1)=0
        E(1)=DEG/2.
        CALL UPHI(2,1)
        NP=2
        DEG=DEI-DPI*DCSTH
        DPGL=DPI-DEI*DCSTH
        DCOSTH=DPGL/DEG
        COSTHE=DCOSTH
        SINTHE=-DSQRT(1.D0-DCOSTH*DCOSTH)
        IQ(2)=0
        E(2)=DEG/2.
        CALL UPHI(3,2)
      END IF
4940  CONTINUE
4941    CONTINUE
        IF((IQ(NP).EQ.0))GO TO 4950
4961    CONTINUE
4970      CALL ELECTR(IRCODE)
          IF((IRCODE.EQ.2))GO TO4962
4950      CALL PHOTON(IRCODE)
          IF((IRCODE.EQ.2))GO TO4962
        GO TO 4961
4962    CONTINUE
        IF((NP.LE.0))GO TO4942
      GO TO 4941
4942  CONTINUE
      RETURN
      END
      SUBROUTINE UPHI(IENTRY,LVL)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      integer*4 IENTRY,LVL
      real*4 CTHET,  RNNO38,  PHI,  CPHI,  A,B,C,  SINPS2,  SINPSI,  US,
     *VS,  SINDEL,COSDEL
      integer*4 IARG,  LPHI,LTHETA,LCTHET,LCPHI
      real*4 xphi,xphi2,yphi,yphi2,rhophi2
      double precision norm
      IARG=21
      IF ((IAUSFL(IARG+1).NE.0)) THEN
        CALL AUSGAB(IARG)
      END IF
      GO TO (4980,4990,5000),IENTRY
      GO TO 5010
4980  CONTINUE
      SINTHE=sin(THETA)
      CTHET=PI5D2-THETA
      COSTHE=sin(CTHET)
4990  CONTINUE
5021  CONTINUE
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
        IF(rhophi2.LE.1)GO TO5022
      GO TO 5021
5022  CONTINUE
      rhophi2 = 1/rhophi2
      cosphi = (xphi2 - yphi2)*rhophi2
      sinphi = 2*xphi*yphi*rhophi2
5000  GO TO (5030,5040,5050),LVL
      GO TO 5010
5030  A=U(NP)
      B=V(NP)
      C=W(NP)
      GO TO 5060
5050  A=U(NP-1)
      B=V(NP-1)
      C=W(NP-1)
5040  X(NP)=X(NP-1)
      Y(NP)=Y(NP-1)
      Z(NP)=Z(NP-1)
      IR(NP)=IR(NP-1)
      WT(NP)=WT(NP-1)
      DNEAR(NP)=DNEAR(NP-1)
      LATCH(NP)=LATCH(NP-1)
5060  SINPS2=A*A+B*B
      IF ((SINPS2.LT.1.0E-20)) THEN
        U(NP)=SINTHE*COSPHI
        V(NP)=SINTHE*SINPHI
c
c known bug corrected - mirko
c
c        W(NP)=COSTHE
c
        W(NP)=C*COSTHE
cnico normalize u,v,w
	norm = u(np)*u(np)+v(np)*v(np)+w(np)*w(np)
	norm = 1./sqrt(norm)
	u(np) = u(np) * norm
	v(np) = v(np) * norm
	w(np) = w(np) * norm
cnico
c
c-------------------------
c
      ELSE
        SINPSI=SQRT(SINPS2)
        US=SINTHE*COSPHI
        VS=SINTHE*SINPHI
        SINDEL=B/SINPSI
        COSDEL=A/SINPSI
        U(NP)=C*COSDEL*US-SINDEL*VS+A*COSTHE
        V(NP)=C*SINDEL*US+COSDEL*VS+B*COSTHE
        W(NP)=-SINPSI*US+C*COSTHE
cnico normalize u,v,w
	norm = u(np)*u(np)+v(np)*v(np)+w(np)*w(np)
	norm = 1./sqrt(norm)
	u(np) = u(np) * norm
	v(np) = v(np) * norm
	w(np) = w(np) * norm
cnico
      END IF
      IARG=22
      IF ((IAUSFL(IARG+1).NE.0)) THEN
        CALL AUSGAB(IARG)
      END IF
      RETURN
5010  WRITE(6,5070)IENTRY,LVL
5070  FORMAT(' STOPPED IN UPHI WITH IENTRY,LVL=',2I6)                           
      STOP
      END
c
      subroutine ausgab ( iarg )
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
c
c
      integer iarg
      real edepr
c          
c-----------------------------------------------------------c
c----------    tracking section ----------------------------c
c
c  keep a running sum of the energy deposited in each region
c
c--------------statistics of regions -----------------------c
c
c mirko testing
c
c       if ((abs(w(np)).gt.1.0  ).or.(abs(u(np)).gt.1.0  ).or.
c     a     (abs(v(np)).gt.1.0  ).or.(abs(x(np)).gt.1000.).or.
c     b     (abs(y(np)).gt.1000.).or.(abs(z(np)).gt.1000.))  then
c       write(*,*) 'ausgab  ', x(np), y(np), z(np), u(np), v(np), w(np)
c     a, np, iarg, edep, idisc, iq(np)
c
c       endif
c
c      if ((iarg.eq.0).and.(edep.ne.0.d0)) write(*,*) 'iarg = 0', edep
c
c end - mirko testing
c
      if (edep.eq.0.d0) goto 10
c
      if (idisc.eq.1000) then 
c
          toteoob = toteoob + edep
          goto 10
c
      endif
c
cnico nico cambia nplan con regoffset
c      if(ir(np).le.nplan) then
      if(ir(np).le.regoffset) then
c
        deodx(ir(np)) = deodx(ir(np)) + edep
        totedeodx = totedeodx + edep
c
      else
c
        totemod = totemod + edep
c
      endif
c
c
      edepr = edep
      if (edepr.ne.0.) call hfill (120,edepr,0.,1.)
c
      if ((index(ir(np))/10000).eq.labeldet) then
c
          if (first_int_label.and.(edep.gt.0.d0)) then
c
             efirst = edep
             xfirst = x(np)
             yfirst = y(np)
             zfirst = z(np)
             first_int_label = .false.
c
          endif
c
          edepr=edep
          if (edepr.ne.0.) call hfill (121,edepr,0.,1.)
c
          totedet = totedet + edep
          elost     =  elost  +   edep
c
          eavg      =  eavg +  edep
          xavg      =  xavg +  edep * x(np)
          yavg      =  yavg +  edep * y(np)
          zavg      =  zavg +  edep * z(np)
c
         if (edep.ge.ecountmax) then
c
             ecountmax = edep
             xcountmax = x(np)
             ycountmax = y(np)
             zcountmax = z(np)
c
         endif
c
      elseif (index(ir(np))/10000.eq.labelbiomass) then 
c
          totebio = totebio + edep
          biosource = biosource + edep
c
          edepr=edep
          if (edepr.ne.0.) call hfill (122,edepr,0.,1.)        
c
      else
c
          toteew = toteew + edep 
c
      endif
c
 10   return
      end
c
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
      subroutine msdist_pI(e0,eloss,tustep_orig,rhof_orig,medium_orig,
     *qel,spin_effects_orig,u0,v0,w0,x0,y0,z0,us,vs,ws,xf,yf,zf,
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
      real*4 blccc,xcccc,z_orig,r,z2,r2,r2max,chia2,chilog,cphi0,
     * cphi, sphi, e_orig, elke_orig, beta2, etap, xi_corr, ms_corr, 
     *epsilon,  temp,  factor,lambda,p2,p2i,q1,  rhophi2,  sint,  
     *sint0,  sint02,  sint0i,sphi0,u2p,ut,vt,wt_orig,xi,xphi,  
     *xphi2,  yphi,  yphi2
      logical find_index,  spin_index
      integer*4 lelke
      blccc = blcc(medium)
      xcccc = xcc(medium)
      e_orig = e0 - 0.5*eloss
      p2 = e_orig*(e_orig + rmt2)
      p2i = 1/p2
      chia2 = xcccc*p2i/(4*blccc)
      beta2 = p2/(p2 + rmsq)
      lambda = tustep_orig*rhof_orig*blccc/beta2
      factor = 1/(1 + 0.9784671*e_orig)
      epsilon= eloss/e0
      epsilon= epsilon/(1-0.5*epsilon)
      temp = 0.25*(1 - factor*(1 - 0.333333*factor))*epsilon**2
      lambda = lambda*(1 + temp)
      IF (( spin_effects_orig )) THEN
        elke_orig = Log(e_orig)
        Lelke=eke1(MEDIUM)*elke_orig+eke0(MEDIUM)
        IF (( lelke .LT. 1 )) THEN
          lelke = 1
          elke_orig = (1 - eke0(medium))/eke1(medium)
        END IF
        IF (( qel .EQ. 0 )) THEN
          etap=etae_ms1(Lelke,MEDIUM)*elke_orig+etae_ms0(Lelke,MEDIUM)
          xi_corr=q1ce_ms1(Lelke,MEDIUM)*elke_orig+
     *q1ce_ms0(Lelke,MEDIUM)
        ELSE
          etap=etap_ms1(Lelke,MEDIUM)*elke_orig+etap_ms0(Lelke,MEDIUM)
          xi_corr=q1cp_ms1(Lelke,MEDIUM)*elke_orig+
     *q1cp_ms0(Lelke,MEDIUM)
        END IF
        ms_corr=blcce1(Lelke,MEDIUM)*elke_orig+blcce0(Lelke,MEDIUM)
      ELSE
        etap = 1
        xi_corr = 1
        ms_corr = 1
      END IF
      chia2 = xcccc*p2i/(4*blccc)*etap
      lambda = lambda/etap/(1+chia2)*ms_corr
      chilog = Log(1 + 1/chia2)
      q1 = 2*chia2*(chilog*(1 + chia2) - 1)
      xi = q1*lambda
      find_index = .true.
      spin_index = .true.
      call mscat(lambda,chia2,xi,elke_orig,beta2,qel,medium,
     * spin_effects_orig,find_index,spin_index, ws,sint)
4481  CONTINUE
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
        IF(rhophi2.LE.1)GO TO4482
      GO TO 4481
4482  CONTINUE
      rhophi2 = 1/rhophi2
      cphi = (xphi2 - yphi2)*rhophi2
      sphi = 2*xphi*yphi*rhophi2
      us = sint*cphi
      vs = sint*sphi
      xi = xi*xi_corr
      IF (( xi .LT. 0.1 )) THEN
        z_orig = 1 - xi*(0.5 - xi*(0.166666667 - 0.041666667*xi))
      ELSE
        z_orig = (1 - Exp(-xi))/xi
      END IF
      r = 0.5*sint
      r2 = r*r
      z2 = z_orig*z_orig
      r2max = 1 - z2
      IF (( r2max .LT. r2 )) THEN
        r2 = r2max
        r = Sqrt(r2)
      END IF
      ut = r*cphi
      vt = r*sphi
      wt_orig = z_orig
      ustep_orig = Sqrt(z2 + r2)*tustep_orig
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
      subroutine mscat(lambda,chia2,q1,elke_orig,beta2,qel,medium_orig
     *,spin_effects_orig,find_index,spin_index, cost,sint)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c

      real*4 lambda, chia2,q1,elke_orig,beta2,cost,sint
      integer*4 qel,medium_orig
      logical spin_effects_orig,find_index,spin_index
      real*4 sprob,explambda,wsum,wprob,xi,rejf,spin_rejection, cosz,sin
     *z,phi,omega2,llmbda,ai,aj,ak,a,u_orig,du,x1,rnno
      integer*4 icount,i,j,k
      save i,j,omega2
      IF ((lambda .LE. 13.8)) THEN
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        sprob = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        explambda = Exp(-lambda)
        IF ((sprob .LT. explambda)) THEN
          cost = 1
          sint = 0
          return
        END IF
        wsum = (1+lambda)*explambda
        IF (( sprob .LT. wsum )) THEN
4010      CONTINUE
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          xi = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          xi = 2*chia2*xi/(1 - xi + chia2)
          cost = 1 - xi
          IF (( spin_effects_orig )) THEN
            rejf = spin_rejection(qel,medium_orig,elke_orig,beta2,
     *       q1,cost,spin_index,.false.)
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            rnno = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            IF (( rnno .GT. rejf )) THEN
              GOTO 4010
            END IF
          END IF
          sint = sqrt(xi*(2 - xi))
          return
        END IF
        IF (( lambda .LE. 1 )) THEN
          wprob = explambda
          wsum = explambda
          cost = 1
          sint = 0
          icount = 0
4021      CONTINUE
            icount = icount + 1
            IF((icount .GT. 20))GO TO4022
            wprob = wprob*lambda/icount
            wsum = wsum + wprob
4030        CONTINUE
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            xi = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            xi = 2*chia2*xi/(1 - xi + chia2)
            cosz = 1 - xi
            IF (( spin_effects_orig )) THEN
              rejf = spin_rejection(qel,medium_orig,elke_orig,beta2
     *         ,q1,cosz,spin_index,.false.)
              IF (( rng_seed .GT. 24 )) THEN
                call ranlux(rng_array)
                rng_seed = 1
              END IF
              rnno = rng_array(rng_seed)
              rng_seed = rng_seed + 1
              IF (( rnno .GT. rejf )) THEN
                GOTO 4030
              END IF
            END IF
            sinz = xi*(2 - xi)
            IF (( sinz .GT. 1.e-20 )) THEN
              sinz = Sqrt(sinz)
              IF (( rng_seed .GT. 24 )) THEN
                call ranlux(rng_array)
                rng_seed = 1
              END IF
              xi = rng_array(rng_seed)
              rng_seed = rng_seed + 1
              phi = xi*6.2831853
              cost = cost*cosz - sint*sinz*Cos(phi)
              sint = Sqrt(Max(0.0,(1-cost)*(1+cost)))
            END IF
            IF((( wsum .GT. sprob)))GO TO4022
          GO TO 4021
4022      CONTINUE
          return
        END IF
      END IF
      IF ((lambda .LE. 1e5 )) THEN
        IF ((find_index)) THEN
          llmbda = log(lambda)
          ai = llmbda*dllambi
          i = ai
          ai = ai - i
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          xi = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF((xi .LT. ai))i = i + 1
          IF (( q1 .LT. 1e-3 )) THEN
            j = 0
          ELSE IF(( q1 .LT. 0.5 )) THEN
            aj = q1*dqmsi
            j = aj
            aj = aj - j
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            xi = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            IF((xi .LT. aj))j = j + 1
          ELSE
            j = 7
          END IF
          IF ((llmbda .LT. 2.2299)) THEN
            omega2 = chia2*(lambda + 4)*(1.347006 + llmbda*( 0.209364 -
     *      llmbda*(0.45525 - llmbda*(0.50142 - 0.081234*llmbda))))
          ELSE
            omega2 = chia2*(lambda + 4)*(-2.77164 + llmbda*(2.94874 - 
     *      llmbda*(0.1535754 - llmbda*0.00552888)))
          END IF
          find_index = .false.
        END IF
4040    CONTINUE
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        xi = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        ak = xi*31
        k = ak
        ak = ak - k
        IF((ak .GT. wms_array(i,j,k)))k = ims_array(i,j,k)
        a = fms_array(i,j,k)
        u_orig = ums_array(i,j,k)
        du = ums_array(i,j,k+1) - u_orig
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        xi = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        IF (( abs(a) .LT. 0.2 )) THEN
          x1 = 0.5*(1-xi)*a
          u_orig = u_orig + xi*du*(1+x1*(1-xi*a))
        ELSE
          u_orig = u_orig - du/a*(1-Sqrt(1+xi*a*(2+a)))
        END IF
        xi = omega2*u_orig/(1 + 0.5*omega2 - u_orig)
        cost = 1 - xi
        IF (( spin_effects_orig )) THEN
          rejf=spin_rejection(qel,medium_orig,elke_orig,beta2,q1,cost,
     *    spin_index,.false.)
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          rnno = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF (( rnno .GT. rejf )) THEN
            GOTO 4040
          END IF
        END IF
        sint = sqrt(xi*(2-xi))
        return
      END IF
      write(6,*) ' '                                                            
      write(6,*) ' *************************************'                       
      write(6,*) ' Maximum step size in mscat exceeded! '                       
      write(6,*) ' Maximum step size initialized: 100000'                       
      write(6,*) ' Present lambda: ',lambda                                     
      write(6,*) ' chia2: ',chia2                                               
      write(6,*) ' q1 elke beta2: ',q1,elke,beta2                               
      write(6,*) ' medium: ',medium_orig                                             
      write(6,*) ' Stopping execution'                                          
      write(6,*) ' *************************************'                       
      stop
      end
      real*4 function spin_rejection(qel,medium_orig,elke_orig,beta2,q1,
     *cost, spin_index,is_single)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      real*4 elke_orig,beta2,q1,cost
      integer*4 qel,medium_orig
      logical spin_index,is_single
      real*4 rnno,ai,qq1,aj,xi,ak
      integer*4 i,j,k
      save i,j
      IF (( spin_index )) THEN
        spin_index = .false.
        IF (( beta2 .GE. b2spin_min )) THEN
          ai = (beta2 - b2spin_min)*dbeta2i
          i = ai
          ai = ai - i
          i = i + 15 + 1
        ELSE IF(( elke_orig .GT. espml )) THEN
          ai = (elke_orig - espml)*dleneri
          i = ai
          ai = ai - i
        ELSE
          i = 0
          ai = 0
        END IF
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        rnno = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        IF((rnno .LT. ai))i = i + 1
        IF (( is_single )) THEN
          j = 0
        ELSE
          qq1 = 2*q1
          qq1 = qq1/(1 + qq1)
          aj = qq1*dqq1i
          j = aj
          aj = aj - j
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          rnno = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF((rnno .LT. aj))j = j + 1
        END IF
      END IF
      xi = Sqrt(0.5*(1-cost))
      ak = xi*31
      k = ak
      ak = ak - k
      spin_rejection = (1-ak)*spin_rej(medium_orig,qel,i,j,k) + 
     *ak*spin_rej(medium_orig,qel,i,j,k+1)
      return
      end
      subroutine sscat(chia2,elke_orig,beta2,qel,medium_orig,
     *spin_effects_orig,cost,sint)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      
      include 'egs4fcomm.for'
c

      real*4 chia2,elke_orig,beta2,cost,sint
      integer*4 qel,medium_orig
      logical spin_effects_orig
      real*4 xi,rnno,rejf,spin_rejection
      logical spin_index
      spin_index = .true.
4050  CONTINUE
      IF (( rng_seed .GT. 24 )) THEN
        call ranlux(rng_array)
        rng_seed = 1
      END IF
      xi = rng_array(rng_seed)
      rng_seed = rng_seed + 1
      xi = 2*chia2*xi/(1 - xi + chia2)
      cost = 1 - xi
      IF (( spin_effects_orig )) THEN
        rejf = spin_rejection(qel,medium_orig,elke_orig,beta2,0.,cost
     *  ,spin_index,.true.)
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        rnno = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        IF((rnno .GT. rejf))goto 4050
      END IF
      sint = sqrt(xi*(2 - xi))
      return
      end
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
      SUBROUTINE ANNIH
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      DOUBLE PRECISION PAVIP,  PESG1,  PESG2
      real*4 AVIP, A, G, T, P, POT, EP0, WSAMP, RNNO01, RNNO02, EP, REJF
     *, ESG1, ESG2, aa, bb, cc, sinpsi, sindel, cosdel, us, vs,cphi,sphi
      integer*4 ibr
      real*4 xphi, xphi2, yphi, yphi2, rhophi2
      NPold = NP
      PAVIP=E(NP)+PRM
      AVIP=PAVIP
      A=AVIP/RM
      G=A-1.0
      T=G-1.0
      P=SQRT(A*T)
      POT=P/T
      EP0=1.0/(A+P)
      WSAMP=LOG((1.0-EP0)/EP0)
      aa = u(np)
      bb = v(np)
      cc = w(np)
      sinpsi = aa*aa + bb*bb
      IF (( sinpsi .GT. 1e-20 )) THEN
        sinpsi = sqrt(sinpsi)
        sindel = bb/sinpsi
        cosdel = aa/sinpsi
      END IF
      IF (( nbr_split .GT. 1 )) THEN
        wt(np) = wt(np)/nbr_split
      END IF
        DO 2811 ibr=1,nbr_split
        IF (( np+1 .GT. 50 )) THEN
          WRITE(6,2820)np+1
2820      FORMAT(//' Stack overflow in ANNIH! np = ',i6, ' Increase $MXS        
     *TACK and try again'//)                                                    
          stop
        END IF
2831    CONTINUE
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          RNNO01 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          EP=EP0*EXP(RNNO01*WSAMP)
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          RNNO02 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          REJF = 1 - (EP*A-1)**2/(EP*(A*A-2))
          IF(((RNNO02 .LE. REJF)))GO TO2832
        GO TO 2831
2832    CONTINUE
        ESG1=AVIP*EP
        PESG1=ESG1
        E(NP)=PESG1
        IQ(NP)=0
        X(np)=X(NPold)
        Y(np)=Y(NPold)
        Z(np)=Z(NPold)
        IR(np)=IR(NPold)
        WT(np)=WT(NPold)
        DNEAR(np)=DNEAR(NPold)
        LATCH(np)=LATCH(NPold)
        COSTHE=MIN(1.0,(ESG1-RM)*POT/ESG1)
        SINTHE=SQRT(1.0-COSTHE*COSTHE)
2841    CONTINUE
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
          IF(rhophi2.LE.1)GO TO2842
        GO TO 2841
2842    CONTINUE
        rhophi2 = 1/rhophi2
        cphi = (xphi2 - yphi2)*rhophi2
        sphi = 2*xphi*yphi*rhophi2
        IF (( sinpsi .GE. 1e-10 )) THEN
          us = sinthe*cphi
          vs = sinthe*sphi
          u(np) = cc*cosdel*us - sindel*vs + aa*costhe
          v(np) = cc*sindel*us + cosdel*vs + bb*costhe
          w(np) = cc*costhe - sinpsi*us
        ELSE
          u(np) = sinthe*cphi
          v(np) = sinthe*sphi
          w(np) = costhe
        END IF
        np = np + 1
        PESG2=PAVIP-PESG1
        esg2 = pesg2
        e(np) = pesg2
        iq(np) = 0
        X(np)=X(NPold)
        Y(np)=Y(NPold)
        Z(np)=Z(NPold)
        IR(np)=IR(NPold)
        WT(np)=WT(NPold)
        DNEAR(np)=DNEAR(NPold)
        LATCH(np)=LATCH(NPold)
        COSTHE=MIN(1.0,(ESG2-RM)*POT/ESG2)
        SINTHE=-SQRT(1.0-COSTHE*COSTHE)
        IF (( sinpsi .GE. 1e-10 )) THEN
          us = sinthe*cphi
          vs = sinthe*sphi
          u(np) = cc*cosdel*us - sindel*vs + aa*costhe
          v(np) = cc*sindel*us + cosdel*vs + bb*costhe
          w(np) = cc*costhe - sinpsi*us
        ELSE
          u(np) = sinthe*cphi
          v(np) = sinthe*sphi
          w(np) = costhe
        END IF
        np = np + 1
2811  CONTINUE
2812  CONTINUE
      np = np-1
      RETURN
      END
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
c
      subroutine howfar
c
c-----------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
c
c
c-----------------------------------------------------------c
c
      real xproxym, yproxym, zproxym, xinterm, yinterm, zinterm
      real xremote, yremote, zremote, xremsave, yremsave, zremsave
      data xproxym,yproxym,zproxym    / 0., 0., 0. /
      data xinterm,yinterm,zinterm    / 0., 0., 0. /
      data xremote,yremote,zremote    / 0., 0., 0. /
      data xremsave,yremsave,zremsave / 0., 0., 0. /
c
      real dstrack 
      integer loopstep, nactual, nremote, nproxym, nremsave
      data dstrack ,loopstep    / 0.0 , 0 /
      data nactual,nremote,nproxym,nremsave / 4*0 /
c
      integer ninterm    
c
c-----------------------------------------------------------------------c
c
      nactual  =  0
      loopstep =  0
      idisc    = 0
c
      nactual = kgeom(x(np),y(np),z(np))
c mirko
c      write (*,*) ' in howfar',x(np), y(np), z(np), u(np), v(np), w(np)
c     a, nactual, iq(np)
c
      if (  nactual.gt.nreg ) then
          idisc = 1000
      else  
cnico          dnear (np) = 0.
          dstrack = amin1(ustep,tracklim)
          xremote = dstrack * u(np) + x(np)
          yremote = dstrack * v(np) + y(np)
          zremote = dstrack * w(np) + z(np)
          nremote = kgeom(xremote,yremote,zremote)
c
          if ( nremote .ne. nactual ) then
             xproxym = x(np)
             yproxym = y(np)
             zproxym = z(np)
             nproxym = nactual
c
   10        continue
                xinterm = 0.5 * (xremote+xproxym)
                yinterm = 0.5 * (yremote+yproxym)
                zinterm = 0.5 * (zremote+zproxym)
                ninterm = kgeom(xinterm,yinterm,zinterm)
c
                if( ninterm.eq.nactual ) then
                   xproxym = xinterm
                   yproxym = yinterm
                   zproxym = zinterm
                else
                   xremote = xinterm
                   yremote = yinterm
                   zremote = zinterm
                   nremote = ninterm
                endif
c
                loopstep = loopstep+1
                dstrack=sqrt((xremote-xproxym)*(xremote-xproxym)
     a                      +(yremote-yproxym)*(yremote-yproxym)
     b                      +(zremote-zproxym)*(zremote-zproxym))
                if( abs(dstrack)-steplim ) 20,10,10
   20        continue
c
             if(loopstep.ge.64) then 
               write(*,*) 'convergence missing '
     a                                ,loopstep,dstrack
c
               write (*,*) nevent, nactual, nremote, x(np), y(np), z(np)
     a,                    xremote, yremote, zremote, dstrack, tracklim
             endif
             if (nremote.le.maxreg) then
                irnew = nremote
                ustep = sqrt ( ( xremote - x(np) ) * ( xremote - x(np) )
     a                      + ( yremote - y(np) ) * ( yremote - y(np) )
     b                      + ( zremote - z(np) ) * ( zremote - z(np) ))
                idisc = 0
cnico
		dnear(np) = 0.	
cnico
             else
	        idisc = 1000
	     endif
          endif
      endif
c mirko
c      write (*,*) 'out howfar',xremote,yremote,zremote,nremote,iq(np)
c     a, idisc
c
      return
c
      end
c
      subroutine hownear(tperp)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
c
c
cale
c
      real*4 tperp
      real r_temp
cale
      if (ir(np).ge.maxreg) then
        tperp=0.0
      elseif (ir(np).gt.regoffset) then
cnico particle inside the collimator
       	if(modlabel.eq.1) then 
c
c     moduli esagonali
c       
      	   if(med(ir(np)).eq.indair) then
c		tperp = dmin1(dnear_temp,xair-dnear_temp)
		tperp = dnear_temp
      	   elseif(med(ir(np)).eq.indlead) then
		tperp = dmin1(dnear_temp,xlead-dnear_temp)
	   else
		write(*,*) 'Error in hownear'
                write (*,*)'med(np)',med(np)
                write (*,*)'iregtype',iregtype(ir(np))
                write (*,*)'ir(np)',ir(np)
                write (*,*)'x y z',x(np),y(np),z(np)
           endif
c
c    fine moduli esagonali       
c
       	elseif(modlabel.eq.2) then 
c
c     moduli rettangolari
c
      	   if(med(ir(np)).eq.indair) then
		tperp = dmin1(xtemp,ytemp,(xair-xtemp),(yair-ytemp))
      	   elseif(med(ir(np)).eq.indlead) then
		tperp = dmin1((xtemp-xair),(ytemp-yair),
     a		         (xcell-xtemp),(ycell-ytemp))
	   endif
c
c    fine moduli rettangolari       
c
	endif
      elseif (ir(np).gt.regoffset_det1mod) then
cnico particle inside the detector
        tperp = amin1( abs(x(np)-xinf(ir(np))) , abs(x(np)-xsup(ir(np)))
     a               , abs(y(np)-yinf(ir(np))) , abs(y(np)-ysup(ir(np)))
     b              , abs(z(np)-zinit(ir(np))), abs(z(np)-zend(ir(np))))
      elseif (iregtype(ir(np)).eq.1)then
        tperp = amin1( abs(x(np)-xinf(ir(np))) , abs(x(np)-xsup(ir(np)))
     a               , abs(y(np)-yinf(ir(np))) , abs(y(np)-ysup(ir(np)))
     b              , abs(z(np)-zinit(ir(np))), abs(z(np)-zend(ir(np))))
c
      elseif (iregtype(ir(np)).eq.2)then
           r_temp=radius(ir(np))-sqrt(((xsup(ir(np))-x(np))*
     a                    (xsup(ir(np))-x(np)))+((y(np)-yinf(ir(np)))*
     b       (y(np)-yinf(ir(np)))))
             tperp=amin1(abs(z(np)-zinit(ir(np))),
     a                  abs(z(np)-zend(ir(np))),r_temp)

      elseif (iregtype(ir(np)).eq.3)then
           r_temp=radius(ir(np))-sqrt(((zend(ir(np))-z(np))*
     a                   (zend(ir(np))-z(np)))+((y(np)-yinf(ir(np)))*
     b       (y(np)-yinf(ir(np)))))
             tperp=amin1(abs(x(np)-xinf(ir(np))),
     a                  abs(x(np)-xsup(ir(np))),r_temp)

cale      endif
      else
	write(*,*) 'Error in hownear ',ir(np),iregtype(ir(np))
 
      endif

      return
c
      end
      SUBROUTINE PHOTON(IRCODE)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      integer*4 IRCODE
      DOUBLE PRECISION PEIG
      real*4 EIG, RNNO35, DPMFP, GMFPR0, GMFP, COHFAC, RNNO37, XXX, X2
     *, Q2, CSQTHE, REJF, RNNORJ, RNNO36, GBR1, GBR2, T
      integer*4 IARG, IDR, IRL, LGLE, LXXX
      IRCODE=1
      PEIG=E(NP)
      EIG=PEIG
      IRL=IR(NP)
      MEDIUM=MED(IRL)
      IF ((EIG.LE.PCUT(IRL))) THEN
        GO TO 4860
      END IF
4870  CONTINUE
4871    CONTINUE
        GLE=LOG(EIG)
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        RNNO35 = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        IF ((RNNO35.EQ.0.0)) THEN
          RNNO35=1.E-30
        END IF
        DPMFP=-LOG(RNNO35)
        IROLD=IR(NP)
4880    CONTINUE
4881      CONTINUE
          IF ((MEDIUM.NE.0)) THEN
            LGLE=GE1(MEDIUM)*GLE+GE0(MEDIUM)
            GMFPR0=GMFP1(LGLE,MEDIUM)*GLE+GMFP0(LGLE,MEDIUM)
          END IF
4890      CONTINUE
4891        CONTINUE
            IF ((MEDIUM.EQ.0)) THEN
              TSTEP=VACDST
            ELSE
              RHOF=RHOR(IRL)/RHO(MEDIUM)
              GMFP=GMFPR0/RHOF
              IF ((IRAYLR(IRL).EQ.1)) THEN
                COHFAC=COHE1(LGLE,MEDIUM)*GLE+COHE0(LGLE,MEDIUM)
                GMFP=GMFP*COHFAC
              END IF
              TSTEP=GMFP*DPMFP
            END IF
            IRNEW=IR(NP)
            IDISC=0
            USTEP=TSTEP
            TUSTEP=USTEP
            IF (( ustep .GT. dnear(np) .OR. wt(np) .LE. 0 )) THEN
              call howfar
            END IF
            IF ((IDISC.GT.0)) THEN
              GO TO 4900
            END IF
            VSTEP=USTEP
            TVSTEP=VSTEP
            EDEP=PZERO
            IARG=0
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
            X(NP)=X(NP)+U(NP)*USTEP
            Y(NP)=Y(NP)+V(NP)*USTEP
            Z(NP)=Z(NP)+W(NP)*USTEP
            DNEAR(NP)=DNEAR(NP)-USTEP
            IF ((MEDIUM.NE.0)) THEN
              DPMFP=MAX(0.,DPMFP-USTEP/GMFP)
            END IF
            IROLD=IR(NP)
            MEDOLD=MEDIUM
            IF ((IRNEW.NE.IROLD)) THEN
              IR(NP)=IRNEW
              IRL=IRNEW
              MEDIUM=MED(IRL)
            END IF
            IARG=5
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
            IF ((EIG.LE.PCUT(IRL))) THEN
              GO TO 4860
            END IF
            IF((IDISC.LT.0))GO TO 4900
            IF((MEDIUM.NE.MEDOLD))GO TO 4892
            IF ((MEDIUM.NE.0.AND.DPMFP.LE.1.E-5)) THEN
              GO TO 4882
            END IF
          GO TO 4891
4892      CONTINUE
        GO TO 4881
4882    CONTINUE
        IF ((IRAYLR(IRL).EQ.1)) THEN
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          RNNO37 = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF ((RNNO37.LE.(1.0-COHFAC))) THEN
            IARG=23
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
            NPold = NP
4910        CONTINUE
4911          CONTINUE
              IF (( rng_seed .GT. 24 )) THEN
                call ranlux(rng_array)
                rng_seed = 1
              END IF
              XXX = rng_array(rng_seed)
              rng_seed = rng_seed + 1
              LXXX=RCO1(MEDIUM)*XXX+RCO0(MEDIUM)
              X2=RSCT1(LXXX,MEDIUM)*XXX+RSCT0(LXXX,MEDIUM)
              Q2=X2*RMSQ/(20.60744*20.60744)
              COSTHE=1.-Q2/(2.*E(NP)*E(NP))
              IF((ABS(COSTHE).GT.1.0))GO TO 4910
              CSQTHE=COSTHE*COSTHE
              REJF=(1.0+CSQTHE)/2.0
              IF (( rng_seed .GT. 24 )) THEN
                call ranlux(rng_array)
                rng_seed = 1
              END IF
              RNNORJ = rng_array(rng_seed)
              rng_seed = rng_seed + 1
              IF(((RNNORJ .LE. REJF)))GO TO4912
            GO TO 4911
4912        CONTINUE
            SINTHE=SQRT(1.0-CSQTHE)
            CALL UPHI(2,1)
            IARG=24
            IF ((IAUSFL(IARG+1).NE.0)) THEN
              CALL AUSGAB(IARG)
            END IF
            GOTO 4870
          END IF
        END IF
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        RNNO36 = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        GBR1=GBR11(LGLE,MEDIUM)*GLE+GBR10(LGLE,MEDIUM)
        IF (((RNNO36.LE.GBR1).AND.(E(NP).GT.RMT2) )) THEN
          IARG=15
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          CALL PAIR
          IARG=16
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          IF (( iq(np) .NE. 0 )) THEN
            GO TO 4872
          ELSE
            goto 4920
          END IF
        END IF
        GBR2=GBR21(LGLE,MEDIUM)*GLE+GBR20(LGLE,MEDIUM)
        IF ((RNNO36.LT.GBR2)) THEN
          IARG=17
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          CALL COMPT
          IARG=18
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          IF((IQ(NP).NE.0))GO TO 4872
        ELSE
          IARG=19
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          CALL PHOTO
          IF ((NP.EQ.0)) THEN
            IRCODE=2
            RETURN
          END IF
          IARG=20
          IF ((IAUSFL(IARG+1).NE.0)) THEN
            CALL AUSGAB(IARG)
          END IF
          IF((IQ(NP) .NE. 0))GO TO 4872
        END IF
4920    PEIG=E(NP)
        EIG=PEIG
        IF((EIG.LT.PCUT(IRL)))GO TO 4860
      GO TO 4871
4872  CONTINUE
      RETURN
4860  IF ((EIG.GT.AP(MEDIUM))) THEN
        IDR=1
      ELSE
        IDR=2
      END IF
      EDEP=PEIG
      IARG=IDR
      IF ((IAUSFL(IARG+1).NE.0)) THEN
        CALL AUSGAB(IARG)
      END IF
      IRCODE=2
      NP=NP-1
      RETURN
4900  EDEP=PEIG
      IARG=3
      IF ((IAUSFL(IARG+1).NE.0)) THEN
        CALL AUSGAB(IARG)
      END IF
      IRCODE=2
      NP=NP-1
      RETURN
      END
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
      subroutine evsum
c
c---------------------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
c
      integer mh
c
c----- energy loss in all regions
c
      nhit         = 0
      multiplicity = 0
c
      if (biosource.gt.0.) call hfill(123,biosource,0.,1.)
      if (elost.gt.0.) call hfill(124,elost,0.,1.)
c
      epixel=0.0
      xpixel=0.0
      ypixel=0.0
      zpixel=0.0
c
      epixmax=0.0
      xpixmax=0.0
      ypixmax=0.0
      zpixmax=0.0
c
cnico nico cambia nplan con regoffset
c      do 10 mh=1,nplan
      do 10 mh=1,regoffset
c
         if ((index(mh)/10000).eq.labeldet) then
c
             if ((nhit.le.maxhit).and.(deodx(mh).gt.trigger)) then
c
                 nhit  = nhit + 1
                 hiteng(nhit)     = deodx(mh)
                 indexcount(nhit) = index(mh)
c
             endif
c
             if ((deodx(mh).gt.esnr).and.(nhit.le.maxhit)) then 
c
		 multiplicity = multiplicity +1
                 epixel=epixel+deodx(mh)
                 xpixel=xpixel+deodx(mh)*(xinf(mh)+xsup(mh))/2.
                 ypixel=ypixel+deodx(mh)*(yinf(mh)+ysup(mh))/2.
                 zpixel=zpixel+deodx(mh)*(zinit(mh)+zend(mh))/2.
c
                 if (deodx(mh).ge.epixmax) then
c
                    epixmax = deodx(mh)
                    xpixmax = (xinf(mh)+xsup(mh))/2.
                    ypixmax = (yinf(mh)+ysup(mh))/2.
                    zpixmax = (zinit(mh)+zend(mh))/2.
c
                 endif
c
             endif
c
         endif
c
 10   continue
c
c---------------------------------------------------------------------c
c
      if ( (elost.gt.treshold).and.(nhit.ge.1) ) then
c
           ilsx = ilsx +1
c
           call tuplesetup
           call hitorder
c
      endif
c
      return
c
      end
c
      subroutine hitorder
c
c-----------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
c
      integer khit, jhit, loop
      real eng_max
c
c----  copy to temporary array 
c
      do 100 khit = 1 , maxhit
c
      temp_eng ( khit ) = hiteng (khit )
      indexord ( khit ) = 0     
c 
  100 continue
c
c----   ordering twofold loop 
c
      do 200 khit = 1,nhit
c
          eng_max = -100.0
          loop = 0
c
          do 220  jhit = 1, nhit
c
              if ( eng_max .lt. temp_eng(jhit) ) then
c
                 eng_max =  temp_eng(jhit)
                 loop    =  jhit
c
              endif                  
c
  220     continue
c
c----      end of inner loop
c
          if ( loop.gt.0) then
c
             indexord (khit ) =   loop
             temp_eng (loop ) = - 100.0
c
          endif
c
c---   end of outer loop
c
  200 continue
c
      return
      END
c
      subroutine tuplesetup
c
c---------------------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
c
      real egaussian, sgaussian_x, sgaussian_y, tuptempvarx, tuptempvary
c
c---------------------------------------------------------------------c
c
      egaussian = 0.0
      call norran ( egaussian)
c
      sgaussian_x = 0.0
      call norran ( sgaussian_x)
c
      sgaussian_y = 0.0
      call norran ( sgaussian_y)
c
      xtup3(1) = nactsx
      xtup3(2) = nevent
c
      xtup3(3)  = xin
      xtup3(4)  = yin
      xtup3(5)  = zin
      xtup3(6)  = uin
      xtup3(7)  = vin
      xtup3(8)  = win
c
      xavg=xavg/eavg
      yavg=yavg/eavg
      zavg=zavg/eavg          
c
      xtup3( 9)  =  eavg
      xtup3(10)  =  xavg
      xtup3(11)  =  yavg
      xtup3(12)  =  zavg
cnico
cc      write(outtxt,100) xavg,yavg
cnico
c
cc      xtup3(13)  =  ecountmax
cc      xtup3(14)  =  xcountmax
cc      xtup3(15)  =  ycountmax
cc      xtup3(16)  =  zcountmax
c
cc      xtup3(17)  =  efirst
cc      xtup3(18)  =  xfirst
cc      xtup3(19)  =  yfirst
cc      xtup3(20)  =  zfirst
c
      xtup3(13)  = amax1( 0. , eavg + spren * egaussian)
c
      if (ispres.eq.1) then
         tuptempvarx = xavg + sgaussian_x * (zdetend - zavg) / 3
         tuptempvary = yavg + sgaussian_y * (zdetend - zavg) / 3
         if ((tuptempvarx.ge.xdetinf).and.(tuptempvarx.le.xdetsup).and.
     a      (tuptempvary.ge.ydetinf).and.(tuptempvary.le.ydetsup)) then
                xtup3(14)  = tuptempvarx
                xtup3(15)  = tuptempvary
         else
                xtup3(14)  = xdetinf - 1.
                xtup3(15)  = ydetinf - 1.
         endif
      else
         xtup3(14)  =  xavg
         xtup3(15)  =  yavg
      endif
c
cc      xtup3(24)  = amax1( 0. , ecountmax + spren * egaussian)
c
      if (ispres.eq.1) then
         tuptempvarx = xcountmax + sgaussian_x * (zdetend - zcountmax)/3
         tuptempvary = ycountmax + sgaussian_y * (zdetend - zcountmax)/3
         if ((tuptempvarx.ge.xdetinf).and.(tuptempvarx.le.xdetsup).and.
     a      (tuptempvary.ge.ydetinf).and.(tuptempvary.le.ydetsup)) then
cc                xtup3(25)  = tuptempvarx
cc                xtup3(26)  = tuptempvary
         else
cc                xtup3(25)  = xdetinf - 1.
cc                xtup3(26)  = ydetinf - 1.
         endif
      else
cc         xtup3(25)  =  xcountmax
cc         xtup3(26)  =  ycountmax
      endif
c
cc      xtup3(27)  = amax1( 0. , efirst + spren * egaussian)
c
      if (ispres.eq.1) then
         tuptempvarx =  xfirst + sgaussian_x * (zdetend - zfirst) / 3
         tuptempvary =  yfirst + sgaussian_y * (zdetend - zfirst) / 3
         if ((tuptempvarx.ge.xdetinf).and.(tuptempvarx.le.xdetsup).and.
     a      (tuptempvary.ge.ydetinf).and.(tuptempvary.le.ydetsup)) then
cc                xtup3(28)  = tuptempvarx
cc                xtup3(29)  = tuptempvary
         else
cc                xtup3(28)  = xdetinf - 1.
cc                xtup3(29)  = ydetinf - 1.
         endif
      else
cc         xtup3(28)  =  xfirst
cc         xtup3(29)  =  yfirst
      endif
c
      xpixel=xpixel/epixel
      ypixel=ypixel/epixel
      zpixel=zpixel/epixel
c
      xtup3(16)  =  epixel
      xtup3(17)  =  amax1( 0. , epixel + spren * egaussian)
      xtup3(18)  =  xpixel
      xtup3(19)  =  ypixel
c      xtup3(34)  =  zpixel
c
cc      xtup3(34)  =  epixmax
cc      xtup3(35)  =  amax1( 0. , epixmax + spren * egaussian)
cc      xtup3(36)  =  xpixmax
cc      xtup3(37)  =  ypixmax
c      xtup3(37)  =  zpixmax          
c
cc      xtup3(38)  = xavg - (xin+uin*(zavg-zin)/win)
cc      xtup3(39)  = yavg - (yin+vin*(zavg-zin)/win)
c
cc      xtup3(40)  = multiplicity
c
      call hfn  (nt_hist,xtup3)
c
      return
c
  100 FORMAT(2(F8.4,1X))
      END
c
