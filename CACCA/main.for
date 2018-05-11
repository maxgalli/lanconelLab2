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
