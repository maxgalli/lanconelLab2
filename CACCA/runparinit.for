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
