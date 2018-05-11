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


