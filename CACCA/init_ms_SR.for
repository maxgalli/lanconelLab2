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
