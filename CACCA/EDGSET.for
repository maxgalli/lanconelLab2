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
