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
