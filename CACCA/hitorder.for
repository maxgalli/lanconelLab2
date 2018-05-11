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
