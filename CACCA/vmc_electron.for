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
