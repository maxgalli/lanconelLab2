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


