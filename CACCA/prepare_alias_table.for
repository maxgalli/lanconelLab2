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
