      subroutine sscat(chia2,elke_orig,beta2,qel,medium_orig,
     *spin_effects_orig,cost,sint)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      
      include 'egs4fcomm.for'
c

      real*4 chia2,elke_orig,beta2,cost,sint
      integer*4 qel,medium_orig
      logical spin_effects_orig
      real*4 xi,rnno,rejf,spin_rejection
      logical spin_index
      spin_index = .true.
4050  CONTINUE
      IF (( rng_seed .GT. 24 )) THEN
        call ranlux(rng_array)
        rng_seed = 1
      END IF
      xi = rng_array(rng_seed)
      rng_seed = rng_seed + 1
      xi = 2*chia2*xi/(1 - xi + chia2)
      cost = 1 - xi
      IF (( spin_effects_orig )) THEN
        rejf = spin_rejection(qel,medium_orig,elke_orig,beta2,0.,cost
     *  ,spin_index,.true.)
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        rnno = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        IF((rnno .GT. rejf))goto 4050
      END IF
      sint = sqrt(xi*(2 - xi))
      return
      end
