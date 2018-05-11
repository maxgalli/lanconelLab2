      subroutine mscat(lambda,chia2,q1,elke_orig,beta2,qel,medium_orig
     *,spin_effects_orig,find_index,spin_index, cost,sint)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c

      real*4 lambda, chia2,q1,elke_orig,beta2,cost,sint
      integer*4 qel,medium_orig
      logical spin_effects_orig,find_index,spin_index
      real*4 sprob,explambda,wsum,wprob,xi,rejf,spin_rejection, cosz,sin
     *z,phi,omega2,llmbda,ai,aj,ak,a,u_orig,du,x1,rnno
      integer*4 icount,i,j,k
      save i,j,omega2
      IF ((lambda .LE. 13.8)) THEN
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        sprob = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        explambda = Exp(-lambda)
        IF ((sprob .LT. explambda)) THEN
          cost = 1
          sint = 0
          return
        END IF
        wsum = (1+lambda)*explambda
        IF (( sprob .LT. wsum )) THEN
4010      CONTINUE
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          xi = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          xi = 2*chia2*xi/(1 - xi + chia2)
          cost = 1 - xi
          IF (( spin_effects_orig )) THEN
            rejf = spin_rejection(qel,medium_orig,elke_orig,beta2,
     *       q1,cost,spin_index,.false.)
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            rnno = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            IF (( rnno .GT. rejf )) THEN
              GOTO 4010
            END IF
          END IF
          sint = sqrt(xi*(2 - xi))
          return
        END IF
        IF (( lambda .LE. 1 )) THEN
          wprob = explambda
          wsum = explambda
          cost = 1
          sint = 0
          icount = 0
4021      CONTINUE
            icount = icount + 1
            IF((icount .GT. 20))GO TO4022
            wprob = wprob*lambda/icount
            wsum = wsum + wprob
4030        CONTINUE
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            xi = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            xi = 2*chia2*xi/(1 - xi + chia2)
            cosz = 1 - xi
            IF (( spin_effects_orig )) THEN
              rejf = spin_rejection(qel,medium_orig,elke_orig,beta2
     *         ,q1,cosz,spin_index,.false.)
              IF (( rng_seed .GT. 24 )) THEN
                call ranlux(rng_array)
                rng_seed = 1
              END IF
              rnno = rng_array(rng_seed)
              rng_seed = rng_seed + 1
              IF (( rnno .GT. rejf )) THEN
                GOTO 4030
              END IF
            END IF
            sinz = xi*(2 - xi)
            IF (( sinz .GT. 1.e-20 )) THEN
              sinz = Sqrt(sinz)
              IF (( rng_seed .GT. 24 )) THEN
                call ranlux(rng_array)
                rng_seed = 1
              END IF
              xi = rng_array(rng_seed)
              rng_seed = rng_seed + 1
              phi = xi*6.2831853
              cost = cost*cosz - sint*sinz*Cos(phi)
              sint = Sqrt(Max(0.0,(1-cost)*(1+cost)))
            END IF
            IF((( wsum .GT. sprob)))GO TO4022
          GO TO 4021
4022      CONTINUE
          return
        END IF
      END IF
      IF ((lambda .LE. 1e5 )) THEN
        IF ((find_index)) THEN
          llmbda = log(lambda)
          ai = llmbda*dllambi
          i = ai
          ai = ai - i
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          xi = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF((xi .LT. ai))i = i + 1
          IF (( q1 .LT. 1e-3 )) THEN
            j = 0
          ELSE IF(( q1 .LT. 0.5 )) THEN
            aj = q1*dqmsi
            j = aj
            aj = aj - j
            IF (( rng_seed .GT. 24 )) THEN
              call ranlux(rng_array)
              rng_seed = 1
            END IF
            xi = rng_array(rng_seed)
            rng_seed = rng_seed + 1
            IF((xi .LT. aj))j = j + 1
          ELSE
            j = 7
          END IF
          IF ((llmbda .LT. 2.2299)) THEN
            omega2 = chia2*(lambda + 4)*(1.347006 + llmbda*( 0.209364 -
     *      llmbda*(0.45525 - llmbda*(0.50142 - 0.081234*llmbda))))
          ELSE
            omega2 = chia2*(lambda + 4)*(-2.77164 + llmbda*(2.94874 - 
     *      llmbda*(0.1535754 - llmbda*0.00552888)))
          END IF
          find_index = .false.
        END IF
4040    CONTINUE
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        xi = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        ak = xi*31
        k = ak
        ak = ak - k
        IF((ak .GT. wms_array(i,j,k)))k = ims_array(i,j,k)
        a = fms_array(i,j,k)
        u_orig = ums_array(i,j,k)
        du = ums_array(i,j,k+1) - u_orig
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        xi = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        IF (( abs(a) .LT. 0.2 )) THEN
          x1 = 0.5*(1-xi)*a
          u_orig = u_orig + xi*du*(1+x1*(1-xi*a))
        ELSE
          u_orig = u_orig - du/a*(1-Sqrt(1+xi*a*(2+a)))
        END IF
        xi = omega2*u_orig/(1 + 0.5*omega2 - u_orig)
        cost = 1 - xi
        IF (( spin_effects_orig )) THEN
          rejf=spin_rejection(qel,medium_orig,elke_orig,beta2,q1,cost,
     *    spin_index,.false.)
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          rnno = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF (( rnno .GT. rejf )) THEN
            GOTO 4040
          END IF
        END IF
        sint = sqrt(xi*(2-xi))
        return
      END IF
      write(6,*) ' '                                                            
      write(6,*) ' *************************************'                       
      write(6,*) ' Maximum step size in mscat exceeded! '                       
      write(6,*) ' Maximum step size initialized: 100000'                       
      write(6,*) ' Present lambda: ',lambda                                     
      write(6,*) ' chia2: ',chia2                                               
      write(6,*) ' q1 elke beta2: ',q1,elke,beta2                               
      write(6,*) ' medium: ',medium_orig                                             
      write(6,*) ' Stopping execution'                                          
      write(6,*) ' *************************************'                       
      stop
      end
      real*4 function spin_rejection(qel,medium_orig,elke_orig,beta2,q1,
     *cost, spin_index,is_single)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      real*4 elke_orig,beta2,q1,cost
      integer*4 qel,medium_orig
      logical spin_index,is_single
      real*4 rnno,ai,qq1,aj,xi,ak
      integer*4 i,j,k
      save i,j
      IF (( spin_index )) THEN
        spin_index = .false.
        IF (( beta2 .GE. b2spin_min )) THEN
          ai = (beta2 - b2spin_min)*dbeta2i
          i = ai
          ai = ai - i
          i = i + 15 + 1
        ELSE IF(( elke_orig .GT. espml )) THEN
          ai = (elke_orig - espml)*dleneri
          i = ai
          ai = ai - i
        ELSE
          i = 0
          ai = 0
        END IF
        IF (( rng_seed .GT. 24 )) THEN
          call ranlux(rng_array)
          rng_seed = 1
        END IF
        rnno = rng_array(rng_seed)
        rng_seed = rng_seed + 1
        IF((rnno .LT. ai))i = i + 1
        IF (( is_single )) THEN
          j = 0
        ELSE
          qq1 = 2*q1
          qq1 = qq1/(1 + qq1)
          aj = qq1*dqq1i
          j = aj
          aj = aj - j
          IF (( rng_seed .GT. 24 )) THEN
            call ranlux(rng_array)
            rng_seed = 1
          END IF
          rnno = rng_array(rng_seed)
          rng_seed = rng_seed + 1
          IF((rnno .LT. aj))j = j + 1
        END IF
      END IF
      xi = Sqrt(0.5*(1-cost))
      ak = xi*31
      k = ak
      ak = ak - k
      spin_rejection = (1-ak)*spin_rej(medium_orig,qel,i,j,k) + 
     *ak*spin_rej(medium_orig,qel,i,j,k+1)
      return
      end
