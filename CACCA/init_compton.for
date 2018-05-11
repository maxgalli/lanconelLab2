      subroutine init_compton
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
cale
      include 'egs4comm.for'

      character*80 incohfilename

      character*80 totincohfilename

cale

c
      integer*4 i,j,iz,nsh
      real*4 aux,pztot
      real*4 aux_erf,erf
        DO 3711 j=1,numgeom
        medium = med(j)
        IF (( medium .GT. 0 .AND. medium .LE. nmed)) THEN
          IF (( ibcmp(j) .GT. 0 )) THEN
            goto 3720
          END IF
        END IF
3711  CONTINUE
3712  CONTINUE
      WRITE(6,3730)
3730  FORMAT(' Bound Compton scattering not requested! '/)                      
      return
3720  WRITE(6,3740)
3740  FORMAT(' Bound Compton scattering requested, reading data ......',        
     *$)
cale
      incohfilename='incoh.dat'

      totincohfilename=pathname(1:pathlength) // incohfilename
cale     write(*,*)'cale:',totincohfilename
      open(unit=78,file=totincohfilename,form='formatted',status='old')
cale
        DO 3751 j=1,18
        read(78,*)
3751  CONTINUE
3752  CONTINUE
        DO 3761 j=1,1538
        read(78,*) iz_array(j),shn_array(j),ne_array(j), Jo_array(j),be_
     *  array(j)
        Jo_array(j) = Jo_array(j)*137.
        be_array(j) = be_array(j)*1e-6/rm
        aux_erf = 0.70710678119*(1+0.3*Jo_array(j))
        erfJo_array(j) = 0.82436063535*(Erf(aux_erf)-1)
3761  CONTINUE
3762  CONTINUE
      WRITE(6,3770)
3770  FORMAT(' Done')                                                           
      WRITE(6,3780)
3780  FORMAT(/' Initializing Bound Compton scattering ......')                  
        DO 3791 medium=1,nmed
        pztot = 0
        nsh = 0
          DO 3801 i=1,nne(medium)
          iz = int(zelem(medium,i))
            DO 3811 j=1,1538
            IF (( iz .EQ. iz_array(j) )) THEN
              nsh = nsh + 1
              IF (( nsh .GT. 50 )) THEN
                WRITE(6,3820)medium,50
3820            FORMAT(/' For medium ',i3,' the number of shells is > ',        
     *          i4,'!')                                                         
                WRITE(6,3830)
3830            FORMAT(' Increase the parameter $MXMDSH! ')                     
                stop
              END IF
              shell_array(nsh,medium) = j
              aux = pz(medium,i)*ne_array(j)
              eno_array(nsh,medium) = aux
              pztot = pztot + aux
            END IF
3811      CONTINUE
3812      CONTINUE
3801    CONTINUE
3802    CONTINUE
        IF (( nsh .EQ. 0 )) THEN
          WRITE(6,3840)medium
3840      FORMAT(' Medium ',i3,' has zero shells! ')                            
          stop
        END IF
        n_shell(medium) = nsh
        WRITE(6,3850)medium,nsh
3850    FORMAT(' Medium ',i3,' has ',i3,' shells: ')                            
          DO 3861 i=1,nsh
          j = shell_array(i,medium)
          eno_array(i,medium) = eno_array(i,medium)/pztot
cale          WRITE(6,3870)i,j,shn_array(j),eno_array(i,medium), Jo_array(j)
cale     *    ,be_array(j)*rm*1000.
3870      FORMAT(i3,i4,i3,f9.5,e10.3,f10.3)
3861    CONTINUE
3862    CONTINUE
3791  CONTINUE
3792  CONTINUE
      WRITE(6,3880)
3880  FORMAT('...... Done.'/)                                                   
        DO 3891 j=1,numgeom
        medium = med(j)
        IF (( medium .GT. 0 .AND. medium .LE. nmed)) THEN
          IF (( iedgfl(j) .GT. 0 .AND. iedgfl(j) .LE. 100 )) THEN
            goto 3900
          END IF
        END IF
3891  CONTINUE
3892  CONTINUE
      WRITE(6,3910)
3910  FORMAT(' In subroutine init_compton: ',/ '   fluorescence not set         
     *but relaxation data are required for ',/ '   bound Compton scatter        
     *ing. ',/ '   calling EDGSET. '//)                                         
      iedgfl(1) = 1
      call edgset(1,1)
      iedgfl(1) = 0
3900  CONTINUE
      return
      end
