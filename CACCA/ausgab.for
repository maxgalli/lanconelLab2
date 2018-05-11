c
      subroutine ausgab ( iarg )
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
c
c
      integer iarg
      real edepr
c          
c-----------------------------------------------------------c
c----------    tracking section ----------------------------c
c
c  keep a running sum of the energy deposited in each region
c
c--------------statistics of regions -----------------------c
c
c mirko testing
c
c       if ((abs(w(np)).gt.1.0  ).or.(abs(u(np)).gt.1.0  ).or.
c     a     (abs(v(np)).gt.1.0  ).or.(abs(x(np)).gt.1000.).or.
c     b     (abs(y(np)).gt.1000.).or.(abs(z(np)).gt.1000.))  then
c       write(*,*) 'ausgab  ', x(np), y(np), z(np), u(np), v(np), w(np)
c     a, np, iarg, edep, idisc, iq(np)
c
c       endif
c
c      if ((iarg.eq.0).and.(edep.ne.0.d0)) write(*,*) 'iarg = 0', edep
c
c end - mirko testing
c
      if (edep.eq.0.d0) goto 10
c
      if (idisc.eq.1000) then 
c
          toteoob = toteoob + edep
          goto 10
c
      endif
c
cnico nico cambia nplan con regoffset
c      if(ir(np).le.nplan) then
      if(ir(np).le.regoffset) then
c
        deodx(ir(np)) = deodx(ir(np)) + edep
        totedeodx = totedeodx + edep
c
      else
c
        totemod = totemod + edep
c
      endif
c
c
      edepr = edep
      if (edepr.ne.0.) call hfill (120,edepr,0.,1.)
c
      if ((index(ir(np))/10000).eq.labeldet) then
c
          if (first_int_label.and.(edep.gt.0.d0)) then
c
             efirst = edep
             xfirst = x(np)
             yfirst = y(np)
             zfirst = z(np)
             first_int_label = .false.
c
          endif
c
          edepr=edep
          if (edepr.ne.0.) call hfill (121,edepr,0.,1.)
c
          totedet = totedet + edep
          elost     =  elost  +   edep
c
          eavg      =  eavg +  edep
          xavg      =  xavg +  edep * x(np)
          yavg      =  yavg +  edep * y(np)
          zavg      =  zavg +  edep * z(np)
c
         if (edep.ge.ecountmax) then
c
             ecountmax = edep
             xcountmax = x(np)
             ycountmax = y(np)
             zcountmax = z(np)
c
         endif
c
      elseif (index(ir(np))/10000.eq.labelbiomass) then 
c
          totebio = totebio + edep
          biosource = biosource + edep
c
          edepr=edep
          if (edepr.ne.0.) call hfill (122,edepr,0.,1.)        
c
      else
c
          toteew = toteew + edep 
c
      endif
c
 10   return
      end
c
