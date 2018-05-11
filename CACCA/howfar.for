c
      subroutine howfar
c
c-----------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
c
c
c-----------------------------------------------------------c
c
      real xproxym, yproxym, zproxym, xinterm, yinterm, zinterm
      real xremote, yremote, zremote, xremsave, yremsave, zremsave
      data xproxym,yproxym,zproxym    / 0., 0., 0. /
      data xinterm,yinterm,zinterm    / 0., 0., 0. /
      data xremote,yremote,zremote    / 0., 0., 0. /
      data xremsave,yremsave,zremsave / 0., 0., 0. /
c
      real dstrack 
      integer loopstep, nactual, nremote, nproxym, nremsave
      data dstrack ,loopstep    / 0.0 , 0 /
      data nactual,nremote,nproxym,nremsave / 4*0 /
c
      integer ninterm    
c
c-----------------------------------------------------------------------c
c
      nactual  =  0
      loopstep =  0
      idisc    = 0
c
      nactual = kgeom(x(np),y(np),z(np))
c mirko
c      write (*,*) ' in howfar',x(np), y(np), z(np), u(np), v(np), w(np)
c     a, nactual, iq(np)
c
      if (  nactual.gt.nreg ) then
          idisc = 1000
      else  
cnico          dnear (np) = 0.
          dstrack = amin1(ustep,tracklim)
          xremote = dstrack * u(np) + x(np)
          yremote = dstrack * v(np) + y(np)
          zremote = dstrack * w(np) + z(np)
          nremote = kgeom(xremote,yremote,zremote)
c
          if ( nremote .ne. nactual ) then
             xproxym = x(np)
             yproxym = y(np)
             zproxym = z(np)
             nproxym = nactual
c
   10        continue
                xinterm = 0.5 * (xremote+xproxym)
                yinterm = 0.5 * (yremote+yproxym)
                zinterm = 0.5 * (zremote+zproxym)
                ninterm = kgeom(xinterm,yinterm,zinterm)
c
                if( ninterm.eq.nactual ) then
                   xproxym = xinterm
                   yproxym = yinterm
                   zproxym = zinterm
                else
                   xremote = xinterm
                   yremote = yinterm
                   zremote = zinterm
                   nremote = ninterm
                endif
c
                loopstep = loopstep+1
                dstrack=sqrt((xremote-xproxym)*(xremote-xproxym)
     a                      +(yremote-yproxym)*(yremote-yproxym)
     b                      +(zremote-zproxym)*(zremote-zproxym))
                if( abs(dstrack)-steplim ) 20,10,10
   20        continue
c
             if(loopstep.ge.64) then 
               write(*,*) 'convergence missing '
     a                                ,loopstep,dstrack
c
               write (*,*) nevent, nactual, nremote, x(np), y(np), z(np)
     a,                    xremote, yremote, zremote, dstrack, tracklim
             endif
             if (nremote.le.maxreg) then
                irnew = nremote
                ustep = sqrt ( ( xremote - x(np) ) * ( xremote - x(np) )
     a                      + ( yremote - y(np) ) * ( yremote - y(np) )
     b                      + ( zremote - z(np) ) * ( zremote - z(np) ))
                idisc = 0
cnico
		dnear(np) = 0.	
cnico
             else
	        idisc = 1000
	     endif
          endif
      endif
c mirko
c      write (*,*) 'out howfar',xremote,yremote,zremote,nremote,iq(np)
c     a, idisc
c
      return
c
      end
c
