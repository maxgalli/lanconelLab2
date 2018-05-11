      subroutine hownear(tperp)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
c
c
cale
c
      real*4 tperp
      real r_temp
cale
      if (ir(np).ge.maxreg) then
        tperp=0.0
      elseif (ir(np).gt.regoffset) then
cnico particle inside the collimator
       	if(modlabel.eq.1) then 
c
c     moduli esagonali
c       
      	   if(med(ir(np)).eq.indair) then
c		tperp = dmin1(dnear_temp,xair-dnear_temp)
		tperp = dnear_temp
      	   elseif(med(ir(np)).eq.indlead) then
		tperp = dmin1(dnear_temp,xlead-dnear_temp)
	   else
		write(*,*) 'Error in hownear'
                write (*,*)'med(np)',med(np)
                write (*,*)'iregtype',iregtype(ir(np))
                write (*,*)'ir(np)',ir(np)
                write (*,*)'x y z',x(np),y(np),z(np)
           endif
c
c    fine moduli esagonali       
c
       	elseif(modlabel.eq.2) then 
c
c     moduli rettangolari
c
      	   if(med(ir(np)).eq.indair) then
		tperp = dmin1(xtemp,ytemp,(xair-xtemp),(yair-ytemp))
      	   elseif(med(ir(np)).eq.indlead) then
		tperp = dmin1((xtemp-xair),(ytemp-yair),
     a		         (xcell-xtemp),(ycell-ytemp))
	   endif
c
c    fine moduli rettangolari       
c
	endif
      elseif (ir(np).gt.regoffset_det1mod) then
cnico particle inside the detector
        tperp = amin1( abs(x(np)-xinf(ir(np))) , abs(x(np)-xsup(ir(np)))
     a               , abs(y(np)-yinf(ir(np))) , abs(y(np)-ysup(ir(np)))
     b              , abs(z(np)-zinit(ir(np))), abs(z(np)-zend(ir(np))))
      elseif (iregtype(ir(np)).eq.1)then
        tperp = amin1( abs(x(np)-xinf(ir(np))) , abs(x(np)-xsup(ir(np)))
     a               , abs(y(np)-yinf(ir(np))) , abs(y(np)-ysup(ir(np)))
     b              , abs(z(np)-zinit(ir(np))), abs(z(np)-zend(ir(np))))
c
      elseif (iregtype(ir(np)).eq.2)then
           r_temp=radius(ir(np))-sqrt(((xsup(ir(np))-x(np))*
     a                    (xsup(ir(np))-x(np)))+((y(np)-yinf(ir(np)))*
     b       (y(np)-yinf(ir(np)))))
             tperp=amin1(abs(z(np)-zinit(ir(np))),
     a                  abs(z(np)-zend(ir(np))),r_temp)

      elseif (iregtype(ir(np)).eq.3)then
           r_temp=radius(ir(np))-sqrt(((zend(ir(np))-z(np))*
     a                   (zend(ir(np))-z(np)))+((y(np)-yinf(ir(np)))*
     b       (y(np)-yinf(ir(np)))))
             tperp=amin1(abs(x(np)-xinf(ir(np))),
     a                  abs(x(np)-xsup(ir(np))),r_temp)

cale      endif
      else
	write(*,*) 'Error in hownear ',ir(np),iregtype(ir(np))
 
      endif

      return
c
      end
