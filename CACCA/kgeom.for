       function kgeom(xpos,ypos,zpos)
c-------
c-----------------------------------------------------------c
c
      implicit none

       include 'egs4comm.for'
c
c
      integer k,iu1cell,iu2cell,iu3cell
      integer iu1temp,iu2temp,iu3temp,ind1temp
      integer ind2temp,ind3temp,ind4temp,indcell_temp
      real xpos,ypos,zpos,u2,u3
      real xpostemp,ypostemp,zpostemp
      doubleprecision xrel,yrel,xtemp0,ytemp0
      integer rawcell,colcell,cell
c
c-----------------------------------------------------------c
c
c          test for system  boundaries
c
c          cycle testing slabs
c
c           write (*,*) 'entro kgeom'
	   kgeom  = maxreg + 1 
c
c-------------- test if the position is inside the apparatus
c
           if ( (zmin.gt.zpos ) .or. (zmax.lt.zpos ) .or.
     a          (xmin.gt.xpos ) .or. (xmax.lt.xpos ) .or.
     b          (ymin.gt.ypos ) .or. (ymax.lt.ypos ) )then
c
                 kgeom  = maxreg +1
c           write (*,*) 'out of limits'

c
                 return
c
           endif
c	   
c--------------
c

cnico
	   dnear_temp = 0.
cnico	   
             if (zinitcoll.lt.zpos ) then
	        if (zendcoll.ge.zpos) then
	           if (xinfcoll.lt.xpos) then
	              if (xsupcoll.gt.xpos) then
	                 if (yinfcoll.lt.ypos) then
	                    if (ysupcoll.gt.ypos) then
ccc
ccc	particle inside the collimator
ccc
              if(modlabel.eq.1) then 		
c
c       moduli esagonali
c
		xrel = xpos - xinfcoll
		yrel = ypos - yinfcoll
		xtemp = 0.5 * xrel 
		ytemp = 0.8660254 * yrel
		u2 = xtemp + ytemp
		u3 = -xtemp + ytemp + u3offset
c	remind: xcellinv = 1/xcell_2
		iu1cell = dint(xrel*xcellinv)
		iu2cell = dint(u2*xcellinv)
		iu3cell = dint(u3*xcellinv) 
		iu1temp = mod(iu1cell,3)
		iu2temp = mod(iu2cell,3)
		iu3temp = mod(iu3cell,3)
		ind1temp = iu1cell + 1
		ind2temp = iu3cell - iu2cell
		ind3temp = abs(ind2temp)
		ind4temp = abs(ind2temp+1)
		if(label_hex.eq.0) then
		indcell_temp = regoffset + iu1cell * numcelly 
     a	        + 2*((iu2cell+1)-int((ind1temp+1)*.5)) + mod(ind4temp,2)
     c		+ 2*mod(ind1temp,2)*mod(ind3temp,2)
     		else
		indcell_temp = regoffset + iu1cell * numcelly 
     a	        + 2*((iu2cell+1)-int((ind1temp+1)*.5)) + mod(ind3temp,2)
     c		+ 2*mod(ind1temp,2)*mod(ind4temp,2)
		endif
c
	if(label_hex1.eq.0) then
		if(iu1temp.eq.0) then
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.0) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.0) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.1) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		elseif(iu1temp.eq.1) then     
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.1) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.0) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
			  return
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.0) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
			  return
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		else
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.0) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.1) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.0) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		endif
c     fine ciclo label_hex1 = 0 
	else if(label_hex1.eq.1) then 
		if(iu1temp.eq.0) then
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.1) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.1) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.2) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		elseif(iu1temp.eq.1) then     
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.2) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.1) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
			  return
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.1) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
			  return
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		else
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.1) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.2) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.2) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.1) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		endif
c     fine ciclo label_hex1 = 1 
	else
		if(iu1temp.eq.0) then
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.2) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.2) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.0) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		elseif(iu1temp.eq.1) then     
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.0) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.2) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
			  return
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.2) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
			  return
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		else
		   if(iu2temp.eq.0) then
		      if(iu3temp.eq.2) then
		        if(u3-iu3cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xair
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u3+iu3cell*xcell_2
c		write(*,*) 'A x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.0) then
		        if(u3-iu3cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2-xlead
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u3-iu3cell*xcell_2
c		write(*,*) 'D x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   elseif(iu2temp.eq.1) then
		      if(iu3temp.eq.0) then
		        if(xrel-iu1cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xair
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-xrel+iu1cell*xcell_2
c		write(*,*) 'C x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(u2-iu2cell*xcell_2.ge.xair) then
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xair
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom =  indcell_temp
			    dnear_temp = xair-u2+iu2cell*xcell_2
c		write(*,*) 'B x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   else
		      if(iu3temp.eq.2) then
		        if(xrel-iu1cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2-xlead
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = xrel-iu1cell*xcell_2
c		write(*,*) 'F x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      elseif(iu3temp.eq.1) then
		        if(u2-iu2cell*xcell_2.ge.xlead) then
			    kgeom =  indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2-xlead
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return
			else
			    kgeom = numcell + indcell_temp
			    dnear_temp = u2-iu2cell*xcell_2
c		write(*,*) 'E x y z u2 u3 iu1cell iu2cell iu3cell xtemp 
c     a	ytemp kgeom med' , xpos,ypos,zpos,
c     b	u2,u3,iu1cell,iu2cell,iu3cell,xtemp,ytemp,kgeom ,med(kgeom)
			  return			  
			endif  
		      endif
		   endif
		endif
	endif 
c     fine ciclo label_hex1 = 2 

c
c    fine moduli esagonali       
c
              elseif(modlabel.eq.2) then 
c
c     moduli rettangolari
c       
		xrel = xpos - xinfcoll
		yrel = ypos - yinfcoll
		xtemp0 = xrel * xcellinv
		ytemp0 = yrel * ycellinv
		rawcell = dint (xtemp0)
		colcell = dint (ytemp0)
		xtemp = xrel - rawcell * xcell
		ytemp = yrel - colcell * ycell
		cell = rawcell + 1 + colcell * numcellx
		if   (xtemp.gt.xair) then
			kgeom = regoffset + cell + numcell
			return
		elseif (ytemp.gt.yair) then
			kgeom = regoffset + cell + numcell + numcell
			return
		endif
		kgeom = regoffset + cell
		return
c
c    fine moduli rettangolari       
c
	      endif
                            endif
                         endif
                      endif
                   endif
                endif
             endif
ccc
ccc	particle outside the collimator
ccc
             if (zinit_det1mod.lt.zpos ) then
	        if (zend_det1mod.ge.zpos) then
	           if (xinf_det1mod.lt.xpos) then
	              if (xsup_det1mod.gt.xpos) then
	                 if (yinf_det1mod.lt.ypos) then
	                    if (ysup_det1mod.gt.ypos) then
ccc
ccc	particle inside the detector 
ccc
			xrel = xpos - xinf_det1mod
			yrel = ypos - yinf_det1mod
			xtemp0 = xrel * xcellinv_det1mod
			ytemp0 = yrel * ycellinv_det1mod
			rawcell = dint (xtemp0)
			colcell = dint (ytemp0)
			xtemp = xrel - rawcell * xcell_det1mod
			ytemp = yrel - colcell * ycell_det1mod
			cell = rawcell + 1 + colcell * numcellx_det1mod
cnico aggiunge x dead zone
		if   (xtemp.gt.xair_det1mod) then
			kgeom = regoffset_det1mod + cell + numcell_det1mod
			return
		elseif (ytemp.gt.yair_det1mod) then
			kgeom = regoffset_det1mod + cell + numcell_det1mod 
     a                          + numcell_det1mod
			return
		endif
cnico
			kgeom = regoffset_det1mod + cell
			return
                            endif
                         endif
                      endif
                   endif
                endif
             endif
ccc
ccc	particle outside the detector
ccc

             k = 0
c
   10        k=k+1
c
c-----     test for zpos
c
cale modifica per regioni cilindriche 09 04 2001
c              write (*,*) 'iregtype', iregtype(k)
      if(iregtype(k).eq.1)then
c                              write (*,*) 'inside parall1'
c      
cale: regioni parallelepipedi  
             if (zinit(k).lt.zpos ) then
	        if (zend(k).ge.zpos) then
	           if (xinf(k).lt.xpos) then
	              if (xsup(k).ge.xpos) then
	                 if (yinf(k).lt.ypos) then
	                    if (ysup(k).ge.ypos) then
c-----
c
c   this is the actual volume of detector space
c   and exit the slabs loop
c
                              kgeom = k
c                   
c                              write (*,*) 'inside parall2'
                              return
c
                            endif
                         endif
                      endif
                   endif
                endif
             endif
c
      elseif(iregtype(k).eq.2)then
c b
cale: regioni cilindriche parallele asse z
c123456 
       if(zinit(k).lt.zpos)then
       if(zend(k).ge.zpos)then
       if((((xsup(k)-xpos)*(xsup(k)-xpos))+((ypos-yinf(k))*       
     a (ypos-yinf(k)))).le.raysqrt(k))then
c        write(*,*) ((xpos-xsup(k))**2)+(ypos-yinf(k))**2
c         write(*,*) raysqrt(k)
                                                  kgeom = k
c
                                                  return
                                            endif
                                        endif
                                  endif
      else
c c
cale: regioni cilindriche parallele asse x

      if(xinf(k).lt.xpos)then
      if(xsup(k).ge.xpos)then
      if((((zpos-zend(k))*(zpos-zend(k)))+(ypos-yinf(k))*
     a (ypos-yinf(k))).le.raysqrt(k))then
                                                  kgeom = k
c
                                                  return
                                            endif
                                        endif
                                  endif
         endif
cale: end if iniziale
c
c       try next slab
c-------
c
          if( nplan - k ) 20,20,10
c
   20     continue
ccc	else
c
c     no region found : reject immediate
c
ccc          kgeom   =  maxreg +1
c
ccc          return
c
       end
