      subroutine rungeominit
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
c
c
      integer  k, j
      integer itmpidm,numslab,cont_cell,mod_tempx,mod_tempy
      real tempinf, tempsup
      doubleprecision xdet1,ydet1,xdet2,ydet2
      integer numcelly_det1mod,numcelly_det2mod
c
      read(ltydat,80,err=41,end=46) nplan
c
      if(nplan.ge.maxreg) then
         write(*,*) ' error in region number ' , nplan,' > ',maxreg
      endif
c
      do 10 k=1,nplan
         read(ltydat,100,err=41,end=46) indmat(k), zinit(k),zend(k), 
     + xinf(k),xsup(k),yinf(k),ysup(k),radius(k),index(k),iregtype(k)
c
      write(*,100) indmat(k),
     + zinit(k),zend(k),xinf(k),xsup(k),yinf(k),ysup(k),radius(k)
     + ,index(k),iregtype(k)
c
c ale iregtype,radius 09 04 2001 
c
   10 continue
c
      nreg = nplan
      indmat(maxreg+1) = 0
      xmax  = -10000.0
      xmin  = +10000.0
      ymax  = -10000.0
      ymin  = +10000.0
      zmin  = +10000.0
      zmax  = -10000.0
c
      numslab = 0
c
c------------------------------------------------------c
c
      write (*,*) ' '
      write (*,*) ' DETECTORS '
      write (*,*) ' '
c
      do 30 k=1,nplan
c
c---  sort the geometrical limits of the detector
c
      tempinf  = amin1(zinit(k),zend(k))
      tempsup  = amax1(zinit(k),zend(k))
      zinit(k) = tempinf
      zend(k)  = tempsup
      tempinf  = amin1(xinf(k),xsup(k))
      tempsup  = amax1(xinf(k),xsup(k))
      xinf(k)  = tempinf
      xsup(k)  = tempsup
      tempinf  = amin1(yinf(k),ysup(k))
      tempsup  = amax1(yinf(k),ysup(k))
      yinf(k)  = tempinf
      ysup(k)  = tempsup
c
c-------
c
      if ( (index(k)/10000).eq.labeldet) then 
c
         write(*,100) indmat(k), zinit(k), zend(k), xinf(k), xsup(k)
     +, yinf(k), ysup(k), index(k), iregtype(k)
c
         zdetend = zend(k)
         ydetinf = yinf(k)
         xdetinf = xinf(k)
         ydetsup = ysup(k)
         xdetsup = xsup(k)
c
      endif
c
c ale iregtype, raysqrt, radius 11 04 2001
c
      raysqrt(k)= radius(k)**2
   30 continue
c
      write(*,60) nplan
c
c-------------------------------------------------------c
c
      do 20 k=1,nplan
c
c---  set geometrical limits of apparatus
c---
c
      if(zinit(k).le.zmin)     zmin  = zinit(k)
      if(zend(k).ge.zmax)      zmax  = zend (k)
c
      if(xinf(k).le.xmin)      xmin  = xinf(k)
      if(xsup(k).ge.xmax)      xmax  = xsup(k)
      if(yinf(k).le.ymin)      ymin  = yinf(k)
      if(ysup(k).ge.ymax)      ymax  = ysup(k)
c
      med(k)=indmat(k)
      itmpidm = index(k)/10000
c
      write(*,70) k,indmat(k),(media(j,indmat(k)),j=1,24), med(k),
     +pcut(indmat(k)),ecut(indmat(k)), zinit(k),zend(k),xinf(k),xsup(k),
     +yinf(k),ysup(k), index(k),itmpidm
 
c
   20 continue
c
c----
c
      tracklim  = sqrt ((zmax-zmin)*(zmax-zmin)+
     a                   (xmax-xmin)*(xmax-xmin)+
     b                   (ymax-ymin)*(ymax-ymin))
      steplim = 1.e-04
 
      write(*,*) ' GEOM_LIMIT ',zmin,zmax,xmin,xmax,ymin,ymax
      write(*,*) ' TRACKLIMIT ',tracklim
      write(*,*) ' STEPLIM    ',steplim
      write(*,*) ' '
cc
cale
ccc   NICO aggiunge x geometria modulare detector
ccc 
       read(ltydat,95,err=41,end=46) modlabel,zinit_det1mod,
     a zend_det1mod,xinf_det1mod,xsup_det1mod,yinf_det1mod,ysup_det1mod
     b ,xcell_det1mod,ycell_det1mod,xair_det1mod,yair_det1mod,
     c indlead_det1mod,indair_det1mod
c
      xdet1 = xsup_det1mod - xinf_det1mod
      ydet1 = ysup_det1mod - yinf_det1mod
      xcellinv_det1mod = 1./xcell_det1mod
      ycellinv_det1mod = 1./ycell_det1mod
      numcellx_det1mod = dint (xdet1 / xcell_det1mod)
      numcelly_det1mod = dint (ydet1 / ycell_det1mod)
      numcell_det1mod = numcellx_det1mod * numcelly_det1mod
      regoffset_det1mod = nreg
cnico cambia x dead zone      nreg = nreg + numcell_det1mod
      nreg = nreg + numcell_det1mod * 3
      write(*,*) 'xdet1mod ',xcell_det1mod,xdet1/xcell_det1mod
      write(*,*) 'xdet1 ydet1 ncellxdet1 ncellydet1 regoffsdet1 nreg'
     a,xdet1,ydet1,numcellx_det1mod,numcelly_det1mod,regoffset_det1mod
     b,nreg
      cont_cell = 0
      do 11 k = regoffset_det1mod+1,regoffset_det1mod+numcell_det1mod
     	    med(k) = indair_det1mod
	    indmat(k) = indair_det1mod
	    index(k) = labeldet * 10000
cnico aggiunge x dead zone
     	    med(k+numcell_det1mod) = indlead_det1mod
	    indmat(k+numcell_det1mod) = indlead_det1mod
	    index(k+numcell_det1mod) = 0
     	    med(k+numcell_det1mod+numcell_det1mod) = indlead_det1mod
	    indmat(k+numcell_det1mod+numcell_det1mod) = indlead_det1mod
	    index(k+numcell_det1mod+numcell_det1mod) = 0
cnico
	    mod_tempx = mod(cont_cell,numcellx_det1mod)
	    mod_tempy = int(cont_cell/numcellx_det1mod)
	    xinf(k) = xinf_det1mod + mod_tempx * xcell_det1mod
	    xsup(k) = xinf(k) + xcell_det1mod
	    yinf(k) = yinf_det1mod + mod_tempy * ycell_det1mod
	    ysup(k) = yinf(k) + ycell_det1mod
	    zinit(k) = zinit_det1mod
	    zend(k) = zend_det1mod
	    cont_cell = cont_cell + 1
            if(zend(k).gt.zdetend) zdetend = zend(k)
            if(yinf(k).lt.ydetinf) ydetinf = yinf(k)
            if(xinf(k).lt.xdetinf) xdetinf = xinf(k)
            if(ysup(k).gt.ydetsup) ydetsup = ysup(k)
            if(xsup(k).gt.xdetsup) xdetsup = xsup(k)
   11 continue	    
ccc
ccc   NICO aggiunge x geometria modulare collimatore
ccc 
       read(ltydat,95,err=41,end=46) modlabel,zinitcoll,zendcoll,
     a xinfcoll,xsupcoll,yinfcoll,ysupcoll,xcell,ycell,xair,yair,
     b indlead,indair
c
       if(modlabel.eq.1) then 
c
c     moduli esagonali
c       
c
       septa = xcell
       xcoll = xsupcoll - xinfcoll
       ycoll = ysupcoll - yinfcoll
c       xair = sqrt(3.) * xair / 2
       xcell = 2. * (xair + septa / 2.)
       xcell_2 = xcell / 2.
       xlead = xcell_2 - xair
       xcellinv = 1./xcell_2
       ycell = xcell * 2. / sqrt(3.)
       numcellx = dint (xcoll / xcell_2) +1
       numcelly = dint (ycoll / (ycell/4) + 2)
       numcell = numcellx * numcelly
       regoffset = nreg 
       u3offset = (int((numcellx + 1)/2) ) * xcell_2
       u3celloffset = int((numcellx + 1)/2) + 1
       nreg = nreg + 2 * numcell
       write(*,*) 'zinitcoll zendcoll, xinfcoll xsupcoll, yinfcoll
     a ysupcoll xcell ycell xair yair indlead indair',zinitcoll,zendcoll
     b ,xinfcoll,xsupcoll,yinfcoll,ysupcoll,xcell,ycell,xair,yair,
     c indlead,indair
       write(*,*) 'xcoll ycoll numcellex numcelley regoffset nreg'
     a ,xcoll,ycoll,numcellx,numcelly,regoffset,nreg
       write(*,*) 'u3off u3celloff xcell_2', u3offset,
     a u3celloffset,xcell_2 
	if(mod(u3celloffset,2).eq.0) then
		label_hex = 0
	else
		label_hex = 1
	endif
	label_hex1 = mod(int((2*u3celloffset-1)/2),3)
       write(*,*) 'label_0 label_1',label_hex,label_hex1
       do 12 k = 1,numcellx
          do 12 j = 1,numcelly
		iu1 = k
		iu2 = int((j+k)/2)
		iu3 = int((j-k+2*u3celloffset-1)/2)
		if(label_hex.eq.0) then
		icell_temp = regoffset + (iu1-1)*numcelly + 
     a		2*(iu2-int((iu1+1)/2)) + mod(abs(iu3-iu2+1),2) + 
     b		2*mod(iu1,2)*mod(abs(iu3-iu2),2)
     		else
		icell_temp = regoffset + (iu1-1)*numcelly + 
     a		2*(iu2-int((iu1+1)/2)) + mod(abs(iu3-iu2),2) + 
     b		2*mod(iu1,2)*mod(abs(iu3-iu2+1),2)
		endif
		med(icell_temp) = indair
		indmat(icell_temp) = indair
		med(icell_temp+numcell) = indlead
		indmat(icell_temp+numcell) = indlead
		index(icell_temp) = 0
		index(icell_temp+numcell) = 0
   12  continue
c
c    fine moduli esagonali       
c
       elseif(modlabel.eq.2) then 
c
c     moduli rettangolari
c       
      xcoll = xsupcoll - xinfcoll
      ycoll = ysupcoll - yinfcoll
      xcellinv = 1./xcell
      ycellinv = 1./ycell
      numcellx = dint (xcoll / xcell)
      numcelly = dint (ycoll / ycell)
      numcell = numcellx * numcelly
      regoffset = nreg
      nreg = nreg + 3 * numcell
      do 13 k = regoffset+1,regoffset+numcell
     	    med(k) = indair
	    indmat(k) = indair
	    med(k+numcell) = indlead
	    indmat(k+numcell) = indlead
	    med(k+numcell+numcell) = indlead
	    indmat(k+numcell+numcell) = indlead
   13 continue	    
       write(*,*) 'zinitcoll zendcoll, xinfcoll xsupcoll, yinfcoll
     a ysupcoll xcell ycell xair yair indlead indair',zinitcoll,zendcoll
     b ,xinfcoll,xsupcoll,yinfcoll,ysupcoll,xcell,ycell,xair,yair,
     c indlead,indair
       write(*,*) 'xcoll ycoll numcellex numcelley regoffset nreg'
     a ,xcoll,ycoll,numcellx,numcelly,regoffset,nreg
       write(*,*) 'u3off u3celloff xcell_2', u3offset,
     a u3celloffset,xcell_2 
c
c    fine moduli rettangolari       
c
       elseif(modlabel.eq.0) then 
         regoffset = nreg
       write(*,*) 'zinitcoll zendcoll, xinfcoll xsupcoll, yinfcoll
     a ysupcoll xcell ycell xair yair indlead indair',zinitcoll,zendcoll
     b ,xinfcoll,xsupcoll,yinfcoll,ysupcoll,xcell,ycell,xair,yair,
     c indlead,indair
       write(*,*) 'xcoll ycoll numcellex numcelley regoffset nreg'
     a ,xcoll,ycoll,numcellx,numcelly,regoffset,nreg
       write(*,*) 'u3off u3celloff xcell_2', u3offset,
     a u3celloffset,xcell_2 
       endif
ccc
ccc  fine aggiunta x geometria modulare collimatore
ccc	
cale
c------------------------------------------------------c
c
      return
c
c------------------------------------------------------c
c
   41 write(*,*) ' ERROR IN GEOMETRY FILE - RUNGEOMINIT '
   46 write(*,*) ' END IN GEOMETRY FILE - RUNGEOMINIT '
      close(ltydat)
c
      return
c
c-----------------------------------------------------c
c
   50 FORMAT(1X,I9,1X,6F10.3)
   60 FORMAT(/,2X,' NUMERO PIANI ',I10,/)
   70 FORMAT( I5,' INDMED =',I3,', MED =',2X,24A1,2X,I3,/,5X
     +,' CUTS _ P & E ',2X,E16.4,2X,E16.4,/, 5X ,' COUNT SIZE  '
     +, 6(1X,F 10.3),/,5X,' COUNT INDEX  ', I10 , ' INDX/10000 ', I10,/)
   80 FORMAT(I10)
   90 FORMAT(I10,6F10.0,I10)
   95 FORMAT(I2,10F10.4,2I2)
c
c ale iregtype radius format 09 04 2001
c
  100 FORMAT(I10,7F10.4,I10,I10)
c
c ale format iregtype radius09 04 2001
c
      END
c
