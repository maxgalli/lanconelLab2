      subroutine evsum
c
c---------------------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
c
      integer mh
c
c----- energy loss in all regions
c
      nhit         = 0
      multiplicity = 0
c
      if (biosource.gt.0.) call hfill(123,biosource,0.,1.)
      if (elost.gt.0.) call hfill(124,elost,0.,1.)
c
      epixel=0.0
      xpixel=0.0
      ypixel=0.0
      zpixel=0.0
c
      epixmax=0.0
      xpixmax=0.0
      ypixmax=0.0
      zpixmax=0.0
c
cnico nico cambia nplan con regoffset
c      do 10 mh=1,nplan
      do 10 mh=1,regoffset
c
         if ((index(mh)/10000).eq.labeldet) then
c
             if ((nhit.le.maxhit).and.(deodx(mh).gt.trigger)) then
c
                 nhit  = nhit + 1
                 hiteng(nhit)     = deodx(mh)
                 indexcount(nhit) = index(mh)
c
             endif
c
             if ((deodx(mh).gt.esnr).and.(nhit.le.maxhit)) then 
c
		 multiplicity = multiplicity +1
                 epixel=epixel+deodx(mh)
                 xpixel=xpixel+deodx(mh)*(xinf(mh)+xsup(mh))/2.
                 ypixel=ypixel+deodx(mh)*(yinf(mh)+ysup(mh))/2.
                 zpixel=zpixel+deodx(mh)*(zinit(mh)+zend(mh))/2.
c
                 if (deodx(mh).ge.epixmax) then
c
                    epixmax = deodx(mh)
                    xpixmax = (xinf(mh)+xsup(mh))/2.
                    ypixmax = (yinf(mh)+ysup(mh))/2.
                    zpixmax = (zinit(mh)+zend(mh))/2.
c
                 endif
c
             endif
c
         endif
c
 10   continue
c
c---------------------------------------------------------------------c
c
      if ( (elost.gt.treshold).and.(nhit.ge.1) ) then
c
           ilsx = ilsx +1
c
           call tuplesetup
           call hitorder
c
      endif
c
      return
c
      end
c
