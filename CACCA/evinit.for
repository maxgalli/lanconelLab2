      subroutine evinit
c
      implicit none
      include 'egs4comm.for'
c
      integer mh, nt
c
      first_int_label = .true.
c
      nhit         = 0
      elost        = 0.0
      biosource    = 0.0
      multiplicity = 0       
c
c---  set o zero energy loss in all regions
c
cnico nico cambia nplan con regoffset
c      do 10 mh=1,nplan
      do 10 mh=1,regoffset
            deodx(mh)=0.d0
   10 continue
c
      do 20 nt = 1,size_of_ntup
            xtup3(nt) =  0.0
   20 continue
c
      eavg=0.0
      xavg=0.0
      yavg=0.0
      zavg=0.0
c
      ecountmax=0.0
      xcountmax=0.0
      ycountmax=0.0
      zcountmax=0.0
c
      efirst=0.0
      xfirst=0.0
      yfirst=0.0
      zfirst=0.0
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
      return
c
      end
c
