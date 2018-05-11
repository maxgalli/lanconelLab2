      subroutine tuplesetup
c
c---------------------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
c
      real egaussian, sgaussian_x, sgaussian_y, tuptempvarx, tuptempvary
c
c---------------------------------------------------------------------c
c
      egaussian = 0.0
      call norran ( egaussian)
c
      sgaussian_x = 0.0
      call norran ( sgaussian_x)
c
      sgaussian_y = 0.0
      call norran ( sgaussian_y)
c
      xtup3(1) = nactsx
      xtup3(2) = nevent
c
      xtup3(3)  = xin
      xtup3(4)  = yin
      xtup3(5)  = zin
      xtup3(6)  = uin
      xtup3(7)  = vin
      xtup3(8)  = win
c
      xavg=xavg/eavg
      yavg=yavg/eavg
      zavg=zavg/eavg          
c
      xtup3( 9)  =  eavg
      xtup3(10)  =  xavg
      xtup3(11)  =  yavg
      xtup3(12)  =  zavg
cnico
cc      write(outtxt,100) xavg,yavg
cnico
c
cc      xtup3(13)  =  ecountmax
cc      xtup3(14)  =  xcountmax
cc      xtup3(15)  =  ycountmax
cc      xtup3(16)  =  zcountmax
c
cc      xtup3(17)  =  efirst
cc      xtup3(18)  =  xfirst
cc      xtup3(19)  =  yfirst
cc      xtup3(20)  =  zfirst
c
      xtup3(13)  = amax1( 0. , eavg + spren * egaussian)
c
      if (ispres.eq.1) then
         tuptempvarx = xavg + sgaussian_x * (zdetend - zavg) / 3
         tuptempvary = yavg + sgaussian_y * (zdetend - zavg) / 3
         if ((tuptempvarx.ge.xdetinf).and.(tuptempvarx.le.xdetsup).and.
     a      (tuptempvary.ge.ydetinf).and.(tuptempvary.le.ydetsup)) then
                xtup3(14)  = tuptempvarx
                xtup3(15)  = tuptempvary
         else
                xtup3(14)  = xdetinf - 1.
                xtup3(15)  = ydetinf - 1.
         endif
      else
         xtup3(14)  =  xavg
         xtup3(15)  =  yavg
      endif
c
cc      xtup3(24)  = amax1( 0. , ecountmax + spren * egaussian)
c
      if (ispres.eq.1) then
         tuptempvarx = xcountmax + sgaussian_x * (zdetend - zcountmax)/3
         tuptempvary = ycountmax + sgaussian_y * (zdetend - zcountmax)/3
         if ((tuptempvarx.ge.xdetinf).and.(tuptempvarx.le.xdetsup).and.
     a      (tuptempvary.ge.ydetinf).and.(tuptempvary.le.ydetsup)) then
cc                xtup3(25)  = tuptempvarx
cc                xtup3(26)  = tuptempvary
         else
cc                xtup3(25)  = xdetinf - 1.
cc                xtup3(26)  = ydetinf - 1.
         endif
      else
cc         xtup3(25)  =  xcountmax
cc         xtup3(26)  =  ycountmax
      endif
c
cc      xtup3(27)  = amax1( 0. , efirst + spren * egaussian)
c
      if (ispres.eq.1) then
         tuptempvarx =  xfirst + sgaussian_x * (zdetend - zfirst) / 3
         tuptempvary =  yfirst + sgaussian_y * (zdetend - zfirst) / 3
         if ((tuptempvarx.ge.xdetinf).and.(tuptempvarx.le.xdetsup).and.
     a      (tuptempvary.ge.ydetinf).and.(tuptempvary.le.ydetsup)) then
cc                xtup3(28)  = tuptempvarx
cc                xtup3(29)  = tuptempvary
         else
cc                xtup3(28)  = xdetinf - 1.
cc                xtup3(29)  = ydetinf - 1.
         endif
      else
cc         xtup3(28)  =  xfirst
cc         xtup3(29)  =  yfirst
      endif
c
      xpixel=xpixel/epixel
      ypixel=ypixel/epixel
      zpixel=zpixel/epixel
c
      xtup3(16)  =  epixel
      xtup3(17)  =  amax1( 0. , epixel + spren * egaussian)
      xtup3(18)  =  xpixel
      xtup3(19)  =  ypixel
c      xtup3(34)  =  zpixel
c
cc      xtup3(34)  =  epixmax
cc      xtup3(35)  =  amax1( 0. , epixmax + spren * egaussian)
cc      xtup3(36)  =  xpixmax
cc      xtup3(37)  =  ypixmax
c      xtup3(37)  =  zpixmax          
c
cc      xtup3(38)  = xavg - (xin+uin*(zavg-zin)/win)
cc      xtup3(39)  = yavg - (yin+vin*(zavg-zin)/win)
c
cc      xtup3(40)  = multiplicity
c
      call hfn  (nt_hist,xtup3)
c
      return
c
  100 FORMAT(2(F8.4,1X))
      END
c
