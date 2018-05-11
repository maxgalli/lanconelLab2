
c
      data totim,tleft,tmed,tlim,tevent / 5*0. /
c
c-------------------------------------------------------------------c
c---------  i.o unit definition   ----------------------------------c
c
c      ltyin     =  standard  input       unit = 15
c      ltydat    =  geometry  file        unit = 30
c      ltyfil    =  histogram file        unit = 20
c      kmpi      =  material   for hatch  unit = 12
c      hislun    =  hbook  output         unit = 40
c      lxdata    =  binary output         unit = 80
c      kmpo      =  echo unit for hatch   unit =  8
c
c-------------------------------------------------------------------c
c-------------------------------------------------------------------c
c
      data ltyin, ltyfil, ltydat, hislun, lxdata / 15, 20, 30, 40, 80  /
c
c------------------------------------------------------------------c
c---  set total kinetic energy  -----------------------------------c
c    
      data ichar, labeldet, labelbiomass, enkin, evkin, ncases, itctx
     a, trigger, treshold, spren, esnr, maxreg, nreg , ispres, zdetend
     b, ydetinf, ydetsup, xdetinf, xdetsup / 3*0 , 2*0.d0 , 2*0 
     c, 4*0.0,  MAXR , 2*0 , 5*0.0 /
c
      data xmin, xmax, ymin, ymax, zmin, zmax, maxhit, tracklim, steplim
     a      / 6*0. , 99 , 1.e30 , 1.e-4 /
c
      data zinitcoll,zendcoll,xinfcoll,xsupcoll,yinfcoll,ysupcoll
     a,     xcoll,ycoll,xcell,ycell,xair,yair,xcellinv,ycellinv
     b,     numcellx,numcelly,numcell,regoffset / 14*0.d0, 4*0 /
c
      data zinit_det1mod,zend_det1mod,xinf_det1mod,xsup_det1mod,
     a     yinf_det1mod,ysup_det1mod,xcell_det1mod,ycell_det1mod,
     b     zinit_det2mod,zend_det2mod,xinf_det2mod,xsup_det2mod,
     c     yinf_det2mod,ysup_det2mod,xcell_det2mod,ycell_det2mod
     d / 16*0.d0 /
c
      data  isxnum, isxtype, isxreg, isxev, sx_cent, sy_cent, sz_cent
     a, rsx, rsxsq, rsxsqinv, rsy, rsysq, rsysqinv, rsz, rszsq, rszsqinv
     b, isxangtype, isxgeomtype / 0 , 60*0 , 240*0. , 40*0 /
      data dnear_temp / 0./
      data  thetaangsx, phiangsx, cthinf, cthsup / 80*0. /
      data  sx_inf, sx_sup, sy_inf, sy_sup, sz_inf, sz_sup / 120*0. /
c
c---------------------------------------------------------------------c
c data per histin
c
      data size_of_ntup / 19  /
c
c--
c
      data nt_hist,nt_cycle, num_of_ntup    /
     a         900 ,        1  ,  10000     /
c
c-- run data
c
      data nt_tags( 1) /'sxnum'/
      data nt_tags( 2) /'history'/
c
c-- source accepted point
c
      data nt_tags( 3) /'xorig'/
      data nt_tags( 4) /'yorig'/
      data nt_tags( 5) /'zorig'/
c
c-- zone counter
c
      data nt_tags( 6) /'a1orig'/
      data nt_tags( 7) /'a2orig'/
      data nt_tags( 8) /'a3orig'/
c
c-- physical total energy released ESLAB 
c
      data nt_tags(  9) /'eslab'/
      data nt_tags( 10) /'xslab'/
      data nt_tags( 11) /'yslab'/
      data nt_tags( 12) /'zslab'/
c
c-- physical total energy released with ECMAX 
c
cc      data nt_tags( 13) /'ecmax'/
cc      data nt_tags( 14) /'xcmax'/
cc      data nt_tags( 15) /'ycmax'/
cc      data nt_tags( 16) /'zcmax'/
c
c-- physical total energy released at first interaction in detector
c
cc      data nt_tags( 17) /'efirst'/
cc      data nt_tags( 18) /'xfirst'/
cc      data nt_tags( 19) /'yfirst'/
cc      data nt_tags( 20) /'zfirst'/
c
c-- detector response for crystal ESLAB ( energy and spatial resolution )
c
      data nt_tags( 13) /'eslabesr'/
      data nt_tags( 14) /'xslabesr'/
      data nt_tags( 15) /'yslabesr'/
c
c-- detector response for crystal ECMAX ( energy and spatial resolution ) 
c
cc      data nt_tags( 24) /'ecmaxesr'/
cc      data nt_tags( 25) /'xcmaxesr'/
cc      data nt_tags( 26) /'ycmaxesr'/
c
c-- detector response for crystal EFIRST ( energy and spatial resolution ) 
c
cc      data nt_tags( 27) /'efirsesr'/
cc      data nt_tags( 28) /'xfirsesr'/
cc      data nt_tags( 29) /'yfirsesr'/
c
c-- total energy released with EPIXEL 
c
      data nt_tags( 16) /'epixel'/
      data nt_tags( 17) /'erpixel'/
      data nt_tags( 18) /'xpixel'/
      data nt_tags( 19) /'ypixel'/
c
c-- total energy released with EPIXMAX 
c
cc      data nt_tags( 34) /'epixmax'/
cc      data nt_tags( 35) /'erpixmax'/
cc      data nt_tags( 36) /'xpixmax'/
cc      data nt_tags( 37) /'ypixmax'/
c
c-- BIOMASS and OTHER
c
cc      data nt_tags( 38) /'radsqx'/
cc      data nt_tags( 39) /'radsqy'/
cc      data nt_tags( 40) /'multipl'/
c
c-----------------------------------------------------------c
c
cc      data totup3,xtup3 / 40,  40 *0. /
      data totup3,xtup3 / 19,  19 *0. /
c
      data xray_memdir / 'GAMMARAY'/
c
      data histo_count,plot_count / 2 * 0 /
      data id_hist,bd_hist/ 40*1 /
c
      data maxmem /  3999000 /
c
c---------------------------------------------------------------------c
c
      data  iloop, nevent, ilsx, nactsx / 4 * 0 /
c
c---------------------------------------------------------------------c
c
      data nhit, elost, multiplicity / 0 , 0.0 , 0 /
      data eavg, xavg, yavg, zavg / 4*0.0 /
      data ecountmax, xcountmax, ycountmax, zcountmax / 4*0.0 /
      data efirst, xfirst, yfirst, zfirst / 4*0.0 /
      data epixel, xpixel, ypixel, zpixel / 4*0.0 /
      data epixmax, xpixmax, ypixmax, zpixmax / 4*0.0 /
      data hiteng,indexcount,indexord / 100*0.00, 200*0 /
      data first_int_label / .true. /
c
c---------------------------
c
      data deodx / MAXR1*0.d0 /
c
c--------------------------
c
      data temp_eng  /100 * 0.0  /
c
c---------------------------
c
      data totein, toteoob, totebio, totedet, toteew, totedeodx, totemod  
     a / 7 * 0.d0  /
c
c---------------------------
