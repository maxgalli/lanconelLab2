cale
      integer MAXR,MAXR1
      parameter(MAXR=120000,MAXR1=120001)
cale

      real totim,tleft,tmed,tlim,tevent 
c
      common / iounit / ltyin, ltydat, ltyfil, hislun, lxdata, outtxt
      integer ltyin, ltydat, ltyfil, hislun, lxdata, outtxt
c
      common / runpar / enkin, trigger, treshold, esnr, evkin, spren
     a, zdetend, ydetinf, ydetsup, xdetinf, xdetsup, ichar, labeldet
     b, labelbiomass, ncases, itctx, maxreg, nreg, ispres

      integer ichar, labeldet, labelbiomass, ncases, itctx
      integer maxreg, nreg, ispres
      doubleprecision enkin, trigger, treshold, esnr, evkin
      real spren, zdetend, ydetinf, ydetsup, xdetinf, xdetsup
c
      common /plngeom/ nplan, indmat(MAXR1), index(MAXR1), zinit(MAXR1)
     a, zend(MAXR1), xinf(MAXR1), xsup(MAXR1), yinf(MAXR1),ysup(MAXR1)
     b, iregtype(MAXR1), raysqrt(MAXR1), radius(MAXR1)
      integer nplan, indmat, index, iregtype
      real zinit, zend, xinf, xsup, yinf, ysup, raysqrt, radius
c
      common / limgeom / xmin, xmax, ymin, ymax, zmin, zmax
     a, maxhit, tracklim, steplim
      doubleprecision xmin, xmax, ymin, ymax, zmin, zmax
      real tracklim, steplim
      integer maxhit
c
      common /newgeom  / zinitcoll,zendcoll,xinfcoll,xsupcoll,yinfcoll,
     a  ysupcoll,xcell,ycell,xair,yair,xcellinv,ycellinv,dnear_temp,
     b	xcell_2,xlead,septa,zinit_det1mod,zend_det1mod,xinf_det1mod,
     c	xsup_det1mod,yinf_det1mod,ysup_det1mod,xcell_det1mod,
     d  ycell_det1mod,xair_det1mod,yair_det1mod,xcellinv_det1mod,
     e  ycellinv_det1mod,xtemp,ytemp,indlead_det1mod,indair_det1mod,
     f  zinit_det2mod,zend_det2mod,xinf_det2mod,xsup_det2mod,
     g  yinf_det2mod,ysup_det2mod,xcell_det2mod,ycell_det2mod,
     h  xair_det2mod,yair_det2mod,xcellinv_det2mod,u3offset,
     i  ycellinv_det2mod,indlead_det2mod,indair_det2mod,
     l  regoffset_det2mod,numcell_det2mod,numcellx_det2mod,
     m  numcellx,numcelly,numcell,regoffset,label_hex,
     n  label_hex1,u3celloffset,iu1,iu2,iu3,icell_temp,indair,indlead,       
     o	regoffset_det1mod,numcell_det1mod,numcellx_det1mod,modlabel
      doubleprecision    zinitcoll,zendcoll,xinfcoll,xsupcoll, 
     a       yinfcoll,ysupcoll,xcell,ycell,xair,yair,xcellinv,
     b	     ycellinv,xcell_2,u3offset,xlead,xcoll,ycoll,septa,
     c       dnear_temp,zinit_det1mod,zend_det1mod,xinf_det1mod,
     d       xsup_det1mod,yinf_det1mod,ysup_det1mod,ycell_det1mod,
     e       xair_det1mod,yair_det1mod,xcellinv_det1mod,
     f       ycellinv_det1mod,xcell_det1mod,xtemp,ytemp,
     g       zinit_det2mod,zend_det2mod,xinf_det2mod,xsup_det2mod,
     h       yinf_det2mod,ysup_det2mod,ycell_det2mod,xair_det2mod,
     i       yair_det2mod,xcellinv_det2mod,xcell_det2mod,
     l       ycellinv_det2mod
      integer            numcellx,numcelly,numcell,regoffset,
     a	     label_hex,label_hex1,u3celloffset,icell_temp,iu1,iu2,
     b       iu3,indair,indlead,indlead_det1mod,indair_det1mod,
     c       regoffset_det1mod,numcell_det1mod,numcellx_det1mod,modlabel
     d      ,indlead_det2mod,indair_det2mod,regoffset_det2mod,
     e       numcell_det2mod,numcellx_det2mod
c
      common / sxgeom /  rsx(20),rsxsq(20),rsxsqinv(20),rsy(20),rsz(20)
     a, rsysq(20),rsysqinv(20),rszsq(20),rszsqinv(20),esource(20)
     b, sx_inf(20),sx_sup(20),sy_inf(20),sy_sup(20),sz_inf(20)
     c, sz_sup(20),thetaangsx(20),phiangsx(20),cthsup(20),cthinf(20)
     d, sx_cent(20),sy_cent(20),sz_cent(20),isxangtype(20)
     e, isxgeomtype(20),isxtype(20), isxreg(20), isxev(20), isxnum 
      integer isxnum, isxtype, isxreg, isxev, isxangtype, isxgeomtype
      real sx_cent, sy_cent, sz_cent
      doubleprecision rsx, rsxsq, rsxsqinv
      doubleprecision rsy, rsysq, rsysqinv, rsz, rszsq, rszsqinv, sx_inf
      doubleprecision sx_sup, sy_inf, sy_sup, sz_inf, sz_sup
      doubleprecision thetaangsx, phiangsx, cthsup, cthinf
      doubleprecision esource
c
      integer kgeom
c
      common /sxdata/ein, xin, yin, zin, uin, vin, win, wtin, irin, iqin
      integer iqin, irin
      doubleprecision ein
      real xin, yin, zin, uin, vin, win, wtin
      real eireal
c
c---------------------------------------------------------------------c
c                                                                     c
c     histogramming program histogramming  definitions                c
c                                                                     c
c---  blank common for hbook4 + paw  ---------------------------------c
c
      common /  pawc  /  hmemor(4100000)
      real               hmemor
c
c
      common /mem_parameters/ maxmem,histo_count,plot_count,
     a                        num_of_ntup,size_of_ntup
      integer maxmem,histo_count,plot_count,num_of_ntup,size_of_ntup
c
      common/x_strings/xray_filnam,xray_memdir
      character*64   xray_filnam
      character*64   xray_memdir
c
      common /ntuples_tag/ totup3, xtup3(19), nt_tags(19)
      character*10         nt_tags
      integer              totup3
      real                 xtup3
c
      common/h_idents/id_hist(20),id_cycle(20),hist_ident(20)
     a,               bd_hist(20),bd_cycle(20),bdhs_ident(20)
     b,               nt_hist,nt_cycle,nt_ident
cc
      integer id_hist,id_cycle,hist_ident
      integer bd_hist,bd_cycle,bdhs_ident
      integer nt_hist,nt_cycle,nt_ident
cc
      common / altro / istat
      integer  istat
      common / hist_fil_nam / histofilename,pathname,pathlength
      character*80 histofilename
      character*80 pathname
      integer pathlength
cale
c---------------------------------------------------------------------c
c
      common / rundata / iloop, nevent, ilsx, nactsx
      integer iloop, nevent, ilsx, nactsx
c
      common  /dethits / nhit,multiplicity,elost,biosource
     a,                  eavg,xavg,yavg,zavg
     b,                  ecountmax,xcountmax,ycountmax,zcountmax
     c,                  efirst, xfirst, yfirst, zfirst
     d,                  epixel,xpixel,ypixel,zpixel
     e,                  epixmax,xpixmax,ypixmax,zpixmax
     f,                  first_int_label
     g,                  hiteng(100),indexcount(100),indexord(100)
      integer nhit, multiplicity, indexord, indexcount
      real elost, biosource, eavg, xavg, yavg, zavg
      real ecountmax, xcountmax, ycountmax, zcountmax
      real efirst, xfirst, yfirst, zfirst
      real epixel, xpixel, ypixel, zpixel
      real epixmax, xpixmax, ypixmax, zpixmax, hiteng
      logical first_int_label
c
c--------------------------------------------------------------------c
c
      common /profile/ deodx(MAXR1)
      doubleprecision deodx
c
c--------------------------------------------------------------------c
c
      common/ gensxpl/ cth,sth,phibis,thetarandom,dnorm,xi0,yi0,zi0
      real xi0,yi0,zi0
      doubleprecision cth,sth,phibis,thetarandom,dnorm
c
c---------------------------------------------------------------------c
c
      common / temp_e / temp_eng(100)
      doubleprecision temp_eng
c
c----------------------------
c
      common / tot_ene / totein, toteoob, totebio, totedet, toteew
     a, totedeodx, totemod
      doubleprecision totein, toteoob, totebio, totedet, toteew
     a, totedeodx, totemod
c
c----------------------------
c
