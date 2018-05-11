      subroutine runsxinit
c
c--------------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
      integer k
c
c--------------------------------------------------------------c
c
      write (*,*) ' '
      read(ltydat,50,err=20,end=20) isxnum
c
      if (isxnum.gt.20) then
        write(*,*)' ISXNUM is TOO large', isxnum
        isxnum = 20
        write(*,*)' ISXNUM HAS LIMITED TO', isxnum
c
      else
        write(*,*)' ISXNUM = ', isxnum
      endif
c
c-------
c
      do 10 k=1,isxnum
c
      read(ltydat,60,err=20,end=25) isxtype(k),isxev(k), sx_cent(k),
     + sy_cent(k),sz_cent(k),rsx(k),rsy(k),rsz(k),esource(k),
     + thetaangsx(k), phiangsx(k)
c
      rsxsq(k) = rsx(k)*rsx(k)
      rsysq(k) = rsy(k)*rsy(k)
      rszsq(k) = rsz(k)*rsz(k)
      if(rsx(k).eq.0.) then
         rsxsqinv(k) = 0.
      else
         rsxsqinv(k) = 1./rsxsq(k)
      endif
      if(rsy(k).eq.0.) then
         rsysqinv(k) = 0.
      else
         rsysqinv(k) = 1./rsysq(k)
      endif
      if(rsz(k).eq.0.) then
         rszsqinv(k) = 0.
      else
         rszsqinv(k) = 1./rszsq(k)
      endif
c
c----    
c
       isxreg(k) = kgeom(sx_cent(k),sy_cent(k),sz_cent(k))
c
c---- 
c
      write(*,30)k,isxtype(k),isxreg(k),isxev(k), sx_cent(k),
     + sy_cent(k),sz_cent(k),rsx(k),rsxsq(k),rsy(k),rsysq(k),
     + rsz(k),rszsq(k),esource(k), thetaangsx(k), phiangsx(k) 
c
c---- find out the angular and geometric source type 
c
c------------- isxangtype:
c
c------------- 1:cone-beam source (isotropic)
c------------- 2:fan-beam source (isotropic on a plane)
c
c
c------------- isxgeomtype:
c
c------------- 1:elliptic source
c------------- 2:cilindrical source //z axis
c ale-190401-- 3:cilindrical source //x axis
c ale-190401-- 4:box source
c
c
      isxangtype(k)  = int(isxtype(k)/10)
      isxgeomtype(k) = mod(isxtype(k),10)
      write(*,*) 'Source ang type:',isxangtype(k),
     +	 '  Source geom type:',isxgeomtype(k) 
      if (isxangtype(k).gt.2) then
         write(*,*)' Source ang type is not correct', isxangtype(k)
      endif 
      if (isxgeomtype(k).gt.4) then
         write(*,*)' Source geom type is not correct', isxgeomtype(k)
      endif 
c
      sx_inf(k) = sx_cent(k) - rsx(k)
      sx_sup(k) = 2. * rsx(k)
c
      sy_inf(k) = sy_cent(k) - rsy(k)
      sy_sup(k) = 2. * rsy(k)
c
      sz_inf(k) = sz_cent(k) - rsz(k)
      sz_sup(k) = 2. * rsz(k)
c
      write(*,40) k ,isxev(k),sx_inf(k),sx_sup(k)
     +, sy_inf(k),sy_sup(k), sz_inf(k),sz_sup(k),esource(k)
ccc      read(ltydat,70,err=20,end=25) thetaangsx(k), phiangsx(k)
c
c------  convert phi and theta in radiants
c
      phiangsx(k) = phiangsx(k) * pi / 180.
      thetaangsx(k) = thetaangsx(k) * pi / 180.
      write(*,*) ' '
      write(*,*)' PHI ANG ',phiangsx(k) ,' THETA ANG ',thetaangsx(k)
c
      cthsup(k) = 1 - cos (thetaangsx(k))
      cthinf(k) = 1.    
      write(*,*)' CTHSUP ', cthsup(k),' CTHINF ', cthinf(k)
      write(*,*) ' '
c
   10 continue
c
c
c--------------------------------------------------------c
c      
      return
c
c--------------------------------------------------------c
c
   20 write(*,*) ' ERROR IN GEOMETRY FILE - RUNSXINIT '
   25 write(*,*) ' END IN GEOMETRY FILE - RUNSXINIT '
      close(ltydat)
c
      return
c
   30 FORMAT( /,1X,'SX PARM,INPUT ' , 3I5,I11,/,1X,10(1x,F9.3))
   40 FORMAT( 1X,'SX PARM,TRASF ' ,I5,2X,I11,/,8(1x,F10.3))
   50 FORMAT(I10)
   60 FORMAT(I10,I10,11F10.0)
   70 FORMAT(2F10.0)
c
       END
c
