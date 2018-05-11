      subroutine source
c
c-----------------------------------------------------------c
c
      implicit none
      include 'egs4comm.for'
      include 'egs4fcomm.for'
c
c-----------------------------------------------------------c
c      sets the kinetic and total energy for monoenergetic
c      gamma  decays
c
      call ranlux(rng_array)
      rng_seed = 1
c
cnico      evkin = enkin
      evkin = esource(nactsx) + abs(float(ichar)) * rm
c
cnico      ein    = enkin
      ein    = evkin
c
      eireal = ein
c
      call hfill (100,eireal,0.,1.)
c
      totein  =  totein + ein
c
c------------- if different geometrical sources
c
c------------- isxgeomtype:
c------------- 1:elliptic source
c------------- 2:cilindrical source
c------------- 3:box source
c
c
      if(isxgeomtype(nactsx).eq.1) then
c
c------------- elliptic source
c	     
   11     continue
          iloop = iloop + 1
c
          xi0 = sx_sup(nactsx)*rng_array(rng_seed)+sx_inf(nactsx)
          rng_seed = rng_seed + 1
          yi0 = sy_sup(nactsx)*rng_array(rng_seed)+sy_inf(nactsx)
          rng_seed = rng_seed + 1
          zi0 = sz_sup(nactsx)*rng_array(rng_seed)+sz_inf(nactsx)
          rng_seed = rng_seed + 1
c
c------------- check if the point(xi0,yi0,zi0) is inside the sphere
c
          if ( (((xi0-sx_cent(nactsx))*(xi0-sx_cent(nactsx)))
     +	                            * rsxsqinv(nactsx)
     +         + ((yi0-sy_cent(nactsx))*(yi0-sy_cent(nactsx)))
     +	                            * rsysqinv(nactsx)
     +         + ((zi0-sz_cent(nactsx))*(zi0-sz_cent(nactsx)))
     +	                           * rszsqinv(nactsx)) - 1 )    21,21,99
c
 99       IF (( rng_seed .GT. 20 )) THEN
               call ranlux(rng_array)
               rng_seed = 1
          END IF
          goto 11
c
c------------------------------------------------
c
   21     nevent = nevent + 1
c
      elseif(isxgeomtype(nactsx).eq.2) then
c
c------------- cilindrical source
c
   12     continue
          iloop = iloop + 1
c
          xi0 = sx_sup(nactsx)*rng_array(rng_seed)+sx_inf(nactsx)
          rng_seed = rng_seed + 1
          yi0 = sy_sup(nactsx)*rng_array(rng_seed)+sy_inf(nactsx)
          rng_seed = rng_seed + 1
          zi0 = sz_sup(nactsx)*rng_array(rng_seed)+sz_inf(nactsx)
          rng_seed = rng_seed + 1
c
c------------- check if the point(xi0,yi0,zi0) is inside the disk
c
          if ( (((xi0-sx_cent(nactsx))*(xi0-sx_cent(nactsx)))
     +	                            * rsxsqinv(nactsx)
     +         + ((yi0-sy_cent(nactsx))*(yi0-sy_cent(nactsx))) 
     +	                            * rsysqinv(nactsx)) - 1 )    22,22,98
c
 98       IF (( rng_seed .GT. 20 )) THEN
               call ranlux(rng_array)
               rng_seed = 1
          END IF
          goto 12
c
c------------------------------------------------
c
   22     nevent = nevent + 1
c
c ale aggiunge per avere anche sorgenti //asse x. 19 04 2001
c 
      elseif(isxgeomtype(nactsx).eq.3) then
c
c------------- cilindrical source // x axis
c
  212     continue
          iloop = iloop + 1
c
          xi0 = sx_sup(nactsx)*rng_array(rng_seed)+sx_inf(nactsx)
          rng_seed = rng_seed + 1
          yi0 = sy_sup(nactsx)*rng_array(rng_seed)+sy_inf(nactsx)
          rng_seed = rng_seed + 1
          zi0 = sz_sup(nactsx)*rng_array(rng_seed)+sz_inf(nactsx)
          rng_seed = rng_seed + 1
c
c------------- check if the projected point(0,yi0,zi0) is inside the disk
c
          if ( (((zi0-sz_cent(nactsx))*(zi0-sz_cent(nactsx)))
     +	                         * rszsqinv(nactsx)
     +         + ((yi0-sy_cent(nactsx))*(yi0-sy_cent(nactsx)))
     +	                         * rsysqinv(nactsx)) - 1 )    222,222,97
c
 97       IF (( rng_seed .GT. 20 )) THEN
               call ranlux(rng_array)
               rng_seed = 1
          END IF
          goto 212	 
c
c------------------------------------------------
c
  222     nevent = nevent + 1
c
c ale cambia box source da 3 a 4. 19 04 2001
c
      elseif(isxgeomtype(nactsx).eq.4) then	
c
c------------- box source
c
          iloop = iloop + 1
c
          xi0 = sx_sup(nactsx)*rng_array(rng_seed)+sx_inf(nactsx)
          rng_seed = rng_seed + 1
          yi0 = sy_sup(nactsx)*rng_array(rng_seed)+sy_inf(nactsx)
          rng_seed = rng_seed + 1
          zi0 = sz_sup(nactsx)*rng_array(rng_seed)+sz_inf(nactsx)
          rng_seed = rng_seed + 1
c
c------------------------------------------------
c
          nevent = nevent + 1
c
      endif
c
c----------
c
      xin      =   xi0
      yin      =   yi0
      zin      =   zi0
c
      call hfill (101,xin,0.,1.)
      call hfill (102,yin,0.,1.)
      call hfill (103,zin,0.,1.)
      call hfill (200,xin,yin,1.)
      call hfill (201,xin,zin,1.)
      call hfill (202,yin,zin,1.)
cc
c      call hfill (300,xi,yi,zi)
cc
c------------- parallel-beam source
c
      if(thetaangsx(nactsx).eq.0.) then
c
          cth = 1.
          sth = 0.
          phibis = 0.
          goto 31
c
      endif
c
c------------- if different angular sources
c
c------------- isxangtype:
c------------- 1:cone-beam source (isotropic)
c------------- 2:fan-beam source (isotropic on a plane)
c
c
      IF (( rng_seed .GT. 18 )) THEN
              call ranlux(rng_array)
              rng_seed = 5
      END IF
c
      if(isxangtype(nactsx).eq.1) then
c
c------------- cone-beam source (isotropic)
c
          cth = cthinf(nactsx) - cthsup(nactsx) * rng_array(rng_seed)
          rng_seed = rng_seed + 1
          sth = sqrt ( 1.- cth*cth )
c
          phibis = - pi +  twopi * rng_array(rng_seed)
          rng_seed = rng_seed + 1
c
      elseif(isxangtype(nactsx).eq.2) then
c
c------------- fan-beam source (isotropic on a plane)
c
          thetarandom=rng_array(rng_seed)*thetaangsx(nactsx)
          rng_seed = rng_seed + 1
          cth = cos(thetarandom)
          sth = sin(thetarandom)
c
          if (rng_array(rng_seed).ge.0.5) then
                phibis = phiangsx(nactsx)
          else
                phibis = phiangsx(nactsx) + pi
          endif
          rng_seed = rng_seed + 1
c
      endif
c
c-----------------------------------------------------------
c
 31   win = cth
      uin = sth * cos (phibis)
      vin = sth * sin (phibis)             
c
      dnorm = sqrt(uin*uin+vin*vin+win*win)

      win = win / dnorm
      uin = uin / dnorm
      vin = vin / dnorm
c
      call hfill (110,uin,0.,1.)
      call hfill (111,vin,0.,1.)
      call hfill (112,win,0.,1.)
c
c-----------------------------------------------------------
c
      irin =  kgeom(xin,yin,zin)
      wtin =  1.00000
c
      return
      end
c
c-----------------------------------------------------
c
