      subroutine init_nist_brems
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none


      real*4 energy_array(57),x_array(30), cs_array(57,30,100)
      real*4 xi_array(30)
      real*8 x_gauss(64),w_gauss(64)
      integer*4 nmix,kmix,i,n,k,j,ii
      integer*4 ngauss,i_gauss
      integer*4 ifirst,ilast,nener,neke,leil
      real*4 cs(57,30),ee(57),ele(57)
      real*4 csx(30),afx(30),bfx(30),cfx(30),dfx(30)
      real*4 cse(57),afe(57),bfe(57),cfe(57),dfe(57)
      real*4 Z_orig,sumA
      real*4 emin,xi,res,spline,eil,ei,beta2,aux,sigb,sigt,ebr1,ebr2
      real*4 sigee,sigep,sige,si_esig,si1_esig,si_ebr1,si1_ebr1,ededx, s
     *ig_bhabha,si_psig,si1_psig,si_pbr1,si1_pbr1,si_pbr2,si1_pbr2
      integer*4 iz
      real*4 ple,qle,x_orig,f,error,max_error,x_max_error,f_max_error
      integer*4 ndat,k_max_error
      real*4 amu
      parameter (amu = 1660.5655)
cale
      include 'egs4comm.for'
      include 'egs4fcomm.for'

      character*80 bremfilename

      character*80 totbremfilename

      bremfilename='nist_brems.dat'

      totbremfilename=pathname(1:pathlength) // bremfilename
cale     write(*,*)'cale:',totbremfilename
      open(unit=76,file=totbremfilename,form='formatted',status='old')
cale
      read(76,*)
      read(76,*) nmix,kmix
      IF (( kmix .NE. 30 .OR. nmix .NE. 57)) THEN
        write(6,*) ' init_nist_brems: wrong data file!'                         
        stop
      END IF
      read(76,*) energy_array
      read(76,*) x_array
      read(76,*)
        DO 5081 i=1,100
        read(76,*) ((cs_array(n,k,i),n=1,nmix),k=1,kmix)
5081  CONTINUE
5082  CONTINUE
      close(76)
        DO 5091 k=1,kmix
        xi_array(k)=Log(1-x_array(k)+1e-6)
5091  CONTINUE
5092  CONTINUE
      ngauss = 64
      call gauss_legendre(0d0,1d0,x_gauss,w_gauss,ngauss)
      write(6,*)
      write(6,*) ' Using NIST brems cross sections! ' 
      write(6,*) ' Initializing brems data for media '                          
      write(6,*)
        DO 5101 medium=1,nmed
        log_ap(medium) = log(ap(medium))
cale        write(6,*) ' Initializing brems data for medium ',medium,'...'          
        emin = max(ae(medium) - rm, ap(medium))
          DO 5111 i=1,nmix
          IF((energy_array(i) .GE. emin))GO TO5112
5111    CONTINUE
5112    CONTINUE
        ifirst = i
          DO 5121 i=nmix,1,-1
          IF((energy_array(i) .LT. ue(medium) - rm))GO TO5122
5121    CONTINUE
5122    CONTINUE
        ilast = i+1
        IF (( ifirst .LT. 1 .OR. ilast .GT. nmix )) THEN
          write(6,*) ' init_nist_brems: data available only for '               
          write(6,*) energy_array(i),' <= E <= ',energy_array(nmix)             
          write(6,*) ' will use spline interpolations to get cross '            
          write(6,*) ' sections beyond the available data but this may'         
          write(6,*) ' produce nonsense!'                                       
          IF((ifirst .LT. 1))ifirst=1
          IF((ilast .GT. nmix))ilast = nmix
        END IF
          DO 5131 i=ifirst,ilast
          ii = i+1 - ifirst
          ee(ii) = energy_array(i)
          ele(ii) = log(ee(ii))
          sumA = 0
            DO 5141 j=1,NNE(medium)
            sumA = sumA + pz(medium,j)*wa(medium,j)
5141      CONTINUE
5142      CONTINUE
          sumA = sumA*amu
            DO 5151 k=1,kmix
            cs(ii,k) = 0
              DO 5161 j=1,NNE(medium)
              Z_orig = zelem(medium,j)
              iz = int(Z_orig+0.1)
              Z_orig = Z_orig*Z_orig/sumA
              cs(ii,k) = cs(ii,k) + pz(medium,j)*Z_orig*cs_array(i,k,iz)
5161        CONTINUE
5162        CONTINUE
            csx(k) = Log(cs(ii,k))
5151      CONTINUE
5152      CONTINUE
          call set_spline(xi_array,csx,afx,bfx,cfx,dfx,kmix)
          cse(ii) = 0
          aux = Log(ee(ii)/ap(medium))
            DO 5171 i_gauss=1,ngauss
            xi = log(1 - ap(medium)/ee(ii)*exp(x_gauss(i_gauss)*aux)+1e-
     *      6)
            res = spline(xi,xi_array,afx,bfx,cfx,dfx,kmix)
            cse(ii) = cse(ii) + w_gauss(i_gauss)*exp(res)
5171      CONTINUE
5172      CONTINUE
5131    CONTINUE
5132    CONTINUE
        nener = ilast - ifirst + 1
        call set_spline(ele,cse,afe,bfe,cfe,dfe,nener)
        neke = meke(medium)
        sigee = 1E-15
        sigep = 1E-15
          DO 5181 i=1,neke
          eil = (float(i) - eke0(medium))/eke1(medium)
          ei = exp(eil)
          leil = i
          beta2 = ei*(ei+2*rm)/(ei+rm)**2
          IF (( ei .LE. ap(medium) )) THEN
            sigb = 1e-30
          ELSE
            sigb = spline(eil,ele,afe,bfe,cfe,dfe,nener)
            sigb = sigb*log(ei/ap(medium))/beta2*rho(medium)
          END IF
          sigt=esig1(Leil,MEDIUM)*eil+esig0(Leil,MEDIUM)
          ebr1=ebr11(Leil,MEDIUM)*eil+ebr10(Leil,MEDIUM)
          IF((sigt .LT. 0))sigt = 0
          IF((ebr1 .GT. 1))ebr1 = 1
          IF((ebr1 .LT. 0))ebr1 = 0
          IF (( i .GT. 1 )) THEN
            si_esig = si1_esig
            si_ebr1 = si1_ebr1
            si1_esig = sigt*(1 - ebr1) + sigb
            si1_ebr1 = sigb/si1_esig
            esig1(i-1,medium) = (si1_esig - si_esig)*eke1(medium)
            esig0(i-1,medium) = si1_esig - esig1(i-1,medium)*eil
            ebr11(i-1,medium) = (si1_ebr1 - si_ebr1)*eke1(medium)
            ebr10(i-1,medium) = si1_ebr1 - ebr11(i-1,medium)*eil
          ELSE
            si1_esig = sigt*(1 - ebr1) + sigb
            si1_ebr1 = sigb/si1_esig
          END IF
          sigt=psig1(Leil,MEDIUM)*eil+psig0(Leil,MEDIUM)
          ebr1=pbr11(Leil,MEDIUM)*eil+pbr10(Leil,MEDIUM)
          ebr2=pbr21(Leil,MEDIUM)*eil+pbr20(Leil,MEDIUM)
          IF((sigt .LT. 0))sigt = 0
          IF((ebr1 .GT. 1))ebr1 = 1
          IF((ebr1 .LT. 0))ebr1 = 0
          IF((ebr2 .GT. 1))ebr2 = 1
          IF((ebr2 .LT. 0))ebr2 = 0
          sig_bhabha = sigt*(ebr2 - ebr1)
          IF((sig_bhabha .LT. 0))sig_bhabha = 0
          IF (( i .GT. 1 )) THEN
            si_psig = si1_psig
            si_pbr1 = si1_pbr1
            si_pbr2 = si1_pbr2
            si1_psig = sigt*(1 - ebr1) + sigb
            si1_pbr1 = sigb/si1_psig
            si1_pbr2 = (sigb + sig_bhabha)/si1_psig
            psig1(i-1,medium) = (si1_psig - si_psig)*eke1(medium)
            psig0(i-1,medium) = si1_psig - psig1(i-1,medium)*eil
            pbr11(i-1,medium) = (si1_pbr1 - si_pbr1)*eke1(medium)
            pbr10(i-1,medium) = si1_pbr1 - pbr11(i-1,medium)*eil
            pbr21(i-1,medium) = (si1_pbr2 - si_pbr2)*eke1(medium)
            pbr20(i-1,medium) = si1_pbr2 - pbr20(i-1,medium)*eil
          ELSE
            si1_psig = sigt*(1 - ebr1) + sigb
            si1_pbr1 = sigb/si1_psig
            si1_pbr2 = (sigb + sig_bhabha)/si1_psig
          END IF
          ededx=ededx1(Leil,MEDIUM)*eil+ededx0(Leil,MEDIUM)
          sige = si1_esig/ededx
          IF((sige .GT. sigee))sigee = sige
          ededx=pdedx1(Leil,MEDIUM)*eil+pdedx0(Leil,MEDIUM)
          sige = si1_psig/ededx
          IF((sige .GT. sigep))sigep = sige
5181    CONTINUE
5182    CONTINUE
        esig1(neke,medium) = esig1(neke-1,medium)
        esig0(neke,medium) = esig0(neke-1,medium)
        ebr11(neke,medium) = ebr11(neke-1,medium)
        ebr10(neke,medium) = ebr10(neke-1,medium)
        psig1(neke,medium) = psig1(neke-1,medium)
        psig0(neke,medium) = psig0(neke-1,medium)
        pbr11(neke,medium) = pbr11(neke-1,medium)
        pbr10(neke,medium) = pbr10(neke-1,medium)
        pbr21(neke,medium) = pbr21(neke-1,medium)
        pbr20(neke,medium) = pbr20(neke-1,medium)
cale        write(6,*) ' Max. new cross sections per energy loss: ',sigee,si        
cale     *  gep
        esig_e(medium) = sigee
        psig_e(medium) = sigep
        IF((sigee .GT. esige_max))esige_max = sigee
        IF((sigep .GT. psige_max))psige_max = sigep
        nb_emin(medium) = energy_array(ifirst)
        IF (( nb_emin(medium) .LE. ap(medium) )) THEN
          nb_emin(medium) = energy_array(ifirst+1)
        END IF
        nb_emax(medium) = energy_array(ilast)
        nb_lemin(medium) = log(nb_emin(medium))
        nb_lemax(medium) = log(nb_emax(medium))
        nb_dle(medium) = (nb_lemax(medium) - nb_lemin(medium))/(100-1)
        nb_dlei(medium) = 1/nb_dle(medium)
        eil = nb_lemin(medium) - nb_dle(medium)
          DO 5191 i=1,100
          eil = eil + nb_dle(medium)
          ei = exp(eil)
            DO 5201 ii=1,nener
            IF((ei .LT. ee(ii)))GO TO5202
5201      CONTINUE
5202      CONTINUE
          ii = ii-1
          IF((ii .GT. nener-1))ii = nener-1
          ple = (eil - ele(ii))/(ele(ii+1)-ele(ii))
          qle = 1 - ple
            DO 5211 k=1,30
            csx(k) = log(qle*cs(ii,k) + ple*cs(ii+1,k))
5211      CONTINUE
5212      CONTINUE
          call set_spline(xi_array,csx,afx,bfx,cfx,dfx,kmix)
          x_orig = ap(medium)/ei
          aux = -log(x_orig)
          xi = log(1 - x_orig+1e-6)
          res = spline(xi,xi_array,afx,bfx,cfx,dfx,kmix)
          nb_xdata(0,i,medium) = 0
          nb_fdata(0,i,medium) = exp(res)
            DO 5221 k=1,30
            IF((x_array(k) .GT. x_orig))GO TO5222
5221      CONTINUE
5222      CONTINUE
          IF((k .GT. 30))k = 30
          ndat = 0
            DO 5231 j=k+1,30-1
            ndat = ndat+1
            nb_xdata(ndat,i,medium) = log(x_array(j)/x_orig)/aux
            nb_fdata(ndat,i,medium) = exp(csx(j))
5231      CONTINUE
5232      CONTINUE
          ndat = ndat+1
          nb_xdata(ndat,i,medium) = 1
          nb_fdata(ndat,i,medium) = exp(csx(30))
          IF((ndat .EQ. 50))goto 5240
5251      CONTINUE
            x_max_error = 0
            f_max_error = 0
            k_max_error = 0
            max_error = 0
              DO 5261 k=0,ndat-1
              x_orig = 0.5*(nb_xdata(k,i,medium) +
     a         nb_xdata(k+1,i,medium))
              f = 0.5*(nb_fdata(k,i,medium) + nb_fdata(k+1,i,medium))
              xi = log(1 - ap(medium)/ei*exp(x_orig*aux)+1e-6)
              res = spline(xi,xi_array,afx,bfx,cfx,dfx,kmix)
              res = exp(res)
              error = abs(1-f/res)
              IF (( error .GT. max_error )) THEN
                x_max_error = x_orig
                f_max_error = res
                max_error = error
                k_max_error = k
              END IF
5261        CONTINUE
5262        CONTINUE
            ndat = ndat+1
              DO 5271 k=ndat,k_max_error+2,-1
              nb_xdata(k,i,medium) = nb_xdata(k-1,i,medium)
              nb_fdata(k,i,medium) = nb_fdata(k-1,i,medium)
5271        CONTINUE
5272        CONTINUE
            nb_xdata(k_max_error+1,i,medium) = x_max_error
            nb_fdata(k_max_error+1,i,medium) = f_max_error
            IF(((ndat .EQ. 50)))GO TO5252
          GO TO 5251
5252      CONTINUE
5240      call prepare_alias_table(50,nb_xdata(0,i,medium), nb_fdata(0,i
     *    ,medium),nb_wdata(1,i,medium),nb_idata(1,i,medium))
5191    CONTINUE
5192    CONTINUE
5101  CONTINUE
5102  CONTINUE
      write(6,*)
      write(6,*)
      return
      end
