      subroutine init_spin
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c

      real*4 eta_array(0:1,0: 31), c_array(0:1,0: 31),g_array(0:1,0: 31)
     *, earray(0: 31),tmp_array(0: 31), sum_Z2,sum_Z,sum_A,sum_pz,Z_orig
     *,tmp,Z23,g_m,g_r,sig,dedx, dum1,dum2,dum3,aux_o,tau,tauc,beta2,eta
     *,gamma,fmax, eil,e_orig,si1e,si2e,si1p,si2p,aae,etap, elarray
     *(0: 31),farray(0: 31), af(0: 31),bf(0: 31),cf(0: 31), df(0: 31),
     * spline
      integer*4 medium_orig,iq_orig,i,j,k,i_ele,iii,iZ,iiZ,n_ener,n_q,
     *n_point,je,neke, ndata,leil,length

      character*80 spin_file
cale
      include 'egs4comm.for'

      character*80 prespin_file
cale
      character*6 string
      integer*4 lnblnk1
      real*4 fine,TF_constant
      parameter (fine=137.03604,TF_constant=0.88534138)
cale      call getenv('HEN_HOUSE',spin_file)
      prespin_file='spinms/z000'
      spin_file=pathname(1:pathlength)//prespin_file
cale      write(*,*)'cale:',spin_file
cale
      length = lnblnk1(spin_file)
c
c mirko - modified spin_files location
c
cale      IF (( spin_file(length:length) .NE. 's' )) THEN                           
cale        length = length + 1
cale        spin_file(length:length) = 's'                                          
cale      END IF
cale      spin_file(length+1:length+11) = 'pinms/z000'                             
cale      length = lnblnk1(spin_file)
c
c mirko - end of modified spin_files location
c
cale


        write(6,*) '  Initializing spin data for media '

        DO 4091 medium_orig=1,NMED
cale        write(6,'(a,i4,a,$)') '  Initializing spin data for
cale     * medium_orig ', medium_orig, ' ..................... '
          DO 4101 iq_orig=0,1
            DO 4111 i=0, 31
            eta_array(iq_orig,i)=0
            c_array(iq_orig,i)=0
            g_array(iq_orig,i)=0
              DO 4121 j=0,15
                DO 4131 k=0,31
                spin_rej(medium_orig,iq_orig,i,j,k) = 0
4131          CONTINUE
4132          CONTINUE
4121        CONTINUE
4122        CONTINUE
4111      CONTINUE
4112      CONTINUE
4101    CONTINUE
4102    CONTINUE
        sum_Z2=0
        sum_A=0
        sum_pz=0
        sum_Z=0
          DO 4141 i_ele=1,NNE(medium_orig)
          Z_orig = ZELEM(medium_orig,i_ele)
          iZ = int(Z_orig+0.5)
          tmp = PZ(medium_orig,i_ele)*Z_orig*(Z_orig+1)
          iii = iZ/100
          spin_file(length-2:length-2) = char(iii+48)
          iiZ = iZ - iii*100
          iii = iiZ/10
          spin_file(length-1:length-1) = char(iii+48)
          iiZ = iiZ - 10*iii
          spin_file(length:length) = char(iiZ+48)
          open(61,file=spin_file,status='old',err=4150)   
          read(61,*) espin_min,espin_max,b2spin_min,b2spin_max
          read(61,*) n_ener,n_q,n_point
          IF (( n_ener .NE. 15 .OR. n_q .NE. 15 .OR. n_point .NE. 31)) 
     *    THEN
            write(6,*) ' Wrong spin file for Z = ',iZ                           
            stop
          END IF
          sum_Z2 = sum_Z2 + tmp
          sum_Z = sum_Z + PZ(medium_orig,i_ele)*Z_orig
          sum_A = sum_A + PZ(medium_orig,i_ele)*WA(medium_orig,i_ele)
          sum_pz = sum_pz + PZ(medium_orig,i_ele)
          Z23 = Z_orig**0.6666667
            DO 4161 iq_orig=0,1
            read(61,*)
            read(61,*)
              DO 4171 i=0, 31
              read(61,'(a,g14.6)') string,earray(i)                             
              read(61,*) dum1,dum2,dum3,aux_o
              eta_array(iq_orig,i)=eta_array(iq_orig,i)+tmp*Log
     * (Z23*aux_o)
              tau = earray(i)/511.0034
              beta2 = tau*(tau+2)/(tau+1)**2
              eta = Z23/(fine*TF_constant)**2*aux_o/4/tau/(tau+2)
              c_array(iq_orig,i)=c_array(iq_orig,i)+ tmp*(Log(1+1/eta)
     *        -1/(1+eta))*dum1*dum3
              g_array(iq_orig,i)=g_array(iq_orig,i)+tmp*dum2
                DO 4181 j=0,15
                read(61,*) tmp_array
                  DO 4191 k=0,31
                  spin_rej(medium_orig,iq_orig,i,j,k) = spin_rej(
     *            medium_orig,iq_orig,i,j,k)+ tmp*tmp_array(k)
4191            CONTINUE
4192            CONTINUE
4181          CONTINUE
4182          CONTINUE
4171        CONTINUE
4172        CONTINUE
4161      CONTINUE
4162      CONTINUE
          close(61)
4141    CONTINUE
4142    CONTINUE
          DO 4201 iq_orig=0,1
            DO 4211 i=0, 31
              DO 4221 j=0,15
              fmax = 0
                DO 4231 k=0,31
                IF (( spin_rej(medium_orig,iq_orig,i,j,k) .GT. fmax ))
     * THEN
                  fmax = spin_rej(medium_orig,iq_orig,i,j,k)
                END IF
4231          CONTINUE
4232          CONTINUE
                DO 4241 k=0,31
                spin_rej(medium_orig,iq_orig,i,j,k) = spin_rej(medium
     *          _orig,iq_orig,i,j,k)/fmax
4241          CONTINUE
4242          CONTINUE
4221        CONTINUE
4222        CONTINUE
4211      CONTINUE
4212      CONTINUE
4201    CONTINUE
4202    CONTINUE
          DO 4251 i=0, 31
          tau = earray(i)/511.0034
          beta2 = tau*(tau+2)/(tau+1)**2
            DO 4261 iq_orig=0,1
            aux_o = Exp(eta_array(iq_orig,i)/sum_Z2)/(fine*TF_constant)
     * **2
            eta_array(iq_orig,i) = 0.26112447*aux_o*blcc(medium_orig)
     * /xcc(medium_orig)
            eta = aux_o/4/tau/(tau+2)
            gamma = 3*(1+eta)*(Log(1+1/eta)*(1+2*eta)-2)/ (Log(1+1/eta)
     *      *(1+eta)-1)
            g_array(iq_orig,i) = g_array(iq_orig,i)/sum_Z2/gamma
            c_array(iq_orig,i) = c_array(iq_orig,i)/sum_Z2/(Log(1+1/eta)
     *      -1/(1+eta))
4261      CONTINUE
4262      CONTINUE
4251    CONTINUE
4252    CONTINUE
        espin_min = espin_min/1000
        espin_max = espin_max/1000
        dlener = Log(espin_max/espin_min)/15
        dleneri = 1/dlener
        espml = Log(espin_min)
        dbeta2 = (b2spin_max-b2spin_min)/15
        dbeta2i = 1/dbeta2
        eil = (1 - eke0(medium_orig))/eke1(medium_orig)
        e_orig = Exp(eil)
        IF (( e_orig .LE. espin_min )) THEN
          si1e = eta_array(0,0)
          si1p = eta_array(1,0)
        ELSE
          IF (( e_orig .LE. espin_max )) THEN
            aae = (eil-espml)*dleneri
            je = aae
            aae = aae - je
          ELSE
            tau = e_orig/0.5110034
            beta2 = tau*(tau+2)/(tau+1)**2
            aae = (beta2 - b2spin_min)*dbeta2i
            je = aae
            aae = aae - je
            je = je + 15 + 1
          END IF
          si1e = (1-aae)*eta_array(0,je) + aae*eta_array(0,je+1)
          si1p = (1-aae)*eta_array(1,je) + aae*eta_array(1,je+1)
        END IF
        neke = meke(medium_orig)
          DO 4271 i=1,neke - 1
          eil = (i+1 - eke0(medium_orig))/eke1(medium_orig)
          e_orig = Exp(eil)
          IF (( e_orig .LE. espin_min )) THEN
            si2e = eta_array(0,0)
            si2p = eta_array(1,0)
          ELSE
            IF (( e_orig .LE. espin_max )) THEN
              aae = (eil-espml)*dleneri
              je = aae
              aae = aae - je
            ELSE
              tau = e_orig/0.5110034
              beta2 = tau*(tau+2)/(tau+1)**2
              aae = (beta2 - b2spin_min)*dbeta2i
              je = aae
              aae = aae - je
              je = je + 15 + 1
            END IF
            si2e = (1-aae)*eta_array(0,je) + aae*eta_array(0,je+1)
            si2p = (1-aae)*eta_array(1,je) + aae*eta_array(1,je+1)
          END IF
          etae_ms1(i,medium_orig) = (si2e - si1e)*eke1(medium_orig)
          etae_ms0(i,medium_orig) = si2e - etae_ms1(i,medium_orig)*eil
          etap_ms1(i,medium_orig) = (si2p - si1p)*eke1(medium_orig)
          etap_ms0(i,medium_orig) = si2p - etap_ms1(i,medium_orig)*eil
          si1e = si2e
          si1p = si2p
4271    CONTINUE
4272    CONTINUE
        etae_ms1(neke,medium_orig) = etae_ms1(neke-1,medium_orig)
        etae_ms0(neke,medium_orig) = etae_ms0(neke-1,medium_orig)
        etap_ms1(neke,medium_orig) = etap_ms1(neke-1,medium_orig)
        etap_ms0(neke,medium_orig) = etap_ms0(neke-1,medium_orig)
          DO 4281 i=0,15
          elarray(i) = Log(earray(i)/1000)
          farray(i) = c_array(0,i)
4281    CONTINUE
4282    CONTINUE
          DO 4291 i=15+1, 31-1
          elarray(i) = Log(earray(i+1)/1000)
          farray(i) = c_array(0,i+1)
4291    CONTINUE
4292    CONTINUE
        ndata =  31+1
        IF (( ue(medium_orig) .GT. 1e5 )) THEN
          elarray(ndata-1) = Log(ue(medium_orig))
        ELSE
          elarray(ndata-1) = Log(1e5)
        END IF
        farray(ndata-1) = 1
        call set_spline(elarray,farray,af,bf,cf,df,ndata)
        eil = (1 - eke0(medium_orig))/eke1(medium_orig)
        si1e = spline(eil,elarray,af,bf,cf,df,ndata)
          DO 4301 i=1,neke-1
          eil = (i+1 - eke0(medium_orig))/eke1(medium_orig)
          si2e = spline(eil,elarray,af,bf,cf,df,ndata)
          q1ce_ms1(i,medium_orig) = (si2e - si1e)*eke1(medium_orig)
          q1ce_ms0(i,medium_orig) = si2e - q1ce_ms1(i,medium_orig)*eil
          si1e = si2e
4301    CONTINUE
4302    CONTINUE
        q1ce_ms1(neke,medium_orig) = q1ce_ms1(neke-1,medium_orig)
        q1ce_ms0(neke,medium_orig) = q1ce_ms0(neke-1,medium_orig)
          DO 4311 i=0,15
          farray(i) = c_array(1,i)
4311    CONTINUE
4312    CONTINUE
          DO 4321 i=15+1, 31-1
          farray(i) = c_array(1,i+1)
4321    CONTINUE
4322    CONTINUE
        call set_spline(elarray,farray,af,bf,cf,df,ndata)
        eil = (1 - eke0(medium_orig))/eke1(medium_orig)
        si1e = spline(eil,elarray,af,bf,cf,df,ndata)
          DO 4331 i=1,neke-1
          eil = (i+1 - eke0(medium_orig))/eke1(medium_orig)
          si2e = spline(eil,elarray,af,bf,cf,df,ndata)
          q1cp_ms1(i,medium_orig) = (si2e - si1e)*eke1(medium_orig)
          q1cp_ms0(i,medium_orig) = si2e - q1cp_ms1(i,medium_orig)*eil
          si1e = si2e
4331    CONTINUE
4332    CONTINUE
        q1cp_ms1(neke,medium_orig) = q1cp_ms1(neke-1,medium_orig)
        q1cp_ms0(neke,medium_orig) = q1cp_ms0(neke-1,medium_orig)
          DO 4341 i=0,15
          farray(i) = g_array(0,i)
4341    CONTINUE
4342    CONTINUE
          DO 4351 i=15+1, 31-1
          farray(i) = g_array(0,i+1)
4351    CONTINUE
4352    CONTINUE
        call set_spline(elarray,farray,af,bf,cf,df,ndata)
        eil = (1 - eke0(medium_orig))/eke1(medium_orig)
        si1e = spline(eil,elarray,af,bf,cf,df,ndata)
          DO 4361 i=1,neke-1
          eil = (i+1 - eke0(medium_orig))/eke1(medium_orig)
          si2e = spline(eil,elarray,af,bf,cf,df,ndata)
          q2ce_ms1(i,medium_orig) = (si2e - si1e)*eke1(medium_orig)
          q2ce_ms0(i,medium_orig) = si2e - q2ce_ms1(i,medium_orig)*eil
          si1e = si2e
4361    CONTINUE
4362    CONTINUE
        q2ce_ms1(neke,medium_orig) = q2ce_ms1(neke-1,medium_orig)
        q2ce_ms0(neke,medium_orig) = q2ce_ms0(neke-1,medium_orig)
          DO 4371 i=0,15
          farray(i) = g_array(1,i)
4371    CONTINUE
4372    CONTINUE
          DO 4381 i=15+1, 31-1
          farray(i) = g_array(1,i+1)
4381    CONTINUE
4382    CONTINUE
        call set_spline(elarray,farray,af,bf,cf,df,ndata)
        eil = (1 - eke0(medium_orig))/eke1(medium_orig)
        si1e = spline(eil,elarray,af,bf,cf,df,ndata)
          DO 4391 i=1,neke-1
          eil = (i+1 - eke0(medium_orig))/eke1(medium_orig)
          si2e = spline(eil,elarray,af,bf,cf,df,ndata)
          q2cp_ms1(i,medium_orig) = (si2e - si1e)*eke1(medium_orig)
          q2cp_ms0(i,medium_orig) = si2e - q2cp_ms1(i,medium_orig)*eil
          si1e = si2e
4391    CONTINUE
4392    CONTINUE
        q2cp_ms1(neke,medium_orig) = q2cp_ms1(neke-1,medium_orig)
        q2cp_ms0(neke,medium_orig) = q2cp_ms0(neke-1,medium_orig)
        tauc = te(medium_orig)/0.5110034
        si1e = 1
          DO 4401 i=1,neke-1
          eil = (i+1 - eke0(medium_orig))/eke1(medium_orig)
          e_orig = Exp(eil)
          leil=i+1
          tau=e_orig/0.5110034
          IF (( tau .GT. 2*tauc )) THEN
            sig=esig1(Leil,medium_orig)*eil+esig0(Leil,medium_orig)
            dedx=ededx1(Leil,medium_orig)*eil+ededx0(Leil,medium_orig)
            sig = sig/dedx
            IF (( sig .GT. 1e-6 )) THEN
              etap=etae_ms1(Leil,medium_orig)*eil+etae_ms0(Leil,
     *  medium_orig)
              eta = 0.25*etap*xcc(medium_orig)/blcc(medium_orig)/
     *  tau/(tau+2)
              g_r = (1+2*eta)*Log(1+1/eta)-2
              g_m = Log(0.5*tau/tauc)+ (1+((tau+2)/(tau+1))**2)*Log(2*(t
     *        au-tauc+2)/(tau+4))- 0.25*(tau+2)*(tau+2+2*(2*tau+1)/(tau+
     *        1)**2)* Log((tau+4)*(tau-tauc)/tau/(tau-tauc+2))+ 0.5*(tau
     *        -2*tauc)*(tau+2)*(1/(tau-tauc)-1/(tau+1)**2)
              IF (( g_m .LT. g_r )) THEN
                g_m = g_m/g_r
              ELSE
                g_m = 1
              END IF
              si2e = 1 - g_m*sum_Z/sum_Z2
            ELSE
              si2e = 1
            END IF
          ELSE
            si2e = 1
          END IF
          blcce1(i,medium_orig) = (si2e - si1e)*eke1(medium_orig)
          blcce0(i,medium_orig) = si2e - blcce1(i,medium_orig)*eil
          si1e = si2e
4401    CONTINUE
4402    CONTINUE
        blcce1(neke,medium_orig) = blcce1(neke-1,medium_orig)
        blcce0(neke,medium_orig) = blcce0(neke-1,medium_orig)
cale        write(6,'(a)') ' done'                                                  
4091  CONTINUE
4092  CONTINUE
      return
4150  write(6,*) ' ******************** Error in init_spin *************        
     ******* '                                                                  
      write(6,'(a,a)') '  could not open file ',spin_file                       
      write(6,*) ' terminating execution '                                      
      write(6,*) ' *****************************************************        
     ********'                                                                  
      stop
      end
