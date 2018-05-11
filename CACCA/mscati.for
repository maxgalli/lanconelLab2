      subroutine mscati
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c
      real*4 ededx,ei,eil,eip1,eip1l,si,sip1,eke_orig,elke_orig,aux,
     *ecutmn,tstbm,tstbmn
      real*4 p2,beta2,dedx0,ekef,elkef,estepx,ektmp,elktmp,chi_a2
      integer*4 i,leil,leip1l,neke,lelke,lelkef,lelktmp
      real*4 sigee,sigep,sig
      IF (( bca_algorithm .EQ. 0 )) THEN
        exact_bca = .true.
      ELSE
        exact_bca = .false.
      END IF
      IF (( estepe .LE. 0 .OR. estepe .GE. 1)) THEN
        estepe = 0.25
      END IF
      IF (( ximax .LE. 0 .OR. ximax .GE. 1 )) THEN
        IF (( exact_bca )) THEN
          ximax = 0.5
        ELSE
          ximax = 0.5
        END IF
      END IF
      IF ((transport_algorithm .NE. 0 .AND. transport_algorithm .NE. 1 .
     *AND. transport_algorithm .NE. 2 )) THEN
        transport_algorithm = 0
      END IF
      IF (( skindepth_for_bca .LE. 1 )) THEN
        IF (( transport_algorithm .EQ. 1 .AND. .NOT.exact_bca )) THEN
          write(6,*) ' PRESTA calculates default min. step-size for BCA:        
     * '                                                                        
          ecutmn = 1e30
            DO 3931 i=1,numgeom
            IF (( med(i) .GT. 0 .AND. med(i) .LE. nmed )) THEN
              ecutmn = Min(ecutmn,ecut(i))
            END IF
3931      CONTINUE
3932      CONTINUE
          write(6,*) '     minimum ECUT found: ',ecutmn                         
          tstbmn = 1e30
            DO 3941 medium=1,nmed
            tstbm = (ecutmn-0.5110034)*(ecutmn+0.5110034)/ecutmn**2
            tstbm = blcc(medium)*tstbm*(ecutmn/xcc(medium))**2
            tstbm = Log(tstbm/Log(tstbm))
            tstbmn = Min(tstbmn,tstbm)
3941      CONTINUE
3942      CONTINUE
          write(6,*) '     default BLCMIN is: ',tstbmn                          
          skindepth_for_bca = Exp(tstbmn)
          write(6,*) '     this corresponds to ',skindepth_for_bca, ' el        
     *astic MFPs '                                                              
        ELSE
          skindepth_for_bca = 3
        END IF
      END IF
      call init_ms_SR
        DO 3951 medium=1,nmed
        blcc(medium) = 1.16699413758864573*blcc(medium)
        xcc(medium) = xcc(medium)**2
3951  CONTINUE
3952  CONTINUE
      IF (( spin_effects )) THEN
        call init_spin
      END IF
      write(6,*)
      esige_max = 0
      psige_max = 0
        DO 3961 medium=1,nmed
        sigee = 1E-15
        sigep = 1E-15
        neke = meke(medium)
          DO 3971 i=1,neke
          ei = exp((float(i) - eke0(medium))/eke1(medium))
          eil = log(ei)
          leil = i
          ededx=ededx1(Leil,MEDIUM)*eil+ededx0(Leil,MEDIUM)
          sig=esig1(Leil,MEDIUM)*eil+esig0(Leil,MEDIUM)
          sig = sig/ededx
          IF((sig .GT. sigee))sigee = sig
          ededx=pdedx1(Leil,MEDIUM)*eil+pdedx0(Leil,MEDIUM)
          sig=psig1(Leil,MEDIUM)*eil+psig0(Leil,MEDIUM)
          sig = sig/ededx
          IF((sig .GT. sigep))sigep = sig
3971    CONTINUE
3972    CONTINUE
cale        write(6,*) ' Medium ',medium,' sige = ',sigee,sigep                     
        esig_e(medium) = sigee
        psig_e(medium) = sigep
        IF((sigee .GT. esige_max))esige_max = sigee
        IF((sigep .GT. psige_max))psige_max = sigep
3961  CONTINUE
3962  CONTINUE
      write(6,*)
      write(6,*) ' Initializing tmxs for estepe = ',estepe,' and ximax =        
     * ',ximax                                                                  
      write(6,*)
      rm = 0.5110034
        DO 3981 medium=1,nmed
        ei = exp((1 - eke0(medium))/eke1(medium))
        eil = log(ei)
        leil = 1
        E_array(1,medium) = ei
        expeke1(medium) = Exp(1./eke1(medium))-1
        range_ep(0,1,medium) = 0
        range_ep(1,1,medium) = 0
        neke = meke(medium)
          DO 3991 i=1,neke - 1
          eip1 = exp((float(i + 1) - eke0(medium))/eke1(medium))
          E_array(i+1,medium) = eip1
          eke_orig = 0.5*(eip1+ei)
          elke_orig = Log(eke_orig)
          Lelke=eke1(MEDIUM)*elke_orig+eke0(MEDIUM)
          ededx=pdedx1(Lelke,MEDIUM)*elke_orig+pdedx0(Lelke,MEDIUM)
          aux = pdedx1(i,medium)/ededx
          range_ep(1,i+1,medium) = range_ep(1,i,medium) + (eip1-ei)/eded
     *    x*(1+aux*(1+2*aux)*((eip1-ei)/eke_orig)**2/24)
          ededx=ededx1(Lelke,MEDIUM)*elke_orig+ededx0(Lelke,MEDIUM)
          aux = ededx1(i,medium)/ededx
          range_ep(0,i+1,medium) = range_ep(0,i,medium) + (eip1-ei)/eded
     *    x*(1+aux*(1+2*aux)*((eip1-ei)/eke_orig)**2/24)
          ei = eip1
3991    CONTINUE
3992    CONTINUE
        eil = (1 - eke0(medium))/eke1(medium)
        ei = Exp(eil)
        leil = 1
        p2 = ei*(ei+2*rm)
        beta2 = p2/(p2+rm*rm)
        chi_a2 = Xcc(medium)/(4*p2*blcc(medium))
        dedx0=ededx1(Leil,MEDIUM)*eil+ededx0(Leil,MEDIUM)
        estepx = 2*p2*beta2*dedx0/ei/Xcc(medium)/(Log(1+1./chi_a2)*(1+ch
     *  i_a2)-1)
        estepx = estepx*ximax
        IF (( estepx .GT. estepe )) THEN
          estepx = estepe
        END IF
        si = estepx*ei/dedx0
          DO 4001 i=1,neke - 1
          elke_orig = (i + 1 - eke0(medium))/eke1(medium)
          eke_orig = Exp(elke_orig)
          lelke = i+1
          p2 = eke_orig*(eke_orig+2*rm)
          beta2 = p2/(p2+rm*rm)
          chi_a2 = Xcc(medium)/(4*p2*blcc(medium))
          ededx=ededx1(Lelke,MEDIUM)*elke_orig+ededx0(Lelke,MEDIUM)
          estepx = 2*p2*beta2*ededx/eke_orig/ Xcc(medium)/
     *    (Log(1+1./chi_a2)*(
     *    1+chi_a2)-1)
          estepx = estepx*ximax
          IF (( estepx .GT. estepe )) THEN
            estepx = estepe
          END IF
          ekef = (1-estepx)*eke_orig
          IF (( ekef .LE. E_array(1,medium) )) THEN
            sip1 = (E_array(1,medium) - ekef)/dedx0
            ekef = E_array(1,medium)
            elkef = (1 - eke0(medium))/eke1(medium)
            lelkef = 1
          ELSE
            elkef = Log(ekef)
            Lelkef=eke1(MEDIUM)*elkef+eke0(MEDIUM)
            leip1l = lelkef + 1
            eip1l = (leip1l - eke0(medium))/eke1(medium)
            eip1 = E_array(leip1l,medium)
            aux = (eip1 - ekef)/eip1
            elktmp = 0.5*(elkef+eip1l+0.25*aux*aux*(1+aux*(1+0.875*aux))
     *      )
            ektmp = 0.5*(ekef+eip1)
            lelktmp = lelkef
            ededx=ededx1(Lelktmp,MEDIUM)*elktmp+ededx0(Lelktmp,MEDIUM)
            aux = ededx1(lelktmp,medium)/ededx
            sip1 = (eip1 - ekef)/ededx*( 1+aux*(1+2*aux)*((eip1-ekef)/ek
     *      tmp)**2/24)
          END IF
          sip1 = sip1 + range_ep(0,lelke,medium) - range_ep(0,lelkef+1,m
     *    edium)
          tmxs1(i,medium) = (sip1 - si)*eke1(medium)
          tmxs0(i,medium) = sip1 - tmxs1(i,medium)*elke_orig
          si = sip1
4001    CONTINUE
4002    CONTINUE
        tmxs0(neke,medium) = tmxs0(neke - 1,medium)
        tmxs1(neke,medium) = tmxs1(neke - 1,medium)
3981  CONTINUE
3982  CONTINUE
      return
      end
