c
cale
c      integer MAXR,MAXR1
c      parameter(MAXR=120000,MAXR1=120001)
cale
      common / egs_vr1 / e_max_rr(120001)
      common / egs_vr2 / i_do_rr(120001)
      integer*4 i_do_rr
      real*4 e_max_rr
c
      common / egs_vr / prob_RR, nbr_split, i_play_RR, i_survived_RR
     *, n_RR_warning
      integer*4 nbr_split,i_play_RR,i_survived_RR,n_RR_warning
      real*4 prob_RR
c
cnico aggiunge questi common
c
      common/ compton_data2 / ibcmp(120001)
      integer*4 ibcmp
c
      common / EDGE2 / IEDGFL(120001),IPHTER(120001)
      integer*4 IEDGFL,  IPHTER
c
      COMMON/BOUNDS/ECUT(120001),PCUT(120001),VACDST
      real*4 ECUT,  PCUT,  VACDST
c
      COMMON/MISC/ DUNIT, KMPI, KMPO
      COMMON/MISC2/ RHOR(120001)
      COMMON/MISC3/ MED(120001)
      COMMON/MISC4/ IRAYLR(120001)
      real*4 DUNIT,  RHOR
      integer*4 MED,  IRAYLR,  KMPI,  KMPO
c
      common/ET_control2/ smaxir(120001)
      common/ET_control3/ estepr(120001)
      real*4 smaxir,  estepr
c
c      integer numgeom
c      data numgeom / 120001 /
      common/numgeom/numgeom
      integer numgeom
c
      COMMON/QDEBUG/QDEBUG
      LOGICAL QDEBUG
c
      COMMON/STACK/ E(50),X(50),Y(50),Z(50),U(50),V(50),W(50),DNEAR(50),
     *WT(50),IQ(50),IR(50),LATCH(50), LATCHI,NP,NPold
      DOUBLE PRECISION E
      real*4 X,Y,Z,  U,V,W,  DNEAR,  WT
      integer*4 IQ,  IR,  LATCH,  LATCHI, NP,  NPold
c
      COMMON/UPHIOT/THETA,SINTHE,COSTHE,SINPHI, COSPHI,PI,TWOPI,PI5D2
      real*4 THETA,  SINTHE,  COSTHE,  SINPHI,  COSPHI,  PI,TWOPI,PI5D2
c
      COMMON/USEFUL/PZERO,PRM,PRMT2,RM,MEDIUM,MEDOLD
      DOUBLE PRECISION PZERO,  PRM,  PRMT2
      real*4 RM
      integer*4 MEDIUM,  MEDOLD
c
      common/randomm/ rng_array(24), seeds(24), rng_seed
      real*4 rng_array
      integer*4 rng_seed,  seeds
c
      COMMON/THRESH/RMT2,RMSQ, AP(21),AE(21),UP(21),UE(21),TE(21)
     *,THMOLL(21)
      real*4 RMT2,  RMSQ,  AP,  AE,  UP,  UE,  TE,  THMOLL
c
      COMMON / BREMPR / DL1(8,21), DL2(8,21), DL3(8,21), DL4(8,21)
     *, DL5(8,21), DL6(8,21), ALPHI(2,21), BPAR(2,21), DELPOS(2,21)
     *, WA(21,50), PZ(21,50), ZELEM(21,50), RHOZ(21,50), PWR2I(50)
     *, DELCM(21), ZBRANG(21), NNE(21), IBRDST, IPRDST, ibr_nist
     *, ASYM(21,50,2)
      CHARACTER*4 ASYM
      real*4 DL1,DL2,DL3,DL4,DL5,DL6,   ALPHI,  BPAR,  DELPOS,  WA,  PZ,
     *  ZELEM,  RHOZ,  PWR2I,  DELCM,  ZBRANG
      integer*4 NNE,  IBRDST,  IPRDST,  ibr_nist
c
      COMMON/EPCONT/EDEP,TSTEP,TUSTEP,USTEP,TVSTEP,VSTEP, RHOF,EOLD,ENEW
     *,EKE,ELKE,GLE,E_RANGE, IDISC,IROLD,IRNEW,IAUSFL(28)
      DOUBLE PRECISION EDEP
      real*4 TSTEP,  TUSTEP,  USTEP,  VSTEP,  TVSTEP,  RHOF,  EOLD,  ENE
     *W,  EKE,  ELKE,  GLE,  E_RANGE
      integer*4 IDISC,  IROLD,  IRNEW,  IAUSFL
c
      common / nist_brems / nb_fdata(0:50,100,21), nb_xdata(0:50,100,21)
     *, nb_wdata(50,100,21), nb_idata(50,100,21), nb_emin(21)
     *, nb_emax(21), nb_lemin(21), nb_lemax(21), nb_dle(21), nb_dlei(21)
     *, log_ap(21)
      real*4 nb_fdata,nb_xdata,nb_wdata,nb_emin,nb_emax,nb_lemin,nb_lema
     *x, nb_dle,nb_dlei,log_ap
      integer*4 nb_idata
c
      common / compton_data / iz_array(1538), be_array(1538)
     *, Jo_array(1538), erfJo_array(1538), ne_array(1538)
     *, shn_array(1538), shell_array(50,21), eno_array(50,21)
     *, n_shell(21)
      integer*4 iz_array,ne_array,shn_array, shell_array,n_shell
      real*4 be_array,Jo_array,erfJo_array,eno_array
c
      COMMON / ELECIN / esig_e(21), psig_e(21), esige_max, psige_max
     *, range_ep(0:1,150,21), E_array(150,21), etae_ms0(150,21)
     *, etae_ms1(150,21), etap_ms0(150,21), etap_ms1(150,21)
     *, q1ce_ms0(150,21), q1ce_ms1(150,21), q1cp_ms0(150,21)
     *, q1cp_ms1(150,21), q2ce_ms0(150,21), q2ce_ms1(150,21)
     *, q2cp_ms0(150,21), q2cp_ms1(150,21), blcce0(150,21)
     *, blcce1(150,21), EKE0(21), EKE1(21), XR0(21), TEFF0(21), BLCC(21)
     *, XCC(21), ESIG0(150,21), ESIG1(150,21), PSIG0(150,21)
     *, PSIG1(150,21), EDEDX0(150,21), EDEDX1(150,21), PDEDX0(150,21)
     *, PDEDX1(150,21), EBR10(150,21), EBR11(150,21), PBR10(150,21)
     *, PBR11(150,21), PBR20(150,21), PBR21(150,21), TMXS0(150,21)
     *, TMXS1(150,21), expeke1(21), IUNRST(21), EPSTFL(21), IAPRIM(21)
      real*4 esig_e, psig_e, esige_max, psige_max, range_ep, E_array, 
     * etae_ms0,etae_ms1,  etap_ms0,etap_ms1,  q1ce_ms0,q1ce_ms1,
     *q1cp_ms0,q1cp_ms1,  q2ce_ms0,q2ce_ms1,  q2cp_ms0,q2cp_ms1, blcce0
     *,blcce1, expeke1,EKE0,EKE1, XR0, TEFF0, BLCC, XCC, ESIG0,ESIG1,  
     *PSIG0,PSIG1,  EDEDX0,EDEDX1,  PDEDX0,PDEDX1,  EBR10,EBR11,
     * PBR10,PBR11,  PBR20,PBR21,  TMXS0,TMXS1
      integer*4 IUNRST,  EPSTFL,  IAPRIM
c
      COMMON / MEDIA / RLC(21), RLDU(21), RHO(21), MSGE(21), MGE(21)
     *, MSEKE(21), MEKE(21), MLEKE(21), MCMFP(21), MRANGE(21)
     *, IRAYLM(21), MEDIA(24,21), NMED
      CHARACTER*4 MEDIA
      real*4 RLC,  RLDU,  RHO
      integer*4 MSGE,  MGE,  MSEKE, MEKE,  MLEKE, MCMFP, MRANGE, IRAYLM,
     * NMED
c
      COMMON / PHOTIN / EBINDA(21), GE0(21), GE1(21), GMFP0(200,21)
     *, GMFP1(200,21), GBR10(200,21), GBR11(200,21), GBR20(200,21)
     *, GBR21(200,21), RCO0(21), RCO1(21), RSCT0(100,21), RSCT1(100,21)
     *, COHE0(200,21), COHE1(200,21), MPGEM(1,21), NGR(21)
      real*4 EBINDA,  GE0,GE1,  GMFP0,GMFP1,  GBR10,GBR11,  GBR20,GBR21,
     *  RCO0,RCO1,  RSCT0,RSCT1,  COHE0,COHE1
      integer*4 MPGEM,  NGR
c
      COMMON/EDGE/binding_energies(6,100), interaction_prob(6,100), rela
     *xation_prob(39,100), edge_energies(16,100), edge_number(100), edge
     *_a(16,100), edge_b(16,100), edge_c(16,100), edge_d(16,100)
      real*4 binding_energies,  interaction_prob,    relaxation_prob,  e
     *dge_energies,  edge_a,edge_b,edge_c,edge_d
      integer*4 edge_number
c
      COMMON/UPHIIN/SINC0,SINC1,SIN0(1002),SIN1(1002)
      real*4 SINC0,SINC1,SIN0,SIN1
c
      common/ET_control/ estepe, ximax, skindepth_for_bca
     *, transport_algorithm, bca_algorithm, exact_bca, spin_effects
      real*4 estepe,  ximax, skindepth_for_bca
      integer*4 transport_algorithm, bca_algorithm
      logical exact_bca,  spin_effects
c
      common/CH_steps/ count_pII_steps,count_all_steps
      real*8 count_pII_steps,count_all_steps
c
      common/spin_data/ spin_rej(5,0:1,0: 31,0:15,0:31), espin_min,espin
     *_max,espml,b2spin_min,b2spin_max, dbeta2,dbeta2i,dlener,dleneri,dq
     *q1,dqq1i
      real*4 spin_rej,espin_min,espin_max,espml,b2spin_min,b2spin_max, d
     *beta2,dbeta2i,dlener,dleneri,dqq1,dqq1i
c
      common/ms_data/ ums_array(0:63,0:7,0:31), fms_array(0:63,0:7,0:31)
     *, wms_array(0:63,0:7,0:31), ims_array(0:63,0:7,0:31), llammin,llam
     *max,dllamb,dllambi,dqms,dqmsi
      real*4 ums_array,fms_array,wms_array, llammin,llammax,dllamb,dllam
     *bi,dqms,dqmsi
      integer*2 ims_array
c
      common/slate/isl(40)
      integer isl
