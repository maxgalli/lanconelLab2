      subroutine fix_brems
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
c
      include 'egs4fcomm.for'
c

      integer*4 medium_orig,i
      real*4 Zt,Zb,Zf,Zg,Zv,fmax1,fmax2,Zi,pi_orig,fc,xi,aux,XSIF,FCOULC
        DO 3691 medium_orig=1,nmed
        log_ap(medium_orig) = log(ap(medium_orig))
        Zt = 0
        Zb = 0
        Zf = 0
          DO 3701 i=1,NNE(medium_orig)
          Zi = ZELEM(medium_orig,i)
          pi_orig = PZ(medium_orig,i)
          fc = FCOULC(Zi)
          xi = XSIF(Zi)
          aux = pi_orig*Zi*(Zi + xi)
          Zt = Zt + aux
          Zb = Zb - aux*Log(Zi)/3
          Zf = Zf + aux*fc
3701    CONTINUE
3702    CONTINUE
        Zv = (Zb - Zf)/Zt
        Zg = Zb/Zt
        fmax1 = 2*(20.863 + 4*Zg) - 2*(20.029 + 4*Zg)/3
        fmax2 = 2*(20.863 + 4*Zv) - 2*(20.029 + 4*Zv)/3
        dl1(1,medium_orig) = (20.863 + 4*Zg)/fmax1
        dl2(1,medium_orig) = -3.242/fmax1
        dl3(1,medium_orig) = 0.625/fmax1
        dl4(1,medium_orig) = (21.12+4*Zg)/fmax1
        dl5(1,medium_orig) = -4.184/fmax1
        dl6(1,medium_orig) = 0.952
        dl1(2,medium_orig) = (20.029+4*Zg)/fmax1
        dl2(2,medium_orig) = -1.93/fmax1
        dl3(2,medium_orig) = -0.086/fmax1
        dl4(2,medium_orig) = (21.12+4*Zg)/fmax1
        dl5(2,medium_orig) = -4.184/fmax1
        dl6(2,medium_orig) = 0.952
        dl1(3,medium_orig) = (20.863 + 4*Zv)/fmax2
        dl2(3,medium_orig) = -3.242/fmax2
        dl3(3,medium_orig) = 0.625/fmax2
        dl4(3,medium_orig) = (21.12+4*Zv)/fmax2
        dl5(3,medium_orig) = -4.184/fmax2
        dl6(3,medium_orig) = 0.952
        dl1(4,medium_orig) = (20.029+4*Zv)/fmax2
        dl2(4,medium_orig) = -1.93/fmax2
        dl3(4,medium_orig) = -0.086/fmax2
        dl4(4,medium_orig) = (21.12+4*Zv)/fmax2
        dl5(4,medium_orig) = -4.184/fmax2
        dl6(4,medium_orig) = 0.952
        dl1(5,medium_orig) = (3*(20.863 + 4*Zg) - (20.029 + 4*Zg))
        dl2(5,medium_orig) = (3*(-3.242) - (-1.930))
        dl3(5,medium_orig) = (3*(0.625)-(-0.086))
        dl4(5,medium_orig) = (2*21.12+8*Zg)
        dl5(5,medium_orig) = (2*(-4.184))
        dl6(5,medium_orig) = 0.952
        dl1(6,medium_orig) = (3*(20.863 + 4*Zg) + (20.029 + 4*Zg))
        dl2(6,medium_orig) = (3*(-3.242) + (-1.930))
        dl3(6,medium_orig) = (3*0.625+(-0.086))
        dl4(6,medium_orig) = (4*21.12+16*Zg)
        dl5(6,medium_orig) = (4*(-4.184))
        dl6(6,medium_orig) = 0.952
        dl1(7,medium_orig) = (3*(20.863 + 4*Zv) - (20.029 + 4*Zv))
        dl2(7,medium_orig) = (3*(-3.242) - (-1.930))
        dl3(7,medium_orig) = (3*(0.625)-(-0.086))
        dl4(7,medium_orig) = (2*21.12+8*Zv)
        dl5(7,medium_orig) = (2*(-4.184))
        dl6(7,medium_orig) = 0.952
        dl1(8,medium_orig) = (3*(20.863 + 4*Zv) + (20.029 + 4*Zv))
        dl2(8,medium_orig) = (3*(-3.242) + (-1.930))
        dl3(8,medium_orig) = (3*0.625+(-0.086))
        dl4(8,medium_orig) = (4*21.12+16*Zv)
        dl5(8,medium_orig) = (4*(-4.184))
        dl6(8,medium_orig) = 0.952
        bpar(2,medium_orig) = dl1(7,medium_orig)/(3*dl1(8,medium_orig)
     *  + dl1(7,medium_orig))
        bpar(1,medium_orig) = 12*dl1(8,medium_orig)/(3*dl1(8, 
     *  medium_orig)+ dl1(7,medium_orig))
3691  CONTINUE
3692  CONTINUE
      return
      end
      real*4 function FCOULC(Z)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      real*4 Z
      real*4 fine,asq
      data fine/137.03604/
      asq = Z/fine
      asq = asq*asq
      FCOULC = asq*(1.0/(1.0+ASQ)+0.20206+ASQ*(-0.0369+ASQ*(0.0083+ASQ*(
     *-0.002))))
      return
      end
      real*4 function XSIF(Z)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      real*4 Z
      integer*4 iZ
      real*4 alrad(4),alradp(4),a1440,a183,FCOULC
      data alrad/5.31,4.79,4.74,4.71/
      data alradp/6.144,5.621,5.805,5.924/
      data a1440/1194.0/,A183/184.15/
      IF (( Z .LE. 4 )) THEN
        iZ = Z
        xsif = alradp(iZ)/(alrad(iZ) - FCOULC(Z))
      ELSE
        xsif = Log(A1440*Z**(-0.666667))/(Log(A183*Z**(-0.33333))-FCOULC
     *  (Z))
      END IF
      return
      end
