      subroutine gauss_legendre(x1,x2,x,w,n)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      integer*4 n
      real*8 x1,x2,x(n),w(n)
      real*8 eps,Pi
      parameter (eps = 3.D-14,Pi=3.141592654D0)
      integer*4 i,m,j
      real*8 xm,xl,z,z1,p1,p2,p3,pp
      m = (n + 1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
        DO 5331 i=1,m
        z=cos(Pi*(i-.25d0)/(n+.5d0))
5341    CONTINUE
          p1=1.d0
          p2=0.d0
            DO 5351 j=1,n
            p3 = p2
            p2 = p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
5351      CONTINUE
5352      CONTINUE
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
          IF(((abs(z-z1) .LT. eps)))GO TO5342
        GO TO 5341
5342    CONTINUE
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
5331  CONTINUE
5332  CONTINUE
      return
      end
