      subroutine set_spline(x,f,a,b,c,d,n)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      integer*4 n
      real*4 x(n),f(n),a(n),b(n),c(n),d(n)
      integer*4 m1,m2,m,mr
      real*4 s,r
      m1 = 2
      m2 = n-1
      s = 0
        DO 4411 m=1,m2
        d(m) = x(m+1) - x(m)
        r = (f(m+1) - f(m))/d(m)
        c(m) = r - s
        s = r
4411  CONTINUE
4412  CONTINUE
      s=0
      r=0
      c(1)=0
      c(n)=0
        DO 4421 m=m1,m2
        c(m) = c(m) + r*c(m-1)
        b(m) = 2*(x(m-1) - x(m+1)) - r*s
        s = d(m)
        r = s/b(m)
4421  CONTINUE
4422  CONTINUE
      mr = m2
        DO 4431 m=m1,m2
        c(mr) = (d(mr)*c(mr+1) - c(mr))/b(mr)
        mr = mr - 1
4431  CONTINUE
4432  CONTINUE
        DO 4441 m=1,m2
        s = d(m)
        r = c(m+1) - c(m)
        d(m) = r/s
        c(m) = 3*c(m)
        b(m) = (f(m+1)-f(m))/s - (c(m)+r)*s
        a(m) = f(m)
4441  CONTINUE
4442  CONTINUE
      return
      end
      real*4 function spline(s,x,a,b,c,d,n)
      ! Copyright National Research Council of Canada, 2000.
      ! All rights reserved.
      implicit none
      integer*4 n
      real*4 s,x(n),a(n),b(n),c(n),d(n)
      integer*4 m_lower,m_upper,direction,m,ml,mu,mav
      real*4 q
      IF (( x(1) .GT. x(n) )) THEN
        direction = 1
        m_lower = n
        m_upper = 0
      ELSE
        direction = 0
        m_lower = 0
        m_upper = n
      END IF
      IF (( s .GE. x(m_upper + direction) )) THEN
        m = m_upper + 2*direction - 1
      ELSE IF(( s .LE. x(m_lower+1-direction) )) THEN
        m = m_lower - 2*direction + 1
      ELSE
        ml = m_lower
        mu = m_upper
4451    IF(iabs(mu-ml).LE.1)GO TO 4452
          mav = (ml+mu)/2
          IF (( s .LT. x(mav) )) THEN
            mu = mav
          ELSE
            ml = mav
          END IF
        GO TO 4451
4452    CONTINUE
        m = mu + direction - 1
      END IF
      q = s - x(m)
      spline = a(m) + q*(b(m) + q*(c(m) + q*d(m)))
      return
      end
