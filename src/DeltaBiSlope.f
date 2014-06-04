CCCCCCCCCCCCCCCCCCCC
C     TARGET DELTA1 EQUATION
CCCCCCCCCCCCCCCCCCCC

      real*8 function target1slope(delta1,gamma1,beta1,sigma1,
     &     p,tau,x,xdim)
      integer xdim, i
      real*8 delta1, gamma1(xdim), beta1(xdim), sigma1,
     &     p, tau, x(xdim)
      real*8 pnrm

      real*8 quan, lp

      quan = dot_product(gamma1, x)
      lp = dot_product(beta1, x)

      target1slope=tau-p*pnrm(quan-delta1 - lp,0.d0,sigma1,1,0)-
     &     (1-p)*pnrm(quan-delta1+lp,0.d0,sigma1,1,0)

      return
      end

C------------------------------
C     myzero1: bisection method to solve for delta1
C------------------------------
      real*8 function myzero1slope(gamma1,beta1, sigma1, p,tau, x, xdim)
      integer xdim, i
      real*8 delta1, gamma1(xdim), beta1(xdim), sigma1, p,
     &     tau, x(xdim)
      real*8 target1slope

      real*8 a, b, fa, fb, c, fc, tol, prevstep, newstep
      logical success
      real*8 dx, t1, cb, t2, pp, q
      integer maxit

      tol = 0.00001
      maxit = 40

      a = -100
      b = 100

      fa = target1slope(a, gamma1, beta1, sigma1, p, tau, x, xdim)
      fb = target1slope(b, gamma1, beta1, sigma1, p, tau, x, xdim)
      c = a
      fc = fa

C     First test if root is an endpoint
      if (abs(fa) .le. tol) then
         myzero1slope = a
         return
      end if

      if (abs(fb) .le. tol) then
         myzero1slope = b
         return
      end if

      if (fa * fb .ge. 0) print*, 'root is not included solving D1.'

      do while (maxit .gt. 0)
         c = (a + b)/2.d0

         fc = target1slope(c, gamma1, beta1, sigma1, p, tau, x, xdim)

         if (abs(fc) .le. tol) then
            myzero1slope = c
            return
         end if

         if (fa * fc .gt. 0) then
            a = c
            fa = fc
         else
            b = c
            fb = fc
         end if

         maxit = maxit - 1
      end do

      print*, 'maximum iteration for bisection reached solving D1'

      myzero1slope = c
      return
      end

C------------------------------
C     Target function for Delta2 given other parameters
C------------------------------

      real*8 function target2slope(d2, d1, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim)
      integer xdim, i
      real*8 d2, gamma1(xdim), beta1(xdim), sigma1
      real*8 gamma2(xdim), beta2sp(xdim), sigma21, sigma20, sigma21sp
      real*8 betay, betay0, betaysp, p, tau, x(xdim), d1
      real*8 p1, p2, pnrm

      sigma20 = sigma21 * exp(sigma21sp)
      betay0 = betay + betaysp

      lp1 = dot_product(beta1, x)
      lp2 = dot_product(beta2sp, x)
      quan2 = dot_product(gamma2, x)

      if (betay .ne. 0) then
         p1=pnrm(((-d2+ quan2
     &        )/betay-(d1+lp1))/
     &        sigma1/sqrt(
     &        sigma21**2/
     &        sigma1**2/betay**2+1),
     &        0.d0,1.d0,1,0)
      else
         p1=pnrm((quan2-d2
     &        )/sigma21, 0.d0, 1.d0, 1, 0)
      end if

      if (betay0 .ne. 0) then
         p2=pnrm(((-d2+quan2 - lp2
     &        )/betay0-(d1-lp1))/
     &        sigma1/sqrt(
     &        sigma20**2/
     &        sigma1**2/betay0**2+1),
     &        0.d0,1.d0,1,0)
      else
         p2=pnrm((quan2-d2-lp2
     &        )/sigma20, 0.d0, 1.d0, 1, 0)
      end if

      if (betay < 0) then
         p1 = 1 - p1
      end if

      if (betay0 < 0) then
         p2 = 1 - p2
      end if

      target2slope = tau - p * p1 - (1 - p) * p2

      return
      end

C------------------------------
C     biscetion method to solve for Delta2
C------------------------------

      real*8 function myzero2slope(gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim, d1)
      integer xdim, i
      real*8 d2, gamma1(xdim), beta1(xdim), sigma1
      real*8 gamma2(xdim), beta2sp(xdim), sigma21, sigma20, sigma21sp
      real*8 betay, betay0, betaysp, p, tau, x(xdim), d1
      real*8 p1, p2, pnrm
      real*8 a, b, fa, fb, c, fc, tol, prevstep, newstep
      real*8 target2slope
      logical success
      real*8 dx, t1, cb, t2, pp, q
      integer maxit

      tol = 0.00001
      maxit = 40

      a = -100
      b = 100

      fa = target2slope(a, d1, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim)
      fb = target2slope(b, d1, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim)

C     First test if root is an endpoint
      if (fa .eq. 0.d0) then
         myzero2slope = a
         return
      end if

      if (fb .eq. 0.d0) then
         myzero2slope = b
         return
      end if

      if (fa * fb .ge. 0) then
         print*, 'root is not included for D2.'
         print*, fa, fb, d1, gamma1, beta1, sigma1, gamma2, beta2sp,
     &        sigma21, sigma21sp, betay, betaysp, p, tau,x, xdim
      end if

      do while (maxit .gt. 0)
         c = (a + b)/2.d0

         fc = target2slope(c, d1, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim)

         if (abs(fc) .le. tol) then
            myzero2slope = c
            return
         end if

         if (fa * fc .gt. 0) then
            a = c
            fa = fc
         else
            b = c
            fb = fc
         end if

         maxit = maxit - 1
      end do

      print*, 'maximum iteration for bisection reached for D2'

      myzero2slope = c
      return
      end

C------------------------------
C     mydelta2bise : bisection method to solve delta1, delta2 for all X
C------------------------------

      SUBROUTINE mydelta2biseslope(x, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, n, xdim, delta)
      implicit none
      integer xdim, n
      real*8 x(n, xdim), gamma1(xdim), beta1(xdim),sigma1
      real*8 gamma2(xdim), beta2sp(xdim), sigma21, sigma21sp
      real*8 betay, betaysp
      real*8 p, tau, delta(n, 2)
      real*8 myzero2slope, myzero1slope
      integer i

      do i = 1, n
         delta(i,1) = myzero1slope(gamma1,beta1,sigma1, p,
     &        tau, x(i,:), xdim)
         delta(i,2) = myzero2slope(gamma1, beta1, sigma1,
     &        gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &        p, tau, x(i,:), xdim, delta(i, 1))
      end do
      return
      end
