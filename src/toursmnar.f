C------------------------------
C     Target function for Delta2 given other parameters
C------------------------------

      real*8 function targeteqn2fmnar(d2, d1, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim)
      integer xdim, i
      real*8 d2, gamma1(xdim), beta1, sigma1
      real*8 gamma2(xdim), beta2sp, sigma21, sigma20, sigma21sp
      real*8 betay, betay0, betaysp, p, tau, x(xdim), d1
      real*8 myzeroin1
      real*8 p1, p2, pnrm
      real*8 targeteqn2fmnar

      sigma20 = sigma21 * exp(sigma21sp)
      betay0 = 1

      lp1 = beta1
      lp2 = beta2sp
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

      p2=pnrm(((quan2 - lp2
     &     )/betay0-(d1-lp1))/
     &     sigma1/sqrt(
     &     sigma20**2/
     &     sigma1**2/betay0**2+1),
     &     0.d0,1.d0,1,0)

      if (betay < 0) then
         p1 = 1 - p1
      end if

      targeteqn2fmnar = tau - p * p1 - (1 - p) * p2

      return
      end

C------------------------------
C     biscetion method to solve for Delta2
C------------------------------

      real*8 function myzero2mnar(gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim, d1)
      integer xdim, i
      real*8 d2, gamma1(xdim), beta1, sigma1
      real*8 gamma2(xdim), beta2sp, sigma21, sigma20, sigma21sp
      real*8 betay, betay0, betaysp, p, tau, x(xdim), d1
      real*8 p1, p2, pnrm
      real*8 a, b, fa, fb, c, fc, tol, prevstep, newstep
      real*8 targeteqn2fmnar
      logical success
      real*8 dx, t1, cb, t2, pp, q
      integer maxit

      tol = 0.00001
      maxit = 40

      a = -100
      b = 100

      fa = targeteqn2fmnar(a, d1, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim)
      fb = targeteqn2fmnar(b, d1, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim)

C     First test if root is an endpoint
      if (fa .eq. 0.d0) then
         myzero2mnar = a
         return
      end if

      if (fb .eq. 0.d0) then
         myzero2mnar = b
         return
      end if

      if (fa * fb .ge. 0) then
         print*, 'root is not included for D2.'
         print*, fa, fb, d1, gamma1, beta1, sigma1, gamma2, beta2sp,
     &        sigma21, sigma21sp, betay, betaysp, p, tau,x, xdim
      end if

      do while (maxit .gt. 0)
         c = (a + b)/2.d0

         fc = targeteqn2fmnar(c, d1, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim)

         if (abs(fc) .le. tol) then
            myzero2mnar = c
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

      myzero2mnar = c
      return
      end

C------------------------------
C     mydelta2bisemnar : bisection method to solve delta1, delta2 for all X for TOursMNAR
C------------------------------

      SUBROUTINE mydelta2bisemnar(x, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, n, xdim, delta)
      implicit none
      integer xdim, n
      real*8 x(n, xdim), gamma1(xdim), beta1,sigma1
      real*8 gamma2(xdim), beta2sp, sigma21, sigma21sp
      real*8 betay, betaysp
      real*8 p, tau, delta(n, 2)
      real*8 myzero2mnar, myzero1
      integer i

      do i = 1, n
         delta(i,1) = myzero1(gamma1,beta1,sigma1, p,
     &        tau, x(i,:), xdim)
         delta(i,2) = myzero2mnar(gamma1, beta1, sigma1,
     &        gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &        p, tau, x(i,:), xdim, delta(i, 1))
      end do
      return
      end
