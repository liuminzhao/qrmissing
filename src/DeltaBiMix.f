c===========================================================
c$$$
C$$$  Time-stamp: <liuminzhao 03/25/2014 14:26:21>
c$$$  2014/01/16 Bayesian MCMC for QRMissing Bivariate mixture of normals
c$$$  make use of targetunimix, myzero1mix, mydelta1bisemix function
c===========================================================

C------------------------------
C     Target function for Delta2 given other parameters
C     q = (gamma1(xdim), beta1(xdim), gamma2(xdim), beta2sp(xdim),
C     mu1(K), sigma1(k), mu2(K), sigma2(K), omega11(K), omega10(K),
C     omega21(K), omega20sp(K), betay(1), betaysp(1), pi)
C     external: tau, y(n, 2), X(n, p), R(n), K
C     temporery: G1(n), G2(n)
C
C     targetbimix has nothing to do with y, R, G1, G2
C------------------------------

      real*8 function targetbimix(d2, d1, gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, x, xdim, K)
      integer xdim, i, K, j
      real*8 d2, d1, gamma1(xdim), beta1
      real*8 gamma2(xdim), beta2sp
      real*8 mu1(K), sigma1(K), mu2(K), sigma2(K)
      real*8 omega11(K), omega10(K), omega21(K), omega20sp(K)
      real*8 betay, betay0, betaysp, p, tau, x(xdim)
      real*8 myzero1mix
      real*8 p1, p2, pnrm, ans1, ans2

      betay0 = betay + betaysp

      lp1 = beta1
      lp2 = beta2sp
      quan2 = dot_product(gamma2, x)

      ans1 = 0.d0
      ans2 = 0.d0

      do i = 1, K
         do j = 1, K
            if (betay .ne. 0) then
               p1=pnrm(((-d2+ quan2 - mu2(j)
     &              )/betay-(d1+lp1 + mu1(i)))/
     &              sqrt(
     &              sigma2(j)**2/betay**2+ sigma1(i)**2),
     &              0.d0,1.d0,1,0)
               if (betay .le. 0) p1 = 1 - p1
            else
               p1=pnrm((quan2-d2 - mu2(j)
     &              )/sigma2(j), 0.d0, 1.d0, 1, 0)
            end if
            ans1 = ans1 + omega11(i)*omega21(j)*p1

            if (betay0 .ne. 0) then
               p2=pnrm(((-d2+quan2 - lp2 - mu2(j)
     &              )/betay0-(d1-lp1+mu1(i)))/
     &              sqrt(
     &              sigma2(j)**2/betay0**2+sigma1(i)**2),
     &              0.d0,1.d0,1,0)
               if (betay0 .le. 0) p2 = 1 - p2
            else
               p2=pnrm((quan2-d2-lp2 - mu2(j)
     &              )/sigma2(j), 0.d0, 1.d0, 1, 0)
            end if

            ans2 = ans2 + omega10(i)*omega20sp(j)*p2

         end do
      end do

      targetbimix = tau - p * p1 - (1 - p) * p2

      return
      end


C------------------------------
C     biscetion method to solve for Delta2 for mixture of normals
C------------------------------

      real*8 function myzero2mix(gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, x, xdim, d1, K)
      integer xdim, i, K
      real*8 d2, gamma1(xdim), beta1
      real*8 gamma2(xdim), beta2sp
      real*8 mu1(K), sigma1(K), mu2(K), sigma2(K)
      real*8 omega11(K), omega10(K), omega21(K), omega20sp(K)
      real*8 betay, betay0, betaysp, p, tau, x(xdim), d1
      real*8 p1, p2, pnrm
      real*8 a, b, fa, fb, c, fc, tol, prevstep, newstep
      real*8 targetbimix, myzero2mix
      logical success
      real*8 dx, t1, cb, t2, pp, q
      integer maxit

      tol = 0.00001
      maxit = 40

      a = -100
      b = 100

      fa = targetbimix(a, d1, gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, x, xdim, K)
      fb = targetbimix(b, d1, gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, x, xdim, K)

C     First test if root is an endpoint
      if (fa .eq. 0.d0) then
         myzero2mix = a
         return
      end if

      if (fb .eq. 0.d0) then
         myzero2mix = b
         return
      end if

      if (fa * fb .ge. 0) then
         print*, 'root is not included for D2.'
         print*, fa, fb, d1, gamma1, beta1, sigma1, gamma2, beta2sp,
     &        sigma21, sigma21sp, betay, betaysp, p, tau,x, xdim
      end if

      do while (maxit .gt. 0)
         c = (a + b)/2.d0

         fc = targetbimix(c, d1, gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, x, xdim, K)

         if (abs(fc) .le. tol) then
            myzero2mix = c
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

      myzero2mix = c
      return
      end

C------------------------------
C     mydelta2bisemix : bisection method to solve delta1, delta2 for all X
C     for mixture of normals bivariate case
C------------------------------

      SUBROUTINE mydelta2bisemix(x,gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, n, xdim, delta, K)
      implicit none
      integer xdim, n, K
      real*8 x(n, xdim), gamma1(xdim), beta1
      real*8 gamma2(xdim), beta2sp
      real*8 mu1(K), sigma1(K), mu2(K), sigma2(K)
      real*8 omega11(K), omega10(K), omega21(K), omega20sp(K)
      real*8 betay, betay0, betaysp
      real*8 p, tau, delta(n, 2)
      real*8 myzero2mix, myzero1mix
      integer i

      do i = 1, n
         delta(i,1) = myzero1mix(gamma1,beta1, K, mu1, sigma1,
     &        omega11, omega10, p,tau, x(i, :), xdim)
         delta(i,2) = myzero2mix(gamma1, beta1, gamma2,
     &        beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &        omega21, omega20sp, betay, betaysp,
     &        p, tau, x(i,:), xdim, delta(i, 1), K)

      end do
      return
      end
