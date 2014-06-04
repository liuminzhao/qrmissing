C------------------------------
C     TARGET DELTA EQUATION for univariate  with mixture normals
C------------------------------

      real*8 function target1mixslope(delta1,gamma1,beta1,
     &     G, mu, sigma, omega1, omega0, p,tau,x,xdim)
      integer xdim, i, G
      real*8 delta1, gamma1(xdim), beta1(xdim), mu(G), sigma(G),
     &     omega1(G), omega0(G), p, tau, x(xdim)
      real*8 pnrm

      real*8 quan, lp

      quan = dot_product(gamma1, x)
      lp = dot_product(beta1, x)

      target1mixslope = 0

      do i = 1, G
         target1mixslope = target1mixslope+p*omega1(i)*pnrm(quan-delta1
     &        - lp - mu(i),0.d0,sigma(i),1,0) +
     &     (1-p)*omega0(i)*pnrm(quan-delta1+lp-mu(i),0.d0,sigma(i),1,0)
      end do
      target1mixslope = tau - target1mixslope

      return
      end

C------------------------------
C     myzero1mix: bisection method to solve for delta1 for mixture normal
C------------------------------
      real*8 function myzero1mixslope(gamma1,beta1, G, mu, sigma,
     &     omega1, omega0, p,tau, x, xdim)
      integer xdim, i, G
      real*8 delta1, gamma1(xdim), beta1(xdim), mu(G), sigma(G),
     &     omega1(G), omega0(G), p, tau, x(xdim)
      real*8 target1mixslope

      real*8 a, b, fa, fb, c, fc, tol, prevstep, newstep
      logical success
      real*8 dx, t1, cb, t2, pp, q
      integer maxit

      tol = 0.00001
      maxit = 40

      a = -100
      b = 100

      fa = target1mixslope(a,gamma1,beta1,
     &     G, mu, sigma, omega1, omega0, p,tau,x,xdim)
      fb = target1mixslope(b,gamma1,beta1,
     &     G, mu, sigma, omega1, omega0, p,tau,x,xdim)

      c = a
      fc = fa

C     First test if root is an endpoint
      if (abs(fa) .le. tol) then
         myzero1mixslope = a
         return
      end if

      if (abs(fb) .le. tol) then
         myzero1mixslope = b
         return
      end if

      if (fa * fb .ge. 0) print*, 'root is not included solving D1.'

      do while (maxit .gt. 0)
         c = (a + b)/2.d0

         fc = target1mixslope(c,gamma1,beta1,
     &     G, mu, sigma, omega1, omega0, p,tau,x,xdim)

         if (abs(fc) .le. tol) then
            myzero1mixslope = c
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

      print*, 'maximum iteration for bisection reached solving D1 mix'

      myzero1mixslope = c

      return
      end


C------------------------------
C     target2mixslope
C------------------------------

      real*8 function target2mixslope(d2, d1, gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, x, xdim, K)
      integer xdim, i, K, j
      real*8 d2, d1, gamma1(xdim), beta1(xdim)
      real*8 gamma2(xdim), beta2sp(xdim)
      real*8 mu1(K), sigma1(K), mu2(K), sigma2(K)
      real*8 omega11(K), omega10(K), omega21(K), omega20sp(K)
      real*8 betay, betay0, betaysp, p, tau, x(xdim)
      real*8 p1, p2, pnrm, ans1, ans2, tol

      tol = 0.00001

      betay0 = betay + betaysp

      lp1 = dot_product(beta1, x)
      lp2 = dot_product(beta2sp, x)
      quan2 = dot_product(gamma2, x)

      ans1 = 0.d0
      ans2 = 0.d0

      do i = 1, K
         do j = 1, K
            if (abs(betay) .ge. tol) then
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

            if (abs(betay0) .ge. tol) then
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

      target2mixslope = tau - p * ans1 - (1 - p) * ans2

      return
      end

C------------------------------
C     biscetion method to solve for Delta2 for mixture of normals
C------------------------------

      real*8 function myzero2mixslope(gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, x, xdim, d1, K)
      integer xdim, i, K
      real*8 d2, gamma1(xdim), beta1(xdim)
      real*8 gamma2(xdim), beta2sp(xdim)
      real*8 mu1(K), sigma1(K), mu2(K), sigma2(K)
      real*8 omega11(K), omega10(K), omega21(K), omega20sp(K)
      real*8 betay, betay0, betaysp, p, tau, x(xdim), d1
      real*8 p1, p2, pnrm
      real*8 a, b, fa, fb, c, fc, tol, prevstep, newstep
      real*8 target2mixslope
      logical success
      real*8 dx, t1, cb, t2, pp, q
      integer maxit

      tol = 0.00001
      maxit = 40

      a = -100
      b = 100

      fa = target2mixslope(a, d1, gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, x, xdim, K)
      fb = target2mixslope(b, d1, gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, x, xdim, K)

C     First test if root is an endpoint
      if (fa .eq. 0.d0) then
         myzero2mixslope = a
         return
      end if

      if (fb .eq. 0.d0) then
         myzero2mixslope = b
         return
      end if

      if (fa * fb .ge. 0) then
         print*, 'root is not included for D2.'
c$$$         print*, fa, fb, d1, gamma1, beta1, sigma1, gamma2, beta2sp,
c$$$     &        sigma2, betay, betaysp, p, tau,x, xdim, omega11,
c$$$     &        omega10, omega21, omega20sp, mu1, mu2
      end if

      do while (maxit .gt. 0)
         c = (a + b)/2.d0

         fc = target2mixslope(c, d1, gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, x, xdim, K)

         if (abs(fc) .le. tol) then
            myzero2mixslope = c
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

      myzero2mixslope = c
      return
      end


C------------------------------
C     mydelta2bisemix : bisection method to solve delta1, delta2 for all X
C     for mixture of normals bivariate case
C------------------------------

      SUBROUTINE mydelta2bisemixslope(x,gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, n, xdim, delta, K)
      implicit none
      integer xdim, n, K
      real*8 x(n, xdim), gamma1(xdim), beta1(xdim)
      real*8 gamma2(xdim), beta2sp(xdim)
      real*8 mu1(K), sigma1(K), mu2(K), sigma2(K)
      real*8 omega11(K), omega10(K), omega21(K), omega20sp(K)
      real*8 betay, betay0, betaysp
      real*8 p, tau, delta(n, 2)
      real*8 myzero2mixslope, myzero1mixslope
      integer i

      do i = 1, n
         delta(i,1) = myzero1mixslope(gamma1,beta1, K, mu1, sigma1,
     &        omega11, omega10, p,tau, x(i, :), xdim)
         delta(i,2) = myzero2mixslope(gamma1, beta1, gamma2,
     &        beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &        omega21, omega20sp, betay, betaysp,
     &        p, tau, x(i,:), xdim, delta(i, 1), K)

      end do
      return
      end
