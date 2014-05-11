C------------------------------
C     update G1, G2
C------------------------------

      subroutine updategmnar(x,gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, n, xdim, delta, K, y, r, g1, g2)
      implicit none
      integer xdim, n, K
      real*8 x(n, xdim), gamma1(xdim), beta1
      real*8 gamma2(xdim), beta2sp
      real*8 mu1(K), sigma1(K), mu2(K), sigma2(K)
      real*8 omega11(K), omega10(K), omega21(K), omega20sp(K)
      real*8 betay, betay0, betaysp
      real*8 p, tau, delta(n, 2)
      real*8 myzero2mix, myzero1mix
      integer i, j
      real*8 y(n, 2)
      integer g1(n), g2(n), r(n)

C     Tmp
      real*8 dd(n, 2), prob1(k), prob2(k), dnrm
      integer rcat

      call mydelta2bisemixmnar(x,gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega11, omega10,
     &     omega21, omega20sp, betay, betaysp,
     &     p, tau, n, xdim, dd, K)

      do i = 1, n
         if (r(i) .eq. 1) then
            do j = 1, k
               prob1(j) = omega11(j)*dnrm(y(i,1), dd(i,1) +
     &              beta1 + mu1(j), sigma1(j), 0)
               prob2(j) = omega21(j)*dnrm(y(i,2), dd(i,2) +
     &              betay*y(i,1) + mu2(j) , sigma2(j), 0)
            end do
         else
            do j = 1, k
               prob1(j) = omega10(j)*dnrm(y(i,1), dd(i,1) -
     &              beta1 + mu1(j), sigma1(j), 0)
               prob2 = omega20sp
            end do
         end if

         g1(i) = rcat(prob1, k)
         g2(i) = rcat(prob2, k)

      end do

      return
      end
