      subroutine myldd(x, alpha, k)
      implicit none
      integer k
      real*8 x(k), alpha(k), tmp(k*2)
      real*8 mylddirichlet
      real*8 ans

      print*, size(x)
      ans = mylddirichlet(x, alpha, k)
      print*, tmp
      print*, exp(x(1:3))/(sum(exp(x(1:3))) + 1)

      return
      end


C------------------------------
C     myrcat
C------------------------------

      integer function myrcat(n, prob)

      integer myrcat, n
      real*8 prob(n)

      integer i
      real*8 r, myrunif, cumprob

      r = myrunif(0.d0, 1.d0)
      cumprob = 0

      do i = 1, n
         cumprob = cumprob + prob(i)
         if (r < cumprob) then
            myrcat = i
            return
         end if
      end do

      return
      end

C------------------------------
C     U
C------------------------------

c$$$      real*8 function u(y, x, r, g1, g2,
c$$$     &     gamma1, beta1, gamma2, beta2sp, mu1, sigma1,
c$$$     &     mu2, sigma2, omega1, omega2, betay, pi,
c$$$     &     n, xdim, k, tau)

      real*8 function u(y, x, r, g1, g2, q, n, xdim, k, tau)

      integer n, xdim, k
      integer r(n)
      integer g1(n), g2(n)
      real*8 y(n, 2), x(n, xdim)
      real*8 q(xdim*2 + k*6 + 2)
      real*8 tau

C     internal parameters
      real*8 gamma1(xdim), beta1, gamma2(xdim), beta2sp
      real*8 mu1(k), sigma1(k), mu2(k), sigma2(k)
      real*8 omega1(k), omega2(k), betay, pi

C     internal variables
      real*8 dd(n, 2), ll
      integer i

C     external functions
      real*8 dnrm, u

C     calculate negative loglikelihood

      gamma1 = q(1:xdim)
      beta1 = q(xdim + 1)
      gamma2 = q((xdim + 2): (xdim * 2 + 1))
      beta2sp = q(xdim*2 + 2)
      mu1 = q((xdim*2 + 3):(xdim*2 + K + 2))
      sigma1 = exp(q((xdim*2 + K + 3):(xdim*2 + K*2 + 2)))
      mu2 = q((xdim*2 + K*2+ 3):(xdim*2 + K*3 + 2))
      sigma2 = exp(q((xdim*2 + K*3 + 3):(xdim*2 + K*4 + 2)))
      omega1(1:(K-1)) = exp(q((xdim*2 + K*4 + 3):(xdim*2 +
     &     K*5 + 1)))/(sum(exp(q((xdim*2 + K*4 + 3):(xdim*2 +
     &     K*5 + 1)))) + 1)
      omega1(k) = 1/(sum(exp(q((xdim*2 + K*4 + 3):
     &     (xdim*2 + K*5 + 1)))) + 1)
      omega2(1:(K-1)) = exp(q((xdim*2 + K*5 + 2):(xdim*2 + K*6)))/
     &     (sum(exp(q((xdim*2 + K*5 + 2):(xdim*2 + K*6)))) + 1)
      omega2(K) = 1/(sum(exp(q((xdim*2 + K*5 + 2):
     &     (xdim*2 + K*6)))) + 1)
      betay = q((xdim*2 + K*6 + 1))
      pi = exp(q(xdim*2 + K*6 + 2))/(1 + exp(q(xdim*2 + K*6 + 2)))

c$$$      print*, gamma1, beta1, gamma2, beta2sp, mu1
c$$$      print*, mu2, sigma2, omega1, omega2, betay, pi

      call mydelta2bisemix(x,gamma1, beta1, gamma2,
     &     beta2sp, mu1, sigma1, mu2, sigma2, omega1, omega1,
     &     omega2, omega2, betay, 0.d0,
     &     pi, tau, n, xdim, dd, k)

      ll = 0

      do i = 1, n
         if (r(i) .eq. 1) then
            ll = ll + dnrm(y(i, 1), dd(i, 1) + beta1 + mu1(g1(i)),
     &           sigma1(g1(i)), 1) + log(omega1(g1(i)))
            ll = ll + dnrm(y(i, 2), dd(i, 2) + betay*y(i, 1) +
     &           mu2(g2(i)), sigma2(g2(i)), 1) + log(omega2(g2(i)))
         else
            ll = ll + dnrm(y(i, 1), dd(i, 1) - beta1 + mu1(g1(i)),
     &           sigma1(g1(i)), 1) + log(omega1(g1(i)))
         end if
      end do

      u = -ll

      return
      end


C------------------------------
C     dU
C------------------------------

      subroutine du(y, x, r, g1, g2, q, n, xdim, k, tau, mydu)

      integer n, xdim, k
      integer r(n)
      integer g1(n), g2(n)
      real*8 y(n, 2), x(n, xdim)
      real*8 q(xdim*2 + k*6 + 2)
      real*8 tau
      real*8 mydu(xdim*2 + k*6 + 2)

C     internal parameters
      real*8 gamma1(xdim), beta1, gamma2(xdim), beta2sp
      real*8 mu1(k), sigma1(k), mu2(k), sigma2(k)
      real*8 omega1(k), omega2(k), betay, pi

C     internal variables
      real*8 dd(n, 2), ll
      integer i, l, j
      real*8 e, old, qc(xdim*2 + k*6 + 2), new

C     external functions
      real*8 dnrm, u

      e = 0.003
      old = u(y, x, r, g1, g2, q, n, xdim, k, tau)
      l = xdim*2 + k*6 + 2

      do i = 1, l
         do j = 1, l
            qc(j) = q(j)
         end do
         qc(i) = qc(i) + e
         new = u(y, x, r, g1, g2, qc, n, xdim, k, tau)
         mydu(i) = (new - old)/e
      end do

      return
      end



C------------------------------
C     HMC
C------------------------------

      subroutine hmc(y, x, r, n, xdim, k, tau,
     &     nscan, nburn, ndisp, nskip, nsave,
     &     initial, qsave)

      implicit none
      integer n, xdim, k
      integer r(n)
      real*8 y(n, 2), x(n, xdim), tau
      integer nscan, nburn, ndisp, nskip, nsave
      real*8 qsave(nsave, 2*xdim + k*6 + 2)
      real*8 initial(2*xdim + k*6 + 2)

C     internal
      integer i, j, lenq, iscan, l
      real*8 newq(2*xdim + k*6 + 2)
      real*8 q(2*xdim + k*6 + 2)
      real*8 p(2*xdim + k*6 + 2)
      real*8 currentp(2*xdim + k*6 + 2)
      real*8 epsilon
      integer g1(n), g2(n)
      real*8 currentu, currentk, proposedu, proposedk
      real*8 prob1(k), prob2(k)
      integer skipcount, dispcount, isave
      real*8 tmpdu(2*xdim + k*6 + 2), tmp

C     internal parameters
      real*8 gamma1(xdim), beta1, gamma2(xdim), beta2sp
      real*8 mu1(k), sigma1(k), mu2(k), sigma2(k)
      real*8 omega1(k), omega2(k), betay, pi, dd(n, 2)


C     function
      real*8 u, myrnorm, myrunif, myrcat, dnrm

C     initial
      lenq = 2*xdim + k*6 + 2
      epsilon = 0.01
      l = 15
      do i = 1, n
         g1(i) = 1
         g2(i) = 1
      end do
      skipcount = 0
      dispcount = 0

      do i = 1, lenq
         q(i) = initial(i)
      end do

C     start roll
      do iscan = 1, nscan

         do i = 1, lenq
            newq(i) = q(i)
            p(i) = myrnorm(0.d0, 1.d0)
            currentp(i) = p(i)
         end do

         call du(y, x, r, g1, g2,
     &        q, n, xdim, k, tau, tmpdu)

         do i = 1, lenq
            p(i) = p(i) - epsilon * tmpdu(i)/2
         end do
c$$$         print*, 'after first u'
         do i = 1, l
            do j = 1, lenq
               newq(j) = newq(j) + epsilon * p(j)
            end do
            call du(y, x, r, g1, g2,
     &           newq, n, xdim, k, tau, tmpdu)
            if (i .ne. l) then
               do j = 1, lenq
                  p(j) = p(j) - epsilon * tmpdu(j)
               end do
            end if
         end do
c$$$         print*, 'after first du'
         call du(y, x, r, g1, g2,
     &        newq, n, xdim, k, tau, tmpdu)
         do i = 1, lenq
            p(i) = p(i) - epsilon * tmpdu(i)/2
            p(i) = -p(i)
         end do

         currentu = u(y, x, r, g1, g2,
     &        q, n, xdim, k, tau)/2
         proposedu = u(y, x, r, g1, g2,
     &        newq, n, xdim, k, tau)/2
         currentk = sum(currentp ** 2)/2
         proposedk = sum(p ** 2)/2

C     accept
         if (myrunif(0.d0, 1.d0) < exp(currentu - proposedu +
     &        currentk - proposedk)) then
            do i = 1, lenq
               q(i) = newq(i)
            end do
         end if
         print*, 'before update'
C     UPDATE G1, G2
         gamma1 = q(1:xdim)
         beta1 = q(xdim + 1)
         gamma2 = q((xdim + 2): (xdim * 2 + 1))
         beta2sp = q(xdim*2 + 2)
         mu1 = q((xdim*2 + 3):(xdim*2 + K + 2))
         sigma1 = exp(q((xdim*2 + K + 3):(xdim*2 + K*2 + 2)))
         mu2 = q((xdim*2 + K*2+ 3):(xdim*2 + K*3 + 2))
         sigma2 = exp(q((xdim*2 + K*3 + 3):(xdim*2 + K*4 + 2)))
         omega1(1:(K-1)) = exp(q((xdim*2 + K*4 + 3):(xdim*2 +
     &        K*5 + 1)))/(sum(exp(q((xdim*2 + K*4 + 3):(xdim*2 +
     &        K*5 + 1)))) + 1)
         omega1(k) = 1/(sum(exp(q((xdim*2 + K*4 + 3):
     &        (xdim*2 + K*5 + 1)))) + 1)
         omega2(1:(K-1)) = exp(q((xdim*2 + K*5 + 2):(xdim*2 + K*6)))/
     &        (sum(exp(q((xdim*2 + K*5 + 2):(xdim*2 + K*6)))) + 1)
         omega2(K) = 1/(sum(exp(q((xdim*2 + K*5 + 2):
     &        (xdim*2 + K*6)))) + 1)
         betay = q((xdim*2 + K*6 + 1))
         pi = exp(q(xdim*2 + K*6 + 2))/(1 + exp(q(xdim*2 + K*6 + 2)))

         call mydelta2bisemix(x,gamma1, beta1, gamma2,
     &        beta2sp, mu1, sigma1, mu2, sigma2, omega1, omega1,
     &        omega2, omega2, betay, 0.d0,
     &        pi, tau, n, xdim, dd, k)

         do i = 1, n
            if (r(i) .eq. 1) then
               do j = 1, k
                  prob1(j) = omega1(j) * dnrm(y(i, 1),
     &                 dd(i, 1) + beta1 + mu1(j),
     &                 sigma1(j), 0)
                  prob2(j) = omega2(j) * dnrm(y(i, 2),
     &                 dd(i, 2) + betay*y(i, 1) +
     &                 mu2(j), sigma2(j), 0)
               end do
            else
               do j = 1, k
                  prob1(j) = omega1(j) * dnrm(y(i, 1),
     &                 dd(i, 1) - beta1 + mu1(j),
     &                 sigma1(j), 0)
                  prob2(j) = omega2(j)
               end do
               g1(i) = myrcat(k, prob1)
               g2(i) = myrcat(k, prob2)
            end if
         end do

C     save

         if (iscan .gt. nburn) then
            skipcount = skipcount + 1
            if (skipcount .ge. nskip) then
               isave = isave + 1
               dispcount = dispcount + 1
               do i = 1, lenq
                  qsave(isave, i) = q(i)
               end do
               skipcount = 0
               if (dispcount .ge. ndisp) then
                  dispcount = 0
                  print*, isave, 'out of', nsave, '\n'
               end if
            end if
         end if

      end do

      return
      end
