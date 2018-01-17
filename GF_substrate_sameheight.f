      subroutine gf_substrate_sameheight(height, rho, dist, xk, eps_sub,
     $                        GRxx, GRyy, GRzz, GRxz, GRzx, istat) 
                         ! Numerical calculation of the half-space
                         ! electromagnetic Green's tensor
                         ! by discretization of the Sommerfeld integrals
                         !
                         ! This program uses intrinsic Bessel functions
c--------------------------------------------------------------------------
      implicit none         

      complex*16, intent(out):: GRxx, GRyy, GRzz, GRxz, GRzx
      real*8    , intent(in) :: rho, height, xk, dist
      integer*4 , intent(out):: istat
      complex*16, intent(in) :: eps_sub

      integer*8 N, N1, N2, N3, k
!      parameter(N=1e5, N1=2*(N/2), N2=2*(N/2), N3=2*(N/2))
      real*8 y_max, div
      parameter(div = 10.0d0)

      real*8 ru, rz, two, three, four, half, pi, twopi
      real*8 w, J0, J1, J0Hat, J1Hat
      real*8 y, y2, y3, y4, dy1, dy2, dy3
      real*8 Rhat, Zhat_p, z_p
      real*8 rpe, r_p
      real*8 zh_p, zh_p_2
      real*8 RJ01p, RJ0p, RJ1p
c$$$      real*8 DBESJ0, DBESJ1

      complex*16 cu, cz, eps, eps1, kzeps
      complex*16 cpe, R, Rpr
      complex*16 J01p, J0p, J1p
      complex*16 s1xx,s1yy,s1xz,s1zz,s2xx,s2yy,s2xz,s2zz
      complex*16 sum1xz,sum2xz
      complex*16 sum1xx_a,sum1xx_b,sum1xx_c
      complex*16 sum1yy_a,sum1yy_b,sum1yy_c
      complex*16 sum1zz_a,sum1zz_b,sum1zz_c
      complex*16 sum2xx_a,sum2xx_b,sum2xx_c
      complex*16 sum2yy_a,sum2yy_b,sum2yy_c
      complex*16 sum2zz_a,sum2zz_b,sum2zz_c
      complex*16 Gxx1a,Gxx1b,Gxx2a1,Gxx2a2,Gxx2ah,Gxx2bh,Gxx2b1,Gxx2b2
      complex*16 Gyy1a,Gyy1b,Gyy2a1,Gyy2a2,Gyy2ah,Gyy2bh,Gyy2b1,Gyy2b2
      complex*16 Gzz1,Gzz2b1,Gzz2b2,Gzz2
      complex*16 Gxz1,Gxz2a1,Gxz2a2,Gxz2
c
c------ Mathematical constants ------------------------------------
c
      rz    = 0.0d0
      ru    = 1.0d0
      cz    = (0.0d0,0.0d0)
      cu    = (0.0d0,1.0d0)
      two   = 2.0d0
      three = 3.0d0
      four  = 4.0d0
      half  = 0.5d0
      pi    = two*dasin(1.0d0)
      twopi = two*pi
      
      istat = 0

      if(height .le. rz) then
         write(*,*) 'GF: incorrect parameter [z<=0]; z=', height
         istat = -2
         return
      end if

!      if(rho .le. rz) then
!         write(*,*) 'GF: incorrect parameter: [rho<=0]; rho=', rho
!         istat = -3
!         return
!      end if

!      if(eps_sub .lt. ru) then
!         write(*,*) 'GF: incorrect parameter [eps_sub<1]; eps_sub=', 
!     $                                                    eps_sub
!         istat = -4
!         return
!      end if

!      if(N .lt. 10) then
!         write(*,*) 'GF: incorrect parameter: [N<10]; N=', N
!         istat = -5
!         return
!      end if

!      if(y_max .le. rz) then
!         write(*,*) 'GF: incorrect parameter: [y_max<=0]; y_max=', 
!     $                                                    y_max
!         istat = -6
!         return
!      end if

c------------------------------------------------------- 
c--------------DIELECTRIC SUBSTRATE--------------------- 
c------------------------------------------------------- 
 
!      IF (eps_sub .eq. (2.25,0.01)) THEN 
!      
      if (rho .ge. rz .and. rho .le. dist/0.6d0) then !40
        N1 = 1.0d2
        N2 = 1.0d4
        N3 = 1.0d4
        y_max = 5.0d2
      end if
      
      if (rho .gt. dist/0.6d0 .and. rho .le. dist/0.06d0) then !400
        N1 = 1.0d3
        N2 = 1.0d4
        N3 = 1.0d5
        y_max = 5.0d2
      end if
      
      if (rho .gt. dist/0.06d0 .and. rho .le. dist/0.024d0) then !1000
        N1 = 3.0d3
        N2 = 1.0d4
        N3 = 3.0d5
        y_max = 3.0d2
      end if
      
      if (rho .gt. dist/0.024d0 .and. rho .le. dist*1.0d2) then !2400
        N1 = 2.0d4
        N2 = 2.0d4
        N3 = 5.0d5
        y_max = 3.0d2
      end if
      
      if (rho .gt. dist*1.0d2 .and. rho .le. dist*5.0d2) then !12000
        N1 = 5.0d5
        N2 = 1.0d5
        N3 = 3.0d6
        y_max = 3.0d2
      end if
      
      if (rho .gt. dist*5.0d2 .and. rho .le. dist*1.0d3) then !24000
        N1 = 2.0d6
        N2 = 2.0d5
        N3 = 1.0d6
        y_max = 2.5d2
      end if
!      
!      END IF
      
c------------------------------------------------------- 
c------------------METAL SUBSTRATE---------------------- 
c------------------------------------------------------- 

      if (rho .ge. rz .and. rho .le. dist/0.6d0) then !40
        N1 = 1.0d3
        N2 = 1.0d5
        N3 = 1.0d5
        y_max = 5.0d2
      end if
      
      if (rho .gt. dist/0.6d0 .and. rho .le. dist/0.06d0) then !400
        N1 = 1.0d4
        N2 = 1.0d5
        N3 = 1.0d6
        y_max = 5.0d2
      end if
      
      if (rho .gt. dist/0.06d0 .and. rho .le. dist/0.024d0) then !1000
        N1 = 3.0d4
        N2 = 1.0d5
        N3 = 3.0d6
        y_max = 3.0d2
      end if
      
      if (rho .gt. dist/0.024d0 .and. rho .le. dist*1.0d2) then !2400
        N1 = 2.0d5
        N2 = 2.0d5
        N3 = 5.0d6
        y_max = 3.0d2
      end if
      
      if (rho .gt. dist*1.0d2 .and. rho .le. dist*5.0d2) then !12000
        N1 = 5.0d6
        N2 = 1.0d6
        N3 = 3.0d7
        y_max = 3.0d2
      end if
      
      if (rho .gt. dist*5.0d2 .and. rho .le. dist*1.0d3) then !24000
        N1 = 2.0d7
        N2 = 2.0d6
        N3 = 1.0d7
        y_max = 2.5d2
      end if
c
c------- Input-dependent constants --------------------------
c
      z_p = two*height
      r_p = dsqrt(rho*rho + z_p*z_p)
      zh_p   = z_p/r_p
      zh_p_2 = zh_p*zh_p

      dy1 = ru/dfloat(N1)
      dy2 = div/dfloat(N2)
      dy3 = (y_max-div) / dfloat(N3)
      
      eps  = eps_sub
      eps1 = eps - ru

      Zhat_p = xk*z_p       ! Used for G^R
      Rhat   = xk*rho     

C$$$  J0Hat = DBESJ0(Rhat)
C$$$  J1Hat = DBESJ1(Rhat)/Rhat
      CALL BJ01(Rhat,J0Hat,J1Hat)
c
c-------- End input-dependent constants block ------------------------
c

c
c-------- Integrals over the finite interval (0,1) ---------
c                 dy = dy1 = 1/N1   

      s2xx = -J1Hat 
      s2yy =  rz    
      s2zz =  J0Hat 
      s2xz =  rz    

      kzeps = cdsqrt(eps1 + ru)
      R   = (ru-kzeps) /(ru+kzeps)
      Rpr = (kzeps-eps)/(kzeps+eps)
      cpe = half*cdexp(cu*Zhat_p)
      s1xx = s2xx + cpe*R    
      s1yy = s2yy + cpe*Rpr  
      s1zz = s2zz             
      s1xz = s2xz             

      sum1xx_a = cz
      sum1yy_a = cz
      sum1zz_a = cz
      sum1xx_b = cz
      sum1yy_b = cz
      sum1xz   = cz
      do 11,k=1,N1-1,2
         y     = k*dy1
         y2    = y*y
         y3    = ru-y2
         y4    = y*y3
         kzeps = cdsqrt(eps1 + y2)
         R     = (y-kzeps)/(y+kzeps)
         Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
         cpe   = cdexp(cu*y*Zhat_p)
         w     = Rhat*dsqrt(y3)
c$$$        J0 = DBESJ0(w)
c$$$        J1 = DBESJ1(w)/w
            CALL BJ01(W,J0,J1)
         J01p  = (J0-J1)*cpe
         J1p   = J1*cpe
         J0p   = J0*cpe
         sum1xx_a = sum1xx_a + J1p  * R
         sum1yy_a = sum1yy_a + J1p  * Rpr * y2
         sum1zz_a = sum1zz_a + J0p  * Rpr * y3
         sum1xx_b = sum1xx_b + J01p * Rpr * y2
         sum1yy_b = sum1yy_b + J01p * R
         sum1xz   = sum1xz   + J1p  * Rpr * y4
 11   continue

      sum2xx_a = cz
      sum2yy_a = cz
      sum2zz_a = cz
      sum2xx_b = cz
      sum2yy_b = cz
      sum2xz   = cz
      do 12, k=2,N1-2,2
         y     = k*dy1
         y2    = y*y
         y3    = ru-y2
         y4    = y*y3
         kzeps = cdsqrt(eps1 + y2)
         R     = (y-kzeps)/(y+kzeps)
         Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
         cpe   = cdexp(cu*y*Zhat_p)
         w     = Rhat*dsqrt(y3)
c$$         J0 = DBESJ0(w)
c$$         J1 = DBESJ1(w)/w
            CALL BJ01(W,J0,J1)
         J01p  = (J0-J1)*cpe
         J1p   = J1*cpe
         J0p   = J0*cpe
         sum2xx_a = sum2xx_a + J1p  * R             
         sum2yy_a = sum2yy_a + J1p  * Rpr * y2
         sum2zz_a = sum2zz_a + J0p  * Rpr * y3
         sum2xx_b = sum2xx_b + J01p * Rpr * y2
         sum2yy_b = sum2yy_b + J01p * R
         sum2xz   = sum2xz   + J1p  * Rpr * y4
 12   continue

      Gxx1a = dy1*(two*sum2xx_a + four*sum1xx_a + s1xx)/three
      Gyy1a = dy1*(two*sum2yy_a + four*sum1yy_a + s1yy)/three
      Gzz1  = dy1*(two*sum2zz_a + four*sum1zz_a + s1zz)/three
      Gxz1  = dy1*(two*sum2xz   + four*sum1xz   + s1xz)/three
        
      s2xx = rz
      s2yy = J1Hat-J0Hat

      kzeps = cdsqrt(eps1 + ru)
      Rpr   = (kzeps-eps)/(kzeps+eps)
      R     = (ru-kzeps)/(ru+kzeps)
      cpe   = half*cdexp(cu*Zhat_p)
      s1xx  = s2xx + cpe*Rpr
      s1yy  = s2yy + cpe*R

      Gxx1b = dy1*(two*sum2xx_b + four*sum1xx_b + s1xx)/three
      Gyy1b = dy1*(two*sum2yy_b + four*sum1yy_b + s1yy)/three
c
c-------- Integrals over the infinite interval (0,Inf) ---------
c-------- The integlas is broken in two parts ------------------
c-------- a) First, compute over the interval (0,div)
c                 dy = dy2 = div/N2
c
      s2xx = -J1Hat
      s2yy =  rz
      s2xz =  rz

      y     = div
      y2    = y*y
      rpe   = dexp(-y*Zhat_p)
      kzeps = cdsqrt(y2-eps1)
      R     = (y-kzeps)/(y+kzeps)
      Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
      w     = Rhat*dsqrt(ru+y2)
c$$$     J1 = DBESJ1(w)/w
         CALL BJ01(W,J0,J1)
      s1xx  = s2xx + rpe*J1*R
      s1yy  = s2yy + y2*rpe*J1*Rpr
      s1xz  = s2xz + y*(y2+ru)*rpe*J1*Rpr

      sum1xx_a = cz
      sum1yy_a = cz
      sum1xx_b = cz
      sum1yy_b = cz
      sum1zz_b = cz
      sum1xz   = cz
      do 21, k = 1, N2-1, 2
         y      = k*dy2
         y2     = y*y
         y3     = ru+y2
         y4     = y*y3
         rpe    = dexp(-y*Zhat_p)
         kzeps  = cdsqrt(y2-eps1)
         R      = (y-kzeps)/(y+kzeps)
         Rpr    = (kzeps-eps*y)/(kzeps+eps*y)
         w      = Rhat*dsqrt(y3)
c$$$         J0 = DBESJ0(w)
c$$$         J1 = DBESJ1(w)/w
             CALL BJ01(W,J0,J1)
         RJ01p  = (J0-J1)*rpe
         RJ1p   = J1*rpe
         RJ0p   = J0*rpe
         sum1xx_a = sum1xx_a + RJ1p  * R
         sum1yy_a = sum1yy_a + RJ1p  * Rpr * y2
         sum1xx_b = sum1xx_b + RJ01p * Rpr * y2
         sum1yy_b = sum1yy_b + RJ01p * R
         sum1zz_b = sum1zz_b + RJ0p  * Rpr * y3
         sum1xz   = sum1xz   + RJ1p  * Rpr * y4 
 21   continue

      sum2xx_a = cz
      sum2yy_a = cz
      sum2xx_b = cz
      sum2yy_b = cz
      sum2zz_b = cz
      sum2xz   = cz
      do 22, k = 2, N2-2, 2
         y      = k*dy2
         y2     = y*y
         y3     = ru+y2
         y4     = y*y3
         rpe    = dexp(-y*Zhat_p)
         kzeps  = cdsqrt(y2-eps1)
         R      = (y-kzeps)/(y+kzeps)
         Rpr    = (kzeps-eps*y)/(kzeps+eps*y)
         w      = Rhat*dsqrt(y3)
c$$$         J0 = DBESJ0(w)
c$$$         J1 = DBESJ1(w)/w
             CALL BJ01(W,J0,J1)
         RJ01p  = (J0-J1)*rpe
         RJ1p   = J1*rpe
         RJ0p   = J0*rpe
         sum2xx_a = sum2xx_a + RJ1p  * R
         sum2yy_a = sum2yy_a + RJ1p  * Rpr * y2
         sum2xx_b = sum2xx_b + RJ01p * Rpr * y2
         sum2yy_b = sum2yy_b + RJ01p * R
         sum2zz_b = sum2zz_b + RJ0p  * Rpr * y3
         sum2xz   = sum2xz   + RJ1p  * Rpr * y4
 22   continue

      Gxx2a1 = dy2*(two*sum2xx_a + four*sum1xx_a + s1xx)/three
      Gyy2a1 = dy2*(two*sum2yy_a + four*sum1yy_a + s1yy)/three
      Gxz2a1 = dy2*(two*sum2xz   + four*sum1xz   + s1xz)/three
c
c-------- Integrals over the infinite interval (0,Inf) ---------
c-------- The integlas is broken in two parts ------------------
c-------- b) Second part, compute over the interval (div,y_max)
c                 dy = dy3 = (y_max-div)/N3
c
      y     = div
      y2    = y*y
      rpe   = dexp(-y*Zhat_p)
      kzeps = cdsqrt(y2-eps1)
      R     = (y-kzeps)/(y+kzeps)
      Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
      w     = Rhat*dsqrt(ru+y2)
c$$$     J1 = DBESJ1(w)/w
         CALL BJ01(W,J0,J1)
      s2xx  = rpe*J1*R
      s2yy  = y2*rpe*J1*Rpr
      s2xz  = y*(y2+ru)*rpe*J1*Rpr

      y     = y_max
      y2    = y*y
      rpe   = dexp(-y*Zhat_p)
      kzeps = cdsqrt(y2-eps1)
      R     = (y-kzeps)/(y+kzeps)
      Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
      w     = Rhat*dsqrt(ru+y2)
c$$$     J1 = DBESJ1(w)/w
         CALL BJ01(W,J0,J1)
      s1xx  = s2xx + rpe*J1*R
      s1yy  = s2yy + y2*rpe*J1*Rpr
      s1xz  = s2xz + y*(y2+ru)*rpe*J1*Rpr

      sum1xx_a = cz
      sum1yy_a = cz
      sum1xz   = cz
      sum1xx_c = cz
      sum1yy_c = cz
      sum1zz_c = cz
      do 31, k = 1, N3-1, 2
         y      = div + k*dy3
         y2     = y*y
         y3     = ru+y2
         y4     = y*y3
         rpe    = dexp(-y*Zhat_p)
         kzeps  = cdsqrt(y2-eps1)
         R      = (y-kzeps)/(y+kzeps)
         Rpr    = (kzeps-eps*y)/(kzeps+eps*y)
         w      = Rhat*dsqrt(y3)
c$$$         J0 = DBESJ0(w)
c$$$         J1 = DBESJ1(w)/w
             CALL BJ01(W,J0,J1)
         RJ01p  = (J0-J1)*rpe
         RJ1p   = J1*rpe
         RJ0p   = J0*rpe
         sum1xx_a = sum1xx_a + RJ1p  * R
         sum1yy_a = sum1yy_a + RJ1p  * Rpr * y2
         sum1xz   = sum1xz   + RJ1p  * Rpr * y4
         sum1xx_c = sum1xx_c + RJ01p * Rpr * y2
         sum1yy_c = sum1yy_c + RJ01p * R
         sum1zz_c = sum1zz_c + RJ0p  * Rpr * y3
 31   continue
        
      sum2xx_a = cz
      sum2yy_a = cz
      sum2xz   = cz
      sum2xx_c = cz
      sum2yy_c = cz
      sum2zz_c = cz
      do 32, k = 2, N3-2, 2
         y      = div + k*dy3
         y2     = y*y
         y3     = ru+y2
         y4     = y*y3
         rpe    = dexp(-y*Zhat_p)
         kzeps  = cdsqrt(y2-eps1)
         R      = (y-kzeps)/(y+kzeps)
         Rpr    = (kzeps-eps*y)/(kzeps+eps*y)
         w      = Rhat*dsqrt(y3)
c$$$         J0 = DBESJ0(w)
c$$$         J1 = DBESJ1(w)/w
             CALL BJ01(W,J0,J1)
         RJ01p  = (J0-J1)*rpe
         RJ1p   = J1*rpe
         RJ0p   = J0*rpe
         sum2xx_a = sum2xx_a + RJ1p  * R
         sum2yy_a = sum2yy_a + RJ1p  * Rpr * y2
         sum2xz   = sum2xz   + RJ1p  * Rpr * y4
         sum2xx_c = sum2xx_c + RJ01p * Rpr * y2
         sum2yy_c = sum2yy_c + RJ01p * R  
         sum2zz_c = sum2zz_c + RJ0p  * Rpr * y3
 32   continue

      Gxx2a2 = dy3*(two*sum2xx_a + four*sum1xx_a + s1xx)/three
      Gyy2a2 = dy3*(two*sum2yy_a + four*sum1yy_a + s1yy)/three
      Gxz2a2 = dy3*(two*sum2xz   + four*sum1xz   + s1xz)/three
      Gxx2ah = Gxx2a1 + Gxx2a2
      Gyy2ah = Gyy2a1 + Gyy2a2
      Gxz2   = Gxz2a1 + Gxz2a2

      s2xx = rz
      s2yy = J1Hat - J0Hat
      s2zz = J0Hat

      y     = div
      y2    = y*y
      rpe   = dexp(-y*Zhat_p)
      kzeps = cdsqrt(y2-eps1)
      Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
      R     = (y-kzeps)/(y+kzeps)
      w     = Rhat*dsqrt(ru+y2)
c$$$     J0 = DBESJ0(w)
c$$$     J1 = DBESJ1(w)/w
         CALL BJ01(W,J0,J1)
      s1xx  = s2xx + y2*rpe*Rpr*(J0-J1)
      s1yy  = s2yy + rpe*R*(J0-J1)
      s1zz  = s2zz + (y2+ru)*rpe*Rpr*J0

      Gxx2b1 = dy2*(two*sum2xx_b + four*sum1xx_b + s1xx)/three
      Gyy2b1 = dy2*(two*sum2yy_b + four*sum1yy_b + s1yy)/three
      Gzz2b1 = dy2*(two*sum2zz_b + four*sum1zz_b + s1zz)/three

      y     = div
      y2    = y*y
      rpe   = dexp(-y*Zhat_p)
      kzeps = cdsqrt(y2-eps1)
      Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
      R     = (y-kzeps)/(y+kzeps)
      w     = Rhat*dsqrt(ru+y2)
c$$$     J0 = DBESJ0(w)
c$$$     J1 = DBESJ1(w)/w
         CALL BJ01(W,J0,J1)
      s2xx  = y2*rpe*Rpr*(J0-J1)
      s2yy  = rpe*R*(J0-J1)
      s2zz  = (y2+ru)*rpe*Rpr*J0

      y     = y_max
      y2    = y*y
      rpe   = dexp(-y*Zhat_p)
      kzeps = cdsqrt(y2-eps1)
      Rpr   = (kzeps-eps*y)/(kzeps+eps*y)
      R     = (y-kzeps)/(y+kzeps)
      w     = Rhat*dsqrt(ru+y2)
c$$$     J0 = DBESJ0(w)
c$$$     J1 = DBESJ1(w)/w
         CALL BJ01(W,J0,J1)
      s1xx  = s2xx + y2*rpe*Rpr*(J0-J1)
      s1yy  = s2yy + rpe*R*(J0-J1)
      s1zz  = s2zz + (y2+ru)*rpe*Rpr*J0

      Gxx2b2 = dy3*(two*sum2xx_c + four*sum1xx_c + s1xx)/three
      Gyy2b2 = dy3*(two*sum2yy_c + four*sum1yy_c + s1yy)/three
      Gzz2b2 = dy3*(two*sum2zz_c + four*sum1zz_c + s1zz)/three

      Gxx2bh = Gxx2b1 + Gxx2b2
      Gyy2bh = Gyy2b1 + Gyy2b2
      Gzz2   = Gzz2b1 + Gzz2b2

c     Dimensionless G^R/k**3
      GRxx = cu*(Gxx1a + Gxx1b) + Gxx2ah - Gxx2bh
      GRyy = cu*(Gyy1a + Gyy1b) + Gyy2bh - Gyy2ah
      GRzz = -cu*Gzz1 - Gzz2
      GRxz = -Rhat*(Gxz1 + Gxz2)
      GRzx = -GRxz

      return
      end


        SUBROUTINE BJ01(X,BJ0,BJ1)
C       ============================================================
C        Computes BJ0 = J0(x) and BJ1 = J1(x)/x
C        WARNING: do not use this subroutine in other applications.
C                 It will not return correct values for the case X=0.
C                 Also, the variable BJ1 contains J1(x)/x, not
C                 the Bessel function J1(x) itself
C       ============================================================
C
C       See S.Zhang, J.Jin, Computation of Special Functions, 
C                           Wiley, 1996, Sec.5.2. 
C       This subroutine is a modification of the subroutine JY01B
C       which is described in that section.
        IMPLICIT NONE

        REAL*8 PI,X,T,T2,BJ0,BJ1,A0,P0,Q0,P1,Q1,TA0,TA1
        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           BJ0=1.0D0
           BJ1=0.5D0 
           RETURN
        ELSE IF (X.LE.4.0D0) THEN
           T=X/4.0D0
           T2=T*T
           BJ0=((((((-.5014415D-3*T2+.76771853D-2)*T2
     &         -.0709253492D0)*T2+.4443584263D0)*T2
     &         -1.7777560599D0)*T2+3.9999973021D0)
     &         *T2-3.9999998721D0)*T2+1.0D0
           BJ1=T*(((((((-.1289769D-3*T2+.22069155D-2)
     &         *T2-.0236616773D0)*T2+.1777582922D0)*T2
     &         -.8888839649D0)*T2+2.6666660544D0)*T2
     &         -3.9999999710D0)*T2+1.9999999998D0)
        ELSE
           T=4.0D0/X
           T2=T*T
           A0=DSQRT(2.0D0/(PI*X))
           P0=((((-.9285D-5*T2+.43506D-4)*T2-.122226D-3)*T2
     &        +.434725D-3)*T2-.4394275D-2)*T2+.999999997D0
           Q0=T*(((((.8099D-5*T2-.35614D-4)*T2+.85844D-4)*T2
     &        -.218024D-3)*T2+.1144106D-2)*T2-.031249995D0)
           TA0=X-.25D0*PI
           BJ0=A0*(P0*DCOS(TA0)-Q0*DSIN(TA0))
           P1=((((.10632D-4*T2-.50363D-4)*T2+.145575D-3)*T2
     &        -.559487D-3)*T2+.7323931D-2)*T2+1.000000004D0
           Q1=T*(((((-.9173D-5*T2+.40658D-4)*T2-.99941D-4)*T2
     &        +.266891D-3)*T2-.1601836D-2)*T2+.093749994D0)
           TA1=X-.75D0*PI
           BJ1=A0*(P1*DCOS(TA1)-Q1*DSIN(TA1))
        ENDIF
        BJ1 = BJ1/X

        RETURN
        END
