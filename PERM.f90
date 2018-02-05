SUBROUTINE PERMITTIVITY (mat,lambda,eps_ext,eps) 

      implicit none
      integer*4, intent(in) :: mat
      real*8,    intent(in) :: eps_ext, lambda
      complex*16, intent(out) :: eps
      integer*8 i
      real*8 lambda_p, eps_0, g_inf 
      real*8 ru, pi, twopi
      real*8 omega_1, gamma_1, omega_p, f_1, omega
      real*8 r1, r2, r1_old, r2_old
      real*8 l_min, l_max, l, l_old
      complex*16 cu
      
      real*8 Am1, Am2, Am3, Br1, Br2, Br3, En1, En2, En3
               
      cu = (0.0d0,1.0d0) 
      ru = 1.0d0 
      pi = 4.0d0*datan(ru) 
      twopi = 2*pi 

!c------------------------SILVER--------------------------------------
      if (mat .eq. 1) then 
           
      lambda_p = 136.1
      g_inf = .0019
      eps_0 = 5.
      eps = eps_0 * ru - ( lambda / lambda_p )**2 / ( ru + cu *  g_inf * lambda / lambda_p)
      
 ! lambda_p = 109.5
    ! eps_0 = 5.45
    ! g_inf = .00485
    ! eps = eps_0 * ru - 
    !$ 0.73 * ( lambda / lambda_p )**2 / ( ru +  
    !$                              cu *  g_inf * lambda / lambda_p )
     
      end if
!c------------------------GOLD----------------------------------------
      if (mat .eq. 2) then 
      
      open(unit=81, file='johnson.gold.dat', status='old')

        l_min = 10000.
        l_max = 0.
        read(81,*,end=13) l, r1, r2
        if(l .lt. l_min) l_min = l
        if(l .gt. l_max) l_max = l
        i = 1
        l_old = l
        r1_old = r1
        r2_old = r2

12	continue
	read(81,*,end=13) l, r1, r2
        if(l .lt. l_min) l_min = l
        if(l .gt. l_max) l_max = l
	i = i + 1
        if(l .eq. lambda) goto 13
        if( (lambda-l_old)*(lambda-l) .le. 0.) then
           r1 = r1_old + (lambda-l_old)*(r1-r1_old)/(l-l_old)
           r2 = r2_old + (lambda-l_old)*(r2-r2_old)/(l-l_old)
           goto 13
        end if
        l_old = l
        r1_old = r1
        r2_old = r2
        goto 12
13      continue
        close (unit=81)
        
        eps = (r1 + cu*r2)**2

      end if
!c------------------------TiN 800 C------------------------------------
      if (mat .eq. 3) then 
           
!      g_inf = .1795
!      eps_0 = 4.855
!      f_1 = 3.2907
!      omega_p = 7.9308
!      omega_1 = 4.2196
!      gamma_1 = 2.0341
!      omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
!
!      eps = eps_0*ru - omega_p**2/(omega**2+cu*omega*g_inf) +
!     $ f_1*omega_1**2 / (omega_1**2-omega**2-cu*omega*gamma_1)

      eps_0 = 3.6684
      Am1 = 47.324
      Am2 = 83.573
      Am3 = 3.193
      En1 = 0.0d0
      En2 = 4.7808
      En3 = 1.9289
      Br1 = 0.28058
      Br2 = 1.9188
      Br3 = 0.98543
      
      omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
      
      eps = eps_0 + Am1/(En1**2 - omega**2 - cu*Br1*omega) + Am2/(En2**2 - omega**2 - cu*Br2*omega) + Am3/(En3**2 - omega**2 - cu*Br3*omega)
      
      end if
!c-----------------------------ZrN---------------------------------------
      if (mat .eq. 4) then 

      g_inf = 0.5192
      eps_0 = 3.4656
      f_1 = 2.4509
      omega_p = 8.018
      omega_1 = 5.48
      gamma_1 = 1.7369 
      omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
                  
      eps = eps_0*ru - omega_p**2/(omega**2+cu*omega*g_inf) + f_1*omega_1**2 / (omega_1**2-omega**2-cu*omega*gamma_1)
     
      end if
!c-----------------------------AZO---------------------------------------
      if (mat .eq. 5) then 

      g_inf = 0.04486
      eps_0 = 3.5402
      
      f_1 = 0.5095
      omega_p = 1.7473
      omega_1 = 4.2942
      gamma_1 = 0.1017 
      omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
                  
      eps = eps_0*ru - omega_p**2/(omega**2+cu*omega*g_inf) + f_1*omega_1**2 / (omega_1**2-omega**2-cu*omega*gamma_1)
     
      end if
!c-----------------------------GZO---------------------------------------
      if (mat .eq. 6) then 

      g_inf = 0.1229
      eps_0 = 3.2257
      
      f_1 = 0.3859
      omega_p = 1.9895
      omega_1 = 4.050
      gamma_1 = 0.0924 
      omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
                  
      eps = eps_0*ru - omega_p**2/(omega**2+cu*omega*g_inf) + f_1*omega_1**2 / (omega_1**2-omega**2-cu*omega*gamma_1)
     
      end if
!c-----------------------------ITO---------------------------------------
      if (mat .eq. 7) then 

      g_inf = 0.155
      eps_0 = 3.528
      
      f_1 = 0.3884
      omega_p = 1.78
      omega_1 = 4.210
      gamma_1 = 0.0919 
      omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
                  
      eps = eps_0*ru - omega_p**2/(omega**2+cu*omega*g_inf) + f_1*omega_1**2 / (omega_1**2-omega**2-cu*omega*gamma_1)
     
      end if
!c-----------------------------CUSTOM------------------------------------
      if (mat .eq. 8) then 
      
!      eps = ru - (lambda/314)**2 / ( ru +  
!     $                           cu *  0.014 * lambda/ 314)
!     
      eps = (2.25,0.01)
      
      end if
      
      eps = eps / eps_ext
      
      RETURN

END SUBROUTINE