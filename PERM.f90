SUBROUTINE PERMITTIVITY_ED (mat,lambda,eps_ext,eps) 

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
      
      real*8 eps_inf, G_D
      real*8 omega_L1, omega_01
      real*8 omega_L2, gamma_2, omega_02 
      
      complex*16 cu
      
      real*8 Am1, Am2, Am3, Br1, Br2, Br3, En1, En2, En3
               
      cu = (0.0d0,1.0d0) 
      ru = 1.0d0 
      pi = 4.0d0*datan(ru) 
      twopi = 2*pi 

!c------------------------SILVER--------------------------------------
      if (mat .eq. 1) then 
           
      !lambda_p = 136.1
      !g_inf = .0019
      !eps_0 = 5.
      !eps = eps_0 * ru - ( lambda / lambda_p )**2 / ( ru + cu *  g_inf * lambda / lambda_p)
      
    ! lambda_p = 109.5
    ! eps_0 = 5.45
    ! g_inf = .00485
    ! eps = eps_0 * ru - 0.73 * ( lambda / lambda_p )**2 / ( ru + cu *  g_inf * lambda / lambda_p )
     open(unit=81, file='johnson.silver.dat', status='old')
        
        l_min = 10000.
        l_max = 0.
        read(81,*,end=14) l, r1, r2
        if(l .lt. l_min) l_min = l
        if(l .gt. l_max) l_max = l
        i = 1
        l_old = l
        r1_old = r1
        r2_old = r2

15	continue
	read(81,*,end=14) l, r1, r2
        if(l .lt. l_min) l_min = l
        if(l .gt. l_max) l_max = l
	i = i + 1
        if(l .eq. lambda) goto 14
        if( (lambda-l_old)*(lambda-l) .le. 0.) then
           r1 = r1_old + (lambda-l_old)*(r1-r1_old)/(l-l_old)
           r2 = r2_old + (lambda-l_old)*(r2-r2_old)/(l-l_old)
           goto 14
        end if
        l_old = l
        r1_old = r1
        r2_old = r2
        goto 15
14      continue
        close (unit=81)
        
        eps = (r1 + cu*r2)**2
                   
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
!c-----------------------------Si------------------------------------
      if (mat .eq. 8) then 

            open(unit=81, file='palik.silicon.dat', status='old')
        
        l_min = 10000.
        l_max = 0.
        read(81,*,end=17) l, r1, r2
        if(l .lt. l_min) l_min = l
        if(l .gt. l_max) l_max = l
        i = 1
        l_old = l
        r1_old = r1
        r2_old = r2

16	continue
	read(81,*,end=17) l, r1, r2
        if(l .lt. l_min) l_min = l
        if(l .gt. l_max) l_max = l
	i = i + 1
        if(l .eq. lambda) goto 17
        if( (lambda-l_old)*(lambda-l) .le. 0.) then
           r1 = r1_old + (lambda-l_old)*(r1-r1_old)/(l-l_old)
           r2 = r2_old + (lambda-l_old)*(r2-r2_old)/(l-l_old)
           goto 17
        end if
        l_old = l
        r1_old = r1
        r2_old = r2
        goto 16
17      continue
        close (unit=81)
        
        eps = (r1 + cu*r2)**2
      
      end if
      
      !c------------------------T=023--------------------------------------
	  if (mat .eq. 23) then
		eps_inf = 4.74
		omega_p = 7.01
		G_D		= 0.23
		omega_L1= 36.7
		gamma_1 = 1.45
		omega_01= 4.10
		omega_L2= 3.88
		gamma_2 = 1.33
		omega_02= 1.85
        omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
        
        eps = eps_inf - omega_p**2/(omega**2 + cu*G_D*omega) + omega_L1/(omega_01**2-omega**2-cu*gamma_1*omega) + omega_L2/(omega_02**2-omega**2-cu*gamma_2*omega)
	  end if
!c------------------------T=200--------------------------------------
	  if (mat .eq. 200) then
		eps_inf = 4.35
		omega_p = 7.06
		G_D		= 0.28
		omega_L1= 39.0
		gamma_1 = 1.15
		omega_01= 4.01
		omega_L2= 4.33
		gamma_2 = 1.38
		omega_02= 1.92
        omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
        
        eps = eps_inf - omega_p**2/(omega**2 + cu*G_D*omega) + omega_L1/(omega_01**2-omega**2-cu*gamma_1*omega) + omega_L2/(omega_02**2-omega**2-cu*gamma_2*omega)
	  end if	  
!c------------------------T=400--------------------------------------
	  if (mat .eq. 400) then
		eps_inf = 3.89
		omega_p = 7.20
		G_D		= 0.37
		omega_L1= 54.4
		gamma_1 = 1.40
		omega_01= 4.24
		omega_L2= 3.61
		gamma_2 = 1.46
		omega_02= 1.95
        omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
        
        eps = eps_inf - omega_p**2/(omega**2 + cu*G_D*omega) + omega_L1/(omega_01**2-omega**2-cu*gamma_1*omega) + omega_L2/(omega_02**2-omega**2-cu*gamma_2*omega)
	  end if      
!c------------------------T=500--------------------------------------
	  if (mat .eq. 500) then
		eps_inf = 3.58
		omega_p = 7.25
		G_D		= 0.40
		omega_L1= 66.3
		gamma_1 = 1.54
		omega_01= 4.39
		omega_L2= 3.58
		gamma_2 = 1.54
		omega_02= 1.94
        omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
        
        eps = eps_inf - omega_p**2/(omega**2 + cu*G_D*omega) + omega_L1/(omega_01**2-omega**2-cu*gamma_1*omega) + omega_L2/(omega_02**2-omega**2-cu*gamma_2*omega)
	  end if	  
!c------------------------T=600--------------------------------------
	  if (mat .eq. 600) then
		eps_inf = 3.37
		omega_p = 7.31
		G_D		= 0.45
		omega_L1= 76.8
		gamma_1 = 1.70
		omega_01= 4.50
		omega_L2= 2.91
		gamma_2 = 1.44
		omega_02= 1.96
        omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
        
        eps = eps_inf - omega_p**2/(omega**2 + cu*G_D*omega) + omega_L1/(omega_01**2-omega**2-cu*gamma_1*omega) + omega_L2/(omega_02**2-omega**2-cu*gamma_2*omega)
	  end if
!c------------------------T=700--------------------------------------
	  if (mat .eq. 700) then
		eps_inf = 3.14
		omega_p = 7.37
		G_D		= 0.50
		omega_L1= 92.2
		gamma_1 = 1.94
		omega_01= 4.67
		omega_L2= 2.22
		gamma_2 = 1.27
		omega_02= 1.94
        omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
        
        eps = eps_inf - omega_p**2/(omega**2 + cu*G_D*omega) + omega_L1/(omega_01**2-omega**2-cu*gamma_1*omega) + omega_L2/(omega_02**2-omega**2-cu*gamma_2*omega)
	  end if
!c------------------------T=800--------------------------------------
	  if (mat .eq. 800) then
		eps_inf = 2.67
		omega_p = 7.39
		G_D		= 0.58
		omega_L1= 114.9
		gamma_1 = 2.13
		omega_01= 4.90
		omega_L2= 1.15
		gamma_2 = 0.92
		omega_02= 2.00
        omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
        
        eps = eps_inf - omega_p**2/(omega**2 + cu*G_D*omega) + omega_L1/(omega_01**2-omega**2-cu*gamma_1*omega) + omega_L2/(omega_02**2-omega**2-cu*gamma_2*omega)
	  end if
!c------------------------T=900--------------------------------------
	  if (mat .eq. 900) then
		eps_inf = 2.08
		omega_p = 7.17
		G_D		= 0.66
		omega_L1= 122.6
		gamma_1 = 2.13
		omega_01= 4.96
		omega_L2= 0.46
		gamma_2 = 0.58
		omega_02= 2.01
        omega = 4.135667*1.0d-15*2.99792458*1.0d8/(lambda*1.0d-9)
        
        eps = eps_inf - omega_p**2/(omega**2 + cu*G_D*omega) + omega_L1/(omega_01**2-omega**2-cu*gamma_1*omega) + omega_L2/(omega_02**2-omega**2-cu*gamma_2*omega)
	  end if
!c------------------------FINAL NORMALIZATION--------------------------------------
      eps = eps / eps_ext
      
      RETURN

END SUBROUTINE