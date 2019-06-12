    PROGRAM SPECTRUMS_INFTY
    USE IFLPORT
    implicit none      
    
    character*100 outputfile, mat
      
    integer*4 i, j, alpha, beta, nl, il
    integer*4 N, N_row_col, nc
      
    real*8 eps_ext
    real*8 li, l_min, l_max, dl
    real*8 theta_pol_e, phi_pol_e, theta_pol_m, phi_pol_m, theta_prop, phi_prop
    real*8 ru, rz, pi, twopi
    real*8 wv, wv_0, wv_1, wv_2, wv_3
    real*8 xij, yij, rij, rij_1, rij_2, rij_3, d_alpha_beta 
    real*8 rnijx, rnijy, rnij_alpha, rnij_beta, rnij_alpha_beta
    real*8 t1_g, t2_g, t3_g, arg
    real*8 t1_c, t2_c
    real*8 prop(3), pol_e(3), pol_m(3)
    real*8 time_begin, time_end
    real*8 rad
    real*8 sum_ext_e, Q_ext_e
    real*8 sum_ext_m, Q_ext_m
    real*8 sum_ext_em, Q_ext_em
    real*8 sigma_ext_em, sigma_ext_e, sigma_ext_m
    real*8 h_x, h_y
    
    complex*16 zet_e, zet_m, inf_p_x, inf_p_y, inf_m_x, inf_m_y
    complex*16 cu, cz, phase_fact
    complex*16 eps
    complex*16 G_tenzor_xx, G_tenzor_yy
    
    real*8,       allocatable, dimension(:)     :: x, y
    complex*16,   allocatable, dimension(:,:)   :: sim_G_tenzor
!c
!c------------------ Setting up parameters ------------------------ 
!c  
    parameter( cz = (0.0d0,0.0d0) )
    parameter( cu = (0.0d0,1.0d0) )
    parameter( rz = 0.0d0 )
    parameter( ru = 1.0d0 )
    parameter( pi = 4.0d0*datan(ru) )
    parameter( twopi = 2.0d0*pi )
!------------------------------------------ 
    parameter(h_x = 580)
    parameter(h_y = 480)
    parameter(rad = 65)
    parameter(mat = "Si")
    parameter(l_min = 420)
    parameter(l_max = 620)
    parameter(nl = 801)
    parameter(theta_pol_e = 90 * pi / 180.)
    parameter(phi_pol_e   = 0  * pi / 180.)
    parameter(theta_prop  = 10  * pi / 180.)
    parameter(phi_prop    = 90  * pi / 180.)
    parameter(outputfile = "580_480_y10.out")
    parameter(N_row_col = 10001)
    parameter(eps_ext = 1)
!c      
!c--------------------Checking for input errors-------------------------      
!c
    if (sin(theta_pol_e)*cos(phi_pol_e) * sin(theta_prop)*cos(phi_prop) + &
    sin(theta_pol_e)*sin(phi_pol_e) * sin(theta_prop)*sin(phi_prop) + &
    cos(theta_pol_e)*cos(theta_prop) .gt. 1.E-10) goto 205     
!c
!c------------------ Starting program------------------------ 
!c     
    CALL CPU_TIME (time_begin)    
    open(unit=70,file=outputfile, status='unknown')
!c 
!c------------------ Creating the mesh------------------------ 
!c    
    N = N_row_col**2
    allocate(x(1:N),y(1:N))
    do j = 1, N_row_col
        do i = 1, N_row_col
            x(i + N_row_col*(j-1)) = (i-1) * h_x 
            y(i + N_row_col*(j-1)) = (j-1) * h_y
        enddo
    enddo
    nc = (N-1)/2 + 1
!c      
!c-----------------------Allocating workspace-----------------------------      
!c      
    allocate (sim_G_tenzor(N,2))    
    sim_G_tenzor = cz 
    
    prop(1) = dsin(theta_prop)*dcos(phi_prop) 
    prop(2) = dsin(theta_prop)*dsin(phi_prop)
    prop(3) = dcos(theta_prop)
       
    pol_e(1) = dsin(theta_pol_e)*dcos(phi_pol_e)
    pol_e(2) = dsin(theta_pol_e)*dsin(phi_pol_e)
    pol_e(3) = dcos(theta_pol_e)
       
    pol_m(1) = prop(2)*pol_e(3) - prop(3)*pol_e(2)
    pol_m(2) = prop(3)*pol_e(1) - prop(1)*pol_e(3)
    pol_m(3) = prop(1)*pol_e(2) - prop(2)*pol_e(1)
!c   
!c----------------------------------------------------------------------- 
!c---------------------START LOOP OVER LAMBDA-----------------------------
    do 15 il = 1,nl 

        dl = (l_max - l_min) / dfloat(nl-1)
        li =  l_min + dl * dfloat(il-1)      
        
        print*,li
        
        wv = sqrt(eps_ext) * twopi / li
        wv_0 = ru 
        wv_1 = wv 
        wv_2 = wv*wv 
        wv_3 = wv_2*wv   
 
        call permittivity   (mat, li, eps_ext, eps)
        call susceptibility (eps, rad, wv, zet_e, zet_m)
!c
!c----------- Creating matrix A --------------------------------- 
!c           
        G_tenzor_xx = cz
        G_tenzor_yy = cz

            do j=1,N
               
                if (j .eq. nc) goto 20
                
                xij = x(j) - x(nc) 
                yij = y(j) - y(nc) 
                rij_2 = xij**2 + yij**2 
                rij = dsqrt(rij_2)
                rij_1 = rij 
                rij_3 = rij * rij_2 
                rnijx = xij/rij 
                rnijy = yij/rij 
                arg = wv * rij 
                phase_fact = cdexp(cu*arg)
          
                do 2, alpha=1,2
                    
                    beta = alpha
                    if (alpha .eq. 1) rnij_alpha = rnijx 
                    if (alpha .eq. 2) rnij_alpha = rnijy 
                    if (beta .eq. 1) rnij_beta = rnijx 
                    if (beta .eq. 2) rnij_beta = rnijy 
                    rnij_alpha_beta = rnij_alpha*rnij_beta 
                        
                    t1_g = wv_2 * (ru - 1.0*rnij_alpha_beta)/rij_1 
                    t2_g = wv_1 * (ru - 3.0*rnij_alpha_beta)/rij_2 
                    t3_g = wv_0 * (ru - 3.0*rnij_alpha_beta)/rij_3 

                    if (alpha .eq. 1 .and. beta .eq. 1) G_tenzor_xx = G_tenzor_xx + (t1_g + cu*t2_g - t3_g)*phase_fact
                    if (alpha .eq. 2 .and. beta .eq. 2) G_tenzor_yy = G_tenzor_yy + (t1_g + cu*t2_g - t3_g)*phase_fact
                        
                    if (alpha .eq. 1 .and. beta .eq. 1) sim_G_tenzor(j, 1) = G_tenzor_xx
                    if (alpha .eq. 2 .and. beta .eq. 2) sim_G_tenzor(j, 2) = G_tenzor_yy
                    
2               continue
20              continue
          
            end do

        !calculation average G tenzor
            
            G_tenzor_xx = cz
            G_tenzor_yy = cz
            
            j = 0.8*N
            
            do i = j, N
                G_tenzor_xx = G_tenzor_xx + sim_G_tenzor(i, 1)
                G_tenzor_yy = G_tenzor_yy + sim_G_tenzor(i, 2)
            end do
            
            G_tenzor_xx = dreal(G_tenzor_xx)/(N-j+1) + cu*dimag(G_tenzor_xx)/(N-j+1)
            G_tenzor_yy = dreal(G_tenzor_yy)/(N-j+1) + cu*dimag(G_tenzor_yy)/(N-j+1)
            
        !end of calculation average G tenzor
!c 
!c-------- End of costructing matrix A  -------------------------------- 
!c
!c 
!c--------- Infinite lattice --------------------------------------- 
!c
        inf_p_x = pol_e(1)/(zet_e - G_tenzor_xx)
        inf_p_y = pol_e(2)/(zet_e - G_tenzor_yy)
        inf_m_x = pol_m(1)/(zet_m - G_tenzor_xx)
        inf_m_y = pol_m(2)/(zet_m - G_tenzor_yy) 

        sum_ext_em = rz    
        sum_ext_e  = rz
        sum_ext_m  = rz

        sum_ext_e  = aimag(inf_p_x * pol_e(1) + inf_p_y * pol_e(2))
        sum_ext_m  = aimag(inf_m_x * pol_m(1) + inf_m_y * pol_m(2))
        sum_ext_em = sum_ext_e + sum_ext_m
            
        !sigma_ext_em = 4 * wv_1 * sum_ext_em / (rad_mean**2.)
        !sigma_ext_e = 4 * wv_1 * sum_ext_e / (rad_mean**2.)
        !sigma_ext_m = 4 * wv_1 * sum_ext_m / (rad_mean**2.)
 
        sigma_ext_em = 4 *pi* wv_1 * sum_ext_em * 1e-6
        sigma_ext_e  = 4 *pi* wv_1 * sum_ext_e  * 1e-6
        sigma_ext_m  = 4 *pi* wv_1 * sum_ext_m  * 1e-6
        
        !Q_ext_em = 4.*wv_1*pi*aimag(inf_p_x * pol_e(1) + inf_p_y * pol_e(2))/rad_mean**2.
        !Q_ext_e = 4.*wv_1*pi*aimag(inf_m_x * pol_m(1) + inf_m_y * pol_m(2))/rad_mean**2.
        !Q_ext_m = 4.*wv_1*pi*aimag(inf_p_x * pol_e(1) + inf_p_y * pol_e(2) + inf_m_x * pol_m(1) + inf_m_y * pol_m(2))/rad_mean**2.
        
        !write(70,'(6F7.3)') li, sngl(sigma_ext), sngl(Q_ext_em), sngl(Q_ext_e), sngl(Q_ext_m)
        write(70,'(15F10.3)') li, sigma_ext_em, real(1.0e6*zet_e), real(1.0e6*zet_m), real(1.0e6*G_tenzor_xx), real(1.0e6*G_tenzor_yy) !, real(1.0e-6/zet_et(1)),imag(1.0e-6/zet_et(1)), real(1.0e-6/zet_mt(1)),imag(1.0e-6/zet_mt(1))
!c---------------------------------------------------------------------                  
!c---------------------------------------------------------------------      
15  end do 
!c---------------------END LOOP OVER OMEGA-----------------------------  
!c--------------------------------------------------------------------- 
!c       
    deallocate(sim_G_tenzor)
    
    close(70)
    close(80)
 
    CALL CPU_TIME (time_end)
    write(*,*) 'estimated_time=', time_end-time_begin, 'sec'
!c-----------------------------------------------------------------------      
 500  stop 
!c-------------------------ERRORS----------------------------------------  
  
 205  continue
      write(*,*) 'E and k are not orthogonal; exiting'
      stop          
     
    end