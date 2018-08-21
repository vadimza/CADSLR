    PROGRAM SPECTRUMS_ED_MD
    USE IFLPORT
    implicit none      
    
    character*100 inputfile, outputfile, coordfile, excitation, nar, mati
    character*100 tmpstr
      
    logical exist_flag_work, result
      
    integer*4 N, Na, i, j, alpha, beta, nl, ja
    integer*4 il, ia
    integer*4 lwork, liwork, info, istat, cnt
      
    real*8 aspect, b1, b2, dist, eps_ext, fill, rad_mean
    real*8 li, l_min, l_max, dl
    real*8 theta_pol_e, phi_pol_e, theta_pol_m, phi_pol_m, theta_prop, phi_prop
    real*8 ru, rz, pi, twopi, mu0, eps0, mueps, mueps_sq, epsmu_sq
    real*8 wv, wv_0, wv_1, wv_2, wv_3
    real*8 xij, yij, zij, rij, rij_1, rij_2, rij_3, d_alpha_beta 
    real*8 rnijx, rnijy, rnijz, rnij_alpha, rnij_beta, rnij_alpha_beta
    real*8 t1_g, t2_g, t3_g, arg
    real*8 t1_c, t2_c
    real*8 prop(3), pol_e(3), pol_m(3)
    real*8 time_begin, time_end
    real*8 arg_inc
    real*8 sum_ext_e, Q_ext_e
    real*8 sum_ext_m, Q_ext_m
    real*8 sum_ext_em, Q_ext_em
    real*8 sum_scat_e, Q_scat_e
    real*8 sum_scat_m, Q_scat_m
    real*8 sum_scat_em, Q_scat_em    
    
      
    complex*16 cu,ncu, cz, phase_fact
    complex*16 eps
      
    integer*4,    allocatable, dimension(:)     :: iwork
    real*8,       allocatable, dimension(:)     :: x, y, z, rad
    complex*16,   allocatable, dimension(:)     :: work, zet_e, zet_m
    complex*16,                dimension(3)     :: zet_et, zet_mt
    
    complex*16,   allocatable, dimension(:)     :: aed, amm
    complex*16,   allocatable, dimension(:,:)   :: a, rhs      
    complex*16,   allocatable, dimension(:,:)   :: E_field, H_field
    complex*16,   allocatable, dimension(:,:)   :: a_e, a_m, a_em 
    complex*16,   allocatable, dimension(:,:)   :: rhs_ED_MD
    complex*16,   allocatable, dimension(:,:)   :: G_tenzor, C_tenzor
    character*10, allocatable, dimension(:)     :: mat
!c
!c------------------ Reading parameters from input file ------------------------ 
!c  
    CALL CPU_TIME (time_begin)
    
    cnt = command_argument_count ()
    if(cnt .GT. 0) then
        CALL GET_COMMAND_ARGUMENT (1, nar)
        read (nar,*) inputfile
        inputfile = trim(adjustl(inputfile))
    endif
      
    open(unit=70,file=inputfile,status='old')

    do
        
        read(70, *, end=10) tmpstr
        tmpstr = trim(adjustl(tmpstr))

        if (tmpstr .eq. "end_of_options" )     goto 10

        if (tmpstr .eq. "number_of_particles") read(70,*)   N
        if (tmpstr .eq. "lambda_1")            read(70,*)   l_min
        if (tmpstr .eq. "lambda_2")            read(70,*)   l_max
        if (tmpstr .eq. "lambda_n")            read(70,*)   nl
        if (tmpstr .eq. "pol_E_theta_phi")     read(70,*)   theta_pol_e, phi_pol_e
        if (tmpstr .eq. "prop_k_theta_phi")    read(70,*)   theta_prop, phi_prop
        if (tmpstr .eq. "excitation_type")     read(70,*)   excitation
        if (tmpstr .eq. "host_medium")         read(70,*)   eps_ext
        if (tmpstr .eq. "structure_file")      read(70,*)   coordfile
        if (tmpstr .eq. "output_file")         read(70,*)   outputfile
        
    end do
      
10  close(70)
    
    coordfile  = trim(adjustl(coordfile))
    outputfile = trim(adjustl(outputfile))
    excitation = trim(adjustl(excitation))
!c 
!c------------------ Reading array from coordfile file ------------------------ 
!c 
    allocate(x(1:N),y(1:N),z(1:N),rad(1:N),mat(1:N))
    
    open(unit=70,file=coordfile,status='old')      
    
    do i = 1, N
        
        read(70,*) x(i), y(i), z(i), rad(i), mat(i)
        
    enddo
!c
!c-----------Calculating rad_mean-----------------------
!c
    rad_mean = 0
    
    do i = 1, N
        
        rad_mean = rad_mean + rad(i)
        
    enddo
    
    rad_mean = rad_mean/N
!c 
!c--------- Constants -------------------------- 
!c 
    cz = (0.0d0,0.0d0) 
    cu = (0.0d0,1.0d0) 
    ncu= (0.0d0,-1.0d0) 
    rz = 0.0d0 
    ru = 1.0d0 
    pi = 4.0d0*datan(ru) 
    twopi = 2.0d0*pi
    mu0   = 4.0d0*pi * 1.0d-7
    eps0  = 8.85d0   * 1.0d-12
    mueps = mu0 / eps0
    mueps_sq = sqrt(   mueps)
    epsmu_sq = sqrt(1./mueps)
!c
!c--------------------Converting degrees to radians---------------------
!c
    theta_pol_e = theta_pol_e * pi / 180.
    phi_pol_e   = phi_pol_e   * pi / 180.
    theta_pol_m = theta_pol_m * pi / 180.
    phi_pol_m   = phi_pol_m   * pi / 180.
    theta_prop  = theta_prop  * pi / 180.
    phi_prop    = phi_prop    * pi / 180.
!c      
!c--------------------Checking for input errors-------------------------      
!c
    if (sin(theta_pol_e)*cos(phi_pol_e) * sin(theta_prop)*cos(phi_prop) + sin(theta_pol_e)*sin(phi_pol_e) * sin(theta_prop)*sin(phi_prop) + cos(theta_pol_e)*cos(theta_prop) .gt. 1.E-10) goto 205
!c      
!c-----------------------Allocating workspace-----------------------------      
!c      
    Na = 3*N
  
    allocate(E_field(1:Na,1:1), H_field(1:Na,1:1))
    allocate(a_em(1:2*Na,1:2*Na))
       
    allocate(zet_e(1:Na))
    allocate(zet_m(1:Na))
    
      
    allocate(rhs_ED_MD(1:2*Na,1:1))
    allocate(a_e(1:Na,1:Na), a_m(1:Na,1:Na))
    allocate(G_tenzor(1:Na,1:Na), C_tenzor(1:Na,1:Na))
       
    prop(1) = dsin(theta_prop)*dcos(phi_prop) 
    prop(2) = dsin(theta_prop)*dsin(phi_prop)
    prop(3) = dcos(theta_prop)
       
    pol_e(1)  = dsin(theta_pol_e)*dcos(phi_pol_e)
    pol_e(2)  = dsin(theta_pol_e)*dsin(phi_pol_e)
    pol_e(3)  = dcos(theta_pol_e)
       
    pol_m(1)  = prop(2)*pol_e(3) - prop(3)*pol_e(2)
    pol_m(2)  = prop(3)*pol_e(1) - prop(1)*pol_e(3)
    pol_m(3)  = prop(1)*pol_e(2) - prop(2)*pol_e(1)
!c   
!c----------------------------------------------------------------------- 
!c---------------------START LOOP OVER LAMBDA-----------------------------
    open(unit=70,file=outputfile, status='unknown')
    
    do 15 il = 1,nl
        
        dl = (l_max - l_min) / dfloat(nl-1)
        li =  l_min + dl * dfloat(il-1)      
        print*,li

        wv = sqrt(eps_ext) * twopi / li
        wv_0 = ru 
        wv_1 = wv 
        wv_2 = wv*wv 
        wv_3 = wv_2*wv
        
        do i = 1,N
            
            call permittivity  (mat(i), li, eps_ext, eps)
            call susceptibility (eps, rad(i), wv, zet_et, zet_mt)
            zet_e(3*i-2) = zet_et(1)
            zet_e(3*i-1) = zet_et(2)
            zet_e(3*i  ) = zet_et(3)
            zet_m(3*i-2) = zet_mt(1)
            zet_m(3*i-1) = zet_mt(2)
            zet_m(3*i  ) = zet_mt(3)
            
        enddo
!c
!c----------- Creating matrix A --------------------------------- 
!c           
        a_e = cz 
        G_tenzor = cz
        C_tenzor = cz
      
        do 1, i=1,N-1
            
            do 1, j=i+1,N
      
                xij = x(j) - x(i) 
                yij = y(j) - y(i) 
                zij = z(j) - z(i) 
                rij_2 = xij**2 + yij**2 + zij**2 
                rij = dsqrt(rij_2) !inter-center distance
                rij_1 = rij 
                rij_3 = rij * rij_2 
                rnijx = xij/rij 
                rnijy = yij/rij 
                rnijz = zij/rij 
                arg = wv * rij 
                phase_fact = cdexp(cu*arg)
          
                do 2, alpha=1,3           !row
                    
                    ia = 3*(i-1) + alpha 
                    if (alpha .eq. 1) rnij_alpha = rnijx 
                    if (alpha .eq. 2) rnij_alpha = rnijy 
                    if (alpha .eq. 3) rnij_alpha = rnijz 
            
                    do 3, beta=1,3            !column
                    
                        ja = 3*(j-1) + beta 
                        if (beta .eq. 1) rnij_beta = rnijx 
                        if (beta .eq. 2) rnij_beta = rnijy 
                        if (beta .eq. 3) rnij_beta = rnijz 
                        rnij_alpha_beta = rnij_alpha*rnij_beta 
                        
                        if(alpha .eq. beta) then !Unit tenzor 
                            d_alpha_beta = ru 
                        else 
                            d_alpha_beta = rz 
                        end if 
                        
                        t1_g = wv_2 * (d_alpha_beta - 1.0*rnij_alpha_beta)/rij_1 
                        t2_g = wv_1 * (d_alpha_beta - 3.0*rnij_alpha_beta)/rij_2 
                        t3_g = wv_0 * (d_alpha_beta - 3.0*rnij_alpha_beta)/rij_3 

                        G_tenzor(ja,ia) = (t1_g + cu*t2_g - t3_g)*phase_fact
              
3                   continue
                    
2               continue
          
                do 21, alpha=1,3           !row
                
                    ia = 3*(i-1) + alpha 
                    
                    do 31, beta=1,3            !column
                
                        ja = 3*(j-1) + beta 
                        rnij_alpha_beta = rz
              
                        if (alpha .eq. 1 .and. beta .eq. 2) rnij_alpha_beta =        rnijz
                        if (alpha .eq. 1 .and. beta .eq. 3) rnij_alpha_beta = -1.0 * rnijy
                        if (alpha .eq. 2 .and. beta .eq. 1) rnij_alpha_beta = -1.0 * rnijz
                        if (alpha .eq. 2 .and. beta .eq. 3) rnij_alpha_beta =        rnijx
                        if (alpha .eq. 3 .and. beta .eq. 1) rnij_alpha_beta =        rnijy
                        if (alpha .eq. 3 .and. beta .eq. 2) rnij_alpha_beta = -1.0 * rnijx
              
                        t1_c = wv_2 * rnij_alpha_beta / rij_1
                        t2_c = wv_1 * rnij_alpha_beta / rij_2
                        C_tenzor(ja,ia) = (t1_c + cu*t2_c) * phase_fact 

31                  continue
                    
21              continue
 
1       continue 
          
        do 4, i=2,Na
            
            do 4, j=1,i-1
                
                G_tenzor(j, i) = G_tenzor(i, j)
                C_tenzor(j, i) = C_tenzor(i, j)
                
4       continue
 
        do 5, i=1,Na
            
            do 5, j=1,Na
                
                if (i .eq. j) then        !for diagonal elements a_e matrix
                    
                    a_e(j,i) = zet_e(i)! - G_tenzor(j,i)
                    a_m(j,i) = zet_m(i)! - C_tenzor(j,i)
                  
                else                      !for off diagonal elements G_ED matrix
                    
                    G_tenzor(j,i) = -1.0 * G_tenzor(j,i)
                    
                end if
                
5       continue
!c 
!c-------- Costructing ED MD matrix A  -------------------------------- 
!c
        a_em = cz
    
        do i = 1, 3 * N 
        
            do j = 1, 3 * N
            
                if(i .eq. j) a_em(    i,     j) =  a_e(i, j)                        !left up part           diagonal elements
                if(i .eq. j) a_em(3*N+i, 3*N+j) =  a_m(i, j)                        !right down part        diagonal elements
            
                if(i .ne. j) a_em(    i,     j) =  G_tenzor(i, j)                    !left up part       off diagonal elements
                if(i .ne. j) a_em(3*N+i, 3*N+j) =  G_tenzor(i, j)                    !right down part    off diagonal elements
                if(i .ne. j) a_em(    i, 3*N+j) = -1.0 * epsmu_sq * C_tenzor(i, j)   !left down part     off diagonal elements 
                if(i .ne. j) a_em(3*N+i,     j) =  mueps_sq* C_tenzor(i, j)   !right up part      off diagonal elements
        
            end do
             
        end do
!c 
!c-------- End of costructing ED MD matrix A  -------------------------------- 
!c
    
!c 
!c-------- Costructing right-hand side  -------------------------------- 
!c  
        E_field = cz
        H_field = cz
  
        if (excitation .eq. 'plane') then           !Plane wave
        
            do i=1,N
                
                arg_inc = prop(1) * x(i) + prop(2) * y(i) + prop(3) * z(i)
                E_field(3*i-2,1) =            pol_e(1) * cdexp(cu * wv_1 * arg_inc)
                E_field(3*i-1,1) =            pol_e(2) * cdexp(cu * wv_1 * arg_inc)
                E_field(3*i  ,1) =            pol_e(3) * cdexp(cu * wv_1 * arg_inc)
                H_field(3*i-2,1) = epsmu_sq * pol_m(1) * cdexp(cu * wv_1 * arg_inc)
                H_field(3*i-1,1) = epsmu_sq * pol_m(2) * cdexp(cu * wv_1 * arg_inc)
                H_field(3*i  ,1) = epsmu_sq * pol_m(3) * cdexp(cu * wv_1 * arg_inc)

            end do
   
        endif
        
        if (excitation .eq. 'tip') then           !Tip
        
            E_field(1,1) =            pol_e(1)  
            E_field(2,1) =            pol_e(2) 
            E_field(3,1) =            pol_e(3)
            H_field(1,1) = epsmu_sq * pol_m(1)  
            H_field(2,1) = epsmu_sq * pol_m(2) 
            H_field(3,1) = epsmu_sq * pol_m(3)
        
        endif    

        do i = 1, Na
            
            rhs_ED_MD(   i,1) = E_field(i,1)
            rhs_ED_MD(Na+i,1) = H_field(i,1)
        
        end do
!c 
!c--------- Solving equations --------------------------------------- 
!c  
        liwork = 2*Na 
        allocate(iwork(liwork)) 
        
        allocate(work(1)) 
        call zsysv('U',2*Na,1,a_em,2*Na,iwork,rhs_ED_MD,2*Na,work,-1,info) 
        lwork = work(1) 
        deallocate(work) 
        allocate(work(lwork))
        
        call zsysv('U',2*Na,1,a_em,2*Na,iwork,rhs_ED_MD,2*Na,work,lwork,info) 
        
        if(info .ne. 0) then 
            write(*,*) 'Matrix inversion has failed; exiting' 
            write(*,*) 'info=', info 
            stop 
        end if
        
        deallocate(work) 
        deallocate(iwork)        
!c 
!c---------- Calculating spectrums------------------------------- 
!c  
!________________________ext_________________________        
        sum_ext_em = rz    
        sum_ext_e = rz
        sum_ext_m = rz

        do i=1,Na

            sum_ext_em = sum_ext_em + aimag(rhs_ED_MD(i,1)*dconjg(E_field(i,1))) + mueps*aimag(rhs_ED_MD(Na+i,1)*dconjg(H_field(i,1)))
            sum_ext_e = sum_ext_e + aimag(rhs_ED_MD(i,1) * dconjg(E_field(i,1))) 
            sum_ext_m = sum_ext_m + mueps*aimag(rhs_ED_MD(Na+i,1)*dconjg(H_field(i,1)))
            
        end do
          
        Q_ext_em = 4.*wv_1*sum_ext_em/(float(N)*rad_mean**2.)
        Q_ext_e = 4.*wv_1*sum_ext_e/(float(N)*rad_mean**2.)
        Q_ext_m = 4.*wv_1*sum_ext_m/(float(N)*rad_mean**2)
!_________________________________________________        
        
        
        
!________________________scat_________________________        
        sum_scat_em = rz    
        sum_scat_e = rz
        sum_scat_m = rz        
        
        
        allocate(aed(1:Na))
        allocate(amm(1:Na))
        
        do i = 1, N
            aed(3*i-2) = (2.*cu*wv_3/3. + zet_e(3*i-2)) * rhs_ED_MD(3*i-2,1)
            aed(3*i-1) = (2.*cu*wv_3/3. + zet_e(3*i-1)) * rhs_ED_MD(3*i-1,1)
            aed(3*i  ) = (2.*cu*wv_3/3. + zet_e(3*i  )) * rhs_ED_MD(3*i,1)
            
            amm(3*i-2) = (2.*cu*wv_3/3. + zet_m(3*i-2)) * rhs_ED_MD(3*N+3*i-2,1)
            amm(3*i-1) = (2.*cu*wv_3/3. + zet_m(3*i-1)) * rhs_ED_MD(3*N+3*i-1,1)
            amm(3*i  ) = (2.*cu*wv_3/3. + zet_m(3*i  )) * rhs_ED_MD(3*N+3*i,1)
        end do
        
        do i=1,Na
            
            sum_scat_e = sum_scat_e + aimag(dconjg(rhs_ED_MD(i,1))    * aed(i) - dconjg(rhs_ED_MD(i,1))    * E_field(i,1))  
            sum_scat_m = sum_scat_e + aimag(dconjg(rhs_ED_MD(Na+i,1)) * amm(i) - dconjg(rhs_ED_MD(Na+i,1)) * H_field(i,1))  
            
        end do
        
        Q_scat_e =         sum_scat_e * (4.*wv_1/(float(N)*rad_mean**2.))       !sigma_s^e + P      from MERCHIERS 27a
        Q_scat_m = mueps * sum_scat_e * (4.*wv_1/(float(N)*rad_mean**2.))       !sigma_s^m + Q      from MERCHIERS 27b
        
        Q_scat_em = Q_scat_e + Q_scat_m !sigma_s^e + P + sigma_s^m + Q
!_________________________________________________        

        
        
        write(70,*) sngl(Q_ext_em),sngl(Q_scat_em) !, sngl(Q_ext_e), sngl(Q_ext_m)
!c---------------------------------------------------------------------                  
!c---------------------------------------------------------------------      
15  end do 
!c---------------------END LOOP OVER OMEGA-----------------------------  
!c--------------------------------------------------------------------- 
!c       
    deallocate(a_em)
    deallocate(rhs_ED_MD)
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