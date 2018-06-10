      program Spectrums
      USE IFLPORT
      implicit none      
    
      character*300 dir_work, name_lambda, name_GR, name_dip, dir_tab
      character*100 N_name, b_name, dist_name, eps_name, pol_name 
      character*100 aspect_name, height_name, prop_name, particle_name, sub_name
      character*100 name_omega, core_name, theta_name, phi_name
     
      character*400 out_spec
      character*100 CLARG                                                           !command line argument
      character*15 TMPSTR, TMPSTR_FULL                                              !temporary string
      
      logical exist_flag_work, result
      
      integer*4 N, Na, Nline, i, j, k, alpha, beta, nl, il, ia, ja
      integer*4 nshape, naxis, nw, iw, excitation
      integer*4 substrate, mat_core_1, mat_core_2, mat_shell_1, mat_shell_2, mat_sub
      integer*4 lwork, liwork, info, istat
      integer*4 i_N_name, i_b_name, i_dist_name, i_eps_name
      integer*4 i_aspect_name, i_height_name, i_omega_name, i_core_name
      integer*4 i_dir_work, i_dir_tab
      integer*4 types
      integer :: output_data
      
      real*8 aspect, b, dist, height, eps_ext, fill, a_eff
      real*8 li, omega_min, omega_max, wi, dw
      real*8 theta_pol_ED, phi_pol_ED, theta_pol_MD, phi_pol_MD, theta_prop, phi_prop
      real*8 ru, rz, pi, twopi, mu0, eps0, mueps, mueps_sq, epsmu_sq
      real*8 wv, wv_0, wv_1, wv_2, wv_3
      real*8 xij, yij, zij, rij, rij_1, rij_2, rij_3, d_alpha_beta 
      real*8 rnijx, rnijy, rnijz, rnij_alpha, rnij_beta, rnij_alpha_beta
      real*8 t1_g, t2_g, t3_g, arg
      real*8 t1_c, t2_c
      real*8 prop(3), pol_elec(3), pol_magn(3)
      real*8 time_begin, time_end
      real*8 da, d2a
      real*8 arg_inc, arg_ref
      real*8 sum_e_ED, sum_a_ED, sum_s_ED, qe_ED, qa_ED, qs_ED, yabs, delta
      real*8 sum_e_MD, sum_a_MD, sum_s_MD, qe_MD, qa_MD, qs_MD
      real*8 sum_e_ED_MD, sum_a_ED_MD, sum_s_ED_MD, qe_ED_MD, qa_ED_MD, qs_ED_MD
      
      complex*16 cu, cz, phase_fact, ref
      complex*16 eps_core_1,eps_core_2, eps_shell_1, eps_shell_2, eps_sub, eps_1, eps_2
      complex*16 GRxx, GRyy, GRzz, GRxz, GRzx
      complex*16 dn0, cdx, cdy, cdz
      complex*16 sum_x, sum_y, sum_z
      
      integer*4,    allocatable, dimension(:)     :: iwork
      real*8,       allocatable, dimension(:)     :: x,y,z
      complex*16,   allocatable, dimension(:)     :: work, zet_ED, zet_MD
      complex*16,                dimension(3)     :: zet_ED_1, zet_MD_1,zet_ED_2, zet_MD_2
      complex*16,   allocatable, dimension(:,:)   :: a, rhs
      complex*16,   allocatable, dimension(:,:,:) :: GR      
      complex*16,   allocatable, dimension(:,:)   :: E_field, H_field, field_ED_MD
      complex*16,   allocatable, dimension(:,:)   :: a_ED, a_MD, a_ED_MD 
      complex*16,   allocatable, dimension(:,:)   :: rhs_ED_MD
      complex*16,   allocatable, dimension(:,:)   :: G_tenzor, C_tenzor
      
        interface
            subroutine coordinates (types, dist, height, x, y, z, N, Nline)    
                real*8,       allocatable, dimension(:) :: x, y, z
                integer*4 types, Nline, N
                real*8 dist, height
            end subroutine coordinates
        end interface
       
!c 
!c------------------ Reading parameters from 'inp.par' file ------------------------ 
!c  
      CALL CPU_TIME (time_begin)
      
      !++++++++++++++++++++++++++++++++
      !take parameters from command line
      !copy of input file
      !__________________
      !sample for commande line argument
      !all parameters from input.txt one by one
      !C:\...\dipole_app.exe C:\...\ C:\...\ 38.115 1. 900 838 0 1. 1 2 0.9 5.4 4501 90 90 0 0 1 2 2 3 2.25 0
      
      open(unit=70,file="inp.par",status='old')

      do
        read(70, *, end=10) TMPSTR_FULL

        if (TMPSTR_FULL(1:4) .eq. "#END" ) goto 10

        if (TMPSTR_FULL(2:8) .eq. "WORKING") read(70,*)         dir_work
        if (TMPSTR_FULL(2:4) .eq. "TAB") read(70,*)             dir_tab  
        if (TMPSTR_FULL(2:8) .eq. "SHORTER") read(70,*)         b
        if (TMPSTR_FULL(2:8) .eq. "FILLING") read(70,*)         fill
        if (TMPSTR_FULL(2:10).eq. "STRUCTURE") read(70,*)       types
        if (TMPSTR_FULL(2:7) .eq. "NUMBER") read(70,*)          N
        if (TMPSTR_FULL(2:14).eq. "INTERPARTICLE") read(70,*)   dist
        if (TMPSTR_FULL(2:7) .eq. "HEIGHT") read(70,*)          height
        if (TMPSTR_FULL(2:7) .eq. "ASPECT") read(70,*)          aspect
        if (TMPSTR_FULL(2:9) .eq. "SPHEROID") read(70,*)        nshape
        if (TMPSTR_FULL(2:5) .eq. "AXIS") read(70,*)            naxis
        if (TMPSTR_FULL(2:8) .eq. "OMEGA_1") read(70,*)         omega_min,  omega_max,  nw
        if (TMPSTR_FULL(2:7) .eq. "POL(E)") read(70,*)          theta_pol_ED,  phi_pol_ED
        if (TMPSTR_FULL(2:8) .eq. "PROP(K)") read(70,*)         theta_prop, phi_prop
        if (TMPSTR_FULL(2:11).eq. "EXCITATION") read(70,*)      excitation
        if (TMPSTR_FULL(2:5) .eq. "CORE")  read(70,*)           mat_core_1, mat_core_2
        if (TMPSTR_FULL(2:6) .eq. "SHELL") read(70,*)           mat_shell_1, mat_shell_2
        if (TMPSTR_FULL(2:9) .eq. "MATERIAL") read(70,*)        mat_sub
        if (TMPSTR_FULL(2:8) .eq. "EPSILON") read(70,*)         eps_ext
        if (TMPSTR_FULL(2:10).eq. "SUBSTRATE") read(70,*)       substrate
        !if (TMPSTR_FULL(2:5) .eq. "TYPE") read(70,*)            output_data
      end do
      
10    close(70)
!c
!c----------------------Creating output file-----------------------
!c            
      write (N_name, '(I4)') N 
      N_name = trim(adjustl(N_name))
      write (b_name, '(I4)') int(b) 
      b_name = trim(adjustl(b_name)) 
      write (height_name, '(I4)') int(height) 
      height_name = trim(adjustl(height_name))
      write (dist_name, '(I4)') int(dist) 
      dist_name = trim(adjustl(dist_name))
      write (eps_name, '(I4)') int(eps_ext*100) 
      eps_name = trim(adjustl(eps_name))
      write (aspect_name, '(I4)') int(aspect*10) 
      aspect_name = trim(adjustl(aspect_name))
            
      write (theta_name, '(I4)') int(theta_pol_ED) 
      theta_name = trim(adjustl(theta_name))
      write (phi_name, '(I4)') int(phi_pol_ED) 
      phi_name = trim(adjustl(phi_name))

      write (pol_name, 104) theta_name,phi_name 
  104 format(a<2>,'_',a<2>)   
      pol_name = trim(adjustl(pol_name))
      
      write (theta_name, '(I4)') int(theta_prop) 
      theta_name = trim(adjustl(theta_name))
      write (phi_name, '(I4)') int(phi_prop) 
      phi_name = trim(adjustl(phi_name))

      write (prop_name, 105) theta_name,phi_name 
  105 format(a<2>,'_',a<2>)   
      prop_name = trim(adjustl(prop_name))
            
      if (nshape .eq. 1 .and. aspect .ne. 1.) particle_name = 'PR'
      if (nshape .eq. 2 .and. aspect .ne. 1.) particle_name = 'OB'
      if (aspect .eq. 1.)                     particle_name = 'SP'
            
      if (mat_core_1 .eq. 1) core_name = 'AG'
      if (mat_core_1 .eq. 2) core_name = 'AU'
      if (mat_core_1 .eq. 3) core_name = 'TiN'
      if (mat_core_1 .eq. 4) core_name = 'ZrN'
      if (mat_core_1 .eq. 5) core_name = 'AZO'
      if (mat_core_1 .eq. 6) core_name = 'GZO'
      if (mat_core_1 .eq. 7) core_name = 'ITO'
      if (mat_core_1 .eq. 8) core_name = 'SI'
            
      if (substrate .eq. 0) sub_name = 'free'
      if (substrate .eq. 1) sub_name = 'subs'
            
      i_N_name = len_trim(N_name)
      i_b_name = len_trim(b_name)             
      i_dist_name = len_trim(dist_name) 
      i_height_name = len_trim(height_name)
      i_aspect_name = len_trim(aspect_name)
      i_eps_name = len_trim(eps_name)
      i_dir_work = len_trim(dir_work)
      i_core_name = len_trim(core_name)

      IF (aspect .ne. 1.) THEN
          write(out_spec, 102) dir_work, core_name, particle_name, aspect_name, N_name, b_name, dist_name, pol_name, prop_name, eps_name, sub_name
 102      format(a<i_dir_work>, a<i_core_name>, a<2>, '0', a<i_aspect_name>, '_N=', a<i_N_name>, '_b=', a<i_b_name>, '_h=',a<i_dist_name>, '_pol=',a<5>,'_prop=', a<1>,'_eps=', a<i_eps_name>, '_', a<4>,'.dat')
          ELSE 
              write(out_spec, 103) dir_work, core_name, particle_name, N_name, b_name, dist_name, pol_name, prop_name, eps_name, sub_name
 103          format(a<i_dir_work>,'2D_', a<i_core_name>, '_',a<2>, '_N=', a<i_N_name>, '_b=', a<i_b_name>, '_h=',a<i_dist_name>, '_pol=',a<5>,'_prop=', a<1>, '_eps=', a<i_eps_name>, '_', a<4>, '.dat')
      END IF
         
      i_dir_work = len_trim(dir_work)
      i_dir_tab  = len_trim(dir_tab)
!c 
!c--------- Constants -------------------------- 
!c 
      cz = (0.0d0,0.0d0) 
      cu = (0.0d0,1.0d0) 
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
!c-----------Calculating a_eff-----------------------
!c
      if (aspect .ne. 1) then
          if (nshape .eq. 1) a_eff = b / aspect ** (1./3.)
          if (nshape .eq. 2) a_eff = b / aspect ** (2./3.)
      end if
      if (aspect .eq. 1) a_eff = b     
!c
!c--------------------Converting degrees to radians---------------------
!c
      theta_pol_ED = theta_pol_ED * pi / 180.
      phi_pol_ED   = phi_pol_ED   * pi / 180.
      theta_pol_MD = theta_pol_MD * pi / 180.
      phi_pol_MD   = phi_pol_MD   * pi / 180.
      theta_prop   = theta_prop   * pi / 180.
      phi_prop     = phi_prop     * pi / 180.
!c      
!c--------------------Checking for input errors-------------------------      
!c
      if (nshape .ne. 1 .and. nshape .ne. 2) goto 202
      if (naxis .ne. 1 .and. naxis .ne. 2 .and. naxis .ne. 3) goto 203
      if (sin(theta_pol_ED)*cos(phi_pol_ED) * sin(theta_prop)*cos(phi_prop) + sin(theta_pol_ED)*sin(phi_pol_ED) * sin(theta_prop)*sin(phi_prop) + cos(theta_pol_ED)*cos(theta_prop) .gt. 1.E-10) goto 205
      
      inquire(DIRECTORY=dir_work, EXIST=exist_flag_work)
      if(.not. exist_flag_work) goto 206
!c 
!c------------------- Writing coordinates --------------- 
!c          
!c 
!c------------------- for single chain---------------       
!c      
      IF (types .eq. 1) THEN
         Nline = N
         allocate(x(1:N),y(1:N),z(1:N))
      ENDIF
!c
!c----------------------full lattice------------------------------
!c         0 0 0 0 0
!c         0 0 0 0 0
!c         0 0 0 0 0
!c         0 0 0 0 0
!c         0 0 0 0 0
!c-----------------------------------------------------------------
!c
      IF (types .eq. 2) THEN
          Nline = N
          N = N**2.
          allocate(x(1:N),y(1:N),z(1:N)) 
      ENDIF
!c
!c----------------------cut lattice------------------------------
!c         0 0 0 0 0
!c         0   0   0
!c         0 0 0 0 0
!c         0   0   0
!c         0 0 0 0 0
!c-----------------------------------------------------------------
!c
      IF (types .eq. 3) THEN
          Nline = N
          N = N**2 - ((N-1)/2)**2
          allocate(x(1:N),y(1:N),z(1:N))
      ENDIF
!c
!c----------------------cut lattice------------------------------
!c         0   0   0
!c         0   0   0
!c         0   0   0
!c         0   0   0
!c         0   0   0
!c-----------------------------------------------------------------
!c
      IF (types .eq. 4) THEN
          Nline = N
          N = (N**2 + N)/2
          allocate(x(1:N),y(1:N),z(1:N))
      ENDIF
!c
!c----------------------hex------------------------------
!c           0 0 0 0 0
!c            0 0 0 0 0
!c           0 0 0 0 0
!c            0 0 0 0 0
!c           0 0 0 0 0
!c-----------------------------------------------------------------
!c
      IF (types .eq. 5) THEN
          Nline = N
          N = N**2 
          allocate(x(1:N),y(1:N),z(1:N)) 
      ENDIF  
      
      x = 0.
      y = 0.
      z = 0.
      
      call coordinates(types, dist, height, x, y, z, N, Nline)

      write (*,*) '--------------------------------------------------'
      write (*,*) '2D lattice, NxN=', N
      write (*,*) 'b=', sngl(b)
      write (*,*) 'period=', sngl(dist)
      write (*,*) '--------------------------------------------------'
      !pause
!c      
!c-----------------------Allocating workspace-----------------------------      
!c      
       Na = 3*N
  
       allocate(E_field(1:Na,1:1), H_field(1:Na,1:1))
       allocate(a_ED_MD(1:2*Na,1:2*Na))
       
       allocate(zet_ED(1:Na))
       allocate(zet_MD(1:Na))
       
       allocate(field_ED_MD(1:2*Na,1:1))
       allocate(rhs_ED_MD(1:2*Na,1:1))
       allocate(a_ED(1:Na,1:Na), a_MD(1:Na,1:Na))
       allocate(G_tenzor(1:Na,1:Na), C_tenzor(1:Na,1:Na))
       
       prop(1) = dsin(theta_prop)*dcos(phi_prop) 
       prop(2) = dsin(theta_prop)*dsin(phi_prop)
       prop(3) = dcos(theta_prop)
       
       pol_elec(1)  = dsin(theta_pol_ED)*dcos(phi_pol_ED)
       pol_elec(2)  = dsin(theta_pol_ED)*dsin(phi_pol_ED)
       pol_elec(3)  = dcos(theta_pol_ED)
       
       pol_magn(1)  = prop(2)*pol_elec(3) - prop(3)*pol_elec(2)
       pol_magn(2)  = prop(3)*pol_elec(1) - prop(1)*pol_elec(3)
       pol_magn(3)  = prop(1)*pol_elec(2) - prop(2)*pol_elec(1)
    
      open(unit=70,file=out_spec, status='unknown')          
!c   
!c----------------------------------------------------------------------- 
!c---------------------START LOOP OVER OMEGA-----------------------------             
      do 15 iw = 1,nw
!c----------------------------------------------------------------------- 
!c-----------------------------------------------------------------------       
 
          dw = (omega_max - omega_min) / dfloat(nw-1)
          wi =  omega_min + dw * dfloat(iw-1)
      
          li = twopi * 3.0d2 / wi
      
          write(*,*) wi,li
      
          write (name_omega, '(f5.3)') wi 
          name_omega = trim(adjustl(name_omega)) 
          i_omega_name = len_trim(name_omega)

          wv = sqrt(eps_ext) * twopi / li
          wv_0 = ru 
          wv_1 = wv 
          wv_2 = wv*wv 
          wv_3 = wv_2*wv 
       
          call permittivity (mat_shell_1,li,eps_ext,eps_shell_1)
          call permittivity (mat_core_1, li,eps_ext,eps_core_1)
            eps_1 = eps_shell_1 * (eps_core_1 + 2.*eps_shell_1 + 2.*fill*(eps_core_1-eps_shell_1))/ (eps_core_1 + 2.*eps_shell_1 - fill*(eps_core_1-eps_shell_1))
            
          call permittivity (mat_shell_2,li,eps_ext,eps_shell_2)
          call permittivity (mat_core_2, li,eps_ext,eps_core_2)
            eps_2 = eps_shell_2 * (eps_core_2 + 2.*eps_shell_2 + 2.*fill*(eps_core_2-eps_shell_2))/ (eps_core_2 + 2.*eps_shell_2 - fill*(eps_core_2-eps_shell_2))
                  
          call susceptibility_single(eps_1, N, nshape, naxis, b, aspect, wv, zet_ED_1, zet_MD_1) 
          call susceptibility_single(eps_2, N, nshape, naxis, b, aspect, wv, zet_ED_2, zet_MD_2) 
          
          if (N .eq. 1) then
              zet_ED = zet_ED_1
              zet_MD = zet_MD_1
              else
                do i = 1, Na, 6
                    zet_ED (i)   = zet_ED_1(1)
                    zet_ED (i+1) = zet_ED_1(2)
                    zet_ED (i+2) = zet_ED_1(3)
                    zet_MD (i)   = zet_MD_1(1)
                    zet_MD (i+1) = zet_MD_1(2)
                    zet_MD (i+2) = zet_MD_1(3)
                if(i+3 .lt. Na) then
                    zet_ED (i+3) = zet_ED_2(1)
                    zet_ED (i+4) = zet_ED_2(2)
                    zet_ED (i+5) = zet_ED_2(3)
                    zet_MD (i+3) = zet_MD_2(1)
                    zet_MD (i+4) = zet_MD_2(2)
                    zet_MD (i+5) = zet_MD_2(3)
                end if
                end do
          end if
!c
!c----------- Creating matrix A --------------------------------- 
!c  
!__________________________________________electric part
          
          a_ED = cz 
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
          
!_for G_tenzor 
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
              
3         continue 
2         continue
          
!_for C_tenzor
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
31        continue 
21        continue
 
1     continue 
          
      do 4, i=2,Na 
          do 4, j=1,i-1 
              G_tenzor(j, i) = G_tenzor(i, j)
              C_tenzor(j, i) = C_tenzor(i, j)
 4    continue
 
      do 5, i=1,Na
          do 5, j=1,Na 
              if (i .eq. j) then        !for diagonal elements a_ED matrix
                  a_ED(j,i) = zet_ED(i)! - G_tenzor(j,i)
                  a_MD(j,i) = zet_MD(i)! - C_tenzor(j,i)
                  
              else                      !for off diagonal elements G_ED matrix
                  G_tenzor(j,i) = -1.0 * G_tenzor(j,i)
              end if
5     continue
       
        ! !print C_tenzor
        !open(unit=71,file="C_tenzor.txt",status='UNKNOWN')
        !do i = 1, 3 * N
        !        write (71, '(100(g12.5,2x))') (C_tenzor(i, j), j = 1 , 3 * N)
        !        write (71, "")
        !end do
        !close (unit=71) 
        !
        ! !print G_tenzor
        !open(unit=71,file="G_tenzor.txt",status='UNKNOWN')
        !do i = 1, 3 * N
        !        write (71, '(100(g12.5,2x))') (G_tenzor(i, j), j = 1 , 3 * N)
        !        write (71, "")
        !end do
        !close (unit=71)
!c 
!c-------- Costructing ED MD matrix A  -------------------------------- 
!c
    a_ED_MD = cz
    
    do i = 1, 3 * N 
        
        do j = 1, 3 * N
            
            if(i .eq. j) a_ED_MD(    i,     j) =  a_ED(i, j)                        !left up part           diagonal elements
            if(i .eq. j) a_ED_MD(3*N+i, 3*N+j) =  a_MD(i, j)                        !right down part        diagonal elements
            
            if(i .ne. j) a_ED_MD(    i,     j) =  G_tenzor(i, j)                    !left up part       off diagonal elements
            if(i .ne. j) a_ED_MD(3*N+i, 3*N+j) =  G_tenzor(i, j)                    !right down part    off diagonal elements
            if(i .ne. j) a_ED_MD(    i, 3*N+j) = -1.0 * epsmu_sq * C_tenzor(i, j)   !left down part     off diagonal elements 
            if(i .ne. j) a_ED_MD(3*N+i,     j) =  mueps_sq* C_tenzor(i, j)   !right up part      off diagonal elements
        
        end do
             
    end do
    
    !!print ED MD matrix
    !    open(unit=71,file="matrix_a_ED_MD.txt",status='UNKNOWN')
    !    do i = 1, 6 * N
    !            write (71, '(100(g12.5,2x))') (a_ED_MD(i, j), j = 1 , 6 * N)
    !    end do
    !    close (unit=71)     
!c 
!c-------- End of costructing ED MD matrix A  -------------------------------- 
!c
    
!c 
!c-------- Costructing right-hand side  -------------------------------- 
!c  

    E_field = cz
    H_field = cz
  
    if (excitation .eq. 1) then           !Plane wave
        
        do i=1,N
            arg_inc = prop(1) * x(i) + prop(2) * y(i) + prop(3) * z(i)
            !_________________________________for electric part 
            E_field(3*i-2,1) =            pol_elec(1) * cdexp(cu * wv_1 * arg_inc)
            E_field(3*i-1,1) =            pol_elec(2) * cdexp(cu * wv_1 * arg_inc)
            E_field(3*i  ,1) =            pol_elec(3) * cdexp(cu * wv_1 * arg_inc)
            !_________________________________for magnetic part 
            H_field(3*i-2,1) = epsmu_sq * pol_magn(1) * cdexp(cu * wv_1 * arg_inc)
            H_field(3*i-1,1) = epsmu_sq * pol_magn(2) * cdexp(cu * wv_1 * arg_inc)
            H_field(3*i  ,1) = epsmu_sq * pol_magn(3) * cdexp(cu * wv_1 * arg_inc)
        end do
   
    else                                  !Tip
        
        !_________________________________for electric part 
        E_field(1,1) =            pol_elec(1)  
        E_field(2,1) =            pol_elec(2) 
        E_field(3,1) =            pol_elec(3)
        !_________________________________for magnetic part
        H_field(1,1) = epsmu_sq * pol_magn(1)  
        H_field(2,1) = epsmu_sq * pol_magn(2) 
        H_field(3,1) = epsmu_sq * pol_magn(3)
        
    END IF    

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
      call zsysv('U',2*Na,1,a_ED_MD,2*Na,iwork,rhs_ED_MD,2*Na,work,-1,info) 
      lwork = work(1) 
      deallocate(work) 
      allocate(work(lwork))
      
      call zsysv('U',2*Na,1,a_ED_MD,2*Na,iwork,rhs_ED_MD,2*Na,work,lwork,info) 
 
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
      sum_e_ED_MD = rz
      sum_a_ED_MD = rz
      sum_s_ED_MD = rz
      
      sum_e_ED = rz
      sum_a_ED = rz
      sum_s_ED = rz
      
      sum_e_MD = rz
      sum_a_MD = rz
      sum_s_MD = rz
     
      do i=1,Na
      !ED+MD
      sum_e_ED_MD = sum_e_ED_MD + aimag(rhs_ED_MD(i,1)*dconjg(E_field(i,1))) + mueps*aimag(rhs_ED_MD(Na+i,1)*dconjg(H_field(i,1)))
      !ED
      sum_e_ED = sum_e_ED + aimag(rhs_ED_MD(i,1) * dconjg(E_field(i,1))) 
      !MD
      sum_e_MD = sum_e_MD + mueps*aimag(rhs_ED_MD(Na+i,1)*dconjg(H_field(i,1)))
      end do      
     
      !if (nshape .eq. 1 .and. aspect .ne. ru) delta = -dimag(aspect * (eps_ED + 2.*ru) / (eps_ext * b**3 * (eps_ED - ru)))
     
      !if (nshape .eq. 2 .and. aspect .ne. ru) delta = -dimag(aspect**2 * (eps_ED + 2.*ru) / (eps_ext * b**3 * (eps_ED - ru)))
      
      qe_ED_MD = 4.*wv_1*sum_e_ED_MD/(float(N)*a_eff**2.)
      
      qe_ED = 4.*wv_1*sum_e_ED/(float(N)*a_eff**2.)
      
      qe_MD = 4.*wv_1*sum_e_MD/(float(N)*a_eff**2)
      
      !ED_MD
      write(70,*)  sngl(li), sngl(wi), sngl(qe_ED_MD), sngl(qe_ED), sngl(qe_MD)

!c---------------------------------------------------------------------                  
!c---------------------------------------------------------------------      
 15   end do 
!c---------------------END LOOP OVER OMEGA-----------------------------  
!c--------------------------------------------------------------------- 
!c       
      deallocate(a_ED_MD)
      deallocate(rhs_ED_MD)
      close(70)
      close(80)
 
      CALL CPU_TIME (time_end)
      write(*,*) 'estimated_time=', time_end-time_begin, 'sec'
!c-----------------------------------------------------------------------      
 500  stop 
!c-------------------------ERRORS----------------------------------------  
  
 202  continue
      write(*,*) 'Invalid nshape; exiting'
      stop
 
 203  continue
      write(*,*) 'Invalid naxis; exiting'
      stop     

 205  continue
      write(*,*) 'E and k are not orthogonal; exiting'
      stop          

 206  continue
      i_dir_work = len_trim(dir_work)
      write(*,*) 'The directory to store dipole moments'
      write(*,300) dir_work
 300  format(a<i_dir_work>, ' does not exist or is unaccessable.')
      stop
      
 209  continue
      write(*,*) 'Unable to open binary; exiting'
      stop

 211  continue
      write(*,*) 'Unable to store dipole moments; exiting'
      stop
      
      end