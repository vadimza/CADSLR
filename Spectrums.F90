      program Spectrums
      USE IFLPORT
      implicit none
      
      character*300 dir_work, name_lambda, name_GR, name_dip, dir_tab
      character*100 N_name, b_name, dist_name, eps_name, pol_name 
      character*100 aspect_name, height_name, prop_name, particle_name, sub_name
      character*100 name_omega, core_name, theta_name, phi_name
     
      character*400 out_spec
      character*100 CLARG !command line argument
      character*15 TMPSTR, TMPSTR_FULL !temporary string
      
      logical exist_flag_work, result
      
      integer*4 N, Na, Nline, i, j, k, alpha, beta, nl, il, ia, ja
      integer*4 nshape, naxis, nw, iw, excitation
      integer*4 substrate, mat_core, mat_shell, mat_sub
      integer*4 lwork, liwork, info, istat
      integer*4 i_N_name, i_b_name, i_dist_name, i_eps_name
      integer*4 i_aspect_name, i_height_name, i_omega_name, i_core_name
      integer*4 i_dir_work, i_dir_tab
      integer*4 types
      
      real*8 aspect, b, dist, height, eps_ext, fill, a_eff
      real*8 li, omega_min, omega_max, wi, dw
      real*8 theta_pol, phi_pol, theta_prop, phi_prop
      real*8 ru, rz, pi, twopi
      real*8 wv, wv_0, wv_1, wv_2, wv_3
      real*8 xij, yij, zij, rij, rij_1, rij_2, rij_3, d_alpha_beta 
      real*8 rnijx, rnijy, rnijz, rnij_alpha, rnij_beta, rnij_alpha_beta
      real*8 t1, t2, t3, arg
      real*8 prop(3), pol(3)
      real*8 time_begin, time_end
      real*8 da, d2a
      real*8 arg_inc, arg_ref
      real*8 sum_e, sum_a, sum_s, qe, qa, qs, yabs, delta
      
      complex*16 cu, cz, phase_fact, ref
      complex*16 eps_core, eps_shell, eps_sub, eps_eff
      complex*16 GRxx, GRyy, GRzz, GRxz, GRzx
      complex*16 dn0, cdx, cdy, cdz
      complex*16 sum_x, sum_y, sum_z
      
      integer*4,    allocatable, dimension(:)   :: iwork
      real*8,       allocatable, dimension(:)   :: x,y,z
      complex*16,   allocatable, dimension(:)   :: work, zet
      complex*16,   allocatable, dimension(:,:) :: a, rhs, E
      complex*16,   allocatable, dimension(:,:,:) :: GR     
      
!c 
!c------------------ Reading parameters from 'inp.par' file ------------------------ 
!c  
      CALL CPU_TIME (time_begin)
      
      !++++++++++++++++++++++++++++++++
      !take parameters from command line
      !copy of input file
      !__________________
!#WORKING DIRECTORY
!C:\git_repositories\CADSLR\results\
!#TAB DIRECTORY
!C:\git_repositories\CADSLR\results\
!#SHORTER SEMIAXIS OR RADIUS [NM]	                                        !	IF ASPECT RATIO = 1, THEN THIS IS THE RADIUS OF A SPHERE
!38.115	
!#FILLING FACTOR					
!1.
!#TYPE OF STRUCTURE														!1 - CHAIN, 2 - SQUARE, 3 - CUT FIRST, 4 - CUT ROW, 5 - HEX
!1
!#NUMBER OF PARTICLES													! INTERPARTICLE DISTANCE IN CASE OF 1D CHAIN OR NUMBER OF PARTICLES IN ROW(COLUMN)
!9	
!#INTERPARTICLE DISTANCE
!838	
!#HEIGHT			                                                        ! 	DISTANCE BETWEEN CENTER OF PARTICLE AND SURFACE OF SUBSTRATE
!0
!#ASPECT RATIO			                                                !	0 < ASPECT RATIO <= 1; 1 - SPHERE
!.5
!#SPHEROID TYPE			                                                ! 	1 - PROLATE, 2 - OBLATE [ONLY IF ASPECT RATIO < 1]
!1
!#AXIS OF SYMMETRY		                                                ! 	1 - X, 2 - Y, 3 - Z	[ONLY IF ASPECT RATIO < 1]
!2
!#OMEGA_1 [rad/fs], OMEGA_2 [rad/fs], NUMBER OF POINTS					
!.9	5.4	4501
!#POL(E) THETA, PHI [DEG]	                                                ! 	ANGLE (POL^Z) & ANGLE (POL_PROJECTION^X)      {X [90;0], Y [90;90], Z [0;0]}
!90	90	
!#PROP(K) THETA, PHI [DEG]	                                                !	ANGLE (PROP^Z) & ANGLE (PROP_PROJECTION^X)    {X [90;0], Y [90;90], Z [0;0]}
!0	0	
!#EXCITATION TYPE			                                        !	0 - TIP / 1 - PLANE WAVE
!1	
!#CORE MATERIAL 			                                                !	1 - Ag / 2 - Au / 3 - TiN 800 / 4 - ZrN / 5 - AZO / 6 - GZO / 7 - ITO / 8 - CUSTOM
!2
!#SHELL MATERIAL			                                                !	1 - Ag / 2 - Au / 3 - TiN 800 / 4 - ZrN / 5 - AZO / 6 - GZO / 7 - ITO / 8 - CUSTOM
!2
!#MATERIAL OF SUBSTRATE		                                                !	1 - Ag / 2 - Au / 3 - TiN 800 / 4 - ZrN / 5 - AZO / 6 - GZO / 7 - ITO / 8 - CUSTOM
!2
!#EPSILON OF HOST MEDIUM		                                                !	1. - VACUUM / 1.78 - WATER
!2.25
!#SUBSTRATE 			                                                !	0-OFF / 1-ON
!0
!#END
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
        if (TMPSTR_FULL(2:5) .eq. "TYPE") read(70,*)            types
        if (TMPSTR_FULL(2:7) .eq. "NUMBER") read(70,*)          N
        if (TMPSTR_FULL(2:14) .eq."INTERPARTICLE") read(70,*)   dist
        if (TMPSTR_FULL(2:7) .eq. "HEIGHT") read(70,*)          height
        if (TMPSTR_FULL(2:7) .eq. "ASPECT") read(70,*)          aspect
        if (TMPSTR_FULL(2:9) .eq. "SPHEROID") read(70,*)        nshape
        if (TMPSTR_FULL(2:5) .eq. "AXIS") read(70,*)            naxis
        if (TMPSTR_FULL(2:8) .eq. "OMEGA_1") read(70,*)         omega_min,  omega_max,  nw
        if (TMPSTR_FULL(2:7) .eq. "POL(E)") read(70,*)          theta_pol,  phi_pol
        if (TMPSTR_FULL(2:8) .eq. "PROP(K)") read(70,*)         theta_prop, phi_prop
        if (TMPSTR_FULL(2:11) .eq."EXCITATION") read(70,*)      excitation
        if (TMPSTR_FULL(2:6) .eq. "CORE") read(70,*)            mat_core
        if (TMPSTR_FULL(2:7) .eq. "SHELL") read(70,*)           mat_shell
        if (TMPSTR_FULL(2:9) .eq. "MATERIAL") read(70,*)        mat_sub
        if (TMPSTR_FULL(2:8) .eq. "EPSILON") read(70,*)         eps_ext
        if (TMPSTR_FULL(2:10) .eq."SUBSTRATE") read(70,*)       substrate
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
            
      write (theta_name, '(I4)') int(theta_pol) 
      theta_name = trim(adjustl(theta_name))
      write (phi_name, '(I4)') int(phi_pol) 
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
            
      if (mat_core .eq. 1) core_name = 'AG'
      if (mat_core .eq. 2) core_name = 'AU'
      if (mat_core .eq. 3) core_name = 'TiN'
      if (mat_core .eq. 4) core_name = 'ZrN'
      if (mat_core .eq. 5) core_name = 'AZO'
      if (mat_core .eq. 6) core_name = 'GZO'
      if (mat_core .eq. 7) core_name = 'ITO'
      if (mat_core .eq. 8) core_name = 'CUS'
            
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
!c--------- Mathematical constants -------------------------- 
!c 
      cz = (0.0d0,0.0d0) 
      cu = (0.0d0,1.0d0) 
      rz = 0.0d0 
      ru = 1.0d0 
      pi = 4.0d0*datan(ru) 
      twopi = 2.0d0*pi
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
      theta_pol = theta_pol * pi / 180.
      phi_pol = phi_pol * pi / 180.
      theta_prop = theta_prop * pi / 180.
      phi_prop = phi_prop * pi / 180.
!c      
!c--------------------Checking for input errors-------------------------      
!c
      if (nshape .ne. 1 .and. nshape .ne. 2) goto 202
      if (naxis .ne. 1 .and. naxis .ne. 2 .and. naxis .ne. 3) goto 203
      if (sin(theta_pol)*cos(phi_pol) * sin(theta_prop)*cos(phi_prop) + sin(theta_pol)*sin(phi_pol) * sin(theta_prop)*sin(phi_prop) + cos(theta_pol)*cos(theta_prop) .gt. 1.E-10) goto 205
      
      inquire(DIRECTORY=dir_work, EXIST=exist_flag_work)
      if(.not. exist_flag_work) goto 206
      
      if (substrate .ne. 0 .and. substrate .ne. 1) goto 207      

   

!c 
!c------------------- Writing coordinates for single chain--------------- 
!c   
      write (*,*) '--------------------------------------------------'
      write (*,*) '2D lattice, NxN=', N,'x',N
      write (*,*) 'b=', sngl(b)
      write (*,*) 'period=', sngl(dist)
      write (*,*) '--------------------------------------------------'

     
!c 
!c------------------- for single chain---------------       
!c      
      if (types .eq. 1) THEN
      allocate(x(1:N),y(1:N),z(1:N))
        do j = 1, N
          do i = 1, N 
              x(i + N*(j-1)) = (i-1) * dist 
              y(i + N*(j-1)) = (j-1) * dist
              z(i + N*(j-1)) = height
          end do
        end do      
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
          N = N**2
          allocate(x(1:N),y(1:N),z(1:N)) 
          
          do j = 1, N
              do i = 1, N 
                  x(i + N*(j-1)) = (i-1) * dist 
                  y(i + N*(j-1)) = (j-1) * dist
                  z(i + N*(j-1)) = rz
              end do
          end do
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
          N = N**2 - ((N-1)/2)**2
          allocate(x(1:N),y(1:N),z(1:N))
          k = 1
          
          do j = 1, N, 2
              do i = 1, N
                  x(k) = (i-1) * dist 
                  y(k) = (j-1) * dist
                  z(k) = rz
                  k = k + 1
              end do
          end do

          do j = 2, N-1, 2
              do i = 1, N, 2 
                  x(k) = (i-1) * dist 
                  y(k) = (j-1) * dist
                  z(k) = rz
                  k = k + 1
              end do
          end do
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
          N = (N**2 + N)/2
          allocate(x(1:N),y(1:N),z(1:N))
          k = 1
          
          do j = 1, N, 2
              do i = 1, N
                  x(k) = (i-1) * dist 
                  y(k) = (j-1) * dist
                  z(k) = rz
                  k = k + 1
              end do
          end do
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
        N = N**2 
        allocate(x(1:N),y(1:N),z(1:N)) 
        k =1
        
        do i = 1, N, 2 
            do j = 1, N 
                x(k) = (j-1) * dist
                y(k) = (i-1) * 0.5 * dist * 3**0.5
                z(k) = rz 
                k = k + 1
            end do 
        end do 

        do i = 2, N-1, 2 
            do j = 1, N 
                x(k) = (j-1) * dist + 0.5 * dist 
                y(k) = (i-1) * 0.5 * dist * 3**0.5 
                z(k) = rz
                k = k +1
            end do 
        end do 
      ENDIF
    
!c      
!c-----------------------Allocating workspace-----------------------------      
!c      

       Na = 3*N
       
       allocate(a(1:Na,1:Na))
       allocate(zet(1:Na))
       allocate(E(1:Na,1:1),rhs(1:Na,1:1))
       allocate(GR(1:N,1:N,1:5)) 
  
       prop(1) = dsin(theta_prop)*dcos(phi_prop) 
       prop(2) = dsin(theta_prop)*dsin(phi_prop)
       prop(3) = dcos(theta_prop)
       pol(1)  = dsin(theta_pol)*dcos(phi_pol)
       pol(2)  = dsin(theta_pol)*dsin(phi_pol)
       pol(3)  = dcos(theta_pol)

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
       
          call permittivity (mat_shell,li,eps_ext,eps_shell)
          call permittivity (mat_core, li,eps_ext,eps_core)
      
          if (fill .ne. ru) then
              eps_eff = eps_shell * (eps_core + 2.*eps_shell + 2.*fill*(eps_core-eps_shell))/ (eps_core + 2.*eps_shell - fill*(eps_core-eps_shell))
          else 
              eps_eff = eps_core
          end if

          call susceptibility(eps_eff,N,nshape,naxis,b,aspect,wv,zet)
!c
!c----------- Creating matrix A --------------------------------- 
!c  
          a = cz 
          GR = cz
       
          IF (substrate .eq. 1 ) THEN
      
          call permittivity (mat_sub,li,eps_ext,eps_sub)
      
          write(name_GR, 121) dir_tab, name_omega 
 121      format(a<i_dir_tab>,'\','GR_', a<i_omega_name>, '.dat') 

          do i = 1,N-1
              do j = i+1,N
                  xij = x(j) - x(i) 
                  yij = y(j) - y(i) 
                  zij = z(j) - z(i) 
                  rij_2 = xij**2 + yij**2 + zij**2 
                  rij = dsqrt(rij_2) 
            
                  CALL GF_substrate_sameheight(height,rij,dist,wv_1,eps_sub, GRxx, GRyy, GRzz, GRxz, GRzx, istat)
            
                  GR(i,j,1) = GRxx * wv_3
                  GR(i,j,2) = GRyy * wv_3
                  GR(i,j,3) = GRzz * wv_3
                  GR(i,j,4) = GRxz * wv_3
                  GR(i,j,5) = GRzx * wv_3
              end do
          end do
          
          END IF
      
          close(10)
      
      do 1, i=1,N-1
      do 1, j=i+1,N
      
          xij = x(j) - x(i) 
          yij = y(j) - y(i) 
          zij = z(j) - z(i) 
          rij_2 = xij**2 + yij**2 + zij**2 
          rij = dsqrt(rij_2) 
          rij_1 = rij 
          rij_3 = rij * rij_2 
          rnijx = xij/rij 
          rnijy = yij/rij 
          rnijz = zij/rij 
          arg = wv * rij 
          phase_fact = cdexp(cu*arg)
      
          do 2, alpha=1,3 
              ia = 3*(i-1) + alpha 
              if (alpha .eq. 1) rnij_alpha = rnijx 
              if (alpha .eq. 2) rnij_alpha = rnijy 
              if (alpha .eq. 3) rnij_alpha = rnijz 
            
          do 3, beta=1,3 
              ja = 3*(j-1) + beta 
              if (beta .eq. 1) rnij_beta = rnijx 
              if (beta .eq. 2) rnij_beta = rnijy 
              if (beta .eq. 3) rnij_beta = rnijz 
              rnij_alpha_beta = rnij_alpha*rnij_beta 
              if(alpha .eq. beta) then 
                 d_alpha_beta = ru 
              else 
                 d_alpha_beta = rz 
              end if 
              t1 = wv_2 * (d_alpha_beta - 1.0*rnij_alpha_beta)/rij_1 
              t2 = wv_1 * (d_alpha_beta - 3.0*rnij_alpha_beta)/rij_2 
              t3 = wv_0 * (d_alpha_beta - 3.0*rnij_alpha_beta)/rij_3 
       
              a(ja,ia) = (t1 + cu*t2 - t3)*phase_fact

              IF (substrate .eq. 1 ) THEN
      
                  if (alpha .eq. 1 .and. beta .eq. 1)  a(ja,ia) = a(ja,ia) + GR(i,j,1)
                  if (alpha .eq. 2 .and. beta .eq. 2)  a(ja,ia) = a(ja,ia) + GR(i,j,2)
                  if (alpha .eq. 3 .and. beta .eq. 3)  a(ja,ia) = a(ja,ia) + GR(i,j,3)
                  if (alpha .eq. 1 .and. beta .eq. 3)  a(ja,ia) = a(ja,ia) + GR(i,j,4)
                  if (alpha .eq. 3 .and. beta .eq. 1)  a(ja,ia) = a(ja,ia) + GR(i,j,5)
          
              END IF
                 
 3        continue 
 2        continue
 
 1    continue 
       
      do 4, i=2,Na 
          do 4, j=1,i-1 
              a(j,i) = a(i,j) 
 4    continue
 
      do 5, i=1,Na
          do 5, j=1,Na 
              if (i .eq. j) then 
                  a(j,i) = zet(i) - a(j,i)
              else 
                  a(j,i) =        - a(j,i) 
              end if
 5    continue
 
      IF (substrate .eq. 1) THEN

      do 6, i=1,N
      
      do 7, alpha=1,3 
      ja = 3*(i-1) + alpha 
      
        if (alpha .eq. 1) a(ja,ja) = a(ja,ja) - GR(1,1,1)
        if (alpha .eq. 2) a(ja,ja) = a(ja,ja) - GR(1,1,2)
        if (alpha .eq. 3) a(ja,ja) = a(ja,ja) - GR(1,1,3)
     
 7    continue
      
 6    continue
 
      END IF
      
!      sum_x = cz
!      sum_y = cz
!      sum_z = cz
!      
!      do i = 4,Na
!      
!      sum_x = sum_x + a(i,1)
!      sum_y = sum_y + a(i,2)
!      sum_z = sum_z + a(i,3)
!      
!      end do
!      
!      sum_x = - 2*b**3 * sum_x
!      sum_y = - 2*b**3 * sum_y
!      sum_z = - 2*b**3 * sum_z
!    
!      write(80,*) sngl(li), sngl(wi), sngl(dreal(zet(1))*b**3),
!     $  sngl(dreal(sum_y)), sngl(dimag(zet(1)*b**3 - sum_y))
!      
!      goto 15
!c 
!c-------- Costructing right-hand side  -------------------------------- 
!c  
      E = cz
      
      IF (excitation .eq. 1) THEN           !Plane wave
      
      do i=1,N
      
      arg_inc = prop(1) * x(i) + prop(2) * y(i) + prop(3) * z(i)
           
      E(3*i-2,1) = pol(1) * cdexp(cu * wv_1 * arg_inc)
      E(3*i-1,1) = pol(2) * cdexp(cu * wv_1 * arg_inc)
      E(3*i  ,1) = pol(3) * cdexp(cu * wv_1 * arg_inc)
            
      if (substrate .eq. 1) then
      
      ref = (cdsqrt(eps_sub) - sqrt(eps_ext)) / (sqrt(eps_ext) + cdsqrt(eps_sub))
      
      arg_ref = prop(1) * x(i) + prop(2) * y(i) + prop(3) * (z(i) - height)
      
      E(3*i-2,1) = E(3*i-2,1) + ref*pol(1) * cdexp(cu * wv_1 * arg_ref)
      E(3*i-1,1) = E(3*i-1,1) + ref*pol(2) * cdexp(cu * wv_1 * arg_ref)
      E(3*i  ,1) = E(3*i  ,1) + ref*pol(3) * cdexp(cu * wv_1 * arg_ref)
      
      end if
      
      end do
            
      ELSE                                  !Tip
      
      E(1,1) = pol(1)  
      E(2,1) = pol(2) 
      E(3,1) = pol(3)
      
      if (substrate .eq. 1) then
      
      ref = (cdsqrt(eps_sub) - sqrt(eps_ext)) / (sqrt(eps_ext) + cdsqrt(eps_sub))
     
      E(1,1) = E(1,1) + ref*pol(1)
      E(2,1) = E(2,1) + ref*pol(2)
      E(3,1) = E(3,1) + ref*pol(3)
      
      end if
      
      END IF
            
      rhs = E      
!c 
!c--------- Solving equations --------------------------------------- 
!c  
      liwork = Na 
      allocate(iwork(liwork)) 
 
      allocate(work(1)) 
      call zsysv('U',Na,1,a,Na,iwork,rhs,Na,work,-1,info) 
      lwork = work(1) 
      deallocate(work) 
      allocate(work(lwork))
      
      call zsysv('U',Na,1,a,Na,iwork,rhs,Na,work,lwork,info) 
 
      if(info .ne. 0) then 
         write(*,*) 'Matrix inversion has failed; exiting' 
         write(*,*) 'info=', info 
         stop 
      end if
       
      deallocate(work) 
      deallocate(iwork)
!c 
!c---------- Saving out dipole moments------------------------------- 
!c  
!      write(name_dip, 122) dir_work, name_omega 
! 122  format(a<i_dir_work>,'\','dip_', a<i_omega_name>, 'nm.dat') 
!      
!      open(unit=20, file=name_dip, status='unknown', err=211)
!      
!      dn0 = rhs(1,1)*dconjg(rhs(1,1))+rhs(2,1)*dconjg(rhs(2,1))
!     $ + rhs(3,1)*dconjg(rhs(3,1))
!      
!      do 8, i=1,N 
!      cdx = rhs(3*i-2,1)*dconjg(rhs(3*i-2,1))
!      cdy = rhs(3*i-1,1)*dconjg(rhs(3*i-1,1))
!      cdz = rhs(3*i  ,1)*dconjg(rhs(3*i  ,1))
!      d2a = cdx + cdy + cdz
!      d2a = d2a / dn0
!      cdx = cdx / dn0
!      cdy = cdy / dn0
!      cdz = cdz / dn0
!      write(20,*)  i, sngl(d2a), sngl(cdabs(cdx)),  
!     $                          sngl(cdabs(cdy)), sngl(cdabs(cdz))
!
! 8    continue 
! 
!      close(20)
!
!      write(70,*)  sngl(li), sngl(d2a)!, 10d0*dlog10(d2a)
!c 
!c---------- Calculating spectrums------------------------------- 
!c  
      sum_e = rz
      sum_a = rz
      sum_s = rz
     
      do i=1,Na
      
      sum_e = sum_e + aimag(rhs(i,1) * dconjg(E(i,1))) 
      sum_a = sum_a + rhs(i,1) * dconjg(rhs(i,1))
      sum_s = sum_s + rhs(i,1) * dconjg(rhs(i,1))
      
      end do
      
      if (aspect .eq. 1) delta = -dimag((eps_eff + 2.*ru) / (eps_ext * b**3 * (eps_eff - ru)))
     
      if (nshape .eq. 1 .and. aspect .ne. ru) delta = -dimag(aspect * (eps_eff + 2.*ru) / (eps_ext * b**3 * (eps_eff - ru)))
     
      if (nshape .eq. 2 .and. aspect .ne. ru) delta = -dimag(aspect**2 * (eps_eff + 2.*ru) / (eps_ext * b**3 * (eps_eff - ru)))
      
      qe = 4.*wv_1*sum_e/(float(N)*a_eff**2)
      qa = 4.*wv_1*delta*sum_a/(float(N)*a_eff**2)
      qs = qe - qa
      
      write(70,*)  sngl(li), sngl(wi), sngl(qe), sngl(qa), sngl(qs)
!c---------------------------------------------------------------------                  
!c---------------------------------------------------------------------      
 15   end do 
!c---------------------END LOOP OVER OMEGA-----------------------------  
!c--------------------------------------------------------------------- 
!c       
      deallocate(a)
      deallocate(rhs)
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
      
 207  continue
      write(*,*) 'Variable for substrate is invalid; exiting' 
      stop  
      
 209  continue
      write(*,*) 'Unable to open binary; exiting'
      stop

 211  continue
      write(*,*) 'Unable to store dipole moments; exiting'
      stop
      
      end