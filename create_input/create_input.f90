    PROGRAM CREATE_INPUT

    implicit none
    
    character*50 dist_name, dist_namex, dist_namey
    character*50 dist_namexx, dist_nameyy
    character*50 sig_pos_name, sig_siz_name
    character*50 n_name, ndel_name
    character*50 chain, excitation_type
    character*50 input_file, output_file, structure_file
    
    integer*4 Nline, Npatch, N 
    integer*4 Ndel1, Ndel2, indel
    integer*4 xx1, xx2, yy1, yy2
    integer*4 sig_pos1, sig_pos2, isig_pos
    integer*4 sig_siz1, sig_siz2, isig_siz
    integer*8 distx1, distx2, idistx
    integer*8 disty1, disty2, idisty
    integer*8 distxx1, distxx2, idistxx
    integer*8 distyy1, distyy2, idistyy
    integer*8 lambda1, lambda2, lambdan
    integer*8 theta_E1, theta_E2, phi_E1, phi_E2
    integer*8 theta_k1, theta_k2, phi_k1, phi_k2
    integer :: i, j, id, k, ii, jj
    integer*8 i_dist_name, i_dist_namex, i_dist_namey
    integer*8 i_dist_namexx, i_dist_nameyy    
    integer*8 i_sig_pos_name, i_sig_siz_name
    integer*8 i_n_name, i_ndel_name
    
    real*8 host_medium, rz, rand
    real*8 radius
    
    integer*4,    allocatable, dimension(:) :: n_dis_ind
    real*8,       allocatable, dimension(:) :: x,y,z,rad
    character*10, allocatable, dimension(:) :: mat
    
    parameter (Nline = 30)               ! particles in line
    parameter (Npatch = 10)             ! number of patches
    parameter (chain='2D')          ! type of chain
    parameter (radius = 50)             ! radius, nm
    
    parameter (Ndel1 = 0)               ! delete particles
    parameter (Ndel2 = 0)
    parameter (sig_pos1 = 0)            ! positional disorder
    parameter (sig_pos2 = 0)
    parameter (sig_siz1 = 0)            ! size disorder
    parameter (sig_siz2 = 0)
    
    parameter (distx1 = 400)            ! period X
    parameter (distx2 = 700)
    parameter (disty1 = 400)            ! period Y
    parameter (disty2 = 700)
    
    parameter (distxx1 = 2000)            ! period X large
    parameter (distxx2 = 10000)
    parameter (distyy1 = 2000)            ! period Y large
    parameter (distyy2 = 10000)

    parameter (lambda1 = 400)           ! lambda, nm
    parameter (lambda2 = 700)
    parameter (lambdan = 301)
    
    parameter (theta_E1 = 90)
    parameter (theta_E2 = 90)
    parameter (phi_E1 = 0)
    parameter (phi_E2 = 0)
    parameter (theta_k1 = 0)
    parameter (theta_k2 = 0)
    parameter (phi_k1 = 0)
    parameter (phi_k2 = 0)
    parameter (excitation_type='plane')
    
    parameter (host_medium = 1)

    call init_random_seed

    rz = 0.0d0
    id = 0
    
    open(unit=80,file='run.bat',status='unknown')!,position='append')
    
    do indel = ndel1, ndel2, 200
    
    do idistx = distx1,distx2,10
    do idisty = disty1,disty2,10
        
    !do idistxx = distxx1,distxx2,500
    !do idistyy = distyy1,distyy2,500        

    do isig_pos = sig_pos1,sig_pos2,50
    do isig_siz = sig_siz1,sig_siz2,5
        
        id = id + 1
        
        write (n_name, '(I4)') Nline 
        n_name   = trim(adjustl(n_name))
        i_n_name = len_trim(n_name)        
        
        write (dist_namex, '(I4)') idistx 
        dist_namex   = trim(adjustl(dist_namex))
        i_dist_namex = len_trim(dist_namex)
        
        write (dist_namey, '(I4)') idisty 
        dist_namey   = trim(adjustl(dist_namey))
        i_dist_namey = len_trim(dist_namey)
        
        write (dist_namexx, '(I5)') idistxx 
        dist_namexx   = trim(adjustl(dist_namexx))
        i_dist_namexx = len_trim(dist_namexx)
        
        write (dist_nameyy, '(I5)') idistyy 
        dist_nameyy   = trim(adjustl(dist_nameyy))
        i_dist_nameyy = len_trim(dist_nameyy)
        
        write (sig_pos_name, '(I4)') isig_pos 
        sig_pos_name   = trim(adjustl(sig_pos_name))
        i_sig_pos_name = len_trim(sig_pos_name)
        
        write (sig_siz_name, '(I4)') isig_siz 
        sig_siz_name   = trim(adjustl(sig_siz_name))
        i_sig_siz_name = len_trim(sig_siz_name)
        
        write (ndel_name, '(I4)') indel 
        ndel_name   = trim(adjustl(ndel_name))
        i_ndel_name = len_trim(ndel_name)                
        
!        i_dist_name = i_dist_namex + i_dist_namey + i_dist_namexx + i_dist_nameyy + 7
!        
!        write(dist_name, 100) dist_namex, dist_namey, dist_namexx, dist_nameyy
!100     format('x',a<i_dist_namex>,'_y',a<i_dist_namey>,'_X',a<i_dist_namexx>,'_Y',a<i_dist_nameyy>)

        i_dist_name = i_n_name + i_dist_namex + i_dist_namey + 2
        
        write(dist_name, 100) n_name, dist_namex, dist_namey
100     format(a<i_n_name>,'_',a<i_dist_namex>,'_',a<i_dist_namey>)


        if (chain .eq. '1D') then
            N = Nline
            allocate(x(1:N),y(1:N),z(1:N),rad(1:N))
            rad = radius
            do i = 1, Nline 
                x(i) = (i-1) * idistx 
                y(i) = rz
                z(i) = rz
            enddo
        endif
        
        if (chain .eq. '1D_dis') then
            N = Nline
            allocate(x(1:N),y(1:N),z(1:N),rad(1:N))
            rad = radius
            do i = 1, Nline
                call RANDOM_NUMBER(rand) 
                x(i) = (i-1) * idistx + 2 * (rand - 0.5) * isig_pos
                call RANDOM_NUMBER(rand)
                y(i) = 2 * (rand - 0.5) * isig_pos
                z(i) = rz
            enddo
        endif

        if (chain .eq. '2D') then
            N = Nline**2
            allocate(x(1:N),y(1:N),z(1:N),rad(1:N))
            rad = radius
            do j = 1, Nline
                do i = 1, Nline 
                    x(i + Nline*(j-1)) = (i-1) * idistx 
                    y(i + Nline*(j-1)) = (j-1) * idisty
                    z(i + Nline*(j-1)) = rz
                enddo
            enddo
        endif
        
        if (chain .eq. '2D_H') then
            N = Nline**2 * Npatch**2
            allocate(x(1:N),y(1:N),z(1:N),rad(1:N))
            rad = radius
            do jj = 1, Npatch
                do ii = 1, Npatch
                    do j = 1, Nline
                        do i = 1, Nline
                            x(i + Nline*(j-1) + Npatch*Nline**2*(jj-1) + Nline**2*(ii-1)) = (i-1) * idistx + (ii-1) * idistxx
                            y(i + Nline*(j-1) + Npatch*Nline**2*(jj-1) + Nline**2*(ii-1)) = (j-1) * idisty + (jj-1) * idistyy
                            z(i + Nline*(j-1) + Npatch*Nline**2*(jj-1) + Nline**2*(ii-1)) = rz
                        enddo
                    enddo
                enddo
            enddo
        endif

        if (chain .eq. '2D_del') then
            N = Nline**2 - indel
            allocate(n_dis_ind(1:N))
            do i = 1, N
11              call RANDOM_NUMBER(rand)
                n_dis_ind(i) = int(Nline**2*rand)+1
                if (i .ne. 1) then
                    do j = 1, i-1
                        if (n_dis_ind(i) .eq. 0) goto 11
                        if (n_dis_ind(i) .eq. n_dis_ind(j)) goto 11
                    enddo
                endif
                goto 12
12          enddo

            allocate(x(1:N),y(1:N),z(1:N),rad(1:N))
            rad = radius
            do j = 1, Nline
                do i = 1, Nline
                    do k = 1,N
                        if (i + Nline*(j-1) .eq. n_dis_ind(k)) then
                            x(k) = (i-1) * idistx 
                            y(k) = (j-1) * idisty
                            z(k) = rz
                        endif
                    enddo
                enddo
            enddo
        endif        
        
        if (chain .eq. '2D_dis') then
            N = Nline**2
            allocate(x(1:N),y(1:N),z(1:N),rad(1:N))
            do j = 1, Nline
                do i = 1, Nline 
                    call RANDOM_NUMBER(rand) 
                    x(i + Nline*(j-1)) = (i-1) * idistx + 2 * (rand - 0.5) * isig_pos
                    call RANDOM_NUMBER(rand) 
                    y(i + Nline*(j-1)) = (j-1) * idisty + 2 * (rand - 0.5) * isig_pos
                    z(i + Nline*(j-1)) = rz
                    call RANDOM_NUMBER(rand)
                    rad(i + Nline*(j-1)) = radius + 2 * (rand - 0.5) * isig_siz
                enddo
            enddo
        endif            
        
        write(structure_file, 101) dist_name
101     format(a<i_dist_name>,'.coord')
        
        write(input_file, 102) dist_name
102     format(a<i_dist_name>,'.inp')

        write(output_file, 103) dist_name
103     format(a<i_dist_name>,'.out')

!        write(structure_file, 101) sig_pos_name
!101     format(a<i_sig_pos_name>,'.coord')
!        
!        write(input_file, 102) sig_pos_name
!102     format(a<i_sig_pos_name>,'.inp')
!
!        write(output_file, 103) sig_pos_name
!103     format(a<i_sig_pos_name>,'.out')
      
        open(unit=70,file=structure_file,status='unknown')
    
        do i = 1,N
            write(70,*) sngl(x(i)), sngl(y(i)), sngl(z(i)), sngl(rad(i)), 'Si'
        enddo
        
        close(unit=70)
        deallocate(x,y,z,rad)
        !,n_dis_ind)
        
        open(unit=70,file=input_file,status='unknown')
    
        write(70,*) 'number_of_particles'
        write(70,*) N
        
        write(70,*) 'lambda_1'
        write(70,*) lambda1

        write(70,*) 'lambda_2'
        write(70,*) lambda2
        
        write(70,*) 'lambda_n'
        write(70,*) lambdan        
        
        write(70,*) 'pol_E_theta_phi'
        write(70,*) theta_E1, phi_E1
        
        write(70,*) 'prop_k_theta_phi'
        write(70,*) theta_k1, phi_k1        
        
        write(70,*) 'excitation_type'
        write(70,*) excitation_type
        
        write(70,*) 'host_medium'
        write(70,*) host_medium
        
        write(70,*) 'structure_file'
        write(70,*) structure_file
        
        write(70,*) 'output_file'
        write(70,*) output_file        
        
        write(70,*) 'end_of_options'
        
        close(unit=70)
        
!        if (mod(id,8) .eq. 0)  write(80,*) 'START /WAIT ','CADSLR.exe ', input_file
!        if (mod(id,8) .ne. 0)  write(80,*) 'START ','CADSLR.exe ', input_file
        write(80,*) './cadslr.out ', input_file
        
    enddo
    enddo
    enddo
    enddo    
    enddo
    !enddo
    !enddo
    
    close(80)
    
    END
    
    SUBROUTINE init_random_seed()
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
       
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
          
    CALL SYSTEM_CLOCK(COUNT=clock)
          
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
          
    DEALLOCATE(seed)
    END SUBROUTINE

!!c 
!!c------------------- for single chain---------------       
!!c      
!      IF (types .eq. 1) THEN
!         Nline = N
!         allocate(x(1:N),y(1:N),z(1:N))
!      ENDIF
!!c
!!c----------------------full lattice------------------------------
!!c         0 0 0 0 0
!!c         0 0 0 0 0
!!c         0 0 0 0 0
!!c         0 0 0 0 0
!!c         0 0 0 0 0
!!c-----------------------------------------------------------------
!!c
!      IF (types .eq. 2) THEN
!          Nline = N
!          N = N**2.
!          allocate(x(1:N),y(1:N),z(1:N)) 
!      ENDIF
!!c
!!c----------------------cut lattice------------------------------
!!c         0 0 0 0 0
!!c         0   0   0
!!c         0 0 0 0 0
!!c         0   0   0
!!c         0 0 0 0 0
!!c-----------------------------------------------------------------
!!c
!      IF (types .eq. 3) THEN
!          Nline = N
!          N = N**2 - ((N-1)/2)**2
!          allocate(x(1:N),y(1:N),z(1:N))
!      ENDIF
!!c
!!c----------------------cut lattice------------------------------
!!c         0   0   0
!!c         0   0   0
!!c         0   0   0
!!c         0   0   0
!!c         0   0   0
!!c-----------------------------------------------------------------
!!c
!      IF (types .eq. 4) THEN
!          Nline = N
!          N = (N**2 + N)/2
!          allocate(x(1:N),y(1:N),z(1:N))
!      ENDIF
!!c
!!c----------------------hex------------------------------
!!c           0 0 0 0 0
!!c            0 0 0 0 0
!!c           0 0 0 0 0
!!c            0 0 0 0 0
!!c           0 0 0 0 0
!!c-----------------------------------------------------------------
!!c
!      IF (types .eq. 5) THEN
!          Nline = N
!          N = N**2 
!          allocate(x(1:N),y(1:N),z(1:N)) 
!      ENDIF  
!      
!!      x = 0.
!!      y = 0.
!!      z = 0.
!
!
!    
