    PROGRAM CREATE_INPUT
!------------------------------------------------------------------------------
!   DEFINING VARIABLES
!------------------------------------------------------------------------------
    implicit none
!   CHARACTERS    
    character*50 dist_name, dist_namex, dist_namey
    character*50 dist_namexx, dist_nameyy
    character*50 sig_pos_name, sig_siz_name
    character*50 phiE_name, phik_name, thetak_name
    character*50 n_namex, n_namey, ndel_name
    character*50 chain, excitation_type
    character*50 input_file, output_file, structure_file
!   INTEGERS    
    integer*4 Nline1, Nline2, N 
    integer*4 Ndel1, Ndel2, indel
    integer*4 xx1, xx2, yy1, yy2
    integer*4 sig_pos1, sig_pos2, isig_pos
    integer*4 sig_siz1, sig_siz2, isig_siz
    integer*8 distx1, distx2, idistx
    integer*8 disty1, disty2, idisty
    integer*8 lambda1, lambda2, lambdan
    integer*8 theta_E1, theta_E2, phi_E1, phi_E2, iphiE
    integer*8 theta_k1, theta_k2, phi_k1, phi_k2, iphik, ithetak
    integer :: i, j, id, k, ii, jj
    integer*8 i_dist_name, i_dist_namex, i_dist_namey
    integer*8 i_dist_namexx, i_dist_nameyy    
    integer*8 i_sig_pos_name, i_sig_siz_name
    integer*8 i_n_namex, i_n_namey, i_ndel_name
    integer*8 i_phiE_name, i_phik_name, i_thetak_name
!   REALS    
    real*8 host_medium, rz, rand
    real*8 radius
!   ARRAYS    
    integer*4,    allocatable, dimension(:) :: n_dis_ind
    real*8,       allocatable, dimension(:) :: x, y, z, rad
    character*10, allocatable, dimension(:) :: mat
!------------------------------------------------------------------------------
!   CONSTANT PARAMETERS    
!------------------------------------------------------------------------------
!   ARRAY GEOMETRY
    parameter (chain='2D')          ! type of chain
    parameter (Nline1 = 10)               ! particles in line X
    parameter (Nline2 = 10)               ! particles in line Y
    parameter (distx1 = 540)            ! period X
    parameter (distx2 = 540)
    parameter (disty1 = 400)            ! period Y
    parameter (disty2 = 500)
!   PARTICLES
    parameter (radius = 65)             ! radius, nm
!   ILLUMINATION
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
!   DEFECTS
    parameter (Ndel1 = 0)               ! delete particles
    parameter (Ndel2 = 0)
    parameter (sig_pos1 = 0)            ! positional disorder
    parameter (sig_pos2 = 0)
    parameter (sig_siz1 = 0)            ! size disorder
    parameter (sig_siz2 = 0)
!   SURROUNDING MEDIUM
    parameter (host_medium = 1)
!------------------------------------------------------------------------------
!   START    
!------------------------------------------------------------------------------
    call init_random_seed

    rz = 0.0d0
    id = 0
    
    open(unit=80,file='run.bat',status='unknown')!,position='append')
    
!    do indel = ndel1, ndel2, 200
    
    do idistx = distx1,distx2,5
    do idisty = disty1,disty2,5
    do iphik  = phi_k1,phi_k2,1
    do ithetak = theta_k1,theta_k2,2
    
!    do iphiE  = phi_E1,phi_E2,2
        
!    do isig_pos = sig_pos1,sig_pos2,50
!    do isig_siz = sig_siz1,sig_siz2,5
        
        id = id + 1
        
        write (n_namex, '(I4)') Nline1 
        n_namex   = trim(adjustl(n_namex))
        i_n_namex = len_trim(n_namex)
        
        write (n_namey, '(I4)') Nline2 
        n_namey   = trim(adjustl(n_namey))
        i_n_namey = len_trim(n_namey)
                
        write (dist_namex, '(I4)') idistx 
        dist_namex   = trim(adjustl(dist_namex))
        i_dist_namex = len_trim(dist_namex)
        
        write (dist_namey, '(I4)') idisty 
        dist_namey   = trim(adjustl(dist_namey))
        i_dist_namey = len_trim(dist_namey)
        
!        write (dist_namexx, '(I5)') idistxx 
!        dist_namexx   = trim(adjustl(dist_namexx))
!        i_dist_namexx = len_trim(dist_namexx)
!        
!        write (dist_nameyy, '(I5)') idistyy 
!        dist_nameyy   = trim(adjustl(dist_nameyy))
!        i_dist_nameyy = len_trim(dist_nameyy)
        
        write (sig_pos_name, '(I4)') isig_pos 
        sig_pos_name   = trim(adjustl(sig_pos_name))
        i_sig_pos_name = len_trim(sig_pos_name)
        
        write (sig_siz_name, '(I4)') isig_siz 
        sig_siz_name   = trim(adjustl(sig_siz_name))
        i_sig_siz_name = len_trim(sig_siz_name)
        
        write (ndel_name, '(I4)') indel 
        ndel_name   = trim(adjustl(ndel_name))
        i_ndel_name = len_trim(ndel_name)
        
        write (phiE_name, '(I4)') iphiE 
        phiE_name   = trim(adjustl(phiE_name))
        i_phiE_name = len_trim(phiE_name)
        
        !write (phik_name, '(I4)') iphik 
        !phik_name   = trim(adjustl(phik_name))
        !i_phik_name = len_trim(phik_name)
        !
        !write (thetak_name, '(I4)') ithetak 
        !thetak_name   = trim(adjustl(thetak_name))
        !i_thetak_name = len_trim(thetak_name)                         
        
        i_dist_name = i_n_namex + i_n_namey + i_dist_namex + i_dist_namey + 3
        
        write(dist_name, 100) n_namex, n_namey, dist_namex, dist_namey
100     format(a<i_n_namex>,'_',a<i_n_namey>,'_',a<i_dist_namex>,'_',a<i_dist_namey>)

!        if (chain .eq. '1D') then
!            N = Nline
!            allocate(x(1:N),y(1:N),z(1:N),rad(1:N))
!            rad = radius
!            do i = 1, Nline 
!                x(i) = (i-1) * idistx 
!                y(i) = rz
!                z(i) = rz
!            enddo
!        endif
        
!        if (chain .eq. '1D_dis') then
!            N = Nline
!            allocate(x(1:N),y(1:N),z(1:N),rad(1:N))
!            rad = radius
!            do i = 1, Nline
!                call RANDOM_NUMBER(rand) 
!                x(i) = (i-1) * idistx + 2 * (rand - 0.5) * isig_pos
!                call RANDOM_NUMBER(rand)
!                y(i) = 2 * (rand - 0.5) * isig_pos
!                z(i) = rz
!            enddo
!        endif

        if (chain .eq. '2D') then
            N = Nline1*Nline2
            allocate(x(1:N),y(1:N),z(1:N),rad(1:N))
            rad = radius
            do j = 1, Nline2
                do i = 1, Nline1 
                    x(i + Nline1*(j-1)) = (i-1) * idistx 
                    y(i + Nline1*(j-1)) = (j-1) * idisty
                    z(i + Nline1*(j-1)) = rz
                enddo
            enddo
        endif
        
!        if (chain .eq. '2D_del') then
!            N = Nline**2 - indel
!            allocate(n_dis_ind(1:N))
!            do i = 1, N
!11              call RANDOM_NUMBER(rand)
!                n_dis_ind(i) = int(Nline**2*rand)+1
!                if (i .ne. 1) then
!                    do j = 1, i-1
!                        if (n_dis_ind(i) .eq. 0) goto 11
!                        if (n_dis_ind(i) .eq. n_dis_ind(j)) goto 11
!                    enddo
!                endif
!                goto 12
!12          enddo
!
!            allocate(x(1:N),y(1:N),z(1:N),rad(1:N))
!            rad = radius
!            do j = 1, Nline
!                do i = 1, Nline
!                    do k = 1,N
!                        if (i + Nline*(j-1) .eq. n_dis_ind(k)) then
!                            x(k) = (i-1) * idistx 
!                            y(k) = (j-1) * idisty
!                            z(k) = rz
!                        endif
!                    enddo
!                enddo
!            enddo
!        endif        
!        
!        if (chain .eq. '2D_dis') then
!            N = Nline**2
!            allocate(x(1:N),y(1:N),z(1:N),rad(1:N))
!            do j = 1, Nline
!                do i = 1, Nline 
!                    call RANDOM_NUMBER(rand) 
!                    x(i + Nline*(j-1)) = (i-1) * idistx + 2 * (rand - 0.5) * isig_pos
!                    call RANDOM_NUMBER(rand) 
!                    y(i + Nline*(j-1)) = (j-1) * idisty + 2 * (rand - 0.5) * isig_pos
!                    z(i + Nline*(j-1)) = rz
!                    call RANDOM_NUMBER(rand)
!                    rad(i + Nline*(j-1)) = radius + 2 * (rand - 0.5) * isig_siz
!                enddo
!            enddo
!        endif            
        
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
        write(70,*) theta_E1+ithetak, iphiE
        
        write(70,*) 'prop_k_theta_phi'
        write(70,*) ithetak, iphik        
        
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
        
        if (mod(id,4) .eq. 0)  write(80,*) 'START /WAIT ','CADSLR.exe ', input_file
        if (mod(id,4) .ne. 0)  write(80,*) 'START ','CADSLR.exe ', input_file
        !write(80,*) './cadslr.out ', input_file
        
    enddo
    enddo
    enddo
    enddo
!    enddo
!    enddo    
!    enddo
    
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
