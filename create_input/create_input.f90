    PROGRAM CREATE_INPUT
    use iflport
    implicit none
    
    character*4 dist_name
    character*20 chain, excitation_type
    character*20 input_file, output_file, structure_file
    
    integer*4 Nline, N
    integer*8 dist1, dist2, idist
    integer*8 lambda1, lambda2, lambdan
    integer*8 theta_E1, theta_E2, phi_E1, phi_E2
    integer*8 theta_k1, theta_k2, phi_k1, phi_k2
    integer :: i, j
    integer*8 i_dist_name
    
    real*8 host_medium, rz
    real*8 radius
    
    real*8,       allocatable, dimension(:) :: x,y,z,rad
    character*10, allocatable, dimension(:) :: mat
    
    parameter (Nline = 20)
    parameter (chain='2D')
    parameter (radius = 80)
    parameter (dist1 = 400)
    parameter (dist2 = 410)
    parameter (lambda1 = 400)
    parameter (lambda2 = 800)
    parameter (lambdan = 401)
    parameter (theta_E1 = 90)
    parameter (theta_E2 = 90)
    parameter (phi_E1 = 0)
    parameter (phi_E2 = 0)
    parameter (theta_k1 = 0)
    parameter (theta_k2 = 0)
    parameter (phi_k1 = 0)
    parameter (phi_k2 = 0)
    parameter (excitation_type='plane')
    parameter (host_medium = 1.0)
    
    rz = 0.0d0
    
    open(unit=80,file='run.bat',status='unknown')
    
    do idist = dist1,dist2
        
        write (dist_name, '(I4)') idist 
        dist_name   = trim(adjustl(dist_name))
        i_dist_name = len_trim(dist_name)
        
        if (chain .eq. '2D') then
            N = Nline**2
            allocate(x(1:N),y(1:N),z(1:N),rad(1:N))
            rad = radius
            do j = 1, Nline
                do i = 1, Nline 
                    x(i + Nline*(j-1)) = (i-1) * idist 
                    y(i + Nline*(j-1)) = (j-1) * idist
                    z(i + Nline*(j-1)) = rz
                enddo
            enddo
        endif
        
        write(structure_file, 100) dist_name
100     format(a<i_dist_name>,'.coord')
        
        write(input_file, 101) dist_name
101     format(a<i_dist_name>,'.inp')

        write(output_file, 102) dist_name
102     format(a<i_dist_name>,'.out')
        
        open(unit=70,file=structure_file,status='unknown')
    
        do i = 1,N-1,2
        
            write(70,*) sngl(x(i)), sngl(y(i)), sngl(z(i)), sngl(rad(i)), 'Au'
            write(70,*) sngl(x(i+1)), sngl(y(i+1)), sngl(z(i+1)), sngl(rad(i+1)), 'Si' 
        
        enddo
        
        close(unit=70)
        deallocate(x,y,z,rad)
        
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
        
        if (mod(idist,4) .eq. 0)  write(80,*) 'START /WAIT ','CADSLR.exe ', input_file
        if (mod(idist,4) .ne. 0)  write(80,*) 'START ','CADSLR.exe ', input_file
        
    enddo
    
    close(80)
    
    END

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
!    !c 
!!c------------------- 1D CHAIN--------------- 
!!c         0 0 0 0 0
!!c-----------------------------------------------------------------
!!c 
!      IF (types .eq. 1) THEN
!        do j = 1, Nline
!              x(j) = (j-1) * dist 
!              y(j ) = rz!(j-1) * dist
!              z(j ) = height
!        end do      
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
!          do j = 1, Nline
!              do i = 1, Nline 
!                  x(i + Nline*(j-1)) = (i-1) * dist 
!                  y(i + Nline*(j-1)) = (j-1) * dist
!                  z(i + Nline*(j-1)) = height
!              end do
!          end do
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
!
!          k = 1
!          
!          do j = 1, Nline, 2
!              do i = 1, Nline
!                  x(k) = (i-1) * dist 
!                  y(k) = (j-1) * dist
!                  z(k) = height
!                  k = k + 1
!              end do
!          end do
!
!          do j = 2, Nline-1, 2
!              do i = 1, Nline, 2 
!                  x(k) = (i-1) * dist 
!                  y(k) = (j-1) * dist
!                  z(k) = height
!                  k = k + 1
!              end do
!          end do
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
!          k = 1
!          
!          do j = 1, Nline, 2
!              do i = 1, Nline
!                  x(k) = (i-1) * dist 
!                  y(k) = (j-1) * dist
!                  z(k) = height
!                  k = k + 1
!              end do
!          end do
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
!        k =1
!        
!        do i = 1, Nline, 2 
!            do j = 1, Nline 
!                x(k) = (j-1) * dist
!                y(k) = (i-1) * 0.5 * dist * 3**0.5
!                z(k) = height 
!                k = k + 1
!            end do 
!        end do 
!
!        do i = 2, Nline-1, 2 
!            do j = 1, Nline 
!                x(k) = (j-1) * dist + 0.5 * dist 
!                y(k) = (i-1) * 0.5 * dist * 3**0.5 
!                z(k) = height
!                k = k +1
!            end do 
!        end do 
!      ENDIF