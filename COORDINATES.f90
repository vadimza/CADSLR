SUBROUTINE COORDINATES (types, dist, height, x, y, z, N, Nline)    
    implicit none
    integer*4 types, Nline, N
    real*8 dist, height, rz
    real*8, dimension(:) :: x,y,z
    integer :: i, j, k
       
    rz = 0.0d0
    
!c 
!c------------------- 2d chain --------------- 
!c         0 0 0 0 0
!c-----------------------------------------------------------------
!c 
      
      IF (types .eq. 1) THEN
        do j = 1, Nline
              x(j) = (j-1) * dist 
              y(j ) = rz!(j-1) * dist
              z(j ) = height
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
          do j = 1, Nline
              do i = 1, Nline 
                  x(i + Nline*(j-1)) = (i-1) * dist 
                  y(i + Nline*(j-1)) = (j-1) * dist
                  z(i + Nline*(j-1)) = height
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

          k = 1
          
          do j = 1, Nline, 2
              do i = 1, Nline
                  x(k) = (i-1) * dist 
                  y(k) = (j-1) * dist
                  z(k) = height
                  k = k + 1
              end do
          end do

          do j = 2, Nline-1, 2
              do i = 1, Nline, 2 
                  x(k) = (i-1) * dist 
                  y(k) = (j-1) * dist
                  z(k) = height
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

          k = 1
          
          do j = 1, Nline, 2
              do i = 1, Nline
                  x(k) = (i-1) * dist 
                  y(k) = (j-1) * dist
                  z(k) = height
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

        k =1
        
        do i = 1, Nline, 2 
            do j = 1, Nline 
                x(k) = (j-1) * dist
                y(k) = (i-1) * 0.5 * dist * 3**0.5
                z(k) = height 
                k = k + 1
            end do 
        end do 

        do i = 2, Nline-1, 2 
            do j = 1, Nline 
                x(k) = (j-1) * dist + 0.5 * dist 
                y(k) = (i-1) * 0.5 * dist * 3**0.5 
                z(k) = height
                k = k +1
            end do 
        end do 
      ENDIF

return
END SUBROUTINE