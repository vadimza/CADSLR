SUBROUTINE SUSCEPTIBILITY (eps, N, nshape, naxis, b, aspect, wv, zet) 
!
!   CALCULATION OF THE RECIPROCAL POLARIZABILITY ZET FOR PARTICLES OF SPHERICAL
!   AND SPHEROIDAL SHAPES. NOTE THAT THE INPUT REQUIRES NORMALIZED EPSILON WITH
!   RESPECT TO THE SURROUNDING MEDIUM. CALCULATIONS ARE BASED ON:
!
!   [1] Alexander Moroz, "Depolarization field of spheroidal particles," 
!                         J. Opt. Soc. Am. B 26, 517-527 (2009)
!   
    implicit none
    integer*4, intent(in) :: N, nshape, naxis 
    real*8,    intent(in) :: aspect, b, wv
    complex*16,intent(in) :: eps
    complex*16,intent(out) :: zet(3*N)
    integer*4 i
    complex*16 cu
    real*8 D_axis, D_ort
    real*8 wv_1, wv_2, wv_3
    real*8 ecc_1, ecc_2, ge 
    real*8 ru, pi, twopi
    real*8 L_axis, L_ort
      
    cu = (0.0d0,1.0d0) 
    ru = 1.0d0 
    pi = 4.0d0*datan(ru) 
    twopi = 2*pi 
      
    wv_1 = wv
    wv_2 = wv_1**2
    wv_3 = wv_1**3        

    !------------------------------PROLATE-----------------------------------
    IF (nshape .eq. 1 .and. aspect .ne. ru) THEN
      
        ecc_2 = 1 - aspect**2
        ecc_1 = sqrt(ecc_2)      

        L_axis = (1-ecc_2)*(0.5*dlog((1+ecc_1)/(1-ecc_1))/ecc_1-ru)/ecc_2
        L_ort = (1 - L_axis)/2.
      
        D_axis = 3*(1+ecc_2)/(1-ecc_2)*L_axis/4 + 1
        D_ort = 0.5*aspect*(1.5*dlog((1+ecc_1)/(1-ecc_1))/ecc_1 - D_axis)

        zet = 3 * aspect * (L_ort + ru / (eps - ru)) / b**3 - wv_2 / b * D_ort - 2. * cu * wv_3 / 3.
      
        do i=0,N-1
                
            zet(naxis+3*i) = 3 * aspect * (L_axis + ru / (eps - ru)) / b**3 - wv_2 / (b / aspect) * D_axis - 2. * cu * wv_3 / 3.
      
        end do
      
    END IF
    !------------------------------OBLATE------------------------------------
    IF (nshape .eq. 2 .and. aspect .ne. ru) THEN
      
        ecc_2 = 1 - aspect**2
        ecc_1 = sqrt(ecc_2)  
        ge = sqrt((1-ecc_2)/ecc_2)
      
        L_ort = 0.5*ge*(pi/2.-datan(ge))/ecc_2 - ge**2/2.
        L_axis = 1 - 2.*L_ort
      
        D_axis = 3*(1-2*ecc_2)*L_axis/4 + 1
        D_ort = 0.5*(3*ge*dasin(ecc_1) - D_axis)/aspect
  
        zet = 3 * aspect**2 * (L_ort + ru / (eps - ru)) / b**3 - wv_2 / (b / aspect) * D_ort - 2. * cu *wv_3 / 3.

        do i=0,N-1
      
            zet(naxis+3*i) = 3 * aspect**2 * (L_axis + ru / (eps - ru)) / b**3 - wv_2 / b * D_axis - 2. * cu * wv_3 / 3.
      
        end do
       
    END IF
    !------------------------------SPHERE------------------------------------
    IF (aspect .eq. ru) THEN

        zet = (eps + 2.*ru) / (eps - ru) / b**3 - wv_2 / b - 2. * cu * wv_3 / 3.
      
    END IF
      
    RETURN
      
END SUBROUTINE