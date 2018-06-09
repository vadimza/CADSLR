SUBROUTINE SUSCEPTIBILITY (eps, N, nshape, naxis, b, aspect, wv, zet_ED, zet_MD) 
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
    complex*16,intent(out) :: zet_ED(3*N), zet_MD(3*N)
    integer*4 i
    complex*16 cu
    complex :: m
    real*8 D_axis, D_ort
    real*8 wv_1, wv_2, wv_3
    real*8 ecc_1, ecc_2, ge 
    real*8 ru, pi, twopi
    real*8 L_axis, L_ort
    
    complex :: psi_kr_0, psi_kr_1
    complex :: psi_mkr_0, psi_mkr_1
    complex :: xi_kr_0, xi_kr_1
    
    cu = (0.0d0,1.0d0) 
    ru = 1.0d0 
    pi = 4.0d0*datan(ru) 
    twopi = 2*pi 
      
    wv_1 = wv
    wv_2 = wv_1**2
    wv_3 = wv_1**3        

    !!------------------------------PROLATE-----------------------------------
    !IF (nshape .eq. 1 .and. aspect .ne. ru) THEN
    !  
    !    ecc_2 = 1 - aspect**2
    !    ecc_1 = sqrt(ecc_2)      
    !
    !    L_axis = (1-ecc_2)*(0.5*dlog((1+ecc_1)/(1-ecc_1))/ecc_1-ru)/ecc_2
    !    L_ort = (1 - L_axis)/2.
    !  
    !    D_axis = 3*(1+ecc_2)/(1-ecc_2)*L_axis/4 + 1
    !    D_ort = 0.5*aspect*(1.5*dlog((1+ecc_1)/(1-ecc_1))/ecc_1 - D_axis)
    !
    !    zet_ED = 3 * aspect * (L_ort + ru / (eps - ru)) / b**3 - wv_2 / b * D_ort - 2. * cu * wv_3 / 3.
    !  
    !    do i=0,N-1
    !            
    !        zet_ED(naxis+3*i) = 3 * aspect * (L_axis + ru / (eps - ru)) / b**3 - wv_2 / (b / aspect) * D_axis - 2. * cu * wv_3 / 3.
    !  
    !    end do
    !  
    !END IF
    !!------------------------------OBLATE------------------------------------
    !IF (nshape .eq. 2 .and. aspect .ne. ru) THEN
    !  
    !    ecc_2 = 1 - aspect**2
    !    ecc_1 = sqrt(ecc_2)  
    !    ge = sqrt((1-ecc_2)/ecc_2)
    !  
    !    L_ort = 0.5*ge*(pi/2.-datan(ge))/ecc_2 - ge**2/2.
    !    L_axis = 1 - 2.*L_ort
    !  
    !    D_axis = 3*(1-2*ecc_2)*L_axis/4 + 1
    !    D_ort = 0.5*(3*ge*dasin(ecc_1) - D_axis)/aspect
    !
    !    zet_ED = 3 * aspect**2 * (L_ort + ru / (eps - ru)) / b**3 - wv_2 / (b / aspect) * D_ort - 2. * cu *wv_3 / 3.
    !
    !    do i=0,N-1
    !  
    !        zet_ED(naxis+3*i) = 3 * aspect**2 * (L_axis + ru / (eps - ru)) / b**3 - wv_2 / b * D_axis - 2. * cu * wv_3 / 3.
    !  
    !    end do
    !   
    !END IF
    !------------------------------SPHERE------------------------------------
    IF (aspect .eq. ru) THEN
        
      m = cdsqrt(eps)
      psi_mkr_0 = cdsin(m*wv_1*b)/(m*wv_1*b) - cdcos(m*wv_1*b)
      psi_mkr_1 = cdcos(m*wv_1*b)/(m*wv_1*b) - cdsin(m*wv_1*b)/((m*wv_1*b)**2) + cdsin(m*wv_1*b)
     
      psi_kr_0 = dsin(wv_1*b)/(wv_1*b) - dcos(wv_1*b)
      psi_kr_1 = dcos(wv_1*b)/(wv_1*b) - dsin(wv_1*b)/((wv_1*b)**2.) + dsin(wv_1*b)
     
      xi_kr_0 = psi_kr_0 - cu * (dcos(wv_1*b)/(wv_1*b) + dsin(wv_1*b))
      xi_kr_1 = psi_kr_1 + cu * (dsin(wv_1*b)/(wv_1*b) + dcos(wv_1*b)/((wv_1*b)**2) - dcos(wv_1*b))
      
      zet_ED = - 2. * cu * wv_3 / 3. * (m*psi_mkr_0*xi_kr_1 -   xi_kr_0*psi_mkr_1) / (m*psi_mkr_0*psi_kr_1 -   psi_kr_0*psi_mkr_1)
      zet_MD = - 2. * cu * wv_3 / 3. * (  psi_mkr_0*xi_kr_1 - m*xi_kr_0*psi_mkr_1) / (  psi_mkr_0*psi_kr_1 - m*psi_kr_0*psi_mkr_1)
        
    END IF
      
    RETURN
      
END SUBROUTINE