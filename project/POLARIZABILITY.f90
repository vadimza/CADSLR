    SUBROUTINE POLARIZABILITY (eps, rad, wv, zet_e, zet_m) 

    implicit none
    real*8,    intent(in)  :: rad, wv
    complex*16,intent(in)  :: eps
    complex*16,intent(out) :: zet_e(3), zet_m(3)

    real*8 wv_1, wv_2, wv_3
    real*8 ru, pi, twopi
    complex*16 cu, m
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
        
    m = cdsqrt(eps)
    psi_mkr_0 = cdsin(m*wv_1*rad)/(m*wv_1*rad) - cdcos(m*wv_1*rad)
    psi_mkr_1 = cdcos(m*wv_1*rad)/(m*wv_1*rad) - cdsin(m*wv_1*rad)/((m*wv_1*rad)**2) + cdsin(m*wv_1*rad)
     
    psi_kr_0 = dsin(wv_1*rad)/(wv_1*rad) - dcos(wv_1*rad)
    psi_kr_1 = dcos(wv_1*rad)/(wv_1*rad) - dsin(wv_1*rad)/((wv_1*rad)**2.) + dsin(wv_1*rad)
     
    xi_kr_0 = psi_kr_0 - cu * (dcos(wv_1*rad)/(wv_1*rad) + dsin(wv_1*rad))
    xi_kr_1 = psi_kr_1 + cu * (dsin(wv_1*rad)/(wv_1*rad) + dcos(wv_1*rad)/((wv_1*rad)**2) - dcos(wv_1*rad))
      
    zet_e = - 2. * cu * wv_3 / 3. * (m*psi_mkr_0*xi_kr_1 -   xi_kr_0*psi_mkr_1) / &
    (m*psi_mkr_0*psi_kr_1 -   psi_kr_0*psi_mkr_1)
    zet_m = - 2. * cu * wv_3 / 3. * (  psi_mkr_0*xi_kr_1 - m*xi_kr_0*psi_mkr_1) / &
    (  psi_mkr_0*psi_kr_1 - m*psi_kr_0*psi_mkr_1)
      
    RETURN
      
    END SUBROUTINE