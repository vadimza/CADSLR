 SUBROUTINE BJ01(X,BJ0,BJ1)
!C       ============================================================
!C        Computes BJ0 = J0(x) and BJ1 = J1(x)/x
!C        WARNING: do not use this subroutine in other applications.
!C                 It will not return correct values for the case X=0.
!C                 Also, the variable BJ1 contains J1(x)/x, not
!C                 the Bessel function J1(x) itself
!C       ============================================================
!C
!C       See S.Zhang, J.Jin, Computation of Special Functions, 
!C                           Wiley, 1996, Sec.5.2. 
!C       This subroutine is a modification of the subroutine JY01B
!C       which is described in that section.
        IMPLICIT NONE

        REAL*8 PI,X,T,T2,BJ0,BJ1,A0,P0,Q0,P1,Q1,TA0,TA1
        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           BJ0=1.0D0
           BJ1=0.5D0 
           RETURN
        ELSE IF (X.LE.4.0D0) THEN
           T=X/4.0D0
           T2=T*T
           BJ0=((((((-.5014415D-3*T2+.76771853D-2)*T2-.0709253492D0)*T2+.4443584263D0)*T2-1.7777560599D0)*T2+3.9999973021D0)*T2-3.9999998721D0)*T2+1.0D0
           BJ1=T*(((((((-.1289769D-3*T2+.22069155D-2)*T2-.0236616773D0)*T2+.1777582922D0)*T2-.8888839649D0)*T2+2.6666660544D0)*T2-3.9999999710D0)*T2+1.9999999998D0)
        ELSE
           T=4.0D0/X
           T2=T*T
           A0=DSQRT(2.0D0/(PI*X))
           P0=((((-.9285D-5*T2+.43506D-4)*T2-.122226D-3)*T2+.434725D-3)*T2-.4394275D-2)*T2+.999999997D0
           Q0=T*(((((.8099D-5*T2-.35614D-4)*T2+.85844D-4)*T2-.218024D-3)*T2+.1144106D-2)*T2-.031249995D0)
           TA0=X-.25D0*PI
           BJ0=A0*(P0*DCOS(TA0)-Q0*DSIN(TA0))
           P1=((((.10632D-4*T2-.50363D-4)*T2+.145575D-3)*T2-.559487D-3)*T2+.7323931D-2)*T2+1.000000004D0
           Q1=T*(((((-.9173D-5*T2+.40658D-4)*T2-.99941D-4)*T2+.266891D-3)*T2-.1601836D-2)*T2+.093749994D0)
           TA1=X-.75D0*PI
           BJ1=A0*(P1*DCOS(TA1)-Q1*DSIN(TA1))
        ENDIF
        BJ1 = BJ1/X

        RETURN
        END SUBROUTINE