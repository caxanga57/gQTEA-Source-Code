SUBROUTINE SR024
IMPLICIT NONE
! Numerov method to solve one dimensional shrodinger equation
! for hormonic oscillatior
! IRA N. LEVINE, QUANTUM CHEMISTRY, seventh edition, pg 74

DOUBLE PRECISION :: mXr,Ea,Gr1,Gr2,Gr3,Phi1,Phi2,Phi3,Xr1,Xr2,Xr3, &
					step,tmp1,tmp2,ss
INTEGER :: i,Nstep

OPEN(90,FILE='Phi.dat',ACTION='WRITE',STATUS='REPLACE')


WRITE(*,*) 'Enter minimum Xr,Ea,step,Nstep in this order'
READ(*,*) mXr,Ea,step,Nstep

ss = step*step
Xr1 = mXr
Xr2 = mXr+step
Xr3 = mXr+2*step

Phi1 = 0.0
Phi2 = 0.0001

DO i=1,Nstep
	!Equation 4.72
	Gr1 = Xr1*Xr1 - 2*Ea
	Gr2 = Xr2*Xr2 - 2*Ea
	Gr3 = Xr3*Xr3 - 2*Ea
	!Equation 4.67
	tmp1 = 2*Phi2 - Phi1 + 5.0*Gr2*Phi2*ss/6.0 + Gr1*Phi1*ss/12.0
	tmp2 = 1.0 - Gr3*ss/12.0

	Phi3 = tmp1/tmp2

	WRITE(90,*) Xr3,Phi3
	
	Xr1 = Xr2
	Xr2 = Xr1+step
	Xr3 = Xr1+2.0*step
	
	Phi1 = Phi2
	Phi2 = Phi3
	   
END DO

CLOSE(90)



END SUBROUTINE SR024