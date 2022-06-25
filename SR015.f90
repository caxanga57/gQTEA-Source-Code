	SUBROUTINE SR015()
	IMPLICIT NONE
! 	Ref. http://en.wikipedia.org/wiki/Fractional_coordinates	
! 	Subroutine used to convert fractional coordinates to cartesian coordinates.
! 	alfa, beta, gamma are the angles in degrees of the unit cell.
! 	alfarad,betarad,gammarad are the angles in radians of the unit cell.
! 	ahat,bhat, and chat are the coordinates in cartesian coordinates
! 	xhat, yhat, and zhat are the coordinates in fractional coordinates.
! 	a, b, and c are the edges of the unit cell.
! 	nat is the number of atoms in the unit cell.
! 	fractional.dat is an input file: 
! 	first line: nat (total number of atoms)
! 	second line: alfa, beta, and gamma in degrees
! 	third line: a, b, and c edges of the unit cell.
! 	The remaing lines contain the atomic symbols and the fractional coordinates.
! 	exemple of the fractional.dat input:
! 	 117
! 	 88.311 116.213 91.565
! 	 17.619 17.805 7.374
! 	
! 	Si 0.1806 0.1782 0.7480
! 	Si 0.6735 0.6738 0.7583
! 	Si 0.8263 0.1745 0.5610
! 	Si 0.3208 0.6769 0.5705
! 	Si 0.8248 0.8337 0.5587


	REAL, PARAMETER :: PI=3.1415926
	INTEGER :: nat, i
	REAL    :: alfa,beta,gamma,alfarad,betarad,gammarad,a,b,c, &
			   ahat,bhat,chat,xhat,yhat,zhat,v,vsq
	CHARACTER(LEN=2) atmsymbl
	
	
	
	OPEN(10,FILE='fractional.dat',ERR=100,STATUS='OLD')
	OPEN(11,FILE='cartesian.xyz',ERR=101,STATUS='REPLACE')
	
	READ(10,*) nat
	READ(10,*) alfa,beta,gamma
		alfarad  = (alfa*PI)/180.0
		betarad  = (beta*PI)/180.0
		gammarad = (gamma*PI)/180.0
	READ(10,*) a,b,c
	
	vsq = 1.0 - cos(alfarad)*cos(alfarad)- &
	        cos(betarad)*cos(betarad) - &
	        cos(gammarad)*cos(gammarad)+ &
	        2*cos(alfarad)*cos(betarad)*cos(gammarad)
	v = SQRT(vsq)
	
	WRITE(11,'(1X,I4)') nat
	WRITE(11,*)
	
	DO i=1,nat
		READ(10,*) atmsymbl,xhat,yhat,zhat 
		ahat = a*xhat + b*cos(gammarad)*yhat + c*cos(betarad)*zhat
		bhat = b*sin(gammarad)*yhat+c*(cos(alfarad)-cos(betarad)*cos(gammarad))*zhat/sin(gammarad)
		chat = c*v*zhat/sin(gammarad)
		WRITE(11,'(1X,A,3F10.4)') atmsymbl, ahat,bhat,chat
	END DO
	CLOSE(10)
	CLOSE(11)
100 WRITE(*,*) 'Impossible to open the file fractional.dat'
101 WRITE(*,*) 'Impossible to open the file cartesian.xyz'

	END SUBROUTINE SR015
	
	