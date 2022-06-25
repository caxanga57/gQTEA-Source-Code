SUBROUTINE SR026()
IMPLICIT NONE

!This subroutine calculates the speed of the rover molecule
!from TRAJECTORY file of the CPMD program

LOGICAL :: EXISTE
INTEGER :: atomIdx,totNumAtoms,frame,I,J,IO,oldFrame,newFrame,nFrame
REAL, ALLOCATABLE :: posXYZ(:,:),velXYZ(:,:) 
REAL :: innerProd,speedAU,speedMS
	
INQUIRE(FILE='TRAJECTORY',EXIST=EXISTE)
IF(EXISTE) THEN
	OPEN(10,FILE='TRAJECTORY',ACTION='READ',STATUS='OLD')
ELSE
	WRITE(*,*)
	WRITE(*,'(A)') 'The TRAJECTORY file is not in the same folder as the gqtea program.'
	WRITE(*,*)
	WRITE(*,'(A)') 'TYPE RETURN TO EXIT'
	READ(*,*)
	STOP
END IF

WRITE(*,'(A)') 'ENTER ATOM INDEX TO CALCULATE ITS SPEED'
READ(*,*) atomIdx		

OPEN(11,FILE='roverSpeed.dat',ACTION='WRITE',STATUS='REPLACE')	

READ(10,*) oldFrame
READ(10,*) newFrame
totNumAtoms = 1
DO WHILE(newFrame == oldFrame)
	totNumAtoms = totNumAtoms +1
	READ(10,*) newFrame
END DO
WRITE(*,*) 'THE SYSTEM CONTAIN ',totNumAtoms,' ATOMS'
REWIND(10)
IO = 0
frame=1
ALLOCATE(posXYZ(totNumAtoms,3),velXYZ(totNumAtoms,3))
READ(10,*) frame,(posXYZ(1,j),j=1,3),(velXYZ(1,j),j=1,3)
nFrame=1
DO WHILE(IO == 0)
	DO I=2,totNumAtoms
		READ(10,*) frame,(posXYZ(I,J),J=1,3),(velXYZ(I,J),J=1,3)
	END DO
	innerProd = velXYZ(atomIdx,1)*velXYZ(atomIdx,1)+&
	            velXYZ(atomIdx,2)*velXYZ(atomIdx,2)+&
				velXYZ(atomIdx,3)*velXYZ(atomIdx,3)
	speedAU = SQRT(innerProd)        !speed in atomic unit
	speedMS = speedAU*2187691.2541   !speed in meter/second
	WRITE(11,'(1X,I5,F14.7,ES14.7)') nFrame, speedAU, speedMS	
	READ(10,*,IOSTAT=IO) frame,(posXYZ(1,J),J=1,3),(velXYZ(1,J),J=1,3)
	nFrame = nFrame + 1
END DO

CLOSE(UNIT=10)
CLOSE(UNIT=11)
END SUBROUTINE SR026