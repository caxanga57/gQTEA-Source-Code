   SUBROUTINE SR003()
   IMPLICIT NONE
   !Generate a trajectory xyz for vmd from cpmd program
   
   
   LOGICAL :: EXISTE
   INTEGER :: natoms,IO,i
   CHARACTER(LEN=10) :: dammy,atom
   REAL :: cx,cy,cz
   
	
	INQUIRE(FILE='TRAJEC.xyz',EXIST=EXISTE)
	IF(EXISTE) THEN
	   OPEN(90,FILE='TRAJEC.xyz',ACTION='READ',STATUS='OLD')
	ELSE
	WRITE(*,*)
	WRITE(*,'(A)') 'SORRY! YOU HAVE TO PROVIDE THE TRAJEC.xyz FILE FROM CMPD'
	WRITE(*,*)
	WRITE(*,*) 'TYPE RETURN TO EXIT'
	READ(*,*)
	STOP
	END IF
	OPEN(91,FILE='vmd.xyz',ACTION='WRITE',STATUS='REPLACE')		
    READ(90,*) natoms
    READ(90,*) dammy
    IO = 0
    DO
      WRITE(91,'(I6)') natoms
      WRITE(91,*)    
      DO i=1,natoms
         READ(90,*) atom,cx,cy,cz
         WRITE(91,'(A3,3F14.6)') atom,cx,cy,cz
      END DO
      READ(90,*,IOSTAT=IO) natoms
      IF(IO /= 0) THEN
         EXIT
      END IF
      READ(90,*) dammy         
   END DO
    
   END SUBROUTINE SR003

