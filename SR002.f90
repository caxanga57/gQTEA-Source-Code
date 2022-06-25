   SUBROUTINE SR002()
   IMPLICIT NONE
    !Generate trajectory xyz for vmd from Quantum Espresso output
    
    	
    REAL, PARAMETER   :: bohr=0.52917720859
    LOGICAL :: existe	
    INTEGER  :: i,natoms, err, IO,trajoption,nstep      
    CHARACTER(LEN=40) :: infile1, infile2, dammy
    CHARACTER(LEN=6) :: atom
    CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: atoms     
    REAL :: cx, cy, cz
   
   
    WRITE(*,*) 'GIVE THE NUMBER OF atoms IN THE SYSTEM'
	READ(*,*) natoms
	WRITE(*,*) 'GIVE THE OUTPUT *.out FROM QE'
	READ(*,*) infile1
    WRITE(*,*) 'GIVE THE OUTPUT *.pos FROM QE'
    READ(*,*) infile2
      
	ALLOCATE(atoms(natoms),STAT=err)
	IF (err /= 0) PRINT*, 'MEMORY CAN NOT BE ALLOCATED FOR atoms ARRAY'
	
    OPEN(90,FILE= infile1,ACTION='READ',STATUS='OLD')
    OPEN(91,FILE= infile2,ACTION='READ',STATUS='OLD')
    OPEN(92,FILE='vmd.xyz',ACTION='WRITE',STATUS='REPLACE')
    IO = 0  
	DO WHILE(IO == 0)
	   READ(90,*,IOSTAT=IO) dammy
      IF (dammy == 'Scaled') THEN
		   DO i=1,natoms
			   READ(90,*) dammy
            atoms(i)= dammy
		   END DO
		   EXIT
		END IF
	END DO     
    IO=0
    READ(91,*) dammy
	 nstep = 1	
    DO WHILE(IO == 0)
      WRITE(*,'(A,I7)') 'WRITING FRAME ', nstep
      WRITE(92,FMT='(I6)') natoms
      WRITE(92,*) 
        
      DO i=1,natoms
        	READ(91,*) cx, cy, cz
         WRITE(92,'(A3,3F14.6)') atoms(i), cx*bohr,cy*bohr,cz*bohr
		END DO
		READ(91,*,IOSTAT=IO) dammy
        nstep = nstep + 1
    END DO
        
    CLOSE(UNIT=90)
    CLOSE(UNIT=91)
    CLOSE(UNIT=92)
    WRITE(*,*)
    WRITE(*,'(A)') 'TRAJECTORY WAS CREATED SUCCESSFULLY!'
    WRITE(*,*)
	
	END SUBROUTINE SR002