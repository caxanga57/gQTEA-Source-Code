SUBROUTINE SR017()
IMPLICIT NONE
	!This soubroutine convert the velocity quantum espresso 4.3 file into 
	!cpmd dipole file format to be used as input for trajec2atoms.x and fourier.x programs
	!to calculate the power-spectrum.

	INTEGER :: natoms,steppos,stepvel,step=1,&
	           IOS1,IOS2,sampling, i
	CHARACTER(LEN=40) :: infile_pos, infile_vel
	DOUBLE PRECISION :: atomx,atomy,atomz,velx,vely,velz
	
	WRITE(*,*) 'WRITE THE ATOMIC POSTION FILE NAME FROM QE 4.3'
	WRITE(*,*) 
	READ(*,*) infile_pos
	WRITE(*,*)
	WRITE(*,*) 'WRITE THE VELOCITIES FILE NAME FROM QE 4.3'
	WRITE(*,*)
	READ(*,*) infile_vel
	WRITE(*,*) 
	WRITE(*,*) 'WRITE THE TOTAL NUMBER OF ATOMS IN THE CELL OR SUPERCELL'
	WRITE(*,*)
	READ(*,*) natoms
	WRITE(*,*) 
	WRITE(*,*) 'GIVE THE SAMPLING USED TO SAVE TRAJECTORY'
	WRITE(*,*)
	READ(*,*) sampling
		
	OPEN(UNIT=91,FILE=infile_pos,STATUS='OLD')
	OPEN(UNIT=92,FILE=infile_vel,STATUS='OLD')
	OPEN(UNIT=93,FILE='TRAJECTORY',STATUS='REPLACE')
	
	READ(91,*) steppos
	READ(92,*) stepvel
	IOS1 = 0
	IOS2 = 0
	DO WHILE(IOS1 == 0 .OR. IOS2 == 0)
		DO i=1,natoms
			READ(91,*) atomx,atomy,atomz
			READ(92,*) velx,vely,velz
			WRITE(93,'(1X,I6,6F24.14)') step, atomx,atomy,atomz,velx,vely,velz
	 	END DO
		step = step + sampling
		READ(91,*,IOSTAT=IOS1) steppos
		READ(92,*,IOSTAT=IOS2) stepvel
	END DO
END SUBROUTINE SR017