	SUBROUTINE SR014()
	USE SR014_MOD
	IMPLICIT NONE

	
	INTEGER :: option
100	WRITE(*,'(A)') '1 - Select frames from TRAJEC.xyz and generate gaussian job file'
	WRITE(*,'(A)') '2 - Only select frames from TRAJEC.xyz file'
	WRITE(*,'(A)') '3 - Select frames from TRAJECTOTY file from CPMD program and'
	WRITE(*,'(A)') '    generate input files for SPECTRA and SHTDDFT calculations'
	WRITE(*,*) 
	READ(*,*) option
	SELECT CASE(option)
		CASE(1)
			CALL INPUTG09()
		CASE(2)
			CALL SELECTFRAMES()
		CASE(3)
			CALL SR025()
		CASE DEFAULT
			WRITE(*,*) 'THIS OPTION IS NOT IMPLEMENTED'
			WRITE(*,*) 
			GO TO 100
	END SELECT
	END SUBROUTINE SR014