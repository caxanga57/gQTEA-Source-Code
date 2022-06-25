	SUBROUTINE SR023()
	USE SR023MOD
	IMPLICIT NONE

	INTEGER :: option
	WRITE(*,'(A)') '1 - Box solvation solute'
	WRITE(*,'(A)') '2 - Solvation a single solute molecule by a mixture of solvente'
	WRITE(*,'(A)') '3 - Mixture of solvente'
	WRITE(*,*) 
100	READ(*,*) option
	SELECT CASE(option)
		CASE(1)
			CALL BoxSS()
		CASE(2)
			CALL SELECTFRAMES()
		CASE DEFAULT
			WRITE(*,*) 'THIS OPTION IS NOT IMPLEMENTED YET'			
			WRITE(*,*) 'CHOSE 1 OR 2'
			GO TO 100
	END SELECT
	
	
	END SUBROUTINE SR023