    SUBROUTINE SR022()
    IMPLICIT NONE
    !Compute the VELOCITY AUTOCORRELATION FUNCTION - VAF
	!and diffusion coefficient using VAF

	CHARACTER(LEN=8) :: option                     
                                                                      	
	WRITE(*,*) 'SELECT ATOMS FOR VAF? [yes|no]'
900	READ(*,*) option
	SELECT CASE(option)
		CASE('yes')
			CALL SR022A()
		CASE('no')
			CALL SR022B()
		CASE DEFAULT
			WRITE(*,*)
			WRITE(*,*) 'YOU MUST TYPE yes or no'
		GO TO 900
	END SELECT
	
    END SUBROUTINE SR022
