   SUBROUTINE SR008()
   IMPLICIT NONE
   !Subroutine to select one kind of RDF (Radial Distribution Function) to be calculated
   
   INTEGER :: option
   
100 WRITE(*,*)

	WRITE(*,*) ' 1 - RADIAL DISTRIBUTION FUNCTION g(r) WITH PBC'
	WRITE(*,*) 
	WRITE(*,*) ' 2 - RADIAL DISTRIBUTION FUNCTION g(r) WITHOUT PBC'
  
    READ(*,*) option
   
    SELECT CASE(option)
      CASE(1)
		CALL SR009()
      CASE(2)
		CALL SR010()
      CASE DEFAULT 
        WRITE(*,*) 
        WRITE(*,*) 'THIS IS NOT AN VALID OPTION'
        GO TO 100
    END SELECT  
   
  END SUBROUTINE SR008
   
   