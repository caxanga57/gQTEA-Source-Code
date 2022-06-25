   SUBROUTINE SR004()
   IMPLICIT NONE
   !Choose a subroutine to calculate bond lenght, angle, 
   !dihedral angle, and rover molecule speed in collisional simulation
   
   INTEGER :: option
   
998 WRITE(*,*)
    WRITE(*,*) 'CALCULATION OPTIONS'
    WRITE(*,*) 
	WRITE(*,*) ' 1 - BOND LENGHT CALCULATION'
	WRITE(*,*) 
	WRITE(*,*) ' 2 - INERATOMIC ANGLE CALCULATION'
	WRITE(*,*)
	WRITE(*,*) ' 3 - DIHEDRAL ANGLE CALCULATION'
	WRITE(*,*)
	WRITE(*,*) ' 4 - ROVER MOLECULE SPEED IN COLLISIONAL SIMULATION IN CPMD'
  
	READ(*,*) option
   
	SELECT CASE(option)
		CASE(1)
			CALL SR005()
		CASE(2)
			CALL SR006()
		CASE(3)
			CALL SR007()
		CASE(4)
			CALL SR026()
		CASE DEFAULT
			WRITE(*,*)
			WRITE(*,*) 'THIS IS NOT AN VALID OPTION'
			GO TO 998
	END SELECT  
   
  END SUBROUTINE SR004
   
   