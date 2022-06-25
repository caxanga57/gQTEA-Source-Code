   SUBROUTINE SR001()
   
   IMPLICIT NONE
   
   INTEGER :: option
   
999 WRITE(*,*)
    WRITE(*,*) 'CHOOSE AN OPTION'
    WRITE(*,*) 
    WRITE(*,*) ' 1 - TO CREATE vmd.xyz FILE FROM QUNTUM ESPRESSO OUTPUT FILES'
    WRITE(*,*)    
    WRITE(*,*) ' 2 - TO CREATE vmd.xyz FILE FROM TRAJEC.xyz OUTPUT FROM CPMD'
    READ(*,*) option
   
    SELECT CASE(option)
       CASE(1)
          CALL SR002()
       CASE(2)
          CALL SR003()
       CASE DEFAULT
          WRITE(*,*)
          WRITE(*,*) 'THIS IS NOT A VALID OPTION'
          GO TO 999
    END SELECT
    
  END SUBROUTINE SR001
