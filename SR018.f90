   SUBROUTINE SR018()   
   IMPLICIT NONE
   
!    Input is the vmd.xyz file
!    Outputs: hbonds.dat (it contains the distance of the hydrogen bonds
!    between selected atoms).
!    hangles.dat (it contains the angles between selected atoms involved 
!    in the hbond)
   
   INTEGER :: option
   
996 WRITE(*,*)
    WRITE(*,*)
    WRITE(*,'(A)')'CHOOSE AN OPTION'
    WRITE(*,*)
    WRITE(*,'(A)') '  1 - HYDROGEN BOND LENGHT CALCULATION (caution: it needs vmd.xyz file)'
    WRITE(*,*)
    WRITE(*,'(A)') '  2 - HYDORGEN BOND ANGLE CALCULATION (caution: it needs tmphbond.dat file)'
    WRITE(*,*)
    READ(*,*) option 
   
    SELECT CASE(option)
      CASE(1)
         CALL SR019()
      CASE(2)
         CALL SR020()
      CASE DEFAULT
         WRITE(*,*)
         WRITE(*,'(A)') 'THIS IS NOT A VALID OPTION'
         GO TO 996
    END SELECT
   
   
   END SUBROUTINE SR018