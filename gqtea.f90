	PROGRAM gqtea
	IMPLICIT NONE
! 	Version 0.4
! 	Implemented in October 01, 2011
! 	Last revision in September 02, 2013
! 	Last revision in September 04, 2013
!   Last revision in May 14, 2015 by Valter 
!   Last revision in June 23, 2016 by A. J. Camargo 

!   Subroutine SR001 : Selects one kind of trajectory to convert to vmd.xyz (VMD) file
! 	Subroutine SR002 : Converts trajectory (QE4) to vmd.xyz (VMD)
! 	Subroutine SR003 : Converts TRAJEC.xyz (CPMD) to vmd.xyz (VMD)
!   Subroutine SR004 : Selects bond lenght, angle or dihedral angle calculation
!   Subroutine SR005 : Computes bond lenght, bond lenght distribution
!                      function, and bond lenght free energy
!   Subroutine SR006 : Computes bond angle, bond angle distribution 
!                      function, and bond angle free energy
! 	Subroutine SR007 : Computes dihedral angle and dihedral andgle distribution function, etc
! 	Subroutine SR008 : Selects the Radial Distribution Function: 
!                      g(r) with PBC or g(r) without PBC
! 	Subroutine SR009 : Computes g(r) with PBC
! 	Subroutine SR010 : Computes g(r) without PBC
! 	Subroutine SR011 : Computes the difusion coefficient for a single 
!                      molecule using the molecular center of mass.
! 	Subroutine SR012 : Computes the mean residence time (MRT)
! 	Subroutine SR013 : Edit the *.evp file from cp.x to be opened on origin.
! 	Subroutine SR014 : Selects frames and generate gaussian, spectra, and SHTDDFT input files.
! 	Subroutine SR015 : Converts fractional coordinates to cartesian coordinates.
!   Subroutine SR016 : Creats cpmd input file using cartesian coordinate file, 
!                      i.e., in xyz format
! 	Subroutine SR017 : Converts the *.vel file from quantum espresso 4.3 into 
! 	                   cpmd dipole file format to be used as input 
!                      on trajec2atoms.x and fourier.x programs
! 	                   to calculate the power-spectrum.
!   Subroutine SR018 : Selects hydrogen bond lenght or hydrogen bond angle calculations
!   Subroutine SR019 : Computes hydrogen bond lenghts and its distribution function
!   Subroutine SR020 : Computes hydrogen bond angles and its distribution function
!   Subroutine SR021 : Computes d-TST rate constant
!	Subroutine SR022 : Selects the type of VAF (Velocity Autocorrelation Functionto) 
!                      to be computed
!	Subroutine SR022A: Computes the VAF with selected atoms from the system
!	Subroutine SR022B: Computes the VAF with all system's atoms
!	Subroutine SR023 : Creats a solvation box
!	Subroutine SR024 : Numerov method to solve monodimensional Schrodinger 
!                      equation for hamonic oscillator
!	Subroutine SR025 : Selects frames from TRAJECTORY file from CPMD and generates input files for
!                      SPECTRA and SHTDDFT dynamics calculations.
	
    INTEGER :: option
INTEGER, DIMENSION(8) :: values_in,values_out
REAL :: start,finish
WRITE(*,*)
WRITE(*,*) 'gqtea4 program: version 4'
WRITE(*,*) 'This program is brought you by gQTEA group from UEG univerty.'
WRITE(*,*)
CALL date_and_time(VALUES=values_in)
CALL cpu_time(start)
WRITE(*,*) 'AUTHORS:'
WRITE(*,*) 'A. J. CAMARGO      ajc@ueg.br'
WRITE(*,*) 'V. H. C. SILVA     fatioleg@ueg.br'
WRITE(*,*) 'S. S. OLIVEIRA     solemar@ueg.br'
WRITE(*,*) 'H. B. NAPOLITANO   hamilton@ueg.br'
	
  
100 WRITE(*,*)  
    WRITE(*,'(A)') 'Select an option below'
    WRITE(*,*)
    WRITE(*,'(A)') ' 1 - Creates vmd.xyz file used by gqtea'
    WRITE(*,'(A)') ' 2 - Bond lenght, angles, dihedral angles, and rover molecule speed'
    WRITE(*,'(A)') ' 3 - Radial distribution function'
    WRITE(*,'(A)') ' 4 - Diffusion coefficiente using the centre of molecular mass'
    WRITE(*,'(A)') ' 5 - Mean residence time (MRT)'
    WRITE(*,'(A)') ' 6 - Edit the file *.evp from cp.x'
    WRITE(*,'(A)') ' 7 - Selects frames and generates gaussian, spectra, and SHTDDFT input files'
    WRITE(*,'(A)') ' 8 - Converts fractional to cartesian coordinates'
    WRITE(*,'(A)') ' 9 - Creates cpmd input file'
    WRITE(*,'(A)') '10 - Converts *.vel (QE4) to cpmd dipole file format used by fourier program'
    WRITE(*,'(A)') '11 - Hydrogen bond analysis'
    WRITE(*,'(A)') '12 - Deformed transition state theory d-TST'
    WRITE(*,'(A)') '13 - Velocity autocorrelation function - VAF'
    WRITE(*,'(A)') '14 - Solvation box'
    WRITE(*,'(A)') '15 - Numerov method for harmonic oscillator: eductional purpose'
	

    READ(*,*) option
  
    SELECT CASE(option)
	    CASE(1)
	        CALL SR001()
	    CASE(2)
	        CALL SR004()
	    CASE(3)
	        CALL SR008()
	    CASE(4)
	    	CALL SR011()
	    CASE(5)
	    	CALL SR012()
	    CASE(6)
	    	CALL SR013()
	    CASE(7)
	    	CALL SR014()
	    CASE(8) 
	    	CALL SR015()
	    CASE(9)
	    	CALL SR027()
	    CASE(10)
	    	CALL SR017()
	    CASE(11)
	        CALL SR018()
	    CASE(12)
	        CALL SR021()
	    CASE(13)
	    	CALL SR022()
!		CASE(14)
!			CALL SR023()
	    CASE(15)
	    	CALL SR024()
	    CASE DEFAULT
	        WRITE(*,*)
	    	WRITE(*,*)'<<< THIS OPTION IS NOT IMPLEMENTED. TRY AGAIN. >>>'
	    	WRITE(*,*)
	    	GO TO 100
    END SELECT
    CALL cpu_time(finish)
    CALL date_and_time(VALUES=values_out)
    CALL DATE_TIME(values_in,values_out)
    WRITE(*,*)
    PRINT '("CPU TIME = ",F14.3," seconds.")',finish - start
    PRINT '("WALL CLOCK TIME ",I5," YEARS ",I5," MONTHS ",I5," DAYS ",I5," HOURS ",I5," MIN ",I5," SEC ", I5," MSEC")', &
    values_out(1)-values_in(1),values_out(2)-values_in(2),values_out(3)-values_in(3),&
    values_out(5)-values_in(5),values_out(6)-values_in(6),values_out(7)-values_in(7),values_out(8)-values_in(8)
    WRITE(*,*)
    WRITE(*,*)'PRESS ENTER TO EXIT'
    READ(*,*)

    END PROGRAM gqtea
