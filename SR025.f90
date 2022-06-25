    SUBROUTINE SR025()
    IMPLICIT NONE
    !Subroutine to select frames from TRAJECTORY file from CPMD program
	!and generate input files for SPECTRA and SHTDDFT dynamics calculations
	
	!EXISTE : Logical variable
	!Nsteps : Counter for the MOLECULAR DYNAMICS steps
	!Natoms : Counter for the number of atoms
	!I,J,K,L  : Counters
	!Nframes: Counter for the number of frames
	!StartFrame: Initial frame to be slected
	!ERRV : Status variable
	!dammy,nfi: Auxiliary variables
	!NumSFrame: number of slected frames to be extracted
	!PTABLE: vector containing the elements of the periodic table
	!ELMTS: total number of elements in the vector PTABLE
	!NELMTS: vector containing the number of each element of the PTABLE. 
	!        Exemple: for PTABLE=(/'H','C','O','N','Cl'/) and NELMT=(/3,0,4,5,1/) means that 
	!		there are 3 H, 0 C, 4 O, 5 N, 1 Cl.
	!NstatesSHTDDFT: number of excited state to be used in the SHTDDFT dynamics
	!ForceSTATE: excited state to start the SHTDDFT dynamics
	
	LOGICAL :: EXISTE
	INTEGER :: nfi,Natoms,dammy,I,J,K,L,Nsteps,IO,&        
				Nframes,StartFrame,ERRV,NumSFrame,Iframe,&
				charge,ELMTS=7,Nstates,Cutoff,DUAL,NstatesSHTDDFT,&
				ForceSTATE
	REAL :: LatA, LatB, LatC, CosA, CosB, CosC
	CHARACTER(LEN=3) :: ext,symb
	CHARACTER(LEN=2),DIMENSION(7) :: PTABLE=(/'H ','C ','O ','N ','F ','Cl','Br'/)
	CHARACTER(LEN=1),DIMENSION(7) :: ANG=(/'S','P','P','P','P','P','D'/)
	INTEGER, DIMENSION(7) ::NELMTS=(/0,0,0,0,0,0,0/)
	CHARACTER(LEN=40) :: outfile1,outfile2,outfile3,outfile4,InputF
	REAL, ALLOCATABLE :: Cxyz(:,:),Vxyz(:,:)
	CHARACTER(LEN=3), ALLOCATABLE :: SYMBOLS(:)
	REAL, PARAMETER :: Bohr=0.52917720859
                                                                      			
	INQUIRE(FILE='TRAJECTORY',EXIST=EXISTE)
	IF(EXISTE) THEN
		OPEN(90,FILE='TRAJECTORY',ACTION='READ',STATUS='OLD')
	ELSE
		WRITE(*,*)
		WRITE(*,'(A)') 'THE TRAJECTORY FILE FROM CPMD WAS NOT FOUND IN THIS DIRECTORY'
		STOP
	END IF
	INQUIRE(FILE='GEOMETRY.xyz',EXIST=EXISTE)
	IF(EXISTE) THEN
	   OPEN(93,FILE='GEOMETRY.xyz',ACTION='READ',STATUS='OLD')
	ELSE
	   WRITE(*,*)
	   WRITE(*,'(A)') 'THE GEOMETRY.xyz FILE FROM CPMD WAS NOT FOUND IN THIS DIRECTORY'
	   STOP
	END IF
	WRITE(*,*) 'ENTER THE PREFIX FOR THE INPUT FILE NAMES: e.g. VitC'
	READ(*,*) InputF
	READ(90,*) nfi
	dammy=nfi
	IO = 0
	Nframes = 1
!Loop to get the number of frames
	DO WHILE(IO == 0)
		IF(nfi /= dammy) THEN
			dammy = nfi
			Nframes = Nframes + 1
		END IF
		READ(90,*,IOSTAT=IO) nfi
	END DO
	REWIND(90)
	READ(90,*) dammy
	nfi=dammy
	Natoms=0
!Loop to get the number of atoms in each frame
	DO WHILE( nfi == dammy)
		Natoms = Natoms+1
		READ(90,*) nfi
	END DO

	WRITE(*,*) ' NUMBER OF ATOMS IN THE SYSTEM  : ', Natoms	
	WRITE(*,*) ' NUMBER OF FRAMES IN THE SYSTEM : ', Nframes
		
	WRITE(*,*) 'ENTER THE START FRAME TO BE SELECTED FROM THE TRAJECTORY FILE'	
	READ(*,*) StartFrame
	WRITE(*,*) ' ENTER THE NUMBER OF FRAMES TO BE SKIPPED BETWEEN TWO SELECTED FRAMES'
	READ(*,*) Nsteps
	
	NumSFrame = INT((Nframes - StartFrame)/Nsteps)
	
	WRITE(*,*) ' IT WILL BE SELECTED ',NumSFrame,' FRAMES'
	REWIND(90)
	

	
	ALLOCATE(Cxyz(Natoms,3),Vxyz(Natoms,3),SYMBOLS(Natoms))

	READ(93,*) SYMBOLS(1)
	READ(93,*) SYMBOLS(2)
	DO I=1,Natoms
		READ(93,*) SYMBOLS(I)
	END DO
	REWIND(93)

	DO I=1,StartFrame
		DO J=1,Natoms
			READ(90,*) nfi
		END DO
	END DO
!>>>>>>>>>>>>>> GETTING THE NUMBER OF EACH DIFERENT ATOMS IN THE SYSTEM <<<<<<<<<<<<<<<<<<<<<
		DO I=1,ELMTS
			DO J=1,Natoms
				IF(PTABLE(I) == SYMBOLS(J)) THEN
					NELMTS(I) = NELMTS(I) + 1
				END IF
			END DO
			IF(NELMTS(I) /= 0) THEN
				WRITE(*,'(A,I3,A,A)') 'THERE ARE ',NELMTS(I),' ',PTABLE(I)
			END IF
		END DO
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	
	WRITE(*,*) 'ENTER THE NUMBER OF STATES TO BE CALCULATED IN THE ELECTRONIC SPECTRA'
	READ(*,*)   Nstates
	WRITE(*,*) 'ENTER THE NUMBER OF STATES TO BE USED IN THE SHTDDFT DYNAMICS'
	READ(*,*)   NstatesSHTDDFT
	WRITE(*,*) 'ENTER THE INITIAL STATE TO START THE SHTDDFT DYNAMICS'
	READ(*,*)   ForceSTATE
	WRITE(*,*) 'ENTER THE CHARGE ON THE SYSTEM'
	READ(*,*)   charge
	WRITE(*,*) 'ENTER THE CELL PARAMETER: LATTICE A, B,AND C; COSA,COSB,COSC IN THE SAME LINE'
	READ(*,*)   LatA, LatB, LatC, CosA, CosB, CosC
	WRITE(*,*) 'ENTER THE CUTOFF ENERGY IN RY'
	READ(*,*)   Cutoff
	WRITE(*,*) 'ENTER THE DUAL FOR CHARGE DENSITY EXPANSION IN PLANEWAVE'
	READ(*,*)   DUAL

! >>>>>>>>>>>>>>>> Initiating the extraction of selected frames and creating the input files. <<<<<<<<<<<<<<<<<<<<<<<<<	
	WRITE(*,*) 'Initiating the extraction of selected frames and creating the input files....'
	DO  Iframe= 1, NumSFrame

		IF(Iframe < 10 ) THEN
			WRITE(ext,'(''00'',I1)') Iframe
			outfile1 = 'GEOMETRY-'//ext//' '
			outfile2 = 'GEOMETRY-'//ext//'.xyz'
			outfile3 = TRIM(InputF)//'-SPECTRA-'//ext//'.inp'
			outfile4 = TRIM(InputF)//'-SHTDDFT-'//ext//'.inp'

		ELSE IF(Iframe > 99) THEN
			WRITE(ext,'(I3)') Iframe
			outfile1 = 'GEOMETRY-'//ext//' '
			outfile2 = 'GEOMETRY-'//ext//'.xyz'
			outfile3 = TRIM(InputF)//'-SPECTRA-'//ext//'.inp'
			outfile4 = TRIM(InputF)//'-SHTDDFT-'//ext//'.inp'

		ELSE
		    WRITE(ext,'(''0'',I2)') Iframe
		    outfile1 = 'GEOMETRY-'//ext//' '
			outfile2 = 'GEOMETRY-'//ext//'.xyz'
			outfile3 = TRIM(InputF)//'-SPECTRA-'//ext//'.inp'
			outfile4 = TRIM(InputF)//'-SHTDDFT-'//ext//'.inp'

		END IF
		
		OPEN(91,FILE=outfile1,ACTION='WRITE',STATUS='REPLACE')
		OPEN(92,FILE=outfile2,ACTION='WRITE',STATUS='REPLACE')
		OPEN(94,FILE=outfile3,ACTION='WRITE',STATUS='REPLACE')
		OPEN(95,FILE=outfile4,ACTION='WRITE',STATUS='REPLACE')
		
		DO J = 1, Natoms
			READ(90,*) nfi,(Cxyz(J,K),K=1,3), (Vxyz(J,L),L=1,3)
		END DO
		WRITE(92,'(I4)') Natoms
		WRITE(92,'(A,I5)') 'GEOMETRY ', Iframe
		DO J = 1, Natoms
			WRITE(91,*) (Cxyz(J,K),K=1,3), (Vxyz(J,L),L=1,3)
			WRITE(92,*) SYMBOLS(J),Bohr*Cxyz(J,1),Bohr*Cxyz(J,2),Bohr*Cxyz(J,3)
		END DO
		DO L=1,Nsteps
			DO J=1,Natoms
				READ(90,*,IOSTAT=IO) nfi
				IF(IO /= 0) THEN
					GOTO 999
				END IF
			END DO
		END DO
		
		
		! ---------------------- SPECTRA INPUT GENERATION ---------------------------
		WRITE(94,'(A)')   '&INFO'
		WRITE(95,'(A)')   '&INFO'
        WRITE(94,'(2XA)') 'CPMD INPUT FILE FOR SPECTRA CALCULATION GENERATED BY gqtea PROGRAM'
        WRITE(95,'(2XA)') 'CPMD INPUT FILE FOR SHTDDFT DYNAMICS GENERATED BY gqtea PROGRAM'
		WRITE(94,'(A)')   '&END'
		WRITE(95,'(A)')   '&END'
		WRITE(94,'(A)')
		WRITE(94,'(A)')
		WRITE(94,'(A)')   '&CPMD'
		WRITE(95,'(A)')   '&CPMD'
		WRITE(94,'(A)')
		WRITE(95,'(A)')
		WRITE(94,'(2XA)')  'ELECTRONIC SPECTRA'
		WRITE(95,'(2XA)')  'MOLECULAR DYNAMICS BO'
		WRITE(94,'(2XA)')  'DIAGONALIZATION LANCZOS'
		WRITE(95,'(2XA)')  'TDDFT'
		WRITE(94,'(2XA)')  'COMPRESS WRITE32' 
		WRITE(95,'(2XA)')  'RESTART  COORDINATES VELOCITIES GEOFILE  LINRES LATEST' 
		WRITE(94,'(A)')
		WRITE(95,'(A)')
		WRITE(94,'(2XA)')  'MEMORY BIG'
		WRITE(95,'(2XA)')  'MEMORY BIG'
		WRITE(94,'(A)')
		WRITE(95,'(2XA)')	   'STORE'
		WRITE(95,'(2XA)')    '  200'
		WRITE(95,'(2XA)') 'TIMESTEP'
		WRITE(95,'(2XA)') '  25'
		WRITE(95,'(2XA)') 'MAXSTEP'
		WRITE(95,'(2XA)') '  10000'
		WRITE(95,'(A)')
		WRITE(95,'(2XA)') 'TRAJECTORY XYZ'
		WRITE(95,'(A)')   ' 5'
		WRITE(95,'(2XA)') 'NOSE IONS'
		WRITE(95,'(2XA)') '  500.0 2000.0'		
		WRITE(95,'(A)')
		WRITE(94,'(A)')    '&END'
		WRITE(95,'(A)')    '&END'
		WRITE(94,'(A)')
		WRITE(95,'(A)')
		WRITE(94,'(A)')    '&TDDFT'
		WRITE(95,'(A)')    '&TDDFT'
		WRITE(94,'(A)')		
		WRITE(94,'(2XA)')  'STATES SINGLET'
		WRITE(95,'(2XA)')  'STATES SINGLET'
		WRITE(94,'(I5)') Nstates
		WRITE(95,'(I5)') NstatesSHTDDFT
		WRITE(94,'(A)')
		WRITE(95,'(2XA)')  'T-SHTDDFT'
		WRITE(95,'(A)')  
		WRITE(95,'(2XA)')    'FORCE STATE'
		WRITE(95,'(I5)')    ForceSTATE		
		WRITE(94,'(2XA)')    'TAMM-DANCOFF'
		WRITE(95,'(2XA)')    'TAMM-DANCOFF'
		WRITE(94,'(2XA)')    'DAVIDSON PARAMETER'
		WRITE(95,'(2XA)')    'DAVIDSON PARAMETER'
		WRITE(94,'(2XA)')    '  150 1.D-7 50'
		WRITE(95,'(2XA)')    '  150 1.D-7 50'
		WRITE(94,'(A)')
		WRITE(95,'(A)')
		WRITE(94,'(A)')  '&END'
		WRITE(95,'(A)')  '&END'
		WRITE(94,'(A)')
		WRITE(95,'(A)')
		WRITE(94,'(A)')    '&DFT'
		WRITE(95,'(A)')    '&DFT'
		WRITE(94,'(2XA)')  'NEWCODE'
		WRITE(95,'(2XA)')  'NEWCODE'
		WRITE(94,'(2XA)')  'FUNCTIONAL PBE'
		WRITE(95,'(2XA)')  'FUNCTIONAL PBE'
		WRITE(94,'(A)')    '&END'
		WRITE(95,'(A)')    '&END'
		WRITE(94,'(A)')
		WRITE(95,'(A)')
		WRITE(94,'(A)')    '&SYSTEM'
		WRITE(95,'(A)')    '&SYSTEM'
		WRITE(94,'(2XA)')  'CHARGE'
		WRITE(95,'(2XA)')  'CHARGE'
		WRITE(94,'(I5)')  charge
		WRITE(95,'(I5)')  charge
		WRITE(94,'(2XA)')  'SYMMETRY'
		WRITE(95,'(2XA)')  'SYMMETRY'
		WRITE(94,'(2XA)')  '  1'
		WRITE(95,'(2XA)')  '  1'
		WRITE(94,'(2XA)')   'ANGSTROM'
		WRITE(95,'(2XA)')   'ANGSTROM'
		WRITE(94,'(2XA)')   'CELL'
		WRITE(95,'(2XA)')   'CELL'
		WRITE(94,'(2X,6F7.2)') LatA,LatB/LatA,LatC/LatA, CosA,CosB, CosC
		WRITE(95,'(2X,6F7.2)') LatA,LatB/LatA,LatC/LatA, CosA,CosB, CosC
		WRITE(94,'(A)')
		WRITE(95,'(A)')
		WRITE(94,'(2XA)')    'CUTOFF'
		WRITE(95,'(2XA)')    'CUTOFF'
		WRITE(94,'(I5)')  cutoff
		WRITE(95,'(I5)')  cutoff
		WRITE(94,'(2XA)')    'DUAL'
		WRITE(95,'(2XA)')    'DUAL'
		WRITE(94,'(I4)')  DUAL
		WRITE(95,'(I4)')  DUAL
		WRITE(94,'(A)')
		WRITE(95,'(A)')
		WRITE(94,'(A)')    '&END'
		WRITE(95,'(A)')    '&END'
		WRITE(94,'(A)')
		WRITE(95,'(A)')
		WRITE(94,'(A)')    '&ATOMS'
		WRITE(95,'(A)')    '&ATOMS'
		
		DO I=1,ELMTS
			IF(NELMTS(I)/= 0) THEN
				WRITE(94,'(A)')'*'//TRIM(PTABLE(I))//'_MT_PBE.psp'
				WRITE(95,'(A)')'*'//TRIM(PTABLE(I))//'_MT_PBE.psp'
				WRITE(94,'(A,A)')'LMAX=',ANG(I)
				WRITE(94,'(I3)') NELMTS(I)
				WRITE(95,'(A,A)')'LMAX=',ANG(I)
				WRITE(95,'(I3)') NELMTS(I)

				DO J=1,Natoms
					IF(PTABLE(I) == SYMBOLS(J)) THEN
						WRITE(94,'(3F14.5)') Bohr*Cxyz(J,1),Bohr*Cxyz(J,2),Bohr*Cxyz(J,3)
						WRITE(95,'(3F14.5)') Bohr*Cxyz(J,1),Bohr*Cxyz(J,2),Bohr*Cxyz(J,3)
					END IF
				END DO
			ELSE 
				CYCLE
			END IF
			WRITE(94,*)
			WRITE(95,*)
		END DO
		WRITE(94,'(A)')'&END'
		WRITE(95,'(A)')'&END'

		CLOSE(UNIT=91)		   
		CLOSE(UNIT=92)
		CLOSE(UNIT=94)
		CLOSE(UNIT=95)
	END DO

999 CLOSE(UNIT=90)
	CLOSE(UNIT=93)
	
	END SUBROUTINE SR025

