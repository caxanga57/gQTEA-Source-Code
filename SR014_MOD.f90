	MODULE SR014_MOD
	IMPLICIT NONE
	!Subroutine used to select frames from vmd.xyz file
	!and it generate gaussian input file
	!EXISTE : Logical variable
	!StartFrame : Start frame to start to extract from vmd.xyz file
	!EndFrame : Last frame to be extracted from vmd.xyz file
	!Sampling : Space between two sampling
	!I,J,M,N,P,Q,K : Counters
	!Natoms : Number of atoms
	!NFrames: Number of frames in the trajectory
	!Mult : Multiplicity, Mult=2*S+1
	!NSatoms: Number of selected atoms
	!chk : Checkpoint file name
	!Mem : Memory for gaussian allocation
	!nproc : Number of process to be used in gaussian calculation
	!Route : Section route from gaussian input file
	!Atoms : Vector to store element symbols
	!Option: Char variable for option
	!Charge: Charge on molecular system
	!Cxyz  : Matrix with atoms coordinates
	!Satoms: Vector with selected atoms indexes
	
	
	CONTAINS
	SUBROUTINE INPUTG09()
	
	LOGICAL :: EXISTE
	INTEGER :: StartFrame,Sampling,EndFrame,I,J,M,N,P,Q,K, &
	           Natoms,err,IO,NFrames,Mult,NSatoms
    CHARACTER(LEN=3) :: Atom,option 
    CHARACTER(LEN=2) :: charge
    CHARACTER(LEN=60) :: chk,Mem,nproc,Route
    REAL :: cx,cy,cz
	CHARACTER(LEN=3), ALLOCATABLE :: Atoms(:)
	REAL, ALLOCATABLE :: Cxyz(:,:)
	INTEGER, ALLOCATABLE :: Satoms(:) 
!Verify if exist vmd.xyz file	
	INQUIRE(FILE='vmd.xyz',EXIST=EXISTE)
	IF(EXISTE) THEN
	   OPEN(UNIT=90,FILE='vmd.xyz',ACTION='READ',STATUS='OLD')
	ELSE
	   WRITE(*,*)
	   WRITE(*,'(A)') 'SORRY! YOU HAVE TO PROVIDE THE vmd.xyz FILE'
	   STOP
	END IF
	
	WRITE(*,*)
	WRITE(*,'(A)') ' ANALYZING vmd.xyz FILE....PLEASE WAIT....THIS CAN TAKE A WHILE..!'
	WRITE(*,*)
	
	READ(90,*) Natoms
	NFrames=0	
	IO=0
	DO
	   IF(IO /= 0) THEN
	      EXIT
	   END IF	
	   DO I=1,Natoms
	      READ(90,*) Atom,cx,cy,cz
	   END DO
	   NFrames=NFrames + 1
	   READ(90,*,IOSTAT=IO) Natoms
	END DO
	WRITE(*,'(A,I7,A)') 'THERE ARE ',NFrames, ' FRAMES IN vmd.xyz FILE'
	WRITE(*,'(A,I7,A)') 'THERE ARE ',Natoms, ' ATOMS IN EACH FRAME'
	WRITE(*,*)
	WRITE(*,'(A)')'GIVE INITIAL FRAME'
	READ(*,*) StartFrame
	WRITE(*,*)
	WRITE(*,'(A)') 'GIVE Sampling FRAMES'
	READ(*,*) Sampling
	WRITE(*,*)
	WRITE(*,'(A)')'GIVE END FRAME'
	READ(*,*) EndFrame
	WRITE(*,*)
	WRITE(*,'(A)') 'ENTER THE CHK FILE NAME'
	READ(*,*) chk
	WRITE(*,'(A)') 'ENTER THE AMOUT OF Mem TO BE USED: e.g. 100MB'
	READ(*,*) Mem
	WRITE(*,'(A)') 'ENTER THE NUMBER OF PROCESSORS TO BE USED'
	READ(*,*) nproc
	WRITE(*,'(A)') 'ENTER THE SECTION ROUTE ENCLOSED BY DOUBLE QUOTES'
	READ(*,*) Route
	WRITE(*,'(A)') 'ENTER THE CHARGE AND Mult IN THIS ORDER'
	READ(*,*) charge,Mult
	WRITE(*,'(A)') 'WOULD YOU LIKE TO SETUP WIBERG BOND INDEX: yes or no'
	READ(*,*) option


	OPEN(UNIT=91,FILE='new-vmd.xyz',ACTION='READWRITE',STATUS='REPLACE')
	OPEN(UNIT=92,FILE='Ginput.gjf',ACTION='WRITE',STATUS='REPLACE')
	
	REWIND(UNIT=90)
	
	ALLOCATE(Atoms(Natoms),STAT=err)
	IF(err /= 0) PRINT*, 'Mem CAN NOT BE ALLOCATED FOR atoms ARRAY'

	ALLOCATE(Cxyz(Natoms,3),STAT=err)
	IF(err /= 0) PRINT*, 'MEMORY CAN NOT BE ALLOCATED FOR Cxyz ARRAY'

	READ(90,*) Natoms
	N = StartFrame
	DO I=1,NFrames			
		DO J = 1,Natoms
			READ(90,*) Atoms(J),(Cxyz(J,M),M = 1,3)
		END DO
		 
		IF(I == N) THEN
			WRITE(91,'(I7)') Natoms
			WRITE(91,*)
			
		    WRITE(92,'(A)') '--Link1--'
			WRITE(92,'(2A)') '%chk=',chk
			WRITE(92,'(2A)') '%Mem=',Mem
			WRITE(92,'(2A)') '%nprocshared=',nproc
			WRITE(92,'(A)')  Route
			WRITE(92,*)
			WRITE(92,'(A,I7)') 'GAUSSIAN CALCULATION ON FRAME ',n
			WRITE(92,*)
			WRITE(92,'(A,I2)') charge, Mult
			 
			DO P = 1,Natoms
				WRITE(91,'(A3,3F14.7)') atoms(p),(Cxyz(p,q),q = 1,3)
				WRITE(92,'(A3,3F14.7)') atoms(p),(Cxyz(p,q),q = 1,3)				
			END DO
			WRITE(92,*)
			IF(option == 'yes') THEN
			    WRITE(92,'(A)') '$nbo bndidx $end'
			    WRITE(92,*)
			END IF
			WRITE(*,'(A,I7)')'EXTRACTING FRAME ', n			
			N = N + Sampling
		END IF
		IF(I > EndFrame) THEN
		   EXIT
		END IF
		IF(I < NFrames) THEN
		   READ(90,*) Natoms
		END IF
	END DO
!=============================================================
	WRITE(*,'(A)') 'WOULD YOU LIKE TO SELECT ATOMS FROM SELECTED FRAMES(yes/no)'
	READ(*,*) option
	IF(option == 'yes') THEN
	    IO = 0
	    N = StartFrame
	    REWIND(UNIT=91)
	    OPEN(93,FILE='Molecule.gjf',ACTION='WRITE',STATUS='REPLACE')
	    OPEN(94,FILE='Solvent.gjf',ACTION='WRITE',STATUS='REPLACE')
	    
 	    READ(91,*) Natoms
	    WRITE(*,'(A)') 'HOW MANY ATOMS WOULD YOU LIKE TO SELECT'
	    READ(*,*) NSatoms
	    ALLOCATE(Satoms(NSatoms))
	    DO I=1,NSatoms
	        WRITE(*,'(A,I4)') 'TYPE THE INDEX OF THE ATOM SELECTED',I
	        READ(*,*) Satoms(I)
	    END DO
	    DO WHILE(IO == 0)
		        WRITE(93,'(A)') '--Link1--'
			    WRITE(93,'(2A)') '%chk=',chk
			    WRITE(93,'(2A)') '%Mem=',Mem
			    WRITE(93,'(2A)') '%nprocshared=',nproc
			    WRITE(93,'(A)')  Route
			    WRITE(93,*)
			    WRITE(93,'(A,I7)') 'GAUSSIAN CALCULATION ON FRAME ',N
			    WRITE(93,*)
			    WRITE(93,'(A,I2)') charge, Mult
		        WRITE(94,'(A)') '--Link1--'
			    WRITE(94,'(2A)') '%chk=',chk
			    WRITE(94,'(2A)') '%Mem=',Mem
			    WRITE(94,'(2A)') '%nprocshared=',nproc
			    WRITE(94,'(A)')  Route
			    WRITE(94,*)
			    WRITE(94,'(A,I7)') 'GAUSSIAN CALCULATION ON FRAME ',N
			    WRITE(94,*)
			    WRITE(94,'(A,I2)') charge, Mult
		        DO J = 1,Natoms
			        READ(91,*) atom,cx,cy,cz
			        DO K=1,NSatoms
			            IF(Satoms(K) == J) THEN
			                WRITE(93,'(A3,3F14.7)') atom,cx,cy,cz
			                GO TO 900
			            END IF
			        END DO
			        WRITE(94,'(A3,3F14.7)') Atom,cx,cy,cz
900			        CYCLE			        
		        END DO
                READ(91,*,IOSTAT=IO) Natoms
            	WRITE(93,*)
            	WRITE(94,*)
            	N = N + Sampling
	    END DO
	END IF
	
    
	CLOSE(90)
	CLOSE(91)
    CLOSE(92)
    CLOSE(93)
    CLOSE(94)
	END SUBROUTINE INPUTG09
	
!Subroutine to select frames

	SUBROUTINE SELECTFRAMES()
	
	LOGICAL :: EXISTE
	INTEGER :: StartFrame,Sampling,EndFrame,I,J,M,N,P, &
	           Natoms,err,IO,NFrames
    CHARACTER(LEN=3) :: Atom
    REAL :: cx,cy,cz
	CHARACTER(LEN=3), ALLOCATABLE :: Atoms(:)
	REAL, ALLOCATABLE :: Cxyz(:,:)

!Verify if exist vmd.xyz file	
	INQUIRE(FILE='vmd.xyz',EXIST=EXISTE)
	IF(EXISTE) THEN
	   OPEN(UNIT=90,FILE='vmd.xyz',ACTION='READ',STATUS='OLD')
	ELSE
	   WRITE(*,*)
	   WRITE(*,'(A)') 'SORRY! YOU HAVE TO PROVIDE THE vmd.xyz FILE'
	   STOP
	END IF
	
	WRITE(*,*)
	WRITE(*,'(A)') ' ANALYZING vmd.xyz FILE....PLEASE WAIT....THIS CAN TAKE A WHILE..!'
	WRITE(*,*)
	
	READ(90,*) Natoms
	NFrames=0	
	IO=0
	DO
	   IF(IO /= 0) THEN
	      EXIT
	   END IF	
	   DO I=1,Natoms
	      READ(90,*) Atom,cx,cy,cz
	   END DO
	   NFrames=NFrames + 1
	   READ(90,*,IOSTAT=IO) Natoms
	END DO
	WRITE(*,'(A,I7,A)') 'THERE ARE ',NFrames, ' FRAMES IN vmd.xyz FILE'
	WRITE(*,'(A,I7,A)') 'THERE ARE ',Natoms, ' ATOMS IN EACH FRAME'
	WRITE(*,*)
	WRITE(*,'(A)')'ENTER THE START FRAME TO BE SELECTED'
	READ(*,*) StartFrame
	WRITE(*,*)
	WRITE(*,'(A)') 'ENTER THE SAMPLING FRAMES'
	READ(*,*) Sampling
	WRITE(*,*)
	WRITE(*,'(A)')'ENTER THE FINAL FRAME TO BE SELECTED'
	READ(*,*) EndFrame
	WRITE(*,*)

	OPEN(UNIT=91,FILE='new-vmd.xyz',ACTION='READWRITE',STATUS='REPLACE')	
	REWIND(UNIT=90)	
	ALLOCATE(Atoms(Natoms),STAT=err)
	IF(err /= 0) PRINT*, 'MEMORY CAN NOT BE ALLOCATED FOR ATOMS ARRAY'

	ALLOCATE(Cxyz(Natoms,3),STAT=err)
	IF(err /= 0) PRINT*, 'MEMORY CAN NOT BE ALLOCATED FOR Cxyz ARRAY'

	READ(90,*) Natoms
	N = StartFrame
	DO I=1,NFrames			
		DO J = 1,Natoms
			READ(90,*) Atoms(J),(Cxyz(J,M),M = 1,3)
		END DO		 
		IF(I == N) THEN
			WRITE(91,'(I7)') Natoms
			WRITE(91,*)			 
			DO P = 1,Natoms
				WRITE(91,'(A3,3F14.7)') Atoms(P),(Cxyz(P,M),M = 1,3)
			END DO
			WRITE(*,'(A,I7)')'EXTRACTING FRAME ', N			
			N = N + Sampling
		END IF
		IF(I > EndFrame) THEN
		   EXIT
		END IF
		IF(I < NFrames) THEN
		   READ(90,*) Natoms
		END IF
	END DO
    
	CLOSE(90)
	CLOSE(91)
	
	END SUBROUTINE SELECTFRAMES
	
	END MODULE SR014_MOD
	
	
	
	