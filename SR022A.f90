   SUBROUTINE SR022A()
    IMPLICIT NONE
	!10:47 AM Friday, June 24, 2016
    !Compute the VELOCITY AUTOCORRELATION FUNCTION  (VAF)
	!and diffusion coefficient for a set of SELECTED ATOMS
	
	!EXISTE : Logical variable
	!Nsteps : Counter for the MOLECULAR DYNAMICS steps
	!Natoms : Counter for the number of atoms
	!I,J,K,M: Counters
	!NumVAF : Counter for the number of VAF to be calculated for average
	!Nframes: Counter for the frame number
	!StartFrame: Initial frame for VAF calculation
	!ERRV : Status variable
	!Cv : Autocorrelation vector
	!T0 : Matrix w/ velocity components at time origin
	!V  : Matrix w/ velocity components
	!D  : Diffusion coefficient
	!aux1,aux2,aux3,AUX,dammy,nfi: Auxiliary variables
	!SA : Vector of selected Atoms
	!NSA: Number of selected atoms
	!option : Character variable to select atoms
	
	LOGICAL :: EXISTE
	INTEGER :: nfi,Natoms,dammy,I,J,K,L,M,Nsteps,NumVAF,IO,&        
				Nframes,StartFrame,ERRV,NSA
	CHARACTER(LEN=8) :: option                     
	REAL :: aux1,aux2,aux3,AUX,D  

	INTEGER, ALLOCATABLE :: SA(:)
	REAL, ALLOCATABLE :: Cv(:),T0(:,:,:,:),V(:,:,:,:)   
                                                                      			
	INQUIRE(FILE='TRAJECTORY',EXIST=EXISTE)
	IF(EXISTE) THEN
	OPEN(90,FILE='TRAJECTORY',ACTION='READ',STATUS='OLD')
	ELSE
	   WRITE(*,*)
	   WRITE(*,'(A)') 'SORRY! YOU HAVE TO PROVIDE THE TRAJECTORY FILE FROM CPMD PROGRAM'
	   STOP
	END IF
	READ(90,*) nfi
	dammy=nfi
	IO = 0
	Nframes = 1
!Loop to count the number of frames
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
!Loop to count the number of atoms in the first frame
	DO WHILE( nfi == dammy)
		Natoms = Natoms+1
		READ(90,*) nfi
	END DO
	WRITE(*,'(A)') 'ENTER THE START FRAME'	
	READ(*,*) StartFrame
	WRITE(*,'(A,I6)') ' NUMBER OF ATOMS IN THE SYSTEM  : ', Natoms	
	WRITE(*,'(A,I6)') ' NUMBER OF FRAMES IN THE SYSTEM : ', Nframes
	WRITE(*,'(A)') ' ENTER THE NUMBER OF SEQUENTIAL FRAMES TO BE USED IN EACH VAF CALCULATION:'
	READ(*,*) Nsteps
	AUX=REAL((Nframes-StartFrame))/REAL(Nsteps)
	WRITE(*,'(A,I5,A)') 'YOU CAN USE UP TO ',CEILING(AUX)-1 ,' VAF CALCULATION FOR AVERAGE'
	WRITE(*,'(A)') 'ENTER THE NUMBER OF VAF CALCULATION TO BE USED FOR AVERAGE?'
	READ(*,*) NumVAF
	WRITE(*,'(A)') 'ENTER THE NUMBER OF ATOMS TO BE SELECTED'
	READ(*,*) NSA
	
	OPEN(91,FILE='VAF.dat',ACTION='READWRITE',STATUS='REPLACE')
	OPEN(92,FILE='Diffusion.txt',ACTION='READWRITE',STATUS='REPLACE')

	REWIND(90)
	ALLOCATE(T0(NumVAF,1,Natoms,3),Cv(Nsteps),&
			V(NumVAF,Nsteps,Natoms,3),SA(NSA),STAT=ERRV)
    IF (ERRV /= 0) THEN 
		WRITE(*,'(A)') 'ERROR IN ALLOCATION'
		STOP 
	END IF
!Loop to read trajectory from CPMD
	WRITE(*,'(A)') ' READING TRAJECTORY FILE FROM CPMD ....'
	DO K=1,NumVAF
		DO I=1,Nsteps
			DO L=1,Natoms		
				READ(90,*) dammy,aux1,aux2,aux3,(V(K,I,L,J),J=1,3)
			END DO
		END DO
	END DO

! Loop to initialize Cv(I)
	DO I=1,Nsteps
		Cv(I) = 0.0
	END DO
!Loop to select atoms
	DO I=1,NSA
		WRITE(*,'(A,I5)') 'ENTER INDEX OF ATOM',I
		READ(*,*) SA(I)
	END DO

!Loop to create T0 matrix w/ the velocity components at time origin for each VAF	
	DO K=1,NumVAF
		DO L=1,Natoms
			DO J=1,3
			T0(K,1,L,J) = V(K,1,L,J)
			END DO
		END DO
	END DO
!Loop to calculate VAF for selected atoms 
	DO K=1,NumVAF
		DO I=1,Nsteps
			DO L=1,Natoms
				DO M=1,NSA				
					IF (SA(M) == L) THEN
						DO J=1,3				
							Cv(I)= Cv(I)+V(K,I,L,J)*T0(K,1,L,J)
						END DO
					END IF
				END DO
			END DO
		END DO
	END DO
!Scaling by Natoms and NumVAF
	AUX=NumVAF*NSA
	DO I=1,Nsteps
		Cv(I)=Cv(I)/AUX
	END DO
!Scaling by Cv(1)
	AUX=Cv(1)
	DO I=1,Nsteps
		Cv(I)=Cv(I)/AUX
	END DO
!Calculating diffusion coeficient D 
	D = SUM(Cv)/6.0
!Print results
	DO I=1,Nsteps
		WRITE(91,'(I6,F12.5)') I,Cv(I)
	END DO
	WRITE(92,FMT='(A,F10.5)') 'D = ',D
	
	DEALLOCATE(T0,Cv,V,SA)
	CLOSE(UNIT=90)		   
    CLOSE(UNIT=91)
    CLOSE(UNIT=92)
		
	END SUBROUTINE SR022A
	
