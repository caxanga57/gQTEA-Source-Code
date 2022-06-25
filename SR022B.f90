    SUBROUTINE SR022B()
    IMPLICIT NONE
    !Compute the VELOCITY AUTOCORRELATION FUNCTION (VAF)
	!and diffusion coefficient using VAF for ALL ATOMS in the system
	
	!EXISTE : Logical variable
	!Nsteps : Counter for the MOLECULAR DYNAMICS steps
	!Natoms : Counter for the number of atoms
	!I,J,K  : Counters
	!NumVAF : Counter for the number of VAF to be calculated for average
	!Nframes: Counter of frames
	!StartFrame: Initial frame for VAF calculation
	!ERRV : Status variable
	!Cv : Autocorrelation vector
	!T0 : Matrix w/ velocity components at time origin
	!V  : Matrix w/ velocity components
	!D  : Diffusion coefficient
	!aux1,aux2,aux3,AUX,dammy,nfi: Auxiliary variables
	
	LOGICAL :: EXISTE
	INTEGER :: nfi,Natoms,dammy,I,J,K,L,Nsteps,NumVAF,IO,&        
				Nframes,StartFrame,ERRV                     
	REAL :: aux1,aux2,aux3,AUX,D  
	REAL, ALLOCATABLE :: Cv(:),T0(:,:,:,:),V(:,:,:,:)   
                                                                      			
	INQUIRE(FILE='TRAJECTORY',EXIST=EXISTE)
	IF(EXISTE) THEN
	   OPEN(90,FILE='TRAJECTORY',ACTION='READ',STATUS='OLD')
	ELSE
	   WRITE(*,*)
	   WRITE(*,'(A)') 'SORRY! YOU HAVE TO PROVIDE THE TRAJECTORY FILE FROM CPMD'
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
!Loop to count the number of atoms in each frame
	DO WHILE( nfi == dammy)
		Natoms = Natoms+1
		READ(90,*) nfi
	END DO
	WRITE(*,*) 'ENTER THE START FRAME'	
	READ(*,*) StartFrame
	WRITE(*,*) ' NUMBER OF ATOMS IN THE SYSTEM  : ', Natoms	
	WRITE(*,*) ' NUMBER OF FRAMES IN THE SYSTEM : ', Nframes
	WRITE(*,*) ' ENTER THE NUMBER OF STEPS TO BE USED IN THE VAF CALCULATION:'
	READ(*,*) Nsteps
	AUX=REAL((Nframes-StartFrame))/REAL(Nsteps)
	WRITE(*,*) 'YOU CAN USE UP TO ',CEILING(AUX)-1 ,' VAF CALCULATION FOR AVERAGE'
	WRITE(*,*) 'ENTER THE NUMBER OF VAF CALCULATION TO BE USED FOR AVERAGE?'
	READ(*,*) NumVAF 
	OPEN(91,FILE='VAF.dat',ACTION='READWRITE',STATUS='REPLACE')
	OPEN(92,FILE='Diffusion.txt',ACTION='READWRITE',STATUS='REPLACE')

	REWIND(90)
	ALLOCATE(T0(NumVAF,1,Natoms,3),Cv(Nsteps),V(NumVAF,Nsteps,Natoms,3),STAT=ERRV)
    IF (ERRV /= 0) THEN 
		WRITE(*,*) 'ERROR IN ALLOCATION'
		STOP 
	END IF
! Loop to initialize the autocorrelation vector Cv(I)
	DO I=1,Nsteps
		Cv(I) = 0.0
	END DO
!Loop to read trajectory from CPMD
	DO K=1,NumVAF
		DO I=1,Nsteps
			DO L=1,Natoms		
				READ(90,*) dammy,aux1,aux2,aux3,(V(K,I,L,J),J=1,3)
			END DO
		END DO
	END DO
	
!Loop to create T0 matrix w/ the velocity components at time origin	
	DO K=1,NumVAF
		DO L=1,Natoms
			DO J=1,3
			T0(K,1,L,J) = V(K,1,L,J)
			END DO
		END DO
	END DO
!Loop to calculate VAF,i.e. dot product of velocity 
	DO K=1,NumVAF
		DO I=1,Nsteps
			DO L=1,Natoms
				DO J=1,3				
				Cv(I)= Cv(I)+V(K,I,L,J)*T0(K,1,L,J)
				END DO
			END DO
		END DO
	END DO
!Scaling by Natoms and NumVAF
	AUX=NumVAF*Natoms
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
		WRITE(91,FMT='(I5,F12.5)') I,Cv(I)
	END DO
	WRITE(92,FMT='(A,F10.5)') 'D = ',D
	CLOSE(UNIT=90)		   
    CLOSE(UNIT=91)
    CLOSE(UNIT=92)
	
	END SUBROUTINE SR022B

