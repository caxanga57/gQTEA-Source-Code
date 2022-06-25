   SUBROUTINE SR020()
   IMPLICIT NONE
!  This subroutine is used to compute the hydrogen bond angles and its distribution function
	
	REAL, PARAMETER :: PI = 3.1415926,atufs=0.02418884326505, &
	                   ps=1000.0, R=0.001987204
	LOGICAL :: EXISTE
	INTEGER :: idx1,idx2,idx3,nframes, i, &
			   IO, bin, maxbin,j,iprint
	REAL :: soma,stdhbondangle, averagehbondangle, varhbondangle, &
			normaa,normab,innerprod,cosalfa,alfarad,alfadegree,deltatheta,centerbin,&
			dt,timefs,dammy,largesthbondangle,lowesthbondangle,ax,ay,az,bx,by,bz,&
			innerr,upperr,sumHISTbin,T,FreeEnergy
	CHARACTER(LEN=14) :: atom1,atom2,atom3,changle
	REAL,DIMENSION(:), ALLOCATABLE :: bondangle
	INTEGER,DIMENSION(:),ALLOCATABLE :: HIST
	REAL, DIMENSION(3,3) :: cxyz
	CHARACTER(LEN=6), DIMENSION(3) :: atoms
		
	WRITE(*,*)
	WRITE(*,'(A)') 'GIVE DELTA THETA TO CALCULATE BOND ANGLE DISTRIBUTION FUNCTION'
	READ(*,*) deltatheta
	WRITE(*,'(A)') 'GIVE dt (THIS IS THE TIME STEP IN ATOMIC UNIT USED IN QE INPUT FILE)'
	READ(*,*) dt
	WRITE(*,'(A)') 'GIVE THE SAMPLING INTERVAL (THIS IS THE iprint USED IN QE INPUT FILE)'
	READ(*,*) iprint
	WRITE(*,'(A)') 'GIVE INNER RADIUS OF THE SHELL'	
	READ(*,*) innerr
	WRITE(*,'(A)') 'GIVE OUTER RADIUS OF THE SHELL'	
	READ(*,*) upperr
	WRITE(*,'(A)') 'GIVE THE SIMULATION TEMPERATURE'	
	READ(*,*) T

	INQUIRE(FILE='tmphbond.dat',EXIST=EXISTE)
	IF(EXISTE) THEN
	   OPEN(90,FILE='tmphbond.dat',ACTION='READ',STATUS='OLD')
	ELSE
	   WRITE(*,*)
	   WRITE(*,'(A)') 'SORRY! YOU HAVE TO PROVIDE THE tmphbond.dat FILE'
	   STOP
	END IF		
	OPEN(91,FILE='hbondangle.dat',ACTION='READWRITE',STATUS='REPLACE')
	OPEN(92,FILE='hbondangledfc.dat',ACTION='WRITE',STATUS='REPLACE')
	OPEN(93,FILE='hbondangle.out',ACTION='WRITE',STATUS='REPLACE')
	OPEN(94,FILE='hangleFreeEnergy.dat',ACTION='WRITE',STATUS='REPLACE')			
			
   DO i=1,3
		READ(90,*) atoms(i), (cxyz(i,j),j=1,3)
	END DO
	 		
	IO = 0
	nframes = 0
	DO WHILE( IO == 0)
      
      ax = cxyz(1,1)-cxyz(2,1)
      ay = cxyz(1,2)-cxyz(2,2)
      az = cxyz(1,3)-cxyz(2,3)
      
      bx = cxyz(3,1)-cxyz(2,1)
      by = cxyz(3,2)-cxyz(2,2)
      bz = cxyz(3,3)-cxyz(2,3)
      			
		innerprod = ax*bx+ay*by+az*bz
		normaa = SQRT(ax*ax+ay*ay+az*az)
		normab = SQRT(bx*bx+by*by+bz*bz)
		IF(normaa > innerr .AND. normaa < upperr) THEN
		   cosalfa   = innerprod/(normaa*normab)
		   alfarad = ACOS(cosalfa)
		   alfadegree = alfarad*180.0/PI
		   timefs = REAL((nframes+1))*dt*REAL(iprint)*atufs
		   WRITE(91,'(F14.7,F14.7)') timefs, alfadegree
		   nframes = nframes+1
		END IF
	   DO i=1,3
			READ(90,*,IOSTAT=IO) atoms(i), (cxyz(i,j),j=1,3)
		END DO
	END DO
	atom1=atoms(1)
	atom2=atoms(2)
	atom3=atoms(3)
		
	REWIND(UNIT=91)
	
	maxbin = INT(200.0/deltatheta)+1
		
	ALLOCATE(HIST(maxbin))
	ALLOCATE(bondangle(nframes))
		
	DO bin=1,maxbin
		HIST(bin) = 0.0
	END DO
		
	DO i=1,nframes
		READ(91,*) dammy, bondangle(i)
		bin = INT(bondangle(i)/deltatheta)
		IF(bin <= maxbin) THEN
 			HIST(bin) = HIST(bin)+1
 		END IF
	END DO
	largesthbondangle = MAXVAL(bondangle)
	lowesthbondangle = MINVAL(bondangle)
	soma = 0
	DO i=1,nframes
		soma = soma + bondangle(i)
	END DO
	averagehbondangle = soma/REAL(nframes)
		
	DO i=1,nframes
		bondangle(i) = bondangle(i) - averagehbondangle
	END DO
		
	DO i=1,nframes
		bondangle(i) = bondangle(i)*bondangle(i)
	END DO
	sumHISTbin=SUM(HIST)	
	centerbin = deltatheta/2
	DO bin = 1,maxbin
		WRITE(92,'(F6.2,I7,F14.7)') centerbin, HIST(bin),(REAL(HIST(bin))/sumHISTbin)*100
		IF(HIST(bin) /= 0) THEN
		    FreeEnergy = -R*T*log(HIST(bin)/sumHISTbin) !Free energy using probability distribution function
			WRITE(94,'(F7.3,F14.4)') centerbin, FreeEnergy     		
		END IF
		centerbin = centerbin + deltatheta
	END DO
		
	varhbondangle = SUM(bondangle)/REAL(nframes-1)
	stdhbondangle = SQRT(SUM(bondangle)/REAL(nframes-1))
	
   changle=TRIM(atom1)//'-'//TRIM(atom2)//'-'//TRIM(atom3)
   
   WRITE(93,'(A)') 'tmphbond.dat (input, must be given)'
   WRITE(93,'(A)') 'hbondangle.dat (output,column 1: simulation time in fs; column 2: bondangle(degrees)'
   WRITE(93,'(A)') 'hbonddfc (output, column 1: atomic angles in degrees,'
   WRITE(93,'(A)') 'column 2: distribution function for the hydrogen bond angles'
   WRITE(93,*) 
   WRITE(93,'(A)')       '              SUMMARY OF THE CALCULATION'
   WRITE(93,*) 
	WRITE(93,'(A,A14)')   '       SELECTED ANGLE:    ',changle
	WRITE(93,'(A,I14)')   'NUMBER OF FRAMES USED: ',nframes
	WRITE(93,'(A,F14.7)') '    THE LARGEST ANGLE: ',largesthbondangle
	WRITE(93,'(A,F14.7)') '     THE LOWEST ANGLE: ',lowesthbondangle	
	WRITE(93,'(A,F14.7)') '        AVERAGE ANGLE: ',averagehbondangle
	WRITE(93,'(A,F14.7)') '             VARIANCE: ',varhbondangle
	WRITE(93,'(A,F14.7)') '   STANDARD DEVIATION: ',stdhbondangle	

	DEALLOCATE(HIST)
	DEALLOCATE(bondangle)

   CLOSE(UNIT=90)		   
   CLOSE(UNIT=91)
   CLOSE(UNIT=92)
   CLOSE(UNIT=93)
   
   END SUBROUTINE SR020
