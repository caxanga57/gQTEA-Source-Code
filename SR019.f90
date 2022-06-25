   SUBROUTINE SR019()
   IMPLICIT NONE
   
!    This subroutine calculates the hydrogen bond lenght, hydrogen bond distribution function
!    hydrogen bond free energy using as input the vmd.xyz file

!    vmd.xyz (input, must be given)
!    tmphbond.dat (output used during the calculation)
!    hbond.dat (output,column 1: simulation time in fs; column 2: bondlenght
!    between first and second atom indexes given; column 3: bondlenght
!    between second and third atom indices given).
!    hbonddfc (output, column 1: atomic distance in Angstrom,
!    column 2: distribution function for first and second atom indexes;
!    column 3: distribution function for second and third atom indexes).
!    vbond=R2(O--H) - R1(O-H); O--H is Hbond and O-H is OH bond
   
   REAL, PARAMETER :: atufs=0.02418884326505,R=0.001987204
	LOGICAL :: EXISTE 
   INTEGER :: natoms,IO,j,i,nframes,sampling,m,&
              maxbin,bin
   INTEGER,DIMENSION(3) :: idx
   CHARACTER(LEN=14) :: atom,chidx,out1,out2
   CHARACTER(LEN=7),DIMENSION(3) :: chatom
   REAL :: cxx,cyy,czz,maxr,deltar,dt,innerr,&
           upperr,deltax,deltay,deltaz,deltaxsq,&
           deltaysq,deltazsq,bondx1,bondx2,rtime,&
           soma1,soma2,averagehbondx1,averagehbondx2,&
           largesthbondx1,largesthbondx2,varhbondx1,&
           varhbondx2,stdhbondx1,stdhbondx2,centerbin,&
           lowesthbondx1,lowesthbondx2,T,FreeEnergy,HFreeEnergy, &
           sumvHISTminus,sumvHISTplus
           
   REAL,DIMENSION(2) :: sumHISTbin
   REAL,DIMENSION(3) :: cx,cy,cz
   REAL,DIMENSION(:,:),ALLOCATABLE :: bondx
   INTEGER,DIMENSION(:,:),ALLOCATABLE :: HIST
   REAL,DIMENSION(:), ALLOCATABLE :: vbond, tmp,vHISTplus,vHISTminus

	INQUIRE(FILE='vmd.xyz',EXIST=EXISTE)
	IF(EXISTE) THEN
	   OPEN(90,FILE='vmd.xyz',ACTION='READ',STATUS='OLD')
	ELSE
	   WRITE(*,*)
	   WRITE(*,'(A)') 'SORRY! YOU HAVE TO PROVIDE THE tmphbond.dat FILE'
	   STOP
	END IF		
   
   WRITE(*,'(A)') 'GIVE 3 ATOMIC INDEXES. CAUTION: THE ORDER IS IMPORTANT'
   WRITE(*,'(A)')
   WRITE(*,'(A)') 'e.g. if you have O1-H2---O3, type idx3 idx2 idx1 in this order.'
   WRITE(*,'(A)') 'If you have O1---H2-O3, type idx1 idx2 idx3 in this order.'
   WRITE(*,'(A)') 'O1---H2 stand for hydrogen bond between idx1 and idx2.'

   READ(*,*) idx(1),idx(2),idx(3)
	WRITE(*,'(A)') 'GIVE THE MAX BOND LENGHT TO CALCULATE THE BOND LENGHT DISTRIBUTION FUNCTION'
	READ(*,*) maxr
	WRITE(*,'(A)') 'GIVE THE WIDTH OF THE HISTOGRAM BINS IN ANGSTROM'
	READ(*,*) deltar
	WRITE(*,'(A)') 'GIVE dt (THIS IS THE TIME STEP IN ATOMIC UNIT USED IN QE INPUT FILE)'
	READ(*,*) dt
	WRITE(*,'(A)') 'GIVE THE SAMPLING INTERVAL (THIS IS THE iprint USED IN QE INPUT FILE)'
	READ(*,*) sampling   
	WRITE(*,'(A)') 'GIVE INNER RADIUS OF THE SHELL'	
	READ(*,*) innerr
	WRITE(*,'(A)') 'GIVE OUTER RADIUS OF THE SHELL'	
	READ(*,*) upperr
	WRITE(*,'(A)') 'GIVE THE SIMULATION TEMPERATURE IN KELVIN'
	READ(*,*) T


	OPEN(91,FILE='tmphbond.dat',ACTION='READWRITE',STATUS='REPLACE')
	OPEN(92,FILE='hbond.dat',ACTION='READWRITE',STATUS='REPLACE')
	OPEN(93,FILE='hbonddfc.dat',ACTION='WRITE',STATUS='REPLACE')
	OPEN(94,FILE='hbond.out',ACTION='WRITE',STATUS='REPLACE')
	OPEN(95,FILE='hbonddfc-Norm.dat',ACTION='WRITE',STATUS='REPLACE')
	OPEN(96,FILE='hbondFreeEnergy.dat',ACTION='WRITE',STATUS='REPLACE')
	OPEN(97,FILE='O-H_FreeEnergy.dat',ACTION='WRITE',STATUS='REPLACE')
    OPEN(98,FILE='vHtransfer.dat',ACTION='WRITE',STATUS='REPLACE')
    OPEN(99,FILE='vHtransferdfc.dat',ACTION='WRITE',STATUS='REPLACE')
    OPEN(100,FILE='vHtransferFreeEnergy.dat',ACTION='WRITE',STATUS='REPLACE')
     
    READ(90,*) natoms
    IO = 0
    nframes=0
    DO WHILE(IO == 0)
      DO i=1,natoms
         READ(90,*) atom,cxx,cyy,czz
         IF(i == idx(1)) THEN
			   IF(i < 10) THEN
			      WRITE(chidx,'(I1)') i
			   ELSE IF(i >= 10 .AND. i < 100) THEN
			      WRITE(chidx,'(I2)') i
			   ELSE IF(i >= 100 .AND. i < 1000) THEN
			      WRITE(chidx,'(I3)') i
			   ELSE
			      WRITE(chidx,'(I4)') i
		      END IF
			   chatom(1) = TRIM(atom)//TRIM(chidx)
			   cx(1) = cxx
			   cy(1) = cyy
			   cz(1) = czz
			ELSE IF(i == idx(2)) THEN   
			   IF(i < 10) THEN
			      WRITE(chidx,'(I1)') i
			   ELSE IF(i >= 10 .AND. i < 100) THEN
			      WRITE(chidx,'(I2)') i
			   ELSE IF(i >= 100 .AND. i < 1000) THEN
			      WRITE(chidx,'(I3)') i
			   ELSE
			      WRITE(chidx,'(I4)') i
		      END IF
			   chatom(2) = TRIM(atom)//TRIM(chidx)
			   cx(2) = cxx
			   cy(2) = cyy
			   cz(2) = czz			   
			ELSE IF(i == idx(3)) THEN   
			   IF(i < 10) THEN
			      WRITE(chidx,'(I1)') i
			   ELSE IF(i >= 10 .AND. i < 100) THEN
			      WRITE(chidx,'(I2)') i
			   ELSE IF(i >= 100 .AND. i < 1000) THEN
			      WRITE(chidx,'(I3)') i
			   ELSE
			      WRITE(chidx,'(I4)') i
		      END IF
			   chatom(3) = TRIM(atom)//TRIM(chidx)
			   cx(3) = cxx
			   cy(3) = cyy
			   cz(3) = czz	
		   END IF			   			   
      END DO
 
		WRITE(91,'(A6,3F14.7)') chatom(1),cx(1),cy(1),cz(1)			
		WRITE(91,'(A6,3F14.7)') chatom(2),cx(2),cy(2),cz(2)			
		WRITE(91,'(A6,3F14.7)') chatom(3),cx(3),cy(3),cz(3)
		nframes = nframes + 1			
        WRITE(91,*)
        READ(90,*,IOSTAT=IO) natoms           
    END DO
    REWIND(UNIT=91)
    ALLOCATE(bondx(nframes,2))
    ALLOCATE(vbond(nframes))
    DO i=1,nframes
      DO j=1,3
         READ(91,*) chatom(j),cx(j),cy(j),cz(j)
      END DO
      deltax=cx(2)-cx(1)
      deltay=cy(2)-cy(1)
      deltaz=cz(2)-cz(1)
      
      deltaxsq=deltax*deltax
      deltaysq=deltay*deltay
      deltazsq=deltaz*deltaz
      bondx1=SQRT(deltaxsq+deltaysq+deltazsq)
            
      deltax=cx(3)-cx(2)
      deltay=cy(3)-cy(2)
      deltaz=cz(3)-cz(2)
      
      deltaxsq=deltax*deltax
      deltaysq=deltay*deltay
      deltazsq=deltaz*deltaz
      bondx2=SQRT(deltaxsq+deltaysq+deltazsq)
      
      IF(bondx1 > innerr .AND. bondx1 < upperr) THEN
         bondx(i,1) = bondx1
         bondx(i,2) = bondx2
      ELSE
         bondx(i,1) = 0.0
         bondx(i,2) = 0.0
      END IF      
   END DO

   m = 0
   DO i=1,nframes
      IF(bondx(i,1) /= 0) THEN
         m = m + 1                                             !vbond(i) measure H transfer
         rtime = REAL(m)*dt*REAL(sampling)*atufs               !bondx(i,1) is the hydrogen bond O--H
         vbond(i) = bondx(i,1) - bondx(i,2)                    !bondx(i,2) is the O-H bond.
         WRITE(92,'(F14.4,2F14.4)') rtime,bondx(i,1),bondx(i,2)
         WRITE(98,'(F14.4,F14.4)') rtime,vbond(i)
         
      END IF 
   END DO
!//-------------------------------
		maxbin = INT(maxr/deltar)+1
		ALLOCATE(HIST(maxbin,2))
		ALLOCATE(vHISTplus(maxbin))
		ALLOCATE(vHISTminus(maxbin))
		DO bin = 1,maxbin
			HIST(bin,1) = 0.0
			HIST(bin,2) = 0.0
		END DO
		DO bin = 1,maxbin
		    vHISTplus(bin) = 0.0
		    vHISTminus(bin) = 0.0
		END DO
		
		DO i=1,nframes
			bin = INT(bondx(i,1)/deltar)+1
			IF(bondx(i,1) /= 0.0 .AND. bin <= maxbin) THEN
				HIST(bin,1) = HIST(bin,1)+1
			END IF
			bin = INT(bondx(i,2)/deltar)+1
			IF(bondx(i,2) /= 0.0 .AND. bin <= maxbin) THEN
				HIST(bin,2) = HIST(bin,2)+1
			END IF
		END DO
		
		DO i=1,nframes
		    IF(vbond(i) > 0.0) THEN
		        bin=INT(vbond(i)/deltar)+1
			    IF(bin <= maxbin) THEN
		            vHISTplus(bin) = vHISTplus(bin)+1
		        END IF
		    ELSE IF(vbond(i) < 0.0) THEN
		        bin=INT(ABS(vbond(i))/deltar)+1
			    IF(bin <= maxbin) THEN
		           vHISTminus(bin) = vHISTminus(bin)+1
		        END IF
		    END IF
		END DO
		
		soma1 = 0
		soma2 = 0
		
		DO i=1,nframes
			soma1 = soma1 + bondx(i,1)
			soma2 = soma2 + bondx(i,2)
		END DO
		
		averagehbondx1 = soma1/REAL(m)
		averagehbondx2 = soma2/REAL(m)
		
		largesthbondx1 = 0.0
		largesthbondx2 = 0.0
		
		DO i=1,nframes
		   largesthbondx1= MAX(largesthbondx1,bondx(i,1))
		   largesthbondx2= MAX(largesthbondx2,bondx(i,2))
		END DO
		
		lowesthbondx1 = 1000.0
		lowesthbondx2 = 1000.0
		
		DO i=1,nframes
		   IF(bondx(i,1) /= 0.0) THEN
		      lowesthbondx1 = MIN(lowesthbondx1,bondx(i,1))
		   END IF
		   IF(bondx(i,2) /= 0.0) THEN		   
		      lowesthbondx2 = MIN(lowesthbondx2,bondx(i,2))
		   END IF
		END DO
		DO i=1,nframes
		   IF(bondx(i,1) /= 0.0) THEN
		      bondx(i,1) = bondx(i,1) - averagehbondx1
		   END IF
		   IF(bondx(i,2) /= 0.0) THEN
		      bondx(i,2) = bondx(i,2) - averagehbondx2
		   END IF			
		END DO
		
		
		DO i=1,nframes
		   bondx(i,1) = bondx(i,1)*bondx(i,1)
		   bondx(i,2) = bondx(i,2)*bondx(i,2)
		END DO
		
		soma1 = 0.0
		soma2 = 0.0
		
		DO i=1,nframes
			soma1 = soma1 + bondx(i,1)
			soma2 = soma2 + bondx(i,2)
		END DO
		
		varhbondx1 = soma1/REAL(m-1)
		varhbondx2 = soma2/REAL(m-1)

		stdhbondx1 = SQRT(varhbondx1)
		stdhbondx2 = SQRT(varhbondx2)
		
		sumHISTbin = SUM(HIST,DIM=1)
				
		centerbin = - deltar/2
		ALLOCATE(tmp(maxbin))
		DO bin = 1,maxbin
			tmp(bin) = centerbin
		    centerbin = centerbin - deltar
		END DO
		
		sumvHISTminus = ABS(SUM(vHISTminus))
		sumvHISTplus  = SUM(vHISTplus)
		
		DO bin = 1, maxbin
		    WRITE(99,'(F7.3,2F14.7)') tmp(maxbin-bin+1),(vHISTminus(maxbin-bin+1)/(sumvHISTminus+sumvHISTplus))*100
		    IF(vHISTminus(maxbin-bin+1) /= 0.0) THEN
		        HFreeEnergy = -R*T*log(vHISTminus(maxbin-bin+1)/(sumvHISTminus+sumvHISTplus))
		        WRITE(100,'(F7.3,F14.7)') tmp(maxbin-bin+1),HFreeEnergy
		    END IF		    
		END DO
		
		DO bin = 1,maxbin
		    WRITE(99,'(F7.3,2F14.7)') ABS(tmp(bin)),(vHISTplus(bin)/(sumvHISTminus+sumvHISTplus))*100
		    IF(vHISTplus(bin) /= 0.0) THEN
		        HFreeEnergy = -R*T*log(vHISTplus(bin)/(sumvHISTminus+sumvHISTplus))
		        WRITE(100,'(F7.3,F14.7)') ABS(tmp(bin)),HFreeEnergy
		    END IF		    
		    		    
		END DO
		
		centerbin = deltar/2		
		DO bin = 1,maxbin
			WRITE(93,'(F7.3,2I12)') centerbin,HIST(bin,1),HIST(bin,2)
			
			WRITE(95,'(F7.3,2F14.4)') centerbin,(HIST(bin,1)/sumHISTbin(1))*100, &
			                        (HIST(bin,2)/sumHISTbin(2))*100
			IF(HIST(bin,1) /= 0) THEN
			    FreeEnergy = -R*T*log(HIST(bin,1)/sumHISTbin(1)) !Free energy using probability distribution function

			    WRITE(96,'(F7.3,2F14.4)') centerbin, FreeEnergy     !Hydorgen bond Free energy calculation		
			END IF
			
			IF(HIST(bin,2) /= 0) THEN
			    FreeEnergy = -R*T*log(HIST(bin,2)/sumHISTbin(2)) !Free energy using probability distribution function

			    WRITE(97,'(F7.3,2F14.4)') centerbin, FreeEnergy      	!X-H bond Free energy calculation		
			END IF		
			centerbin = centerbin + deltar
		END DO
		out1=TRIM(chatom(1))//'-'//TRIM(chatom(2))
		out2=TRIM(chatom(2))//'-'//TRIM(chatom(3))

      WRITE(94,'(A)') 'vmd.xyz (input, must be given)'
      WRITE(94,'(A)') 'tmphbond.dat (output used during the calculation)'
      WRITE(94,'(A)') 'hbond.dat (output,column 1: simulation time in fs; column 2: bondlenght'
      WRITE(94,'(A)') 'between first and second atom indexes given; column 3: bondlenght'
      WRITE(94,'(A)') 'between second and third atom indexes given).'
      WRITE(94,'(A)') 'hbonddfc (output, column 1: atomic distance in Angstrom,'
      WRITE(94,'(A)') 'column 2: distribution function for first and second atom indexes;'
      WRITE(94,'(A)') 'column 3: distribution function for second and third atom indexes).'
      WRITE(94,*)
      WRITE(94,'(A)')        '                               SUMMARY'
		WRITE(94,*)				
		WRITE(94,'(A,2A14)')   '           Selected atoms:      ',out1,out2					
		WRITE(94,'(A,2F14.7)') ' The largest bond lenghts: ',largesthbondx1,largesthbondx2
		WRITE(94,'(A,2F14.7)') '   The lowest bond lenght: ',lowesthbondx1,lowesthbondx2
		WRITE(94,'(A,2F14.7)') ' Average of bond lenghts : ',averagehbondx1,averagehbondx2
		WRITE(94,'(A,2F14.7)') '                 Variance: ',varhbondx1,varhbondx2
		WRITE(94,'(A,2F14.7)') '       Standard deviation: ',stdhbondx1,stdhbondx2
		WRITE(94,'(A,I14)')    '          Frames analized: ',nframes
		WRITE(94,'(A,I14)')    'Frames with hydrogen bond: ',m
		WRITE(94,'(A,F14.4)')  '      Simulation time(fs): ',REAL(nframes)*dt*REAL(sampling)*atufs		
		WRITE(94,'(A,F14.4)')  '   Hydrogen bond time(fs): ',REAL(m)*dt*REAL(sampling)*atufs

      DEALLOCATE(bondx)
      DEALLOCATE(HIST)
      DEALLOCATE(vHISTplus)
      DEALLOCATE(vHISTminus)
      DEALLOCATE(tmp)
      CLOSE(UNIT=90)
      CLOSE(UNIT=91)
      CLOSE(UNIT=92)
      CLOSE(UNIT=93)
      CLOSE(UNIT=94)
      CLOSE(UNIT=98)
      CLOSE(UNIT=99)
      CLOSE(UNIT=100)
      
   END SUBROUTINE SR019