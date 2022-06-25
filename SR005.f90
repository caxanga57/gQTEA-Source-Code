   SUBROUTINE SR005()
   IMPLICIT NONE
!  This subroutine calculates bond lenghts, average bond lenght, bond lenght distribution function and
!  Free Energy from probabilyt distribution function. The result is given in kcal/K.mol
!  Bond lenght is given in angstroms	
!  G=-R*T*ln(P(r)), where P(r) is the probability distribution function P(r)= (number of frames per bin)/(total of frames)
!  of the bond lenghts as a function of r; R is Boltzmann constant and T is the simulation temperature
!  J Mol Model (2011) 17:2159–2168 DOI 10.1007/s00894-010-0939-6
!  Science 275, 817 (1997) DOI: 10.1126/science.275.5301.817
			
	REAL, PARAMETER :: atufs=0.02418884326505, ps=1000.0, R=0.001987204 
	LOGICAL :: EXISTE
 	INTEGER :: idx1,idx2, tnatm, nstep, i, IO,maxbin,bin,&
	           nbond,ibond,iprint
	REAL :: cx,cy,cz,cx1,cy1,cz1,cx2,cy2,cz2,delx,dely,delz,rijsq,rij,soma,std, &
			averagebond, var,maxr,delr,centerbin,largestbond,lowestbond,&
			dammy,dt,timeps, T,HISTbinNormalized,FreeEnergy,sumHISTbin
			
	CHARACTER(LEN=3) :: atype, atype1,atype2,ext
	CHARACTER(LEN=7) :: out1,out2,chidx1,chidx2
	REAL,DIMENSION(:), ALLOCATABLE :: vecbond, HIST
	CHARACTER(LEN=40) :: outfile1, outfile2,outfile3,outfile4
	
	INQUIRE(FILE='vmd.xyz',EXIST=EXISTE)
	IF(EXISTE) THEN
	   OPEN(10,FILE='vmd.xyz',ACTION='READ',STATUS='OLD')
	ELSE
	WRITE(*,*)
	WRITE(*,'(A)') 'SORRY! YOU HAVE TO PROVIDE THE vmd.xyz FILE'
	WRITE(*,*)
	WRITE(*,'(A)') 'TYPE RETURN TO EXIT'
	READ(*,*)
	STOP
	END IF


	WRITE(*,'(A)') 'HOW MANY BOND LENGHT DO YOU WANT TO CALCULATE?'
	READ(*,*) nbond		
	WRITE(*,'(A)') 'GIVE THE MAX BOND LENGHT TO CALCULATE THE BOND LENGHT DISTRIBUTION FUNCTION'
	READ(*,*) maxr
	WRITE(*,'(A)') 'GIVE THE WIDTH OF THE HISTOGRAM BINS IN ANGSTROM'
	READ(*,*) delr
	WRITE(*,'(A)') 'GIVE dt (THIS IS THE TIME STEP IN ATOMIC UNIT USED IN QE INPUT FILE)'
	READ(*,*) dt
	WRITE(*,'(A)') 'GIVE THE SAMPLING INTERVAL (THIS IS THE iprint USED IN QE INPUT FILE)'
	READ(*,*) iprint
	WRITE(*,'(A)') 'GIVE THE SIMULATION TEMPERATURE IN KELVIN'
	READ(*,*) T
		
	DO ibond=1,nbond
		WRITE(*,'(A,I3)') 'GIVE TWO ATOMS INDEXES TO CALCULATE BOND LENGHT',ibond
		READ(*,*) idx1, idx2
        WRITE(*,*)
		IF(ibond < 10 ) THEN
			WRITE(ext,'(''00'',I1)') ibond
			outfile1 = 'bond-'//ext//'.dat'
			outfile2 = 'bonddfc-'//ext//'.dat'
			outfile3 = 'bond-out-'//ext//'.dat'
			outfile4 = 'bondFreeEnergy-'//ext//'.dat'

		ELSE IF(ibond > 99) THEN
			WRITE(ext,'(I3)') ibond
			outfile1 = 'bond-'//ext//'.dat'
			outfile2 = 'bonddfc-'//ext//'.dat'
			outfile3 = 'bond-out-'//ext//'.dat'
			outfile4 = 'bondFreeEnergy-'//ext//'.dat'

		ELSE
		    WRITE(ext,'(''0'',I2)') ibond
		    outfile1 = 'bond-'//ext//'.dat'
			outfile2 = 'bonddfc-'//ext//'.dat'
			outfile3 = 'bond-out-'//ext//'.dat'
			outfile4 = 'bondFreeEnergy-'//ext//'.dat'

		END IF		

		OPEN(11,FILE=outfile1,ACTION='READWRITE',STATUS='REPLACE')
		OPEN(12,FILE=outfile2,ACTION='WRITE',STATUS='REPLACE')
		OPEN(13,FILE=outfile3,ACTION='WRITE',STATUS='REPLACE')		
		OPEN(14,FILE=outfile4,ACTION='WRITE',STATUS='REPLACE')		
		
		READ(10,*) tnatm
		IO = 0
		nstep = 0
		DO WHILE( IO == 0)
			DO i=1,tnatm
				READ(10,*) atype,cx,cy,cz
				IF(i == idx1) THEN
					atype1 = atype
					cx1 = cx
					cy1 = cy
					cz1 = cz
				ELSE IF(i == idx2) THEN
					atype2 = atype
					cx2 = cx
					cy2 = cy
					cz2 = cz
				END IF
			END DO
			delx = cx2 - cx1
			dely = cy2 - cy1
			delz = cz2 - cz1
			rijsq = delx*delx+dely*dely+delz*delz
			rij =SQRT(rijsq)
			timeps = REAL((nstep+1))*dt*REAL(iprint)*atufs/ps
			WRITE(11,'(1X,F14.7,F14.7)') timeps, rij
			nstep = nstep+1
			READ(10,*,IOSTAT=IO) tnatm		
		END DO

		REWIND(UNIT=11)
		maxbin = INT(maxr/delr)+1
		ALLOCATE(HIST(maxbin))
		ALLOCATE(vecbond(nstep))
		DO bin = 1,maxbin
			HIST(bin) = 0.0
		END DO
		
		DO i=1,nstep
			READ(11,*) dammy,vecbond(i)
			bin = INT(vecbond(i)/delr)+1
			IF(bin .LE. maxbin) THEN
				HIST(bin) = HIST(bin)+1
			END IF
		END DO
		soma = 0
		DO i=1,nstep
			soma = soma + vecbond(i)
		END DO
		averagebond = soma/REAL(nstep)
		largestbond= MAXVAL(vecbond)
		lowestbond = MINVAL(vecbond)
		
		DO i=1,nstep
			vecbond(i) = vecbond(i) - averagebond
		END DO
		
		DO i=1,nstep
			vecbond(i) = vecbond(i)*vecbond(i)
		END DO
		
		var = sum(vecbond)/REAL(nstep-1)
		std = SQRT(sum(vecbond)/REAL(nstep-1))
		
		centerbin = delr/2
		sumHISTbin = SUM(HIST)
		
		DO bin = 1,maxbin
		    HISTbinNormalized=(HIST(bin)/sumHISTbin)*100
			WRITE(12,'(F7.3,2F14.3)') centerbin,HIST(bin),HISTbinNormalized
			IF(HIST(bin) /= 0) THEN
			    FreeEnergy = -R*T*log(HIST(bin)/sumHISTbin) !Free energy using probability distribution function
			    WRITE(14,'(F7.3,F14.3)') centerbin, FreeEnergy     			
			END IF
			centerbin = centerbin + delr
		END DO
	IF(idx1 < 10) THEN
    	WRITE(chidx1,'(I1)') idx1	
    	ELSE IF(idx1 >= 100) THEN
    	WRITE(chidx1,'(I3)') idx1
    	ELSE 
    	WRITE(chidx1,'(I2)') idx1
    END IF
	IF(idx2 < 10) THEN
    	WRITE(chidx2,'(I1)') idx2	
    	ELSE IF(idx2 >= 100) THEN
    	WRITE(chidx2,'(I3)') idx2
    	ELSE 
    	WRITE(chidx2,'(I2)') idx2
    END IF		
    out1 = TRIM(atype1)//TRIM(chidx1)
    out2 = TRIM(atype2)//TRIM(chidx2)
				
		WRITE(13,'(1X,A,2A5)')   '                 Selected atoms: ',out1,out2						
		WRITE(13,'(1X,A,I14)')   '         Total number of frames: ',nstep
		WRITE(13,'(1X,A,F14.7)') '        The largest bond lenght: ',largestbond
		WRITE(13,'(1X,A,F14.7)') '         The lowest bond lenght: ',lowestbond
		WRITE(13,'(1X,A,F14.7)') '           Bond lenght average : ',averagebond
		WRITE(13,'(1X,A,F14.7)') '          Bond lenght variance : ', var
		WRITE(13,'(1X,A,F14.7)') 'Bond lenght standard deviation : ', std

		REWIND(UNIT=10)
		CLOSE(UNIT=11)
		CLOSE(UNIT=12)
		CLOSE(UNIT=13)
		DEALLOCATE(HIST)
		DEALLOCATE(vecbond)
	END DO
	CLOSE(UNIT=10)
 END SUBROUTINE SR005