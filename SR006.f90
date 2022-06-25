   SUBROUTINE SR006()
   IMPLICIT NONE
    !This subroutine calculates angle, angle distribution function, etc
    
    
	REAL, PARAMETER :: PI = 3.1415926,atufs=0.02418884326505, ps=1000.0, &
	                   R=0.001987204  !Kcal/KMol
	LOGICAL :: EXISTE
	INTEGER :: idx1,idx2,idx3, tnatm, nstep, i, &
			   IO, bin, maxbin,j,err,nangles,iangle,iprint
	REAL :: cx1,cy1,cz1,cx2,cy2,cz2,cx3,cy3,cz3, &
			soma,std, averageangle, var,vecx1,vecy1,vecz1, &
			vecx2,vecy2,vecz2,normavec1,normavec2,innerprod,&
			cosalfa,alfarad,alfadegree,delteta,centerbin,&
			dt,timeps,dammy,largestangle,lowestangle,T,maxHISTbin,&
			Lnormalized1,Lnormalized2,FreeEnergy
	CHARACTER(LEN=3) :: atype, atype1,atype2,atype3,ext, &
					chidx1,chidx2,chidx3
	CHARACTER(LEN=5) :: out1, out2, out3
	CHARACTER(LEN=30) :: outfile1,outfile2,outfile3,outfile4
	REAL,DIMENSION(:), ALLOCATABLE :: vecangle, HIST
	REAL, DIMENSION(:,:), ALLOCATABLE :: mxyz
	CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: vecatmsbl

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
	
		
	WRITE(*,'(A)') 'HOW MANY ANGLES DO YOU WANT TO CALCULATE?'
	READ(*,*)  nangles
	WRITE(*,'(A)')
	WRITE(*,'(A)') 'GIVE DELTA THETA TO CALCULATE BOND ANGLE DISTRIBUTION FUNCTION'
	READ(*,*) delteta
	WRITE(*,'(A)') 'GIVE dt (THIS IS THE TIME STEP IN ATOMIC UNIT USED IN QE INPUT FILE)'
	READ(*,*) dt
	WRITE(*,'(A)') 'GIVE THE SAMPLING INTERVAL (THIS IS THE iprint USED IN QE INPUT FILE)'
	READ(*,*) iprint
	WRITE(*,'(A)') 'GIVE THE SIMULATION TEMPERATURE IN KELVIN'
	READ(*,*) T
			
 DO iangle=1,nangles
	
	WRITE(*,'(A,I3)') 'GIVE 3 ATOMIC INDEXES FOR ANGLE ', iangle
	READ(*,*) idx1,idx2, idx3
	
	IF(iangle < 10 ) THEN
	  WRITE(ext,'(''00'',I1)') iangle
	  outfile1 = 'angle-'//ext//'.dat'
	  outfile2 = 'angledfc-'//ext//'.dat'
	  outfile3 = 'angle-out-'//ext//'.dat'
      outfile4 = 'angleFreeEnergy-'//ext//'.dat'

	ELSE IF(iangle > 99) THEN
	  WRITE(ext,'(I3)') iangle
	  outfile1 = 'angle-'//ext//'.dat'
	  outfile2 = 'angledfc-'//ext//'.dat'
	  outfile3 = 'angle-out-'//ext//'.dat'
      outfile4 = 'angleFreeEnergy-'//ext//'.dat'

	ELSE
	  WRITE(ext,'(''0'',I2)') iangle
	  outfile1 = 'angle-'//ext//'.dat'
	  outfile2 = 'angledfc-'//ext//'.dat'
	  outfile3 = 'angle-out-'//ext//'.dat'
      outfile4 = 'angleFreeEnergy-'//ext//'.dat'
	  
	END IF
			
	 OPEN(11,FILE=outfile1,ACTION='READWRITE',STATUS='REPLACE')
	 OPEN(12,FILE=outfile2,ACTION='WRITE',STATUS='REPLACE')
	 OPEN(13,FILE=outfile3,ACTION='WRITE',STATUS='REPLACE')			
     OPEN(14,FILE=outfile4,ACTION='WRITE',STATUS='REPLACE')		
			
	 READ(10,*) tnatm
	 ALLOCATE(vecatmsbl(tnatm),STAT=err)
	 IF(err .NE. 0) PRINT*, 'Impossible to allocate memory for vecatmsbl array'
	 ALLOCATE(mxyz(tnatm,3),STAT=err)
	 IF(err .NE. 0) PRINT*, 'Impossible to allocate memory for mxyz array'
		
	 IO = 0
	 nstep = 0
	 DO WHILE( IO == 0)
		DO i=1,tnatm
			READ(10,*) vecatmsbl(i), (mxyz(i,j),j=1,3)
		END DO
			
		atype1 = vecatmsbl(idx1)
		cx1 = mxyz(idx1,1)
		cy1 = mxyz(idx1,2)
		cz1 = mxyz(idx1,3)

		atype2 = vecatmsbl(idx2)
		cx2 = mxyz(idx2,1)
		cy2 = mxyz(idx2,2)
		cz2 = mxyz(idx2,3)

		atype3 = vecatmsbl(idx3)
		cx3 = mxyz(idx3,1)
		cy3 = mxyz(idx3,2)
		cz3 = mxyz(idx3,3)

		vecx1 = cx1 - cx2
		vecy1 = cy1 - cy2
		vecz1 = cz1 - cz2
		vecx2 = cx3 - cx2
		vecy2 = cy3 - cy2
		vecz2 = cz3 - cz2
			
		innerprod = vecx1*vecx2+vecy1*vecy2+vecz1*vecz2
		normavec1 = SQRT(vecx1*vecx1+vecy1*vecy1+vecz1*vecz1)
		normavec2 = SQRT(vecx2*vecx2+vecy2*vecy2+vecz2*vecz2)
		cosalfa   = innerprod/(normavec1*normavec2)
		alfarad = ACOS(cosalfa)
		alfadegree = alfarad*180.0/PI
		timeps = REAL((nstep+1))*dt*REAL(iprint)*atufs/ps		
		WRITE(11,'(F14.7,F14.7)') timeps, alfadegree
		nstep = nstep+1
		READ(10,*,IOSTAT=IO) tnatm		
	END DO
		
	REWIND(UNIT=11)
	maxbin = INT(180.0/delteta)+1
		
	ALLOCATE(HIST(maxbin))
	ALLOCATE(vecangle(nstep))
		
	DO bin=1,maxbin
		HIST(bin) = 0.0
	END DO
		
	DO i=1,nstep
		READ(11,*) dammy, vecangle(i)
		bin = INT(vecangle(i)/delteta)
		IF(bin <= maxbin) THEN
 			HIST(bin) = HIST(bin)+1
 		END IF
	END DO
	largestangle = MAXVAL(vecangle)
	lowestangle = MINVAL(vecangle)
	soma = 0
	DO i=1,nstep
		soma = soma + vecangle(i)
	END DO
	averageangle = soma/REAL(nstep)
		
	DO i=1,nstep
		vecangle(i) = vecangle(i) - averageangle
	END DO
		
	DO i=1,nstep
		vecangle(i) = vecangle(i)*vecangle(i)
	END DO
	
	maxHISTbin = MAXVAL(HIST)	
	centerbin = delteta/2
	
	DO bin = 1,maxbin
    Lnormalized1=HIST(bin)/maxHISTbin
    Lnormalized2=Lnormalized1*100
	WRITE(12,'(F7.3,3F14.3)') centerbin,HIST(bin),Lnormalized1,Lnormalized2

! 		WRITE(12,*) centerbin, HIST(bin)
	IF(HIST(bin) /= 0) THEN
	    FreeEnergy = -R*T*log(HIST(bin)/maxHISTbin)
		WRITE(14,'(F7.3,F14.3)') centerbin, FreeEnergy      	!Free energy calculation		
	END IF
	centerbin = centerbin + delteta
	END DO
		
	var = sum(vecangle)/REAL(nstep-1)
	std = SQRT(sum(vecangle)/REAL(nstep-1))
	
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
	IF(idx3 < 10) THEN
    	WRITE(chidx3,'(I1)') idx3	
    	ELSE IF(idx3 >= 100) THEN
    	WRITE(chidx3,'(I3)') idx3
    	ELSE 
    	WRITE(chidx3,'(I2)') idx3
    END IF
    out1 = TRIM(atype1)//TRIM(chidx1)
    out2 = TRIM(atype2)//TRIM(chidx2)
    out3 = TRIM(atype3)//TRIM(chidx3)
	WRITE(13,'(A,3A5)')   '          SELECTED ATOMS: ', out1,out2,out3
	WRITE(13,'(A,I14)')   '   NUMBER OF FRAMES USED: ',nstep
	WRITE(13,'(A,F14.7)') '       THE LARGEST ANGLE: ',largestangle
	WRITE(13,'(A,F14.7)') '        THE LOWEST ANGLE: ',lowestangle	
	WRITE(13,'(A,F14.7)') '           ANGLE AVERAGE: ',averageangle
	WRITE(13,'(A,F14.7)') '          ANGLE VARIANCE: ',var
	WRITE(13,'(A,F14.7)') 'ANGLE STANDARD DEVIATION: ',std	
	REWIND(UNIT=10)
	CLOSE(UNIT=11)
	CLOSE(UNIT=12)
	CLOSE(UNIT=13)
	DEALLOCATE(HIST)
	DEALLOCATE(vecangle)
	DEALLOCATE(vecatmsbl)
	DEALLOCATE(mxyz)
   END DO
   CLOSE(UNIT=10)		
   END SUBROUTINE SR006
 