   SUBROUTINE SR007()
	IMPLICIT NONE
    !This subroutine calculates dihedral angle, dihedral angle distribution function, etc
    
    
	REAL, PARAMETER :: PI = 3.1415926,atufs=0.02418884326505,ps=1000.0
	LOGICAL :: EXISTE		
	INTEGER :: idx1,idx2,idx3,idx4,natoms,IO,nstep,bin,maxbin, &
					err,i,k,j,idihedral,ndihedral,iprint
	CHARACTER(LEN=3) :: atype1,atype2,atype3,atype4,ext
	CHARACTER(LEN=40) :: outfile1,outfile2,outfile3,outfile4
	CHARACTER(LEN=7)  :: out1,out2,out3,out4,chidx1,chidx2,chidx3,chidx4
	REAL :: cx1,cy1,cz1,cx2,cy2,cz2,cx3,cy3,cz3, &
	      cx4,cy4,cz4,ax,ay,az,bx,by,bz,cx,cy,cz,normb,&
		crosbcx,crosbcy,crosbcz,crosabx,crosaby,crosabz, &
		atan2x,atan2y, phi,deltatheta,phidegree,vardihedral,stddihedral,&
		averagedihedral,timeps,dt,dammy,largestangle,lowestangle,&
		incrdihedral
	REAL, DIMENSION(:), ALLOCATABLE :: vecdihedral,dihedral
	REAL, DIMENSION(:,:), ALLOCATABLE :: cxyz,tmp,HIST
	CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: atoms

	INQUIRE(FILE='vmd.xyz',EXIST=EXISTE)
	IF(EXISTE) THEN
	   OPEN(90,FILE='vmd.xyz',ACTION='READ',STATUS='OLD')
	ELSE
	WRITE(*,*)
	WRITE(*,'(A)') 'SORRY! YOU HAVE TO PROVIDE THE vmd.xyz FILE'
	WRITE(*,*)
	WRITE(*,'(A)') 'TYPE RETURN TO EXIT'
	READ(*,*)
	STOP
	END IF

			
	WRITE(*,'(A)') 'HOW MANY DIHEDRAL ANGLES DO YOU WANT TO CALCULATE?'
	READ(*,*) ndihedral
	WRITE(*,'(A)') 'GIVE DELTA THETA TO CALCULATE THE DIHEDRAL DISTRIBUTION FUNCTION'
	READ(*,*) deltatheta
	WRITE(*,'(A)') 'GIVE dt (THIS IS THE TIME STEP IN ATOMIC UNIT USED IN QE INPUT FILE)'
	READ(*,*) dt
	WRITE(*,'(A)') 'GIVE THE SAMPLING INTERVAL (THIS IS THE iprint USED IN QE INPUT FILE)'
	READ(*,*) iprint
				
	DO idihedral=1,ndihedral			
	   WRITE(*,'(A,I3)') 'GIVE 4 ATOMIC INDEXES TO CALCULATE THE DIHEDRAL ANGLE ', idihedral
		READ(*,*) idx1,idx2,idx3,idx4
      WRITE(*,*)
		IF(idihedral < 10 ) THEN
			WRITE(ext,'(''00'',I1)') idihedral
			outfile1 = 'dihedral-'//ext//'.dat'
			outfile2 = 'dihedraldfc1-'//ext//'.dat'
         outfile3 = 'dihedral-out-'//ext//'.dat'
         outfile4 = 'dihedraldfc2-'//ext//'.dat'
		ELSE IF(idihedral > 99) THEN
			WRITE(ext,'(I3)') idihedral
			outfile1 = 'dihedral-'//ext//'.dat'
			outfile2 = 'dihedraldfc1-'//ext//'.dat'
            outfile3 = 'dihedral-out-'//ext//'.dat'
            outfile4 = 'dihedraldfc2-'//ext//'.dat'
		ELSE
         WRITE(ext,'(''0'',I2)') idihedral
			outfile1 = 'dihedral-'//ext//'.dat'
			outfile2 = 'dihedraldfc1-'//ext//'.dat'
            outfile3 = 'dihedral-out-'//ext//'.dat'
            outfile4 = 'dihedraldfc2-'//ext//'.dat'
		END IF		
		
		OPEN(91,FILE=outfile1,ACTION='READWRITE',STATUS='REPLACE')
		OPEN(92,FILE=outfile2,ACTION='READWRITE',STATUS='REPLACE')
		OPEN(93,FILE=outfile3,ACTION='WRITE',STATUS='REPLACE')
		OPEN(94,FILE=outfile4,ACTION='WRITE',STATUS='REPLACE')
		
		READ(90,*) natoms
		
		ALLOCATE(atoms(natoms),STAT=err)
		IF(err /= 0) PRINT*, 'MEMORY CAN NOT BE ALLOCATED FOR atoms ARRAY'
		ALLOCATE(cxyz(natoms,3),STAT=err)
		IF(err /= 0) PRINT*, 'MEMORY CAN NOT BE ALLOCATED FOR cxyz ARRAY'
		
		IO = 0
		nstep = 0
		DO WHILE( IO == 0)
			DO i=1,natoms
				READ(90,*) atoms(i), (cxyz(i,j),j=1,3)
			END DO	
			atype1 = atoms(idx1)
			cx1 = cxyz(idx1,1)
			cy1 = cxyz(idx1,2)
			cz1 = cxyz(idx1,3)
			
			atype2 = atoms(idx2)
			cx2 = cxyz(idx2,1)
			cy2 = cxyz(idx2,2)
			cz2 = cxyz(idx2,3)
			
			atype3 = atoms(idx3)
			cx3 = cxyz(idx3,1)
			cy3 = cxyz(idx3,2)
			cz3 = cxyz(idx3,3)
			
			atype4 = atoms(idx4)
			cx4 = cxyz(idx4,1)
			cy4 = cxyz(idx4,2)
			cz4 = cxyz(idx4,3)
			
			ax = cx2 - cx1
			ay = cy2 - cy1
			az = cz2 - cz1
			bx = cx3 - cx2
			by = cy3 - cy2
			bz = cz3 - cz2
			cx = cx4 - cx3
			cy = cy4 - cy3
			cz = cz4 - cz3
			normb = SQRT(bx*bx+by*by+bz*bz)
			
			crosbcx = by*cz - cy*bz
			crosbcy = bz*cx - cz*bx
			crosbcz = bx*cy - cx*by
			
			crosabx = ay*bz - by*az
			crosaby = az*bx - bz*ax
			crosabz = ax*by - bx*ay
			
			atan2x = normb*(ax*crosbcx+ay*crosbcy+az*crosbcz)
			atan2y = crosabx*crosbcx+crosaby*crosbcy+crosabz*crosbcz
			
			phi = atan2(atan2x,atan2y)
			
			phidegree = phi*180.0/PI
			timeps = REAL((nstep+1))*dt*REAL(iprint)*atufs/ps			
			WRITE(91,'(2F14.7)') timeps,phidegree
			nstep = nstep+1
			READ(90,*,IOSTAT=IO) natoms		
		END DO
				
		REWIND(UNIT=91)
		
		maxbin = INT(360.0/deltatheta) 
        ALLOCATE(HIST(maxbin,2),STAT=err)
		IF(err /= 0) PRINT*, 'MEMORY CAN NOT BE ALLOCATED FOR HIST ARRAY'        		
		ALLOCATE(vecdihedral(nstep),STAT=err)
		IF(err /= 0) PRINT*, 'Impossible to allocate memory for vecdihedral array'
		ALLOCATE(dihedral(nstep),STAT=err)
		IF(err /= 0) PRINT*, 'MEMORY CAN NOT BE ALLOCATED FOR DIHEDRAL VECTOR'
		
		incrdihedral = deltatheta/2.0
		DO bin=1,maxbin
			HIST(bin,1) = incrdihedral
			HIST(bin,2) = 0.0
			incrdihedral = incrdihedral+deltatheta
		END DO
		
		DO i=1,nstep
			READ(91,*) dammy,vecdihedral(i)
			dihedral(i) = vecdihedral(i)
		END DO
		
		DO i=1,nstep
		   IF(vecdihedral(i) > 0.0) THEN
				bin=INT(vecdihedral(i)/deltatheta)+1
				HIST(bin,2)=HIST(bin,2)+1
			ELSE IF(vecdihedral(i) < 0.0) THEN
				bin=INT((vecdihedral(i)+360.0)/deltatheta)+1
				HIST(bin,2)=HIST(bin,2)+1
			END IF
		END DO
		k=0
		DO i=1,maxbin
			WRITE(92,'(F6.2,F14.7)') HIST(i,1),HIST(i,2)
			IF(HIST(maxbin-i+1,1) > 180.0) THEN
			   k=k+1
			END IF
		END DO
		ALLOCATE(tmp(k,2))
		DO i=1,k
			tmp(i,1)= -(360.0-HIST(maxbin-i+1,1))
			tmp(i,2)=HIST(maxbin-i+1,2)
		END DO
		DO i=1,k
			WRITE(94,'(F8.2,F16.7)') tmp(k-i+1,1),tmp(k-i+1,2)
		END DO
		DO i=1,maxbin
			IF(HIST(i,1) < 180.0) THEN
			WRITE(94,'(F8.2,F16.7)') HIST(i,1),HIST(i,2)
			END IF
		END DO
		
		averagedihedral = SUM(dihedral)/REAL(nstep)
		largestangle = MAXVAL(dihedral)
		lowestangle  = MINVAL(dihedral)
		DO i=1,nstep
			dihedral(i) = dihedral(i) - averagedihedral
			dihedral(i) = dihedral(i)*dihedral(i)
		END DO
		vardihedral = SUM(dihedral)/REAL(nstep-1)
		stddihedral = SQRT(vardihedral)
				
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
		IF(idx4 < 10) THEN
    	WRITE(chidx4,'(I1)') idx4	
    	ELSE IF(idx4 >= 100) THEN
    	WRITE(chidx4,'(I3)') idx4
    	ELSE 
    	WRITE(chidx4,'(I2)') idx4
    	END IF        
    	out1 = TRIM(atype1)//TRIM(chidx1)
    	out2 = TRIM(atype2)//TRIM(chidx2)
    	out3 = TRIM(atype3)//TRIM(chidx3)
    	out4 = TRIM(atype4)//TRIM(chidx4)
    
    	WRITE(93,'(A)')              '                    SUMMARY OF CALCULATIONS'
    	WRITE(93,*)
	   WRITE(93,'(1X,A,4A5)')       '                            Selected atoms: ',out1,out2,out3,out4
	   WRITE(93,'(1X,A,I14)')       '                    Total number of frames: ',nstep
	   WRITE(93,'(1X,A,F14.7)')     '              Dihedral angle average value: ',averagedihedral
	   WRITE(93,'(1X,A,F14.7)')     '                    Largest dihedral angle: ',largestangle
	   WRITE(93,'(1X,A,F14.7)')     '                     Lowest dihedral angle: ',lowestangle
	   WRITE(93,'(1X,A,F14.7)')     '                   Dihedral angle variance: ',vardihedral
	   WRITE(93,'(1X,A,F14.7)')     '         Dihedral angle standard deviation: ',stddihedral
				
	   CLOSE(UNIT=91)
	   CLOSE(UNIT=92)
	   CLOSE(UNIT=93)
	   REWIND(UNIT=90)
	   DEALLOCATE(vecdihedral)
	   DEALLOCATE(atoms)
	   DEALLOCATE(cxyz)
	   DEALLOCATE(dihedral)
	   DEALLOCATE(HIST)
	   DEALLOCATE(tmp)	
	END DO
	CLOSE(UNIT=90)
	
	END SUBROUTINE SR007