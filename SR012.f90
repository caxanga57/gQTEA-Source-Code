	SUBROUTINE SR012()
	IMPLICIT NONE
	!This subroutine computes the mean residence time (MRT)
	
   REAL, PARAMETER :: atufs=0.02418884326505
	LOGICAL :: EXISTE
	INTEGER :: i,j,n,nframes,IO,cnt, k,idx,natoms,tatominshell,atominshell, &
			   nstep,err,cntnees,natomselected,natomsexcluded,&
			   natominshell
	REAL    :: cx,cy,cz,cx2,cy2,cz2,innerr,outerr,dt, &
			   deltax,deltay,deltaz,lenghtsq,lenght,mrtime,freq,&
			   sampling,soma,soma2,tsimultime,mcn,stime
	CHARACTER(LEN=3) :: atomselected, atom,atom2
	CHARACTER(LEN=6) :: chidx,chidx2,chatom
	INTEGER,DIMENSION(:,:),ALLOCATABLE :: mshell
	INTEGER,DIMENSION(:),ALLOCATABLE :: idxx
	REAL,DIMENSION(:,:), ALLOCATABLE   :: RT,T
	CHARACTER(LEN=6), DIMENSION(:), ALLOCATABLE :: tmpatoms   

	INQUIRE(FILE='vmd.xyz',EXIST=EXISTE)
	IF(EXISTE) THEN
	   OPEN(90,FILE='vmd.xyz',ACTION='READ',STATUS='OLD')
	ELSE
	   WRITE(*,*)
	   WRITE(*,'(A)') 'SORRY! YOU HAVE TO PROVIDE THE vmd.xyz FILE'
	   STOP
	END IF		
	
	
	WRITE(*,'(A)') 'GIVE ATOM INDEX AT THE CENTER OF THE SHELL'
	READ(*,*) idx
	WRITE(*,'(A)') 'GIVE ATOM SYMBOL TO CALCULATE THE RESIDENCE TIME FOR IT'
	READ(*,*) atomselected
	WRITE(*,'(A)') 'GIVE INNER RADIUS OF THE SHELL'	
	READ(*,*) innerr
	WRITE(*,'(A)') 'GIVE OUTER RADIUS OF THE SHELL'	
	READ(*,*) outerr
	WRITE(*,'(A)') 'GIVE THE MOLECULAR DYNAMICS SIMULATION TIME STEP (dt) IN ATU'
	READ(*,*) dt
	WRITE(*,'(A)') 'GIVE THE SAMPLING FREQUENCY (FOR QE IS iprint VALUE)'
	READ(*,*) sampling
   WRITE(*,'(A)') 'GIVE THE NUMBER OF ATOMS TO BE EXCLUDED FROM THE CALCULATION. TYPE ZERO FOR NO ONE'
   READ(*,*) natomsexcluded
   
	OPEN(91,FILE='tmp1.dat',ACTION='READWRITE',STATUS='REPLACE')
	OPEN(92,FILE='tmp2.dat',ACTION='READWRITE',STATUS='REPLACE')
	OPEN(93,FILE='mshell.dat',ACTION='READWRITE',STATUS='REPLACE')     
	OPEN(94,FILE='coordnumber.dat',ACTION='READWRITE',STATUS='REPLACE')    
	OPEN(95,FILE='residencetime.dat',STATUS='REPLACE')     
	OPEN(96,FILE='meanresidencetime.out',ACTION='WRITE',STATUS='REPLACE')
	OPEN(97,FILE='tmpatoms.txt',ACTION='READWRITE',STATUS='REPLACE')
   
   ALLOCATE(idxx(natomsexcluded))
   DO i=1,natomsexcluded
      WRITE(*,'(A,I3,A)') 'WRITE THE INDEX OF THE ATOM ',i, ' TO BE EXCLUDED'
      READ(*,*) idxx(i)
   END DO
		
	READ(90,*) natoms
	IO = 0
	nframes = 0                                   
	DO WHILE(IO == 0)	                 
		DO i=1,natoms                     
			READ(90,*) atom,cx,cy,cz   
			IF(i == idx) THEN
			   IF(i < 10) THEN
			      WRITE(chidx2,'(I1)') i
			   ELSE IF(i >= 10 .AND. i < 100) THEN
			      WRITE(chidx2,'(I2)') i
			   ELSE IF(i >= 100 .AND. i < 1000) THEN
			      WRITE(chidx2,'(I3)') i
			   ELSE
			      WRITE(chidx2,'(I4)') i
		      END IF
			   chatom = TRIM(atom)//TRIM(chidx2)
			   WRITE(91,'(A6,3F14.7)') chatom,cx,cy,cz			
			END IF
		END DO
		nframes = nframes + 1 		
		READ(90,*,IOSTAT=IO) natoms
	END DO
	
	REWIND(UNIT=90)                       
	
	READ(90,*) natoms                      
	natomselected=0
	DO i = 1,natoms                         
		READ(90,*) atom,cx,cy,cz         
			IF(atom == atomselected) THEN
			   DO j=1,natomsexcluded
			      IF(idxx(j) == i) THEN
			         GO TO 998
			      END IF
			   END DO
			   IF(i < 10) THEN
			      WRITE(chidx,'(I1)') i
			   ELSE IF(i >= 10 .AND. i < 100) THEN
			      WRITE(chidx,'(I2)') i
			   ELSE IF(i >= 100 .AND. i < 1000) THEN
			      WRITE(chidx,'(I3)') i
			   ELSE
			      WRITE(chidx,'(I4)') i
		      END IF
			chatom = TRIM(atom)//TRIM(chidx)
			WRITE(92,'(A6,3F14.7)') chatom,cx,cy,cz
			WRITE(97,'(A)') chatom
			natomselected = natomselected + 1         
			END IF
998		CYCLE
	END DO
	WRITE(*,'(A,I7)') 'NUMBER OF FRAMES IN vmd.xyz FILE IS ',nframes
	WRITE(*,'(A,I5,A,A2,A)') 'THERE ARE ',natomselected,' ATOMS OF ',atomselected,' IN EACH FRAME' 
	
	ALLOCATE(tmpatoms(natomselected)) !it contains the selected atoms symbols
	
	REWIND(UNIT=97)
	
	DO i=1,natomselected
      READ(97,*) tmpatoms(i)	   
	END DO
	
	WRITE(93,'(A13,1000A8)') 'TIME(ps)', (tmpatoms(j),j=1,natomselected)
	
	READ(90,*) natoms	                     
	IO = 0
	DO WHILE(IO == 0)              
		DO i=1,natoms               
			READ(90,*) atom,cx,cy,cz
			IF(atom == atomselected) THEN
			   DO j=1,natomsexcluded
			      IF(idxx(j) == i) THEN
			         GO TO 999
			      END IF
			   END DO
			   IF(i < 10) THEN
			      WRITE(chidx,'(I1)') i
			   ELSE IF(i >= 10 .AND. i < 100) THEN
			      WRITE(chidx,'(I2)') i
			   ELSE IF(i >= 100 .AND. i < 1000) THEN
			      WRITE(chidx,'(I3)') i
			   ELSE
			      WRITE(chidx,'(I4)') i
		      END IF
			   chatom = TRIM(atom)//TRIM(chidx)
			   WRITE(92,'(A6,3F14.7)') chatom,cx,cy,cz
			END IF
999	CYCLE
		END DO
		READ(90,*,IOSTAT=IO) natoms
	END DO 	
	
	REWIND(UNIT=91)                
	REWIND(UNIT=92)
	REWIND(UNIT=97)                

	ALLOCATE(mshell(nframes,natomselected),STAT=err)             
	IF(err /= 0) PRINT*, 'MEMORY CAN NOT BE ALLOCATED FOR mshell ARRAY'
	natominshell=0
	DO i = 1,nframes                          
		READ(91,*) atom,cx,cy,cz
		DO j=1,natomselected
			READ(92,*) atom2,cx2,cy2,cz2
			deltax = cx - cx2
			deltay = cy - cy2
			deltaz = cz - cz2
			lenghtsq = deltax*deltax + deltay*deltay + deltaz*deltaz
			lenght = SQRT(lenghtsq)
			IF(lenght > innerr .AND. lenght < outerr) THEN
				mshell(i,j) = 1
				natominshell=natominshell+1 !compute the total number of atoms in shell for all frames
				ELSE 
					mshell(i,j) = 0
			END IF
		END DO
	END DO
	!next loop write in mshell.dat. Col 1 is the time in ps	
	DO i = 1,nframes
	   stime = REAL(i)*REAL(dt)*REAL(sampling)*atufs/1000.0                             
		WRITE(93,'(F10.4,1000I8)') stime, (mshell(i,j),j=1,natomselected)
	END DO
    !The next loop counts the number of moleculs that is in the shell for each frame and store
    !in variable atominshell.
   tatominshell=0                         		
	DO i = 1,nframes 
		atominshell = 0   !Number of atoms in the shell by frame                 
		DO j = 1,natomselected
			IF(mshell(i,j) == 1) THEN
				atominshell = atominshell + 1
			END IF
		END DO
		stime = REAL(i)*REAL(dt)*REAL(sampling)*atufs/1000.0
		WRITE(94,'(F10.4,I6)') stime, atominshell
		tatominshell=tatominshell+atominshell 
	END DO
   mcn = REAL(tatominshell)/REAL(nframes) !compute the mean coordination number
	k = 0
	DO j = 1,natomselected             !This loop counts how many times the molecules go in and go out from the shell
		i = 1                     
		cntnees = 0
		DO WHILE(i <= nframes)			 
			IF(mshell(i,j) == 1) THEN
				cntnees = cntnees + 1
				DO WHILE(mshell(i,j) == 1)
					i=i+1
				END DO
			END IF
			i=i+1
		END DO
		IF(k < cntnees) THEN
			k = cntnees
		END IF
	END DO
	!RT is a matrix conataining the residence time for each selected atom
	ALLOCATE(RT(natomselected,k),STAT = err)
	IF(err .NE. 0) PRINT*, 'MEMORY CAN NOT BE ALLOCATED FOR RT ARRAY'
	
	DO i=1,natomselected                    
		DO j=1,k
			RT(i,j) = 0.0
		END DO
	END DO
	
	DO j=1,natomselected              !This loop counts the time that a molecule stay in the shell
		n=0                            !The results are put in RT matrix.
		nstep=0
		i=1
		DO WHILE(i <= nframes)
			IF(mshell(i,j)==1) THEN				
				DO WHILE(mshell(i,j)==1 .AND. i <= nframes)
					nstep = nstep + 1
					i=i+1
				END DO
				n=n+1
				RT(j,n)= REAL(nstep)*REAL(dt)*REAL(sampling)*atufs  !RT is given in fs
			END IF
			i=i+1
			nstep = 0
		END DO
	END DO
	
	DO i=1,natomselected                  
		WRITE(95,'(A6,400F14.7)') tmpatoms(i),(RT(i,j),j=1,k)
	END DO	
	
	ALLOCATE(T(natomselected,3))
	
	DO i=1,natomselected
		DO j=1,3
			T(i,j)=0.0
		END DO
	END DO
	
	DO i=1,natomselected
		n=0
		DO j=1,k
			IF(RT(i,j) /= 0.0) THEN
				n=n+1
				T(i,1)=T(i,1)+RT(i,j)
			END IF
		END DO
		T(i,2)= n
		IF(n /= 0) THEN
			T(i,3)=T(i,1)/REAL(n)
		END IF
	END DO
	
	WRITE(96,'(A)') 'ATOM      TOTAL TIME(fs)  #EXCHANGE  AVERAGE TIME(fs)'
	
	DO i=1,natomselected                  
		WRITE(96,'(A6,F14.4,I10,F14.4)') tmpatoms(i),T(i,1),INT(T(i,2)),T(i,3)
	END DO
	
	n = 0                              !This loop calculates the Total Mean Residence Time
	soma = 0
	DO i=1,natomselected
		IF(T(i,3) /= 0.0) THEN
			soma = soma + T(i,3)
			n = n + 1
		END IF
	END DO
	mrtime = soma/REAL(n)                   
	
   soma2=0
   DO i=1,natomselected
      soma2 = soma2 + T(i,2)
	END DO
	tsimultime=REAL(nframes)*REAL(dt)*REAL(sampling)*atufs
	freq = soma2/tsimultime
	
	WRITE(96,*)
	WRITE(96,'(A,F10.4,A)') '       SHELL INNER RADIUS = ',innerr, ' angstroms'
	WRITE(96,'(A,F10.4,A)') '       SHELL OUTER RADIUS = ',outerr, ' angstroms'
	WRITE(96,'(A,I10,A)') 'NUMBER OF FRAMES ANALIZED = ',nframes, ' frames'
	WRITE(96,'(A,F10.4,A)') '    TOTAL SIMULATION TIME = ',tsimultime,' fentosecond(fs)'	
	WRITE(96,'(A,F10.4,A)') '      MEAN RESIDENCE TIME = ',mrtime,' fentosecond(fs)'
	WRITE(96,'(A,F10.4,A)') ' MEAN COORDINATION NUMBER = ',mcn,' atoms/frame'
	WRITE(96,'(A,F12.7,A)') '       EXCHANGE FREQUENCY = ',freq,' echange/fs'
	CLOSE(UNIT=90)
	CLOSE(UNIT=91,STATUS='DELETE')
	CLOSE(UNIT=92,STATUS='DELETE')
	CLOSE(UNIT=93)
	CLOSE(UNIT=94)
	CLOSE(UNIT=95)
	CLOSE(UNIT=96)
	CLOSE(UNIT=97,STATUS='DELETE')
	DEALLOCATE(T)
    DEALLOCATE(mshell)
    DEALLOCATE(tmpatoms)
    DEALLOCATE(RT)
    DEALLOCATE(idxx)
END SUBROUTINE SR012