	SUBROUTINE SR011()   
		IMPLICIT NONE
		!This subroutine calculates the Diffusion Coefficent of a single molecule
		!using its Centre of Mass.

		INTEGER :: i, j,k,nstep=0,natmmol,tnatm, IO, n, err
		INTEGER, DIMENSION(:), ALLOCATABLE :: vecidxmol, vecidx
		REAL, DIMENSION(:), ALLOCATABLE :: amass, vecmsd, vecsmsd
		REAL, DIMENSION(:,:), ALLOCATABLE :: mxyz, tmpxyz, st, stt
		REAL, DIMENSION(2,2) :: A, B, ID, CA
		REAL, DIMENSION(2)   :: beta
		CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: vecatmsbl
		CHARACTER(LEN=6) :: keyword
		REAL :: tam=0,cmx=0,cmy=0,cmz=0, cmxold,cmyold,cmzold,step, &
				delcmx, delcmy, delcmz
		
		OPEN(11,FILE='vmd.xyz',STATUS='OLD')   
		OPEN(12,FILE='cmxyz.dat',STATUS='REPLACE')
		OPEN(13,FILE='smsd.dat',STATUS='REPLACE')
		OPEN(14,FILE='dccm.out',STATUS='REPLACE') 
		
		WRITE(*,*)'Write the number of atoms in the molecule'
		READ(*,*) natmmol
		WRITE(*,*) 'Write the step in femtosecond (fs)'
		READ(*,*) step     
		
		ALLOCATE(vecidxmol(natmmol),STAT=err)
		IF(err .NE. 0) WRITE(*,*) 'Impossible to allocate memory for vecidxmol array'
		ALLOCATE(amass(natmmol),STAT=err)
		IF(err .NE. 0) WRITE(*,*) 'Impossible to allcate memory for amass array'
		WRITE(*,*)
		
		DO i=1,natmmol 
			WRITE(*,'(1X,(A),(I3),(A),(I3))') 'Enter the index and mass of atom ',i,' out of ',natmmol                 
			READ(*,*) vecidxmol(i), amass(i)
		END DO
		
		!Calculating the total atomic mass of the molecule
		
		DO i=1,natmmol
			tam = tam + amass(i)
		END DO
		WRITE(*,*) 'The total atomic mass is ', tam

		!Reading the vmd.xyz file and calculating the centre of mass
		
		READ(11,*) tnatm
		WRITE(*,*) 'The total number of atoms is ',tnatm
		
		ALLOCATE(vecatmsbl(tnatm),STAT=err)
		IF(err .NE. 0) WRITE(*,*) 'Impossible to allocate memory for vecatmsbl array'		
		ALLOCATE(mxyz(tnatm,3),STAT=err)
		IF(err .NE. 0) WRITE(*,*) 'Impossible to allocate memory for mxyz matrix'		
		ALLOCATE(tmpxyz(natmmol,3),STAT=err)
		IF(err .NE. 0) WRITE(*,*) 'Impossible to allocate memory for tmpxyz matrix'
		ALLOCATE(vecidx(tnatm),STAT=err)
		IF(err .NE. 0) WRITE(*,*) 'Impossible to allocate memory for vecidx array'


		IO = 0		
		DO WHILE(IO == 0)
			DO i=1,tnatm
				READ(11,*) vecatmsbl(i), (mxyz(i,j),j=1,3)
				vecidx(i)=i
			END DO
			DO i=1,tnatm
				DO j=1,natmmol
					IF( vecidxmol(j) == vecidx(i) ) THEN
!						WRITE(*,*) vecatmsbl(i), (mxyz(i,k),k=1,3)
						cmx = cmx + amass(j)*mxyz(i,1)
						cmy = cmy + amass(j)*mxyz(i,2)
						cmz = cmz + amass(j)*mxyz(i,3)
					END IF	
				END DO
			END DO
			WRITE(12,*) cmx/tam, cmy/tam, cmz/tam
			nstep = nstep + 1	
			cmx = 0
			cmy = 0
			cmz = 0								
			READ(11,*,IOSTAT=IO) tnatm
		END DO		
		REWIND(12)
		
		!Calculating the mean square displacement (MSD) between 2 steps
		
		ALLOCATE(vecmsd(nstep-1),STAT=err)
		IF(err .NE. 0) WRITE(*,*) 'Impossible to allocate memory for vecmsd array'
		ALLOCATE(vecsmsd(nstep-1),STAT=err)
		IF(err .NE. 0) WRITE(*,*) 'Impossible to allocate memory for vecsmsd array'
		
		READ(12,*) cmxold,cmyold,cmzold
		DO i=1,nstep-1
			READ(12,*) cmx,cmy,cmz
			delcmx = cmx-cmxold
			delcmy = cmy-cmyold
			delcmz = cmz-cmzold
			vecmsd(i)= delcmx*delcmx + delcmy*delcmy + delcmz*delcmz
			cmxold = cmx
			cmyold = cmy
			cmzold = cmz
		END DO
		
		!Calculating the partial sum of msd, i.e., vecsmsd
		
		vecsmsd(1)= vecmsd(1)
		DO i=2,nstep-1
			vecsmsd(i)= vecsmsd(i-1) + vecmsd(i)
		END DO
		
		!Calculating the time simulation matrix st
		
		ALLOCATE(st(nstep-1,2),STAT=err)
		IF(err .NE. 0) WRITE(*,*) 'Impossible to allocate memory for st matrix'
		st(1,1) = step
		st(1,2) = 1.0
		DO i=2,nstep-1
			st(i,1) = st(i-1,1)+step
			st(i,2) = 1.0
		END DO
		
		DO i=1,nstep-1
			WRITE(13,*) st(i,1), vecsmsd(i)
		END DO
		
		!Calating the transpose of st, i.e., stt
		
		ALLOCATE(stt(2,nstep-1),STAT=err)
		IF(err .NE. 0) WRITE(*,*) 'Impossible to allocate memory for stt matrix'
		stt = TRANSPOSE(st)
		
		!Multiplying stt*st
		
		A = MATMUL(stt,st)
		CA(1,1)=A(1,1)
		CA(1,2)=A(1,2)
		CA(2,1)=A(2,1)
		CA(2,2)=A(2,2)
		
		!Calculating the inverse of A
		n = 2
		CALL matrixinv(A,B,n)
		
		ID=MATMUL(CA,B)
		WRITE(*,*)
		WRITE(*,*) 'Number of frames ',nstep
		WRITE(14,*) 'Number of frames ',nstep
		WRITE(*,*) 'Total time of simulation ',nstep*step,'fs'
		WRITE(14,*) 'Total time of simulation ',nstep*step,'fs'

		WRITE(*,*)
		WRITE(*,*) 'Product of A*inv(A)'
		WRITE(*,*)
		WRITE(*,*) ID(1,1),ID(1,2)
		WRITE(*,*) ID(2,1),ID(2,2)
		
		!Calculating the regretion vector beta
		
		beta = MATMUL(MATMUL(B,stt),vecsmsd)
		WRITE(*,*)
		WRITE(14,*) 
		WRITE(*,*) 'Regression equation by least square'
		WRITE(14,*) 'Regression equation by least square'
		WRITE(14,*)
		WRITE(*,*)
    	WRITE(*,*) 'vecsmsd=',beta(1),'*t +',beta(2) 
    	WRITE(14,*), 'vecsmsd=',beta(1),'*t +',beta(2) 
		WRITE(14,*)
    	WRITE(*,*)
    	WRITE(*,*) 'DIFFUSION COEFFICIENT D=',beta(1)/6.0,' A**2/fs ou ',0.1*(beta(1)/6.0),'cm**2/s'
    	WRITE(14,*), 'DIFFUSION COEFFICIENT D=',beta(1)/6.0,' A**2/fs ou ',0.1*(beta(1)/6.0),'cm**2/s'
 		
    	CLOSE(UNIT=11)
    	CLOSE(UNIT=12)
    	CLOSE(UNIT=13)
    	CLOSE(UNIT=14)
    			
	END SUBROUTINE SR011
