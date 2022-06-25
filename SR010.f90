	SUBROUTINE SR010()
	IMPLICIT NONE
    !This subroutine calculates the g(r) without PBC
    
 		REAL, PARAMETER :: PI = 3.1415926
 		INTEGER :: tnatm,IO,i,j,idx,maxbin,bin,nstep,cntatm,err
 		REAL, DIMENSION(:,:), ALLOCATABLE :: mxyz
 		REAL, DIMENSION(:),   ALLOCATABLE :: HIST, vecr,vecgr,vecintgr
 		CHARACTER(LEN=3),DIMENSION(:), ALLOCATABLE    :: vecatmsbl
 		REAL :: a,b,c, rijsq,rij,delr, maxr, rho,delx,dely,delz, &
 			    rlower,rupper,const,nideal
 		CHARACTER(LEN=3) :: atype
 	
 		WRITE(*,*)'Entre com a distancia maxima da funcao g(r) '
 		READ(*,*) maxr
 		WRITE(*,*) 'Entre com a largura (delta r) do bin do histograma'
 		READ(*,*) delr
 		WRITE(*,*) 'Entre com o indice do atomo que se deseja calcular a funcao g(r) as seu redor'
 		READ(*,*) idx
 		WRITE(*,*) 'Entre com o simbolo do atomo que se deseja calcular a funcao g(r)'
 		READ(*,*) atype
 		WRITE(*,*) 'Escreva os parametros a, b e c em angstrom da caixa'
 		READ(*,*) a,b,c
 	
		OPEN(10,FILE='vmd.xyz',STATUS='OLD')
		OPEN(11,FILE='gr.dat',STATUS='REPLACE')
 	 	
 		maxbin= INT(maxr/delr)+1
 		WRITE(*,*) 'Numero de bin no histograma: ', maxbin
 	
 		ALLOCATE(HIST(maxbin),STAT=err)
 		IF(err .NE. 0) PRINT*, 'Impossible allocate memory for HIST array'
 		ALLOCATE(vecgr(maxbin),STAT=err)
 		IF(err .NE. 0) PRINT*, 'Impossible allocate memory for vecgr array'
 		ALLOCATE(vecintgr(maxbin),STAT=err)
 		IF(err .NE. 0) PRINT*, 'Impossible allocate memory for vecintgr array'
 		ALLOCATE(vecr(maxbin),STAT=err)
 		IF(err .NE. 0) PRINT*, 'Impossible allocate memory for vecr array'
 	
 		DO i=1,maxbin
 			HIST(i) = 0.0
 			vecgr(i) = 0.0
 			vecintgr(i) = 0.0
 		END DO
 	
 		vecr(1) = delr/2
 		DO i=2,maxbin
 			vecr(i)=vecr(i-1)+delr
 		END DO
 	
 		!Calculating the gas ideal density rho
 		READ(10,*) tnatm
 		ALLOCATE(mxyz(tnatm,3),STAT=err)
 		IF(err .NE. 0) PRINT*, 'Impossible allocate memory for mxyz array'
 		ALLOCATE(vecatmsbl(tnatm),STAT=err)
 		IF(err .NE. 0) PRINT*, 'Impossible allocate memory for vecatmsbl array'
 	 	
 	cntatm = 0
 	DO i=1,tnatm
 		READ(10,*) vecatmsbl(i), (mxyz(i,j),j=1,3)
 		IF(vecatmsbl(i)== atype) THEN 		
 			cntatm = cntatm + 1
 		END IF
 	END DO
 	WRITE(*,*) 'Existem ',cntatm,' atomos de ', atype,' em cada frame'
 	rho = cntatm/(a*b*c)
 	WRITE(*,*) 'A densidade de ',atype,'tratado como gas ideal e de ',rho,'atmos/A**3'
 	
 	REWIND(UNIT=10)
 	
 	
 	READ(10,*) tnatm 
 	nstep = 0	
 	IO=0
 	DO WHILE(IO == 0)
 		!Reading the frame
 		DO i=1,tnatm
 			READ(10,*) vecatmsbl(i), (mxyz(i,j),j=1,3)
 		END DO
 		
 		DO i=1,tnatm
 			IF(vecatmsbl(i) == atype ) THEN
 				delx = mxyz(idx,1)-mxyz(i,1)
 				dely = mxyz(idx,2)-mxyz(i,2)
 				delz = mxyz(idx,3)-mxyz(i,3)
 				rijsq = delx*delx+dely*dely+delz*delz 
 				rij = SQRT(rijsq)
 				bin = INT(rij/delr)+1
 				IF(bin <= maxbin) THEN
 					HIST(bin) = HIST(bin)+1
 				END IF
 			END IF
 		END DO
 		nstep = nstep + 1
 		READ(10,*,IOSTAT=IO) tnatm		
 			
 	END DO
 	
 	!Normalizing the bins
 	const = 4.0*PI*rho/3.0
 	DO bin = 1,maxbin
 		rlower = REAL(bin-1)*delr
 		rupper = rlower + delr
 		nideal = const*(rupper**3-rlower**3)
 		vecgr(bin) = REAL(HIST(bin))/(REAL(nstep)*nideal)
 	END DO
 	
 	!Calculating the g(r) integral
 	vecintgr(1) = HIST(1)/nstep
 	DO bin = 2,maxbin
 		vecintgr(bin) = vecintgr(bin-1) + HIST(bin)/nstep	
 	END DO
 	
 	DO bin = 1,maxbin
 		WRITE(11,*) vecr(bin), vecgr(bin), vecintgr(bin)
 	END DO
		
	END SUBROUTINE SR010
   