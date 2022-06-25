   SUBROUTINE SR009()
   IMPLICIT NONE
   !It calculates the radial distribution function with PBC
   
      REAL, PARAMETER :: PI = 3.1415926
      REAL, DIMENSION(:,:), ALLOCATABLE :: cxyz
      REAL, DIMENSION(:),   ALLOCATABLE :: HIST, vecr,vecgr,vecintgr, cx,cy,cz
      CHARACTER(LEN=3),DIMENSION(:), ALLOCATABLE    :: vecatoms, atmsymbl
      REAL :: a,b,c, rijsq,rij,deltar,maxr,rho,delx,dely,delz,rlower,rupper,const,nideal
      CHARACTER(LEN=3) :: symbol, dammy='W',answersave
      INTEGER :: natoms, natomsexcluded, IO, frame,err1,err2, &
                 i,j,k,idx,maxbin,bin,nstep,cntatoms
      INTEGER, DIMENSION(:), ALLOCATABLE :: idxx
        
        
      WRITE(*,'(A)') 'GIVE MAXIMUM r FOR g(r) FUNCTION'
      READ(*,*) maxr
      WRITE(*,'(A)') 'GIVE THE WIDTH OF THE HISTOGRAM BINS'
      READ(*,*) deltar
      WRITE(*,'(A)') 'GIVE ATOM INDEX TO CALCULATE g(r) AROUND IT'
      READ(*,*) idx
      WRITE(*,'(A)') 'GIVE ATOM SYMBOL TO CALCULATE THE g(r) DISTRIBUTION FUNCTION'
      READ(*,*) symbol
      WRITE(*,'(A)') 'GIVE CELL DIMENSION a, b, AND c IN ANGSTROM'
      READ(*,*) a,b,c
      WRITE(*,'(A)') 'GIVE THE NUMBER OF ATOMS TO BE EXCLUDED FROM THE CALCULATION. TYPE ZERO FOR NO ONE'
      READ(*,*) natomsexcluded
      WRITE(*,'(A)') 'DO YOU WISH TO SAVE THE REPLICATED CONFIGURATION FILE? yes/no'
      READ(*,*) answersave
              
      OPEN(90,FILE='vmd.xyz',ACTION='READ', STATUS='OLD')
      OPEN(91,FILE='tmp.dat',ACTION='READWRITE', STATUS='REPLACE')
      OPEN(92,FILE='grpbc.dat',ACTION='WRITE',STATUS='REPLACE')
      OPEN(93,FILE='replicated.xyz',ACTION='WRITE',STATUS='REPLACE')
        
      READ(90,*) natoms
      ALLOCATE(idxx(natomsexcluded))
      DO i=1,natomsexcluded
         WRITE(*,'(A,I3,A)') 'WRITE THE INDEX OF THE ATOM ',i, ' TO BE EXCLUDED'
         READ(*,*) idxx(i)
      END DO

      ALLOCATE(atmsymbl(natoms), STAT=err1)
      IF(err1 /= 0) PRINT*, 'MEMORY CAN NOT BE ALLOCATED FOR atmsybl ARRAY'
      ALLOCATE(cx(natoms),cy(natoms),cz(natoms),STAT=err2)
      IF(err2 /= 0) PRINT*, 'MEMORY CAN NOT BE ALLOCATED FOR cx,cy, AND cz ARRAYS'
 
      IF(natomsexcluded /= 0) THEN
         WRITE(*,*)
         WRITE(*,'(A,I8)') '....EXCLUDING SELECTED ATOMS FROM THE CALCULATION FOR EACH FRAME....'
         WRITE(*,*)
         WRITE(*,'(A)')    '........THIS CAN TAKE WHILE! BE PATIENTE........'
         WRITE(*,*)
      ELSE
         WRITE(*,*)
         WRITE(*,'(A)') '......CALCULATING THE tmp.dat FILE........'
         WRITE(*,*)
         WRITE(*,'(A)') '......THIS CAN TAKE WHILE! BE PATIENTE......'
         WRITE(*,*)
      END IF
      
      IO = 0
      DO WHILE(IO == 0)

         WRITE(91,"(I5)") natoms
         WRITE(91,*)
         DO i=1,natoms
            READ(90,*) atmsymbl(i), cx(i),cy(i),cz(i)
         END DO
         DO i=1,natoms
            DO k=1,natomsexcluded
               IF(idxx(k) == i) THEN
                  WRITE(91,'(A,3F15.5)') dammy,cx(i),cy(i),cz(i)
                  GO TO 9999
               END IF
            END DO
         WRITE(91,'(A,3F15.5)') atmsymbl(i), cx(i), cy(i), cz(i)
  9999   CYCLE
         END DO
         READ(90,*,IOSTAT=IO) natoms
         
      END DO
      CLOSE(UNIT=90)
      REWIND(UNIT=91)
      maxbin= INT(maxr/deltar)+1
 
      ALLOCATE(HIST(maxbin))
      ALLOCATE(vecgr(maxbin))
      ALLOCATE(vecintgr(maxbin))
      ALLOCATE(vecr(maxbin))
      DO i=1,maxbin
         HIST(i) = 0.0
         vecgr(i) = 0.0
         vecintgr(i) = 0.0
      END DO
 
      vecr(1) = deltar/2
      DO i=2,maxbin
         vecr(i)=vecr(i-1)+deltar
      END DO

      !Calculating the gas ideal density rho
      READ(91,*) natoms
      ALLOCATE(cxyz(natoms*27,3))
      ALLOCATE(vecatoms(natoms*27))
         
      cntatoms = 0
      DO i=1,natoms
         READ(91,*) vecatoms(i), (cxyz(i,j),j=1,3)
         IF(vecatoms(i)== symbol) THEN 
            cntatoms = cntatoms + 1
         END IF
      END DO
      WRITE(*,'(A,I5,A,A2,A)') 'THERE ARE ',cntatoms,' ATOMS OF ', symbol,' IN EACH FRAME'
      WRITE(*,*)
      rho = REAL(cntatoms)/(a*b*c)
      WRITE(*,'(A,A2,A,F10.5,A)') 'THE DENSITY OF ',symbol,'CONSIDERED AS IDEAL GAS IS OF',rho,' ATOMS/A**3'
 
      REWIND(UNIT=91)
      READ(91,*) natoms 
      nstep = 0
      IO=0
      WRITE(*,'(A,I7)') 'CALCULATING THE g(r) DISTRIBUTION FUNCTION ... PLEASE WAIT ...'
      WRITE(*,*)

      DO WHILE(IO == 0)
         DO i=1,natoms
            READ(91,*) vecatoms(i), (cxyz(i,j),j=1,3)
         END DO
         !Replicating the frame
         DO i=1,natoms
            vecatoms(natoms+i)=vecatoms(i)
            cxyz(natoms+i,1)= cxyz(i,1)+a
            cxyz(natoms+i,2)= cxyz(i,2)
            cxyz(natoms+i,3)= cxyz(i,3)
 
            vecatoms(2*natoms+i)=vecatoms(i)
            cxyz(2*natoms+i,1)= cxyz(i,1)-a
            cxyz(2*natoms+i,2)= cxyz(i,2)
            cxyz(2*natoms+i,3)= cxyz(i,3)
 
            vecatoms(3*natoms+i)=vecatoms(i)
            cxyz(3*natoms+i,1)= cxyz(i,1)
            cxyz(3*natoms+i,2)= cxyz(i,2)+b
            cxyz(3*natoms+i,3)= cxyz(i,3)
 
            vecatoms(4*natoms+i)=vecatoms(i)
            cxyz(4*natoms+i,1)= cxyz(i,1)
            cxyz(4*natoms+i,2)= cxyz(i,2)-b
            cxyz(4*natoms+i,3)= cxyz(i,3)
 
            vecatoms(5*natoms+i)=vecatoms(i)
            cxyz(5*natoms+i,1)= cxyz(i,1)
            cxyz(5*natoms+i,2)= cxyz(i,2)
            cxyz(5*natoms+i,3)= cxyz(i,3)+c
 
            vecatoms(6*natoms+i)=vecatoms(i)
            cxyz(6*natoms+i,1)= cxyz(i,1)
            cxyz(6*natoms+i,2)= cxyz(i,2)
            cxyz(6*natoms+i,3)= cxyz(i,3)-c
 
            vecatoms(7*natoms+i)=vecatoms(i)
            cxyz(7*natoms+i,1)= cxyz(i,1)+a
            cxyz(7*natoms+i,2)= cxyz(i,2)+b
            cxyz(7*natoms+i,3)= cxyz(i,3)

            vecatoms(8*natoms+i)=vecatoms(i)
            cxyz(8*natoms+i,1)= cxyz(i,1)+a
            cxyz(8*natoms+i,2)= cxyz(i,2)-b
            cxyz(8*natoms+i,3)= cxyz(i,3)
 
            vecatoms(9*natoms+i)=vecatoms(i)
            cxyz(9*natoms+i,1)= cxyz(i,1)-a
            cxyz(9*natoms+i,2)= cxyz(i,2)+b
            cxyz(9*natoms+i,3)= cxyz(i,3)
 
            vecatoms(10*natoms+i)=vecatoms(i)
            cxyz(10*natoms+i,1)= cxyz(i,1)-a
            cxyz(10*natoms+i,2)= cxyz(i,2)-b
            cxyz(10*natoms+i,3)= cxyz(i,3)
 
            vecatoms(11*natoms+i)=vecatoms(i)
            cxyz(11*natoms+i,1)= cxyz(i,1)+a
            cxyz(11*natoms+i,2)= cxyz(i,2)
            cxyz(11*natoms+i,3)= cxyz(i,3)+c
 
            vecatoms(12*natoms+i)=vecatoms(i)
            cxyz(12*natoms+i,1)= cxyz(i,1)-a
            cxyz(12*natoms+i,2)= cxyz(i,2)
            cxyz(12*natoms+i,3)= cxyz(i,3)+c
 
            vecatoms(13*natoms+i)=vecatoms(i)
            cxyz(13*natoms+i,1)= cxyz(i,1)+a
            cxyz(13*natoms+i,2)= cxyz(i,2)
            cxyz(13*natoms+i,3)= cxyz(i,3)-c
 
            vecatoms(14*natoms+i)=vecatoms(i)
            cxyz(14*natoms+i,1)= cxyz(i,1)-a
            cxyz(14*natoms+i,2)= cxyz(i,2)
            cxyz(14*natoms+i,3)= cxyz(i,3)-c
 
            vecatoms(15*natoms+i)=vecatoms(i)
            cxyz(15*natoms+i,1)= cxyz(i,1)
            cxyz(15*natoms+i,2)= cxyz(i,2)+b
            cxyz(15*natoms+i,3)= cxyz(i,3)+c
 
            vecatoms(16*natoms+i)=vecatoms(i)
            cxyz(16*natoms+i,1)= cxyz(i,1)
            cxyz(16*natoms+i,2)= cxyz(i,2)-b
            cxyz(16*natoms+i,3)= cxyz(i,3)+c
 
            vecatoms(17*natoms+i)=vecatoms(i)
            cxyz(17*natoms+i,1)= cxyz(i,1)
            cxyz(17*natoms+i,2)= cxyz(i,2)+b
            cxyz(17*natoms+i,3)= cxyz(i,3)-c
 
            vecatoms(18*natoms+i)=vecatoms(i)
            cxyz(18*natoms+i,1)= cxyz(i,1)
            cxyz(18*natoms+i,2)= cxyz(i,2)-b
            cxyz(18*natoms+i,3)= cxyz(i,3)-c
 
            vecatoms(19*natoms+i)=vecatoms(i)
            cxyz(19*natoms+i,1)= cxyz(i,1)-a
            cxyz(19*natoms+i,2)= cxyz(i,2)+b
            cxyz(19*natoms+i,3)= cxyz(i,3)+c
 
            vecatoms(20*natoms+i)=vecatoms(i)
            cxyz(20*natoms+i,1)= cxyz(i,1)+a
            cxyz(20*natoms+i,2)= cxyz(i,2)+b
            cxyz(20*natoms+i,3)= cxyz(i,3)+c
 
            vecatoms(21*natoms+i)=vecatoms(i)
            cxyz(21*natoms+i,1)= cxyz(i,1)-a
            cxyz(21*natoms+i,2)= cxyz(i,2)-b
            cxyz(21*natoms+i,3)= cxyz(i,3)+c
 
            vecatoms(22*natoms+i)=vecatoms(i)
            cxyz(22*natoms+i,1)= cxyz(i,1)+a
            cxyz(22*natoms+i,2)= cxyz(i,2)-b
            cxyz(22*natoms+i,3)= cxyz(i,3)+c
 
            vecatoms(23*natoms+i)=vecatoms(i)
            cxyz(23*natoms+i,1)= cxyz(i,1)-a
            cxyz(23*natoms+i,2)= cxyz(i,2)+b
            cxyz(23*natoms+i,3)= cxyz(i,3)-c
 
            vecatoms(24*natoms+i)=vecatoms(i)
            cxyz(24*natoms+i,1)= cxyz(i,1)+a
            cxyz(24*natoms+i,2)= cxyz(i,2)+b
            cxyz(24*natoms+i,3)= cxyz(i,3)-c
 
            vecatoms(25*natoms+i)=vecatoms(i)
            cxyz(25*natoms+i,1)= cxyz(i,1)-a
            cxyz(25*natoms+i,2)= cxyz(i,2)-b
            cxyz(25*natoms+i,3)= cxyz(i,3)-c
 
            vecatoms(26*natoms+i)=vecatoms(i)
            cxyz(26*natoms+i,1)= cxyz(i,1)+a
            cxyz(26*natoms+i,2)= cxyz(i,2)-b
            cxyz(26*natoms+i,3)= cxyz(i,3)-c          
         END DO
      DO i=1,27*natoms
         IF(vecatoms(i) == symbol ) THEN
            delx = cxyz(idx,1)-cxyz(i,1)
            dely = cxyz(idx,2)-cxyz(i,2)
            delz = cxyz(idx,3)-cxyz(i,3)
            rijsq = delx*delx+dely*dely+delz*delz 
            rij = SQRT(rijsq)
            bin = INT(rij/deltar)+1
            IF(bin <= maxbin) THEN
               HIST(bin) = HIST(bin)+1
            END IF
         END IF
      END DO
      IF(answersave == 'yes') THEN
        WRITE(93,'(I4)') 27*natoms
        WRITE(93,*)
        DO i=1,27*natoms
            WRITE(93,'(A2,3F15.7)') vecatoms(i),(cxyz(i,j),j=1,3)    
        END DO
      END IF
           
      nstep = nstep + 1

      READ(91,*,IOSTAT=IO) natoms
      END DO
 
        !Normalizing the bins
        const = 4.0*PI*rho/3.0
        DO bin = 1,maxbin
                rlower = REAL(bin-1)*deltar
                rupper = rlower + deltar
                nideal = const*(rupper**3-rlower**3)
                vecgr(bin) = REAL(HIST(bin))/(REAL(nstep)*nideal)
        END DO
 
        !Calculating the g(r) integral
        vecintgr(1) = HIST(1)/nstep
        DO bin = 2,maxbin
                vecintgr(bin) = vecintgr(bin-1) + HIST(bin)/nstep
        END DO
 
        DO bin = 1,maxbin
                WRITE(92,*) vecr(bin), vecgr(bin), vecintgr(bin)
        END DO
      DEALLOCATE(cxyz)
      DEALLOCATE(vecatoms)
      DEALLOCATE(HIST)
      DEALLOCATE(vecgr)
      DEALLOCATE(vecintgr)
      DEALLOCATE(vecr)
      DEALLOCATE(cx,cy,cz)
      DEALLOCATE(atmsymbl)
      DEALLOCATE(idxx)
      CLOSE(UNIT=91)
      CLOSE(UNIT=92)
      CLOSE(UNIT=93)     

   END SUBROUTINE SR009
