	MODULE SR023MOD
	IMPLICIT NONE
	
	CONTAINS
	SUBROUTINE BoxSS()
	!BoxSS: Box Solvation Solute
	!Box with solute and solvent
	!InFileS : Solute file
	!InFileS1: Solvente 1 file; InFileS2: Solvente 2 file
	!InFileS3: Solvente 3 file
	!nSLVT number of type of solvente to be used.
	!nAtomS, nAtomS1,nAtomS2 number of atoms of the solute and solventes
	!MatGrid is matrix with grid of points where solvente molecules are centered.
	!PBSX,PBSY,PBSZ are periodic box size X,Y, and Z, respectively.

	IMPLICIT NONE
	CHARACTER(LEN=20) :: InFileS,InFileS1,InFileS2,InFileS3
	INTEGER   :: nSLVT,nAtomS,nAtomS1,ERR1,I,J 
	REAL :: LaSLT,LaSLVT1,PBSX,PBSY,PBSZ,mDist
	REAL, ALLOCATABLE :: SLT(:,:),SLVT1(:,:),SLVT2(:,:),SLVT3(:,:), &
						 SYBSLT(:),SYBSLT1(:),SYBSLT2(:),SYBSLT3(:), &
						 MatGrid(:,:)
		
	WRITE(*,*) 'How many type of solvent will be used in solvation box?'
	READ(*,*) nSLVT
	IF(nSLVT == 1) THEN		
		WRITE(*,*) 'Input solute file name'
		READ(*,*) InFileS
		WRITE(*,*) 'Input solvente file name'
		READ(*,*) InFileS1
		OPEN(UNIT=90,FILE=InFileS,ACTION='READ',STATUS='OLD')
		OPEN(UNIT=91,FILE=InFileS1,ACTION='READ',STATUS='OLD')
		OPEN(UNIT=94,FILE='BoxSS.xyz',ACTION='WRITE',STATUS='REPLACE')
		READ(UNIT=90,FMT=*) nAtomS
		READ(UNIT=91,FMT=*) nAtomS1
		ALLOCATE(SLT(nAtomS,3),SLVT1(nAtomS1,3),STAT=ERR1)
		IF (ERR1 /= 0) PRINT*, 'MEMORY CAN NOT BE ALLOCATED'
		DO I=1,nAtomS
			READ(UNIT=90,FMT=*) SYBSLT(I),(SLT(I,J),J=1,3)
		END DO
		DO I=1,nAtomS1
			READ(UNIT=91,FMT=*) SYBSLT1(I),(SLVT1(I,J),J=1,3)
		END DO
		CALL CenterMol(SLT,nAtomS,3)
		CALL CenterMol(SLVT1,nAtomS1,3)
		CALL SmallestBox(SLT,nAtomS,LaSLT)
		WRITE(UNIT=*,FMT=*) 'The smallest box is ', 'a=',LaSLT,'b=',LaSLT,'c=',LaSLT
		WRITE(UNIT=*,FMT=*) 'Input the periodic box size vectors a,b, and c:'
		READ(UNIT=*,FMT=*) PBSX,PBSY,PBSZ
		WRITE(UNIT=*,FMT=*) 'Input the minimum distance between two molecule:'		
		READ(UNIT=*,FMT=*) mDist
		CALL SmallestBox(SLVT1,nAtomS1,LaSLVT1)
		CALL GridXYZ(MatGrid,LaSLT,LaSLVT1,PBSX,PBSY,PBSZ,mDist)
	END IF
	END SUBROUTINE BoxSS
	
	SUBROUTINE CenterMol(ARR1,NL,NC)
	!Subroutine to center molecule on origin
	!ARR1 stands for array
	!NL number of lines
	!NC number of columns
	IMPLICIT NONE
	REAL,DIMENSION(:,:),INTENT(INOUT) :: ARR1
	INTEGER,INTENT(IN) :: NL,NC
	REAL :: SomaX=0,SomaY=0,SomaZ=0
	INTEGER :: I
	DO I=1,NL
		SomaX = SomaX + ARR1(I,1)
		SomaY = SomaY + ARR1(I,2)
		SomaZ = SomaZ + ARR1(I,3)
	END DO
		SomaX = SomaX/REAL(NL)
		SomaY = SomaY/REAL(NL)
		SomaZ = SomaZ/REAL(NL)
	DO I=1,NL
		ARR1(I,1) = ARR1(I,1) - SomaX
		ARR1(I,2) = ARR1(I,2) - SomaY
		ARR1(I,3) = ARR1(I,3) - SomaZ
	END DO
		
	END SUBROUTINE CenterMol
	
	SUBROUTINE SmallestBox(ARR2,nAtom,La)
	!Smalest Cubic box enclosing solute or solvent in angstroms
	!La stands for Lattice a.
	!nAtom stands for number of atoms.
	!Lght stands for Lenght between two atoms inside molecule.
	
	REAL,DIMENSION(:,:),INTENT(IN)::ARR2
	REAL, INTENT(OUT) :: La
	INTEGER, INTENT(IN) :: nAtom
	REAL :: Lght,SQRX,SQRY,SQRZ
	INTEGER :: I,J
	
	La=0
	DO I=1,nAtom
		DO J=I+1,nAtom
			SQRX = (ARR2(I,1)-ARR2(J,1))*(ARR2(I,1)-ARR2(J,1))
			SQRY = (ARR2(I,2)-ARR2(J,2))*(ARR2(I,2)-ARR2(J,2))
			SQRZ = (ARR2(I,3)-ARR2(J,3))*(ARR2(I,3)-ARR2(J,3))
			Lght = SQRT(SQRX + SQRY + SQRZ)
			IF(Lght > La) THEN
				La = Lght
			END IF
		END DO
	END DO
	
	END SUBROUTINE SmallestBox
	
	SUBROUTINE GridXYZ(ARR3,LaSLT,LbSLT,LcSLT, &
				LaSLVT1,LbSLVT1,LcSLVT1,A,B,C,Dist)
	IMPLICIT NONE
	!A,B,C are box solvation size inputed by user
	!RngX range of X to put solvent
	!LbinX length given by LaSLVT1 + Dist
	!Dist minimal distance between two molecule. It should be input by user
	!nbinX number of bins for X coordinate
	
	REAL,ALLOCATABLE,INTENT(OUT) :: ARR3(:,:)
	REAL,INTENT(IN) :: LaSLT,LbSLT,LcSLT,LaSLVT1,LbSLVT1,LcSLVT1,A,B,C,Dist
	REAL :: RngX,RngY,RngZ,LbinX,LbinY,LbinZ,SomaX,SomaY,SomaZ
	INTEGER :: I, nbinX,nbinY,nbinZ
	
	RngX  = (A - LaSLT)/2
	RngY  = (B - LbSLT)/2
	RngZ  = (C - LcSLT)/2
	
	LbinX = LaSLVT1 + Dist
	LbinY = LbSLVT1 + Dist
	LbinZ = LcSLVT1 + Dist
	
	nbinX = FLOOR(RngX/LbinX)
	nbinY = FLOOR(RngY/LbinY)
	nbinZ = FLOOR(RngZ/LbinZ)
	
	SomaX = (LaSLT + LbinX)/2
	SomaY = (LbSLT + LbinY)/2
	SomaZ = (LcSLT + LbinZ)/2
	
	! ALLOCATE(ARR3(nbinX+1,)
	M = 0
	DO I=0,nbinX
		DO J=0,nbinY
			DO K=0,nbinZ
				ARR3(I+M,1) = REAL(I)*SomaX + REAL(I)*LbinX
				ARR3(I+M,2) = REAL(J)*SomaY + REAL(J)*LbinY
				ARR3(I+M,3) = REAL(J)*SomaZ + REAL(K)*LbinZ
				M = J + k +1
			END DO
		END DO
	END DO
	
	END SUBROUTINE GridXYZ
	
	SUBROUTINE ROTX()
	
	END SUBROUTINE ROTX
	
	SUBROUTINE ROTY()
	
	END SUBROUTINE ROTY
	
	SUBROUTINE ROTZ()
	
	END SUBROUTINE ROTZ
	
	SUBROUTINE TRANSL()
	
	END SUBROUTINE TRANSL
	
	
	END MODULE SR023MOD