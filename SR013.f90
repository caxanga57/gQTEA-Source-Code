   SUBROUTINE SR013()
   IMPLICIT NONE
   !This subroutine generates evp.dat file to be opened using origin
      
   INTEGER :: IO,isave,dammy,step
   CHARACTER(LEN=30) :: inputevp
   REAL :: nfi,ekinc,temph,tempp,etot,enthal,econs,econt,dt
   
   WRITE(*,*) 'GIVE INPUT FILE NAME THAT HAS EXTENSION .evp'
   READ(*,*)  inputevp
   WRITE(*,*) 'GIVE SAMPLING INTERVAL (This is iprint keyword in cp.x input)'
   READ(*,*)  isave
   WRITE(*,*) 'GIVE dt (time step) IN atu'
   READ(*,*) dt
   
   OPEN(10,FILE=inputevp,ACTION='READ', STATUS='OLD')
   OPEN(11,FILE='evp.dat',ACTION='WRITE',STATUS='REPLACE')
   
   WRITE(*,*) 'Short Legend and Physical Units in the Output from left to right '
   WRITE(*,*) '--------------------------------------------- '
   WRITE(*,*) 'NFI     [int]          - nfi will be given in ps '
   WRITE(*,*) 'EKINC   [HARTREE A.U.] - kinetic energy of the fictitious electronic dynamics'
   WRITE(*,*) 'TEMPH   [K]            - Temperature of the fictitious cell dynamics'
   WRITE(*,*) 'TEMP    [K]            - Ionic temperature'
   WRITE(*,*) 'ETOT    [HARTREE A.U.] - Scf total energy (Kohn-Sham hamiltonian)'
   WRITE(*,*) 'ENTHAL  [HARTREE A.U.] - Enthalpy ( ETOT + P * V )'
   WRITE(*,*) 'ECONS   [HARTREE A.U.] - Enthalpy + kinetic energy of ions and cell'
   WRITE(*,*) 'ECONT   [HARTREE A.U.] - Constant of motion for the CP lagrangian'
   WRITE(*,*) 'EKIONIC [HARTREE A.U.] - ECONS - ENTAHL '
   
   IO=0
   READ(10,*) nfi,ekinc,temph,tempp,etot,enthal,econs,econt
   nfi = REAL(isave)
   DO WHILE(IO == 0)
      WRITE(11,'(F14.4,F15.5,2F8.1,5F15.5)') nfi*dt*0.02418884326505/1000.0, &
                                         ekinc,temph, tempp, etot,enthal, econs, econt, econs-enthal
      nfi = nfi + isave
      READ(10,*,IOSTAT=IO) dammy,ekinc,temph,tempp,etot,enthal,econs,econt
   END DO
   END SUBROUTINE SR013

	   