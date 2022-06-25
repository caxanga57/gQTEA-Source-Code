      SUBROUTINE SR021()
C
C     Programa desenvolvido por Patricia RP Barreto-2002/2007 e
C     redimenssionado por Valter Carvalho-2013/2015 para c†lculo das correá‰es
C     de tunelamento de Bell 35 e 58 e implementaá∆o da formulaá∆o d-TST
C     Ultima modificacao - 10 de MarÁo de 2016
C     com o calculo da taxa unimolecular no limite de alta pressao
C     calculo da MEP em unidades de mi^1/2 a0
                                    
C
C     Programa para:
C     1. calcular as taxas de reacoes usando TST, com coeficiente de transmissao de
C        1.1 Wigner
C        1.2 Eckart
C        1.3 Tsallis
C        1.4 Bell de 1935 e 1958
C     2. determinar a MEP
C     3 determinar os parametros A,n e Ea de Ahrrenius Modificada
C
      implicit real*8(a-h,o-z)
        real(8) temp(500),ental(500),tk(100)
	real(8) vag(500),vg(500),fac(500)
	real(8) beta1(500), beta2(500),xfim(500)
	real(8) tipo(500,500),sigma(500),veff(500)
	real(8) xm(500,500),xmn(500),xmp2(500)
	real(8) xi(500,500),yi(500,500),zi(500,500)
	real(8) xcm(500),ycm(500),zcm(500),pia(500)
	real(8) fvib(500,500),deg(500,500)
	real(8) stt(500,500),et(500,500)
	real(8) sr(500,500),er(500,500)
	real(8) sv(500,500),ev(500,500)
	real(8) sto(500,500),eto(500,500)
	real(8) zep(500),xif(500),et298(500)
	real(8) cpt(500,500),cpr(500,500)
	real(8) cpv(500,500),cpto(500,500)
	real(8) xmab(500),xmts(500)
	real(8) xmp2rc(500),xmp2pc(500)
	real(8) xmp2ab(500),xmp2ts(500)
	real(8) xifab(500),xifts(500)
	real(8) zpeab(500),zpets(500)
	real(8) piab(500),pits(500)
	real(8) qtab(100,500),eatab(100,500)
	real(8) qtts(100,500),eatts(100,500)
	real(8) qt(100,500),eat(100,500)
	real(8) qr(100,500),ear(100,500)
	real(8) qv(100,500),eav(100,500)
	real(8) qtt(100,500),eatt(100,500),xb1kk(100,500)
	real(8) a(100,500),ea(100,500),xb1kapa(500),xb22tkapa(500)
	real(8) xk(100,500),dxk(100,500),xkk(100,500),dxkapa(100,500)
	real(8) xkk2(100,500),arh(100,100),xkapa(500),dtsxkapa(100,500)
	real(8) coefa(10,10), coefb(10),dexpxkapa(100,500),xbkk(100,500)
	real(8) coef(500,500),xmepi(500),xbkapa(500),dtsxk(100,500)
	real(8) caa(10,10), cbb(10),dexpxk(100,500),xb22tkk(100,500)
	real(8) aa(500,500),bb(500),xx(500)
        real(8) xarh(500),af(500,500),daf(500,500),dtsaf(500,500)
C        real(8) xkappa,xii,xff,xkappa2,xpp,xmp,tk2,bbu,cc,aau
        real(8) xii,xff,xpp,xmp,tk2,bbu,cc,aau,xkappa2(500)
	real(8) s(210),vmep(200,210),vagg(200,210)
	real(8) nti(2), ntf(2)
	real(8) f0,f2n,sm,sp,a1,a2,dt,ter1,ter2,ter3
	integer na(500),nv(500,500),nesp(100)
	character*10 espab(500),espts(500),esp(500)
	character*10 reg1(500),reg2(500),ts(500)
	character*10 prod1(500),prod2(500) 
        character*15 arq1,arq2,arq3
	character*16 arq4,arq5
	character*17 arq6,arq7,arq8,arq9,arq10,arq11
	character*20 arq15
        character larq1*2,larq2*2
	character*10 sai(500),sailn(500)
	
C     PERMEABILIDADE DE BARREIRA DE ECKART

      xkappa(tr1,tr2,delt,x,tk2,pi)=1.0d0-(dcosh(2.0d0*pi*(tr1))
     +           +dcosh(2.0d0*pi*delt))/(dcosh(2.0d0*pi*(tr2))
     +           +dcosh(2.0d0*pi*delt))
	alf(x)=a1t*dsqrt(x)
	bt(x,aa1,aa2)=b1t*dsqrt((1.0d0+x)*aa1-aa2)
	dtt(aa1,aa2)=(1.0d0/pi)*dsqrt(abs(aa1*aa2-pi*pi/4.0d0))

C     DEFINICAO DAS CONSTANTES
C
      xkb=1.3806d-23                                                      CONTANTE DE BOLTZMANN   J/K
	h=6.6257d-34                                                      CONTANTE DE PLANCK
	c=2.99793d10                                                      VELOCIDADE DA LUZ
	press=1.013d+05
	xmc=1.6605402d-27                                                 UNIDADE DE MASSA ATOMICA
	xnav=6.0224d26
	xnav2=6.0224d23                                                   CONSTANTE DE AVOGADRO
	xme=9.108d-31                                                     MASSA DE ELETRON
	r=1.9872d0                                                        CONSTANTE DOS GASES    CAL.K-1.MOL-1
	rj=4186.0d0                                                       TRANSFORMANDO PARA  JOULE
	pi=4.0d0*atan(1.0)
	hb=h/(2.0*pi)                                                     h CORTADO

C
C     IMPRESSAO DE INFORMACOES INICIAIS
C

	write(*,*)"       "
	write(*,*)"       "
	write(*,*)"Programa desenvolvido por Patricia Barreto-2002/2007"
	write(*,*)"e redimenssionado por Valter Carvalho -2013/2015"
	write(*,*)" Ultima modificacao - Outubro de 2015     "
	write(*,*)"       "
	write(*,*)"       "
      write(*,*)"Programa para: "
	write(*,*)"1. calcular as taxas de reacoes usando TST, com: "
	write(*,*)"   1.1 coeficiente de transmissao de Wigner"
	write(*,*)"   1.2 coeficiente de transmissao de Eckart"
	write(*,*)"   1.3 coeficiente de transmissao de Bell-1935/58"
	write(*,*)"   1.4 coeficiente de transmissao de Tsallis"
	write(*,*)"2.determinar a MEP       "
	write(*,*)"3.determinar os parametros A,n e Ea de Ahrrenius "
	write(*,*)"       "
	write(*,*)"       "
C
C     ENTRADA DE DADOS: QUANTIDADE DE ESPECIES E NOME DO ARQUIVO
C
	write(*,*)"       "
	write(*,*)" Entrada de dados      "
	write(*,*)"       "
        write(*,*)"Numero de especies do arquivo reagente"
	read(*,*) nab
	write(*,*)"       "
	write(*,*)"Arquivo que contem os dados para os reagentes (AB)"
	read(*,'(a2)') larq1
	write(*,*)"       "
        write(*,*)"Numero de especies do arquivo TS"
	read(*,*) nts
	write(*,*)"       "
        write(*,*)"arquivo que contem os dados para os TS (TS)"
	read(*,'(a2)') larq2
        nr=nts
	arq1="especie"//larq1//".dat"                                     ARQUIVO DE ENTRADA DE DADOS DOS REAGENTES
        arq2="especie"//larq2//".dat"                                     ARQUIVO DE ENTRADA DE DADOS DA TS
	arq3="reacao"//larq1//larq2//".dat"                               ARQUIVO DE ENTRADA DE DADOS DE AMBOS
C
      open(1,file=arq1,status='old')
      open(2,file=arq2,status='old')
      open(3,file=arq3,status='old')
C
      write(*,*)"      "
      write(*,*)"Iniciando os calculos"
      write(*,*)"      "

C
C     DEFINICAO DA FAIXA DE TEMPERATURA
C
	nt=65
        nti(1)=2
	ntf(1)=15
        nti(2)=15
	ntf(2)=65
	
C     PRIMEIRA FAIXA DE TEMPERATURA
	do i=1,3
         xin=float(i-1)
	   temp(i)=100.0d0*xin

	enddo

C     SEGUNDA FAIXA DE TEMPERATURA
	   temp(4)=250.0d0
	   temp(5)=298.15d0

	   
C     TERCEIRA FAIXA DE TEMPERATURA
	do i=6,10
           xin=float(i-6)
	   temp(i)=300.0d0+50.0d0*xin

	enddo
	
C     QUARTA FAIXA DE TEMPERATURA
	do i=11,65
           xin=float(i-11)
  	   temp(i)=600.0d0+100.0d0*xin
  	enddo

C
C     CALCULO DAS FUNCOES DE PARTICAO REAGENTES E PRODUTOS
C
C     ENTRADA DE DADOS
C
        nesp(1)=nab                                                       NUMERO DE ESPECIES DOS REAGENTES
	nesp(2)=nts                                                       NUMERO DE ESPECIES DA TS
        ient=1
	do i=1,nab
C     TIPO=>µTOMO:1,MOLECULA LINEAR:2, MOLECULA N∆O-LINEAR:3; SIGMA=>N£MERO DE SIMETRIA:TABELADO(PODE USAR 1)
		read(1,*)esp(i),na(i),tipo(1,i),sigma(i)
		
		do j=1,na(i)
C     LEITURA DA GEOMETRIA DE CADA REAGENTE
			read(1,*)xm(j,i),xi(j,i),yi(j,i),zi(j,i)
		enddo
		
		read(1,*)nv(1,i)
		zep(i)=0.0d0
		
		do j=1,nv(1,i)
C     LEITURA DAS FREQ E DEGENERESCENCIA
			read(1,*)fvib(j,i),deg(j,i)
C     ENERGIA DE ZPVE PARA CADA FREQ
			t2ev=0.5d0*xnav*deg(j,i)*h*c*fvib(j,i)
C     SOMA TOTAL DAS ENERGIAS DE ZPVE PARA FREQ EM KILOCALORIAS
			zep(i)=zep(i)+1.0d-03*t2ev/4.184d0
		enddo	
		read(1,*)xmp2kcal,xif(i)
C	        xmp2(i)=xmp2kcal/627.5095                                 ENERGIA DA MOLêCULA EM KCALORIAS
C     ENERGIA DA MOLêCULA EM HARTREE
                xmp2(i)=xmp2kcal
	enddo
	close(1)

C
C     CALCULO DAS FUNCOES DE PARTICAO DE TRANSLACAO
C
        
      do i=1,nesp(1)
	
	  sum=0.0d0
	  
	    do ia=1,na(i)
C     CALCULO MA MASSA MOLECULAR PARA CADA MOLECULA
	      sum=sum+xm(ia,i)
	    enddo
C     MASSA MOLECULAR TOTAL DA MOLECULA
	  xmn(i)=sum
C     MASSA MOLECULAR EM KG
	  xmass=xmn(i)*xmc
C     FUNÄAO DE TRANSLAÄ«O PARA AS 65=nt TEMPERATURAS
          do k=1,nt
C     CALCULO DO VOLUME
            t1=r*temp(k)/press
C     DEFINIÄÂES CONFUSAS
	    t2=2*pi*xmass*xkb*temp(k)/(h*h)
	    qt(k,i)=(2.0*pi*xmc*xmn(i)*xkb*temp(k)/(h*h))**(3./2.)
 	    if(xmn(i).eq.0) qt(k,i)=1d6*xnav
	  enddo
C
C     CALCULO DO CENTRO DE MASSA
C
        sumx=0.0d0
	sumy=0.0d0
	sumz=0.0d0
	
	do ia=1,na(i)
	sumx=sumx+xm(ia,i)*xi(ia,i)
	sumy=sumy+xm(ia,i)*yi(ia,i)
	sumz=sumz+xm(ia,i)*zi(ia,i)
	enddo
	
	xcm(i)=sumx/xmn(i)
	ycm(i)=sumy/xmn(i)
	zcm(i)=sumz/xmn(i)
C
C     CALCULO DOS MOMENTOS DE INERCIA EM RELAÄ«O AO CENTRO DE MASSA
C
	sumx=0.0d0
	sumy=0.0d0
	sumz=0.0d0
	sumx2=0.0d0
	sumy2=0.0d0
	sumz2=0.0d0
	sumxy=0.0d0
	sumxz=0.0d0
	sumyz=0.0d0
	sumxyq=0.0d0
	sumxzq=0.0d0
	sumyzq=0.0d0
	
	 do ia=1,na(i)
	   sumx=sumx+xm(ia,i)*(xi(ia,i))
	   sumy=sumy+xm(ia,i)*(yi(ia,i))
	   sumz=sumz+xm(ia,i)*(zi(ia,i))
	   sumx2=sumx2+xm(ia,i)*(xi(ia,i)-xcm(i))**2
	   sumy2=sumy2+xm(ia,i)*(yi(ia,i)-ycm(i))**2
	   sumz2=sumz2+xm(ia,i)*(zi(ia,i)-zcm(i))**2
	   
	   if(sumx2.ne.0.0d0) sumfim=sumx2
	   if(sumy2.ne.0.0d0) sumfim=sumy2
 	   if(sumz2.ne.0.0d0) sumfim=sumz2
 	   
	   sumxy=sumxy+xm(ia,i)*(xi(ia,i))*(yi(ia,i))
	   sumxz=sumxz+xm(ia,i)*(xi(ia,i))*(zi(ia,i))
	   sumyz=sumyz+xm(ia,i)*(yi(ia,i))*(zi(ia,i))
	   sumxyq=sumxyq+xm(ia,i)*((xi(ia,i))**2.0+(yi(ia,i))**2.0)
	   sumxzq=sumxzq+xm(ia,i)*((xi(ia,i))**2.0+(zi(ia,i))**2.0)
	   sumyzq=sumyzq+xm(ia,i)*((yi(ia,i))**2.0+(zi(ia,i))**2.0)
	   
	enddo
C     COMPONENTES DO MOMENTO DE INERCIA
	xa=sumyzq-(sumy**2.0+sumz**2.0)/xmn(i)
	xb=sumxzq-(sumx**2.0+sumz**2.0)/xmn(i)
	xc=sumxyq-(sumx**2.0+sumy**2.0)/xmn(i)
	xd=sumxy-sumx*sumy/xmn(i)
	xe=sumxz-sumx*sumz/xmn(i)
	xf=sumyz-sumy*sumz/xmn(i)
C     A ESCOLHA DO TIPO IMPACTA NO CALCULO DA FUNÄ«O DE PARTIÄ«O ROTACIONAL
C     PIA S«O OS MOMENTOS DE INERCIA
        if(tipo(1,i).eq.1) pia(i)=0.0d0
        if(tipo(1,i).eq.2) pia(i)=sumfim
	if(tipo(1,i).eq.3) pia(i)=(xa*xb*xc-xa*xf*xf-xc*xd*xd-2*xd*xe*xf-
     =                         xb*xe*xe)
C
C     CALCULO DAS FUNCOES DE PARTICAO DE ROTACAO
C

       xfp=xmc*1.0d-20

       do k=1,nt

C
C     FUNÄ«O DE PARTIÄ«O PARA µTOMOS
C
        if(tipo(1,i).eq.1) then
         qr(k,i)=1.0d0
	  ear(k,i)=0.0d0
C
C     FUNÄ«O DE PARTIÄ«O PARA ESPECIES LINEARES
C
	else 
	 if(tipo(1,i).eq.2) then
	 t1=8*pi*pi*pia(i)*xfp*xkb*temp(k)/(sigma(i)*h*h)
	 qr(k,i)=1.66d-47*((8*pi*pi*pia(i)*xkb*temp(k))/(sigma(i)*h*h))    
C
C     FUNÄ«O DE PARTIÄ«O PARA ESPECIES N«O-LINEARES
C
	else
 	 qr(k,i)=(1.66d-47)*(sqrt(pi)/sigma(i))*(8*pi*pi*pia(i)*
     = 				          xkb*temp(k)/h**2)**(3/2)
	endif
C
	endif

	if(xmn(i).eq.0) qr(k,i)=1
C
C     CALCULO DAS FUNCOES DE PARTICAO DE VIBRACAO
C
	 qv(k,i)=1
         fac(k)=h*c/(xkb*temp(k))
         
	 do j=1,nv(1,i)
	    x=(fac(k)*fvib(j,i))
            qv(k,i)=qv(k,i)*(1-exp(-fac(k)*fvib(j,i)))**(-deg(j,i))
	enddo

	if(xmn(i).eq.0) qv(k,i)=1

C
C     CALCULO FINAL DAS FUNCOES DE PARTICAO
C
	qtt(k,i)=qt(k,i)*qr(k,i)*qv(k,i)

       enddo
      
      enddo
      
      do i=1,nab
	  espab(i)=esp(i)
	  xmab(i)=xmn(i)
	  piab(i)=pia(i)
	  zpeab(i)=zep(i)
	  xmp2ab(i)=xmp2(i)
	  xifab(i)=xif(i)
	  iti=1
	  
	  do it=3,6
	   tk(iti)=temp(it)
	   qtab(iti,i)=qtt(it,i)
	   iti=iti+1
	  enddo
	  
	  do it=8,10,2
	   qtab(iti,i)=qtt(it,i)
	   tk(iti)=temp(it)
	   iti=iti+1
	  enddo
	  
	  do it=11,45,2
	   qtab(iti,i)=qtt(it,i)
	   tk(iti)=temp(it)
	   iti=iti+1
	  enddo
	  
      enddo
C
C     CALCULO DAS FUNCOES DE PARTICAO PARA AS ESTRUTURAS DE TRANSICAO
C
C     ENTRADA DE DADOS
C
      ient=2

	do i=1,nts
          read(2,*)esp(i),na(i),tipo(2,i),sigma(i)
C         write(*,*)esp(i),na(i),tipo(i),sigma(i)
          espts(i)=esp(i)
          
	  do j=1,na(i)
	    read(2,*)xm(j,i),xi(j,i),yi(j,i),zi(j,i)
	  enddo
	read(2,*)nv(2,i)
	  
        zep(i)=0.0d0
        
	do j=1,nv(2,i)
	  read(2,*)fvib(j,i),deg(j,i)
	  t2ev=0.5d0*xnav*deg(j,i)*h*c*fvib(j,i)
	  zep(i)=zep(i)+1.0d-03*t2ev/4.184d0
	enddo
        	
	read(2,*)xmp2kcal,xif(i)
C	xmp2(i)=xmp2kcal/627.5095
	xmp2(i)=xmp2kcal
	enddo
	close(2)
C
C      CALCULO DA ENTALPIA A 298K  TS
C

      do i=1,nesp(1)
	  sum=0.0d0
	  
	  do ia=1,na(i)
	    sum=sum+xm(ia,i)
	  enddo
	  
	  xmn(i)=sum
	  xmass=xmn(i)*xmc
	  
          do k=1,nt
            t1=r*temp(k)/press
	    t2=2*pi*xmass*xkb*temp(k)/(h*h)
	    qt(k,i)=(2.0*pi*xmc*xmn(i)*xkb*temp(k)/(h*h))**(3./2.)
	    if(xmn(i).eq.0) qt(k,i)=1d6*xnav
	  enddo
C
C     CALCULO DO CENTRO DE MASSA TS
C
        sumx=0.0d0
	sumy=0.0d0
	sumz=0.0d0
	
	do ia=1,na(i)
	  sumx=sumx+xm(ia,i)*xi(ia,i)
	  sumy=sumy+xm(ia,i)*yi(ia,i)
	  sumz=sumz+xm(ia,i)*zi(ia,i)
	enddo
	
	xcm(i)=sumx/xmn(i)
	ycm(i)=sumy/xmn(i)
	zcm(i)=sumz/xmn(i)
c	
c     CALCULO DOS MOMENTOS DE INERCIA
c
	sumx=0.0d0
	sumy=0.0d0
	sumz=0.0d0
	sumx2=0.0d0
	sumy2=0.0d0
	sumz2=0.0d0
	sumxy=0.0d0
	sumxz=0.0d0
	sumyz=0.0d0
	sumxyq=0.0d0
	sumxzq=0.0d0
	sumyzq=0.0d0
	
	do ia=1,na(i)
	  sumx=sumx+xm(ia,i)*(xi(ia,i))
	  sumy=sumy+xm(ia,i)*(yi(ia,i))
	  sumz=sumz+xm(ia,i)*(zi(ia,i))
	  sumx2=sumx2+xm(ia,i)*(xi(ia,i)-xcm(i))**2
	  sumy2=sumy2+xm(ia,i)*(yi(ia,i)-ycm(i))**2
	  sumz2=sumz2+xm(ia,i)*(zi(ia,i)-zcm(i))**2
	  
	  if(sumx2.ne.0.0d0) sumfim=sumx2
  	  if(sumy2.ne.0.0d0) sumfim=sumy2
	  if(sumz2.ne.0.0d0) sumfim=sumz2
	  
	  sumxy=sumxy+xm(ia,i)*(xi(ia,i))*(yi(ia,i))
	  sumxz=sumxz+xm(ia,i)*(xi(ia,i))*(zi(ia,i))
	  sumyz=sumyz+xm(ia,i)*(yi(ia,i))*(zi(ia,i))
	  sumxyq=sumxyq+xm(ia,i)*((xi(ia,i))**2.0+(yi(ia,i))**2.0)
	  sumxzq=sumxzq+xm(ia,i)*((xi(ia,i))**2.0+(zi(ia,i))**2.0)
	  sumyzq=sumyzq+xm(ia,i)*((yi(ia,i))**2.0+(zi(ia,i))**2.0)
	enddo
	
	xa=sumyzq-(sumy**2.0+sumz**2.0)/xmn(i)
	xb=sumxzq-(sumx**2.0+sumz**2.0)/xmn(i)
	xc=sumxyq-(sumx**2.0+sumy**2.0)/xmn(i)
	xd=sumxy-sumx*sumy/xmn(i)
	xe=sumxz-sumx*sumz/xmn(i)
	xf=sumyz-sumy*sumz/xmn(i)
        if(tipo(2,i).eq.1) pia(i)=0.0d0
        if(tipo(2,i).eq.2) pia(i)=sumfim
	if(tipo(2,i).eq.3) pia(i)=(xa*xb*xc-xa*xf*xf-xc*xd*xd-2*xd*xe*xf-
     =                         xb*xe*xe)
C
C     CALCULO DAS FUNCOES DE PARTICAO DE ROTACAO TS
C

        xfp=xmc*1.0d-20

        do k=1,nt

C
C     FUNÄ«O DE PARTIÄ«O PARA µTOMOS     TS
C
      if(tipo(2,i).eq.1) then
       qr(k,i)=1.0d0
C
C     FUN«√O DE PARTI«√O PARA ESPECIES LINEARES      TS
C
	else 
	if(tipo(2,i).eq.2) then
	 t1=8*pi*pi*pia(i)*xfp*xkb*temp(k)/(sigma(i)*h*h)
	 qr(k,i)=1.66d-47*((8*pi*pi*pia(i)*xkb*temp(k))/(sigma(i)*h*h))    
C
C     FUN«√O DE PARTI«√O PARA ESPECIES N«O-LINEARES  TS
C
	else
	 x3=xfp**3
	 t1=sqrt(pi*pia(i)*x3)/sigma(i)
	 t2=8*pi*pi*xkb*temp(k)/(h*h)
 	 qr(k,i)=(1.66d-47)*(sqrt(pi)/sigma(i))*(8*pi*pi*pia(i)*
     = 				          xkb*temp(k)/h**2)**(3/2)
	endif
C
	endif

	if(xmn(i).eq.0) qr(k,i)=1
C
C     CALCULO DA FUN«√O PARTI«√O DE VIBRACAO TS
C
	qv(k,i)=1
        fac(k)=h*c/(xkb*temp(k))
        
	do j=1,nv(2,i)
	 x=(fac(k)*fvib(j,i))
         qv(k,i)=qv(k,i)*(1-exp(-fac(k)*fvib(j,i)))**(-deg(j,i))
	enddo

	if(xmn(i).eq.0) eav(k,i)=0
	if(xmn(i).eq.0) qv(k,i)=1

C
C     CALCULO FINAL DAS FUNCOES DE PARTICAO TS
C
	qtt(k,i)=qt(k,i)*qr(k,i)*qv(k,i)
	
       enddo
      enddo

      do i=1,nts
   	espts(i)=esp(i)
	xmts(i)=xmn(i)
	pits(i)=pia(i)
	zpets(i)=zep(i)
	xmp2ts(i)=xmp2(i)
C     TOMANDO O VALOR ABSOLUTO DA FREQU“NCIA NEGATIVA
	xifts(i)=abs(xif(i))
	iti=1
C     REDEFINDINDO OS LIMITES DE TEMPERATURA
	do it=3,6
	  tk(iti)=temp(it)
	  qtts(iti,i)=qtt(it,i)
	  iti=iti+1
	enddo
	
	do it=8,10,2
	  qtts(iti,i)=qtt(it,i)
	  tk(iti)=temp(it)
	  iti=iti+1
	enddo
	
	do it=11,45,2
	 qtts(iti,i)=qtt(it,i)
	 tk(iti)=temp(it)
	 iti=iti+1
	enddo
	
      enddo

      write(*,*)"     "
      write(*,*)"     "
      write(*,*)"    Para reagentes e produtos   "
      write(*,*)" calculo das funcoes de particao: EFETUADO"
      write(*,*)"     "
      write(*,*)"    Para as estruturas de transicao   "
      write(*,*)" calculo das funcoes de particao: EFETUADO"
      write(*,*)"     "
      write(*,*)"     "
      write(*,*)" Inciando proxima etapa dos calculos    "
      write(*,*)"     "
      write(*,*)"  Calculo da taxa de reacao e da MEP"
      write(*,*)"  Ajuste da taxa de reacao na forma de Ahrrenius"
      write(*,*)"     "
      
	do i=1,nr
	  read(3,*)reg1(i),reg2(i),ts(i),prod1(i),prod2(i),xmp2rc(i),
     +         xmp2pc(i),xmepi(i),xfim(i)
	enddo
	close(3)

        nt=24
	do it=1,nt
	 temp(it)=tk(it)
	enddo

	arq10="reacao"//larq1//larq2//"zep.dat"

        open(98,file=arq10,status='unknown')
C
C     CALCULO DA TAXA DE REACAO EM FUNCAO DA TEMPERATURA
C
	do i=1,nr
C     REINDEXANDO OS NOMES DA TS, REAGENTE E PRODUTOS
	  do j=1,nts
	   if (ts(i).eq.espts(j)) n1=j
	  enddo
	  
	  do j=1,nab
	   if (reg1(i).eq.espab(j)) n2=j
	  enddo
	  
	  do j=1,nab
	   if (reg2(i).eq.espab(j)) n3=j
	  enddo
	  
	  do j=1,nab
	   if (prod1(i).eq.espab(j)) n4=j
	  enddo
	  
	  do j=1,nab
  	    if (prod2(i).eq.espab(j)) n5=j
	  enddo

          do k=1,nt
C     -------------PARA REAÄÂES UNIMOLECULARES -------------------
	    if(reg2(i).eq.'ZERO     ') then
C	    write(*,*)'dentro do IF   ',reg2(i)
C     RAZAO DAS FUNLÄOES DE PARTIÄ«O
	      qttg=qtts(k,n1)/qtab(k,n2)
C     DIFERENÄA DOS ZPVE DA TS PELO REAGENTE
              dzpe=xfim(i)*(zpets(n1)-zpeab(n2))
C     ALTURA DA BARREIRA EM KKCAL/MOL
	      ener=(xmp2ts(n1)-xmp2ab(n2))*627.5095d3
C     ALTURA DA BARREIRA COM CORREÄ«O DE ZPVE
	      ea(k,i)=ener+dzpe
C	      ea(k,i)=ener
              write(101,*) xmp2ts(n1), xmp2ab(n2), ea(k,i)
C	    af(k,i)=(1.0d6*xnav*1.0d-03*xkb*temp(k)/h)*qttg
C     FATOR PRE-EXPONENCIAL DA TAXA PARA TST
	      af(k,i)=(1.0d6*1.0d-03*xkb*temp(k)/h)*qttg
C     CALCULO DO ANGULO DE SKEW
              beta1(i)=0.0
              beta2(i)=0.0
C	write(99,*)temp(k),qtab(k,n2),qtab(k,n3),qtts(k,n1),qttg,af(k,i)
C     PREVIS«O DO d
              d=-(1.0d0/3.0d0)*(c*h*xifts(i)*xnav2/(2.0d0*ea(k,i)
     =         *(rj/1000.0d0)))**2

C     PREVIS«O DA TEMPERATURA DE CROSS-OVER
              Tc=(c*h*xifts(i))/(2.0d0*pi*xkb)

C     PREVIS«O DA TEMPERATURA DE VALIDADE DA d-TST

              Td=Tc + d*ea(k,i)*(rj/1000.0d0)/((2.0d0*xkb)*xnav2)
	     else
C	write(*,*)'fora do IF   '
C     -------------PARA REAÄÂES BIMOLECULARES-------------------
C     RAZAO DAS FUNLÄOES DE PARTIÄ«O
  	      qttg=qtts(k,n1)/(qtab(k,n2)*qtab(k,n3))
C     DIFERENÄA DOS ZPVE DA TS PELO REAGENTE
	      dzpe=xfim(i)*(zpets(n1)-(zpeab(n2)+zpeab(n3)))
C	    dzpe2=xfim(i)*(zpets(n1)-(zpeab(n2)+zpeab(n3)))
C     ALTURA DA BARREIRA EM KCAL/MOL
	      ener=(xmp2ts(n1)-(xmp2ab(n2)+xmp2ab(n3)))*627.5095d3
C     ALTURA DA BARREIRA COM CORREÄ«O DE ZPVE EM CAL/MOL
	      ea(k,i)=ener+dzpe
C	      ea(k,i)=ener
C	      write(101,*) ener,dzpe,ea(k,i)
C     PREVIS«O DO d

              d=-(1.0d0/3.0d0)*(c*h*xifts(i)*xnav2/(2.0d0*ea(k,i)
     =         *(rj/1000.0d0)))**2

C     PREVIS«O DA TEMPERATURA DE CROSS-OVER
              Tc=(c*h*xifts(i))/(2.0d0*pi*xkb)
              
C     PREVIS«O DA TEMPERATURA DE VALIDADE DA d-TST

              Td=Tc + d*ea(k,i)*(rj/1000.0d0)/((2.0d0*xkb)*xnav2)

C     COMPLEMENTO DE ROTAÄ«O VIA DISTRIBUIÄ«O TSALLIS TST
C               drotn2=(1.0d0/(2.0d0-d))**(tipo(1,n2)/2.0d0)
C              drotn3=(1.0d0/(2.0d0-d))**(tipo(1,n3)/2.0d0)
C               drotTS=(1.0d0/(2.0d0-d))**(tipo(2,n1)/2.0d0)

C     COMPLEMENTO DE VIBRAÄ«O VIA DISTRIBUIÄ«O TSALLIS TST
C     REAGENTE 1
C             dvibn2=0.0d0
C              if(nv(1,n2).ne. 0) then
C                do l=1,nv(1,n2)
C                   dvibn2=dvibn2+(1.0d0/(1+l*d))
C                enddo
C              else
C                 dvibn2=1.0d0
C              endif
C     REAGENTE 2
C              dvibn3=0.0d0
C              if(nv(1,n3).ne. 0) then
C                do l=1,nv(1,n3)
C                   dvibn3=dvibn3+(1.0d0/(1+l*d))
C                enddo
C              else
C                 dvibn3=1.0d0
C              endif
C     TS
C              dvibts=0.0d0
C              do l=1,nv(2,n1)
C                 dvibts=dvibts+(1.0d0/(1+l*d))
C              enddo
C     COMPLEMENTO DE TRANSLAÄ«O VIA DISTRIBUIÄ«O TSALLIS TST - FUNÄ«O GAMMA
C             zz=(-1.0d0/d)+(-3.0d0/2.0d0)
C              zz1=(-1.0d0/d)
C              gammaln=0.0d0
C              gammaln1=0.0d0
C              do l=1,1000000
C
C                  gammaln=gammaln+dlog((((1.0d0+1.0d0/l)**zz)/
C     =             (1.0d0+zz/l)))
C                  gammaln1=gammaln1+dlog((((1.0d0+1.0d0/l)**zz1)/
C     =             (1.0d0+zz1/l)))
C
C              enddo
C               dtrans=(-1.0d0/d)**(3.0d0/2.0d0)*
C     =          ((1.0d0/zz)/(1.0d0/zz1))*dexp(gammaln-gammaln1)
C                 write(101,*) k,temp(k),ea(k,i),d,d1,zz,zz1,gamma
C	    ea(k,i)=ener+dzpe+(eatts(k,n1)-(eatab(k,n2)+eatab(k,n3)))
C     FATOR PRE-EXPONENCIAL DA TAXA PARA TST
	      af(k,i)=(1.0d6*xnav*1.0d-03*xkb*temp(k)/h)*qttg

C     RAZAO DAS FUNLÄOES DE PARTIÄ«O VIA TSALLIS COM TODAS AS ENTIDADES
C  	      dqttg=drotTS*dvibts*qtts(k,n1)/(drotn2*dvibn2*qtab(k,n2)
C     =           *drotn3*dvibn3*dtrans*qtab(k,n3))
C
C     RAZAO DAS FUNLÄOES DE PARTIÄ«O VIA TSALLIS CONSIDERANDO SOMENTE TS
C  	      dtsqttg=drotTS*dvibts*qtts(k,n1)/(qtab(k,n2)
C     =           *qtab(k,n3))
C     FATOR PRE-EXPONENCIAL DA TAXA PARA TST VIA TSALLIS  COM TODAS AS ENTIDADES
C	      daf(k,i)=(1.0d6*xnav*1.0d-03*xkb*temp(k)/h)*dqttg
C     FATOR PRE-EXPONENCIAL DA TAXA PARA TST VIA TSALLIS  CONSIDERANDO SOMENTE TS
C	      dtsaf(k,i)=(1.0d6*xnav*1.0d-03*xkb*temp(k)/h)*dtsqttg
C     TERMO PARA O CALCULO DO ANGULO DE SKEW 1
              xb=(xmab(n4)*xmab(n3))/(xmab(n2)*xmab(n5))
C     CALCULO DO ANGULO DE SKEW 1
	      beta1(i)=(180./pi)*acos(sqrt(xb))
C     TERMO PARA O CALCULO DO ANGULO DE SKEW 2
	      xb=(xmab(n2)-xmab(n4))*xmts(n1)/(xmab(n3)*xmab(n4))
C     CALCULO DO ANGULO DE SKEW 2
	      beta2(i)=(180./pi)*atan(sqrt(xb))
	
C	write(99,*)temp(k),qtab(k,n2),qtab(k,n3),qtts(k,n1),qttg,af(k,i)
	     endif
C     ------TAXA DA TERORIA DO ESTADO DE TRANSIÄ«O---------------------
	     xk(k,i)=af(k,i)*exp(-ea(k,i)/(r*temp(k)))
C     ------TAXA DA TERORIA DO ESTADO DE TRANSIÄ«O VIA TSALLIS COM TODAS AS ENTIDADES---------
C	     dxk(k,i)=daf(k,i)*(1.0d0-d*ea(k,i)/(r*temp(k)))
C     =        **((1.0+d)/d)
C             dxkapa(k,i)=  dxk(k,i)/xk(k,i)
C
C     ------TAXA DA TERORIA DO ESTADO DE TRANSIÄ«O VIA TSALLIS CONSIDERANDO SOMENTE TS---------
C	dtsxk(k,i)=dtsaf(k,i)*(1.0d0-d*ea(k,i)/(r*temp(k)))**((1.0+d)/d)
C
C             dtsxkapa(k,i)=  dtsxk(k,i)/xk(k,i)
C     ------TAXA DA TERORIA DO ESTADO DE TRANSIÄ«O VIA TSALLIS CONSIDERANDO SOMENTE EULER---------
	dexpxk(k,i)=af(k,i)*(1.0d0-d*ea(k,i)/(r*temp(k)))**((1.0d0)/d)


                   dexpxkapa(k,i)=  dexpxk(k,i)/xk(k,i)
C     ------CORREÄ«O DE TUNELAMENTO DE WIGNER--------------------------
  	     xkapa(k)=1.0+(1./24.)*(h*xifts(n1)*c/(xkb*temp(k)))**2
C     ------CORREÄ«O DE TUNELAMENTO DE BELL 1935--------------------------
             xb1kapa(k)=  ((xkb*temp(k))- (hb)*xifts(n1)*c
     =       *exp(ea(k,i)/(r*temp(k)))
     =      *exp(-(ea(k,i)*(rj/1000.0d0))/(c*(hb)*xifts(i)*xnav2)))
     =       /((xkb*temp(k))- (hb)*xifts(n1)*c)
C     ------CORREÄ«O DE TUNELAMENTO DE BELL 1958--------------------------
  	     xbkapa(k)=1.0/2.0*(h*xifts(n1)*c/(xkb*temp(k)))
     =       /sin(1.0/2.0*(h*xifts(n1)*c/(xkb*temp(k))))
C     ------CORREÄ«O DE TUNELAMENTO DE BELL 1958 COM 2 TERMOS------------
  	     xb22tkapa(k)=1.0/2.0*(h*xifts(n1)*c/(xkb*temp(k)))
     =       /sin(1.0/2.0*(h*xifts(n1)*c/(xkb*temp(k))))
     =       -(exp(ea(k,i)/(r*temp(k)))
     =       *exp(-(2.*pi*ea(k,i)*(rj/1000.0d0))/(c*h*xifts(i)*xnav2)))
     =       /((2.*pi*xkb*temp(k)/(h*xifts(n1)*c))-1.)
C     -----TAXA TST COM CORREÄ«O DE WIGNER-----------------------------
	     xkk(k,i)=xk(k,i)*xkapa(k)
C     -----TAXA TST COM CORREÄ«O DE BELL 1958-----------------------------
	     xbkk(k,i)=xk(k,i)*xbkapa(k)
C     -----TAXA TST COM CORREÄ«O DE BELL 1935-----------------------------
	     xb1kk(k,i)=xk(k,i)*xb1kapa(k)
C     -----TAXA TST COM CORREÄ«O DE BELL 1958 COM 2 TERMOS----------------
	     xb22tkk(k,i)=xk(k,i)*xb22tkapa(k)
C     ENTALPIA DE REAÄ«O COM CORREÄ«O DE ZPVE EM KCALORIAS/MOL
	     ental(i)=1.0d-3*(((xmp2ab(n4)+xmp2ab(n5))-
     *          (xmp2ab(n2)+xmp2ab(n3)))*627.5095d3+
     *          xfim(i)*(zpeab(n4)+zpeab(n5))-(zpeab(n2)+zpeab(n3)))
     

C     ALTURA DA BARREIRA COM CORREÄ«O DE ZPVE
      vg(i)=(xmp2ts(n1)-(xmp2ab(n2)+xmp2ab(n3)))*627.5095
      vag(i)=vg(i)+xfim(i)*(zpets(n1)-(zpeab(n2)+zpeab(n3)))*1d-3

C	write(95,*) temp(k),xk(k,i),xkk(k,i),xkapa
C	pause
C
	  enddo
	enddo

C
C     DETERMINACAO DA MEP
C
	do i=1,nr
	
	   do j=1,nts
	    if (ts(i).eq.espts(j)) n1=j
	   enddo
	   
	   do j=1,nab
	    if (reg1(i).eq.espab(j)) n2=j
	   enddo
	   
	   do j=1,nab
	    if (reg2(i).eq.espab(j)) n3=j
	   enddo
	   
	   do j=1,nab
	    if (prod1(i).eq.espab(j)) n4=j
	   enddo
	   
	   do j=1,nab
	    if (prod2(i).eq.espab(j)) n5=j
           enddo
C     ENERGIA DOS REAGENTES E DOS PRODUTOS
          xmp2r=xmp2ab(n2)+xmp2ab(n3)
	  xmp2p=xmp2ab(n4)+xmp2ab(n5)

	write(98,*)'reag= ',xmp2r
	write(98,*)'RC= ',xmp2rc(n1)
	write(98,*)'TS= ',xmp2ts(n1)
	write(98,*)'PC= ',xmp2pc(n1)
	write(98,*)'prod= ',xmp2p

	if(vg(i).gt.0) then
	aau=((xmp2ab(n4)+xmp2ab(n5))-
     *    (xmp2ab(n2)+xmp2ab(n3)))*627.5095
	aac=((xmp2ab(n4)+xmp2ab(n5))-
     *    (xmp2ab(n2)+xmp2ab(n3)))*627.5095+xfim(i)*
     *    ((zpeab(n4)+zpeab(n5))-(zpeab(n2)+zpeab(n3)))*1.0d-3
	vdu=(xmp2ts(n1)-(xmp2ab(n2)+xmp2ab(n3)))*627.5095
	vdc=(xmp2ts(n1)-(xmp2ab(n2)+xmp2ab(n3)))*627.5095+
     *    xfim(i)*((zpets(n1)-(zpeab(n2)+zpeab(n3))))*1.0d-3   
	veff(i)=vdc
	 else
	if(xmp2r.le.xmp2rc(n1)) then
	   aau=((xmp2ab(n4)+xmp2ab(n5))-
     *       (xmp2ab(n2)+xmp2ab(n3)))*627.5095
	   aac=((xmp2ab(n4)+xmp2ab(n5))-
     *       (xmp2ab(n2)+xmp2ab(n3)))*627.5095+xfim(i)*
     *       ((zpeab(n4)+zpeab(n5))-(zpeab(n2)+zpeab(n3)))*1.0d-3
	   vdu=(xmp2ts(n1)-(xmp2ab(n2)+xmp2ab(n3)))*627.5095
	   vdc=(xmp2ts(n1)-(xmp2ab(n2)+xmp2ab(n3)))*627.5095+
     *       xfim(i)*((zpets(n1)-(zpeab(n2)+zpeab(n3))))*1.0d-3 
         veff(i)=vdc  
	 else
	   if(xmp2p.le.xmp2ts(n1)) then
	      aau=((xmp2ab(n4)+xmp2ab(n5))-xmp2rc(n1))*627.5095
	      aac=((xmp2ab(n4)+xmp2ab(n5))-xmp2rc(n1))*627.5095+xfim(i)*
     *          ((zpeab(n4)+zpeab(n5))-(zpeab(n2)+zpeab(n3)))*1.0d-3
	      vdu=(xmp2ts(n1)-xmp2rc(n1))*627.5095
	      vdc=(xmp2ts(n1)-xmp2rc(n1))*627.5095+
     *          xfim(i)*((zpets(n1)-(zpeab(n2)+zpeab(n3))))*1.0d-3  
            veff(i)=vdc   
	  else
	      aau=(xmp2pc(n1)-xmp2rc(n1))*627.5095
	      aac=(xmp2pc(n1)-xmp2rc(n1))*627.5095+xfim(i)*
     *          ((zpeab(n4)+zpeab(n5))-(zpeab(n2)+zpeab(n3)))*1.0d-3
	      vdu=(xmp2ts(n1)-xmp2rc(n1))*627.5095
	      vdc=(xmp2ts(n1)-xmp2rc(n1))*627.5095+
     *          xfim(i)*((zpets(n1)-(zpeab(n2)+zpeab(n3))))*1.0d-3  
            veff(i)=vdc   
        endif
      endif
	endif
	
	write(*,*)ts(i),vdu,aau,vdu-aau,veff(i)
c	write(*,*)ts(i),vdc,aac,vdc-aac,veff(i)
c	pause

	  bbu=(2.*vdu-aau)+2.*sqrt(vdu*(vdu-aau))
	  bbc=(2.*vdc-aac)+2.*sqrt(vdc*(vdc-aac))

	xmia=xmab(n2)*xmab(n3)/xmts(n1)
      xmib=1.0d0
	fca2a=627.5095d0*xmc*(2.0d0*pi*c*0.529177249d-10)**2
     *      /4.36d-18
	fca2b=(627.5095d0*xmc/4.36d-18)*(2.0d-10*pi*c)**2
	fca2c=(rj/xnav)*(c/(2*pi))**2. 

	xmi=xmib
	fca2=fca2b

	term1=fca2a*xmia*bbu*xifts(n1)**2./(2.*vdu*(vdu-aau))
	term2=fca2b*xmia*bbu*xifts(n1)**2./(2.*vdu*(vdu-aau))
	term3=fca2c*xmia*bbu*xifts(n1)**2./(2.*vdu*(vdu-aau))
	term4=fca2a*xmib*bbu*xifts(n1)**2./(2.*vdu*(vdu-aau))
	term5=fca2b*xmib*bbu*xifts(n1)**2./(2.*vdu*(vdu-aau))
	term6=fca2c*xmib*bbu*xifts(n1)**2./(2.*vdu*(vdu-aau))

      alfa1=sqrt(abs(term1))
      alfa2=sqrt(abs(term2))
      alfa3=sqrt(abs(term3))
      alfa4=sqrt(abs(term4))
      alfa5=sqrt(abs(term5))
      alfa6=sqrt(abs(term6))

	sou1=-(1/alfa1)*(log((aau+bbu)/(abs(bbu-aau))))
	sou2=-(1/alfa2)*(log((aau+bbu)/(abs(bbu-aau))))
	sou3=-(1/alfa3)*(log((aau+bbu)/(abs(bbu-aau))))
	sou4=-(1/alfa4)*(log((aau+bbu)/(abs(bbu-aau))))
	sou5=-(1/alfa5)*(log((aau+bbu)/(abs(bbu-aau))))
	sou6=-(1/alfa6)*(log((aau+bbu)/(abs(bbu-aau))))

	soc1=-(1/alfa1)*(log((aac+bbc)/(abs(bbc-aac))))
	soc2=-(1/alfa2)*(log((aac+bbc)/(abs(bbc-aac))))
	soc3=-(1/alfa3)*(log((aac+bbc)/(abs(bbc-aac))))
	soc4=-(1/alfa4)*(log((aac+bbc)/(abs(bbc-aac))))
	soc5=-(1/alfa5)*(log((aac+bbc)/(abs(bbc-aac))))
	soc6=-(1/alfa6)*(log((aac+bbc)/(abs(bbc-aac))))
   
	xlixo=xmi*bbu*xifts(n1)**2./(2.*vdu*(vdu-aau))
	alfa=sqrt(abs(fca2*xmi*bbu*xifts(n1)**2./(2.*vdu*(vdu-aau))))
	sou=-(1/alfa)*(log((aau+bbu)/(abs(bbu-aau))))
	soc=-(1/alfa)*(log((aac+bbc)/(abs(bbc-aac))))
	cc=xfim(i)*(zpeab(n2)+zpeab(n3))*1.0d-3
	
c	write(*,*)ts(i),vdu,vdc,aau,bbu,sou
c	write(*,*)aac,bbc,soc,alfa,cc
	write(98,*)reg1(i),xfim(i)*zpeab(n2)*1.0e-3
	write(98,*)reg2(i),xfim(i)*zpeab(n3)*1.0e-3
	write(98,*)prod1(i),xfim(i)*zpeab(n4)*1.0e-3
      write(98,*)prod2(i),xfim(i)*zpeab(n5)*1.0e-3
      write(98,*)ts(i),xfim(i)*zpets(n1)*1.0e-3
	dzep1=xfim(i)*(zpets(n1)-(zpeab(n2)+zpeab(n3)))*1.0e-3
	dzep2=xfim(i)*((zpeab(n4)+zpeab(n5))-(zpeab(n2)+zpeab(n3)))*
     &  1.0d-3
	write(98,*)'dzpe(TS-R)=  ',dzep1
	write(98,*)'dzpe(P-R)=  ',dzep2

	do is=1,201
	fis=float(is)
	s(is)=-xmepi(i)*(1.-2.*(fis-1.)/200.0)
	yu=exp(alfa*(s(is)-sou))
	yc=exp(alfa*(s(is)-soc))
	vmep(i,is)=aau*yu/(1+yu)+(bbu*yu/((1+yu)**2))
	vagg(i,is)=aac*yc/(1+yc)+(bbc*yc/((1+yc)**2))+cc

	yu1=exp(alfa1*(s(is)-sou1))
	yu2=exp(alfa2*(s(is)-sou2))
	yu3=exp(alfa3*(s(is)-sou3))
	yu4=exp(alfa4*(s(is)-sou4))
	yu5=exp(alfa5*(s(is)-sou5))
	yu6=exp(alfa6*(s(is)-sou6))

	yc1=exp(alfa1*(s(is)-soc1))
	yc2=exp(alfa2*(s(is)-soc2))
	yc3=exp(alfa3*(s(is)-soc3))
	yc4=exp(alfa4*(s(is)-soc4))
	yc5=exp(alfa5*(s(is)-soc5))
	yc6=exp(alfa6*(s(is)-soc6))
	
	vmep1=aau*yu1/(1+yu1)+(bbu*yu1/((1+yu1)**2))
	vmep2=aau*yu2/(1+yu2)+(bbu*yu2/((1+yu2)**2))
	vmep3=aau*yu3/(1+yu3)+(bbu*yu3/((1+yu3)**2))
	vmep4=aau*yu4/(1+yu4)+(bbu*yu4/((1+yu4)**2))
	vmep5=aau*yu5/(1+yu5)+(bbu*yu5/((1+yu5)**2))
	vmep6=aau*yu6/(1+yu6)+(bbu*yu6/((1+yu6)**2))
	
	vag1=aac*yc1/(1+yc1)+(bbc*yc1/((1+yc1)**2))+cc
	vag2=aac*yc2/(1+yc2)+(bbc*yc2/((1+yc2)**2))+cc
	vag3=aac*yc3/(1+yc3)+(bbc*yc3/((1+yc3)**2))+cc
	vag4=aac*yc4/(1+yc4)+(bbc*yc4/((1+yc4)**2))+cc
	vag5=aac*yc5/(1+yc5)+(bbc*yc5/((1+yc5)**2))+cc
	vag6=aac*yc6/(1+yc6)+(bbc*yc6/((1+yc6)**2))+cc


c	write(97,*)is,s(is),vmep(i,is),vagg(i,is)

c	write(99,'(14(3x, e11.5))')s(is),vmep1,vmep2,vmep3,vmep4,vmep5,
c     =                 vmep6,vag1,vag2,vag3,vag4,vag5,vag6
	enddo

c	write(*,*)'MEP calculada'
c	write(*,*)'iniciar calculo de kappa'

c	pause

c
c     integral numerica de kappa(T) de Eckart
c

      fca1a2=2*pi*4.36d-18/(627.5095*h*c)
	a1=vdc*fca1a2/(xifts(n1))
	a2=(vdc-aac)*fca1a2/(xifts(n1))

	if (a1.gt.a2) then 
	 vmax=vdc
	 att1=a1
	 att2=a2
	else
	 vmax=vdc-aac
	 att1=a2
	 att2=a1
      endif

      a1=abs(att1)
	a2=abs(att2)

	vint=vmax*4.36d-18/(627.5095d0*xkb)
c	vint=4.36d-18/(627.5095d0*xkb)

c	ter3=dtt(a1,a2)
	a1t=(1.0d0/pi)*dsqrt(a1)*(1.0d0/dsqrt(a1)+1.0d0/dsqrt(a2))**(-1.0)
	b1t=(1.0d0/pi)*(1.0d0/dsqrt(a1)+1.0d0/dsqrt(a2))**(-1.0)

c      ter3t=0.5d0*dsqrt((bbc-alfa)/alfa)

c	write(*,*)'ter3   ',ts(i),ter3,ter3t

      xmax=(a1**3*a2**2 - 2*a1**2.5*a2**2.5 + a1**2*a2**3 - 
     -    25538*a1**2*a2*Pi**2 + 25538*a1*a2**2*Pi**2 + 
     -    163047361*a1*Pi**4 + 
     -    326094722*Sqrt(a1)*Sqrt(a2)*Pi**4 + 
     -    163047361*a2*Pi**4)/(51076.*a1**2*a2*Pi**2)

      xmin=(a2**3*a1**2 - 2*a2**2.5*a1**2.5 + a2**2*a1**3 - 
     -    25538*a2**2*a1*Pi**2 + 25538*a2*a1**2*Pi**2 + 
     -    163047361*a2*Pi**4 + 
     -    326094722*Sqrt(a2)*Sqrt(a1)*Pi**4 + 
     -    163047361*a1*Pi**4)/(51076.*a2**2*a1*Pi**2)
    
      xff=xmax
c      write(*,*)a1,a2,vmax,xmax,xmin
c      write(*,*)'ter3= ',ter3
      
	do 2222 it=1,nt

	tk2=temp(it)
c      a1=abs(att1)*fca1a2*xkb*
c	a2=abs(att2)*fca1a2/(xifts(n1))

	ter3=dtt(a1,a2)
c	write(*,*)'ter3= ',ter3

      if(ter3.le.113.0d0) then
c      goto 1111

	nint=5000
	xii=0.0d0
c      xff=xfim(i)
	xinc=float(nint)
	ter1=alf(xii)
	ter2=bt(xii,a1,a2)

c	ter1=0.5d0*dsqrt(xii/alfa)
c	ter2=0.5d0*dsqrt(abs(xii-aac)/alfa)

	ttest1=ter1+ter2
	ttest2=(ter1-ter2)

	tert1=0.5d0*dsqrt(xii/alfa)+0.5d0*dsqrt(abs(xii-aac)/alfa)
	tert1=0.5d0*dsqrt(xii/alfa)-0.5d0*dsqrt(abs(xii-aac)/alfa)

c	write(*,*)ttest1,tert1,ttest2,tert2
c	write(*,*)'f2n',ter1,ter2,ttest1,ttest2,ter3
c	pause

      f0=xkappa(ttest2,ttest1,ter3,xii,tk2,pi)*dexp(-xii*vint/tk2)

	ter1=alf(xff)
	ter2=bt(xff,a1,a2)

c	ter1=0.5d0*dsqrt(xff/alfa)
c	ter2=0.5d0*dsqrt(abs(xff-aac)/alfa)

	ttest1=(ter1+ter2)
	ttest2=abs(ter1-ter2)

	tert1=0.5d0*dsqrt(xii/alfa)+0.5d0*dsqrt(abs(xii-aac)/alfa)
	tert1=0.5d0*dsqrt(xii/alfa)-0.5d0*dsqrt(abs(xii-aac)/alfa)

c	write(*,*)ttest1,tert1,ttest2,tert2
c	write(*,*)'f2n',ter1,ter2,ttest1,ttest2,ter3
c	pause
	f2n=xkappa(ttest2,ttest1,ter3,xff,tk2,pi)*exp(-xff*vint/tk2)

	sm=0.d0
	sp=0.d0
	smk=0.0d0
	spk=0.0d0
	xint1=0.d0
	xhp=(xff-xii)/(2.0d0*xinc)

	do iin=1,2*nint,2
	xmp=(-xii+float(iin)*xhp)
	ter1=alf(xmp)
	ter2=bt(xmp,a1,a2)

c	ter1=0.5d0*dsqrt(xmp/alfa)
c	ter2=0.5d0*dsqrt(abs(xmp-aac)/alfa)

	ttest1=(ter1+ter2)
	ttest2=(ter1-ter2)
	
c	write(97,*)'sm',ter1,ter2,ttest1,ttest2,ter3

	sm=sm+xkappa(ttest2,ttest1,ter3,xmp,tk2,pi)*exp(-xmp*vint/tk2)
	smkk=xkappa(ttest2,ttest1,ter3,xmp,tk2,pi)
	smk=smk+xkappa(ttest2,ttest1,ter3,xmp,tk2,pi)
	xint1=4.0d0*sm
	enddo
	xint3=0.0D0
	do iin=2,2*nint-1,2
	xint2=xint3+xint1
	xpp=(-xii+float(iin)*xhp)
	ter1=alf(xpp)
	ter2=bt(xpp,a1,a2)
	
c	ter1=0.5d0*dsqrt(xpp/alfa)
c	ter2=0.5d0*dsqrt(abs(xpp-aac)/alfa)

	ttest1=(ter1+ter2)
	ttest2=(ter1-ter2)

c	write(97,*)'sp',ter1,ter2,ttest1,ttest2,ter3

  	sp=sp+xkappa(ttest2,ttest1,ter3,xmp,tk2,pi)*exp(-xpp*vint/tk2)
	spk=spk+xkappa(ttest2,ttest1,ter3,xmp,tk2,pi)
	spkk=xkappa(ttest2,ttest1,ter3,xmp,tk2,pi)
	xint1=2.0d0*sp

	enddo
	
	xint=(xhp/3.0d0)*((f0+4.0d0*sm+2.0d0*sp+f2n))
      xkappa2(it)=exp(vint/tk2)*xint*vint/tk2
C     -----TAXA TST COM CORREÄ«O DE ECKART-----------------------------
      xkk2(it,i)=xkappa2(it)*xk(it,i)
      
c	write(97,*)ts(i),tk2,temp(it),xk(it,i),xkk(it,i),xkk2(it,i),
c     +           xkk(it,i)/xk(it,i),xkappa2
c      write(*,*)tk2,xkk(it,i)/xk(it,i),xkappa2

      else

      xkappa2=0.0d0
      xkk2(it,i)=xk(it,i)
c	write(97,*)ts(i),tk2,temp(it),xk(it,i),xkk(it,i),xkk2(it,i),
c     +           xkk(it,i)/xk(it,i),xkappa2
c      write(*,*)tk2,xkk(it,i)/xk(it,i),xkappa2

	endif

c	write(*,*)ts(i),tk2,f0,4*sm,2*sp,f2n,xint,xkappa2
	
c	goto 2222

c1111 xkk2(it,i)=xk(it,i)
c      goto 2222

c      write(*,*)tk2,xkk(it,i)/xk(it,i),xkappa2
	
 2222 continue
      
	enddo
	
      write(*,*)"segunda parte do calculo efetuada"
      write(*,*)"     "
c
c     impressao parcial dos resultados
c
      sai(1)="reg1"
	sai(2)="reg2"
	sai(3)="TS"
	sai(4)="V(tst)"
	sai(5)="VaG"
	sai(6)="Veff"
	sai(7)="Ental"
	sai(8)="beta1"
	sai(9)="beta2"
c
	arq4="reacao"//larq1//larq2//"ke.dat"
      open(4,file=arq4,status='unknown')

      write(4,'(9(3x, A6))')(sai(ii),ii=1,9)
	do i=1,nr
      write(4,'(3(3x, A6),6(3x, e11.5))')reg1(i),reg2(i),ts(i),vg(i),
     *                          vag(i),veff(i),ental(i),beta1(i),
     *                           beta2(i)    
	enddo
      close(4)
	close(10)

C     -----IMPRESS«O DOS DADOS CINêTICOS  de k e lnk
	sai(1)="Temp"
	sai(2)="1000/Temp"
	sai(3)="k(T)"
	sai(4)="kWg*k(T))"
        sai(5)="kEck*k(T)"
        sai(6)="kBell-1958"
        sai(7)="kdTST"
        sai(8)="kBell-1935"
        sai(9)="kBell-1958-2Ter"
        
       	sailn(1)="Temp"
	sailn(2)="1000/Temp"
	sailn(3)="lnk(T)"
	sailn(4)="lnkWg*k"
        sailn(5)="lnkEck*k"
        sailn(6)="lnkBell-58"
        sailn(7)="lnkdTST"
        sailn(8)="lnkBell35"
        sailn(9)="lnkBel582T"
        arq6="reacao"//larq1//larq2//"ktr.dat"
        arq11="reacao"//larq1//larq2//"lnk.dat"
      open(6,file=arq6,status='unknown')
      open(25,file=arq11,status='unknown')
C      write(6,*)(ts(i),i=1,nr)
C      write(25,*)(ts(i),i=1,nr)
      write(6,'(9(3x, A12))')(sai(ii),ii=1,9)
      write(6,*)  "  k x T "
      write(25,'(9(2x, A12))')(sailn(ii),ii=1,9)
      write(25,*) " lnk x T "
	do k=1,nt
      write(6,'(60(4x, e11.5))')temp(k),1.0d3/temp(k),(xk(k,i),i=1,nr),
     +                           (xkk(k,i),i=1,nr),(xkk2(k,i),i=1,nr),
     +                           (xbkk(k,i),i=1,nr),
     +                           (dexpxk(k,i),i=1,nr),
     +                         (xb1kk(k,i),i=1,nr),(xb22tkk(k,i),i=1,nr)
C       write(*,*)   (xkk2(k,i),i=1,nr)
      write(25,'(50(4x, e11.5))')temp(k),1.0d3/temp(k),
     +     (dlog(xk(k,i)),i=1,nr),(dlog(xkk(k,i)),i=1,nr),
     +  (dlog(xkk2(k,i)),i=1,nr),(dlog(xbkk(k,i)),i=1,nr),
     +  (dlog(dexpxk(k,i)),i=1,nr),
     +  (dlog(xb1kk(k,i)),i=1,nr),
     +  (dlog(xb22tkk(k,i)),i=1,nr)
     

	enddo
	close(6)
	close(25)

      sai(1)="s"
	sai(2)="Vmep"
	sai(3)="VaG"
	arq9="reacao"//larq1//larq2//"v_n.dat"
      open(9,file=arq9,status='unknown')
      write(9,*)(ts(i),i=1,nr)
      write(9,'(5(3x, A6))')(sai(ii),ii=1,3)
	do is=1,201
	write(9,'(101(3x, e11.5))')s(is),(vmep(ir,is),ir=1,nr),
     +                            (vagg(ir,is),ir=1,nr)
	enddo
	close(9)

c	
c     Ajuste de Ahrrenius
c
      slnt=0.
	s1t=0.
	st=0.
	stlnt=0.
	slnt2=0.
	slnt1t=0.
	do it=1,nt
	slnt=slnt+dlog(temp(it))
	s1t=s1t+1./temp(it)
	st=st+temp(it)
	stlnt=stlnt+temp(it)*dlog(temp(it))
	slnt2=slnt2+(dlog(temp(it)))**2.
	slnt1t=slnt1t+(dlog(temp(it)))/temp(it)
	enddo
	coefa(1,1)=st
	coefa(1,2)=stlnt
	coefa(1,3)=-nt
	coefa(2,1)=nt
	coefa(2,2)=slnt
	coefa(2,3)=-s1t
	coefa(3,1)=slnt
	coefa(3,2)=slnt2
	coefa(3,3)=-slnt1t
      do ir=1,nr
	do ic=1,3
	slnk=0.
	stlnk=0.
	slntk=0.
	serr=0.
	do it=1,nt
       if(ic.eq.1) xarh(it)=xk(it,ir)
	 if(ic.eq.2) xarh(it)=xkk(it,ir)
	 if(ic.eq.3) xarh(it)=xkk2(it,ir)
      slnk=slnk+dlog(xarh(it))
	stlnk=stlnk+temp(it)*dlog(xarh(it))
	slntk=slntk+dlog(xarh(it))*dlog(temp(it))
	enddo
	coefb(1)=stlnk
      coefb(2)=slnk
	coefb(3)=slntk
	caa(1,1)=coefa(1,1)
	caa(1,2)=coefa(1,2)
	caa(1,3)=coefa(1,3)
	caa(2,1)=coefa(2,1)*coefa(1,1)-coefa(2,1)*coefa(1,1)
	caa(2,2)=coefa(2,1)*coefa(1,2)-coefa(2,2)*coefa(1,1)
	caa(2,3)=coefa(2,1)*coefa(1,3)-coefa(2,3)*coefa(1,1)
	caa(3,1)=coefa(3,1)*coefa(1,1)-coefa(3,1)*coefa(1,1)
	aux32=coefa(3,1)*coefa(1,2)-coefa(3,2)*coefa(1,1)
	aux33=coefa(3,1)*coefa(1,3)-coefa(3,3)*coefa(1,1)
	caa(3,2)=0.
	caa(3,3)=caa(2,2)*aux33-aux32*caa(2,3)
	cbb(1)=coefb(1)
	cbb(2)=coefa(2,1)*coefb(1)-coefa(1,1)*coefb(2)
	aux3=coefa(3,1)*coefb(1)-coefa(1,1)*coefb(3)
	cbb(3)=caa(2,2)*aux3-aux32*cbb(2)
    	inc=(ic-1)*4+2
	inc1=inc+1
	inc2=inc+2
	inc3=inc+3
	arh(ir,1)=ir
	arh(ir,inc2)=r*(cbb(3)/caa(3,3))
	arh(ir,inc1)=(cbb(2)-caa(2,3)*arh(ir,inc2)/r)/caa(2,2)
	arh(ir,inc)=exp((cbb(1)-caa(1,2)*arh(ir,inc1)-
     +             caa(1,3)*arh(ir,inc2)/r)/caa(1,1))

      do it=1,nt
	xc1=arh(ir,inc)*(temp(it)**arh(ir,inc1))*
     -    exp(-arh(ir,inc2)/(r*temp(it)))
	serr=serr+((dlog(xk(it,ir))-dlog(xc1)))**2.
	enddo
	arh(ir,inc3)=sqrt(serr/nt)
	enddo
	enddo
      write(*,*)"     "
      write(*,*)"     "
      write(*,*)"calculo concluido"
      write(*,*)"1. determinado a taxa de reacao via TST "
      write(*,*)"2. aplicacao dos coeficientes de transmissao de: "
      write(*,*)"   2.1 Wigner e   "
      write(*,*)"   2.2 Eckart  "
      write(*,*)"3. determinada a MEP     "
      write(*,*)"4. Ajuste da taxa de reacao conforme Ahrrenius"
      write(*,*)"     "
      write(*,*)"     "
      write(*,*)"Impressao dos resultados finais     "
      write(*,*)"     "
      write(*,*)"     "

      sai(1)="TS"
	sai(2)="A_tst"
	sai(3)="n_tst"
	sai(4)="Ea_tst"
	sai(5)="diff1"
	sai(6)="A_tstW"
	sai(7)="n_tstW"
	sai(8)="Ea_tstW"
	sai(9)="diff2"
	sai(10)="A_tstE"
	sai(11)="n_tstE"
	sai(12)="Ea_tstE"
	sai(13)="diff3"

	arq10="reacao"//larq1//larq2//"arh.dat"

      open(10,file=arq10,status='unknown')
      open(99,file="Properties.dat")
      write(10,'(14(3x, A9))')(sai(ii),ii=1,14)
c	do inc=2,13
	do ir=1,nr
	write(10,'((3x, A6),12(3x, e11.5))')ts(ir),(arh(ir,inc),inc=2,13)
c	write(10,'((3x, A6),12(3x, e11.5))')ts(ir),(arh(ir,inc),ir=1,nr)
	enddo
	close(10)
        write(99,*)"Propridades Reacionais em kcal/mol e cm-1 com ZPVE"
        write(99,*)"       d                 Eo              Freq_Imag
     +            dH                Einv              Ezpe
     +           Tc               RecpTc          Td             RecpTd"
                Einv= 1.0d-3*ea(1,1)- ental(1)
        write(99,'(80(8x, e14.8))')d,1.0d-3*ea(1,1), xifts(1),ental(1),
     +   Einv,1.0d-3*dzpe,Tc,1000.0/Tc,Td,1000.0/Td
        write(99,*)"      "
        write(99,*)"Correcao de Tunelamento"
        write(99,*)"        Temp            1000/Temp          QWigner
     +         QEckart          QBell1958
     +QdTST                 QdTST1935          QBell1958-2Ter"
        do k=1,nt
            write(99,'(60(4x, e14.8))')temp(k),1.0d3/temp(k),xkapa(k)
     +      ,xkappa2(k),xbkapa(k),(dexpxkapa(k,i),i=1,nr),
     +     xb1kapa(k),xb22tkapa(k)
	enddo
        close(99)

	write(*,*)"      "
	write(*,*)"      "
	write(*,*)" Programa terminado     "
	write(*,*)"      "
	write(*,*)"      "

C	pause

      stop
      end
