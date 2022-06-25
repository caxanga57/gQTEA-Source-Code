	SUBROUTINE matrixinv(A,B,n)
 		! subroutina para calcular a inversa de A usando a eliminacao de Gauss-Jordan elimination
 		! A inversa da matriz A(n,n) e' calculada e amazenada na matriz B(n,n)
 		INTEGER :: i,j,k,l,m,irow
 		INTEGER, INTENT(IN) :: n
 		REAL   :: big,dum
 		REAL, DIMENSION(n,n),INTENT(INOUT) :: A    !A= matriz quatrada 2x2 a ser invertida 
 		REAL, DIMENSION(n,n),INTENT(OUT) :: B      !B= inversa de A

                                            !Este loop cria a matriz identidade
 		DO i = 1,n
    		DO j = 1,n
       			b(i,j) = 0.0
 			END DO
 			b(i,i) = 1.0
 		END DO

 		DO i = 1,n                                  ! Este e' o big loop sobre todas as colunas de A(n,n)
			big = a(i,i)
 			DO j = i,n
 				IF (a(j,i).gt.big) THEN
 					big = a(j,i)
 					irow = j
 				END IF
 			END DO
                                  ! troca a linha i com irow para ambas a() e b() matrizes
 			IF (big.gt.a(i,i)) THEN
 				DO k = 1,n
 					dum = a(i,k)                     ! matriz a()
 					a(i,k) = a(irow,k)
 					a(irow,k) = dum
 					dum = b(i,k) ! matrix b()
 					b(i,k) = b(irow,k)
 					b(irow,k) = dum
 				END DO
 			END IF
                                  ! divide todas entradas na linha i de a(i,j) pelo valor a(i,i);
                                  ! a mesma operacao para a matriz identidade 
 			dum = a(i,i)
 			DO j = 1,n
 				a(i,j) = a(i,j)/dum
 				b(i,j) = b(i,j)/dum
 			END DO
                                  ! zera todas entradas na coluna a(j,i); mesma operacao para indent()
 			DO j = i+1,n
 				dum = a(j,i)
 				DO k = 1,n
 					a(j,k) = a(j,k) - dum*a(i,k)
 					b(j,k) = b(j,k) - dum*b(i,k)
 				END DO
 			END DO
 		END DO
                                  ! subtrai um multiplo apropriado da row j da row j-1
 		DO i = 1,n-1
 			DO j = i+1,n
 				dum = a(i,j)
 				DO l = 1,n
 					a(i,l) = a(i,l)-dum*a(j,l)
 					b(i,l) = b(i,l)-dum*b(j,l)
 				END DO
 			END DO
 		END DO

    END SUBROUTINE matrixinv
