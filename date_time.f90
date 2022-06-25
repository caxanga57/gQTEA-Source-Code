	SUBROUTINE DATE_TIME(values_in,values_out)
	IMPLICIT NONE
	!values_*(1) = year
	!values_*(2) = month
	!values_*(3) = day
	!values_*(4) = not used here
	!values_*(5) = hours
	!values_*(6) = minute
	!values_*(7) = seconds
	!values_*(8) = miliseconds
	
	INTEGER,DIMENSION(8),INTENT(INOUT) :: values_in,values_out
		IF(values_out(8) < values_in(8)) THEN
			values_out(8) = values_out(8) + 1000
			values_out(7) = values_out(7) - 1
		END IF

		IF(values_out(7) < values_in(7)) THEN
			values_out(7) = values_out(7) + 60
			values_out(6) = values_out(6) - 1
		END IF
		
		IF(values_out(6) < values_in(6)) THEN
			values_out(6) = values_out(6) + 60
			values_out(5) = values_out(5) - 1
		END IF
		
		IF(values_out(5) < values_in(5)) THEN
			values_out(5) = values_out(5) + 24
			values_out(3) = values_out(3) - 1
		END IF
			
		IF(values_out(3) < values_in(3)) THEN
			values_out(3) = values_out(3) + 30
			values_out(2) = values_out(2) - 1
		END IF
		
		IF(values_out(2) < values_in(2)) THEN
			values_out(2) = values_out(2) + 12
			values_out(1) = values_out(1) - 1
		END IF

	END SUBROUTINE DATE_TIME