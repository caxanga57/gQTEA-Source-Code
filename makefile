.SUFFIXES: $(SUFFIXES) .f90 .f
PROG = gqtea.exe
SRCS = SR014_MOD.f90 gqtea.f90 SR001.f90 SR002.f90 SR003.f90 \
SR004.f90 SR005.f90 SR006.f90 SR007.f90 SR008.f90 SR009.f90 SR010.f90 \
SR011.f90 SR012.f90 SR013.f90 SR014.f90 SR015.f90 SR016.f90 \
SR017.f90 SR018.f90 SR019.f90 SR020.f90 date_time.f90 matrixinv.f90 \
SR021.f SR022.f90 SR022A.f90 SR022B.f90 SR024.f90 SR025.f90 SR026.f90 \
SR027.f90 SR028.f90

OBJS = SR014_MOD.o gqtea.o SR001.o SR002.o SR003.o SR004.o SR005.o SR006.o SR007.o \
SR008.o SR009.o SR010.o SR011.o SR012.o SR013.o  SR014.o SR015.o \
SR016.o SR017.o SR018.o SR019.o SR020.o SR021.o SR022.o SR022A.o \
SR022B.o matrixinv.o date_time.o SR024.o SR025.o SR026.o SR027.o SR028.o
 

%.o: %.mod
	
LIBS = 
FC = gfortran

FCFLAGS = -O3

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) -static -o $@ $(OBJS) $(LIBS)
	
clean:
# 	rm -f $(PROG) $(OBJS) *.mod
#	del $(PROG) $(OBJS) *.mod 
#	rm *.o *.mod *.exe
	del *.o *.mod *.exe
.f90.o:
	$(FC) $(FCFLAGS) -c $*.f90

