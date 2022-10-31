# makefile - elapid - david gallina - university of kassel - 2022
#...............................................................................

# compiler
FC = gfortran
#FC = ifort
# flags
#FFLAGS = -O2 -C -fbounds-check -fopenmp
FFLAGS = -O3 -C -fopenmp 

# objects in source folder

OBJS = constants.o functions.o random.o main.o efol.o potential.o \
       approx.o lbfgs.o explore.o checkpot.o plds1d.o plds2d.o

LIBS = liblapack.a libblas.a

EXEC = ../bin
 
# Rules and Targets
.SUFFIXES: 
.SUFFIXES: .o .f .f90

.f90.o:
	$(FC) $(FFLAGS) -c  $<
.f.o:
	$(FC) $(FFLAGS) -c  $<
	
$(EXEC)/chara: $(OBJS)
	$(FC) $(FFLAGS) -o $(EXEC)/xc *.o $(LIBS)

clean: 
	rm -f $(OBJS) *.o *.mod
	
# dependencies
$(OBJS):     	    $(LIBS)
main.o:           constants.o functions.o random.o potential.o explore.o \
                  checkpot.o
explore.o:        constants.o functions.o random.o potential.o lbfgs.o efol.o \
                  approx.o
efol.o:           constants.o functions.o random.o potential.o
lbfgs.o:          constants.o functions.o random.o potential.o
potential.o:      constants.o functions.o dpls1d.o dpls2d.o
approx.o:         constants.o functions.o potential.o
dpls1d.o:         constants.o functions.o random.o
dpls2d.o:         constants.o functions.o random.o
random.o:         constants.o
functions.o:      constants.o



