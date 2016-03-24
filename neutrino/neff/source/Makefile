CAMBPATH = $(HOME)/Projects/CMBtools/cosmologist.info/camb/Release

CFLAGS = -O3 -ip -heap-arrays -vec-report0 -openmp -fpp -CB -traceback -W0 -WB
FFLAGS = -CB -traceback -openmp -O3 -ip -W0 -WB -fpp -vec_report0

PRJ  = test_const_omegab_zeq_thetas

OBJ  = subroutines.o inifile.o equations.o bbn.o recfast.o expansion.o $(PRJ).o
FC90 = ifort

default: $(PRJ)

$(PRJ): camb $(OBJ)
	$(FC90) -o $(PRJ) $(OBJ) -I$(CAMBPATH) -L$(CAMBPATH) -lcamb_recfast $(FFLAGS)

%.o: %.f90
	$(FC90) $(FFLAGS) -I$(CAMBPATH) -c $*.f90

camb:
cd $(CAMBPATH) && $(MAKE) --file=Makefile_main libcamb_recfast.a
	
clean:
	rm -f *.mod *.o $(PRJ)
