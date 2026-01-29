#OPTIONS+=-CB
OPTIONS+=-DHALOFIND
#OPTIONS+=-DSPEEDTEST

MODFILE:=$(wildcard *.f90)
OBJFILE:=$(addprefix ,$(notdir $(MODFILE:.f90=.o)))

all: main.x
	@echo "done"
main.x: $(OBJFILE)
	@echo "Link files:"
	$(FC) $(XFLAG) $(OPTIONS) $(OBJFILE) -o $@ $(FFTFLAG)

$(OBJFILE): variables.o
parameters.o: Makefile basic_functions.f08
variables.o: parameters.o
# fft.o: parameters.o
# particle_mesh.o: variables.o
initialize.o: variables.o Green.o z_checkpoint.txt
main.o: $(OBJFILE)
#$(OBJFILE): variables.o

parameters.o: parameters.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
%.o: %.f90 Makefile
	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)

clean:
	rm -f *.mod *.o *.out *.err *.x *~
