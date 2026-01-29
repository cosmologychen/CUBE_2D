OPTIONS=-DFFTFINE
OPTIONS+=-Dsigma8
OPTIONS+=-DPID
#OPTIONS+=-CB

MODEFILE:=$(wildcard *.f90)
OBJFILE:=$(addprefix ,$(notdir $(MODEFILE:.f90=.o)))

all: ic.x density_power.x dsp.x DTFE_dsp.x fof.x DT_void.x void.x# fof_merge.x #write_gadget_format.x 
	@echo "made:" $^
ic.x: ic.o parameters.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
fof.x: fof.o parameters.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
DT_void.x: DT_void.o parameters.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
void.x: void.o parameters.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
density_power.x: density_power.o parameters.o powerspectrum.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
dsp.x: parameters.o displacement.o powerspectrum.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
DTFE_dsp.x: parameters.o displacement_DTFE.o powerspectrum.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@

parameters.o: Makefile ../basic_functions.f08

parameters.o: ../parameters.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
ic.o: ic.f90 Makefile parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $(FFTFLAG) $< -o $@
fof.o: fof.f90 Makefile parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $(FFTFLAG) $< -o $@
DT_void.o: DT_void.f90 Makefile parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $(FFTFLAG) $< -o $@
void.o: void.f90 Makefile parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $(FFTFLAG) $< -o $@
density_power.o: density_power.f90 Makefile parameters.o powerspectrum.o
	$(FC) $(OFLAG) $(OPTIONS) $(FFTFLAG) $< -o $@
powerspectrum.o: powerspectrum.f90 Makefile parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $(FFTFLAG) $< -o $@
displacement.o: displacement.f90 Makefile parameters.o powerspectrum.o
	$(FC) $(OFLAG) $(OPTIONS) $(FFTFLAG) $< -o $@
displacement_DTFE.o: displacement_DTFE.f90 Makefile parameters.o powerspectrum.o
	$(FC) $(OFLAG) $(OPTIONS) $(FFTFLAG) $< -o $@

*.o: Makefile

clean:
	rm -f *.mod *.o *.out *.err *.x *~
