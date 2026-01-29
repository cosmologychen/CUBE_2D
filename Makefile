#OPTIONS+=-CB
OPTIONS+=-DHALOFIND
#OPTIONS+=-DSPEEDTEST

MODFILE:=$(wildcard *.f90)
OBJFILE:=$(addprefix ,$(notdir $(MODFILE:.f90=.o)))

# 添加 utilities 目录的构建目标
all: utilities main.x
	@echo "done"

# 构建 utilities 目录下的程序
utilities:
	@echo "Building utilities..."
	$(MAKE) -C utilities

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

# 添加 clean-all 目标清理整个项目
clean:
	rm -f *.mod *.o *.out *.err *.x *~
	$(MAKE) -C utilities clean

# 专门清理 utilities 的目标
clean-utilities:
	$(MAKE) -C utilities clean

.PHONY: all utilities clean clean-utilities