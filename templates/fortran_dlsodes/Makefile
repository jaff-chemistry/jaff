#executable name
exe = test

# test if ifx is present
wres = $(shell which ifx > /dev/null; echo $$?)
ifeq "$(wres)" "0"
	fc = ifx
	switchOPT = -O3 -g -traceback -xHost -fp-model=precise
	switchDBG = -O0 -check all -warn all -save-temps
	switchDBG += -u -traceback -warn nounused -g
	switchDBG += -init=snan,zero,arrays -ftrapuv -check noarg_temp_created
	switchOMP = -qopenmp
	nowarn = -nowarn
else
	fc = gfortran
	switchOPT = -ffree-line-length-none -O3 -g -fallow-argument-mismatch
	switchDBG = -fbacktrace -g -fallow-argument-mismatch
	switchDBG += -fbounds-check -ffree-line-length-none -O0
	nowarn = -w
endif

# test if icx is present
wres_c = $(shell which icx > /dev/null; echo $$?)
ifeq "$(wres_c)" "0"
	cc = icx
	switchCLIB = -lifcore -limf -lsvml -lintlc -ldl
	switchCDBG = -O0 -Wall -debug full
	cswitchOPT = -O3
else
	cc = gcc
	switchCLIB = -lgfortran -lm
	switchCDBG = -O0
	cswitch =
endif

#default switch
fswitch = $(switchOPT)
cswitch = $(cswitchOPT)

#objects
objs = opkda1.o
objs += opkda2.o
objs += opkdmain.o
objs += user_functions.o
objs += chemistry.o
objs += solver.o

obj_main = main.o

#default target
all:	$(objs) $(obj_main)
	$(fc) $(fswitch) $(objs) $(obj_main) -o $(exe) $(lib)

# Module dependencies
# user_functions has no dependencies; chemistry uses it
user_functions.o: user_functions.f90
	$(fc) $(fswitch) -c $< -o $@ $(lib)

chemistry.o: chemistry.f90 user_functions.o
	$(fc) $(fswitch) -c $< -o $@ $(lib)

# solver uses chemistry; main uses both
solver.o: solver.f90 chemistry.o
	$(fc) $(fswitch) -c $< -o $@ $(lib)

main.o: main.f90 chemistry.o solver.o
	$(fc) $(fswitch) -c $< -o $@ $(lib)

#full debug target
debug: fswitch = $(switchDBG)
debug: all

#clean target
clean:
	rm -f *.o *.mod *__genmod.f90 *~ $(exe)

.PHONY: clean

#rule for f90
%.o:%.f90
	$(fc) $(fswitch) -c $^ -o $@ $(lib)

#rule for f
%.o:%.f
	$(fc) $(fswitch) -c $^ -o $@ $(lib) $(nowarn)

