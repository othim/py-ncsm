# the compiler to use
CC  := gcc
CXX := g++-9
XLF := gfortran

SRCFORT_90_FILES := tbi-constants.f90 my_braket.f90
OBJFORT_90_FILES := tbi-constants.o my_braket.o

#On mac
LDFLAGS := -lgfortran -L/usr/local/Cellar/gcc/12.2.0/lib/gcc/12/ 
# On linux
#LDFLAGS := -lgfortran 

FORTFLAGS_90 := -O2 -Wall -g -fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000 -fPIC -std=gnu


so: $(OBJFORT_90_FILES) 
	$(XLF) -shared -o lib.so $(OBJFORT_90_FILES) $(LDFLAGS)

%.o: %.f90
	$(XLF) $(FORTFLAGS_90) -c $< -o $@

.PHONY: clean cleanall
# Clean everything but FORTRAN files
clean:
	rm *.o *.mod
