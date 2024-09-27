# the compiler to use
CC  := g++
CXX := g++-9
XLF := gfortran
#XLF_77 := gfortran -std=legacy -O2 -fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000 -fPIC
TARGET    := run
#H_FILES   := $(shell find ./ -regex [^\#]*\\.h$)
SRC_FILES := $(shell find ./ -regex [^\#]*\\.cpp$)
OBJ_FILES := $(SRC_FILES:.cpp=.o)
SRCFORT_90_FILES := cdbonn_interface.f90
SRCFORT_77_FILES := cdbonn.f
OBJFORT_90_FILES := cdbonn_interface.o
OBJFORT_77_FILES := cdbonn.o

#On mac
#LDFLAGS := -lgfortran -L/usr/local/Cellar/gcc/11.2.0/lib/gcc/11/
# On linux
LDFLAGS := -lgfortran 
CPPFLAGS := -Wall -O2 -ggdb -std=c++11 

FORTFLAGS_90 := -O2 -fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000 -fPIC -std=gnu
FORTFLAGS_77 := -O2 -fdefault-real-8 -fdefault-double-8 -fPIC -static  -std=gnu

$(TARGET): $(OBJ_FILES) $(OBJFORT_90_FILES) $(OBJFORT_77_FILES)
	#$(CC) $^ -o $@ $(LDFLAGS) $(LDLIBS) 

%.o: %.cpp makefile
	$(CC) $(CPPFLAGS) -MMD -MP -c $< -o $@

%.o: %.f90
	$(XLF) $(FORTFLAGS_90) -c $< -o $@
%.o: %.f
	$(XLF) $(FORTFLAGS_77) -c $< -o $@

.PHONY: clean cleanall
# Clean everything but FORTRAN files
clean:
	rm -f $(TARGET) $(OBJ_FILES) $(DEP_FILES) *.o *.exe
# Clean everything
cleanall:
	rm $(OBJFORT_90_FILES) $(OBJFORT_77_FILES)
