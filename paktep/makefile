# ======================================================================
# Declarations
# ======================================================================

# The compiler
FCOMP = gfortran

# flags for compiling
FCFLAGS = -g -O2 -fbackslash

# flags for linking
FLFLAGS = -g -O2 -fbackslash

# flags for debugging
FDFLAGS = -g -Wall -Wextra

# libraries needed for linking 
LDFLAGS = -L${MIEKTILS}/paktep/dierckx -ldierckx -L/Users/michaelselby/lapack-3.10.0 -llapack

# directory containing modules
MODDIR = ${MIEKTILS}/paktep 

# List of executables to be built within the package
TEST_NAME = test
TEST_SRC = $(TEST_NAME).f90
TEST = $(TEST_NAME).exe
TEST_OBJ = $(TEST_NAME).o

PRG_NAME = main
PRG_SRC = $(PRG_NAME).f90
PROGRAM = $(PRG_NAME).exe
PRG_OBJ = $(PRG_NAME).o

# List of subroutines to be built within the package

TEST_MODULES = persistence_length.f90 twist_writhe.f90 geom.f90 spline.f90 bend.f90 readers.f90
TEST_SRCS = readers.f90 bend.f90 spline.f90 geom.f90 twist_writhe.f90 persistence_length.f90 $(TEST_SRC)
TEST_OBJECTS = $(TEST_SRCS:.f90=.o)
TEST_SRC_OBJS =$(TEST_MODULES:.f90=.o)

MODULES = processor.f90 persistence_length.f90 twist_writhe.f90 geom.f90 spline.f90 bend.f90 readers.f90
SRCS = readers.f90 bend.f90 spline.f90 geom.f90 twist_writhe.f90 persistence_length.f90 processor.f90 test.f90 $(PRG_SRC)
OBJECTS = $(SRCS:.f90=.o)
SRC_OBJS =$(MODULES:.f90=.o)

# "make" builds all
test: $(TEST)
all: $(TEST) $(PROGRAM) 


# ======================================================================
# General Rules
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends

$(TEST): $(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $(TEST_OBJ) $(TEST_SRC_OBJS) -I$(MODDIR) $(LDFLAGS)

$(PROGRAM): $(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $(PRG_OBJ) $(SRC_OBJS) -I$(MODDIR) $(LDFLAGS)


# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
 

%.o: %.f08
	$(FCOMP) $(FLFLAGS) -c $< -I$(MODDIR)
 
%.o: %.F08
	$(FCOMP) $(FLFLAGS) -c $< -I$(MODDIR)

%.o: %.f90
	$(FCOMP) $(FLFLAGS) -c $< -I$(MODDIR)

%.o: %.F90
	$(FCOMP) $(FLFLAGS) -c $< -I$(MODDIR)


# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
