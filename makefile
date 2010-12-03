# For each file add the object file that it needs
# to be made.
OBJS = Runge_Kutta.o Cloud.o ShieldedCoulombForce.o DragForce.o ConfinementForce.o RectConfinementForce.o ThermalForce.o DrivingForce.o ThermalForceLocalized.o RotationalForce.o mtrand.o TimeVaryingDragForce.o TimeVaryingThermalForce.o
SRCS = $(OBJS,.o=.cpp)
HDRS = $(OBJS,.o=.h)

# Command name
CMD = DEMON
# Library name
LIBS = libSimulation.a
CFITSIO = /Users/robert/Documents/Jefferson/cfitsio

# Compling Flags
# Switch comments if compiler doesn't
# support OpenMP.
#CXXFLAGS = -O3 -funroll-loops -m64 -Wall -mdynamic-no-pic -march=core2 -fomit-frame-pointer -falign-functions -mfpmath=sse -msse4.1 -fno-stack-protector -I $(CFITSIO)
#CXXFLAGS = -O3 -fopenmp -funroll-loops -m64 -Wall -mdynamic-no-pic -march=core2 -fomit-frame-pointer -falign-functions -mfpmath=sse -msse4.1 -fno-stack-protector -I $(CFITSIO)
CXXFLAGS = -O3 -m64 -Wall -mdynamic-no-pic -march=core2 -fomit-frame-pointer -msse4.1 -fno-stack-protector -I $(CFITSIO)
#CXXFLAGS = -O4 -fopenmp -funroll-loops -m64 -Wall -mdynamic-no-pic -march=core2 -fomit-frame-pointer -falign-functions -mfpmath=sse -msse4.1 -fno-stack-protector -I $(CFITSIO)

# Linking Flags
LDFLAGS = $(CXXFLAGS) -L$(CFITSIO)

# Complier command
#CXX = g++-4.2
#CXX = clang++
CXX = /Developer/usr/bin/llvm-g++-4.2

# The below should never need to be changed.
all: $(CMD) $(LIBS) $(HDRS)

lib: $(LIBS) $(HDRS)

$(CMD): $(LIBS) driver.o
	$(CXX) $(LDFLAGS) -o $(@) $(LIBS) driver.o -lcfitsio

$(LIBS): $(OBJS)
	ar rcs $(@) $(OBJS)

$driver.o: driver.cpp
	$(CXX) -c driver.cpp driver.o

clean:
	-rm $(OBJS) driver.o
	
cleanlib:
	-rm $(LIBS) 

cleanall:
	-rm $(OBJS) driver.o $(CMD) $(LIBS)
