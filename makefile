HDRS = $(wildcard *.h)
SRCS = $(filter-out driver.cpp, $(wildcard *.cpp))
OBJS = $(filter-out driver.o, $(patsubst %.cpp, %.o, $(SRCS)))

# Command name
CMD = DEMON
# Library name
LIBS = libSimulation.a

include makefile.include

# Compling Flags
CXXFLAGS += -funroll-loops -m64 -Wall -mdynamic-no-pic -march=core2 -fomit-frame-pointer -falign-functions -mfpmath=sse -msse4.1 -fno-stack-protector -I $(CFITSIO)

# Linking Flags
LDFLAGS = $(CXXFLAGS) -L$(CFITSIO)

# The below should never need to be changed.
all: $(CMD)

lib: $(LIBS)

$(CMD): $(LIBS) driver.o
	$(CXX) $(LDFLAGS) -o $(@) $(LIBS) driver.o -lcfitsio

$(LIBS): $(OBJS) $(HDRS)
	ar rcs $(@) $(OBJS)

clean:
	-rm $(OBJS) driver.o
	
cleanlib:
	-rm $(LIBS) 

cleanall:
	-rm $(OBJS) driver.o $(CMD) $(LIBS)
