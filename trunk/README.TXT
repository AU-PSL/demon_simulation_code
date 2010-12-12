Requirements:

Intel CPU with SSE2 or higher support.
libcfitsio.a Avalible from http://heasarc.gsfc.nasa.gov/fitsio/

Compiling:

    Download and build libcfitsio.a from the above link.

1. Using Make:

    In the same directory as the DEMON makefile create a file called 
makefile.include" with the following entries.

CFITSIO = "path to the complied cfitsio"
CXX = "command name of the C++ complier"
CXXFLAGS = "-O3 and any other compiler flags needed."

For example compling with g++ where libcfitsio.a has been compiled into 
/foo/cfitsio the makefile.include should appear as

CFITSIO = /foo/cfitsio
CXX = g++
CXXFLAGS = -O3

where -O3 is the flag controlling the optimization level. To compile DEMON type
make in the command line.

1. Using xcode:

    In the xcode project, right click the Simulation project icon under 
"Groups & files". Select "Get Info" and click the build tab. Edit the entry for 
"User Header Search Paths" with the path to the directories containing 
cfitsio.h. Repeat for both release and debug build configurations. In targets 
click on the DEMON target under targets and select get info. Add libcfitsio.a to
linked libraries by clicking the + and choosing the pre compiled library. Click 
build to compile. Compiled program and libSimulation.a can be found in the 
/xcode/build/<Debug or Release>/ where Debug or Release depends on wheather a 
Debug build or Release has been created.
