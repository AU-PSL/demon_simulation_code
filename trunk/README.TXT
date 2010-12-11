Requirements:

Intel CPU with SSE2 or higher support.
libcfitsio.a Avalible from http://heasarc.gsfc.nasa.gov/fitsio/

Compiling:

1. Using Make:

	Download and build libcfitsio.a from the above link. In the makefile change
	
CFITSIO = "PATH/TO/CFITSIO/DIRECTORY"

type "make".


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
