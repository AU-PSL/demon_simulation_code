Requirements:

Intel CPU with SSE2 or higher support.
libcfitsio.a Avalible from http://heasarc.gsfc.nasa.gov/fitsio/

Compiling:
	
	Download and build libcfitsio.a from the above link. In the makefile change
	
CFITSIO = "PATH/TO/CFITSIO/DIRECTORY"

in the xcode project, right click on libcfitsio.a, longnam.h, and fitsio.h. 
Select "Get Info" and choose the path to the repsective files.
	
	To build type "make" or in xcode click build.
