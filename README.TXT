Requirements:

Intel CPU with SSE4.2 or higher support. 
Compiler with C++11 support and C++11 standard library support.
libcfitsio.a Avalible from http://heasarc.gsfc.nasa.gov/fitsio/

Compiling:

    Download and build libcfitsio.a from the above link.

1. Using CMake:

    DEMON requires that the cfitsio library and include paths be set with the
following command.

cmake -DCFITSIO_INC_PATH=/path/to/include -DCFITSIO_LIB_PATH=/path/to/lib

extra build flags maybe included by setting

cmake -DCMAKE_CXX_FLAGS="flags"

A optimized build maybe specifed by

cmake -DCMAKE_BUILD_TYPE=Release

Cmake will generate a make file. The built executable can be found in 

build/bin

2. Using xcode:

    In the xcode project, click on the Simulation project. In the "Build
Settings" Tab, under the "Search Paths" section add the path to the cfitsio 
header files and libcfitsio headers. Under the linking tab check to make sure 
"Other Linking Flags" is set to -lcfitsio.
