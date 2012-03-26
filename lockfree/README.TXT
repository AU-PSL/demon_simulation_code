Requirements:

Intel CPU with SSE4.2 or higher support. 
Compiler with C++11 support and C++11 standard library support.
libcfitsio.a Avalible from http://heasarc.gsfc.nasa.gov/fitsio/

Compiling:

    Download and build libcfitsio.a from the above link.

1. Using Make:

    In the same directory as the DEMON makefile create a file called 
makefile.include" with the following entries.

CFITSIO = "path to the complied cfitsio"
CXX = "command name of the C++ complier"
CXXFLAGS = "-O3 and any other compiler flags needed."
LIBCXX = "C++ standard library."

For example compling with clang using libc++ where libcfitsio.a has been 
compiled into /foo/cfitsio the makefile.include should appear as

CFITSIO = /foo/cfitsio
CXX = clang++
CXXFLAGS = -O3
LIBCXX = libc++

where -O3 is the flag controlling the optimization level. To compile DEMON type
make in the command line. Parallel builds can be activatled by adding -openmp to
to the CXXFLAGS if using an OpenMP compatable compiler.

2. Using xcode:

    In the xcode project, click on the Simulation project. In the "Build
Settings" Tab, under the "Search Paths" section add the path to the cfitsio 
header files and libcfitsio headers. Under the linking tab check to make sure 
"Other Linking Flags" is set to -lcfitsio.
