Requirements:

Intel CPU with SSE4.2 or higher support. 
Compiler with C++11 support and C++11 standard library support.
libcfitsio.a Avalible from http://heasarc.gsfc.nasa.gov/fitsio/

Before installing DEMON, be sure you have the cfitsio library installed from
here: http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html

Also make sure you have CMake installed from here:
http://www.cmake.org/

A step by step installation guide for mac and linux devices will be in section 4

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


3. Multicore Support on Linux:

For parallel processing support on linux, use the -fopenmp flag,

cmake -DCMAKE_CXX_FLAGS="-fopenmp"



4. Step-by-step installation on linux and mac devices

     Install cfitsio, if you use your distribution's software repository,
be sure to keep track of where the library and include paths are located.
Otherwise do the following:

mkdir builds

cd builds

wget ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3360.tar.gz

tar xvzf cfitsio3360.tar.gz

cd cfitsio

./configure

make

cd ..


     Now, make sure that cmake is installed correctly. More than likely,
you can do this with your package manager, but if you insist on doing this
through the command line, do the following:

wget http://www.cmake.org/files/v2.8/cmake-2.8.12.2.tar.gz

tar xvzf cmake-2.8.12.2.tar.gz

cd cmake-2.8.12.2

./bootstrap

make

make install


     Now we can start to work with DEMON. Do the following:

svn checkout http://demonsimulationcode.googlecode.com/svn/trunk/ demonsimulationcode-read-only

cmake -DCFITSIO_INC_PATH=../cfitsio -DCFITSIO_LIB_PATH=../lib

make


     Now DEMON should be on your system, with the executable in the following
directory: ~builds/demonsimulationcode-read-only/build/bin/ . Now we just need
to make this directory callable in any directory on your system. do the following:

pwd

     Now place a path to that directory in your bashrc file. On Linux:

echo export PATH=$PATH:/home/USER/builds/demonsimulationcode-read-only//build/bin >> ~/.bashrc 

     On Mac:

echo export PATH=$PATH:/home/USER/builds/demonsimulationcode-read-only//build/bin >> ~/.PROFILE

     And there you have it, a full copy of DEMON on your machine!

5. Notes on github

    DEMON is now available on github and can be mondified with both svn and git.

6. Looking at Data

    DEMON data will be stored in a fits file. Please use additional software, such as astropy.io to use this data.
