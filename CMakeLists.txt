cmake_minimum_required (VERSION 2.8)

project (demon_simulation_code CXX)

set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build/bin)
set (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build/lib)

find_library (CFITSIO_LIB NAMES cfitsio PATHS ${CFITSIO_LIB_PATH})
include_directories(${CFITSIO_INC_PATH})

find_library(FFTW_LIBRARY NAMES fftw3 fftw)
set(FFTW_LIBRARIES "${FFTW_LIBRARY}")
if(UNIX AND NOT WIN32)
     find_library(FFTW_libm_LIBRARY NAMES m)
     list(APPEND FFTW_LIBRARIES "${FFTW_libm_LIBRARY}")
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -msse4.2")

list (APPEND demon_sources
	CacheOperator.cpp
	CacheOperator.h
	Cloud.cpp
	Cloud.h
	ConfinementForce.cpp
	ConfinementForce.h
	ConfinementForceVoid.cpp
	ConfinementForceVoid.h
	DragForce.cpp
	DragForce.h
	DrivingForce.cpp
	DrivingForce.h
	Force.h
	GravitationalForce.cpp
	GravitationalForce.h
	Integrator.cpp
	Integrator.h
	MagneticForce.cpp
	MagneticForce.h
	Operator.h
	Parallel.h
	RandomNumbers.cpp
	RandomNumbers.h
	RectConfinementForce.cpp
	RectConfinementForce.h
	RotationalForce.cpp
	RotationalForce.h
	Runge_Kutta2.cpp
	Runge_Kutta2.h
	Runge_Kutta4.cpp
	Runge_Kutta4.h
	ShieldedCoulombForce.cpp
	ShieldedCoulombForce.h
	ThermalForce.cpp
	ThermalForce.h
	ThermalForceLocalized.cpp
	ThermalForceLocalized.h
	TimeVaryingDragForce.cpp
	TimeVaryingDragForce.h
	TimeVaryingThermalForce.cpp
	TimeVaryingThermalForce.h
        ElectricForce.cpp
        ElectricForce.h
	VertElectricForce.cpp
	VertElectricForce.h
)

add_library (simulation STATIC ${demon_sources})
add_executable (DEMON driver.cpp)
add_executable (ANGEL ANGEL.cpp)
add_executable (FFTAnalysis FFTAnalysis.cpp)
add_dependencies (DEMON simulation)
add_dependencies (ANGEL simulation)
add_dependencies (FFTAnalysis simulation)
target_link_libraries (DEMON simulation ${CFITSIO_LIB})
target_link_libraries (ANGEL simulation ${CFITSIO_LIB})
target_link_libraries (FFTAnalysis simulation ${CFITSIO_LIB} ${FFTW_LIBRARIES})
