cmake_minimum_required (VERSION 2.8)

project (demon_simulation_code CXX)

set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build/bin)
set (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build/lib)

find_library (CFITSIO_LIB NAMES cfitsio PATHS ${CFITSIO_LIB_PATH})
include_directories(${CFITSIO_INC_PATH})

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -msse4.2 -pedantic")

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
)

add_library (simulation STATIC ${demon_sources})
add_executable (DEMON driver.cpp)
add_dependencies (DEMON simulation)
target_link_libraries (DEMON simulation ${CFITSIO_LIB})
