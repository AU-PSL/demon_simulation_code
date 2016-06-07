/**
* @file  Force.h
* @brief Defines the abstract Force class
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#ifndef FORCE_H
#define FORCE_H

#include "Cloud.h"
#include <vector>

class Force {
public:
	Cloud * const cloud; //!< Cloud object
	
	/**
	* @brief Constructor method
	* @param[in] C Cloud object
	**/
	Force(Cloud * const C) : cloud(C) {} 

	/**
	* @brief Destructor method
	**/
	virtual ~Force() {}

	/**
	* @brief Computes force on particles at first RK substep
	* @param[in] currentTime The current time of the simulation
	**/
	virtual void force1(const double currentTime)=0; // rk substep 1

	/**
	* @brief Computes force on particles at second RK substep
	* @param[in] currentTime The current time of the simulation
	**/
	virtual void force2(const double currentTime)=0; // rk substep 2

	/**
	* @brief Computes force on particles at third RK substep
	* @param[in] currentTime The current time of the simulation
	**/
	virtual void force3(const double currentTime)=0; // rk substep 3

	/**
	* @brief Computes force on particles at fourth RK substep
	* @param[in] currentTime The current time of the simulation
	**/
	virtual void force4(const double currentTime)=0; // rk substep 4

	/**
	* @brief Writes the force information to the fits file
	*
	* @param[in]     file  The name of the fits file
	* @param[in,out] error Error status code
	**/
	virtual void writeForce(fitsfile * const file, int * const error) const=0;	// output force information to file

	/**
	* @brief Reads the force data from a fits file
	*
	* @param[in]     file  The name of the fits file
	* @param[in,out] error Error status code
	**/
	virtual void readForce(fitsfile * const file, int * const error)=0;	// read force information from file
};
	
typedef std::vector<Force *> ForceArray; //!< Vector of Force objects
typedef long force_flags;				 //!< Force flag for  bit-packed FORCES keyword in fits file

//!< Binary assignments for the bit-packed FORCES keyword in Fits file:
enum ForceFlag : force_flags {
	ConfinementForceFlag = 1,          // 000000000000001
	DragForceFlag = 2,                 // 000000000000010
	ShieldedCoulombForceFlag = 4,      // 000000000000100
	GravitationalForceFlag = 8,        // 000000000001000
	ThermalForceFlag = 16,             // 000000000010000
	ThermalForceLocalizedFlag = 32,    // 000000000100000
	DrivingForceFlag = 64,             // 000000001000000
	RotationalForceFlag = 128,         // 000000010000000
	TimeVaryingDragForceFlag = 256,    // 000000100000000
	TimeVaryingThermalForceFlag = 512, // 000001000000000
	MagneticForceFlag = 1024,          // 000010000000000
	ConfinementForceVoidFlag = 2048,   // 000100000000000
    ElectricForceFlag = 4096,          // 001000000000000
	RectConfinementForceFlag = 8192,   // 010000000000000
	VertElectricForceFlag = 16384      // 100000000000000
};

#endif // FORCE_H
