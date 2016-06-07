/**
* @file  TimeVaryingThermalForce.cpp
* @class TimeVaryingThermalForce TimeVaryingThermalForce.h
*
* @brief Computes a random force to model thermal effects
*		 with a scale factor that changes with a defined function
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#include "TimeVaryingThermalForce.h"

void TimeVaryingThermalForce::force1(const double currentTime) {
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force1(currentTime);
}

void TimeVaryingThermalForce::force2(const double currentTime) {
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force2(currentTime);
}

void TimeVaryingThermalForce::force3(const double currentTime) {
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force3(currentTime);
}

void TimeVaryingThermalForce::force4(const double currentTime) {
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force4(currentTime);
}

/**
* @brief Computes strength of thermal force with function c(t) = m*t + b
		 where the thermal force F has form F = c(t)*L
*
* @param[in] currentTime Current simulation time
**/
inline const double TimeVaryingThermalForce::calculateHeatVal(const double currentTime) const {
	return heatValScale*currentTime + heatValOffset;
}

void TimeVaryingThermalForce::writeForce(fitsfile * const file, int * const error) const {
	ThermalForce::writeForce(file, error);
	
	// add flag indicating that the thermal force is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add TimeVaryingThermalForce bit:
		forceFlags |= ThermalForceFlag;
		forceFlags |= TimeVaryingThermalForceFlag;

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, 
                            const_cast<char *> ("Force configuration."), error);
	}

	if (!*error) {
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("heatingValueScale"), heatValScale, 
                           6, const_cast<char *> ("[N/s] (TimeVaryingThermalForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("heatingValueOffset"), heatValOffset, 
                           6, const_cast<char *> ("[N] (TimeVaryingThermalForce)"), error);
	}
}

void TimeVaryingThermalForce::readForce(fitsfile * const file, int * const error) {
	ThermalForce::readForce(file, error);
	
	if (!*error) {
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("heatingValueScale"), &heatValScale, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("heatingValueOffset"), &heatValOffset, NULL, error);
	}
}
