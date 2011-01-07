/*===- TimeVaryingThermalForce.cpp - libSimulation -============================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#include "TimeVaryingThermalForce.h"

TimeVaryingThermalForce::TimeVaryingThermalForce(Cloud * const myCloud, const double scale, const double offset) 
: ThermalForce(myCloud, offset), heatValScale(scale), heatValOffset(offset) {}

//1D:
void TimeVaryingThermalForce::force1_1D(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force1_1D(currentTime);
}

void TimeVaryingThermalForce::force2_1D(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force2_1D(currentTime);
}

void TimeVaryingThermalForce::force3_1D(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force3_1D(currentTime);
}

void TimeVaryingThermalForce::force4_1D(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force4_1D(currentTime);
}

//2D:
void TimeVaryingThermalForce::force1_2D(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force1_2D(currentTime);
}

void TimeVaryingThermalForce::force2_2D(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force2_2D(currentTime);
}

void TimeVaryingThermalForce::force3_2D(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force3_2D(currentTime);
}

void TimeVaryingThermalForce::force4_2D(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force4_2D(currentTime);
}

//3D:
void TimeVaryingThermalForce::force1_3D(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force1_3D(currentTime);
}

void TimeVaryingThermalForce::force2_3D(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force2_3D(currentTime);
}

void TimeVaryingThermalForce::force3_3D(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force3_3D(currentTime);
}

void TimeVaryingThermalForce::force4_3D(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce::force4_3D(currentTime);
}

inline const double TimeVaryingThermalForce::calculateHeatVal(const double currentTime) const
{
	return heatValScale*currentTime + heatValOffset;
}

void TimeVaryingThermalForce::writeForce(fitsfile * const file, int * const error, const int dimension) const
{
	//move to primary HDU:
	if(!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	//add flag indicating that the thermal force is used:
	if(!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		//add ThermalForce bit:
		forceFlags |= TimeVaryingThermalForceFlag; //compound bitwise OR

		if(*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0;                        //clear above error.

		//add or update keyword:
		if(!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if(!*error)
	{
		//file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("heatingValueScale"), heatValScale, 6, const_cast<char *> ("[N/s] (TimeVaryingThermalForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("heatingValueOffset"), heatValOffset, 6, const_cast<char *> ("[N] (TimeVaryingThermalForce)"), error);
	}
}

void TimeVaryingThermalForce::readForce(fitsfile * const file, int * const error, const int dimension)
{
	//move to primary HDU:
	if(!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if(!*error)
	{
		//file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("heatingValueScale"), &heatValScale, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("heatingValueOffset"), &heatValOffset, NULL, error);
	}
}
