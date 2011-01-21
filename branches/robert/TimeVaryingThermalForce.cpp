/*===- TimeVaryingThermalForce.cpp - libSimulation -============================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#include "TimeVaryingThermalForce.h"

TimeVaryingThermalForce1D::TimeVaryingThermalForce1D(Cloud * const myCloud, const double scale, const double offset) 
: ThermalForce1D(myCloud, offset), heatValScale(scale), heatValOffset(offset) {}
TimeVaryingThermalForce2D::TimeVaryingThermalForce2D(Cloud * const myCloud, const double scale, const double offset) 
: ThermalForce2D(myCloud, offset), TimeVaryingThermalForce1D(myCloud, scale, offset) {}
TimeVaryingThermalForce3D::TimeVaryingThermalForce3D(Cloud * const myCloud, const double scale, const double offset) 
: ThermalForce3D(myCloud, offset), TimeVaryingThermalForce2D(myCloud, scale, offset) {}

//1D:
void TimeVaryingThermalForce1D::force1(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce1D::force1(currentTime);
}

void TimeVaryingThermalForce1D::force2(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce1D::force2(currentTime);
}

void TimeVaryingThermalForce1D::force3(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce1D::force3(currentTime);
}

void TimeVaryingThermalForce1D::force4(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce1D::force4(currentTime);
}

//2D:
void TimeVaryingThermalForce2D::force1(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce2D::force1(currentTime);
}

void TimeVaryingThermalForce2D::force2(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce2D::force2(currentTime);
}

void TimeVaryingThermalForce2D::force3(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce2D::force3(currentTime);
}

void TimeVaryingThermalForce2D::force4(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce2D::force4(currentTime);
}

//3D:
void TimeVaryingThermalForce3D::force1(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce3D::force1(currentTime);
}

void TimeVaryingThermalForce3D::force2(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce3D::force2(currentTime);
}

void TimeVaryingThermalForce3D::force3(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce3D::force3(currentTime);
}

void TimeVaryingThermalForce3D::force4(const double currentTime)
{
	heatVal = calculateHeatVal(currentTime);
	ThermalForce3D::force4(currentTime);
}

//calculuateHeatVal:
inline const double TimeVaryingThermalForce1D::calculateHeatVal(const double currentTime) const
{
	return heatValScale*currentTime + heatValOffset;
}

//writeForce:
void TimeVaryingThermalForce1D::writeForce(fitsfile * const file, int * const error) const
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

void TimeVaryingThermalForce1D::readForce(fitsfile * const file, int * const error)
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
