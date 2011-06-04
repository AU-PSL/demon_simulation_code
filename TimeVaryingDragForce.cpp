/*===- TimeVaryingDragForce.cpp - libSimulation -===============================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#include "TimeVaryingDragForce.h"

TimeVaryingDragForce1D::TimeVaryingDragForce1D(Cloud * const myCloud, const double scale, const double offset) 
: DragForce3D(myCloud, -offset), scaleConst(scale), offsetConst(offset) {}
TimeVaryingDragForce2D::TimeVaryingDragForce2D(Cloud * const myCloud, const double scale, const double offset) 
: TimeVaryingDragForce1D(myCloud, scale, offset) {}
TimeVaryingDragForce3D::TimeVaryingDragForce3D(Cloud * const myCloud, const double scale, const double offset) 
: TimeVaryingDragForce2D(myCloud, scale, offset) {}

//1D:
void TimeVaryingDragForce1D::force1(const double currentTime)
{
	dragConst = calculateGamma(currentTime);
	DragForce1D::force1(currentTime);
}

void TimeVaryingDragForce1D::force2(const double currentTime)
{
	dragConst = calculateGamma(currentTime);
	DragForce1D::force2(currentTime);
}

void TimeVaryingDragForce1D::force3(const double currentTime)
{
	dragConst = calculateGamma(currentTime);
	DragForce1D::force3(currentTime);
}

void TimeVaryingDragForce1D::force4(const double currentTime)
{
	dragConst = calculateGamma(currentTime);
	DragForce1D::force4(currentTime);
}

//2D:
void TimeVaryingDragForce2D::force1(const double currentTime)
{
	dragConst = calculateGamma(currentTime);
	DragForce2D::force1(currentTime);
}

void TimeVaryingDragForce2D::force2(const double currentTime)
{
	dragConst = calculateGamma(currentTime);
	DragForce2D::force2(currentTime);
}

void TimeVaryingDragForce2D::force3(const double currentTime)
{
	dragConst = calculateGamma(currentTime);
	DragForce2D::force3(currentTime);
}

void TimeVaryingDragForce2D::force4(const double currentTime)
{
	dragConst = calculateGamma(currentTime);
	DragForce2D::force4(currentTime);
}

//3D:
void TimeVaryingDragForce3D::force1(const double currentTime)
{
	dragConst = calculateGamma(currentTime);
	DragForce3D::force1(currentTime);
}

void TimeVaryingDragForce3D::force2(const double currentTime)
{
	dragConst = calculateGamma(currentTime);
	DragForce3D::force2(currentTime);
}

void TimeVaryingDragForce3D::force3(const double currentTime)
{
	dragConst = calculateGamma(currentTime);
	DragForce3D::force3(currentTime);
}

void TimeVaryingDragForce3D::force4(const double currentTime)
{
	dragConst = calculateGamma(currentTime);
	DragForce3D::force4(currentTime);
}

//calculateGamma:
inline const double TimeVaryingDragForce1D::calculateGamma(const double currentTime) const
{
	return -(scaleConst*currentTime + offsetConst);
}

//writeForce:
void TimeVaryingDragForce1D::writeForce(fitsfile * const file, int * const error) const
{
	//move to primary HDU:
	if (!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	//add flag indicating that the time varying drag force is used:
	if (!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		//add TimeVaryingDragForce bit:
		forceFlags |= TimeVaryingDragForceFlag;

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; //clear above error.

		//add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configureation."), error);
	}

	if (!*error)
	{
		//file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("TVDragScaleConst"), scaleConst, 6, const_cast<char *> ("[s^-2] (TimeVaryingDragForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("TVDragOffsetConst"), offsetConst, 6, const_cast<char *> ("[s^-1] (TimeVaryingDragForce)"), error);
	}
}

//readForce:
void TimeVaryingDragForce1D::readForce(fitsfile * const file, int * const error)
{
	//move to primary HDU:
	if (!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if (!*error)
	{
		//file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("TVDragScaleConst"), &scaleConst, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("TVDragOffsetConst"), &offsetConst, NULL, error);
	}
}
