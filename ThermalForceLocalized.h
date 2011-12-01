/*===- ThermalForceLocalized.h - libSimulation -================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#ifndef THERMALFORCELOCALIZED_H
#define THERMALFORCELOCALIZED_H

#include "Force.h"

class ThermalForceLocalized : public Force {
public:
	ThermalForceLocalized(Cloud * const C, const double thermRed1, const double thermRed2, 
                          const double specifiedRadius);
	~ThermalForceLocalized();

	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

private:
	double heatingRadius, heatVal1, heatVal2; // [m], [N], [N]

    RandCache *evenRandCache, *oddRandCache;
#ifdef DISPATCH_QUEUES
	dispatch_group_t evenRandGroup, oddRandGroup;
	dispatch_queue_t randQueue;
#endif

	void force(const cloud_index currentParticle, const doubleV displacementX, const doubleV displacementY, 
               const RandCache &rc);
};

#endif // THERMALFORCELOCALIZED_H
