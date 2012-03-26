/*===- ThermalForce.h - libSimulation -=========================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#ifndef THERMALFORCE_H
#define THERMALFORCE_H

#include "Force.h"

class ThermalForce : public Force {
public:
	ThermalForce(Cloud * const C, const double redFactor);
	~ThermalForce();

	virtual void force1(const double currentTime); // rk substep 1
	virtual void force2(const double currentTime); // rk substep 2
	virtual void force3(const double currentTime); // rk substep 3
	virtual void force4(const double currentTime); // rk substep 4

	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);

private:
    RandCache *evenRandCache, *oddRandCache;
#ifdef DISPATCH_QUEUES
    dispatch_group_t evenRandGroup, oddRandGroup;
    dispatch_queue_t randQueue;
#endif

	void force(const cloud_index currentParticle, const RandCache &RC);
    static const doubleV randomCos(const RandCache &RC);
    static const doubleV randomSin(const RandCache &RC);
    
protected:
	double heatVal; // [N]
};

#endif // THERMALFORCE_H
