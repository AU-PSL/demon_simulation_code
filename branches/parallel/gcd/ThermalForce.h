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
#include "mtrand.h"	// MT header

class ThermalForce : public Force
{	
public:
	ThermalForce(Cloud * const myCloud, const double redFactor);
	~ThermalForce(); // destructor

// public functions:
	// Note: currentTime parameter is necessary (due to parent class) but unused
	virtual void force1(const double currentTime); // rk substep 1
	virtual void force2(const double currentTime); // rk substep 2
	virtual void force3(const double currentTime); // rk substep 3
	virtual void force4(const double currentTime); // rk substep 4

	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);

private:
// private class
	class RandCache {
	public:
		__m128d r;
		double l, h;
		
		RandCache(const __m128d r_ = _mm_set1_pd(0.0), 
		          const double l_ = 0.0, const double h_ = 0.0) 
		: r(r_), l(l_), h(h_) {}
	};
	
// private variables:
	MTRand mt;
	
	RandCache *evenRandCache, *oddRandCache;
	dispatch_group_t evenRandGroup, oddRandGroup;
	dispatch_queue_t randQueue;
		
// private functions:
	void force(const cloud_index currentParticle, const RandCache &rc);

protected:
// protected variables:
	double heatVal;
};

#endif // THERMALFORCE_H
