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
#include "mtrand.h" // MT header
#include "VectorCompatibility.h"

class ThermalForceLocalized : public Force
{	
public:
	ThermalForceLocalized(Cloud * const myCloud, const double thermRed1, const double thermRed2, const double specifiedRadius);	//overloaded constructor
	~ThermalForceLocalized(); // destructor

// public functions:
	// Note: currentTime parameter is necessary (due to parent class) but unused
	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

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
	double heatingRadius;
	double heatVal1;
	double heatVal2;
	
	RandCache *evenRandCache, *oddRandCache;
	dispatch_group_t evenRandGroup, oddRandGroup;
	dispatch_queue_t randQueue;

// private functions:
	void force(const cloud_index currentParticle, const __m128d displacementX, const __m128d displacementY, const RandCache &rc);
};

#endif // THERMALFORCELOCALIZED_H
