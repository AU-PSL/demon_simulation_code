/*===- Cloud.h - libSimulation -================================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#ifndef CLOUD_H
#define CLOUD_H

#include "fitsio.h"
#include "Parallel.h"
#include "RandomNumbers.h"

class Cloud {	
public:
	Cloud(const cloud_index numPar);
	~Cloud();

	const cloud_index n; // number of elements (particles)
	double * const x, * const y, * const Vx, * const Vy; // current positions and velocities
	double * const charge, * const mass;
	double * const k1, * const k2, * const k3, * const k4; // velocityX (Runge-Kutta) tidbits
	double * const l1, * const l2, * const l3, * const l4; // positionsX (Runge-Kutta) tidbits
	double * const m1, * const m2, * const m3, * const m4; // velocityY (Runge-Kutta) tidbits
	double * const n1, * const n2, * const n3, * const n4; // positionsY (Runge-Kutta) tidbits
	double * const forceX, * const forceY;
	__m128d * const xCache, * const yCache, * const VxCache, * const VyCache;
	
	RandomNumbers rands;
	
	static const double interParticleSpacing;
	static const double electronCharge;
	static const double epsilon0;
	static const double particleRadius;
    static const double dustParticleMassDensity;

	void setPosition(const cloud_index index, const double initialPosX, const double initialPosY) const;
	void setVelocity(const cloud_index index) const;
	void setCharge(const double qMean, const double qSigma);	
	void setMass(const double rMean, const double rSigma);

	void writeCloudSetup(fitsfile * const file, int &error) const;
	void writeTimeStep(fitsfile * const file, int &error, double currentTime) const;
    
	const __m128d getx1_pd(const cloud_index i) const;
	const __m128d getx2_pd(const cloud_index i) const;
	const __m128d getx3_pd(const cloud_index i) const;
	const __m128d getx4_pd(const cloud_index i) const;

	const __m128d getx1r_pd(const cloud_index i) const;
	const __m128d getx2r_pd(const cloud_index i) const;
	const __m128d getx3r_pd(const cloud_index i) const;
	const __m128d getx4r_pd(const cloud_index i) const;
	
	const __m128d gety1_pd(const cloud_index i) const;
	const __m128d gety2_pd(const cloud_index i) const;
	const __m128d gety3_pd(const cloud_index i) const;
	const __m128d gety4_pd(const cloud_index i) const;
    
	const __m128d gety1r_pd(const cloud_index i) const;
	const __m128d gety2r_pd(const cloud_index i) const;
	const __m128d gety3r_pd(const cloud_index i) const;
	const __m128d gety4r_pd(const cloud_index i) const;
    
	const __m128d getVx1_pd(const cloud_index i) const;
	const __m128d getVx2_pd(const cloud_index i) const;
	const __m128d getVx3_pd(const cloud_index i) const;
	const __m128d getVx4_pd(const cloud_index i) const;
	
	const __m128d getVy1_pd(const cloud_index i) const;
	const __m128d getVy2_pd(const cloud_index i) const;
	const __m128d getVy3_pd(const cloud_index i) const;
	const __m128d getVy4_pd(const cloud_index i) const;
	
	static Cloud * const initializeGrid(const cloud_index numParticles,
										const double rMean, const double rSigma,
                                        const double qMean, const double qSigma);
	static Cloud * const initializeFromFile(fitsfile * const file, int &error, 
                                            double * const currentTime);
};

#endif // CLOUD_H
