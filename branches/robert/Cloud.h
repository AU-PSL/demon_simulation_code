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
#include "VectorCompatibility.h"

typedef unsigned int cloud_index;

class Cloud
{
public:
	Cloud(cloud_index numPar); //overloaded constructor
	~Cloud();

// public variables:
	const cloud_index n;                                   //number of elements (particles)
	double * const k1, * const k2, * const k3, * const k4; //velocityX (Runge-Kutta) tidbits
	double * const l1, * const l2, * const l3, * const l4; //positionsX (Runge-Kutta) tidbits
	double * const m1, * const m2, * const m3, * const m4; //velocityY (Runge-Kutta) tidbits
	double * const n1, * const n2, * const n3, * const n4; //positionsY (Runge-Kutta) tidbits
	double * const o1, * const o2, * const o3, * const o4; //velocityZ (Runge-Kutta) tidbits
	double * const p1, * const p2, * const p3, * const p4; //positionsZ (Runge-Kutta) tidbits
	double * const x, * const y, * const z;                //current positions
	double * const Vx, * const Vy, * const Vz;             //current velocities
	double * const charge, * const mass;
	double * const forceX, * const forceY, * const forceZ;
	__m128d * const xCache, * const yCache, * const zCache;
	__m128d * const VxCache, * const VyCache, * const VzCache;
	
	static const double interParticleSpacing;

//public functions:
	//Input: int index, initialPosX, intialPosY, initialPosZ
	//Preconditions: 0 <= index < number of particles
	//Postconditions: x,y,z position of particle #index set to initialPosX,initialPosY,initialPosZ
	void setPosition(const cloud_index index, const double initialPosX, const double initialPosY, const double initialPosZ) const;

	//Input: int index
	//Preconditions: 0 <= index < number of particles
	//Postconditions: velocity vector of particle #index initialized to zero vector
	void setVelocity(const cloud_index index) const;

	//Input: none
	//Preconditions: none
	//Postconditions: charge of each particle randomly set, range 5900 to 6100 electrons
	void setCharge() const;

	//Input: none
	//Preconditions: none
	//Postconditions: mass of each particle set according to radius, density
	void setMass() const;

	//Input: fitsfile *file, int *error
	//Preconditions: fitsfile exists, error = 0
	//Postconditions: initial cloud data, including mass & charge, output to file
	void writeCloudSetup(fitsfile * const file, int * const error) const;
	
	//Input: fitsfile *file, int *error, double currentTime
	//Preconditions: fitsfile exists, error = 0, currentTime > 0, writeCloudSetup has previously been called
	//Postconditions: positions and velocities for current time step output to file
	void writeTimeStep(fitsfile * const file, int * const error, double currentTime) const;
   
	//RK4 substep helper functions: 
	const __m128d getx1_pd(const cloud_index i) const;
	const __m128d getx2_pd(const cloud_index i) const;
	const __m128d getx3_pd(const cloud_index i) const;
	const __m128d getx4_pd(const cloud_index i) const;

	const __m128d getx1r_pd(const unsigned int i) const;
	const __m128d getx2r_pd(const unsigned int i) const;
	const __m128d getx3r_pd(const unsigned int i) const;
	const __m128d getx4r_pd(const unsigned int i) const;
	
	const __m128d gety1_pd(const cloud_index i) const;
	const __m128d gety2_pd(const cloud_index i) const;
	const __m128d gety3_pd(const cloud_index i) const;
	const __m128d gety4_pd(const cloud_index i) const;

	const __m128d gety1r_pd(const unsigned int i) const;
	const __m128d gety2r_pd(const unsigned int i) const;
	const __m128d gety3r_pd(const unsigned int i) const;
	const __m128d gety4r_pd(const unsigned int i) const;

	const __m128d getz1_pd(const cloud_index i) const;
	const __m128d getz2_pd(const cloud_index i) const;
	const __m128d getz3_pd(const cloud_index i) const;
	const __m128d getz4_pd(const cloud_index i) const;

	const __m128d getz1r_pd(const unsigned int i) const;
	const __m128d getz2r_pd(const unsigned int i) const;
	const __m128d getz3r_pd(const unsigned int i) const;
	const __m128d getz4r_pd(const unsigned int i) const;
    
	const __m128d getVx1_pd(const cloud_index i) const;
	const __m128d getVx2_pd(const cloud_index i) const;
	const __m128d getVx3_pd(const cloud_index i) const;
	const __m128d getVx4_pd(const cloud_index i) const;
    
	const __m128d getVy1_pd(const cloud_index i) const;
	const __m128d getVy2_pd(const cloud_index i) const;
	const __m128d getVy3_pd(const cloud_index i) const;
	const __m128d getVy4_pd(const cloud_index i) const;

	const __m128d getVz1_pd(const cloud_index i) const;
	const __m128d getVz2_pd(const cloud_index i) const;
	const __m128d getVz3_pd(const cloud_index i) const;
	const __m128d getVz4_pd(const cloud_index i) const;

//static functions:
	static Cloud * const initializeLine(const cloud_index numParticles);   //1D
	static Cloud * const initializeSquare(const cloud_index numParticles); //2D
	static Cloud * const initializeCube(const cloud_index numParticles);   //3D

	//Input: fitsFile *file, int *error
	//Preconditions: fitsfile exists, error = 0
	//Postconditions: cloud initialized with last time step of fitsfile,
	//	forces and other simulation information extracted as well
	static Cloud * const initializeFromFile(fitsfile * const file, int * const error, double * const currentTime);
};

#endif /* CLOUD_H */
