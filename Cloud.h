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

class Cloud
{
public:
	Cloud(unsigned int numPar, double sizeOfCloud); //overloaded constructor
	~Cloud();

//public variables:
	unsigned int n;            //number of elements (particles)
	double cloudSize;
	double *k1, *k2, *k3, *k4; //velocityX (Runge-Kutta) tidbits
	double *l1, *l2, *l3, *l4; //positionX (Runge-Kutta) tidbits
	double *m1, *m2, *m3, *m4; //velocityY (Runge-Kutta) tidbits
	double *n1, *n2, *n3, *n4; //positionY (Runge-Kutta) tidbits
	double *o1, *o2, *o3, *o4; //velocityZ (Runge-Kutta) tidbits
	double *p1, *p2, *p3, *p4; //positionZ (Rupge-Kutta) tidbits
	double *x, *y, *z;         //current positions
	double *Vx, *Vy, *Vz;      //current velocities
	double *charge, *mass;
	double *forceX, *forceY, *forceZ;
	__m128d *xCache, *yCache, *zCache, *VxCache, *VyCache, *VzCache;

//public functions:
	//Input: int index, initialPosX, intialPosY, initialPosZ
	//Preconditions: 0 <= index < number of particles
	//Postconditions: x,y,z position of particle #index set to initialPosX,initialPosY,initialPosZ
	void setPosition(const unsigned int index, const double initialPosX, const double initialPosY, const double initialPosZ);

	//Input: int index
	//Preconditions: 0 <= index < number of particles
	//Postconditions: velocity vector of particle #index initialized to zero vector
	void setVelocity(const unsigned int index);

	//Input: int index
	//Preconditions: 0 <= index < number of particles
	//Postconditions: charge of particle #index randomly set, range 5900 to 6100 electrons
	void setCharge(const unsigned int index);

	//Input: int index
	//Preconditions: 0 <= index < number of particles
	//Postconditions: mass of particle #index set according to radius, density
	void setMass(const unsigned int index);

	//Input: fitsfile *file, int *error
	//Preconditions: fitsfile exists, error = 0
	//Postconditions: initial cloud data, including mass & charge, output to file
	void writeCloudSetup(fitsfile * const file, int * const error) const;
	
	//Input: fitsfile *file, int *error, double currentTime
	//Preconditions: fitsfile exists, error = 0, currentTime > 0, writeCloudSetup has previously been called
	//Postconditions: positions and velocities for current time step output to file
	void writeTimeStep(fitsfile * const file, int * const error, double currentTime) const;
   
	//RK4 substep helper functions: 
	const __m128d getx1_pd(const unsigned int i) const;
	const __m128d getx2_pd(const unsigned int i) const;
	const __m128d getx3_pd(const unsigned int i) const;
	const __m128d getx4_pd(const unsigned int i) const;
	
	const __m128d gety1_pd(const unsigned int i) const;
	const __m128d gety2_pd(const unsigned int i) const;
	const __m128d gety3_pd(const unsigned int i) const;
	const __m128d gety4_pd(const unsigned int i) const;

	const __m128d getz1_pd(const unsigned int i) const;
	const __m128d getz2_pd(const unsigned int i) const;
	const __m128d getz3_pd(const unsigned int i) const;
	const __m128d getz4_pd(const unsigned int i) const;

	const __m128d getVx1_pd(const unsigned int i) const;
	const __m128d getVx2_pd(const unsigned int i) const;
	const __m128d getVx3_pd(const unsigned int i) const;
	const __m128d getVx4_pd(const unsigned int i) const;
	
	const __m128d getVy1_pd(const unsigned int i) const;
	const __m128d getVy2_pd(const unsigned int i) const;
	const __m128d getVy3_pd(const unsigned int i) const;
	const __m128d getVy4_pd(const unsigned int i) const;

	const __m128d getVz1_pd(const unsigned int i) const;
	const __m128d getVz2_pd(const unsigned int i) const;
	const __m128d getVz3_pd(const unsigned int i) const;
	const __m128d getVz4_pd(const unsigned int i) const;

//static functions:
	//Input: int numParticles, double cloudSize
	//Preconditions: both inputs positive
	//Postconditions: cloud initialized on spatial grid with side length = 2*cloudSize
	static Cloud * const initializeLine(const unsigned int numParticles, const double cloudSize);   //1D
	static Cloud * const initializeSquare(const unsigned int numParticles, const double cloudSize); //2D
	static Cloud * const initializeCube(const unsigned int numParticles, const double cloudSize);   //3D

	//Input: fitsFile *file, int *error
	//Preconditions: fitsfile exists, error = 0
	//Postconditions: cloud initialized with last time step of fitsfile,
	//	forces and other simulation information extracted as well
	static Cloud * const initializeFromFile(fitsfile * const file, int * const error, double * const currentTime);
};

#endif /* CLOUD_H */
