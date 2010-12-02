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
	Cloud(unsigned int numPar, double sizeOfCloud);	//overloaded constructor
	~Cloud();

//public variables:
	unsigned int n;			//number of elements (particles)
	double cloudSize;
	double *k1, *k2, *k3, *k4;	//velocityX (Runge-Kutta) tidbits
	double *l1, *l2, *l3, *l4;	//positionsX (Runge-Kutta) tidbits
	double *m1, *m2, *m3, *m4;	//velocityY (Runge-Kutta) tidbits
	double *n1, *n2, *n3, *n4;	//positionsY (Runge-Kutta) tidbits
	double *x, *y, *Vx, *Vy;	//current positions and velocities=
	double *charge, *mass;
	double *forceX, *forceY;

//public functions:
	//Input: int index
	//Preconditions: 0 <= index < number of particles
	//Postconditions: x,y positions of particle #index randomly set
	void setPosition(const unsigned int index);
	
	//Input: int index, initialPosX, intialPosY
	//Preconditions: 0 <= index < number of particles
	//Postconditions: x,y positions of particle #index set to initialPosX,initialPosY
	void setPosition(const unsigned int index, const double initialPosX, const double initialPosY);

	//Input: int index
	//Preconditions: 0 <= index < number of particles
	//Postconditions: velocity vector of particle #index randomly set
	void setVelocity(const unsigned int index);

	//Input: int index
	//Preconditions: 0 <= index < number of particles
	//Postconditions: charge of particle #index randomly set, range 5900 to 6100 *1.6E-19
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
    
    const __m128d getx1_pd(const unsigned int i) const;
    const __m128d getx2_pd(const unsigned int i) const;
    const __m128d getx3_pd(const unsigned int i) const;
    const __m128d getx4_pd(const unsigned int i) const;
    
    const __m128d gety1_pd(const unsigned int i) const;
    const __m128d gety2_pd(const unsigned int i) const;
    const __m128d gety3_pd(const unsigned int i) const;
    const __m128d gety4_pd(const unsigned int i) const;
    
    const __m128d getVx1_pd(const unsigned int i) const;
    const __m128d getVx2_pd(const unsigned int i) const;
    const __m128d getVx3_pd(const unsigned int i) const;
    const __m128d getVx4_pd(const unsigned int i) const;
    
    const __m128d getVy1_pd(const unsigned int i) const;
    const __m128d getVy2_pd(const unsigned int i) const;
    const __m128d getVy3_pd(const unsigned int i) const;
    const __m128d getVy4_pd(const unsigned int i) const;

//static functions:
	//Input: int numParticles, double cloudSize
	//Preconditions: both inputs positive
	//Postconditions: randomly initialized cloud with radius = cloudSize
	static Cloud * const initializeNew(const unsigned int numParticles, const double cloudSize);

	//Input: int numParticles, double cloudSize
	//Preconditions: both inputs positive
	//Postconditions: cloud initialized on spatial grid with side length = 2*cloudSize
	static Cloud * const initializeGrid(const unsigned int numParticles, const double cloudSize);

	//Input: fitsFile *file, int *error
	//Preconditions: fitsfile exists, error = 0
	//Postconditions: cloud initialized with last time step of fitsfile,
	//	forces and other simulation information extracted as well
	static Cloud * const initializeFromFile(fitsfile * const file, int * const error, double * const currentTime);
};

#endif /* CLOUD_H */
