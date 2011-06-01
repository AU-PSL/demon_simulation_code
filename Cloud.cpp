/*===- Cloud.cpp - libSimulation -==============================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#include "Cloud.h"
#include <cmath>
#include <sstream>

using namespace std;

//typical spacing of physical clouds:
const double Cloud::interParticleSpacing = 0.0003;

Cloud::Cloud(const cloud_index numPar) : n(numPar),
k1(new double[n]), k2(new double[n]), k3(new double[n]), k4(new double[n]),
l1(new double[n]), l2(new double[n]), l3(new double[n]), l4(new double[n]),
m1(new double[n]), m2(new double[n]), m3(new double[n]), m4(new double[n]),
n1(new double[n]), n2(new double[n]), n3(new double[n]), n4(new double[n]),
o1(new double[n]), o2(new double[n]), o3(new double[n]), o4(new double[n]),
p1(new double[n]), p2(new double[n]), p3(new double[n]), p4(new double[n]),
x(new double[n]), y(new double[n]), z(new double[n]),
Vx(new double[n]), Vy(new double[n]), Vz(new double[n]),
charge(new double[n]), mass(new double[n]),
forceX(new double[n]), forceY(new double[n]), forceZ(new double[n]),
xCache(new __m128d[n/2]), yCache(new __m128d[n/2]), zCache(new __m128d[n/2]),
VxCache(new __m128d[n/2]), VyCache(new __m128d[n/2]), VzCache(new __m128d[n/2]) {}

Cloud::~Cloud() 
{
	delete[] k1; delete[] k2; delete[] k3; delete[] k4;
	delete[] l1; delete[] l2; delete[] l3; delete[] l4;
	delete[] m1; delete[] m2; delete[] m3; delete[] m4;
	delete[] n1; delete[] n2; delete[] n3; delete[] n4;
	delete[] o1; delete[] o2; delete[] o3; delete[] o4;
	delete[] p1; delete[] p2; delete[] p3; delete[] p4;
	delete[] x; delete[] y; delete[] z;
	delete[] Vx; delete[] Vy; delete[] Vz;
	delete[] charge; delete[] mass; 
	delete[] forceX; delete[] forceY; delete[] forceZ;
	delete[] xCache; delete[] yCache; delete[] zCache;
	delete[] VxCache; delete[] VxCache; delete[] VzCache;
}

inline void Cloud::setPosition(const cloud_index index, const double xVal, const double yVal, const double zVal) const
{
	x[index] = xVal;
	y[index] = yVal;
	z[index] = zVal;
}

inline void Cloud::setVelocity(const cloud_index index) const
{
	Vx[index] = 0.0;
	Vy[index] = 0.0;
	Vz[index] = 0.0;
}

inline void Cloud::setCharge() const
{
	srand((int)time(NULL));
	for(cloud_index i = 0; i < n; i++)
		charge[i] = (rand()%201 + 5900)*1.6E-19;
}

inline void Cloud::setMass() const
{
	const double radius = 1.45E-6;
	const double particleDensity = 2200.0;
	const double particleMass = (4.0/3.0)*M_PI*radius*radius*radius*particleDensity;
	for(cloud_index i = 0; i < n; i++)
		mass[i] = particleMass;
}

Cloud * const Cloud::initializeLine(const cloud_index numParticles)
{
	Cloud * const cloud = new Cloud(numParticles);

	//Ensure that a given row is centered for both odd and even numbers of particles:
	const double cloudHalfSize = (double)numParticles/2.0*interParticleSpacing
		- ((numParticles%2) ? interParticleSpacing/2.0 : 0.0);

	cloud->setCharge();
	cloud->setMass();

	//initialize as a line centered on the x-axis:
	for(cloud_index i = 0; i < numParticles; i++)
	{
		cloud->setPosition(i,
			cloudHalfSize - (double)(i)*interParticleSpacing, 0.0, 0.0);
		cloud->setVelocity(i);
	}

	return cloud;
}

Cloud * const Cloud::initializeSquare(const cloud_index numParticles)
{
	Cloud * const cloud = new Cloud(numParticles);

	const cloud_index sqrtNumPar = (cloud_index)floor(sqrt(numParticles));

	//Ensure that a given row is centered for both odd and even numbers of particles:
	const double cloudHalfSize = (double)sqrtNumPar/2.0*interParticleSpacing
		- ((sqrtNumPar%2) ? interParticleSpacing/2.0 : 0.0);

	cloud->setCharge();
	cloud->setMass();

	//initialize as a square in the xy-plane with sqrtNumPar particles on each side:
	for(cloud_index i = 0; i < numParticles; i++)
	{
		cloud->setPosition(i,
			cloudHalfSize - (double)(i%sqrtNumPar)*interParticleSpacing,
			cloudHalfSize - (double)(i/sqrtNumPar)*interParticleSpacing, 0.0);
		cloud->setVelocity(i);
	}

	return cloud;
}

Cloud * const Cloud::initializeCube(const cloud_index numParticles)
{
	Cloud * const cloud = new Cloud(numParticles);

	const cloud_index cbrtNumPar = (cloud_index)floor(pow(numParticles,1.0/3.0));

	//Ensure that a given row is centered for both odd and even numbers of particles:
	const double cloudHalfSize = (double)cbrtNumPar/2.0*interParticleSpacing
		- ((cbrtNumPar%2) ? interParticleSpacing/2.0 : 0.0);

	cloud->setCharge();
	cloud->setMass();

	//initialize as a cube with cbrtNumPar particles on each side:
	for(cloud_index i = 0; i < numParticles; i++)
	{
		cloud->setPosition(i,
			cloudHalfSize - (double)(i%cbrtNumPar)*interParticleSpacing,
			cloudHalfSize - (double)((i/cbrtNumPar)%cbrtNumPar)*interParticleSpacing,
			cloudHalfSize - (double)(i/(cbrtNumPar*cbrtNumPar))*interParticleSpacing);
		cloud->setVelocity(i);
	}

	return cloud;
}

Cloud * const Cloud::initializeFromFile(fitsfile * const file, int * const error, double * const currentTime)
{
	int *anyNull = NULL;
	long numParticles = 0;
	long numTimeSteps = 0;

	//move to CLOUD HDU:
	if(!*error)
		fits_movnam_hdu(file, BINARY_TBL, const_cast<char *> ("CLOUD"), 0, error);

	//get number of particles:
	if(!*error)
		fits_get_num_rows(file, &numParticles, error);

	//create cloud:
	Cloud* cloud = new Cloud((cloud_index)numParticles); 

	//read mass and charge information:
	if(!*error)
	{
		//file, column #, starting row, first element, num elements, mass array, pointless pointer, error
		fits_read_col_dbl(file, 1, 1, 1, numParticles, 0.0, cloud->charge, anyNull, error);
		fits_read_col_dbl(file, 2, 1, 1, numParticles, 0.0, cloud->mass, anyNull, error);
	}

	//move to TIME_STEP HDU:
	if(!*error)
		fits_movnam_hdu(file, BINARY_TBL, const_cast<char *> ("TIME_STEP"), 0, error);

	//get number of time steps:
	if(!*error)
		fits_get_num_rows(file, &numTimeSteps, error);

	if(!*error)
	{
		if(currentTime)
			fits_read_col_dbl(file, 1, numTimeSteps, 1, 1, 0.0, currentTime, anyNull, error);

		fits_read_col_dbl(file, 2, numTimeSteps, 1, numParticles, 0.0, cloud->x, anyNull, error);
		fits_read_col_dbl(file, 3, numTimeSteps, 1, numParticles, 0.0, cloud->y, anyNull, error);
		fits_read_col_dbl(file, 4, numTimeSteps, 1, numParticles, 0.0, cloud->z, anyNull, error);
		fits_read_col_dbl(file, 5, numTimeSteps, 1, numParticles, 0.0, cloud->Vx, anyNull, error);
		fits_read_col_dbl(file, 6, numTimeSteps, 1, numParticles, 0.0, cloud->Vy, anyNull, error);
		fits_read_col_dbl(file, 7, numTimeSteps, 1, numParticles, 0.0, cloud->Vz, anyNull, error);
	}

	return cloud;
}
//writeCloudSetup:
void Cloud::writeCloudSetup(fitsfile * const file, int * const error) const
{
	//format number of elements of type double as string, e.g. 1024D
	stringstream numStream;
	numStream << n << "D";
	const string numString = numStream.str();

	char *ttypeCloud[] = {const_cast<char *> ("CHARGE"), const_cast<char *> ("MASS")};
	char *tformCloud[] = {const_cast<char *> ("D"), const_cast<char *> ("D")};
	char *tunitCloud[] = {const_cast<char *> ("C"), const_cast<char *> ("kg")};

	//write mass and charge:
	if (!*error)
		//file, storage type, num rows, num columns, ...
		fits_create_tbl(file, BINARY_TBL, n, 2, ttypeCloud, tformCloud, tunitCloud, "CLOUD", error);
	if(!*error)
	{
		//file, column #, starting row, first element, num elements, mass array, error
		fits_write_col_dbl(file, 1, 1, 1, n, charge, error);
		fits_write_col_dbl(file, 2, 1, 1, n, mass, error);
	}

	//write position and velocity:
	char *ttypeRun[] = {const_cast<char *> ("TIME"),
		const_cast<char *> ("X_POSITION"), const_cast<char *> ("Y_POSITION"), const_cast<char *> ("Z_POSITION"), 
		const_cast<char *> ("X_VELOCITY"), const_cast<char *> ("Y_VELOCITY"), const_cast<char *> ("Z_VELOCITY")};
	char *tformRun[] = {const_cast<char *> ("D"), 
		const_cast<char *> (numString.c_str()), const_cast<char *> (numString.c_str()), const_cast<char *> (numString.c_str()), 
		const_cast<char *> (numString.c_str()), const_cast<char *> (numString.c_str()), const_cast<char *> (numString.c_str())};
	char *tunitRun[] = {const_cast<char *> ("s"),
		const_cast<char *> ("m"), const_cast<char *> ("m"), const_cast<char *> ("m"), 
		const_cast<char *> ("m/s"), const_cast<char *> ("m/s"), const_cast<char *> ("m/s")};

	if(!*error)
		fits_create_tbl(file, BINARY_TBL, 0, 7, ttypeRun, tformRun, tunitRun, "TIME_STEP", error);
		//n.b. num rows automatically incremented. Increment from 0 as opposed to preallocating 
		//to ensure proper output in the event of program interruption.

	if(!*error)
	{
		double time = 0.0;
		fits_write_col_dbl(file, 1, 1, 1, 1, &time, error);
		fits_write_col_dbl(file, 2, 1, 1, n, x, error);
		fits_write_col_dbl(file, 3, 1, 1, n, y, error);
		fits_write_col_dbl(file, 4, 1, 1, n, z, error);
		fits_write_col_dbl(file, 5, 1, 1, n, Vx, error);
		fits_write_col_dbl(file, 6, 1, 1, n, Vy, error);
		fits_write_col_dbl(file, 7, 1, 1, n, Vz, error);
	}

	//write buffer, close file, reopen at same point:
	fits_flush_file(file, error);
}

void Cloud::writeTimeStep(fitsfile * const file, int * const error, double currentTime) const
{
	if (!*error)
	{
		long numRows = 0;
		fits_get_num_rows(file, &numRows, error);
		fits_write_col_dbl(file, 1, ++numRows, 1, 1, &currentTime, error);
		fits_write_col_dbl(file, 2, numRows, 1, n, x, error);
		fits_write_col_dbl(file, 3, numRows, 1, n, y, error);
		fits_write_col_dbl(file, 4, numRows, 1, n, z, error);
		fits_write_col_dbl(file, 5, numRows, 1, n, Vx, error);
		fits_write_col_dbl(file, 6, numRows, 1, n, Vy, error);
		fits_write_col_dbl(file, 7, numRows, 1, n, Vz, error);
	}

	//write buffer, close file, reopen at same point:
	fits_flush_file(file, error);
}

// 4th order Runge-Kutta substep helper methods:

// x-position helper functions -------------------------------------------------
const __m128d Cloud::getx1_pd(const cloud_index i) const
{
	return _mm_load_pd(x + i); //x
}

const __m128d Cloud::getx2_pd(const cloud_index i) const
{
    return xCache[i/2]; //x + l1/2
}

const __m128d Cloud::getx3_pd(const cloud_index i) const
{
    return xCache[i/2]; //x + l2/2
}

const __m128d Cloud::getx4_pd(const cloud_index i) const
{
    return xCache[i/2]; //x + l3
}

const __m128d Cloud::getx1r_pd(const cloud_index i) const
{
	return _mm_loadr_pd(x + i);
}

const __m128d Cloud::getx2r_pd(const cloud_index i) const
{
	const cloud_index j = i/2;
	return _mm_shuffle_pd(xCache[j], xCache[j], _MM_SHUFFLE2(0, 1));
}

const __m128d Cloud::getx3r_pd(const cloud_index i) const
{
	const cloud_index j = i/2;
	return _mm_shuffle_pd(xCache[j], xCache[j], _MM_SHUFFLE2(0, 1));
}

const __m128d Cloud::getx4r_pd(const cloud_index i) const
{
	const cloud_index j = i/2;
	return _mm_shuffle_pd(xCache[j], xCache[j], _MM_SHUFFLE2(0, 1));
}

// y-position helper functions -------------------------------------------------
const __m128d Cloud::gety1_pd(const cloud_index i) const
{
    return _mm_load_pd(y + i); // y
}

const __m128d Cloud::gety2_pd(const cloud_index i) const
{
    return yCache[i/2]; // y + n1/2
}

const __m128d Cloud::gety3_pd(const cloud_index i) const
{
    return yCache[i/2]; //y + n2/2
}

const __m128d Cloud::gety4_pd(const cloud_index i) const
{
    return yCache[i/2]; // y + n3
}

const __m128d Cloud::gety1r_pd(const cloud_index i) const
{
    return _mm_loadr_pd(y + i);
}

const __m128d Cloud::gety2r_pd(const cloud_index i) const
{
    const cloud_index j = i/2;
    return _mm_shuffle_pd(yCache[j], yCache[j], _MM_SHUFFLE2(0, 1));
}

const __m128d Cloud::gety3r_pd(const cloud_index i) const
{
    const cloud_index j = i/2;
    return _mm_shuffle_pd(yCache[j], yCache[j], _MM_SHUFFLE2(0, 1));
}

const __m128d Cloud::gety4r_pd(const cloud_index i) const
{
    const cloud_index j = i/2;
    return _mm_shuffle_pd(yCache[j], yCache[j], _MM_SHUFFLE2(0, 1));
}

// z-position helper functions -------------------------------------------------
const __m128d Cloud::getz1_pd(const cloud_index i) const
{
	return _mm_load_pd(z + i); // z
}

const __m128d Cloud::getz2_pd(const cloud_index i) const
{
    return zCache[i/2]; // z + p1/2
}

const __m128d Cloud::getz3_pd(const cloud_index i) const
{
    return zCache[i/2]; //z + p2/2
}

const __m128d Cloud::getz4_pd(const cloud_index i) const
{
    return zCache[i/2]; // z + p3
}

const __m128d Cloud::getz1r_pd(const cloud_index i) const
{
    return _mm_loadr_pd(z + i);
}

const __m128d Cloud::getz2r_pd(const cloud_index i) const
{
    const cloud_index j = i/2;
    return _mm_shuffle_pd(zCache[j], zCache[j], _MM_SHUFFLE2(0, 1));
}

const __m128d Cloud::getz3r_pd(const cloud_index i) const
{
    const cloud_index j = i/2;
    return _mm_shuffle_pd(zCache[j], zCache[j], _MM_SHUFFLE2(0, 1));
}

const __m128d Cloud::getz4r_pd(const cloud_index i) const
{
    const cloud_index j = i/2;
    return _mm_shuffle_pd(zCache[j], zCache[j], _MM_SHUFFLE2(0, 1));
}

//x-velocity helper functions ------------------------------------------------
const __m128d Cloud::getVx1_pd(const cloud_index i) const
{
    return _mm_load_pd(Vx + i); // Vx
}

const __m128d Cloud::getVx2_pd(const cloud_index i) const
{
    return VxCache[i/2]; // Vx + k1/2
}

const __m128d Cloud::getVx3_pd(const cloud_index i) const
{
    return VxCache[i/2]; // Vx + k2/2
}

const __m128d Cloud::getVx4_pd(const cloud_index i) const
{
    return VxCache[i/2]; // Vx + k3
}

//y-velocity helper functions ------------------------------------------------
const __m128d Cloud::getVy1_pd(const cloud_index i) const
{
    return _mm_load_pd(Vy + i); // Vy
}

const __m128d Cloud::getVy2_pd(const cloud_index i) const
{
    return VyCache[i/2]; // Vy + k1/2
}

const __m128d Cloud::getVy3_pd(const cloud_index i) const
{
    return VyCache[i/2]; // Vy + k2/2
}

const __m128d Cloud::getVy4_pd(const cloud_index i) const
{
    return VyCache[i/2]; // Vy + k3
}

//z-velocity helper functions ------------------------------------------------
const __m128d Cloud::getVz1_pd(const cloud_index i) const
{
    return _mm_load_pd(Vz + i); // Vz
}

const __m128d Cloud::getVz2_pd(const cloud_index i) const
{
    return VzCache[i/2]; // Vz + k1/2
}

const __m128d Cloud::getVz3_pd(const cloud_index i) const
{
    return VzCache[i/2]; // Vz + k2/2
}

const __m128d Cloud::getVz4_pd(const cloud_index i) const
{
    return VzCache[i/2]; // Vz + k3
}
