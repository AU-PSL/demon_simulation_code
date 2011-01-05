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

Cloud::Cloud(unsigned int numPar, double sizeOfCloud) : n(numPar), cloudSize(sizeOfCloud),
k1(new double[n]), k2(new double[n]), k3(new double[n]), k4(new double[n]),
l1(new double[n]), l2(new double[n]), l3(new double[n]), l4(new double[n]),
m1(new double[n]), m2(new double[n]), m3(new double[n]), m4(new double[n]),
n1(new double[n]), n2(new double[n]), n3(new double[n]), n4(new double[n]),
o1(new double[n]), o2(new double[n]), o3(new double[n]), o4(new double[n]),
p1(new double[n]), p2(new double[n]), p3(new double[n]), p4(new double[n]),
x(new double[n]), y(new double[n]), z(new double[n]),
Vx(new double[n]), Vy(new double[n]), Vz(new double[n]),
charge(new double[n]), mass(new double[n]), 
forceX(new double[n]), forceY(new double[n]), forceZ(new double[n]) {}

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
}

inline void Cloud::setPosition(const unsigned int index, const double xVal, const double yVal, const double zVal)
{
	x[index] = xVal;
	y[index] = yVal;
	z[index] = zVal;
}

inline void Cloud::setVelocity(const unsigned int index)
{
	Vx[index] = 0.0;
	Vy[index] = 0.0;
	Vz[index] = 0.0;
}

inline void Cloud::setCharge(const unsigned int index)
{
	charge[index] = (rand()%201 + 5900)*1.6E-19;
}

inline void Cloud::setMass(const unsigned int index)
{
	const double radius = 1.45E-6;
	const double particleDensity = 2200.0;
	mass[index] = (4.0/3.0)*M_PI*radius*radius*radius*particleDensity;
}

Cloud * const Cloud::initializeGrid(const unsigned int numParticles, const double cloudSize, const bool make3D)
{
	Cloud * const cloud = new Cloud(numParticles, cloudSize);

	//seed rand function with time(NULL):
	srand((int)time(NULL));		//needed for Cloud::setCharge

	double tempPosX = cloudSize;	//position of first particle
	double tempPosY = cloudSize;

	if(!make3D)
	{
		const double sqrtNumPar = floor(sqrt(numParticles));
		const double gridUnit = 2.0*cloudSize/sqrtNumPar;	//number of particles per row/column
    
		//initialize dust cloud:
		for(unsigned int i = 0; i < numParticles; i++)
		{
			cloud->setPosition(i, tempPosX, tempPosY, 0.0);
			cloud->setVelocity(i);
			cloud->setCharge(i);
			cloud->setMass(i);

			tempPosX -= gridUnit;
			if(tempPosX <= -cloudSize) //end of row
			{
				tempPosX = cloudSize; //reset
				tempPosY -= gridUnit; //move to next row
			}
		}
	}
	else
	{
		double tempPosZ = cloudSize;	

		const double cubeNumPar = floor(pow(numParticles,(double)(1/3)));
		const double gridUnit = 2.0*cloudSize/cubeNumPar;

		//initialize dust cloud:
		for(unsigned int i = 0; i < numParticles; i++)
		{
			cloud->setPosition(i, tempPosX, tempPosY, tempPosZ);
			cloud->setVelocity(i);
			cloud->setCharge(i);
			cloud->setMass(i);

			tempPosX -= gridUnit;
			if(tempPosX <= -cloudSize) //end of row
			{
				tempPosX = cloudSize; //reset
				tempPosY -= gridUnit; //move to next row

				if(tempPosY <= -cloudSize) //end of plane
				{
					tempPosY = cloudSize; //reset
					tempPosZ -= gridUnit; //move to next plane
				}
			}
		}
	}

	return cloud;
}

//FIXME: Allow initialization from files made with 2D code.
Cloud * const Cloud::initializeFromFile(fitsfile * const file, int * const error, double * const currentTime)
{
	int *anyNull = NULL;
	long numParticles;
	long numTimeSteps;

	//move to CLOUD HDU:
	if(!*error)
		fits_movnam_hdu(file, BINARY_TBL, const_cast<char *> ("CLOUD"), 0, error);

	//get number of particles:
	if(!*error)
		fits_get_num_rows(file, &numParticles, error);

	//create cloud:
	Cloud* cloud = new Cloud((unsigned int)numParticles, 0.0);	//cloudSize not used in this case, so set to zero

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
		if (currentTime)
			fits_read_col_dbl(file, 1, numTimeSteps, 1, 1, 0.0, currentTime, anyNull, error);

		fits_read_col_dbl(file, 2, numTimeSteps, 1, numParticles, 0.0, cloud->x, anyNull, error);
		fits_read_col_dbl(file, 3, numTimeSteps, 1, numParticles, 0.0, cloud->y, anyNull, error);
		fits_read_col_dbl(file, 3, numTimeSteps, 1, numParticles, 0.0, cloud->z, anyNull, error);
		fits_read_col_dbl(file, 4, numTimeSteps, 1, numParticles, 0.0, cloud->Vx, anyNull, error);
		fits_read_col_dbl(file, 5, numTimeSteps, 1, numParticles, 0.0, cloud->Vy, anyNull, error);
		fits_read_col_dbl(file, 5, numTimeSteps, 1, numParticles, 0.0, cloud->Vz, anyNull, error);
	}

	return cloud;
}

void Cloud::writeCloudSetup(fitsfile * const file, int * const error) const
{
	//format number of elements of type double as string, e.g. 1024D
	stringstream numStream;
	numStream << n << "D";
	const string numString = numStream.str();

	char *ttypeCloud[] = {const_cast<char *> ("CHARGE"), const_cast<char *> ("MASS")};
	char *tformCloud[] = {const_cast<char *> ("D"), const_cast<char *> ("D")};
	char *tunitCloud[] = {const_cast<char *> ("C"), const_cast<char *> ("kg")};	

	char *ttypeRun[] = {const_cast<char *> ("TIME"),
		const_cast<char *> ("X_POSITION"), const_cast<char *> ("Y_POSITION"), const_cast<char *> ("Z_POSITION"), 
		const_cast<char *> ("X_VELOCITY"), const_cast<char *> ("Y_VELOCITY"), const_cast<char *> ("Z_VELOCITY")};
	char *tformRun[] = {const_cast<char *> ("D"), 
		const_cast<char *> (numString.c_str()), const_cast<char *> (numString.c_str()), const_cast<char *> (numString.c_str()), 
		const_cast<char *> (numString.c_str()), const_cast<char *> (numString.c_str()), const_cast<char *> (numString.c_str())};
	char *tunitRun[] = {const_cast<char *> ("s"),
		const_cast<char *> ("m"), const_cast<char *> ("m"), const_cast<char *> ("m"), 
		const_cast<char *> ("m/s"), const_cast<char *> ("m/s"), const_cast<char *> ("m/s")};

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
	if (!*error)
		fits_create_tbl(file, BINARY_TBL, 0, 5, ttypeRun, tformRun, tunitRun, "TIME_STEP", error);
		//n.b. num rows automatically incremented.
		// Increment from 0 as opposed to preallocating to ensure
		// proper output in the event of program interruption.
	if (!*error)
	{
		double time = 0.0;
		fits_write_col_dbl(file, 1, 1, 1, 1, &time, error);
		fits_write_col_dbl(file, 2, 1, 1, n, x, error);
		fits_write_col_dbl(file, 3, 1, 1, n, y, error);
		fits_write_col_dbl(file, 3, 1, 1, n, z, error);
		fits_write_col_dbl(file, 4, 1, 1, n, Vx, error);
		fits_write_col_dbl(file, 5, 1, 1, n, Vy, error);
		fits_write_col_dbl(file, 5, 1, 1, n, Vz, error);
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
		fits_write_col_dbl(file, 3, numRows, 1, n, z, error);
		fits_write_col_dbl(file, 4, numRows, 1, n, Vx, error);
		fits_write_col_dbl(file, 5, numRows, 1, n, Vy, error);
		fits_write_col_dbl(file, 5, numRows, 1, n, Vz, error);
	}

	//write buffer, close file, reopen at same point:
	fits_flush_file(file, error);
}

// 4th order Runge-Kutta substep helper methods:

// x-position helper functions -------------------------------------------------
const __m128d Cloud::getx1_pd(const unsigned int i) const	// x
{
    return _mm_load_pd(x + i);
}

const __m128d Cloud::getx2_pd(const unsigned int i) const	// x + l1/2
{
    return _mm_load_pd(x + i) + _mm_load_pd(l1 + i)/_mm_set1_pd(2.0);
}

const __m128d Cloud::getx3_pd(const unsigned int i) const	// x + l2/2
{
    return _mm_load_pd(x + i) + _mm_load_pd(l2 + i)/_mm_set1_pd(2.0);
}

const __m128d Cloud::getx4_pd(const unsigned int i) const	// x + l3
{
    return _mm_load_pd(x + i) + _mm_load_pd(l3 + i);
}

// y-position helper functions -------------------------------------------------
const __m128d Cloud::gety1_pd(const unsigned int i) const	// y
{
    return _mm_load_pd(y + i);
}

const __m128d Cloud::gety2_pd(const unsigned int i) const	// y + n1/2
{
    return _mm_load_pd(y + i) + _mm_load_pd(n1 + i)/_mm_set1_pd(2.0);
}

const __m128d Cloud::gety3_pd(const unsigned int i) const	// y + n2/2
{
    return _mm_load_pd(y + i) + _mm_load_pd(n2 + i)/_mm_set1_pd(2.0);
}

const __m128d Cloud::gety4_pd(const unsigned int i) const	// y + n3
{
    return _mm_load_pd(y + i) + _mm_load_pd(n3 + i);
}

// z-position helper functions -------------------------------------------------
const __m128d Cloud::getz1_pd(const unsigned int i) const	// z
{
    return _mm_load_pd(z + i);
}

const __m128d Cloud::getz2_pd(const unsigned int i) const	// z + p1/2
{
    return _mm_load_pd(z + i) + _mm_load_pd(p1 + i)/_mm_set1_pd(2.0);
}

const __m128d Cloud::getz3_pd(const unsigned int i) const	// z + p2/2
{
    return _mm_load_pd(z + i) + _mm_load_pd(p2 + i)/_mm_set1_pd(2.0);
}

const __m128d Cloud::getz4_pd(const unsigned int i) const	// z + p3
{
    return _mm_load_pd(z + i) + _mm_load_pd(p3 + i);
}

// x-velocity helper functions ------------------------------------------------
const __m128d Cloud::getVx1_pd(const unsigned int i) const	// Vx
{
    return _mm_load_pd(Vx + i);
}

const __m128d Cloud::getVx2_pd(const unsigned int i) const	// Vx + k1/2
{
    return _mm_load_pd(Vx + i) + _mm_load_pd(k1 + i)/_mm_set1_pd(2.0);
}

const __m128d Cloud::getVx3_pd(const unsigned int i) const	// Vx + k2/2
{
    return _mm_load_pd(Vx + i) + _mm_load_pd(k2 + i)/_mm_set1_pd(2.0);
}

const __m128d Cloud::getVx4_pd(const unsigned int i) const	// Vx + k3
{
    return _mm_load_pd(Vx + i) + _mm_load_pd(k3 + i);
}

// y-velocity helper functions ------------------------------------------------
const __m128d Cloud::getVy1_pd(const unsigned int i) const	// Vy
{
    return _mm_load_pd(Vy + i);
}

const __m128d Cloud::getVy2_pd(const unsigned int i) const	// Vy + m1/2
{
    return _mm_load_pd(Vy + i) + _mm_load_pd(m1 + i)/_mm_set1_pd(2.0);
}

const __m128d Cloud::getVy3_pd(const unsigned int i) const	// Vy + m2/2
{
    return _mm_load_pd(Vy + i) + _mm_load_pd(m2 + i)/_mm_set1_pd(2.0);
}

const __m128d Cloud::getVy4_pd(const unsigned int i) const	// Vy + m3
{
    return _mm_load_pd(Vy + i) + _mm_load_pd(m3 + i);
}

// z-velocity helper functions ------------------------------------------------
const __m128d Cloud::getVz1_pd(const unsigned int i) const	// Vz
{
    return _mm_load_pd(Vz + i);
}

const __m128d Cloud::getVz2_pd(const unsigned int i) const	// Vz + m1/2
{
    return _mm_load_pd(Vz + i) + _mm_load_pd(o1 + i)/_mm_set1_pd(2.0);
}

const __m128d Cloud::getVz3_pd(const unsigned int i) const	// Vz + m2/2
{
    return _mm_load_pd(Vz + i) + _mm_load_pd(o2 + i)/_mm_set1_pd(2.0);
}

const __m128d Cloud::getVz4_pd(const unsigned int i) const	// Vz + m3
{
    return _mm_load_pd(Vz + i) + _mm_load_pd(o3 + i);
}
