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
#include <iostream>
#include <sstream>

using namespace std;

Cloud::Cloud(unsigned int numPar, double sizeOfCloud) : n(numPar), cloudSize(sizeOfCloud),
k1(new double[n]), k2(new double[n]), k3(new double[n]), k4(new double[n]),
l1(new double[n]), l2(new double[n]), l3(new double[n]), l4(new double[n]),
m1(new double[n]), m2(new double[n]), m3(new double[n]), m4(new double[n]),
n1(new double[n]), n2(new double[n]), n3(new double[n]), n4(new double[n]),
x(new double[n]), y(new double[n]), Vx(new double[n]), Vy(new double[n]), 
charge(new double[n]), mass(new double[n]), 
forceX(new double[n]), forceY(new double[n]) {}

Cloud::~Cloud() 
{
	delete[] k1; delete[] k2; delete[] k3; delete[] k4;
	delete[] l1; delete[] l2; delete[] l3; delete[] l4;
	delete[] m1; delete[] m2; delete[] m3; delete[] m4;
	delete[] n1; delete[] n2; delete[] n3; delete[] n4;
	delete[] x; delete[] y; delete[] Vx; delete[] Vy;
	delete[] charge; delete[] mass; 
	delete[] forceX; delete[] forceY;
}

inline void Cloud::setPosition(const unsigned int index)
{
	double tempR = arc4random();
	double tempTheta = arc4random();
	double radius = (tempR/RAND_MAX)*cloudSize;
	double theta = (tempTheta/RAND_MAX)*M_2_PI;

	x[index] = radius*cos(theta);
	y[index] = radius*sin(theta);
}

inline void Cloud::setPosition(const unsigned int index, const double xVal, const double yVal)
{
	x[index] = xVal;
	y[index] = yVal;
}

inline void Cloud::setVelocity(const unsigned int index)
{
	double initialVel = 0.0;

	double temp = rand();
	double theta = (temp/RAND_MAX)*(M_2_PI);

	Vx[index] = initialVel*sin(theta);
	Vy[index] = initialVel*cos(theta);
}

inline void Cloud::setCharge(const unsigned int index)
{
	charge[index] = (rand()%201 + 5900)*1.6E-19;
}

inline void Cloud::setMass(const unsigned int index)
{
	const double radius = 1.45E-6;
	const double particleDensity = 2200.0;
	mass[index] = (4.0/3.0)*M_PI*(pow(radius,3))*particleDensity;
}

Cloud *Cloud::initializeNew(const unsigned int numParticles, const double cloudSize)
{
	cout << "\nGenerating original data.\n\n";

	Cloud* cloud = new Cloud(numParticles, cloudSize);

	//initialize dust cloud:
	for(unsigned int i = 0; i < numParticles; i++)
	{
		cloud->setPosition(i);
		cloud->setVelocity(i);
		cloud->setCharge(i);
		cloud->setMass(i);
	}

	return cloud;
}

Cloud *Cloud::initializeGrid(const unsigned int numParticles, const double cloudSize)
{
	cout << "\nInitializing grid.\n\n";

	Cloud* cloud = new Cloud(numParticles, cloudSize);

	const double sqrtNumPar = floor(sqrt(numParticles));
	const double gridUnit = 2.0*cloudSize/sqrtNumPar; //number of particles per row/column
	double tempPosX = cloudSize; //position of first particle
	double tempPosY = cloudSize;

	//initialize dust cloud:
	for(unsigned int i = 0; i < numParticles; i++)
	{
		cloud->setPosition(i, tempPosX, tempPosY);
		cloud->setVelocity(i);
		cloud->setCharge(i);
		cloud->setMass(i);

		tempPosX -= gridUnit;
		if(tempPosX <= -1*cloudSize) //end of row
		{
			tempPosX = cloudSize; //reset
			tempPosY -= gridUnit; //move to next row
		}
	}

	cout << "\nInitialization complete.\n";

	return cloud;
}

Cloud *Cloud::initializeFromFile(fitsfile *file, int *error, double *currentTime)
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
		fits_read_col_dbl(file, 4, numTimeSteps, 1, numParticles, 0.0, cloud->Vx, anyNull, error);
		fits_read_col_dbl(file, 5, numTimeSteps, 1, numParticles, 0.0, cloud->Vy, anyNull, error);
	}

	return cloud;
}

void Cloud::writeCloudSetup(fitsfile *file, int *error) const
{
	//format number of elements of type double as string, e.g. 1024D
	stringstream numStream;
	numStream << n << "D";
	string numString = numStream.str();

	char *ttypeCloud[] = {const_cast<char *> ("CHARGE"), const_cast<char *> ("MASS")};
	char *tformCloud[] = {const_cast<char *> ("D"), const_cast<char *> ("D")};
	char *tunitCloud[] = {const_cast<char *> ("C"), const_cast<char *> ("kg")};	

	char *ttypeRun[] = {const_cast<char *> ("TIME"),
		const_cast<char *> ("X_POSITION"), const_cast<char *> ("Y_POSITION"), 
		const_cast<char *> ("X_VELOCITY"), const_cast<char *> ("Y_VELOCITY")};
	char *tformRun[] = {const_cast<char *> ("D"), 
		const_cast<char *> (numString.c_str()), const_cast<char *> (numString.c_str()), 
		const_cast<char *> (numString.c_str()), const_cast<char *> (numString.c_str())};
	char *tunitRun[] = {const_cast<char *> ("s"),
		const_cast<char *> ("m"), const_cast<char *> ("m"), 
		const_cast<char *> ("m/s"), const_cast<char *> ("m/s")};

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
		fits_write_col_dbl(file, 4, 1, 1, n, Vx, error);
		fits_write_col_dbl(file, 5, 1, 1, n, Vy, error);
	}

	//write buffer, close file, reopen at same point:
	fits_flush_file(file, error);
}

void Cloud::writeTimeStep(fitsfile *file, int *error, double currentTime) const
{
	if (!*error)
	{
		long numRows = 0;
		fits_get_num_rows(file, &numRows, error);
		fits_write_col_dbl(file, 1, ++numRows, 1, 1, &currentTime, error);
		fits_write_col_dbl(file, 2, numRows, 1, n, x, error);
		fits_write_col_dbl(file, 3, numRows, 1, n, y, error);
		fits_write_col_dbl(file, 4, numRows, 1, n, Vx, error);
		fits_write_col_dbl(file, 5, numRows, 1, n, Vy, error);
	}

	//write buffer, close file, reopen at same point:
	fits_flush_file(file, error);
}
