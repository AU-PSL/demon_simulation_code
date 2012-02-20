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

const double Cloud::interParticleSpacing = 0.0003;
const double Cloud::electronCharge = -1.602E-19;
const double Cloud::epsilon0 = 8.8541878E-12;
const double Cloud::dustParticleMassDensity = 2200.0;

Cloud::Cloud(const cloud_index numPar) 
: n(numPar),
x(new double[n]), y(new double[n]), Vx(new double[n]), Vy(new double[n]), 
charge(new double[n]), mass(new double[n]),
k1(new double[n]), k2(new double[n]), k3(new double[n]), k4(new double[n]),
l1(new double[n]), l2(new double[n]), l3(new double[n]), l4(new double[n]),
m1(new double[n]), m2(new double[n]), m3(new double[n]), m4(new double[n]),
n1(new double[n]), n2(new double[n]), n3(new double[n]), n4(new double[n]),
forceX(new double[n]), forceY(new double[n]),
xCache(new doubleV[n/DOUBLE_STRIDE]), yCache(new doubleV[n/DOUBLE_STRIDE]), 
VxCache(new doubleV[n/DOUBLE_STRIDE]), VyCache(new doubleV[n/DOUBLE_STRIDE]) {
#ifdef _OPENMP
	omp_set_num_threads(omp_get_num_procs());
#endif
}

Cloud::~Cloud() {
	delete[] x; delete[] y; delete[] Vx; delete[] Vy;
	delete[] charge; delete[] mass; 
	delete[] k1; delete[] k2; delete[] k3; delete[] k4;
	delete[] l1; delete[] l2; delete[] l3; delete[] l4;
	delete[] m1; delete[] m2; delete[] m3; delete[] m4;
	delete[] n1; delete[] n2; delete[] n3; delete[] n4;
	delete[] forceX; delete[] forceY;
	delete[] xCache; delete[] yCache; 
	delete[] VxCache; delete[] VyCache;
}

// Particle charges are set as a guassian distribution. To set a uniform 
// distribution the charge sigma should be set to zero.
inline void Cloud::initCharge(const double qMean, const double qSigma) {
	std::normal_distribution<double> dist(qMean, qSigma);
	for (cloud_index i = 0; i < n; i++)
		charge[i] = rands.guassian(dist)*electronCharge;
}

// Particle sizes are set as a guassian distribution. To set a uniform 
// distribution the radius sigma should be set to zero.
inline void Cloud::initMass(const double rMean, const double rSigma) {
	const double particleMassConsant = (4.0/3.0)*M_PI*dustParticleMassDensity;
	std::normal_distribution<double> dist(rMean, rSigma);
	for (cloud_index i = 0; i < n; i++) {
		const double r = rands.guassian(dist);
		mass[i] = particleMassConsant*r*r*r;
	}
}

// Generates a cloud on a square grid with zero inital velocity.
Cloud * const Cloud::initializeGrid(const cloud_index numParticles, 
									const double rMean, const double rSigma,
                                    const double qMean, const double qSigma) {
	Cloud * const cloud = new Cloud(numParticles);

	const cloud_index sqrtNumPar = (cloud_index)floor(sqrt(numParticles));
	
	// For even numbers of partciles on a row center the row over the origin.
	const double cloudHalfSize = (double)sqrtNumPar/2.0*interParticleSpacing 
		- ((sqrtNumPar%2) ? 0.0 : interParticleSpacing/2.0);

	cloud->initCharge(qMean, qSigma);
	cloud->initMass(rMean, rSigma);
    BEGIN_PARALLEL_FOR(i, e, numParticles, 1, static)
    	cloud->x[i] = cloudHalfSize - (double)(i%sqrtNumPar)*interParticleSpacing;
	    cloud->y[i] = cloudHalfSize - (double)(i%sqrtNumPar)*interParticleSpacing;
		cloud->Vx[i] = cloud->Vy[i] = 0.0;
    END_PARALLEL_FOR
	return cloud;
}

// Generates a cloud using the last time step of the specified file.
Cloud * const Cloud::initializeFromFile(fitsfile * const file, int &error, 
										double * const currentTime) {
	int anyNull = 0;
	long numParticles = 0;
	long numTimeSteps = 0;

	// move to CLOUD HDU:
	if (!error)
		fits_movnam_hdu(file, BINARY_TBL, const_cast<char *> ("CLOUD"), 0, &error);

	// get number of particles:
	if (!error)
		fits_get_num_rows(file, &numParticles, &error);

	// create cloud:
	Cloud * const cloud = new Cloud((cloud_index)numParticles);

	// read mass information:
	if (!error) {
		// file, column #, starting row, first element, num elements, mass array, pointless pointer, error
		fits_read_col_dbl(file, 1, 1, 1, numParticles, 0.0, cloud->mass, &anyNull, &error);
		fits_read_col_dbl(file, 2, 1, 1, numParticles, 0.0, cloud->charge, &anyNull, &error);
	}

	// move to TIME_STEP HDU:
	if (!error)
		fits_movnam_hdu(file, BINARY_TBL, const_cast<char *> ("TIME_STEP"), 0, &error);

	// get number of time steps:
	if (!error)
		fits_get_num_rows(file, &numTimeSteps, &error);

	if (!error) {
		if (currentTime)
			fits_read_col_dbl(file, 1, numTimeSteps, 1, 1, 0.0, currentTime, &anyNull, &error);

		fits_read_col_dbl(file, 2, numTimeSteps, 1, numParticles, 0.0, cloud->x, &anyNull, &error);
		fits_read_col_dbl(file, 3, numTimeSteps, 1, numParticles, 0.0, cloud->y, &anyNull, &error);
		fits_read_col_dbl(file, 4, numTimeSteps, 1, numParticles, 0.0, cloud->Vx, &anyNull, &error);
		fits_read_col_dbl(file, 5, numTimeSteps, 1, numParticles, 0.0, cloud->Vy, &anyNull, &error);
	}

	return cloud;
}

// Sets up and writes the intial timestep data to a the specified fitsfile.
void Cloud::writeCloudSetup(fitsfile * const file, int &error) const {
	// format number of elements of type double as string, e.g. 1024D
	std::stringstream numStream;
	numStream << n << "D";
	const std::string numString = numStream.str();

	char *ttypeCloud[] = {const_cast<char *> ("MASS"), const_cast<char *> ("CHARGE")};
	char *tformCloud[] = {const_cast<char *> ("D"), const_cast<char *> ("D")};
	char *tunitCloud[] = {const_cast<char *> ("kg"), const_cast<char *> ("C")};	

	char *ttypeRun[] = {const_cast<char *> ("TIME"),
		const_cast<char *> ("X_POSITION"), const_cast<char *> ("Y_POSITION"), 
		const_cast<char *> ("X_VELOCITY"), const_cast<char *> ("Y_VELOCITY")};
	char *tformRun[] = {const_cast<char *> ("D"), 
		const_cast<char *> (numString.c_str()), const_cast<char *> (numString.c_str()), 
		const_cast<char *> (numString.c_str()), const_cast<char *> (numString.c_str())};
	char *tunitRun[] = {const_cast<char *> ("s"),
		const_cast<char *> ("m"), const_cast<char *> ("m"), 
		const_cast<char *> ("m/s"), const_cast<char *> ("m/s")};

	// write mass:
	if (!error)
		// file, storage type, num rows, num columns, ...
		fits_create_tbl(file, BINARY_TBL, (LONGLONG)n, 2, ttypeCloud, tformCloud, tunitCloud, "CLOUD", &error);	
	if (!error) {
		// file, column #, starting row, first element, num elements, mass array, error
		fits_write_col_dbl(file, 1, 1, 1, (LONGLONG)n, mass, &error);
		fits_write_col_dbl(file, 2, 1, 1, (LONGLONG)n, charge, &error);
	}

	// write position and velocity:
	if (!error)
		fits_create_tbl(file, BINARY_TBL, 0, 5, ttypeRun, tformRun, tunitRun, "TIME_STEP", &error);
		// n.b. num rows automatically incremented.
		// Increment from 0 as opposed to preallocating to ensure
		// proper output in the event of program interruption.
	if (!error) {
		double time = 0.0;
		fits_write_col_dbl(file, 1, 1, 1, 1, &time, &error);
		fits_write_col_dbl(file, 2, 1, 1, (LONGLONG)n, x, &error);
		fits_write_col_dbl(file, 3, 1, 1, (LONGLONG)n, y, &error);
		fits_write_col_dbl(file, 4, 1, 1, (LONGLONG)n, Vx, &error);
		fits_write_col_dbl(file, 5, 1, 1, (LONGLONG)n, Vy, &error);
	}

	// write buffer, close file, reopen at same point:
	fits_flush_file(file, &error);
}

// Appends the current cloud data to a fits file. This method requires that the
// current HDU be "TIME_STEP". All data is flushed to the file incase exicution
// is interupted.
void Cloud::writeTimeStep(fitsfile * const file, int &error, double currentTime) const {
	if (!error) {
		long numRows = 0;
		fits_get_num_rows(file, &numRows, &error);
		fits_write_col_dbl(file, 1, ++numRows, 1, 1, &currentTime, &error);
		fits_write_col_dbl(file, 2, numRows, 1, (LONGLONG)n, x, &error);
		fits_write_col_dbl(file, 3, numRows, 1, (LONGLONG)n, y, &error);
		fits_write_col_dbl(file, 4, numRows, 1, (LONGLONG)n, Vx, &error);
		fits_write_col_dbl(file, 5, numRows, 1, (LONGLONG)n, Vy, &error);
	}

	// write buffer, close file, reopen at same point:
	fits_flush_file(file, &error);
}

// 4th order Runge-Kutta subsetp helper methods. 

// X position helper functions -------------------------------------------------
const doubleV Cloud::getx1_pd(const cloud_index i) const {
	return load_pd(x + i); // x
}

const doubleV Cloud::getx2_pd(const cloud_index i) const {
	return xCache[i/DOUBLE_STRIDE]; // x + l1/2
}

const doubleV Cloud::getx3_pd(const cloud_index i) const {
	return xCache[i/DOUBLE_STRIDE]; // x + l2/2
}

const doubleV Cloud::getx4_pd(const cloud_index i) const {
	return xCache[i/DOUBLE_STRIDE]; // x + l3
}

const doubleV Cloud::getx1r_pd(const cloud_index i) const {
	return _mm_loadr_pd(x + i);
}

const doubleV Cloud::getx2r_pd(const cloud_index i) const {
	const cloud_index j = i/DOUBLE_STRIDE;
	return _mm_shuffle_pd(xCache[j], xCache[j], _MM_SHUFFLE2(0, 1));
}

const doubleV Cloud::getx3r_pd(const cloud_index i) const {
	const cloud_index j = i/DOUBLE_STRIDE;
	return _mm_shuffle_pd(xCache[j], xCache[j], _MM_SHUFFLE2(0, 1));
}

const doubleV Cloud::getx4r_pd(const cloud_index i) const {
	const cloud_index j = i/DOUBLE_STRIDE;
	return _mm_shuffle_pd(xCache[j], xCache[j], _MM_SHUFFLE2(0, 1));
}

// Y position helper functions -------------------------------------------------
const doubleV Cloud::gety1_pd(const cloud_index i) const {
	return load_pd(y + i); // y
}

const doubleV Cloud::gety2_pd(const cloud_index i) const {
	return yCache[i/DOUBLE_STRIDE]; // y + n1/2
}

const doubleV Cloud::gety3_pd(const cloud_index i) const {
	return yCache[i/DOUBLE_STRIDE]; // y + n2/2
}

const doubleV Cloud::gety4_pd(const cloud_index i) const {
	return yCache[i/DOUBLE_STRIDE]; // y + n3
}

const doubleV Cloud::gety1r_pd(const cloud_index i) const 
{
	return _mm_loadr_pd(y + i);
}

const doubleV Cloud::gety2r_pd(const cloud_index i) const {
	const cloud_index j = i/DOUBLE_STRIDE;
	return _mm_shuffle_pd(yCache[j], yCache[j], _MM_SHUFFLE2(0, 1));
}

const doubleV Cloud::gety3r_pd(const cloud_index i) const {
	const cloud_index j = i/DOUBLE_STRIDE;
	return _mm_shuffle_pd(yCache[j], yCache[j], _MM_SHUFFLE2(0, 1));
}

const doubleV Cloud::gety4r_pd(const cloud_index i) const {
	const cloud_index j = i/DOUBLE_STRIDE;
	return _mm_shuffle_pd(yCache[j], yCache[j], _MM_SHUFFLE2(0, 1));
}

// Vx position helper functions ------------------------------------------------
const doubleV Cloud::getVx1_pd(const cloud_index i) const {
	return load_pd(Vx + i); // Vx
}

const doubleV Cloud::getVx2_pd(const cloud_index i) const {
	return VxCache[i/DOUBLE_STRIDE]; // Vx + k1/2
}

const doubleV Cloud::getVx3_pd(const cloud_index i) const {
	return VxCache[i/DOUBLE_STRIDE]; // Vx + k2/2
}

const doubleV Cloud::getVx4_pd(const cloud_index i) const {
	return VxCache[i/DOUBLE_STRIDE]; // Vx + k3
}

// Vy position helper functions ------------------------------------------------
const doubleV Cloud::getVy1_pd(const cloud_index i) const {
	return load_pd(Vy + i); // Vy
}

const doubleV Cloud::getVy2_pd(const cloud_index i) const {
	return VyCache[i/DOUBLE_STRIDE]; // Vy + m1/2
}

const doubleV Cloud::getVy3_pd(const cloud_index i) const {
	return VyCache[i/DOUBLE_STRIDE]; // Vy + m2/2
}

const doubleV Cloud::getVy4_pd(const cloud_index i) const {
	return VyCache[i/DOUBLE_STRIDE]; // Vy + m3
}
