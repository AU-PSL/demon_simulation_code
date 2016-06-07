/**
* @file  Cloud.cpp
* @class Cloud Cloud.h
*
* @brief Defines the physical parameters of the dust cloud
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#include "Cloud.h"
#include <cmath>
#include <sstream>

const double Cloud::electronCharge = -1.602E-19;
const double Cloud::epsilon0 = 8.8541878E-12;


/**
* @brief Constructor for the cloud class
* @param[in] numPar The number of particles
**/
Cloud::Cloud(const cloud_index numPar) :
	n(numPar),
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

/**
* @brief Destructor for the cloud class
**/
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

/**
* @brief Sets the charges for the dust particles as a gaussian distribution.
*		   
*
* @param[in] qMean  The average charge in Coulombs
* @param[in] qSigma The standard deviation for the charge in Coulombs
**/
inline void Cloud::initCharge(const double qMean, const double qSigma) {
	std::normal_distribution<double> dist(qMean, qSigma);
	for (cloud_index i = 0; i < n; i++)
		charge[i] = rands.gaussian(dist)*electronCharge;
}

/**
* @brief Sets the masses for the dust particles as a gaussian distribution.
*		   
*
* @param[in] rMean  The average radius in meters
* @param[in] rSigma The standard deviation for the radius in meters
**/
inline void Cloud::initMass(const double rMean, const double rSigma) {
	const double particleMassConstant = (4.0/3.0)*M_PI*dustParticleMassDensity;
	std::normal_distribution<double> dist(rMean, rSigma);
	for (cloud_index i = 0; i < n; i++) {
		const double r = rands.gaussian(dist);
		mass[i] = particleMassConstant*r*r*r;
	}
}

/**
* @brief Sets the initial position for each dust particle.
*		   
*
* @param[in] numParticles  The total number of particles
* @param[in] row_x_particles The number of rows in the x-direction
* @param[in] row_y_particles The number of rows in the y-direction
* @param[in] rMean  The average radius in meters
* @param[in] rSigma The standard deviation for the radius in meters
* @param[in] qMean  The average charge in Coulombs
* @param[in] qSigma The standard deviation for the charge in Coulombs
**/
Cloud * const Cloud::initializeGrid(const cloud_index numParticles,
									cloud_index row_x_particles,
									cloud_index row_y_particles,
									const double rMean, const double rSigma,
                                    const double qMean, const double qSigma) {

	Cloud * const cloud = new Cloud(numParticles);

	const cloud_index sqrtNumPar = (cloud_index)floor(sqrt(numParticles));

	//If only numParticles defined, creates as square of a grid as possible
	if(row_y_particles==0) {
		row_x_particles = sqrtNumPar;
		while(numParticles % row_x_particles != 0) {
			++row_x_particles;
		}
		row_y_particles = numParticles / row_x_particles;
	}

	const double cloudHalfSizeX = ((double)row_x_particles-1.0)/2.0*interParticleSpacing;
	const double cloudHalfSizeY = ((double)row_y_particles-1.0)/2.0*interParticleSpacing;

	cloud->initCharge(qMean, qSigma);
	cloud->initMass(rMean, rSigma);
	
	//Put particles into grid
    for(int i = 0; i < row_x_particles; i++) {
    	for(int j = 0; j < row_y_particles; j++) {
      		cloud->x[i*row_y_particles+j] = cloudHalfSizeX - (i*interParticleSpacing) + justX;
      		cloud->y[i*row_y_particles+j] = cloudHalfSizeY - (j*interParticleSpacing) + justY;
	    }
	}

	//Set particle velocities
    BEGIN_PARALLEL_FOR(l, e, numParticles, 1, static)
    	cloud->Vx[l] = velX;
		cloud->Vy[l] = velY;
    END_PARALLEL_FOR
	
	return cloud;
}


/**
* @brief Generates a cloud using the last time step of the specified file
*		   
*
* @param[in]  file        The name of the fits file
* @param[out] error       The error code (if any) that was produced when opening the file
* @param[in]  currentTime ??UNKNOWN??
**/
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

/**
* @brief Sets up and writes the intial timestep data to a the specified fits file
*		   
*
* @param[in]  file        The name of the fits file
* @param[out] error       The error code (if any) that was produced when opening the file
**/
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


/**
* @brief Appends the current cloud data to a fits file
* @details Appends the current cloud data to a fits file. This method requires that the
*          current HDU be "TIME_STEP". All data is flushed to the file in case execution
*          is interrupted.	   
*
* @param[in]  file        The name of the fits file
* @param[out] error       The error code (if any) that was produced when opening the file
* @param[in]  currentTime ??UNKNOWN??
**/
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

/**
* @brief x-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getx1_pd(const cloud_index i) const {
	return load_pd(x + i); // x
}

/**
* @brief x-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getx2_pd(const cloud_index i) const {
	return xCache[i/DOUBLE_STRIDE]; // x + l1/2
}

/**
* @brief x-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getx3_pd(const cloud_index i) const {
	return xCache[i/DOUBLE_STRIDE]; // x + l2/2
}

/**
* @brief x-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getx4_pd(const cloud_index i) const {
	return xCache[i/DOUBLE_STRIDE]; // x + l3
}

/**
* @brief x-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getx1r_pd(const cloud_index i) const {
	return _mm_loadr_pd(x + i);
}

/**
* @brief x-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getx2r_pd(const cloud_index i) const {
	const cloud_index j = i/DOUBLE_STRIDE;
	return _mm_shuffle_pd(xCache[j], xCache[j], _MM_SHUFFLE2(0, 1));
}

/**
* @brief x-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getx3r_pd(const cloud_index i) const {
	const cloud_index j = i/DOUBLE_STRIDE;
	return _mm_shuffle_pd(xCache[j], xCache[j], _MM_SHUFFLE2(0, 1));
}

/**
* @brief x-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getx4r_pd(const cloud_index i) const {
	const cloud_index j = i/DOUBLE_STRIDE;
	return _mm_shuffle_pd(xCache[j], xCache[j], _MM_SHUFFLE2(0, 1));
}

// Y position helper functions -------------------------------------------------

/**
* @brief y-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::gety1_pd(const cloud_index i) const {
	return load_pd(y + i); // y
}

/**
* @brief y-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::gety2_pd(const cloud_index i) const {
	return yCache[i/DOUBLE_STRIDE]; // y + n1/2
}

/**
* @brief y-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::gety3_pd(const cloud_index i) const {
	return yCache[i/DOUBLE_STRIDE]; // y + n2/2
}

/**
* @brief y-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::gety4_pd(const cloud_index i) const {
	return yCache[i/DOUBLE_STRIDE]; // y + n3
}

/**
* @brief y-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::gety1r_pd(const cloud_index i) const {
	return _mm_loadr_pd(y + i);
}

/**
* @brief y-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::gety2r_pd(const cloud_index i) const {
	const cloud_index j = i/DOUBLE_STRIDE;
	return _mm_shuffle_pd(yCache[j], yCache[j], _MM_SHUFFLE2(0, 1));
}

/**
* @brief y-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::gety3r_pd(const cloud_index i) const {
	const cloud_index j = i/DOUBLE_STRIDE;
	return _mm_shuffle_pd(yCache[j], yCache[j], _MM_SHUFFLE2(0, 1));
}

/**
* @brief y-position helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::gety4r_pd(const cloud_index i) const {
	const cloud_index j = i/DOUBLE_STRIDE;
	return _mm_shuffle_pd(yCache[j], yCache[j], _MM_SHUFFLE2(0, 1));
}

// Vx position helper functions ------------------------------------------------

/**
* @brief x-velocity helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getVx1_pd(const cloud_index i) const {
	return load_pd(Vx + i); // Vx
}

/**
* @brief x-velocity helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getVx2_pd(const cloud_index i) const {
	return VxCache[i/DOUBLE_STRIDE]; // Vx + k1/2
}

/**
* @brief x-velocity helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getVx3_pd(const cloud_index i) const {
	return VxCache[i/DOUBLE_STRIDE]; // Vx + k2/2
}

/**
* @brief x-velocity helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getVx4_pd(const cloud_index i) const {
	return VxCache[i/DOUBLE_STRIDE]; // Vx + k3
}

// Vy position helper functions ------------------------------------------------

/**
* @brief y-velocity helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getVy1_pd(const cloud_index i) const {
	return load_pd(Vy + i); // Vy
}

/**
* @brief y-velocity helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getVy2_pd(const cloud_index i) const {
	return VyCache[i/DOUBLE_STRIDE]; // Vy + m1/2
}

/**
* @brief y-velocity helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getVy3_pd(const cloud_index i) const {
	return VyCache[i/DOUBLE_STRIDE]; // Vy + m2/2
}

/**
* @brief y-velocity helper method for 4th-order Runge-Kutta substep
* @param[in] i ??UNKNOWN??
**/
const doubleV Cloud::getVy4_pd(const cloud_index i) const {
	return VyCache[i/DOUBLE_STRIDE]; // Vy + m3
}
