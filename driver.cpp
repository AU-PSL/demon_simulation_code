/**
* @file  driver.cpp
* @brief Parses input parameters and begins simulation
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#include "ConfinementForceVoid.h"
#include "DrivingForce.h"
#include "MagneticForce.h"
#include "RectConfinementForce.h"
#include "RotationalForce.h"
#include "Runge_Kutta4.h"
#include "ShieldedCoulombForce.h"
#include "ThermalForceLocalized.h"
#include "TimeVaryingDragForce.h"
#include "TimeVaryingThermalForce.h"
#include "ElectricForce.h"
#include "GravitationalForce.h"
#include "VertElectricForce.h"

#include <iostream>
#include <cstdarg>
#include <cassert>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>

void help();
void checkForce(const size_t numChecks, ...);
bool isUnsigned(const char *val);
bool isDouble(const char *val);
bool isOption(const char *val);
template <typename T>
void optionWarning(const char option, const char *name, const T val);
void checkOption(const int argc, char * const argv[], int &optionIndex, 
                 const char option, unsigned numOptions, ...);
void parseCommandLineOptions(int argc, char * const argv[]);
void checkParams(const char* inputFile);
void checkFitsError(const int error, const int lineNumber);
void deleteFitsFile(char * const filename, int &error);
void fitsFileExists(char * const filename, int &error);
void fitsFileCreate(fitsfile **file, char * const fileName, int &error);
void setParticleRows();

using namespace std;
using namespace chrono;

#define clear_line "\33[2K"                 //!< VT100 signal to clear line.
typedef duration<long, ratio<86400> > days; //!< Used to keep track of code execution time


enum clFlagType : int {						
	CI, //!< cloud_index
	D,  //!< double
	F,  //!< file_index
};

typedef int file_index;             //!< Used to keep track of file input arguments

bool Mach = false;                  //!< Perform Mach cone experiment
bool rk4 = true;					//!< Use Runge-Kutta-4 integration
bool pflag = false;                 //!< Parameter flag, true if parameter file used
double voidDecay = 0.4;             //!< Decay constant in ConfinementForceVoid [m^-1]
double magneticFieldStrength = 1.0; //!< Magnitude of B-field in z-direction [T]
double startTime = 0.0;             //!< Start time of simulation [s]
double dataTimeStep = 0.01;         //!< Time step for data output [s]
double simTimeStep =
    dataTimeStep/100.0;             //!< Time step for simulation [s]
double spacing = 0.003;             //!< Inter-particle spacing [m]
double endTime = 5.0;               //!< End time for simulation [s]
double confinementConst = 100;      //!< Strength of ConfinementForce [V/m^2]
double confinementConstX = 100;     //!< Strength of RectConfinementForce in x-direction [V/m^2]
double confinementConstY = 1000;    //!< Strength of RectConfinementForce in y-direction [V/m^2]
double shieldingConstant = 2E4;     //!< Defines inverse of distance where ShieldedCoulombForce kicks in [m^-1]
									//!< corresponds to 10*(ion debye length) [m^-1]
double shieldingConstX = 2E4;       //!< Inverse shielding length in x-direction [m^-1]
double shieldingConstY = 2E4;       //!< Inverse shielding lenght in y-direction [m^-1]
double dragGamma = 10.0;            //!< Dust drag frequency [Hz]
double thermRed = 1E-14;            //!< Thermal reduction factor [N]
double thermRed1 = thermRed;        //!< Outer reduction factor (-L) [N]
double thermScale = 1E-14;          //!< Default for TimeVaryingThermalForce [N/s]
double thermOffset = 0.0;           //!< Default for TimeVaryingThermalForce [N]
double heatRadius = 0.001;          //!< Apply ThermalForce only within this radius [m]
double driveConst = 0.00001;        //!< Used in DrivingForce for waves [m^2]
double waveAmplitude = 1E-13;       //!< Driving wave amplitude (default comparable to other forces throughout cloud) [N]
double waveShift = 0.007;           //!< Driving wave shift [m]
double machSpeed = 0.2;             //!< Firing speed for Mach Cone experiment [m/s]
double massFactor = 100;            //!< mass ultiplier for fired Mach Cone particle
double qMean = 6000.0;              //!< Mean number of charges of the gaussian charge distriburion [c]
double qSigma = 100.0;              //!< Standard deviation of number of charges [c]
double rMean = 1.45E-6;             //!< Mean dust particle radius of the gaussian size distribution [m]
double rSigma = 0.0;                //!< Standard deviation of the dust size distribution [m]
double rmin = 
	spacing* 
	5.0;  							//!< Inner radius of shear layer [m]
double rmax = 
	spacing*10.0; 					//!< Outer radius of shear layer [m]
double rotConst = 1E-15;            //!< Rotational force in shear layer [N]
double dragScale = -1.0;            //!< Scale factor in TimeVaryingDragForce [Hz/s]
double electricFieldStrength = 0;   //!< Stength of ElectricForce [V/m^2]
double plasmaRadius = 1;            //!< Decay constant for ElectricForce [m]
double vertElectricFieldStrength=0; //!< Strength of VertElectricForce [V/m^2]
double verticalDecay = 1;          	//!< Decay constant for VertElectricForce [m]
double gravitationalFieldStrength=0;//!< Strenght of GravitationalForce [m]
double justifyX = 0;                //!< Translation of cloud in x-direction [m]
double justifyY = 0;                //!< Translation of cloud in y-direction [m]
double massDensity = 2200;          //!< Density of cloud particles [kg/m^3]
double velocityX = 0.0;             //!< Initial x-velocity of cloud [m/s]
double velocityY = 0.0;             //!< Initial y-velocity of cloud [m/s]

force_flags usedForces = 0;         //!< Bitpacked forces
cloud_index numParticles = 4;		//!< Number of dust particles
cloud_index row_x_particles = 4;	//!< Number of rows in the x-direction
cloud_index row_y_particles = 0;	//!< Number of rows in the y-direction

file_index continueFileIndex = 0;   //!< Index of argv array that holds the file name of the fitsfile to continue. 
file_index finalsFileIndex = 0;     //!< Index of argv array that holds the file name of the fitsfile to use finals of.
file_index outputFileIndex = 0;     //!< Index of argv array that holds the file name of the fitsfile to output.
file_index inputFileIndex = 0;      //!< Input parameter file

double Cloud::interParticleSpacing = spacing;
double Cloud::dustParticleMassDensity = massDensity;
double Cloud::justX = justifyX;
double Cloud::justY = justifyY;
double Cloud::velX = velocityX;
double Cloud::velY = velocityY;


/**
* @brief Displays help to the console.
*
* @details Display help. The output is white space sensitive to render correctly
*          in an 80 column terminal environment. There should be no tabs.
**/
void help() {
     cout << endl 
          << "                                      DEMON" << endl
          << "        Dynamic Exploration of Microparticle clouds Optimized Numerically" << endl << endl
          << "Options:" << endl << endl
          << " -B 1.0                 set magnitude of B-field in z-direction [T]" << endl
          << " -c noDefault.fits      continue run from file" << endl
          << " -C 100.0               set confinementConst [V/m^2]" << endl
          << " -D -1.0 10.0           use TimeVaryingDragForce; set scale [Hz/s], offset [Hz]" << endl
          << " -d 2200                set dust density [kg/m^-3]" << endl 
          << " -E 100 10              set Electric field strength [V/m]; decay constant [m]" << endl
          << " -e 5.0                 set simulation end time [s]" << endl
          << " -f noDefault.fits      use final positions and velocities from file" << endl
	  	  << " -F 100 1               set verticalElectricField and top and bottom positions" << endl
          << " -g 10.0                set dragGamma (magnitute of drag constant) [Hz]" << endl
	      << " -G 0.0                 set Gravitational field strength [m/s^2]" << endl
          << " -h                     display Help (instead of running)" << endl
          << " -I                     use 2nd order Runge-Kutta integrator" << endl
          << " -k 0 0                 kick the particles in the x;y directions [m/s]" << endl
          << " -i 0.003               set initial inter-particle spacing [m]" << endl
          << " -L 0.001 1E-14 1E-14   use ThermalForceLocalized; set radius [m], in,out" << endl
          << "                        thermal values [N]" << endl
          << " -M 0.2 100             create Mach Cone; set bullet velocity [m/s], mass factor" << endl
          << " -n 8                   set number of particles" << endl
          << " -o 0.01                set the data Output time step [s]" << endl
          << " -O data.fits           set the name of the output file" << endl
          << " -P Parameters.cfg      Read parameters from file" << endl
          << " -p 0 0                 set initial x;y positions [m] of cloud" << endl
          << " -q 6000.0 100.0        set charge mean and sigma [c]" << endl
          << " -R 100.0 1000.0        use RectConfinementForce; set confineConstX,Y [V/m^2]" << endl
          << " -r 1.45E-6 0.0         set mean particle radius and sigma [m]" << endl
          << " -s 2E4                 set coulomb shielding constant [m^-1]" << endl
          << " -S 1E-15 0.005 0.007   use RotationalForce; set strength [N], rmin, rmax [m]" << endl
          << " -t 0.0001              set the simulation time step [s]" << endl
          << " -T 1E-14               use ThermalForce; set thermal reduction factor [N]" << endl
          << " -v 1E-14 0.0           use TimeVaryingThermalForce; set scale [N/s]" << endl
          << "                        and offset [N]" << endl
          << " -V 0.4                 use ConfinementForceVoid; set void decay constant [m^-1]" << endl
          << " -w 1E-13 0.007 0.00001 use DrivingForce; set amplitude [N], shift [m]," << endl
          << "                        driveConst [m^-2]" << endl << endl

          << "Notes: " << endl << endl
          << " Parameters specified above represent the default values and accepted type," << endl
          << "    with the exception of -c and -f, for which there are no default values." << endl
          << " -c appends to file; ignores all force flags (use -f to run with different" << endl
          << "    forces). -c overrides -f if both are specified" << endl
          << " -D uses strengthening drag if scale > 0, weakening drag if scale < 0." << endl
          << " -M is best used by loading up a previous cloud that has reached equilibrium." << endl
          << " -n expects even number, else will add 1 (required for SIMD)." << endl
          << " -S creates a shear layer between rmin = cloudsize/2 and" << endl
          << "    rmax = rmin + cloudsize/5." << endl
          << " -T runs with heat; otherwise, runs cold." << endl
          << " -v uses increases temp if scale > 0, decreasing temp if scale < 0." << endl
          << " -w creates acoustic waves along the x-axis (best with -R)." << endl 
          << " -E is set to 0 0 initially. If you would like to run DEMON with an" << endl
          << "    Electric force, you may also want to turn the Confinement Force to 0." << endl <<endl;
}


/**
* @brief Checks fits file for errors.
*
* @param[in] error      The error code
* @param[in] lineNumber The line number where the error occured
*
**/
void checkFitsError(const int error, const int lineNumber) {
	if (!error)
		return;

	char message[80];
	fits_read_errmsg(message);
	cout << "Error: Fits file error " << error 
	<< " at line number " << lineNumber 
	<< " (driver_2D.cpp)" << endl 
	<< message << endl;
	exit(1);
}

/**
* @brief Deletes existing fits file.
*
* @param[in]  filename The name of the fits file to delete
* @param[out] error    The error code (if any) that is produced when deleting the file
*
**/
void deleteFitsFile(char * const filename, int &error) {
	int exists = 0;
	fits_file_exists(filename, &exists, &error);

	if (exists) {
		cout << "Warning: Removing pre-existing \"" << filename << "\" file." << endl;
		remove(filename);
	}
	checkFitsError(error, __LINE__);
}

/**
* @brief Checks if fits file exists.
*
* @param[in]  filename The name of the fits file to delete
* @param[out] error    The error code (if any) that is produced when deleting the file
*
**/
void fitsFileExists(char * const filename, int &error) {
    int exists = 0;
    fits_file_exists(filename, &exists, &error);
    if (!exists) {
        cout << "Error: Fits file \"" << filename << "\" does not exist." << endl;
        help();
        exit(1);
    }
    checkFitsError(error, __LINE__);
}

 
/**
* @brief Creates a new fits file by deleting an existing one if it has the same name.
*
* @param[in]  file     The fits file to be created
* @param[in]  fileName The name of the fits file to delete
* @param[out] error    The error code (if any) that is produced when deleting the file
*
**/
void fitsFileCreate(fitsfile **file, char * const fileName, int &error) {
	deleteFitsFile(fileName, error);
	fits_create_file(file, fileName, &error);
	checkFitsError(error, __LINE__);
}


/**
* @brief Changes number/arrangement of particles if necesssary
*
* @details If only one number is defined, DEMON will try to arrange the grid
*          as squarely as possible. Otherwise, they will be set as rows and columns.
*		   There must be a multiple of either 4 or 8 particles, depending on AVX support.
**/
void setParticleRows() {
	if(row_y_particles == 0)
		numParticles = row_x_particles;
	else
		numParticles = row_x_particles*row_y_particles;

	if (numParticles < 4) {
		row_x_particles = 2;
		row_y_particles = 2;
		numParticles = row_x_particles * row_y_particles;
		cout << "Warning: -n requires multiples of " << FLOAT_STRIDE << " numbers of particles. Incrementing number of particles to (" 
		<< numParticles << ")." << endl;
	}
	else if(numParticles%FLOAT_STRIDE) {
		while(numParticles%FLOAT_STRIDE) {
			if(row_y_particles == 0) {
				++row_x_particles;
				numParticles = row_x_particles;
			}
			else if(row_y_particles < row_x_particles) {
				++row_y_particles;
				numParticles = row_x_particles * row_y_particles;
			}
			else {
				++row_x_particles;
				numParticles = row_x_particles * row_y_particles;
			}
		}

		cout << "Warning: -n requires multiples of " << FLOAT_STRIDE << " numbers of particles. Incrementing number of particles to (" 
		<< numParticles << ")." << endl;
	}
}


/**
* @brief Parses command line, prepares fits files, and begins simulation
*
* @param[in] argc Number of command line arguments
* @param[in] argv Array of command line arguments
**/
int main (int argc, char * const argv[]) {
	steady_clock::time_point start = steady_clock::now();

	parseCommandLineOptions(argc, argv);

    // All simulations require the folling three forces if subsitutes are not 
    // used.
	if (!(usedForces & TimeVaryingDragForceFlag))
		usedForces |= DragForceFlag;
	if (!(usedForces & RectConfinementForceFlag) && !(usedForces & ConfinementForceVoidFlag))
		usedForces |= ConfinementForceFlag;
	usedForces |= ShieldedCoulombForceFlag;

	fitsfile *file = NULL;
	int error = 0;
	Cloud *cloud;

	if (continueFileIndex) {
        // Create a cloud using a specified fits file. Subsequent time step data
        // will be appended to this fits file.
		fitsFileExists(argv[continueFileIndex], error);

		fits_open_file(&file, argv[continueFileIndex], READWRITE, &error);
		checkFitsError(error, __LINE__);

		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &usedForces, NULL, &error);
		checkFitsError(error, __LINE__);

		cloud = Cloud::initializeFromFile(file, error, &startTime);
		checkFitsError(error, __LINE__);
	} else if (finalsFileIndex) {
        // Create a cloud using the last time step of a specified fits file.
        // Subsequent time step data will be written to a new file.
		fitsFileExists(argv[finalsFileIndex], error);
		fits_open_file(&file, argv[finalsFileIndex], READONLY, &error);
		checkFitsError(error, __LINE__);
        
		cloud = Cloud::initializeFromFile(file, error, NULL);
		checkFitsError(error, __LINE__);
		
 		fits_close_file(file, &error);
		checkFitsError(error, __LINE__);
	} else
		cloud = Cloud::initializeGrid(numParticles, row_x_particles, row_y_particles, rMean, rSigma, qMean, qSigma);

	// Create a new file if we aren't continuing an old one.
	if (!continueFileIndex) {
		fitsFileCreate(&file, outputFileIndex ? argv[outputFileIndex] 
                               : const_cast<char *> ("data.fits"), error);
	
		// create "proper" primary HDU
		// (prevents fits from generating errors when creating binary tables)
		fits_create_img(file, 16, 0, NULL, &error);
		checkFitsError(error, __LINE__);
	}
	
    // Create all forces specified in used forces.
    ForceArray forces;
	if (usedForces & ConfinementForceFlag)
		forces.push_back(new ConfinementForce(cloud, confinementConst));
	if (usedForces & ConfinementForceVoidFlag)
		forces.push_back(new ConfinementForceVoid(cloud, confinementConst, voidDecay));
	if (usedForces & DragForceFlag) 
		forces.push_back(new DragForce(cloud, dragGamma));
	if (usedForces & DrivingForceFlag)
		forces.push_back(new DrivingForce(cloud, driveConst, waveAmplitude, waveShift));
	if (usedForces & MagneticForceFlag)
		forces.push_back(new MagneticForce(cloud, magneticFieldStrength));
	if (usedForces & RectConfinementForceFlag)
		forces.push_back(new RectConfinementForce(cloud, confinementConstX, confinementConstY));
	if (usedForces & RotationalForceFlag)
		forces.push_back(new RotationalForce(cloud, rmin, rmax, rotConst));
    /**
    if (usedForces & ShieldedCoulombForceFlag)
        forces.push_back(new ShieldedCoulombForce(cloud, shieldingConstant));
    **/
    if (usedForces & ShieldedCoulombForceFlag)
		forces.push_back(new ShieldedCoulombForce(cloud, shieldingConstX, shieldingConstY));
	if (usedForces & ThermalForceFlag)
		forces.push_back(new ThermalForce(cloud, thermRed));
	if (usedForces & ThermalForceLocalizedFlag)
		forces.push_back(new ThermalForceLocalized(cloud, thermRed, thermRed1, heatRadius));
	if (usedForces & TimeVaryingDragForceFlag)
		forces.push_back(new TimeVaryingDragForce(cloud, dragScale, dragGamma));
	if (usedForces & TimeVaryingThermalForceFlag)
		forces.push_back(new TimeVaryingThermalForce(cloud, thermScale, thermOffset));
	if (usedForces & ElectricForceFlag)
		forces.push_back(new ElectricForce(cloud, electricFieldStrength, plasmaRadius));
	if (usedForces & GravitationalForceFlag)
		forces.push_back(new GravitationalForce(cloud, gravitationalFieldStrength));
	if (usedForces & VertElectricForceFlag)
		forces.push_back(new VertElectricForce(cloud, vertElectricFieldStrength, verticalDecay));

	
	if (continueFileIndex) { // Initialize forces from old file.
		for (Force * const F : forces)
			F->readForce(file, &error);
		checkFitsError(error, __LINE__);
	} else { // Write force config to new file.
		for (Force * const F : forces)
			F->writeForce(file, &error);
		checkFitsError(error, __LINE__);
	}

	
    // Write initial data.
	if (!continueFileIndex) {
		cloud->writeCloudSetup(file, error);
		checkFitsError(error, __LINE__);
	} else {
		fits_movnam_hdu(file, BINARY_TBL, const_cast<char *> ("TIME_STEP"), 0, &error);
		checkFitsError(error, __LINE__);
	}
	
    // If performing the mach cone experiments alter the first particle. The 
    // first particle is moved to the left of the cloud. Its mass is increased
    // and it is given an inital velocity toward the main cloud.
	if (Mach) {
		cloud->x[0] = -0.75*sqrt((double)cloud->n)*spacing;
		cloud->y[0] = 0.0;
		cloud->Vx[0] = machSpeed;
		cloud->Vy[0] = 0.0;
		cloud->mass[0] *= massFactor;
	}
    
    // Create 2nd or 4th order Runge-Kutta integrator.
    Integrator * const I = rk4 ? new Runge_Kutta4(cloud, forces, simTimeStep, startTime)
                               : new Runge_Kutta2(cloud, forces, simTimeStep, startTime);

	// Run the simulation. Add a blank line to provide space between warnings
    // the completion counter.
    cout << endl;
	while (startTime < endTime) {
		cout << clear_line << "\rCurrent Time: " << I->currentTime << "s (" 
		<< I->currentTime/endTime*100.0 << "% Complete)" << flush;
		
		// Advance simulation to next timestep.
		I->moveParticles(startTime += dataTimeStep);
		cloud->writeTimeStep(file, error, I->currentTime);
	}

	// Close fits file.
	fits_close_file(file, &error);

	// clean up objects:
	for (Force * const F : forces)
		delete F;
	delete cloud;
    delete I;

	// Calculate and display elapsed time.
	const auto totalTime = steady_clock::now() - start;
	const days d = duration_cast<days> (totalTime);
	const hours h = duration_cast<hours> (totalTime - d);
	const minutes m = duration_cast<minutes> (totalTime - d - h);
	const seconds s = duration_cast<seconds> (totalTime - d - h - m);
	
	cout << clear_line << "\rTime elapsed: " 
	<< d.count() << (d.count() == 1 ? " day, " : " days, ") 
	<< h.count() << (h.count() == 1 ? " hour, " : " hours, ") 
	<< m.count() << (m.count() == 1 ? " minute, " : " minutes, ") 
	<< s.count() << (s.count() == 1 ? " second." : " seconds.") << endl;
	
	return 0;
}

/**
* @brief Checks if force is used or conflicts with a previously set force
*
* @param[in] numChecks Number of forces to check
* @param[in] ...       List of forces
**/
void checkForce(const size_t numChecks, ...) {
	va_list arglist;
	va_start(arglist, numChecks);
	
	const char firstOption = (char)va_arg(arglist, int);
	const ForceFlag firstFlag = (ForceFlag)va_arg(arglist, long);
	
	if (usedForces & firstFlag) {
		cout << "Error: option -" << firstOption << " already set." << endl;
		help();
		va_end(arglist);
		exit(1);
	}
	
	for (size_t i = 1; i < numChecks; i++) {
		const char nextOption = (char)va_arg(arglist, int);
		const ForceFlag nextFlag = (ForceFlag)va_arg(arglist, int);
		
		if (usedForces & nextFlag) {
			cout << "Error: option -" << firstOption << " conflicts with option -" << nextOption << endl;
			help();
			va_end(arglist);
			exit(1);
		}
	}
	
	va_end(arglist);
	usedForces |= firstFlag;
}

/**
* @brief Checks if string is a positive integer
*
* @param[in] val The string to check
*
* @return True if val is a positive integer, false if not
**/
bool isUnsigned(const char *val) {
	for (const char *c = val; *c != '\0'; c++)
		if (!isdigit(*c))
			return false;
	return true;
}

/**
* @brief Checks if string uses a decimal or scientific notation.
*
* @param[in] val The string to check
*
* @return True if val uses a decimal or scientific notation, false if it does not
**/
bool isDouble(const char *val) {
	for (const char *c = val; *c != '\0'; c++)
		if (!isdigit(*c) && *c != 'e' && *c != 'E' && *c != '.' && *c != '-')
			return false;
	return true;
}

/**
* @brief Checks if string string has the form "-x".
*
* @param[in] val The string to check
*
* @return True if val has the form "-x"
**/
bool isOption(const char *val) {
	return val[0] == '-' && isalpha(val[1]) && val[2] == '\0';
}

/**
* @brief Warns about incomplete command line parses.
*
* @param[in] option The character for the option that was incomplete
* @param[in] name   The complete name of the option that was incomplete
* @param[in] val 	The default value for the option
**/
template <typename T>
void optionWarning(const char option, const char *name, const T val) {
	cout << "Warning: -" << option << " option incomplete. Using default " 
	<< name << " (" << val << ")." << endl;
}

/**
* @brief Check command-line options and uses defaults if values are missing.
*
* @param[in]     argc        Number of command-line arguments
* @param[in]     argv        Array of command-line arguments
* @param[in,out] optionIndex Current index of option that is being checked
* @param[in] 	 option      Character of current option that is being checked
* @param[in] 	 numOptions  Number of options to check
*
**/
void checkOption(const int argc, char * const argv[], int &optionIndex, const char option, unsigned numOptions, ...) {
	++optionIndex;
	va_list arglist;
	va_start(arglist, numOptions);
	
	for (unsigned int i = 0; i < numOptions; i++) {
		const char *name = va_arg(arglist, char *);
		const clFlagType type = (clFlagType)va_arg(arglist, int);
		void *val = va_arg(arglist, void *);
		
		switch (type) {
			case CI: { // cloud_index argument
				cloud_index *ci = (cloud_index *)val;
				if (optionIndex < argc && isUnsigned(argv[optionIndex]))
					*ci = (cloud_index)atoi(argv[optionIndex++]);
				else
					optionWarning<cloud_index> (option, name, *ci);
				break;
			}
			case D: { // double argument 
				double *d = (double *)val;
				if (optionIndex < argc && !isOption(argv[optionIndex]) && isDouble(argv[optionIndex]))
					*d = atof(argv[optionIndex++]);
				else
					optionWarning<double> (option, name, *d);
				break;
			}
			case F: { // file_index argument
				const char *defaultFileName = va_arg(arglist, char *);
				file_index *fi = (file_index *)val;
				if (optionIndex < argc && !isOption(argv[optionIndex]) && !isDouble(argv[optionIndex]) && !isUnsigned(argv[optionIndex]))
					*fi = optionIndex++;
				else
					optionWarning<const char *> (option, name, defaultFileName);
				break;
			}
			default:
				va_end(arglist);
				assert(false && "Undefined Argument Type");
		}
	}
	va_end(arglist);
}


/**
* @brief Checks parameters from the parameter file.
*
* @param[in] inputFile The name of the parameter file to check
*
**/
void checkParams(const char* inputFile){
//    cout << "check " << inputFile << endl;
    ifstream paramfile;
    string varname, value;
    paramfile.open(inputFile);
    if (!inputFile){
        cout << "Incorrect parameter filepath. Please input an appropriate parameter file." << endl;
    }

    while (paramfile >> varname >> value){

        // Now a bunch of if statements to parse this.
        if (varname == "voidDecay"){
            voidDecay = atof(value.c_str());
        }
        if (varname == "magneticFieldStrength"){
            magneticFieldStrength = atof(value.c_str());
        }
        if (varname == "numParticles"){
            numParticles = atof(value.c_str());
        }
        if (varname == "startTime"){
            startTime = atof(value.c_str());
        }
        if (varname == "dataTimeStep"){
            dataTimeStep = atof(value.c_str());
        }
        if (varname == "simTimeStep"){
            simTimeStep = atof(value.c_str());
        }
        if (varname == "spacing"){
            Cloud::interParticleSpacing = atof(value.c_str());
        }
        if (varname == "endTime"){
            endTime = atof(value.c_str());
        }
        if (varname == "confinementConst"){
            confinementConst = atof(value.c_str());
        }
        if (varname == "confinementConstX"){
            confinementConstX = atof(value.c_str());
        }
        if (varname == "confinementConstY"){
            confinementConstY = atof(value.c_str());
        }
        if (varname == "shieldingConstant"){
            shieldingConstant = atof(value.c_str());
        }
        if (varname == "shieldingConstX"){
            shieldingConstX = atof(value.c_str());
        }
        if (varname == "shieldingConstY"){
            shieldingConstY = atof(value.c_str());
        }
        if (varname == "dragGamma"){
            dragGamma = atof(value.c_str());
        }
        if (varname == "thermRed"){
            thermRed = atof(value.c_str());
        }
        if (varname == "thermRed1"){
            thermRed1 = atof(value.c_str());
        }
        if (varname == "thermScale"){
            thermScale = atof(value.c_str());
        }
        if (varname == "thermOffset"){
            thermOffset = atof(value.c_str());
        }
        if (varname == "heatRadius"){
            heatRadius = atof(value.c_str());
        }
        if (varname == "driveConst"){
            driveConst = atof(value.c_str());
        }
        if (varname == "waveAmplitude"){
            waveAmplitude = atof(value.c_str());
        }
        if (varname == "waveShift"){
            waveShift = atof(value.c_str());
        }
        if (varname == "machSpeed"){
            machSpeed = atof(value.c_str());
        }
        if (varname == "massFactor"){
            massFactor = atof(value.c_str());
        }
        if (varname == "qMean"){
            qMean = atof(value.c_str());
        }
        if (varname == "qSigma"){
            qSigma = atof(value.c_str());
        }
        if (varname == "rMean"){
            rMean = atof(value.c_str());
        }
        if (varname == "rSigma"){
            rSigma = atof(value.c_str());
        }
        if (varname == "rmin"){
            rmin = atof(value.c_str());
        }
        if (varname == "rmax"){
            rmax = atof(value.c_str());
        }
        if (varname == "rotConst"){
            rotConst = atof(value.c_str());
        }
        if (varname == "dragScale"){
            dragScale = atof(value.c_str());
        }
        if (varname == "electricFieldStrength"){
            electricFieldStrength = atof(value.c_str());
        }
        if (varname == "plasmaRadius"){
            plasmaRadius = atof(value.c_str());
        }
        if (varname == "vertElectricFieldStrength"){
            vertElectricFieldStrength = atof(value.c_str());
        }
        if (varname == "verticalDecay"){
            verticalDecay = atof(value.c_str());
        }
        if (varname == "gravitationalFieldStrength"){
            gravitationalFieldStrength = atof(value.c_str());
        }
        if (varname == "justifyX"){
            Cloud::justX = atof(value.c_str());
        }
        if (varname == "justifyY"){
            Cloud::justY = atof(value.c_str());
        }
        if (varname == "massDensity"){
            Cloud::dustParticleMassDensity = atof(value.c_str());
        }
        if (varname == "velocityX"){
            Cloud::velX = atof(value.c_str());
        }
        if (varname == "velocityY"){
            Cloud::velY = atof(value.c_str());
        }
        if (varname == "forceFlags"){
            // Now we need to flip the appropriate force flags
            vector<string> flags;
            int comma;
            while (value.size() > 1){
                comma = value.find(",");
                flags.push_back(value.substr(0,value.find(",")));
                value = value.substr(comma+1,value.size() - comma);
            }

           // Now we need to figure out which forces to use
           for (int i = 0; i < flags.size(); i++){
               if (flags[i] == "B"){
                   checkForce(1, 'B', MagneticForceFlag);
               }
               if (flags[i] == "D"){
                   checkForce(1, 'D', TimeVaryingDragForceFlag);
               }
               if (flags[i] == "L"){
                   checkForce(3, 
                              'L', ThermalForceLocalizedFlag, 
                              'T', ThermalForceFlag, 
                              'v', TimeVaryingThermalForceFlag);

               }
               if (flags[i] == "R"){
                   checkForce(1, 'R', RectConfinementForceFlag);
               }
               if (flags[i] == "S"){
                   checkForce(1, 'S', RotationalForceFlag);
               }
               if (flags[i] == "T"){
                   checkForce(3, 
                              'T', ThermalForceFlag, 
                              'L', ThermalForceLocalizedFlag, 
                              'v', TimeVaryingThermalForceFlag);

               }
               if (flags[i] == "v"){
                   checkForce(3,
                              'v', TimeVaryingThermalForceFlag,
                              'L', ThermalForceLocalizedFlag,
                              'T', ThermalForceFlag);

               }
               if (flags[i] == "V"){
                   checkForce(1, 'V', ConfinementForceVoidFlag);
               }
               if (flags[i] == "w"){
                   checkForce(1, 'w', DrivingForceFlag);
               }
               if (flags[i] == "E"){
                   checkForce(1, 'E', ElectricForceFlag);
               }
               if (flags[i] == "F"){
                   checkForce(1, 'F', VertElectricForceFlag);
               }
               if (flags[i] == "G"){
                   checkForce(1, 'G', GravitationalForceFlag);
               }
           }
        }
    } 
}


/**
* @brief Parses the command-line options
*
* @param[in] argc Number of command-line arguments
* @param[in] argv Array of command-line arguments
*
**/
void parseCommandLineOptions(int argc, char * const argv[]){
    // argv[0] is the name of the executable. The routine checkOption increments
    // the array index internally. If check option is not used, it must be
    // incremented manually.
	for (int i = 1; i < argc;) {
		switch (argv[i][1]) {
	        // All F cases:
	        case 'P': // Read "P"arameter file
                checkOption(argc, argv, i, 'P', 1, 
            	"input file", F, &inputFileIndex, "Params.cfg");
                pflag = true;

                // We now need to parse the config file:
                // Follow flags below to determine how to read everything in
                checkParams(argv[inputFileIndex]);
                break;

	        case 'O': // name "O"utput file:
                checkOption(argc, argv, i, 'O', 1, 
            				"output file", F, &outputFileIndex, "data.fits");
                break;

	        case 'c': // "c"ontinue from file:
                checkOption(argc, argv, i, 'c', 1, 
           				    "continue file", F, &continueFileIndex, "");
                break;

	        case 'f': // use "f"inal positions and velocities from previous run:
                checkOption(argc, argv, i, 'f', 1, 
            				"finals file", F, &finalsFileIndex, "");
                break;

	        // All CI cases
	        case 'n': // set "n"umber of particles:
				checkOption(argc, argv, i, 'n', 2, 
	                        "number of row x particles", CI, &row_x_particles,
	                        "number of row y particles", CI, &row_y_particles);
				setParticleRows();
				break;

            // All D Cases. 
            // Note: if pflag = true, these are not going to be read.
            if (pflag == false) {
			
				case 'B': // set "B"-field:
					checkForce(1, 'B', MagneticForceFlag);
					checkOption(argc, argv, i, 'B', 1, 
	                            "magnetic field", D, &magneticFieldStrength);
					break;
				case 'C': // set "C"onfinementConst:
					checkOption(argc, argv, i, 'C', 1, 
	                            "confinementConst", D, &confinementConst);
					break;
				case 'D': // use TimeVarying"D"ragForce:
					checkForce(1, 'D', TimeVaryingDragForceFlag);
					checkOption(argc, argv, i, 'D', 2, 
	                            "scale factor", D, &dragScale, 
	                            "offset",       D, &dragGamma);
					break;
				case 'e': // set "e"nd time:
					checkOption(argc, argv, i, 'e', 1, 
	                            "end time", D, &endTime);
					break;		
				case 'g': // set "g"amma:
					checkOption(argc, argv, i, 'g', 1, 
	                            "dragGamma", D, &dragGamma);
					break;
				case 'h': // display "h"elp:
					help();
					exit(0);
                case 'I': // use 2nd order "i"ntegrator
                        rk4 = false;
                        i++;
                        break;
				case 'L': // perform "L"ocalized heating experiment:
					checkForce(3, 
	                           'L', ThermalForceLocalizedFlag, 
	                           'T', ThermalForceFlag, 
	                           'v', TimeVaryingThermalForceFlag);
					checkOption(argc, argv, i, 'L', 3, 
	                            "radius",       D, &heatRadius, 
	                            "heat factor1", D, &thermRed, 
	                            "heat factor2", D, &thermRed1);
					break;
				case 'M': // perform "M"ach Cone experiment:
					Mach = true;
					checkOption(argc, argv, i, 'M', 2, 
	                            "velocity", D, &machSpeed, 
	                            "mass",     D, &massFactor);
					break;
				case 'o': // set dataTimeStep, which conrols "o"utput rate:
					checkOption(argc, argv, i, 'o', 1, 
	                            "data time step", D, &dataTimeStep);
					break;
                case 'q':
                    checkOption(argc, argv, i, 'q', 2, 
                				"mean number of charges",  D, &qMean,
                				"number of charges sigma", D, &qSigma);
                    break;
                 case 'r': // set dust "r"adius
                    checkOption(argc, argv, i, 'r', 2, 
                				"mean dust radius", D, &rMean,
                				"dust radius sigma", D, &rSigma);
                    break;
				case 'R': // use "R"ectangular confinement:
					checkForce(1, 'R', RectConfinementForceFlag);
					checkOption(argc, argv, i, 'R', 2, 
	                            "confine constantX", D, &confinementConstX, 
	                            "confine constantY", D, &confinementConstY);
					break;
				case 's': // set "s"hielding constant:
					checkOption(argc, argv, i, 's', 1, 
	                            "shielding constant", D, &shieldingConstant);
					break;
				case 'S': // create rotational "S"hear layer:
					checkForce(1, 'S', RotationalForceFlag);
					checkOption(argc, argv, i, 'S', 3, 
	                            "force constant", D, &rotConst, 
	                            "rmin",           D, &rmin, 
	                            "rmax",           D, &rmax);
					break;
				case 't': // set "t"imestep:
					checkOption(argc, argv, i, 't', 1, "time step", D, &simTimeStep);
					break;
				case 'T': // set "T"emperature reduction factor:
					checkForce(3, 
	                           'T', ThermalForceFlag,
	                           'L', ThermalForceLocalizedFlag, 
	                           'v', TimeVaryingThermalForceFlag);
					checkOption(argc, argv, i, 'T', 1, 
	                            "heat factor", D, &thermRed);
					break;
				case 'v': // use time ""arying thermal force:
					checkForce(3,
	                           'v', TimeVaryingThermalForceFlag,
	                           'L', ThermalForceLocalizedFlag,
	                           'T', ThermalForceFlag);
					checkOption(argc, argv, i, 'v', 2,
	                            "heat value scale",  D, &thermScale,
	                            "heat value offset", D, &thermOffset);
					break;
				case 'V': // use ConfinementForceVoid:
					checkForce(1, 'V', ConfinementForceVoidFlag);
					checkOption(argc, argv, i, 'V', 1, 
	                            "void decay", D, &voidDecay);
					break;
				case 'w': // drive "w"aves:
					checkForce(1, 'w', DrivingForceFlag);
					checkOption(argc, argv, i, 'w', 3, 
	                            "amplitude",        D, &waveAmplitude, 
	                            "wave shift",       D, &waveShift, 
	                            "driving constant", D, &driveConst);
					break;
        		case 'E': // use "E"lectricForce:
        			checkForce(1, 'E', ElectricForceFlag);
        			checkOption(argc, argv, i, 'E', 2, 
                            	"electric field const", D, &electricFieldStrength, 
                            	"plasma radius", D, &plasmaRadius);
        			break;
				case 'F': // use VertElectricForce:
        			checkForce(1, 'F', VertElectricForceFlag);
        			checkOption(argc, argv, i, 'F', 2, 
                            "vert electric field const",  D, &vertElectricFieldStrength, 
                            "vertical decay const",       D, &verticalDecay);
					break;
				case 'G': // use "G"ravitationalForce
					checkForce(1, 'G', GravitationalForceFlag);
					checkOption(argc, argv, i, 'G', 1,
				    			"gravitational field const", D, &gravitationalFieldStrength);
					break;
				case 'd': // set dust mass "d"ensity:
			        checkOption(argc, argv, i, 'd', 1, 
                           		"massDensity", D, &Cloud::dustParticleMassDensity);
                    break;
				case 'i': // set "i"nter-particle spacing:
			        checkOption(argc, argv, i, 'i', 1, 
                            	"spacing", D, &Cloud::interParticleSpacing);
                    break;
                case 'p': // set "p"osition [x,y]:
                    checkOption(argc, argv, i, 'p', 2, 
                				"justify x", D, &Cloud::justX, 
                				"justify y", D, &Cloud::justY);
                    break;
                case 'k': // velocity "kick" [x,y]:
                    checkOption(argc, argv, i, 'k', 2, 
                    			"velocity x", D, &Cloud::velX, 
                    			"velocity y", D, &Cloud::velY);
        			break;
	        }

			default: // Handle unknown options by issuing error.
				cout << "Error: Unknown option " << argv[i] << endl;
				help();
				exit(1);
		}
	}
}
