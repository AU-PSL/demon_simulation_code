/*===- driver_2D.cpp - Driver -=================================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

//Force classes:
#include "ConfinementForce.h"
#include "DragForce.h"
#include "DrivingForce.h"
#include "RectConfinementForce.h"
#include "RotationalForce.h"
#include "ShieldedCoulombForce.h"
#include "ThermalForce.h"
#include "ThermalForceLocalized.h"
#include "TimeVaryingDragForce.h"
#include "TimeVaryingThermalForce.h"

//Runge-Kutta class:
#include "Runge_Kutta.h"

//Other dependencies:
#include <ctime>
#include <iostream>
using namespace std;

#define clear_line "\33[2K" // VT100 signal to clear line.

void help()
{
 //80 cols is ********************************************************************************
     cout << endl 
          << "                                      DEMON" << endl
          << "        Dynamic Exploration of Microparticle clouds Optimized Numerically" << endl << endl
          << "Options:" << endl << endl
          << " -c noDefault.fits      continue run from file" << endl
          << " -C 1E-13               set confinementConst" << endl
          << " -d -1.0 10.0           use TimeVaryingDragForce; set scale, offset" << endl
          << " -D 2                   set number of spatial dimensions" << endl
          << " -e 5.0                 set simulation end time" << endl
          << " -f noDefaut.fits       use final positions and velocities from file" << endl
          << " -g 10.0                set gamma (magnitute of drag constant)" << endl
          << " -h                     display Help (instead of running)" << endl
          << " -L 0.001 1E-14 1E-14   use ThermalForceLocalized; set rad, in/out therm vals" << endl
          << " -M 0.2 100             create Mach Cone; set bullet velocity, mass factor" << endl
          << " -n 10                  set number of particles" << endl
          << " -o 0.01                set the data Output time step" << endl
          << " -O data.fits           set the name of the output file" << endl
          << " -r 0.01                set cloud radius (one-half side length)" << endl
          << " -R 1E-13 1E-12         use RectConfinementForce; set confineConstX,Y" << endl
          << " -s 2E4                 set coulomb shelding constant" << endl
          << " -S 1E-15 0.005 0.007   use RotationalForce; set strength, rmin, rmax" << endl
          << " -t 0.0001              set the simulation time step" << endl
          << " -T 1E-14               use ThermalForce; set thermal reduction factor" << endl
          << " -v 1E-14 0.0           use TimeVaryingThermalForce; set scale and offset" << endl
          << " -w 1E-13 0.007 0.00001 use DrivingForce; set amplitude, shift, driveConst" << endl << endl
          << "Notes: " << endl << endl
          << " Parameters specified above represent the default values and accepted type," << endl
          << "    with the exception of -c and -f, for which there are no default values." << endl
          << " -c appends to file; ignores all force flags (use -f to run with different" << endl
          << "    forces). -c overrides -f if both are specified." << endl
          << " -d uses strengthening drag if scale > 0, weakening drag if scale < 0." << endl
          << " -D may have argument 1, 2, or 3 for 1D, 2D, or 3D clouds, respectively." << endl
          << " -M is best used by loading up a previous cloud that has reached equalibrium." << endl
          << " -n expects even number, else will add 1 (required for SIMD)." << endl
          << " -S creates a shear layer between rmin = cloudsize/2 and" << endl
          << "    rmax = rmin + cloudsize/5." << endl
          << " -T runs with heat; otherwise, runs cold." << endl
          << " -v uses increases temp if scale > 0, decreasing temp if scale < 0." << endl
          << " -w creates acoustic waves along the x-axis (best with -R)." << endl << endl;
}

//check if force is used:
void checkForce(const char option, const long usedForces, const ForceFlag flag)
{
	if(usedForces & flag)
	{
		cout << "Error: -" << option << " already set." << endl;
		help();
		exit(1);
	}
}

//check if using two incompatible forces:
void checkForce(char option1, char option2, long usedForces, ForceFlag flag1, ForceFlag flag2)
{
	checkForce(option1, usedForces, flag1);
	if(usedForces & flag2)
	{
		cout << "Error: -" << option1 << " cannot be used with -" << option2 << endl;
		help();
		exit(1);
	}
}

//check whether character is alphabetical:
inline const bool isCharacter(const char c)
{
	return (c > 'a' && c < 'z') || (c > 'A' && c < 'Z');
}

//check file name exists:
int checkFileOption(const int argc, char * const argv[], int i, const char option,
	const string name, int * const file)
{
	if(i+1 >= argc || argv[i+1][0] == '-')
	{
		cout << "Warning: -" << option << " option incomplete." << endl
			<< name << " missing." << endl;
		help();
		exit(i);
	}
	else
		*file = ++i;

	return i;
}

//check for one command line flag, use default value if absent:
int checkOption(const int argc, char * const argv[], int i, const char option, 
	const string name, double * const value)
{
	if(i+1 >= argc || argv[i+1][0] == '-')
		cout << "Warning: -" << option << " option incomplete." << endl 
			<< "Using default " << name << " (" << *value << ")." << endl;
	else
		*value = atof(argv[++i]);

	return i;
}

//check for one command line flag, use default value if absent (overloaded for int values):
int checkOption(const int argc, char * const argv[], int i, const char option, 
	const string name, int * const value)
{
	if(i+1 >= argc || argv[i+1][0] == '-')
		cout << "Warning: -" << option << " option incomplete." << endl 
			<< "Using default " << name << " (" << *value << ")." << endl;
	else
		*value = atoi(argv[++i]);

	return i;
}

//check for one command line flag, use default value if absent (overloaded for unsigned int values):
int checkOption(const int argc, char * const argv[], int i, const char option, 
	const string name, unsigned int * const value)
{
	if(i+1 >= argc || argv[i+1][0] == '-')
		cout << "Warning: -" << option << " option incomplete." << endl 
			<< "Using default " << name << " (" << *value << ")." << endl;
	else
		*value = atoi(argv[++i]);

	return i;
}

//check for two command line flags, use default values if absent:
int checkOption(const int argc, char * const argv[], int i, const char option, 
	const string name1, double * const value1, 
	const string name2, double * const value2)
{
	if(i+1 >= argc || argv[i+1][0] == '-')
		cout << "Warning: -" << option << " option incomplete." << endl
			<< "Using default "<< name1 << " (" << *value1 << ") and " 
			<< name2 << " (" << *value2 << ")." << endl;
	else
	{
		*value1 = atof(argv[++i]);
		i = checkOption(argc, argv, i, option, name2, value2);
	}

	return i;
}

//check for three command line flags, use default values if absent:
int checkOption(const int argc, char * const argv[], int i, const char option, 
	const string name1, double * const value1, 
	const string name2, double * const value2, 
	const string name3, double * const value3)
{
	if(i+1 >= argc || argv[i+1][0] == '-')
		cout << "Warning: -" << option << " option incomplete." << endl 
			<< "Using default "<< name1 << " (" << *value1 << "), " 
			<< name2 << " (" << *value2 << ") and " 
			<< name3 << " (" << *value3 << ")." << endl;
	else
	{
		*value1 = atof(argv[++i]);
		i = checkOption(argc, argv, i, option, name2, value2, name3, value3);
	}

	return i;
}

//check for two command line flags in the case of a negative argument,
//	use default values if absent:
int checkOptionWithNeg(const int argc, char * const argv[], int i, const char option, 
	const string name1, double * const value1, 
	const string name2, double * const value2)
{
	if(i+1 >= argc || (argv[i+1][0] == '-' && isCharacter(argv[i+1][1])))
		cout << "Warning: -" << option << " option incomplete." << endl 
			<< "Using default " << name1 << " (" << *value1 << ") and " 
			<< name2 << " (" << *value2 << ")." << endl << endl;
	else
	{
		*value1 = atof(argv[++i]);
		i = checkOption(argc, argv, i, option, name2, value2);
	}
	return i;
}

//count number of forces in use:
const unsigned int getNumForces(const long usedForces)
{
	unsigned int i = 0;
	if (usedForces & ConfinementForceFlag)
		++i;
	if (usedForces & DragForceFlag)
		++i;
	if (usedForces & ShieldedCoulombForceFlag)
		++i;
	if (usedForces & RectConfinementForceFlag)
		++i;
	if (usedForces & ThermalForceFlag)
		++i;
	if (usedForces & ThermalForceLocalizedFlag)
		++i;
	if (usedForces & DrivingForceFlag)
		++i;
	if (usedForces & RotationalForceFlag)
		++i;
	if (usedForces & TimeVaryingDragForceFlag)
		++i;
	if (usedForces & TimeVaryingThermalForceFlag)
		++i;
	return i;
}

//check fitsfile for errors:
void checkFitsError(const int error, const int lineNumber)
{
	if(!error)
		return;

	char message[80];
	fits_read_errmsg(message);
	cout << "Error: Fits file error " << error 
		<< " at line number " << lineNumber 
		<< " (driver_2D.cpp)" << endl 
		<< message << endl;
	exit(1);
}

//delete fitsfile:
void deleteFitsFile(char * const filename, int * const error)
{
	//check for pre-existing data file:
	int exists = 0;
	fits_file_exists(filename, &exists, error);

	if(exists)
	{
		cout << "Warning: Removing pre-existing \"" << filename << "\" file." << endl;
		remove(filename); //required by fits, else can't create
	}

	checkFitsError(*error, __LINE__);
}

//check if fits file exists:
void fitsFileExists(char * const filename, int * const error)
{
	int exists = 0;
	fits_file_exists(filename, &exists, error);
	if(exists != 1)
	{
		cout << "Error: Fits file \"" << filename << "\" does not exist." << endl;
		help();
		exit(1);
	}

	checkFitsError(*error, __LINE__);
	cout << "Initializing with fits file \"" << filename << "\"." << endl;
}

int main (int argc, char * const argv[]) 
{
	long timer = time(NULL); //start timer

	//object declarations:
	Cloud *cloud;
	Force **forceArray;                  //new pointer to Force object (will set to array)
	
	//declare variables and set default values:
	bool Mach = false;                   //true -> perform Mach Cone experiment
	double startTime = 0.0;
	double simTimeStep = 0.0001;
	double dataTimeStep = 0.01;
	double endTime = 5;
	double cloudSize = 0.01;             //one-half side length (aka "radius")
	double confinementConst = 1E-13;     //confinementForce
	double confinementConstX = 1E-13;    //RectConfinementForce
	double confinementConstY = 1E-12;    //RectConfinementForce
	double confinementConstZ = 1E-12;    //RectConfinementForce
	double shieldingConstant = 2E4;      //corresponds to 10*(ion debye length)
	double gamma = 10.0;
	double thermRed = 1E-14;             //default thermal reduction factor
	double thermRed1 = thermRed;         //default outer reduction factor (-L)
	double thermScale = 1E-14;           //default for TimeVaryingThermalForce
	double thermOffset = 0.0;            //default for TimeVaryingThermalForce
	double heatRadius = .001;            //apply thermal force only within this radius
	double driveConst = .00001;          //used in DrivingForce.cpp for waves
	double waveAmplitude = 1E-13;        //driving wave amplitude (default comparable to other forces throughout cloud)
	double waveShift = 0.007;            //driving wave shift
	double machSpeed = 0.2;              //firing speed for Mach Cone experiment
	double massFactor = 100;             //mass multiplier for fired Mach Cone particle
	double rmin = cloudSize/2.0;         //inner radius of shear layer
	double rmax = rmin + cloudSize/5.0;  //outer radius of shear layer
	double rotConst = 1E-15;             //rotational force in shear layer
	double dragScale = -1.0;             //used in TimeVaryingDragForce
	int continueFileIndex = 0;           //Index of argv array that holds the file name of the fitsfile to continue. 
	int finalsFileIndex = 0;             //Index of argv array that holds the file name of the fitsfile to use finals of.
	int outputFileIndex = 0;             //Index of argv array that holds the file name of the fitsfile to output.
	int dimension = 2;                   //1D, 2D, or 3D cloud
	long usedForces = 0;                 //bitpacked forces
	unsigned int numParticles = 10;
	unsigned int numForces = 3;

	//process command line flags:
	for(int i = 0; i < argc; i++)
	{
		switch(argv[i][1])
		{
			case 'c': //"c"ontinue from file:
				i = checkFileOption(argc, argv, i, 'c', "Continue file", &continueFileIndex);
				break;
			case 'C': //set "C"onfinementConst:
				i = checkOption(argc, argv, i, 'C', "confinementConst", &confinementConst);
				break;
			case 'd': //use TimeVarying"D"ragForce:
				checkForce('d', usedForces, TimeVaryingDragForceFlag);
				usedForces |= TimeVaryingDragForceFlag;
				//dragScale needs to allow negative numbers:
				i = checkOptionWithNeg(argc, argv, i, 'd', "scale factor", &dragScale, "offset", &gamma);
				break;
			case 'D': //set cloud "D"imension
				i = checkOption(argc, argv, i, 'D', "dimension", &dimension);
				if(dimension != 1 && dimension != 2 && dimension != 3)
				{
					cout << "Error: Invalid spatial dimension.\n";
					help();
					exit(1);
				}
				break;
			case 'e': //set "e"nd time:
				i = checkOption(argc, argv, i, 'e', "end time", &endTime);
				break;
			case 'f': //use "f"inal positions and velocities from previous run:
				i = checkOption(argc, argv, i, 'f', "Finals file", &finalsFileIndex);
				break;
			case 'g': //set "g"amma:
				i = checkOption(argc, argv, i, 'g', "gamma", &gamma);
				break;
			case 'h': //display "h"elp:
				help();
				exit(1);
				break;
			case 'L': //perform "L"ocalized heating experiment:
				checkForce('L', 'T', usedForces, ThermalForceLocalizedFlag, ThermalForceFlag);
				checkForce('L', 'v', usedForces, ThermalForceLocalizedFlag, TimeVaryingThermalForceFlag);
				usedForces |= ThermalForceLocalizedFlag;
				i = checkOption(argc, argv, i, 'L', "radius", &heatRadius, "heat factor1", &thermRed, "heat factor2", &thermRed1);
				break;
			case 'M': //perform "M"ach Cone experiment:
				Mach = true;
				i = checkOption(argc, argv, i, 'M', "velocity", &machSpeed, "mass", &massFactor);
				break;
			case 'n': //set "n"umber of particles:
				i = checkOption(argc, argv, i, 'n', "number of particles", &numParticles);
				if((numParticles % 2) != 0)     //odd
				{
					cout << "Even number of particles required for SIMD." << endl 
						<< "Incrementing number of particles to " << ++numParticles << endl;
				}
				break;
			case 'o': //set dataTimeStep, which conrols "o"utput rate:
				i = checkOption(argc, argv, i, 'o', "data time step", &dataTimeStep);
				break;
			case 'O': //name "O"utput file:
				i = checkOption(argc, argv, i, 'O', "output file", &outputFileIndex);
				break;
			case 'r': //set cloud "r"adius:
				i = checkOption(argc, argv, i, 'r', "cloud size", &cloudSize);
				break;
			case 'R': //use "R"ectangular confinement:
				checkForce('R', usedForces, RectConfinementForceFlag);
				usedForces |= RectConfinementForceFlag;
				i = checkOption(argc, argv, i, 'R', "confine constantX", &confinementConstX, "confine constantY", &confinementConstY, "confine constantY", &confinementConstZ);
				break;
			case 's': //set "s"hielding constant:
				i = checkOption(argc, argv, i, 's', "shielding constant", &shieldingConstant);
				break;
			case 'S': //create rotational "S"hear layer:
				checkForce('S', usedForces, RotationalForceFlag);
				usedForces |= RotationalForceFlag;
				i = checkOption(argc, argv, i, 'S', "force constant", &rotConst, "rmin", &rmin, "rmax", &rmax);
				break;
			case 't': //set "t"imestep:
				i = checkOption(argc, argv, i, 't', "time step", &simTimeStep);
				if(simTimeStep == 0.0) //prevent divide-by-zero error
				{
					cout << "Error: simTimeStep set to 0.0 with -t." << endl 
						<< "Terminating to prevent divide-by-zero." << endl;
					help();
					exit(1);
				}
				break;
			case 'T': //set "T"emperature reduction factor:
				checkForce('T', 'L', usedForces, ThermalForceFlag, ThermalForceLocalizedFlag);
				checkForce('T', 'v', usedForces, ThermalForceFlag, TimeVaryingThermalForceFlag);
				usedForces |= ThermalForceFlag;
				i = checkOption(argc, argv, i, 'T', "heat factor", &thermRed);
				break;
			case 'v': //use time ""arying thermal force:
				checkForce('v', 'T', usedForces, TimeVaryingThermalForceFlag, ThermalForceFlag);
				checkForce('v', 'L', usedForces, TimeVaryingThermalForceFlag, ThermalForceLocalizedFlag);
				usedForces |= TimeVaryingThermalForceFlag;
				i = checkOptionWithNeg(argc, argv, i, 'v', "heat value scale", &thermScale, "heat value offset", &thermOffset);
				break;
			case 'w': //drive "w"aves:
				checkForce('w', usedForces, DrivingForceFlag);
				usedForces |= DrivingForceFlag;
				i = checkOption(argc, argv, i, 'w', "amplitude", &waveAmplitude, "wave shift", &waveShift, "driving constant", &driveConst);
				break;
			default: //default (expressionless)
				break;
		}
	}

	if (!(usedForces & TimeVaryingDragForceFlag))
		usedForces |= DragForceFlag;
	if (!(usedForces & RectConfinementForceFlag))
		usedForces |= ConfinementForceFlag;
	usedForces |= ShieldedCoulombForceFlag;

/*------------------------------------------------------------------------------
 * Initialize cloud:
 -----------------------------------------------------------------------------*/
	cout << "Status: Initializing cloud." << endl;

	//declare fits file and error:
	fitsfile *file;
	int error = 0;

	if(continueFileIndex)
	{
		fitsFileExists(argv[continueFileIndex], &error);
		
		//open file:
		fits_open_file(&file, argv[continueFileIndex], READWRITE, &error); //file pointer, file name (char), read/write, error
		checkFitsError(error, __LINE__);
		
		//use the same forces:
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &usedForces, NULL, &error);
		checkFitsError(error, __LINE__);
		
		//initialize with last time step from file:
		cloud = Cloud::initializeFromFile(file, &error, &startTime);
		checkFitsError(error, __LINE__);
	}
	else if(finalsFileIndex)
	{
		fitsFileExists(argv[finalsFileIndex], &error);

		//open file:
		fits_open_file(&file, argv[finalsFileIndex], READONLY, &error); //file pointer, file name (char), read only, error
		checkFitsError(error, __LINE__);

		//initialize with last time step from file:
		cloud = Cloud::initializeFromFile(file, &error, NULL);
		checkFitsError(error, __LINE__);
		
		//close file:
		fits_close_file(file, &error);
		checkFitsError(error, __LINE__);
	}
	else //initialize new cloud on grid:
	{
		if(dimension == 1)
			cloud = Cloud::initializeLine(numParticles);
		else if(dimension == 2)
			cloud = Cloud::initializeSquare(numParticles);
		else if(dimension == 3)
			cloud = Cloud::initializeCube(numParticles);
		else //this should never happen
		{
			cout << "Error: Impossible case reached when determining dimension "
				"prior to calling Cloud::initializeGrid.\n";
			exit(1);
		}
		
	}

	// Create a new file if we aren't continueing one.
	if (!continueFileIndex)
	{
		if(outputFileIndex) //use specified file name
		{
			deleteFitsFile(argv[outputFileIndex], &error);
			fits_create_file(&file, argv[outputFileIndex], &error);
			checkFitsError(error, __LINE__);
			
			//create "proper" primary HDU
			//	(prevents fits from generating errors when creating binary tables)
			fits_create_img(file, 16, 0, NULL, &error);
			checkFitsError(error, __LINE__);
		}
		else    //use default file name
		{
			deleteFitsFile(const_cast<char *> ("data.fits"), &error);
			fits_create_file(&file, const_cast<char *> ("data.fits"), &error);
			checkFitsError(error, __LINE__);
			
			//create "proper" primary HDU
			//	(prevents fits from generating errors when creating binary tables)
			fits_create_img(file, 16, 0, NULL, &error);
			checkFitsError(error, __LINE__);
		}
	}	
/*------------------------------------------------------------------------------
 * This concludes initialization of cloud.
 -----------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------
 * Initialize array of Force objects:
 -----------------------------------------------------------------------------*/
	cout << "Status: Initializing forces." << endl;

	numForces = getNumForces(usedForces);
	forceArray = new Force*[numForces];
	unsigned int index = 0;
	if(dimension == 1)
	{
		if (usedForces & ConfinementForceFlag)
			forceArray[index++] = new ConfinementForce1D(cloud, confinementConst);
		if (usedForces & DragForceFlag) 
			forceArray[index++] = new DragForce1D(cloud, gamma);
		if (usedForces & ShieldedCoulombForceFlag) 
			forceArray[index++] = new ShieldedCoulombForce1D(cloud, shieldingConstant);
		if (usedForces & RectConfinementForceFlag)
			forceArray[index++] = new RectConfinementForce1D(cloud, confinementConstX);
		if (usedForces & ThermalForceFlag)
			forceArray[index++] = new ThermalForce1D(cloud, thermRed);
		if (usedForces & ThermalForceLocalizedFlag)
			forceArray[index++] = new ThermalForceLocalized1D(cloud, thermRed, thermRed1, heatRadius);
		if (usedForces & DrivingForceFlag)
			forceArray[index++] = new DrivingForce1D(cloud, driveConst, waveAmplitude, waveShift);
		if (usedForces & RotationalForceFlag)
		{
			cout << "Error: RotationalForce not defined for 1D\n";
			exit(1);
		}
		if (usedForces & TimeVaryingDragForceFlag)
			forceArray[index++] = new TimeVaryingDragForce1D(cloud, dragScale, gamma);
		if (usedForces & TimeVaryingThermalForceFlag)
			forceArray[index++] = new TimeVaryingThermalForce1D(cloud, thermScale, thermOffset);
	}
	else if(dimension == 2)
	{
		if (usedForces & ConfinementForceFlag)
			forceArray[index++] = new ConfinementForce2D(cloud, confinementConst);
		if (usedForces & DragForceFlag) 
			forceArray[index++] = new DragForce2D(cloud, gamma);
		if (usedForces & ShieldedCoulombForceFlag) 
			forceArray[index++] = new ShieldedCoulombForce2D(cloud, shieldingConstant);
		if (usedForces & RectConfinementForceFlag)
			forceArray[index++] = new RectConfinementForce2D(cloud, confinementConstX, confinementConstY);
		if (usedForces & ThermalForceFlag)
			forceArray[index++] = new ThermalForce2D(cloud, thermRed);
		if (usedForces & ThermalForceLocalizedFlag)
			forceArray[index++] = new ThermalForceLocalized2D(cloud, thermRed, thermRed1, heatRadius);
		if (usedForces & DrivingForceFlag)
			forceArray[index++] = new DrivingForce2D(cloud, driveConst, waveAmplitude, waveShift);
		if (usedForces & RotationalForceFlag)
			forceArray[index++] = new RotationalForce2D(cloud, rmin, rmax, rotConst);
		if (usedForces & TimeVaryingDragForceFlag)
			forceArray[index++] = new TimeVaryingDragForce2D(cloud, dragScale, gamma);
		if (usedForces & TimeVaryingThermalForceFlag)
			forceArray[index++] = new TimeVaryingThermalForce2D(cloud, thermScale, thermOffset);
	}
	if(dimension == 3)
	{
		if (usedForces & ConfinementForceFlag)
			forceArray[index++] = new ConfinementForce3D(cloud, confinementConst);
		if (usedForces & DragForceFlag) 
			forceArray[index++] = new DragForce3D(cloud, gamma);
		if (usedForces & ShieldedCoulombForceFlag) 
			forceArray[index++] = new ShieldedCoulombForce3D(cloud, shieldingConstant);
		if (usedForces & RectConfinementForceFlag)
			forceArray[index++] = new RectConfinementForce3D(cloud, confinementConstX, confinementConstY, confinementConstZ);
		if (usedForces & ThermalForceFlag)
			forceArray[index++] = new ThermalForce3D(cloud, thermRed);
		if (usedForces & ThermalForceLocalizedFlag)
			forceArray[index++] = new ThermalForceLocalized3D(cloud, thermRed, thermRed1, heatRadius);
		if (usedForces & DrivingForceFlag)
			forceArray[index++] = new DrivingForce3D(cloud, driveConst, waveAmplitude, waveShift);
		if (usedForces & RotationalForceFlag)
			forceArray[index++] = new RotationalForce3D(cloud, rmin, rmax, rotConst);
		if (usedForces & TimeVaryingDragForceFlag)
			forceArray[index++] = new TimeVaryingDragForce3D(cloud, dragScale, gamma);
		if (usedForces & TimeVaryingThermalForceFlag)
			forceArray[index++] = new TimeVaryingThermalForce3D(cloud, thermScale, thermOffset);
	}
	
	if(continueFileIndex) //initialize forces from old file
	{
		for (unsigned int i = 0; i < numForces; i++)
			forceArray[i]->readForce(file, &error);
		checkFitsError(error, __LINE__);
	}
	else                  //write forces to new file
	{
		for (unsigned int i = 0; i < numForces; i++)
			forceArray[i]->writeForce(file, &error);
		checkFitsError(error, __LINE__);
	}
/*------------------------------------------------------------------------------
 * This concludes initialization of Force objects.
 -----------------------------------------------------------------------------*/

	//write initial data:
	if (!continueFileIndex) 
	{
		cloud->writeCloudSetup(file, &error);
		checkFitsError(error, __LINE__);
	}
	else
	{
		fits_movnam_hdu(file, BINARY_TBL, const_cast<char *> ("TIME_STEP"), 0, &error);
		checkFitsError(error, __LINE__);
	}

	if(Mach) //TODO: make Mach 3D
	{
		if(dimension != 2)
		{
			cout << "Error: Mach currently only supported when dimension = 2.\n";
			exit(1);
		}

		//reserve particle 1 for mach experiment
		cloud->x[0] = -2.0*cloudSize;
		cloud->y[0] = 0.0;
		cloud->Vx[0] = machSpeed;
		cloud->Vy[0] = 0.0;
		cloud->mass[0] = massFactor*cloud->mass[0];
	}

/*------------------------------------------------------------------------------
 * Commence Runge-Kutta algorithm:
 -----------------------------------------------------------------------------*/
	cout << "Status: Commencing Runge-Kutta." << endl << endl;

	cout << "Dimension = " << dimension << endl;

	Runge_Kutta rk4(cloud, forceArray, simTimeStep, numForces, startTime);

	//execute simulation for desired length of time:
	if(dimension == 1) //condition farther up than necessary to avoid querying on each time step
	{
		while(startTime < endTime)
		{
			cout << clear_line << "\rCurrent Time: " << rk4.currentTime << "s (" 
				<< rk4.currentTime/endTime*100.0 << "% Complete)" << flush;

			//call Runge-Kutta algorithm:
			rk4.moveParticles(startTime += dataTimeStep);
			//write positions and velocities:
			cloud->writeTimeStep(file, &error, rk4.currentTime);
		}
	}
	if(dimension == 2)
	{
		while(startTime < endTime)
		{
			cout << clear_line << "\rCurrent Time: " << rk4.currentTime << "s (" 
				<< rk4.currentTime/endTime*100.0 << "% Complete)" << flush;

			//call Runge-Kutta algorithm:
			rk4.moveParticles(startTime += dataTimeStep);
			//write positions and velocities:
			cloud->writeTimeStep(file, &error, rk4.currentTime);
		}
	}
	if(dimension == 3)
	{
		while(startTime < endTime)
		{
			cout << clear_line << "\rCurrent Time: " << rk4.currentTime << "s (" 
				<< rk4.currentTime/endTime*100.0 << "% Complete)" << flush;

			//call Runge-Kutta algorithm:
			rk4.moveParticles(startTime += dataTimeStep);
			//write positions and velocities:
			cloud->writeTimeStep(file, &error, rk4.currentTime);
		}
	}
/*------------------------------------------------------------------------------
 * This concludes the Runge-Kutta algorithm.
 -----------------------------------------------------------------------------*/

	//close fits file:
	fits_close_file(file, &error);

	//calculate and display elapsed time:
	long seconds = time(NULL) - timer;
	long minutes = seconds/60;
	long hours = minutes/60;
	long days = hours/24;
	hours -= days*24;
	minutes -= hours*60 + days*1440;
	seconds -= minutes*60 + hours*3600 + days*86400;
	cout << clear_line << "\rTime elapsed: " 
		<< days << (days == 1 ? " day, " : " days, ") 
		<< hours << (hours == 1 ? " hour " : " hours, ") 
		<< minutes << (minutes == 1 ? " minute " : " minutes, ") 
		<< seconds << (seconds == 1 ? " second " : " seconds.") << endl;

	//clean up objects:
	for (unsigned int i = 0; i < numForces; i++)
		delete forceArray[i];

	delete[] forceArray;
//	delete cloud;

	return 0;
}
