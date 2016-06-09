/**
* @file  FFTAnalysis.cpp
* @brief Peforms FFT analysis on DEMON FITS files
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#include <iostream>
#include <cstdarg>
#include <cassert>
#include <cmath>
#include <fstream>
#include <string>
#include <fftw3.h>
#include "fitsio.h"


using namespace std;

#define clear_line "\33[2K"

void help();
bool isUnsigned(const char *val);
bool isDouble(const char *val);
bool isOption(const char *val);
template <typename T>
void optionWarning(const char option, const char *name, const T val);
void checkOption(const int argc, char * const argv[], int &optionIndex, 
                 const char option, unsigned numOptions, ...);
void parseCommandLineOptions(int argc, char * const argv[]);
void checkFitsError(const int error, const int lineNumber);
void averageFFT(fitsfile *file, int hdutype, int numPars, int argc, char * const argv[]);

enum clFlagType : int {
   I, // integer
   D,  // double
   F,  // file_index
};

typedef int file_index;

int numParticles = 1;

file_index continueFileIndex = 0;   //!< Index of argv array that holds the file name of the fitsfile to continue. 
file_index finalsFileIndex = 0;     //!< Index of argv array that holds the file name of the fitsfile to use finals of.
file_index outputFileIndex = 0;     //!< Index of argv array that holds the file name of the fitsfile to output.
file_index inputFileIndex = 0;      //!< Input parameter file

bool performAverageFFT = false;



int main(int argc, char *argv[]) {

   parseCommandLineOptions(argc, argv);

   int error   = 0,
       hdutype = 0;

   fitsfile *file = NULL;

   // Variables to be extracted from the fits file
   long int npars      = 0,
            nrows      = 0,
            ntimesteps = 0;

   // Open input file
   fits_open_file(&file, argv[inputFileIndex], READONLY, &error);
   checkFitsError(error, __LINE__);

   //Get number of particles
   fits_movabs_hdu(file, 2, &hdutype, &error);
   fits_get_num_rows(file, &npars, &error);
   if(npars < numParticles)
      numParticles = npars;

   //If -a argument used, do averageFFT analysis
   if(performAverageFFT)
      averageFFT(file, hdutype, numParticles, argc, argv);


   fits_close_file(file, &error);
}


//---------------------------------------------------------------------------------------------------------------------//

/**
* @brief Displays help to the console.
*
* @details Display help. The output is white space sensitive to render correctly
*          in an 80 column terminal environment. There should be no tabs.
**/
void help() {
     cout << endl 
          << "                                      FFTAnalysis" << endl << endl
          << "Options:" << endl << endl
          << " -a 1                   Compute average FFT of particle trajectories" << endl
          << "                        with specified number of particles" << endl
          << " -h                     display Help (instead of running)" << endl                          
          << " -i data.fits           Set input file" << endl
          << " -o output.dat          Set output file" << endl<< endl;
}

void parseCommandLineOptions(int argc, char * const argv[]) {
    // argv[0] is the name of the executable. The routine checkOption increments
    // the array index internally. If check option is not used, it must be
    // incremented manually.
   for (int i = 1; i < argc;) {
      switch (argv[i][1]) {
         case 'i': 
            checkOption(argc, argv, i, 'i', 1, "input file", F, &inputFileIndex, "data.fits");
            break;

         case 'o': // name "O"utput file:
            checkOption(argc, argv, i, 'o', 1, 
                        "output file", F, &outputFileIndex, "output.dat");
            break;
         case 'a':
            performAverageFFT = true;
            checkOption(argc, argv, i, 'a', 1,
                        "average number", I, &numParticles);
            break;
         case 'h': // display "h"elp:
               help();
               exit(0);
      }
   }
}

void checkOption(const int argc, char * const argv[], int &optionIndex, const char option, unsigned numOptions, ...) {
   ++optionIndex;
   va_list arglist;
   va_start(arglist, numOptions);
   
   for (unsigned int i = 0; i < numOptions; i++) {
      const char *name = va_arg(arglist, char *);
      const clFlagType type = (clFlagType)va_arg(arglist, int);
      void *val = va_arg(arglist, void *);
      
      switch (type) {
         case I: { // cloud_index argument
            int *i = (int *)val;
            if (optionIndex < argc && isUnsigned(argv[optionIndex]))
               *i = (int)atoi(argv[optionIndex++]);
            else
               optionWarning<int> (option, name, *i);
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

bool isUnsigned(const char *val) {
   for (const char *c = val; *c != '\0'; c++)
      if (!isdigit(*c))
         return false;
   return true;
}

bool isDouble(const char *val) {
   for (const char *c = val; *c != '\0'; c++)
      if (!isdigit(*c) && *c != 'e' && *c != 'E' && *c != '.' && *c != '-')
         return false;
   return true;
}

bool isOption(const char *val) {
   return val[0] == '-' && isalpha(val[1]) && val[2] == '\0';
}

template <typename T>
void optionWarning(const char option, const char *name, const T val) {
   cout << "Warning: -" << option << " option incomplete. Using default " 
   << name << " (" << val << ")." << endl;
}

void checkFitsError(const int error, const int lineNumber) {
   if (!error)
      return;

   char message[80];
   fits_read_errmsg(message);
   cout << "Error: Fits file error " << error 
   << " at line number " << lineNumber 
   << " FFTAnalysis.cpp" << endl 
   << message << endl;
   exit(1);
}

void averageFFT(fitsfile *file, int hdutype, int numPars, int argc, char * const argv[]) {
   int error = 0;
   long int ntimesteps;

   // ALlocate memory for FFT
   fftw_complex *in, *out;

   //Allocate memory to store fits file data
   double *cloudx  = new double[numPars];
   double amp = 0.0;

   //Move to the 3rd HDU to get the simulation/cloud information
   fits_movabs_hdu(file, 3, &hdutype, &error);

   //Get the number of timesteps (nrows)
   fits_get_num_rows(file, &ntimesteps, &error);

   // Allocate size of FFT data structures
   in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntimesteps);
   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntimesteps);

   // Create plan to speed up computation
   fftw_plan p;
   p = fftw_plan_dft_1d(ntimesteps, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

   double averageFFT[ntimesteps][2];

   //Note that nt is not indexed starting from zero
   for (long int npart = 0; npart < numPars; npart++) {
      for(long int nt = 1; nt <= ntimesteps; nt++) {
         //Read the file information into our "cloud" arrays
         //fits_read_col_dbl(ffptr,COLUMN,ROW,FIRSTELEMENT,NUMELEMENTS,NULLVALUE,ARRAY,&anynull,&error);
         fits_read_col_dbl(file, 2, nt, 1, numPars, 0.0,  cloudx, NULL, &error);

         in[nt-1][0] = cloudx[npart];
         in[nt-1][1] = 0.0;
      }

      fftw_execute(p);

      for(long int nt = 0; nt < ntimesteps; nt++) {
         averageFFT[nt][0] += out[nt][0];
         averageFFT[nt][1] += out[nt][1];
      }
   }

   //Free the pointer data
   delete[] cloudx;

   fftw_destroy_plan(p);
   fftw_free(in);
   fftw_free(out);

   ofstream outFile (outputFileIndex ? argv[outputFileIndex] 
                               : const_cast<char *> ("output.dat"));

   for(int i = 0; i < ntimesteps/2; i++) {
      
      amp = sqrt(averageFFT[i][0]*averageFFT[i][0] + averageFFT[i][1]*averageFFT[i][1]);

      // Bad fix for bug where random amplitudes are enormous
      if (amp / (double)numPars > 100) {
         averageFFT[i][0] = 0;
         averageFFT[i][1] = 0;
      }
      
      outFile << ((double) i)/((double)ntimesteps * 0.01) << "\t" << amp / ((double) numPars) << "\n";
   }

   outFile.close();
}