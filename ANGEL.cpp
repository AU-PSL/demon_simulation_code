/*===- ANGEL.cpp -=============================================================/
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include <string.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include "fitsio.h"
using namespace std;

int main(int argc, char *argv[])
{
    std::fstream angel;
    angel.open ("com.dat", std::fstream::out);

    fitsfile *fptr; //establishes FITS pointer
    char *val, value[1000], nullstr[]="*";
    char keyword[FLEN_KEYWORD], colname[FLEN_VALUE];
    int hdunum, hdutype, ncols, anynul, dispwidth[1000], status = 0;
    int firstcol, lastcol = 0, linewidth;
    long nrows;

    if (argc != 2) {
      angel << "NO INPUT";
      return(0);
    }

    fits_open_file(&fptr, argv[1], READONLY, &status);
    if (!status) {
      if ( fits_get_hdu_num(fptr, &hdunum) == 1 )
          fits_movabs_hdu(fptr, 2, &hdutype, &status);
       else 
          fits_get_hdu_type(fptr, &hdutype, &status);

        fits_get_num_rows(fptr, &nrows, &status);
        fits_get_num_cols(fptr, &ncols, &status);

         while(lastcol < ncols) {
           linewidth = 0;
           firstcol = lastcol+1;
           for (lastcol = firstcol; lastcol <= ncols; lastcol++) {
              fits_get_col_display_width
                 (fptr, lastcol, &dispwidth[lastcol], &status);
              linewidth += dispwidth[lastcol] + 1;
           }

          if (lastcol > firstcol)lastcol--; 

          val = value; 
          for ( long jj = 1; jj <= nrows && !status; jj++) {
              angel << jj;
              for (int ii = firstcol; ii <= lastcol; ii++)
              {
                  if (fits_read_col_str (fptr,ii,jj, 1, 1, nullstr,
                      &val, &anynul, &status) )
                     break;

                  angel << value;
              }
              angel << endl;
          }
      }
      fits_close_file(fptr, &status);
    } 

      angel.close();

      if (status) fits_report_error(stderr, status);
      return(status);
      return 0;
}
