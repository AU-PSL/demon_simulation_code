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
 {
   freopen ("com.dat","w",stdout); //opens output file for array writing

    fitsfile *fptr; //establishes FITS pointer
    char *val, value[1000], nullstr[]="*";
    char keyword[FLEN_KEYWORD], colname[FLEN_VALUE];
    int status = 0;
    int hdunum, hdutype, ncols, ii, anynul, dispwidth[1000];
    int firstcol, lastcol = 0, linewidth;
    long jj, nrows;

    if (argc != 2) {
      printf("NO INPUT");
      return(0);
    }

    if (!fits_open_file(&fptr, argv[1], READONLY, &status)) //opens FITS file for data
    {
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

          printf("\n");

          val = value; 
          for (jj = 1; jj <= nrows && !status; jj++) {
              printf("%4d ", jj);
              for (ii = firstcol; ii <= lastcol; ii++)
              {
                  if (fits_read_col_str (fptr,ii,jj, 1, 1, nullstr,
                      &val, &anynul, &status) )
                     break;

                  printf("%-*s ",dispwidth[ii], value);
              }
              printf("\n");
          }
      }
      fits_close_file(fptr, &status);
    } 

    if (status) fits_report_error(stderr, status);
    return(status);

   fclose (stdout);
   return 0;
 }

}
