/* ----------------------------------------------------------------------------
   calibrate.c
-------------------------------------------------------------------------------

   Take a filter trace (and possibly atmospheric absorption trace?) and the
   Hayes spectrophotometry of Vega and compute the absolute calibration for the
   bandpass.

   Usage: calibrate trace.txt

-------------------------------------------------------------------------------
   David W. Hogg / hogg@ias.edu / December 1998
---------------------------------------------------------------------------- */

#include <stdio.h>

#define MAXLAM 1000 /* maximum number of wavelength values in traces */

main(argc,argv)
     int argc ;
     char *argv[] ;
{
  /* Declare functions and variables */
  int fclose(),nlam,i ;
  double lamHayes[MAXLAM],lamtrace[MAXLAM],lamflam[MAXLAM],trace[MAXLAM],
    index,lameff ;
  FILE *fp,*fopen() ;

  /* Open and read files */

  /* Compute effective wavelength(s) */

  /* Compute intensity from Vega */
  for(i=0;i<nlam;i+)

  /* Loop over values for the spectral index */

  /* Compute intensity for power law */

  /* Compute and output calibration */
}
