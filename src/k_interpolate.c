#include <stdio.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_interpolate.c
 *
 * Just a linear interpolation of a vector
 *
 * Mike Blanton
 * 2/2001
 *
 */

/* linear interpolate in a given vector; out of bounds --> 0 */
double k_interpolate(double currpos,
										 double vals[],
										 double pos[],
										 IDL_LONG n)
{
	unsigned long i,ip1;
	double sp,currval;

	k_locate(pos, (unsigned long) n, currpos, &i);
	if(i>=n-1 || i<0) return(0.);
	ip1=i+1;
	sp=currpos-pos[i];
	currval=vals[i]+sp*(vals[ip1]-vals[i])/(pos[ip1]-pos[i]);
	return(currval);
} /* end interpolate templates */

/* linear interpolate in a given vector; out of bounds --> 0;
 * assume x-axis is given equally spaced */
double k_interpolate_es(double currpos,
											double vals[],
											double pos[],
											IDL_LONG n)
{
	unsigned long i,ip1;
	double sp,currval;

	i=(unsigned long) floor((double) n*(currpos-pos[0])/(pos[n-1]-pos[0]));
	if(i>=n-1 || i<0) return(0.);
	ip1=i+1;
	sp=currpos-pos[i];
	currval=vals[i]+sp*(vals[ip1]-vals[i])/(pos[ip1]-pos[i]);
	return(currval);
} /* end interpolate templates */
