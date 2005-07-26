#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * ztransform.c: utilities to convert among redshift, lookback time,
 * comoving distance, volume element, enclosed volume, distance
 * modulus, and angular diameter distance.
 *
 * Works up to redshift 10 or so.
 *
 * Originally started circa 1999, Michael Blanton
 *
 */

#define ZRF_NUM 10000
#define ZRF_AMAX 1.
#define ZRF_AMIN 0.1
#define ZRF_C 2.99792e+5 /* in km/s */
#define LN10 2.3025851
#define TOL 1.e-6

static int zrf_initialized=0;
static float zrf_omega0;
static float zrf_omegal0;
static float zrf_z[ZRF_NUM];
static float zrf_t[ZRF_NUM];
static float zrf_r[ZRF_NUM];
static float zrf_dV[ZRF_NUM];
static float zrf_V[ZRF_NUM];
static float zrf_dm[ZRF_NUM];
static float zrf_add[ZRF_NUM];

float k_qromo(float (*func)(float), float a, float b,
	float (*choose)(float(*)(float), float, float, int));
float k_midpnt(float (*func)(float), float a, float b, int n);
void k_locate(float xx[], unsigned long n, float x, unsigned long *j);
void init_zrf(float omega0, float omegal0);
float ztor_flat_integrand(float a);
float z2t_flat_integrand(float a);
float ztor_open_integrand(float a);
float z2t_open_integrand(float a);

/* Converts redshift z into comoving distance r in km/s */
float ztor(float z, 
					 float omega0,
					 float omegal0)
{
	unsigned long i,ip1;
	float sz,r;

	if(!zrf_initialized || omega0!=zrf_omega0 || omegal0!=zrf_omegal0) 
		init_zrf(omega0, omegal0);

	k_locate(zrf_z, ZRF_NUM, z, &i);
	if(i==ZRF_NUM || i==ZRF_NUM-1) i=ZRF_NUM-2;
	if(i==-1) i=0;
	ip1=i+1;
	sz=z-zrf_z[i];
	r=zrf_r[i]+sz*(zrf_r[ip1]-zrf_r[i])/(zrf_z[ip1]-zrf_z[i]);
	return(r);
} /* end ztor */

/* Converts redshift z to comoving volume element dV in h^{-3} Mpc^3 */
float ztodV(float z, 
						 float omega0,
						 float omegal0)
{
	unsigned long i,ip1;
	float sz,dV;

	if(!zrf_initialized || omega0!=zrf_omega0 || omegal0!=zrf_omegal0) 
		init_zrf(omega0, omegal0);

	k_locate(zrf_z, ZRF_NUM, z, &i);
	if(i==ZRF_NUM || i==ZRF_NUM-1) i=ZRF_NUM-2;
	if(i==-1) i=0;
	ip1=i+1;
	sz=z-zrf_z[i];
	dV=zrf_dV[i]+sz*(zrf_dV[ip1]-zrf_dV[i])/(zrf_z[ip1]-zrf_z[i]);
	return(dV);
} /* end ztodV */

/* Converts redshift z to comoving enclosed volume in h^{-3} Mpc^3 */
float ztoV(float z, 
						float omega0,
						float omegal0)
{
	unsigned long i,ip1;
	float sz,V;

	if(!zrf_initialized || omega0!=zrf_omega0 || omegal0!=zrf_omegal0) 
		init_zrf(omega0, omegal0);

	k_locate(zrf_z, ZRF_NUM, z, &i);
	if(i==ZRF_NUM || i==ZRF_NUM-1) i=ZRF_NUM-2;
	if(i==-1) i=0;
	ip1=i+1;
	sz=z-zrf_z[i];
	V=zrf_V[i]+sz*(zrf_V[ip1]-zrf_V[i])/(zrf_z[ip1]-zrf_z[i]);
	return(V);
} /* end ztodV */

/* Converts to redshift z from comoving enclosed volume in h^{-3} Mpc^3 */
float Vtoz(float V, 
						float omega0,
						float omegal0)
{
	unsigned long i,ip1;
	float sV,z;

	if(!zrf_initialized || omega0!=zrf_omega0 || omegal0!=zrf_omegal0) 
		init_zrf(omega0, omegal0);

	k_locate(zrf_V, ZRF_NUM, V, &i);
	if(i==ZRF_NUM || i==ZRF_NUM-1) i=ZRF_NUM-2;
	if(i==-1) i=0;
	ip1=i+1;
	sV=V-zrf_V[i];
	z=zrf_z[i]+sV*(zrf_z[ip1]-zrf_z[i])/(zrf_V[ip1]-zrf_V[i]);
	return(z);
} /* end ztodV */

/* Converts comoving distance r in km/s to redshift z */
float rtoz(float r, 
					 float omega0,
					 float omegal0)
{
	unsigned long i,ip1;
	float sr,z;

	if(!zrf_initialized || omega0!=zrf_omega0 || omegal0!=zrf_omegal0) 
		init_zrf(omega0, omegal0);

	k_locate(zrf_r, ZRF_NUM, r, &i);
	if(i==ZRF_NUM || i==ZRF_NUM-1) i=ZRF_NUM-2;
	if(i==-1) i=0;
	ip1=i+1;
	sr=r-zrf_r[i];
	z=zrf_z[i]+sr*(zrf_z[ip1]-zrf_z[i])/(zrf_r[ip1]-zrf_r[i]);
	return(z);
} /* end rotz */

/* Converts redshift z into distance modulus dm */
float z2dm(float z, 
					 float omega0,
					 float omegal0)
{
	unsigned long i, ip1;
	float sz,dm;

	if(!zrf_initialized || omega0!=zrf_omega0 || omegal0!=zrf_omegal0) 
		init_zrf(omega0, omegal0);

	k_locate(zrf_z, ZRF_NUM, (float)z, &i);
	if(i==ZRF_NUM || i==ZRF_NUM-1) i=ZRF_NUM-2;
	if(i==-1) i=0;
	ip1=i+1;
	sz=z-zrf_z[i];
	dm=(zrf_dm[i]+sz*(zrf_dm[ip1]-zrf_dm[i])/(zrf_z[ip1]-zrf_z[i]));
	return(dm);
} /* end z2dm */
								 
/* Converts redshift z into distance modulus dm */
float dm2z(float dm, 
					 float omega0,
					 float omegal0)
{
	unsigned long i, ip1;
	float sdm,z;

	if(!zrf_initialized || omega0!=zrf_omega0 || omegal0!=zrf_omegal0) 
		init_zrf(omega0, omegal0);

	k_locate(zrf_dm, ZRF_NUM, (float)dm, &i);
	if(i==ZRF_NUM || i==ZRF_NUM-1) i=ZRF_NUM-2;
	if(i==-1) i=0;
	ip1=i+1;
	sdm=dm-zrf_dm[i];
	z=(zrf_z[i]+sdm*(zrf_z[ip1]-zrf_z[i])/(zrf_dm[ip1]-zrf_dm[i]));
	return(z);
} /* end z2dm */
								 
								 
/* Gives the angular diameter distance in km/s at redshift z */
float z2add(float z, 
						float omega0,
						float omegal0)
{
	unsigned long i, ip1;
	float sz,add;

	if(!zrf_initialized || omega0!=zrf_omega0 || omegal0!=zrf_omegal0) 
		init_zrf(omega0, omegal0);

	k_locate(zrf_z, ZRF_NUM, (float)z, &i);
	if(i==ZRF_NUM || i==ZRF_NUM-1) i=ZRF_NUM-2;
	if(i==-1) i=0;
	ip1=i+1;
	sz=z-zrf_z[i];
	add=(zrf_add[i]+sz*(zrf_add[ip1]-zrf_add[i])/(zrf_z[ip1]-zrf_z[i]));
	return(add);
} /* end z2add */

/* returns age of universe at redshift z in h^{-1} Gyrs */
float z2t(float z, 
					float omega0,
					float omegal0)
{
	unsigned long i, ip1;
	float sz,t;

	if(!zrf_initialized || omega0!=zrf_omega0 || omegal0!=zrf_omegal0) 
		init_zrf(omega0, omegal0);

	k_locate(zrf_z, ZRF_NUM, (float)z, &i);
	if(i==ZRF_NUM || i==ZRF_NUM-1) i=ZRF_NUM-2;
	if(i==-1) i=0;
	ip1=i+1;
	sz=z-zrf_z[i];
	t=(zrf_t[i]+sz*(zrf_t[ip1]-zrf_t[i])/(zrf_z[ip1]-zrf_z[i]));
	return(t);
} /* end z2t */

/* returns age of universe at redshift z in h^{-1} Gyrs */
float t2z(float t, 
					float omega0,
					float omegal0)
{
	unsigned long i, ip1;
	float st,z;

	if(!zrf_initialized || omega0!=zrf_omega0 || omegal0!=zrf_omegal0) 
		init_zrf(omega0, omegal0);

	k_locate(zrf_t, ZRF_NUM, (float)t, &i);
	if(i==ZRF_NUM || i==ZRF_NUM-1) i=ZRF_NUM-2;
	if(i==-1) i=0;
	ip1=i+1;
	st=t-zrf_t[i];
	z=(zrf_z[i]+st*(zrf_z[ip1]-zrf_z[i])/(zrf_t[ip1]-zrf_t[i]));
	return(z);
} /* end t2z */

float ztor_open_integrand(float a)
{
	float value;
	
	value=1./sqrt(zrf_omega0*a+(1.-zrf_omega0)*a*a);
	return(value);
} /* end ztor_integrand */

float z2t_open_integrand(float a)
{
	float value;
	
	value=9.78/sqrt(zrf_omega0/a+(1.-zrf_omega0));
	return(value);
} /* end ztor_integrand */

float ztor_flat_integrand(float a)
{
	float a2,value;
	
	a2=a*a;
	value=1./sqrt(zrf_omega0*a+zrf_omegal0*a2*a2);
	return(value);
} /* end ztor_integrand */

float z2t_flat_integrand(float a)
{
	float value;
	
	value=9.78/sqrt(zrf_omega0/a+zrf_omegal0*a*a);
	return(value);
} /* end ztor_integrand */

void init_zrf(float omega0, 
							float omegal0)
{
	float amin,DM,Eint;
	unsigned long i,flag; /* flag=0 for open, 1 for flat */

	if(fabs(omega0+omegal0-1.)<TOL)
		flag=1;
	else if (fabs(omegal0)<TOL && omega0<=1. && omega0>=0.)
		flag=0;
	else {
		printf("omega0=%le and omegal0=%le not allowed in ztransform.c\n",
					 omega0,omegal0);
		exit(1);
	} /* end if..else */

	zrf_initialized=1;
	zrf_omega0=omega0;
	zrf_omegal0=omegal0;
	
	if(flag==0) {
		for(i=0;i<ZRF_NUM;i++) {
			amin=ZRF_AMIN+(ZRF_AMAX-ZRF_AMIN)*((float)i+0.5)/(float)ZRF_NUM;
			zrf_z[i]=1./amin-1.;
			Eint=k_qromo(ztor_open_integrand,amin,1.,k_midpnt);
			zrf_r[i]=ZRF_C*Eint;
			DM=ZRF_C*sinh(sqrt(1.-zrf_omega0)*Eint)/sqrt(1.-zrf_omega0);
			zrf_dm[i]=25.+5.*log10(0.01*DM*(1.+zrf_z[i]));
			zrf_add[i]=DM/(1.+zrf_z[i]);
			zrf_dV[i]=1.e-6*ZRF_C*DM*DM/
				sqrt(omega0*(1.+zrf_z[i])*(1.+zrf_z[i])*(1.+zrf_z[i])
						 +(1.-omega0)*(1.+zrf_z[i])*(1.+zrf_z[i]));
			zrf_V[i]=(1.5*1.e-6*ZRF_C*ZRF_C*ZRF_C/(1.-omega0))*
				(DM/ZRF_C*sqrt(1.+(1.-omega0)*DM*DM/(ZRF_C*ZRF_C))
				 -asinh(sqrt(1.-omega0)*DM/ZRF_C)/sqrt(1.-omega0));
			zrf_t[i]=k_qromo(z2t_open_integrand,0.,amin,k_midpnt);
		} /* end for */
	} else {
		for(i=0;i<ZRF_NUM;i++) {
			amin=ZRF_AMIN+(ZRF_AMAX-ZRF_AMIN)*((float)i+0.5)/(float)ZRF_NUM;
			zrf_z[i]=1./amin-1.;
			zrf_r[i]=ZRF_C*k_qromo(ztor_flat_integrand,amin,1.,k_midpnt);
			zrf_dm[i]=25.+5.*log10(0.01*zrf_r[i]*(1.+zrf_z[i]));
			zrf_dV[i]=1.e-6*ZRF_C*zrf_r[i]*zrf_r[i]/
				sqrt(omega0*(1.+zrf_z[i])*(1.+zrf_z[i])*(1.+zrf_z[i])
						 +omegal0);
			zrf_V[i]=1.e-6*zrf_r[i]*zrf_r[i]*zrf_r[i];
			zrf_add[i]=zrf_r[i]/(1.+zrf_z[i]);
			zrf_t[i]=k_qromo(z2t_flat_integrand,0.,amin,k_midpnt);
		} /* end for */
	} /* end if..else */
} /* end init_zrf */

