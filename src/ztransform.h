/* Converts redshift z into comoving distance r in km/s */
float ztor(float z, float omega0, float omegal0);
/* Converts redshift z to comoving volume element dV in h^{-3} Mpc^3 */
float ztodV(float z, float omega0, float omegal0);
/* Converts redshift z to comoving enclosed volume in h^{-3} Mpc^3 */
float ztoV(float z, float omega0, float omegal0);
/* Converts to redshift z from comoving enclosed volume in h^{-3} Mpc^3 */
float Vtoz(float V, float omega0, float omegal0);
/* Converts comoving distance r in km/s to redshift z */
float rtoz(float r, float omega0, float omegal0);
/* Converts redshift z into distance modulus dm */
float z2dm(float z, float omega0, float omegal0);
/* Converts redshift z into distance modulus dm */
float dm2z(float dm, float omega0, float omegal0);
/* Gives the angular diameter distance in km/s at redshift z */
float z2add(float z, float omega0, float omegal0);
/* returns age of universe at redshift z in h^{-1} Gyrs */
float z2t(float z, float omega0, float omegal0);
