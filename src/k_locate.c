/* 
 * I have altered this routine to be compatible with zero-offset
 * arrays
 */
void k_locate(float xx[], unsigned long n, float x, unsigned long *j)
{
	unsigned long ju,jm,jl;
	int ascnd;

	jl=-1;
	ju=n;
	ascnd=(xx[n-1] > xx[0]);
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x > xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	*j=jl;
}
