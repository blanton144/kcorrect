#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

/* applies evolution to an absolute magnitude */

float k_evolve(float absm, 
               float z, 
               float q0, 
               float q1, 
               float qz0)
{
  return(absm+q0*(1.+q1*(z-qz0))*(z-qz0));
}
