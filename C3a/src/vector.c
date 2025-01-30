#include <stdio.h>
#include <math.h>
#include "vector.h"

double 
l2_norm(
	double *vector, 
	int length
		) 
{
    double sum = 0.0;
    for (int i = 0; i < length; i++) {
        sum += vector[i] * vector[i];
    }
    return sqrt(sum);
}


void
elementwise_addition(
		     double *res,
		     double *v1,
		     double *v2,
		     unsigned int len
	            )
{
    for(int i = 0; i < len; ++i){
	res[i] = v1[i] + v2[i];
    }
}

void
constant_multiplication(
			   double *res,
			   double *v1,
			   double a,
			   unsigned int len
	                  )
{
    for(int i = 0; i < len; ++i){
	res[i] = v1[i] * a;
    }
}

double
dot_product(
	    double *v1,
	    double *v2,
	    unsigned int len
	   )
{
    double result = 0;
    for(int i = 0; i < len; ++i){
	result += v1[i] * v2[i];
    }
    return result;
}
