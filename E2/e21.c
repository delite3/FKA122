#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "attempt.h"

#define PI 3.141592653589

/* **********************************************
 *
 * Perform Monte Carlo integration of the given
 * integral without using importance sampling.
 *
 * Parameters
 * ----------
 *  N - Number of points to sample
 *  k - GSL random number generator object
 *
 * Returns
 * -------
 *  Struct with the result
 *
 * **********************************************/
result_t MC_without_importance_sampling(int N, gsl_rng *k)
{
    result_t result;
    result.integral = 0.0; // Integral result
    result.error = 0.0;    // Error estimate

    double sum = 0.0;
    double sum_of_squares = 0.0;

    // Perform the sampling
    for (int i = 0; i < N; i++) {
        double x = gsl_rng_uniform(k);  // Generate random point in [0, 1]
        double f_x = x * (1 - x);       // Evaluate the function f(x) = x * (1 - x)

        sum += f_x;  // Accumulate the sum of function values
        sum_of_squares += f_x * f_x;  // Accumulate the sum of squares for variance calculation
    }

    // Compute the average of the function values
    result.integral = sum / N;

    // Estimate the error (standard deviation / sqrt(N))
    double variance = (sum_of_squares / N) - (result.integral * result.integral);
    result.error = sqrt(variance / N);

    return result;
}
