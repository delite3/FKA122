#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "attempt.h"

#define PI 3.141592653589

result_t MC_with_importance_sampling(int N, gsl_rng *k)
{
    result_t result;
    result.integral = 0.0;
    result.error = 0.0;

    double sum = 0.0;
    double sum_of_squares = 0.0;

    // Make sure the random number generator is properly initialized
    gsl_rng_env_setup();

    // Loop over N points to perform Monte Carlo integration
    for (int i = 0; i < N; i++) {
        // Generate a random number u in [0, 1) to be used in transformation
        double u = gsl_rng_uniform(k);

        // Use the transformation method to generate x according to p(x) ∝ sin(πx)
        // This ensures that the generated points follow the desired probability density function
        double x = (1.0 / PI) * acos(1.0 - 2.0 * u);

        // Evaluate the function f(x) = x * (1 - x)
        double f_x = x * (1 - x);

        // Calculate the weight function p(x) = (π/2) * sin(π * x)
        double p_x = (PI / 2.0) * sin(PI * x);

        // Ensure p(x) > 0 to prevent division by zero or undefined behavior
        if (p_x > 0.0) {
            // Compute the weighted value f(x) / p(x) for the importance sampling estimate
            double weighted_value = f_x / p_x;

            // Accumulate the sum of weighted values for the integral estimate
            sum += weighted_value;

            // Accumulate the sum of squares of weighted values for variance calculation
            sum_of_squares += weighted_value * weighted_value;
        }
    }

    // After all points have been processed, calculate the integral and error estimates
    if (N > 0) {
        // Compute the Monte Carlo estimate of the integral
        result.integral = sum / N;

        // Compute the variance of the weighted values
        double variance = (sum_of_squares / N) - (result.integral * result.integral);

        // Ensure variance is non-negative (account for potential floating-point errors)
        if (variance < 0.0) {
            variance = 0.0;
        }

        // Compute the standard error of the integral estimate
        result.error = sqrt(variance / N);
    }

    return result;
}
