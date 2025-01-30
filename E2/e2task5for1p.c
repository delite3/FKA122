#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>

#define PI 3.141592653589793

// Define the result_t struct
typedef struct {
    double weight;          // Weight of the current position
    double function_value;  // Value of the function at the current position
    int accepted;           // 1 if the move was accepted, 0 otherwise
} result_t;

/* ***************************************
 *
 * Calculate the normalized weight function w(x, y, z)
 *
 * Parameters
 * ----------
 *  x - current walker position, size = 3
 *
 * Returns
 * -------
 * The weight of the current position
 *
 * ***************************************/
double weight(double *x)
{
    double r_squared = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    return pow(PI, -1.5) * exp(-r_squared);
}

/* ***************************************
 *
 * Calculate the function to be sampled f(x, y, z)
 *
 * Parameters
 * ----------
 *  x - current walker position, size = 3
 *
 * Returns
 * -------
 * The function value
 *
 * ***************************************/
double function(double *x)
{
    double term1 = x[0] * x[0];
    double term2 = term1 * x[1] * x[1];
    double term3 = term2 * x[2] * x[2];
    return term1 + term2 + term3;
}

/* ***************************************
 *
 * Perform a single MCMC step with displacement
 *
 * Parameters
 * ----------
 *  x - current walker position, size = 3
 *  delta - step size
 *  k - GSL random number generator object
 *
 * Returns
 * -------
 * Updates x and returns the result of the step.
 *
 * ***************************************/
result_t MCMC_step_displace_all(double *x, double delta, gsl_rng *k)
{
    result_t result;

    double x_old[3] = {x[0], x[1], x[2]};
    double weight_old = weight(x_old);

    // Propose a new position
    for (int i = 0; i < 3; i++) {
        double r = gsl_rng_uniform(k);
        x[i] = x[i] + delta * (r - 0.5);
    }

    double weight_new = weight(x);
    double acceptance = weight_new / weight_old;
    double r_accept = gsl_rng_uniform(k);

    if (r_accept < acceptance) {
        result.accepted = 1;
    } else {
        result.accepted = 0;
        for (int i = 0; i < 3; i++) {
            x[i] = x_old[i];
        }
    }

    result.weight = weight(x);
    result.function_value = function(x);
    return result;
}

/* ***************************************
 *
 * Main program to compute the integral with dynamic adjustment of delta
 *
 * ***************************************/
int main()
{
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, 12345); // Seed for reproducibility

    double x[3] = {0.0, 0.0, 0.0}; // Initial position
    double delta = 1.0;            // Initial step size
    int steps = 10000;             // Number of steps for integral calculation
    int warmup_steps = 1000;       // Steps to adjust delta
    double target_ratio = 0.5;     // Desired acceptance ratio

    // Warmup phase: Adjust delta dynamically
    for (int i = 0; i < warmup_steps; i++) {
        int accepted_steps = 0;

        for (int j = 0; j < 100; j++) { // Small batch of steps
            result_t result = MCMC_step_displace_all(x, delta, rng);
            if (result.accepted) {
                accepted_steps++;
            }
        }

        double acceptance_ratio = (double)accepted_steps / 100;
        if (acceptance_ratio > target_ratio + 0.05) {
            delta *= 1.1; // Increase delta to reduce acceptance ratio
        } else if (acceptance_ratio < target_ratio - 0.05) {
            delta *= 0.9; // Decrease delta to increase acceptance ratio
        }
        printf("Warmup Step %d: Delta = %f, Acceptance Ratio = %f\n", i, delta, acceptance_ratio);
    }

    // Sampling phase: Perform MCMC and compute the integral
    double integral = 0.0;
    int accepted_steps = 0;

    for (int i = 0; i < steps; i++) {
        result_t result = MCMC_step_displace_all(x, delta, rng);
        if (result.accepted) {
            accepted_steps++;
        }
        integral += result.function_value;
    }

    integral /= steps;

    // Report results
    double final_acceptance_ratio = (double)accepted_steps / steps;
    printf("Final Delta: %f\n", delta);
    printf("Final Acceptance Ratio: %f\n", final_acceptance_ratio);
    printf("Estimated Integral: %f\n", integral);

    gsl_rng_free(rng);
    return 0;
}
