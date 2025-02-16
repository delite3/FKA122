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
    // Compute r_squared = x^2 + y^2 + z^2
    double r_squared = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    // Return normalized weight w(x, y, z) = (1 / (PI^(3/2))) * exp(-r^2)
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
    // Calculate f(x, y, z) = x^2 + x^2*y^2 + x^2*y^2*z^2
    double term1 = x[0] * x[0];
    double term2 = term1 * x[1] * x[1];
    double term3 = term2 * x[2] * x[2];
    return term1 + term2 + term3;
}

/* ***************************************
 *
 * Perform a single MCMC step with displacement
 *
 * - Use gsl_rng_uniform to displace the walker.
 * - The walker should be displaced in the order x, y, z.
 * - The draw for the acceptance condition should be done
 *   after the draws for the displacing of the walker.
 *
 * Parameters
 * ----------
 *  x - current walker position, size = 3
 *  delta - step size
 *  k - GSL random number generator object
 *
 * Returns
 * -------
 * - Updates the x parameter to reflect if the move was
 *   accepted or rejected.
 * - result contains:
 *   - The weight of the "exiting" x parameter.
 *   - The function value of the "exiting" x parameter.
 *   - `accepted` is 1 if the move was accepted, and
 *     0 if it was rejected.
 *
 * ***************************************/
result_t MCMC_step_displace_all(double *x, double delta, gsl_rng *k)
{
    result_t result;

    // Store the current position
    double x_old[3] = {x[0], x[1], x[2]};
    double weight_old = weight(x_old);

    // Propose a new position
    for (int i = 0; i < 3; i++) {
        double r = gsl_rng_uniform(k); // Random number in [0,1]
        x[i] = x[i] + delta * (r - 0.5);
    }

    // Calculate the weight of the new position
    double weight_new = weight(x);

    // Acceptance probability
    double acceptance = weight_new / weight_old;

    // Draw a random number to decide acceptance
    double r_accept = gsl_rng_uniform(k);

    if (r_accept < acceptance) {
        // Accept the move
        result.accepted = 1;
    } else {
        // Reject the move, revert to old position
        result.accepted = 0;
        for (int i = 0; i < 3; i++) {
            x[i] = x_old[i];
        }
    }

    // Store the results
    result.weight = weight(x);              // Weight of the (potentially reverted) position
    result.function_value = function(x);    // Function value at the position

    return result;
}

/* ***************************************
 *
 * Main function to test MCMC sampling
 *
 * ***************************************/
int main()
{
    // Initialize GSL random number generator
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, 12345); // Seed for reproducibility

    // Initial walker position
    double x[3] = {0.0, 0.0, 0.0};
    double delta = 1.0; // Step size

    // Perform MCMC steps
    int steps = 10000;
    double integral = 0.0;

    for (int i = 0; i < steps; i++) {
        result_t result = MCMC_step_displace_all(x, delta, rng);
        integral += result.function_value;
    }

    // Compute the average
    integral /= steps;

    // Print the result
    printf("Estimated integral: %f\n", integral);

    // Free GSL resources
    gsl_rng_free(rng);

    return 0;
}
