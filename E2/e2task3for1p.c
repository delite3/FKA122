#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>

#define PI 3.141592653589

/* Struct to store the result of the Monte Carlo integration */
typedef struct {
    double integral;
    double error;
} result_t;

/* Monte Carlo integration without importance sampling */
result_t MC_without_importance_sampling(int N, gsl_rng *k)
{
    result_t result;
    result.integral = 0.0;
    result.error = 0.0;
    FILE *output = fopen("samples.txt", "w");
    double sum = 0.0;
    double sum_of_squares = 0.0;

    for (int i = 0; i < N; i++) {
        double x = gsl_rng_uniform(k);
        double f_x = x * (1 - x);

        sum += f_x;
        sum_of_squares += f_x * f_x;
        fprintf(output, "%f %f %f %f\n", x, f_x, sum, sum_of_squares);
    }

    result.integral = sum / N;

    double variance = (sum_of_squares / N) - (result.integral * result.integral);
    result.error = sqrt(variance / N);

    

    return result;
}

/* Monte Carlo integration with importance sampling */
result_t MC_with_importance_sampling(int N, gsl_rng *k)
{
    result_t result;
    result.integral = 0.0;
    result.error = 0.0;

    double sum = 0.0;
    double sum_of_squares = 0.0;

    for (int i = 0; i < N; i++) {
        double u = gsl_rng_uniform(k);
        double x = (1.0 / PI) * acos(1.0 - 2.0 * u);
        double f_x = x * (1 - x);
        double p_x = (PI / 2.0) * sin(PI * x);

        double weighted_value = f_x / p_x;
        sum += weighted_value;
        sum_of_squares += weighted_value * weighted_value;
    }

    result.integral = sum / N;

    double variance = (sum_of_squares / N) - (result.integral * result.integral);
    result.error = sqrt(variance / N);

    return result;
}

int main()
{
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

    int N_values[] = {10, 100, 1000, 10000};
    double exact_value = 1.0 / 6.0;

    printf("N\tWithout Importance Sampling\tError\tWith Importance Sampling\tError\n");
    printf("------------------------------------------------------------------------\n");

    for (int i = 0; i < 4; i++) {
        int N = N_values[i];

        result_t no_importance = MC_without_importance_sampling(N, r);
        result_t with_importance = MC_with_importance_sampling(N, r);

        printf("%d\t%.6f\t\t\t%.6f\t%.6f\t\t\t%.6f\n",
               N, no_importance.integral, fabs(no_importance.integral - exact_value),
               with_importance.integral, fabs(with_importance.integral - exact_value));
    }

    gsl_rng_free(r);
    return 0;
}
