#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define KB 1.0         // Boltzmann constant (normalized)
#define J 1.0          // Interaction energy (normalized)
#define T_START 0.1    // Starting temperature
#define T_END 5.0      // Ending temperature
#define T_STEP 0.1     // Temperature step
#define MAX_ITER 10000 // Maximum iterations for convergence
#define TOLERANCE 1e-6 // Convergence tolerance

// Function to compute magnetization m(T)
double find_m(double T) {
    double m = 0.99; // Initial guess close to 1
    double m_new;
    int iter = 0;

    do {
        m_new = tanh(J * m / (KB * T));
        if (fabs(m_new - m) < TOLERANCE) break;
        m = m_new;
        iter++;
    } while (iter < MAX_ITER);

    return m_new;
}

// Function to calculate energy U(T)
double calculate_energy(double m) {
    return -0.5 * J * m * m;
}

// Function to calculate heat capacity C(T)
double calculate_heat_capacity(double U1, double U2, double T_step) {
    return (U2 - U1) / T_step;
}

int main(void) {
    double T, m, U, C, U_prev = 0.0;
    FILE *output = fopen("results_mean_field.csv", "w");
    if (!output) {
        perror("Failed to open file");
        return 1;
    }

    fprintf(output, "T,m,U,C\n");

    for (T = T_START; T <= T_END; T += T_STEP) {
        // Find magnetization m(T)
        m = find_m(T);

        // Calculate U(T)
        U = calculate_energy(m);

        // Calculate C(T) using finite difference
        if (T > T_START) {
            C = calculate_heat_capacity(U_prev, U, T_STEP);
        } else {
            C = 0.0; // Heat capacity undefined at T=0
        }

        // Save results
        fprintf(output, "%.5f,%.5f,%.5f,%.5f\n", T, m, U, C);

        // Update previous U
        U_prev = U;
    }

    fclose(output);
    printf("Results written to results_mean_field.csv\n");
    return 0;
}
