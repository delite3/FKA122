#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h" // Include the header file for the provided functions

void initialize_conditions(double *positions, double *velocities, double *accelerations, double *Q, double *P, unsigned int N, double E_0) {
    // Initialize Q_k and P_k
    for (unsigned int k = 0; k < N; k++) {
        Q[k] = 0.0; // No displacement in any mode initially
        P[k] = 0.0; // No momentum in any mode except k=1
    }
    P[0] = sqrt(2.0 * E_0); // Energy localized in mode k=1

    // Transform Q_k and P_k back to positions and velocities
    transform_to_normal_modes(Q, positions, N);
    transform_to_normal_modes(P, velocities, N);

    // Initialize accelerations
    calculate_acceleration(accelerations, positions, 0.0, N); // alpha = 0 for harmonic case
}

void simulate_motion(double *positions, double *velocities, double *accelerations, double *energies, unsigned int N, double dt, double t_max) {
    unsigned int steps = (unsigned int)(t_max / dt);
    FILE *energy_output = fopen("mode_energies.txt", "w");

    if (!energy_output) {
        perror("Failed to open file for writing");
        exit(EXIT_FAILURE);
    }

    for (unsigned int step = 0; step < steps; step++) {
        double time = step * dt;

        // Perform one velocity Verlet step
        velocity_verlet_one_step(accelerations, positions, velocities, 0.0, dt, N); // alpha = 0 for harmonic case

        // Calculate normal mode energies
        calculate_normal_mode_energies(energies, positions, velocities, N);

        // Write energies to file for the first 5 modes
        fprintf(energy_output, "%f", time); // Write time
        for (unsigned int k = 0; k < 5; k++) {
            fprintf(energy_output, " %f", energies[k]);
        }
        fprintf(energy_output, "\n");
    }

    fclose(energy_output);
}

int main() {
    // Parameters
    const unsigned int N = 32; // Number of particles
    const double E_0 = N;      // Total initial energy
    const double dt = 0.1;     // Time step
    const double t_max = 25000.0; // Maximum simulation time

    // Allocate memory
    double *positions = (double *)calloc(N, sizeof(double));
    double *velocities = (double *)calloc(N, sizeof(double));
    double *accelerations = (double *)calloc(N, sizeof(double));
    double *energies = (double *)calloc(N, sizeof(double));
    double *Q = (double *)calloc(N, sizeof(double));
    double *P = (double *)calloc(N, sizeof(double));

    if (!positions || !velocities || !accelerations || !energies || !Q || !P) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    // Initialize conditions
    initialize_conditions(positions, velocities, accelerations, Q, P, N, E_0);

    // Simulate motion
    simulate_motion(positions, velocities, accelerations, energies, N, dt, t_max);

    // Clean up
    free(positions);
    free(velocities);
    free(accelerations);
    free(energies);
    free(Q);
    free(P);

    printf("Simulation complete. Energies written to mode_energies.txt.\n");
    return 0;
}
