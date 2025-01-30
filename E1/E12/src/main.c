#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"

#define ALPHA 0.1
#define N 32       // Number of particles
#define E_0 32.0   // Initial energy
#define DT 0.1     // Time step
#define T_MAX 1000000 // Maximum simulation time
#define LONGRUN 1

void initialize_conditions(double *positions, double *velocities, double *accelerations, double *Q, double *P) {
    // Initialize Q_k and P_k
    for (unsigned int k = 0; k < N; k++) {
        Q[k] = 0.0;
        P[k] = 0.0; // Only k=1 mode has kinetic energy
    }
    P[0] = sqrt(2.0 * E_0);

    // Transform normal modes to ordinary coordinates
    transform_to_normal_modes(Q, positions, N);

    // Set velocities based on P_k
    transform_to_normal_modes(P, velocities, N);

    // Initialize accelerations
    calculate_acceleration(accelerations, positions, ALPHA, N); 
}

int main() {
    double *positions = calloc(N, sizeof(double));
    double *velocities = calloc(N, sizeof(double));
    double *accelerations = calloc(N, sizeof(double));
    double *Q = calloc(N, sizeof(double));  // Normal mode displacements
    double *P = calloc(N, sizeof(double));  // Normal mode momenta
    double energies[N];                     // Energies of normal modes
    double time = 0.0;
    double energy_sum = 0;
    
    FILE *output = fopen("energies.txt", "w");

    // Initialize system
    initialize_conditions(positions, velocities, accelerations, Q, P);

    unsigned int steps = T_MAX / DT;
    for (unsigned int step = 0; step < steps; step++) {
        // Perform one velocity Verlet step
        velocity_verlet_one_step(accelerations, positions, velocities, ALPHA, DT, N);

        // Calculate normal mode energies and save to file
        if (LONGRUN == 1) {
            if (step % 1000 == 0) {
                calculate_normal_mode_energies(energies, positions, velocities, N);
                
                // Output all normal mode energies
                fprintf(output, "%f", time); // Write the current time
                for (unsigned int k = 0; k < N; k++) {
                    fprintf(output, " %f", energies[k]); // Write the energy of mode k
                    energy_sum += energies[k];
                }
                fprintf(output, " %f", energy_sum); // Write the energy of mode k
                fprintf(output, "\n"); // End the line
                energy_sum = 0;
            }
            
        }
        else {
            calculate_normal_mode_energies(energies, positions, velocities, N);
            fprintf(output, "%f %f %f %f %f %f\n", time, energies[0], energies[1], energies[2], energies[3], energies[4]);
        }

        // Increment time
        time += DT;
    }

    fclose(output);
    free(positions);
    free(velocities);
    free(accelerations);
    free(Q);
    free(P);

    printf("Simulation complete. Results saved to energies.txt\n");
    return 0;
}
