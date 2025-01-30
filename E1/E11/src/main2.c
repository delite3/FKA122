#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"

int main() {
    int num_particles = 3;  // CO2 has three particles
    const unsigned int num_timesteps = 1000;  // Number of time steps
    const double timestep = 0.01;  // Time step in atomic units
    const double kappa = 1.6e3;    // Spring constant in atomic units (N/m)
    const double alat = 10.0;      // Cell size for visualization (arbitrary units)

    double *positions_1D = malloc(num_particles * sizeof(double));
    double *velocities_1D = malloc(num_particles * sizeof(double));
    double *accelerations_1D = malloc(num_particles * sizeof(double));

    if (!positions_1D || !velocities_1D || !accelerations_1D) {
        perror("Memory allocation failed");
        return 1;
    }

    // Initialize values
    positions_1D[0] = 0.01;
    positions_1D[1] = 0.005;
    positions_1D[2] = -0.005;

    for (unsigned int i = 0; i < num_particles; i++) {
        velocities_1D[i] = 0.0;
        accelerations_1D[i] = 0.0;
    }
    // Preallocate memory for 3D positions and velocities (for ASE)
    double **positions_3D = create_2D_array(num_particles, 3);
    double **velocities_3D = create_2D_array(num_particles, 3);

    if (positions_3D == NULL || velocities_3D == NULL) {
        fprintf(stderr, "Failed to allocate memory for 3D positions or velocities.\n");
        return 1;
    }

    // Open the output file for visualization
    FILE *fp = fopen("CO2_simulation.extxyz", "w");
    if (fp == NULL) {
        perror("Failed to open output file.");
        return 1;
    }

    // Simulation loop
    for (unsigned int t = 0; t < num_timesteps; ++t) {
        // Write current 1D positions and velocities to 3D arrays
        for (unsigned int i = 0; i < num_particles; ++i) {
            positions_3D[i][0] = positions_1D[i];  // x-coordinate
            positions_3D[i][1] = 0.0;             // y-coordinate (dummy)
            positions_3D[i][2] = 0.0;             // z-coordinate (dummy)

            velocities_3D[i][0] = velocities_1D[i];  // x-coordinate
            velocities_3D[i][1] = 0.0;               // y-coordinate (dummy)
            velocities_3D[i][2] = 0.0;               // z-coordinate (dummy)
        }

        // Write to the .extxyz file
        write_xyz(fp, "C", positions_3D, velocities_3D, alat, num_particles);

        // Perform one velocity Verlet step
        velocity_verlet_one_step(
            accelerations_1D, positions_1D, velocities_1D, NULL, kappa, timestep
        );
    }

    // Free memory and close the file
    destroy_2D_array(positions_3D, num_particles);
    destroy_2D_array(velocities_3D, num_particles);
    fclose(fp);
    free(positions_1D);
    free(velocities_1D);
    free(accelerations_1D);

    return 0;
}
