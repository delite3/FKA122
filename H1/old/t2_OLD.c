#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"

#define N 4             // 4x4x4 unit cells
#define NATOMS 4*4*4*N  // Total number of atoms
#define DT 0.001        // Default time step
#define STEPS 10000     // Number of simulation steps
#define KB 8.617333262145e-5 // Boltzmann constant in eV/K
#define M 0.00284946023 // mass of aluminium in simulation units

void perturb_positions(double **positions, double **velocities, double latticeparam) {
    double max_displacement = 0.065 * latticeparam; // ±6.5% of lattice spacing
    for (unsigned int i = 0; i < NATOMS; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            double random_shift = (rand() / (double)RAND_MAX) * 2 * max_displacement - max_displacement;
            positions[i][j] += random_shift;
            velocities[i][j] = 0;
        }
    }
}

void run_simulation(double latticeparam, double timestep, const char *filename) {
    double **positions = create_2D_array(NATOMS, 3);
    double **velocities = create_2D_array(NATOMS, 3);
    double **forces = create_2D_array(NATOMS, 3);
    double *masses = calloc(NATOMS, sizeof(double));
    double Epot, Etot;
    
    for (unsigned int i = 0; i < NATOMS; i++) {
        masses[i] = M;
    }
    
    // Initialize FCC lattice
    init_fcc(positions, N, latticeparam);

    // Perturb positions
    perturb_positions(positions, velocities, latticeparam);

    // Open file to save results
    FILE *output = fopen(filename, "w");
    if (!output) {
        perror("Failed to open output file");
        exit(EXIT_FAILURE);
    }
    fprintf(output, "Time E_potential E_kinetic E_total Temperature\n");

    // Simulation loop
    for (unsigned int step = 0; step < STEPS; step++) {
        double time = step * timestep;

        // Get forces and potential energy
        get_forces_AL(forces, positions, N*latticeparam, NATOMS);
        Epot = get_energy_AL(positions, N*latticeparam, NATOMS);

        // Apply velocity Verlet for each dimension
        for (unsigned int dim = 0; dim < 3; dim++) {
            // Extract 1D arrays for the current dimension
            double forces_dim[NATOMS], positions_dim[NATOMS], velocities_dim[NATOMS];
            for (unsigned int i = 0; i < NATOMS; i++) {
                forces_dim[i] = forces[i][dim];
                positions_dim[i] = positions[i][dim];
                velocities_dim[i] = velocities[i][dim];
            }

            // Apply velocity Verlet in the current dimension
            velocity_verlet_one_step(forces_dim, positions_dim, velocities_dim, 0.0, timestep, NATOMS);

            // Update the original arrays with the updated values
            for (unsigned int i = 0; i < NATOMS; i++) {
                forces[i][dim] = forces_dim[i];
                positions[i][dim] = positions_dim[i];
                velocities[i][dim] = velocities_dim[i];
            }
        }
        double Ekin = 0;
        // calculate kinetic energy
        for (unsigned int dim = 0; dim < 3; dim++) {
            // Extract 1D velocity array for the current dimension
            double velocities_dim[NATOMS];
            for (unsigned int i = 0; i < NATOMS; i++) {
                velocities_dim[i] = velocities[i][dim];
                Ekin += calculate_kinetic_energy(velocities_dim, masses); 
            }
        }

        // Total energy
        Etot = Epot + Ekin;

        // Calculate temperature
        double temperature = (2.0 * Ekin) / (3.0 * NATOMS * KB);

        // Write results to file
        fprintf(output, "%f %f %f %f %f\n", time, Epot, Ekin, Etot, temperature);
    }

    fclose(output);

    // Clean up
    destroy_2D_array(positions, NATOMS);
    destroy_2D_array(velocities, NATOMS);
    destroy_2D_array(forces, NATOMS);
}

int run(int argc, char *argv[]) {
    double latticeparam = 4.04; // Lattice constant (Å)

    // Investigate timestep effects
    double timesteps[] = {0.0005, 0.0010, 0.0050}; // Test different timesteps

    char buf[50]; // Adjusted buffer size for safety

    for (unsigned int i = 0; i < sizeof(timesteps) / sizeof(timesteps[0]); i++) {
        // Create a unique filename for each timestep
        snprintf(buf, sizeof(buf), "energy_time_%.4f.txt", timesteps[i]);

        printf("Running simulation with DT = %.5f...\n", timesteps[i]);

        // Call the simulation function and pass the filename
        run_simulation(latticeparam, timesteps[i], buf);
    }

    return 0;
}
