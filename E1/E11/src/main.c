#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"  // Include the header file with function declarations

#define TIMESTEPS 10000  // Number of timesteps
#define DT 0.0001        // Time step [ps]
#define KAPPA 99.8641    // Spring constant [eV/Å^2]

// Function prototypes
void write_positions(FILE *file, double time, double *positions);
void write_energies(FILE *file, double time, double potential, double kinetic, double total);

int main() {
    
    // Define masses (in atomic units)
    double masses[3] = {0.00165820292, 0.00124365219, 0.00165820292};  // C, O, C masses in asu [eV/(ps^2 * Å^2)]

    // Initial positions (in Ångströms)
    double positions[3] = {0.01, 0.005, -0.005};

    // Initial velocities (all zero)
    double velocities[3] = {0.0, 0.0, 0.0};

    // Accelerations (initialized to zero)
    double accelerations[3] = {0.0, 0.0, 0.0};

    // Initialize energies
    double potential_energy = calculate_potential_energy(positions, KAPPA);
    double kinetic_energy = calculate_kinetic_energy(velocities, masses);
    double total_energy = potential_energy + kinetic_energy;

    // Open files to save data
    FILE *pos_file = fopen("positions.txt", "w");
    FILE *energy_file = fopen("energies.txt", "w");
    if (!pos_file || !energy_file) {
        perror("Failed to open output files.");
        exit(EXIT_FAILURE);
    }

    // Time integration loop
    for (int step = 0; step < TIMESTEPS; step++) {
        double time = step * DT;

        // Write positions and energies to files
        write_positions(pos_file, time, positions);
        write_energies(energy_file, time, potential_energy, kinetic_energy, total_energy);

        // Perform a velocity Verlet step
        velocity_verlet_one_step(accelerations, positions, velocities, masses, KAPPA, DT);

        // Update energies
        potential_energy = calculate_potential_energy(positions, KAPPA);
        kinetic_energy = calculate_kinetic_energy(velocities, masses);
        total_energy = potential_energy + kinetic_energy;
    }

    // Close files
    fclose(pos_file);
    fclose(energy_file);

    printf("Simulation complete. Results written to positions.txt and energies.txt\n");
    return 0;
}

// Helper function to write positions to a file
void write_positions(FILE *file, double time, double *positions) {
    fprintf(file, "%lf %lf %lf %lf\n", time, positions[0], positions[1], positions[2]);
}

// Helper function to write energies to a file
void write_energies(FILE *file, double time, double potential, double kinetic, double total) {
    fprintf(file, "%lf %lf %lf %lf\n", time, potential, kinetic, total);
}
