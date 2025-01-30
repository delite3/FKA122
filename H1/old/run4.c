#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"
#include "lattice.h"
#include "potential.h"

#define KB 8.617333262145e-5    // Boltzmann constant in eV/K
#define M 26.0/9649             // Mass in atomic simulation units

void perturb_positions(double **positions, double **velocities, double latticeparam, int N) {
    double max_displacement = 0.065 * latticeparam; // ±6.5% of lattice spacing
    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            double random_shift = (rand() / (double)RAND_MAX) * 2 * max_displacement - max_displacement;
            positions[i][j] += random_shift;
            velocities[i][j] = 0;
        }
    }
}

void velocity_verlet_step(double **positions, double **velocities, double **forces, 
                          double cell_length, double timestep, int N) {
    // Half-step velocity update
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            velocities[i][j] += 0.5 * forces[i][j] * timestep;
        }
    }

    // Full-step position update
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            positions[i][j] += velocities[i][j] * timestep;
        }
    }

    // Recalculate forces
    get_forces_AL(forces, positions, cell_length, N);

    // Half-step velocity update (final)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            velocities[i][j] += 0.5 * forces[i][j] * timestep;
        }
    }
}

double get_kinetic(double **velocities, int n_atoms, double mass){
    double kinetic_energy = 0.0;
    for (int i = 0; i < n_atoms; i++) {
        for (int j = 0; j < 3; j++) {
            kinetic_energy += 0.5 * mass * velocities[i][j] * velocities[i][j];
        }
    }
    return kinetic_energy;
}

double calculate_temperature(double **velocities, int N) {
    double kinetic_energy = get_kinetic(velocities, N, M);

    // Temperature using equipartition theorem
    return (2.0 * kinetic_energy) / ((double)3.0 * N * KB);
}

double calculate_pressure(double temperature, double virial, double volume, int N) {
    return (N * KB * temperature + virial) / volume;
}


void run(int argc, char *argv[]) {
    const int N = 256;              // Number of atoms
    const double a0 = 4.04;         // lattice parameter 
    double cell_length = a0*4;      // Cell length (Å)
    const double timestep = 0.001;  // Time step (ps)
    const int steps = 10000;        // Number of steps
    const double volume = pow(cell_length, 3); // Simulation box volume

    double **positions = create_2D_array(N, 3);
    double **velocities = create_2D_array(N, 3);
    double **forces = create_2D_array(N, 3);

    double potential_energy, virial;

    // Initialize positions, velocities, and forces
    init_fcc(positions, 4, a0);     // FCC lattice with 4x4x4 cells
    perturb_positions(positions, velocities, a0, N);

    get_forces_AL(forces, positions, cell_length, N);

    FILE *output = fopen("output.txt", "w");
    fprintf(output, "Time PotentialEnergy KineticEnergy TotalEnergy Temperature Pressure\n");

    for (int step = 0; step < steps; step++) {
        double time = step * timestep;

        // Perform a velocity Verlet step
        velocity_verlet_step(positions, velocities, forces, cell_length, timestep, N);

        potential_energy = get_energy_AL(positions, cell_length, N);
        
        // Calculate temperature
        double temperature = calculate_temperature(velocities, N);

        // Calculate virial
        virial = get_virial_AL(positions, cell_length, N);

        // Calculate pressure
        double pressure = calculate_pressure(temperature, virial, volume, N);

        // Calculate total energy
        double kinetic_energy = get_kinetic(velocities, N, M);
        double total_energy = potential_energy + kinetic_energy;

        // Save results to file
        fprintf(output, "%f %f %f %f %f %f\n", time, potential_energy, kinetic_energy, total_energy, temperature, pressure);
    }

    fclose(output);

    destroy_2D_array(positions, N);
    destroy_2D_array(velocities, N);
    destroy_2D_array(forces, N);
}
