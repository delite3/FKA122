#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"
#include "lattice.h"
#include "potential.h"

// Constants
#define UNIT 4.0            // Unit Length
#define NATOMS 256          // Number of atoms
#define KB 8.617333262e-5   // Boltzmann constant (eV/K)
#define KT 1.0/760000.0     // AL isothermic compressibility as inverse Bulk Mod
#define EQ 1                // Run with equillibration
#define SIMTIME 20          // simulation time [ps]
#define BURNOFF 5        // Steps to skip for equilibration
#define MAL 26.0/9649.0     // Simulation-Unit mass

// isothermal compressibility for aluminium = 0.01385 GPa = 
#define IC 0.01385 

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

double get_kinetic(double **velocities, int n_atoms, double mass){
    double kinetic_energy = 0.0;
    for (int i = 0; i < n_atoms; i++) {
        for (int j = 0; j < 3; j++) {
            kinetic_energy += 0.5 * mass * velocities[i][j] * velocities[i][j];
        }
    }
    return kinetic_energy;
}

double get_kinetic_atom(double *atom_velocity, double mass){
    double kinetic_energy_atom = 0.0;
    for (int i = 0; i < 3; i++) {
        kinetic_energy_atom += 0.5 * mass * atom_velocity[i] * atom_velocity[i];
    }
    return kinetic_energy_atom;
}

double get_temp(double **velocities, int n_atoms, double mass) {
    
    double E_kin = get_kinetic(velocities, n_atoms, mass);
    
    return E_kin * 2.0/((double)3.0*KB*n_atoms);
}


double get_pressure(double T, double W, double V, int nAtoms){
    return ((nAtoms * KB * T + W) / V) * 1.602e6;   // Convert final result to bars
}

void verlet_step(double** positions, double** velocities, double** forces,
                 int N, double dt,  double m, double L){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < 3; j++){
            velocities[i][j] += dt*forces[i][j]/(2.0*m);
            positions[i][j] += dt*velocities[i][j];
        }
    }
    get_forces_AL(forces,positions, L, N);
    for(int i = 0; i < N; i++)
        for(int j = 0; j < 3; j++)
            velocities[i][j] += dt*forces[i][j]/(2.0*m);

}

double 
get_alphaT(double Teq, double temp, double dt){
    double tau = 100*dt;
    return 1 + 2*(dt/(tau))*(Teq-temp)/(temp); 
}

double 
get_alphaP(double Peq, double pressure, double dt){
    double tau = 300*dt;
    return 1 - KT*(dt/tau)*(Peq - pressure); 
}

void 
equillibriate(double Teq, double Peq, double T, double P, double *cell_length, double** positions, double** velocities, double dt){
    double aT = sqrt(get_alphaT(Teq, T, dt));
    double aP = cbrt(get_alphaP(Peq, P, dt));

    *cell_length = aP* (*cell_length);
    for (unsigned int i = 0; i < NATOMS; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            positions[i][j] *= aP;
            velocities[i][j] *= aT;
        }
    }
}

//################

void run_simulation(double Teq, double Peq, double latticeparam, double dt, const char *filename) {
    FILE *output = fopen(filename, "w");
    double TIME_STEPS = SIMTIME/dt;
    //fprintf(output, "time, potential_energy, kinetic_energy, total_energy, temperature, pressure, lattice_constant, sample_traj_x, sample_traj_y, sample_traj_z\n");

    // Initialize si mulation starting parameters
    double cell_length = UNIT * latticeparam;  // Simulation box length (Å)
    double cell_volume = pow(cell_length, 3); // Initial volume    

    // Allocate memory for positions, velocities, and forces
    double **positions = create_2D_array(NATOMS, 3);
    double **velocities = create_2D_array(NATOMS, 3);
    double **forces = create_2D_array(NATOMS, 3);
    double **sample_trajectory = create_2D_array(TIME_STEPS, 3);

    double kinetic_deviation_sq = 0.0;
    double pressure_sum = 0.0;
    double avg_temp = 0.0;
    
    // Perturbation and initialization
    init_fcc(positions, UNIT, latticeparam); // Initialize FCC lattice
    perturb_positions(positions, velocities, latticeparam); // Add small random displacements

    // Simulation loop
    for (int t = 0; t < TIME_STEPS; t++) {
        // Perform a velocity Verlet step
        verlet_step(positions, velocities, forces, NATOMS, dt, MAL, cell_length);

        cell_volume = pow(cell_length, 3);                                    // Update volume
        double virial = get_virial_AL(positions, cell_length, NATOMS);   // Compute virial
        double temperature = get_temp(velocities, NATOMS, MAL);              // Compute temperature
        double pressure = get_pressure(temperature, virial, cell_volume, NATOMS);            // Compute pressure

        if (EQ == 1){
            equillibriate(Teq, Peq, temperature, pressure, &cell_length, positions, velocities, dt);
        }

        double potential_energy = get_energy_AL(positions, cell_length, NATOMS);
        double kinetic_energy = get_kinetic(velocities, NATOMS, MAL);
        double total_energy = potential_energy + kinetic_energy;

        latticeparam = cell_length / UNIT;
        
        for (int i = 0; i < 3; i++) {
            sample_trajectory[t][i] = positions[0][i]; // Save position of one atom
        }
        
        //Compute fluctuations and average pressure after equilibration:

        double kinetic_mean = kinetic_energy / NATOMS;
        
        for(int i = 0; i < NATOMS; i++){
            double atom_kinetic = 0.0;
            
            atom_kinetic += get_kinetic_atom(velocities[i], MAL);

            double kinetic_deviation = kinetic_mean - atom_kinetic;
            kinetic_deviation_sq += kinetic_deviation*kinetic_deviation;
        }
        pressure_sum += pressure;
        avg_temp += temperature;


        fprintf(output, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                t*dt, potential_energy, kinetic_energy, total_energy, temperature, pressure, latticeparam,
                sample_trajectory[t][0], sample_trajectory[t][1], sample_trajectory[t][2]);
    }

    // Calculate average temperature and pressure
    avg_temp /= (TIME_STEPS - BURNOFF);
    pressure_sum /= (TIME_STEPS - BURNOFF);

    // Calculate heat capacity
    double heat_capacity = (3.0 * NATOMS * KB) / (1.0 - 2.0 * kinetic_deviation_sq / (3.0 * NATOMS * KB * KB * avg_temp * avg_temp));

    printf("Average temperature: %.2f K\n", avg_temp);
    printf("Average pressure: %.2f\n", pressure_sum);
    printf("Heat capacity: %f\n", heat_capacity);

    // Save simulation results to file

    fclose(output);

    // Free allocated memory
    destroy_2D_array(positions, NATOMS);
    destroy_2D_array(velocities, NATOMS);
    destroy_2D_array(forces, NATOMS);
    destroy_2D_array(sample_trajectory, NATOMS);

}


int run(int argc, char *argv[]) {
    double latticeparam = 4.03;               // Lattice parameter (Å)
    double Teq = 793.15;               
    double Peq = 1;               

    // Investigate timestep effects
    //double timesteps[] = {0.001, 0.005, 0.01, 0.015, 0.017}; // Test different timesteps
    double timesteps[] = {0.0005}; // Test different timesteps
    char filename[50]; // 

    for (unsigned int i = 0; i < sizeof(timesteps) / sizeof(timesteps[0]); i++) {
        // Create a unique filename for each timestep
        snprintf(filename, sizeof(filename), "energy_time_%.3f.txt", timesteps[i]);

        printf("Running simulation with DT = %.3f...\n", timesteps[i]);

        // Call the simulation function and pass the filename
        run_simulation(Teq, Peq, latticeparam, timesteps[i], filename);
    }

    return 0;
}
