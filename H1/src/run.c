#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"
#include "lattice.h"
#include "potential.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

// program definitions
#define EQUALIZE 1          // Run with equillibration
#define LIQUIDIZE 1         // Liquize parameters
// Constants
#define UNIT 4.0            // Unit Length
#define NSUPER 256          // Number of supercells
#define MAL 26.0/9649.0     // Simulation-Unit mass 
#define KB 8.617333262e-5   // Boltzmann constant [eV/K]
#define KT 1.0/760000.0     // AL isothermic compressibility as [Bulk Modulus]^-1
#define SIMTIME 10          // simulation time [ps]
#define TMELT 2             // Initial melt time [ps]
#define TSTABLE 3           // Time for simulation to stabilize [ps]
#define PI 3.14159265359


void perturb_positions(double** positions, double**velocities, double latticeparam){
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup ();
    T = gsl_rng_default ;
    r = gsl_rng_alloc (T);
    int seed = time(NULL);
    //int seed = 42;
    gsl_rng_set (r, seed);

    double max_displacement = 0.065 * latticeparam; // ±6.5% of lattice spacing
    for(int i = 0; i < NSUPER; i++){
        for(int j=0; j < 3;j++){
            double ran_pos = gsl_ran_flat(r, -max_displacement, max_displacement);
            positions[i][j] += ran_pos;
            //printf("perturb: %f\n", ran_pos);
            velocities[i][j] = 0.0;
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

double get_temp(double **velocities, int n_atoms, double mass) {
    
    double E_kin = get_kinetic(velocities, n_atoms, mass);
    
    return E_kin * 2.0/((double)3.0*KB*n_atoms);
}


double get_pressure(double T, double W, double V, int nAtoms){
    return ((nAtoms * KB * T + W) / V) * 1602000;   // Convert final result to bars
}

void verlet_step(double** positions, double** velocities, double** forces,
                 int N, double dt,  double m, double L){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < 3; j++){
            velocities[i][j] += dt*forces[i][j]/(2.0*m);
            positions[i][j] += dt*velocities[i][j];
        }
    }
    get_forces_AL(forces, positions, L, N);
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
    for (unsigned int i = 0; i < NSUPER; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            positions[i][j] *= aP;
            velocities[i][j] *= aT;
        }
    }
}

//################

void run_simulation(double Teq, double Peq, double latticeparam, double dt, const char *filename1, const char *filename2) {
    FILE *output1 = fopen(filename1, "w");
    FILE *output2 = fopen(filename2, "w");
    int TIME_STEPS = SIMTIME/dt;
    //fprintf(output, "time, potential_energy, kinetic_energy, total_energy, temperature, pressure, lattice_constant, sample_traj_x, sample_traj_y, sample_traj_z\n");

    // Initialize simulation starting parameters
    double cell_length      = UNIT * latticeparam;  // Simulation box length (Å)
    double cell_volume      = pow(cell_length, 3); // Initial volume    
    double virial           = 0; 
    double temperature      = 0;
    double pressure         = 0;
    double potential_energy = 0;
    double kinetic_energy   = 0;
    double total_energy     = 0;

    // Allocate memory for positions, velocities, and forces
    double **positions = create_2D_array(NSUPER, 3);
    double **velocities = create_2D_array(NSUPER, 3);
    double **forces = create_2D_array(NSUPER, 3);
    
    // Perturbation and initialization
    init_fcc(positions, UNIT, latticeparam); // Initialize FCC lattice
    perturb_positions(positions, velocities, latticeparam); // Add small random displacements
    
    //Task 4 Liquidization
    double TeqO = 1;

    // Task 6 parameters
    int k               = 200;
    double r_max        = cell_length/2;
    double delta_r      = r_max/k;
    double *bin_counts  = calloc(k, sizeof(double));
    double density = NSUPER / cell_volume;

    // Simulation loop
    for (int t = 0; t < TIME_STEPS; t++) {
        
        double ps = t*dt;
        if (LIQUIDIZE == 1 && ps < TMELT){
            TeqO = 2*Teq;
        } else {
            TeqO = Teq;
        }
        
        // Perform a velocity Verlet step
        verlet_step(positions, velocities, forces, NSUPER, dt, MAL, cell_length);

        cell_volume  = pow(cell_length, 3);                                      // Update volume
        virial       = get_virial_AL(positions, cell_length, NSUPER);            // Compute virial
        temperature  = get_temp(velocities, NSUPER, MAL);                        // Compute temperature
        pressure     = get_pressure(temperature, virial, cell_volume, NSUPER);   // Compute pressure

        if (EQUALIZE == 1) equillibriate(TeqO, Peq, temperature, pressure, &cell_length, positions, velocities, dt);

        potential_energy = get_energy_AL(positions, cell_length, NSUPER);
        kinetic_energy   = get_kinetic(velocities, NSUPER, MAL);
        total_energy     = potential_energy + kinetic_energy;

        latticeparam    = cell_length / UNIT;
        r_max           = cell_length/2;
        delta_r         = r_max/k;

        if (t> (TMELT+TSTABLE)/dt){
            density += NSUPER / cell_volume;
            for (int i = 0; i < NSUPER; i++) {
                for (int j = i+1; j < NSUPER; j++) {

                    double dx  = positions[j][0] - positions[i][0];
                    double dy  = positions[j][1] - positions[i][1];
                    double dz  = positions[j][2] - positions[i][2];
                    
                    // Minimum image convention
                    if (dx > r_max)         dx -= cell_length;
                    else if (dx < -r_max)   dx += cell_length;

                    if (dy > r_max)         dy -= cell_length;
                    else if (dy < -r_max)   dy += cell_length;

                    if (dz > r_max)         dz -= cell_length;
                    else if (dz < -r_max)   dz += cell_length;

                    double r = sqrt(dx*dx + dy*dy + dz*dz);

                    int bin = (int)(r / delta_r);
                    if (bin < k) {
                        bin_counts[bin] += 2.0; // each pair contributes to both i an j surroundings
                    }
                }
            }
        }
        // Save simulation results to file
        fprintf(output1, "%f,%f,%f,%f,%f,%f,%f\n",
                t*dt, potential_energy, kinetic_energy, total_energy, temperature, pressure, latticeparam);

        for (int i = 0; i < NSUPER; i++) {
            fprintf(output2, "%f,%f,%f,%f,%f\n",
                    t*dt, cell_length, positions[i][0], positions[i][1], positions[i][2]);}

    }

    density /= (SIMTIME-TMELT-TSTABLE)/dt;

    FILE *fp = fopen("gofr.txt", "w");
    for (int bin = 0; bin < k; bin++) {
        double r_bin = (bin + 0.5)*delta_r; 
        double shell_volume = 4.0 * PI * (r_bin*r_bin)*delta_r;
        double ideal_count = density * shell_volume * NSUPER * (SIMTIME-TMELT-TSTABLE)/dt;
        double g_of_r = bin_counts[bin]/ideal_count;

        fprintf(fp, "%f,%f,%f\n", r_bin, g_of_r, density);
    }
    
    fclose(fp);

    fclose(output1);
    fclose(output2);

    // Free allocated memory
    destroy_2D_array(positions, NSUPER);
    destroy_2D_array(velocities, NSUPER);
    destroy_2D_array(forces, NSUPER);
    
}

int run(int argc, char *argv[]) {
    double latticeparam = 4.03;         // Lattice parameter (Å)
    double Teq = 773.15;                // target temperature [K]
    double Peq = 1;                     // target pressure [bar] 

    if (LIQUIDIZE == 1){
            Teq = 973.15;
    }

    // Investigate timestep effects
    //double timesteps[] = {0.0001, 0.0005, 0.001, 0.005, 0.01}; // Test different timesteps
    double timesteps[] = {0.0005}; // Test different timesteps
    char filename1[50]; // 
    char filename2[50]; // 

    for (unsigned int i = 0; i < sizeof(timesteps) / sizeof(timesteps[0]); i++) {
        // Create a unique filename for each timestep
        if (LIQUIDIZE){
            snprintf(filename1, sizeof(filename1), "Equilibrate_Liquid_dt%.4f.txt", timesteps[i]);
            snprintf(filename2, sizeof(filename2), "Positions_Liquid_dt%.4f.txt", timesteps[i]);
        } else{
            snprintf(filename1, sizeof(filename1), "Equilibrate_Solid_dt%.4f.txt", timesteps[i]);
            snprintf(filename2, sizeof(filename2), "Positions_Solid_dt%.4f.txt", timesteps[i]);
        }

        printf("Running simulation %s...\n", filename1);

        // Call the simulation function and pass the filename
        run_simulation(Teq, Peq, latticeparam, timesteps[i], filename1, filename2);
    }

    return 0;
}
