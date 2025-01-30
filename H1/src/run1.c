#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h" // Include your utility functions here

#define N 4        // 4x4x4 unit cells
#define NATOMS 4*4*4*N // Total number of atoms

void run(int argc, char *argv[]) {
    double latticeparam_start = 3.5; // Start lattice constant (Å)
    double latticeparam_end = 5.0;   // End lattice constant (Å)
    double latticeparam_step = 0.05; // Step size (Å)
    //double *masses = calloc(N, sizeof(double));
    double **positions = create_2D_array(NATOMS, 3);
    FILE *output = fopen("lattice_energy.txt", "w");

    // Loop over lattice constants
    for (double latticeparam = latticeparam_start; latticeparam <= latticeparam_end; latticeparam += latticeparam_step) {
        // Generate FCC lattice
        init_fcc(positions, N, latticeparam);

        // Compute potential energy
        double Epot = get_energy_AL(positions, N * latticeparam, NATOMS);
        
        // Calculate the unit cell volume
        double unit_cell_volume = pow(latticeparam, 3); // Unit cell volume (Å³)

        // Since energy is already for the entire system, calculate energy per unit cell
        double Epot_per_cell = Epot / (N * N * N); // Divide by number of unit cells
        
        // Write lattice constant and energy to file
        fprintf(output, "%f %f\n", unit_cell_volume, Epot_per_cell);
        printf("Lattice constant: %.2f Å, Energy: %.6f eV\n", unit_cell_volume, Epot_per_cell);
    }

    fclose(output);
}
