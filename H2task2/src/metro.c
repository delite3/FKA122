#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define L 10
#define N (L * L * L * 2) // Total number of atoms (2 sublattices)
#define NEIGHBORS 8
#define K_B 1.0 // Boltzmann constant

// Bond energies (in meV)
#define E_AA -436
#define E_BB -113
#define E_AB -294

// Temperatures in Kelvin
#define TEMP1 400
#define TEMP2 600
#define TEMP3 1000

// Simulation parameters
#define N_STEPS 1000000
#define N_EQ 100000

// Function prototypes
void initialize_structure(int *lattice, int ordered);
double calculate_energy(int *lattice, int neighbors[N][NEIGHBORS]);
void create_neighbor_list(int neighbors[N][NEIGHBORS]);
int metropolis_step(int *lattice, int neighbors[N][NEIGHBORS], double T);
void run_simulation(int *lattice, int neighbors[N][NEIGHBORS], double temperatures[]);

int main() {
    int lattice[N];
    int neighbors[N][NEIGHBORS];
    double temperatures[] = {TEMP1, TEMP2, TEMP3};

    // Initialize neighbor list for the BCC structure
    create_neighbor_list(neighbors);

    // Cold start
    printf("Cold Start Simulation:\n");
    initialize_structure(lattice, 1);
    run_simulation(lattice, neighbors, temperatures);

    // Warm start
    printf("\nWarm Start Simulation:\n");
    initialize_structure(lattice, 0);
    run_simulation(lattice, neighbors, temperatures);

    return 0;
}

void 
run_simulation(int *lattice, int neighbors[N][NEIGHBORS], double temperatures[]) {
    double energy; 
    double total_energy;

    for (int t = 0; t < 3; t++) {
        double T = temperatures[t];

        energy = calculate_energy(lattice, neighbors);
        printf("Initial energy at T=%.1f K: %.2f meV\n", T, energy);

        // Equilibration steps
        for (int i = 0; i < N_EQ; i++) {
            metropolis_step(lattice, neighbors, T);
        }

        // Production steps
        total_energy = 0.0;
        for (int i = 0; i < N_STEPS; i++) {
            metropolis_step(lattice, neighbors, T);
            total_energy += calculate_energy(lattice, neighbors);
        }

        // Average energy
        double avg_energy = total_energy / N_STEPS;
        printf("Average energy at T=%.1f K: %.2f meV\n", T, avg_energy);
    }
}

void initialize_structure(int *lattice, int ordered) {
    for (int i = 0; i < N; i++) {
        lattice[i] = (ordered ? (i % 2) : (rand() % 2)); // 0 for A, 1 for B
    }
}

double calculate_energy(int *lattice, int neighbors[N][NEIGHBORS]) {
    double energy = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < NEIGHBORS; j++) {
            int neighbor = neighbors[i][j];
            if (lattice[i] == 0 && lattice[neighbor] == 0) {
                energy += E_AA;
            } else if (lattice[i] == 1 && lattice[neighbor] == 1) {
                energy += E_BB;
            } else {
                energy += E_AB;
            }
        }
    }
    return energy / 2.0; // Each bond is counted twice
}

void create_neighbor_list(int neighbors[N][NEIGHBORS]) {
    int index = 0;
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            for (int z = 0; z < L; z++) {
                // Sublattice a
                neighbors[index][0] = (((x + 1) % L) * L * L + y * L + z) * 2 + 1;
                neighbors[index][1] = (((x - 1 + L) % L) * L * L + y * L + z) * 2 + 1;
                neighbors[index][2] = (x * L * L + ((y + 1) % L) * L + z) * 2 + 1;
                neighbors[index][3] = (x * L * L + ((y - 1 + L) % L) * L + z) * 2 + 1;
                neighbors[index][4] = (x * L * L + y * L + ((z + 1) % L)) * 2 + 1;
                neighbors[index][5] = (x * L * L + y * L + ((z - 1 + L) % L)) * 2 + 1;
                neighbors[index][6] = (((x + 1) % L) * L * L + ((y + 1) % L) * L + ((z + 1) % L)) * 2 + 1;
                neighbors[index][7] = (((x - 1 + L) % L) * L * L + ((y - 1 + L) % L) * L + ((z - 1 + L) % L)) * 2 + 1;

                index++;

                // Sublattice b
                neighbors[index][0] = (((x + 1) % L) * L * L + y * L + z) * 2;
                neighbors[index][1] = (((x - 1 + L) % L) * L * L + y * L + z) * 2;
                neighbors[index][2] = (x * L * L + ((y + 1) % L) * L + z) * 2;
                neighbors[index][3] = (x * L * L + ((y - 1 + L) % L) * L + z) * 2;
                neighbors[index][4] = (x * L * L + y * L + ((z + 1) % L)) * 2;
                neighbors[index][5] = (x * L * L + y * L + ((z - 1 + L) % L)) * 2;
                neighbors[index][6] = (((x + 1) % L) * L * L + ((y + 1) % L) * L + ((z + 1) % L)) * 2;
                neighbors[index][7] = (((x - 1 + L) % L) * L * L + ((y - 1 + L) % L) * L + ((z - 1 + L) % L)) * 2;

                index++;
            }
        }
    }
}

int metropolis_step(int *lattice, int neighbors[N][NEIGHBORS], double T) {
    int i = rand() % N;
    int j = rand() % N;
    if (i == j) return 0;

    int delta_energy = 0;
    for (int k = 0; k < NEIGHBORS; k++) {
        int neighbor = neighbors[i][k];
        if (lattice[i] == 0 && lattice[neighbor] == 0) delta_energy -= E_AA;
        if (lattice[i] == 1 && lattice[neighbor] == 1) delta_energy -= E_BB;
        if (lattice[i] != lattice[neighbor]) delta_energy -= E_AB;
    }

    lattice[i] = 1 - lattice[i]; // Swap atom type

    for (int k = 0; k < NEIGHBORS; k++) {
        int neighbor = neighbors[i][k];
        if (lattice[i] == 0 && lattice[neighbor] == 0) delta_energy += E_AA;
        if (lattice[i] == 1 && lattice[neighbor] == 1) delta_energy += E_BB;
        if (lattice[i] != lattice[neighbor]) delta_energy += E_AB;
    }

    if (delta_energy <= 0 || (rand() / (double)RAND_MAX) < exp(-delta_energy / (K_B * T))) {
        return 1; // Move accepted
    } else {
        lattice[i] = 1 - lattice[i]; // Revert swap
        return 0; // Move rejected
    }
}