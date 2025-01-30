#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define L 10 // Lattice size
#define N (L * L * L) // Total number of atoms
#define k_B 8.617e-5 // Boltzmann constant in eV/K

// Bond energies in meV
#define E_AA -436
#define E_BB -113
#define E_AB -294

// Function to initialize the lattice
void initialize_lattice(int lattice[L][L][L], int ordered) {
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            for (int z = 0; z < L; z++) {
                lattice[x][y][z] = (x + y + z) % 2;
            }
        }
    }

    if (!ordered) {
        for (int i = 0; i < N; i++) {
            int x1 = rand() % L;
            int y1 = rand() % L;
            int z1 = rand() % L;
            int x2 = rand() % L;
            int y2 = rand() % L;
            int z2 = rand() % L;
            int temp = lattice[x1][y1][z1];
            lattice[x1][y1][z1] = lattice[x2][y2][z2];
            lattice[x2][y2][z2] = temp;
        }
    }
}

// Function to calculate the energy of the system
double calculate_energy(int lattice[L][L][L]) {
    double energy = 0;
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            for (int z = 0; z < L; z++) {
                int atom = lattice[x][y][z];
                int neighbors[6][3] = {
                    {(x + 1) % L, y, z},
                    {(x - 1 + L) % L, y, z},
                    {x, (y + 1) % L, z},
                    {x, (y - 1 + L) % L, z},
                    {x, y, (z + 1) % L},
                    {x, y, (z - 1 + L) % L}
                };
                for (int i = 0; i < 6; i++) {
                    int nx = neighbors[i][0];
                    int ny = neighbors[i][1];
                    int nz = neighbors[i][2];
                    int neighbor_atom = lattice[nx][ny][nz];
                    if (atom == neighbor_atom) {
                        energy += (atom == 0) ? E_AA : E_BB;
                    } else {
                        energy += E_AB;
                    }
                }
            }
        }
    }
    return energy / 2; // Each bond is counted twice
}

// Function to perform the Metropolis algorithm
void metropolis(int lattice[L][L][L], double T, int steps) {
    double energy = calculate_energy(lattice);
    for (int step = 0; step < steps; step++) {
        int x1 = rand() % L;
        int y1 = rand() % L;
        int z1 = rand() % L;
        int x2 = rand() % L;
        int y2 = rand() % L;
        int z2 = rand() % L;

        int temp = lattice[x1][y1][z1];
        lattice[x1][y1][z1] = lattice[x2][y2][z2];
        lattice[x2][y2][z2] = temp;

        double new_energy = calculate_energy(lattice);
        double dE = new_energy - energy;

        if (dE < 0 || ((double)rand() / RAND_MAX) < exp(-dE / (k_B * T))) {
            energy = new_energy;
        } else {
            temp = lattice[x1][y1][z1];
            lattice[x1][y1][z1] = lattice[x2][y2][z2];
            lattice[x2][y2][z2] = temp;
        }
    }
}

int run(int argc, char *argv[]) {
    int lattice[L][L][L];
    double T_list[] = {400, 600, 1000};

    for (int i = 0; i < 3; i++) {
        double T = T_list[i];
        initialize_lattice(lattice, 1);
        metropolis(lattice, T, 1000000);
        double energy = calculate_energy(lattice);
        printf("Average energy at T = %.1f K: %.2f meV\n", T, energy);
    }

    return 0;
}
