#include <stdio.h>
#include <math.h>
#include <stdlib.h>


int
run(int argc, char *argv[]){
    double T, 
    m, 
    U, 
    C, 
    U_prev = 0.0;
    FILE *output = fopen("results_mean_field.csv", "w");
    if (!output) {
        perror("Failed to open file");
        return 1;
    }

    fprintf(output, "T,m,U,C\n");

    for (T = T_START; T <= T_END; T += T_STEP) {
        // Find magnetization m(T)
        m = find_m(T);

        // Calculate U(T)
        U = calculate_energy(m);

        // Calculate C(T) using finite difference
        if (T > T_START) {
            C = calculate_heat_capacity(U_prev, U, T_STEP);
        } else {
            C = 0.0; // Heat capacity undefined at T=0
        }

        // Save results
        fprintf(output, "%.5f,%.5f,%.5f,%.5f\n", T, m, U, C);

        // Update previous U
        U_prev = U;
    }

    fclose(output);
    printf("Results written to results_mean_field.csv\n");
    return 0;
}
