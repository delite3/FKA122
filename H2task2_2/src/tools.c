#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <complex.h>

#include "/Users/darre/Desktop/Chalmers/vscode/H1/include/tools.h"


void transform_to_normal_modes(double *positions, double *Q, const unsigned int N)
{
    // Create transformation matrix A (N x N)
    double **A = create_2D_array(N, N);

    // Fill the transformation matrix A
    for (unsigned int k = 0; k < N; k++) {
        for (unsigned int i = 0; i < N; i++) {
            A[k][i] = sqrt(2.0 / (N + 1)) * sin((i + 1) * (k + 1) * M_PI / (N + 1));
        }
    }

    // Multiply matrix A with positions vector
    matrix_vector_multiplication(Q, A, positions, N, N);

    // Clean up
    destroy_2D_array(A, N);
}


void calculate_normal_mode_energies(double *energies, double *positions, double *velocities, const unsigned int N)
{
    double *Q = (double *)calloc(N, sizeof(double));  // Normal mode displacements
    double *P = (double *)calloc(N, sizeof(double));  // Normal mode momenta

    transform_to_normal_modes(positions, Q, N);
    transform_to_normal_modes(velocities, P, N);

    // Compute energies for each mode
    for (unsigned int k = 0; k < N; k++) {
        double omega_k = 2.0 * sin((k + 1) * M_PI / (2.0 * (N + 1))); // Frequency in dimensionless units
        energies[k] = 0.5 * (P[k] * P[k] + omega_k * omega_k * Q[k] * Q[k]);
    }

    // Clean up
    free(P);
    free(Q);
}


void calculate_acceleration(double *accelerations, double *positions, double alpha, const unsigned int N)
{
    // Write to the accelerations vector
for (unsigned int i = 0; i < N; i++)
    {
        double term1 = 0.0, term2 = 0.0;

        if (i == 0) {
            term1 = positions[i + 1] - 2 * positions[i];  // Harmonic term
            term2 = alpha * ((positions[i + 1] - positions[i]) * (positions[i + 1] - positions[i])
                            - (positions[i] * positions[i]));  // Anharmonic term
        } else if (i == N - 1) {
            term1 = positions[i - 1] - 2 * positions[i];  // Harmonic term
            term2 = alpha * ((-positions[i]) * (-positions[i])
                            - (positions[i] - positions[i - 1]) * (positions[i] - positions[i - 1]));  // Anharmonic term
        } else {
            term1 = positions[i + 1] - 2 * positions[i] + positions[i - 1];  // Harmonic term
            term2 = alpha * ((positions[i + 1] - positions[i]) * (positions[i + 1] - positions[i])
                            - (positions[i] - positions[i - 1]) * (positions[i] - positions[i - 1]));  // Anharmonic term
        }

        accelerations[i] = term1 + term2;  // Update acceleration
    }
}

double 
calculate_potential_energy(double *positions, double kappa)
{

    double potential_energy = 0.0;

    // Potential energy for the first and second springs
    potential_energy += 0.5 * kappa * (positions[1] - positions[0]) * (positions[1] - positions[0]);
    potential_energy += 0.5 * kappa * (positions[2] - positions[1]) * (positions[2] - positions[1]);

    return potential_energy;
}

double 
calculate_kinetic_energy(double *velocities, double *masses)
{
    const unsigned int N = 3;
    double kinetic_energy = 0.0;

    for (unsigned int i = 0; i < N; i++)
    {
        kinetic_energy += 0.5 * masses[i] * velocities[i] * velocities[i];
    }

    return kinetic_energy;
}

void velocity_verlet_one_step(double *accelerations, double *positions, double *velocities,
                              double alpha, double timestep, const unsigned int N)
{
    // Step 1: Half-step velocity update
    for (unsigned int i = 0; i < N; i++)
    {
        velocities[i] += 0.5 * accelerations[i] * timestep;
    }

    // Step 2: Full-step position update
    for (unsigned int i = 0; i < N; i++)
    {
        positions[i] += velocities[i] * timestep;
    }

    // Step 3: Recalculate accelerations
    calculate_acceleration(accelerations, positions, alpha, N);

    // Step 4: Final half-step velocity update
    for (unsigned int i = 0; i < N; i++)
    {
        velocities[i] += 0.5 * accelerations[i] * timestep;
    }
}


void
elementwise_addition(
    double *res,     // Output vector to store the results
    double *v1,         // Input vector 1
    double *v2,         // Input vector 2
    unsigned int len   // Number of elements in each vector
)
{
    for (unsigned int i = 0; i < len; i++) {
        res[i] = v1[i] + v2[i];
    }
}


void
elementwise_multiplication(
                           double *res,
                           double *v1,
                           double *v2,
                           unsigned int len
                          )
{
    for (unsigned int i = 0; i < len; i++) {
        res[i] = v1[i] * v2[i];
    }
}

void
addition_with_constant(
                       double *res,
                       double *v,
                       double constant,
                       unsigned int len)
{
    for (unsigned int i = 0; i < len; i++) {
        res[i] = v[i] + constant;
    }
}

void
multiplication_with_constant(
                             double *res,
                             double *v,
                             double constant,
                             unsigned int len)
{
    for (unsigned int i = 0; i < len; i++) {
        res[i] = v[i] * constant;
    }
}

double
dot_product(
            double *v1,
            double *v2,
            unsigned int len
           )
{
    double result = 0.0;
    for (unsigned int i = 0; i < len; i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

double **
create_2D_array(
                unsigned int row_size,
                unsigned int column_size
               )
{
    if (row_size == 0 || column_size == 0) {
        fprintf(stderr, "Invalid dimensions for 2D array: rows=%u, columns=%u\n", row_size, column_size);
        return NULL;
    }

    // Allocate memory for row pointers
    double **array = calloc(row_size,  sizeof(double *));
    if (array == NULL) {
        perror("Failed to allocate memory for row pointers.");
        return NULL;
    }

    // Allocate a single block of memory for all rows
    double *data_block = calloc(row_size * column_size, sizeof(double));
    if (data_block == NULL) {
        perror("Failed to allocate memory for data block.");
        free(array);
        return NULL;
    }

    // Assign row pointers to appropriate offsets in the data block
    for (unsigned int i = 0; i < row_size; ++i) {
        array[i] = &data_block[i * column_size];
    }

    return array;
}


void
destroy_2D_array(
                 double **array,
                 unsigned int row_size
                )
{
    if (array == NULL) {
        return;
    }

    // Free the single block of memory holding all rows
    free(array[0]);  // Free the contiguous data block

    // Free the array of row pointers
    free(array);
}


void
matrix_vector_multiplication(
    double *result,      // Output vector (n x 1)
    double **A,          // Input matrix (n x m)
    double *B,           // Input vector (m x 1)
    unsigned int n,      // Number of rows in A
    unsigned int m       // Number of columns in A (and size of B)
)
{
    // Iterate through each row of the matrix
    for (unsigned int i = 0; i < n; ++i) {
        result[i] = 0; // Initialize result for the i-th row
        for (unsigned int j = 0; j < m; ++j) {
            // Multiply and accumulate: result[i] += A[i][j] * B[j]
            result[i] += A[i][j] * B[j];
        }
    }
}

void
matrix_matrix_multiplication(
                             double **result,
                             double **A,
                             double **B,
                             unsigned int n,
                             unsigned int m,
                             unsigned int k
                            )
{
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < k; j++) {
            result[i][j] = 0.0;
            for (unsigned int l = 0; l < m; l++) {
                result[i][j] += A[i][l] * B[l][j];
            }
        }
    }
}

double
vector_norm(
            double *v1,
            unsigned int len
           )
{
    double sum = 0.0;
    for (unsigned int i = 0; i < len; i++) {
        sum += v1[i] * v1[i];
    }
    return sqrt(sum);
}


void
normalize_vector(
                 double *v1,
                 unsigned int len
                )
{
    double norm = vector_norm(v1, len);
    if (norm == 0.0) return;  // Avoid division by zero
    for (unsigned int i = 0; i < len; i++) {
        v1[i] /= norm;
    }
}

double
average(
        double *v1,
        unsigned int len
       )
{
    if (len == 0) return 0;
    double sum = 0.0;
    for (unsigned int i = 0; i < len; i++) {
        sum += v1[i];
    }
    return sum / len;
}


double
standard_deviation(
                       double *v1,
                       unsigned int len
                  )
{
    if (len == 0) return 0;
    double avg = average(v1, len);
    double sum_sq_diff = 0.0;
    for (unsigned int i = 0; i < len; i++) {
        sum_sq_diff += (v1[i] - avg) * (v1[i] - avg);
    }
    return sqrt(sum_sq_diff / len);
}

double
distance_between_vectors(
                         double *v1,
                         double *v2,
                         unsigned int len
                        )
{
    double sum = 0.0;
    for (unsigned int i = 0; i < len; i++) {
        double diff = v1[i] - v2[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

void
cumulative_integration(
                       double *res,
                       double *v,
                       double dx,
                       unsigned int v_len
                      )
{
    res[0] = 0.0;
    for (unsigned int i = 1; i < v_len; i++) {
        res[i] = res[i - 1] + (v[i - 1] + v[i]) * dx / 2.0;
    }
}


void
write_xyz(
          FILE *fp,
          char *symbol,
          double **positions,
          double **velocities,
          double alat,
          int natoms)
{
    fprintf(fp, "%i\nLattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\" ", natoms, alat, alat, alat);
    fprintf(fp, "Properties=species:S:1:pos:R:3:vel:R:3 pbc=\"T T T\"\n");
    for (int i = 0; i < natoms; ++i) {
        fprintf(fp, "%s %f %f %f %f %f %f\n",
                symbol,
                positions[i][0], positions[i][1], positions[i][2],
                velocities[i][0], velocities[i][1], velocities[i][2]);
    }
}

void fft_freq(
              double *res,
              int n,
              double timestep)
{
    // Total frequency range
    double freq_range = 1.0 / (n * timestep);

    
    //printf("int n: %i\n", n);
    //printf("timestep: %f\n", timestep);

    //printf("freq_range: %f\n", freq_range);

    for (int k = 0; k < n; ++k) {
        if (k < n / 2) {
            // Positive frequencies
            res[k] = 2 * M_PI * k * freq_range; // Multiply by 2π for angular frequency
        } else {
            // Negative frequencies
            res[k] = 2 * M_PI * (k - n) * freq_range; // Multiply by 2π for angular frequency
        }
    
        //printf("k: %i, res[k]: %f\n", k, res[k]);
    }
}

/* Freely given functions */
void
skip_line(FILE *fp)
{
    int c;
    while (c = fgetc(fp), c != '\n' && c != EOF);
}

void 
read_xyz(
    FILE *fp,
    char *symbol,
    double **positions,
    double **velocities,
    double *alat)
{
    int natoms;
    if (fscanf(fp, "%i\nLattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 %lf\" ", &natoms, alat, alat, alat) != 4) {
        fprintf(stderr, "Error reading header line.\n");
        exit(1);
    }

    skip_line(fp);

    for (int i = 0; i < natoms; ++i) {
        if (fscanf(fp, "%s %lf %lf %lf ",
                   symbol, &positions[i][0], &positions[i][1], &positions[i][2]) != 4) {
            fprintf(stderr, "Error reading positions for atom %d.\n", i);
            exit(1);
        }

        if (fscanf(fp, "%lf %lf %lf\n",
                   &velocities[i][0], &velocities[i][1], &velocities[i][2]) != 3) {
            fprintf(stderr, "Error reading velocities for atom %d.\n", i);
            exit(1);
        }
    }
}

void powerspectrum(
           double *res,
           double *signal,
           int n,
                   double timestep)
{
    /* Declaration of variables */
    double *complex_coefficient = malloc(sizeof(double) * 2*n); // array for the complex fft data
    double *data_cp = malloc(sizeof(double) * n);

    /*make copy of data to avoid messing with data in the transform*/
    for (int i = 0; i < n; i++) {
    data_cp[i] = signal[i];
    }

    /* Declare wavetable and workspace for fft */
    gsl_fft_real_wavetable *real;
    gsl_fft_real_workspace *work;

    /* Allocate space for wavetable and workspace for fft */
    work = gsl_fft_real_workspace_alloc(n);
    real = gsl_fft_real_wavetable_alloc(n);

    /* Do the fft*/
    gsl_fft_real_transform(data_cp, 1, n, real, work);

    /* Unpack the output into array with alternating real and imaginary part */
    gsl_fft_halfcomplex_unpack(data_cp, complex_coefficient,1,n);

    /*fill the output powspec_data with the powerspectrum */
    for (int i = 0; i < n; i++) {
    res[i] = (complex_coefficient[2*i]*complex_coefficient[2*i]+complex_coefficient[2*i+1]*complex_coefficient[2*i+1]);
    res[i] *= timestep / n;
    }

    /* Free memory of wavetable and workspace */
    gsl_fft_real_wavetable_free(real);
    gsl_fft_real_workspace_free(work);
    free(complex_coefficient);
    free(data_cp);
}
