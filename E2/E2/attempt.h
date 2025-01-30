#ifndef ATTEMPT_H
#define ATTEMPT_H

#include <stdlib.h>

// Function to dynamically create a 2D array
// Parameters:
//   n - Number of rows
//   m - Number of columns
// Returns:
//   Pointer to a dynamically allocated 2D array
double **create_2D_array(unsigned int n, unsigned int m);

// Function to destroy a dynamically allocated 2D array
// Parameters:
//   array - Pointer to the 2D array
//   n - Number of rows
void destroy_2D_array(double **array, unsigned int n);

// Function for matrix-vector multiplication
// Parameters:
//   result - Vector to store the result
//   A - 2D matrix
//   b - Vector to be multiplied
//   n - Number of rows in the matrix
//   m - Number of columns in the matrix (and size of b)
void matrix_vector_multiplication(double *result, double **A, double *b, unsigned int n, unsigned int m);

// Function to transform to normal modes
// Parameters:
//   positions - Vector of positions
//   Q - Vector to store Q-coordinates
//   N - Number of particles
void transform_to_normal_modes(double *positions, double *Q, const unsigned int N);

// Function to calculate the normal mode energies
// Parameters:
//   energies - Vector to store the energies
//   positions - Vector of positions
//   velocities - Vector of velocities
//   N - Number of particles
void calculate_normal_mode_energies(double *energies, double *positions, double *velocities, const unsigned int N);

#endif // ATTEMPT_H
