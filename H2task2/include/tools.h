#pragma once
#include <stdio.h>

/**
 * Transform to normal modes.
 *
 * Parameters
 * ----------
 *  positions - Vector of positions
 *  Q - Vector of Q-coordinates
 *  N - Number of particles (length of vectors)
 *
*/
void 
transform_to_normal_modes(
                          double *positions, 
                          double *Q, 
                          const unsigned int N
                          );

/**
 * Calculate the normal mode energies
 *
 * Parameters
 * ----------
 *  energies - Vector where energies will be written
 *  positions - Vector of positions
 *  velocities - Vector of positions
 *  N - Number of particles (length of vectors)
 *
*/
void 
calculate_normal_mode_energies(
                               double *energies, 
                               double *positions, 
                               double *velocities, 
                               const unsigned int N
                               );

/**
 * Calculate the acceleration.
*
* Parameters
* ----------
*  accelerations - Vector where accelerations are to be written
*  positions - Vector of positions
*  alpha - Anharmonicity constant
*  N - Number of atoms
*
*/
void
calculate_acceleration(
                       double *accelerations, 
                       double *positions, 
                       double alpha, 
                       const unsigned int N
                       );


/**
* Calculate the potential energy. All vectors are assumed to be of length 3
*
* Parameters
* ----------
* positions - Vector of positions
* kappa - Spring constant
* Returns
* -------
* Potential energy
*
*/

double 
calculate_potential_energy(
                           double *positions,
                           double kappa
                           );


/**
* Calculate the kinetic energy. All vectors are assumed to be of length 3
*
* Parameters
* ----------
* velocities - Vector of velocities
* masses - Vector of masses
*
* Returns
* -------
* Kinetic energy
*
*/
double 
calculate_kinetic_energy(
                         double *velocities, 
                         double *masses
                         );

/**
 * Perform one velocity Verlet step
 *
 * Parameters
 * ----------
 *  accelerations - Vector of accelerations
 *  positions - Vector of positions
 *  velocities - Vector of velocities
 *  alpha - Anharmonicity constant
 *  timestep - Time step
 *  N - Number of atoms
 *
*/
void 
velocity_verlet_one_step(
                         double *accelerations, 
                         double *positions, 
                         double *velocities, 
                         double alpha, 
                         double timestep, 
                         const unsigned int N
                         );


/* **********************************************
 *
 * Add v1 and v2 elementwise
 * results is stored in res.
 * The length of the arrays should be len
 *
 * **********************************************/

void
elementwise_addition(
                     double *res,
                     double *v1,
                     double *v2,
                     unsigned int len
                    );

/* **********************************************
 *
 * Multiply v1 and v2 elementwise
 * results is stored in res.
 * The length of the arrays should be len
 *
 * **********************************************/
void
elementwise_multiplication(
                           double *res,
                           double *v1,
                           double *v2,
                           unsigned int len
                          );

/* **********************************************
 *
 * Add const with every element in v
 * results is stored in res.
 * The length of the arrays should be len
 *
 * **********************************************/
void
addition_with_constant(
                       double *res,
                       double *v,
                       double constant,
                       unsigned int len);

/* **********************************************
 *
 * Multiply const with every element in v
 * results is stored in res.
 * The length of the arrays should be len
 *
 * **********************************************/
void
multiplication_with_constant(
                             double *res,
                             double *v,
                             double constant,
                             unsigned int len);

/* **********************************************
 *
 * Calculate the dot product between
 * v1 and v2 with length len
 *
 * the result is returned as a double
 *
 * **********************************************/
double
dot_product(
            double *v1,
            double *v2,
            unsigned int len
           );


/* **********************************************
 *
 * Allocate the memory to a 2D array
 *
 * Array is the varible in which in the array is
 * stored.
 * The matrix is of size n x m.
 *
 * **********************************************/
double **
create_2D_array(
                unsigned int n,
                unsigned int m
               );

/* **********************************************
 *
 * Deallocate the memory of a 2D array
 *
 * The input n is option and is not necessary.
 * Try to find a solution in create_2D_array
 * that allows you to not use n.
 *
 * **********************************************/
void
destroy_2D_array(
                 double **array,
                 unsigned int n
                );


/* **********************************************
 *
 * Calculates the matrix multiplication:
 *
 * y = A x b.
 *
 * with dimensions
 *
 * n x 1 ::: n x m ::: m x 1
 *
 * The n x 1 matrices y and b are represented as plain vectors
 *
 * **********************************************/
void
matrix_vector_multiplication(
                             double *y,
                             double **A,
                             double *b,
                             unsigned int n,
                             unsigned int m
                            );

/* **********************************************
 *
 * Calculates the matrix multiplication:
 *
 * Y = A x B
 *
 * with dimensions
 *
 * n x k ::: n x m ::: m x k
 *
 * **********************************************/
void
matrix_matrix_multiplication(
                             double **Y,
                             double **A,
                             double **B,
                             unsigned int n,
                             unsigned int m,
                             unsigned int k
                            );

/* **********************************************
 *
 * Calculates the L2 norm of v1 with length
 * len
 *
 * **********************************************/
double
vector_norm(
            double *v1,
            unsigned int len
           );

/* **********************************************
 *
 * Normalizes the vector v1 such that the L2
 * norm of v1 is 1.
 *
 * The length of the vector is len
 *
 * **********************************************/
void
normalize_vector(
                 double *v1,
                 unsigned int len
                );

/* **********************************************
 *
 * Calculates the average of a vector with length
 * len
 *
 * **********************************************/
double
average(
        double *v1,
        unsigned int len
       );

/* **********************************************
 *
 * Calculates the standard deviation of a vector
 * with length len
 *
 * **********************************************/
double
standard_deviation(
                       double *v1,
                       unsigned int len
                  );

/* **********************************************
 *
 * Calculates the Euclidean distance  between
 * two points v1 and v2 with dimension len.
 *
 * **********************************************/
double
distance_between_vectors(
                         double *v1,
                         double *v2,
                         unsigned int len
                        );

/* **********************************************
 *
 * Calculates the cumulative integration of the
 * array v with length len using the
 * trapezoidal rule.
 * result is stored in res which has the same
 * length as v and dx is the step length.
 *
 * **********************************************/
void
cumulative_integration(
                       double *res,
                       double *v,
                       double dx,
                       unsigned int len
                      );


/* **********************************************
 *
 * Writes atoms to extxyz format which can be used
 * to visualize atoms in ase for example. fp is
 * the file-pointer to the file, symbol is the
 * atomic symbol, e.g., "Al", positions are the
 * positions of the atoms, velocities are the
 * velocities of the atoms, alat is the cell
 * length and natoms is the number of atoms.
 *
 * **********************************************/
void
write_xyz(
          FILE *fp,
          char *symbol,
          double **positions,
          double **velocities,
          double alat,
          int natoms);


/* **********************************************
 *
 * Performs a Fourier transform of data.
 * The result is stored in res and signal is the
 * input to be Fourier transformed,
 * n is the length of the signal and timestep is
 * the difference in time between two succesive
 * signal values.
 *
 * **********************************************/
void powerspectrum(
               double *res,
               double *signal,
               int n,
               double timestep);


/* **********************************************
 *
 * Calculates the frequencies for the
 * powerspectrum data
 *
 * The ordering is the same as for numpy, check
 * the numpy documentation/implementation. However,
 * it should be given in angular frequencies, i.e.,
 * multiply by 2*PI
 *
 * check powerspectrum for arguments.
 *
 * **********************************************/
void fft_freq(
          double *res,
          int n,
          double timestep);


/* **********************************************
 *
 * Reads atoms from extxyz format.
 *
 * **********************************************/
void
read_xyz(
          FILE *fp,
          char *symbol,
          double **positions,
          double **velocities,
          double *alat);
