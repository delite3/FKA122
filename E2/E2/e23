#include <math.h>
#include "attempt.h"

/*
 * The following functions are to be used:
 */
/*
double
average(
        double *v,
        unsigned int len
       );
double
variance(
         double *v,
         unsigned int len
        );
*/

/* ***************************************
 *
 * Autocorrelation
 *
 * Parameters
 * ----------
 *  data - the raw data from which the
 *         autocorrelation should be
 *         calculated
 *  data_len     - the len of data
 *  time_lag_ind - the lag (index) at 
 *                 which the autocorrelation
 *                 should be calculated
 *
 * Returns
 * -------
 *  The autocorrelation at a specific time
 *  lag
 *
 * ***************************************/
double autocorrelation(
    double *data,
    int data_len,
    int time_lag_ind)
{
    double mean = average(data, data_len); // Calculate the mean of the data
    double variance_val = variance(data, data_len); // Calculate the variance
    if (variance_val == 0) {
        return 0; // Avoid division by zero
    }

    // Compute the autocorrelation
    double numerator = 0;
    for (int i = 0; i < data_len - time_lag_ind; i++) {
        numerator += (data[i] - mean) * (data[i + time_lag_ind] - mean);
    }
    numerator /= (data_len - time_lag_ind);

    return numerator / variance_val; // Autocorrelation at time lag
}

/* ***************************************
 *
 * Block average
 *
 * Parameters
 * ----------
 *  data       - the raw data from which the
 *               autocorrelation should be
 *               calculated
 *  data_len   - the length of data
 *  block_size - the size of the block
 *
 * Returns
 * -------
 * The statistical inefficiency for a given
 * block size
 *
 * ***************************************/
double block_average(double *data,
                     int data_len,
                     int block_size)
{
    if (block_size > data_len || block_size <= 0) {
        return 0; // Invalid block size
    }

    int num_blocks = data_len / block_size; // Number of blocks
    double *block_means = (double *)malloc(num_blocks * sizeof(double));
    if (!block_means) {
        return 0; // Memory allocation failed
    }

    // Compute the mean for each block
    for (int i = 0; i < num_blocks; i++) {
        block_means[i] = average(&data[i * block_size], block_size);
    }

    // Compute the variance of block means
    double block_variance = variance(block_means, num_blocks);

    // Free allocated memory
    free(block_means);

    if (block_variance == 0) {
        return 0; // Avoid division by zero
    }

    // Calculate statistical inefficiency
    double data_variance = variance(data, data_len);
    return (block_size * block_variance) / data_variance;
}
