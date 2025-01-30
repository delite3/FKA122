#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 
 * Function prototypes
 */
double average(double *v, unsigned int len);
double variance(double *v, unsigned int len);
double autocorrelation(double *data, int data_len, int time_lag);
double find_s_from_autocorr(double *data, int data_len);
double block_average(double *data, int data_len, int block_size);
double find_s_from_block_averaging(double *data, int data_len, double s_auto_hint);

int main(void)
{
    // Number of data points
    int N = 1000000;

    // Allocate memory for data
    double *data = (double *)malloc(N * sizeof(double));
    if (!data) {
        perror("Memory allocation failed");
        return EXIT_FAILURE;
    }

    // Read data from MC.txt
    FILE *fp = fopen("MC.txt", "r");
    if (!fp) {
        perror("Error opening MC.txt");
        free(data);
        return EXIT_FAILURE;
    }

    for (int i = 0; i < N; i++) {
        if (fscanf(fp, "%lf", &data[i]) != 1) {
            fprintf(stderr, "Error reading data at line %d\n", i + 1);
            fclose(fp);
            free(data);
            return EXIT_FAILURE;
        }
    }
    fclose(fp);

    // Compute mean and variance of original data
    double data_mean = average(data, N);
    double data_var = variance(data, N);

    printf("Data length: %d\n", N);
    printf("Mean of data: %lf\n", data_mean);
    printf("Variance of data: %lf\n", data_var);

    // Center the data
    double *data_centered = (double *)malloc(N * sizeof(double));
    if (!data_centered) {
        perror("Memory allocation failed");
        free(data);
        return EXIT_FAILURE;
    }

    for (int i = 0; i < N; i++) {
        data_centered[i] = data[i] - data_mean;
    }

    // Estimate s from autocorrelation on centered data
    double s_auto = find_s_from_autocorr(data_centered, N);
    printf("Estimated s from autocorrelation (centered data): %g\n", s_auto);

    // Estimate s from block averaging on centered data
    // We pass s_auto as a hint for scaling how far we need to go in block sizes
    double s_block = find_s_from_block_averaging(data_centered, N, s_auto);
    printf("Estimated s from block averaging (centered data): %g\n", s_block);

    free(data_centered);
    free(data);

    return 0;
}

/* Compute the average of a vector */
double average(double *v, unsigned int len) {
    double sum = 0.0;
    for (unsigned int i = 0; i < len; i++) {
        sum += v[i];
    }
    return sum / (double)len;
}

/* Compute the variance of a vector */
double variance(double *v, unsigned int len) {
    double mu = average(v, len);
    double var_sum = 0.0;
    for (unsigned int i = 0; i < len; i++) {
        double diff = v[i] - mu;
        var_sum += diff * diff;
    }
    return var_sum / (double)len;
}

/* Compute autocorrelation at a given time lag for centered data */
double autocorrelation(double *data, int data_len, int time_lag) {
    double var = variance(data, data_len);
    if (var == 0.0) {
        return 0.0;
    }

    double sum = 0.0;
    int count = data_len - time_lag;
    for (int i = 0; i < count; i++) {
        sum += data[i] * data[i + time_lag];
    }
    sum /= count;
    return sum / var;
}

/*
 * Find s from autocorrelation by searching for lag k where Phi_k ~ exp(-2).
 * We start at k=1 and go up to a large fraction of data_len if needed.
 */
double find_s_from_autocorr(double *data, int data_len) {
    double target = exp(-2.0); // ~0.135
    // We'll search up to data_len/100 or 100000 lags, whichever is smaller, 
    // to avoid huge computation times. You can adjust this.
    int max_lag = data_len/100;
    if (max_lag > 100000) {
        max_lag = 100000;
    }

    for (int k = 1; k < max_lag; k++) {
        double phi_k = autocorrelation(data, data_len, k);
        if (phi_k <= target) {
            return (double) k;
        }
    }
    // If not found, return the largest searched. You can handle this differently.
    return (double)max_lag;
}

/* Compute statistical inefficiency for a given block size B */
double block_average(double *data, int data_len, int block_size) {
    if (block_size <= 0 || block_size > data_len) {
        return NAN;
    }

    int num_blocks = data_len / block_size;
    double *block_means = (double *)malloc(num_blocks * sizeof(double));
    if (!block_means) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < num_blocks; i++) {
        block_means[i] = average(&data[i * block_size], block_size);
    }

    double block_var = variance(block_means, num_blocks);
    free(block_means);

    double data_var = variance(data, data_len);
    if (data_var == 0.0) {
        return NAN;
    }

    return (block_size * block_var) / data_var;
}

/*
 * Find s from block averaging by increasing block size until s stabilizes.
 * We use s_auto_hint to ensure we test beyond that range.
 * 
 * Strategy:
 * - Double the block size starting from 1 until we reach something like data_len/100.
 * - Keep track of previous s values.
 * - If the relative difference is less than 1% over a few steps, consider it stable.
 */
double find_s_from_block_averaging(double *data, int data_len, double s_auto_hint) {
    double prev_s = -1.0;
    double tolerance = 0.01; // 1% tolerance for stability
    double final_s = -1.0;

    // We will go up to data_len/100 block size, which for N=1,000,000 gives us up to block_size=10,000.
    // This is arbitrary and can be increased if needed. If s is very large, you may need more.
    int max_block_size = data_len / 100;
    if (max_block_size < 1) max_block_size = 1;

    // Start from block_size=1 and double each time until max_block_size is reached
    // or until stabilization is detected.
    for (int B = 1; B <= max_block_size; B *= 2) {
        double s_val = block_average(data, data_len, B);
        if (!isnan(s_val) && s_val > 0.0) {
            // Print diagnostics if you want
            // printf("Block size %d: s = %g\n", B, s_val);

            if (prev_s > 0.0) {
                double rel_diff = fabs(s_val - prev_s) / prev_s;
                if (rel_diff < tolerance && B > s_auto_hint) {
                    // If we are beyond the autocorrelation estimate and stable, stop.
                    final_s = s_val;
                    break;
                }
            }
            prev_s = s_val;
        }
    }

    if (final_s < 0.0) {
        // If we didn't break out due to stabilization, just return the last computed s.
        final_s = prev_s;
    }

    return final_s;
}
