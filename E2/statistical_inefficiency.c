#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Function prototypes */
double average(double *v, unsigned int len);
double variance(double *v, unsigned int len);
double autocorrelation(double *data, int data_len, int time_lag);
double find_s_from_autocorr(double *data, int data_len);
double block_average(double *data, int data_len, int block_size);
double find_s_from_block_averaging(double *data, int data_len);
double moving_average(double *values, int len, int window_size, int idx);

int main(void) {
    int N = 1000000; // Number of data points

    // Allocate memory for the dataset
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

    // Calculate mean and variance of the original data
    double data_mean = average(data, N);
    double data_var = variance(data, N);

    printf("Data length: %d\n", N);
    printf("Mean of data: %.10f\n", data_mean);
    printf("Variance of data: %.10f\n", data_var);

    // Center the data by subtracting the mean
    double *data_centered = (double *)malloc(N * sizeof(double));
    if (!data_centered) {
        perror("Memory allocation failed");
        free(data);
        return EXIT_FAILURE;
    }

    for (int i = 0; i < N; i++) {
        data_centered[i] = data[i] - data_mean;
    }

    // Estimate s using the autocorrelation method
    double s_auto = find_s_from_autocorr(data_centered, N);
    printf("Estimated s from autocorrelation (centered data): %.2f\n", s_auto);

    // Estimate s using the block averaging method
    double s_block = find_s_from_block_averaging(data_centered, N);
    printf("Estimated s from block averaging (centered data): %.2f\n", s_block);

    // Free allocated memory
    free(data_centered);
    free(data);

    return 0;
}

/* Compute the average of an array */
double average(double *v, unsigned int len) {
    double sum = 0.0;
    for (unsigned int i = 0; i < len; i++) {
        sum += v[i];
    }
    return sum / len;
}

/* Compute the variance of an array */
double variance(double *v, unsigned int len) {
    double avg = average(v, len);
    double var_sum = 0.0;
    for (unsigned int i = 0; i < len; i++) {
        double diff = v[i] - avg;
        var_sum += diff * diff;
    }
    return var_sum / len;
}

/* Compute the autocorrelation at a given time lag for centered data */
double autocorrelation(double *data, int data_len, int time_lag) {
    double var = variance(data, data_len);
    if (var == 0.0) {
        return 0.0;
    }

    double sum = 0.0;
    for (int i = 0; i < data_len - time_lag; i++) {
        sum += data[i] * data[i + time_lag];
    }
    sum /= (data_len - time_lag);
    return sum / var;
}

/* Find s using the autocorrelation method with extended lag range */
double find_s_from_autocorr(double *data, int data_len) {
    double target = exp(-2.0); // Threshold for statistical inefficiency
    double negligible_threshold = 1e-5; // Early stopping for negligible correlation
    int max_lag = data_len / 20; // Extend the lag range

    double s = 0.5; // Start with 0.5 (for lag 0)

    for (int lag = 1; lag < max_lag; lag++) {
        double phi_k = autocorrelation(data, data_len, lag);
        if (phi_k < target || phi_k < negligible_threshold) {
            break;
        }
        s += 2.0 * phi_k; // Sum contributions for each lag
    }
    return s;
}

/* Compute statistical inefficiency for a given block size */
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

/* Smooth values using a simple moving average */
double moving_average(double *values, int len, int window_size, int idx) {
    double sum = 0.0;
    int count = 0;
    for (int i = idx - window_size / 2; i <= idx + window_size / 2; i++) {
        if (i >= 0 && i < len) {
            sum += values[i];
            count++;
        }
    }
    return count > 0 ? sum / count : values[idx];
}

/* Find s from block averaging with smoother stabilization */
double find_s_from_block_averaging(double *data, int data_len) {
    int max_block_size = data_len / 10; // Larger max block size for better results
    double tolerance = 0.005; // 0.5% tolerance for stabilization
    int step = 40; // Smaller step size for smoother results
    int stabilization_window = 20; // Larger stabilization window

    double *s_values = (double *)malloc(max_block_size / step * sizeof(double));
    if (!s_values) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    int idx = 0;
    double prev_s = -1.0;
    double final_s = -1.0;
    int stabilization_count = 0;

    for (int block_size = step; block_size <= max_block_size; block_size += step) {
        double s_val = block_average(data, data_len, block_size);
        s_values[idx++] = s_val;

        if (!isnan(s_val) && s_val > 0.0) {
            printf("Block size %d: s = %.5f\n", block_size, s_val);

            // Apply a moving average to smooth results
            double smoothed_s = moving_average(s_values, idx, 5, idx - 1);

            if (prev_s > 0.0) {
                double rel_diff = fabs(smoothed_s - prev_s) / prev_s;
                if (rel_diff < tolerance) {
                    stabilization_count++;
                    if (stabilization_count >= stabilization_window) {
                        final_s = smoothed_s;
                        break;
                    }
                } else {
                    stabilization_count = 0; // Reset stabilization count
                }
            }
            prev_s = smoothed_s;
        }
    }

    free(s_values);
    return final_s > 0.0 ? final_s : prev_s;
}
