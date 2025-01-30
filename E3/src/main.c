#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <string.h>

#define PI 3.141592653589

// Struct to store the result
typedef struct {
    double position;
    double velocity; 
} result_t;

result_t BD3(double initial_position, double initial_velocity, double w0, double dt, double eta, double kB, double mass, double T, gsl_rng *k)
{
    result_t result;
    
    // Constants
    double c0 = exp(-eta * dt);     // Damping factor
    double vth = sqrt(kB * T / mass); // Thermal velocity
    
    // Initialize position and velocity
    double x = initial_position;
    double v = initial_velocity;
    
    // Calculate the random Gaussian numbers G1 and G2
    double G1 = gsl_ran_gaussian(k, 1.0); // G1 ~ N(0,1)
    double G2 = gsl_ran_gaussian(k, 1.0); // G2 ~ N(0,1)
    
    // Compute the acceleration due to external force: F_ext = -m * w0^2 * x
    double a = -w0 * w0 * x;
    
    // Step 1: Update velocity (half-step)
    double v_half = 0.5 * a * dt + sqrt(c0) * v + vth * sqrt(1.0 - c0) * G1;
    
    // Step 2: Update position
    x = x + v_half * dt;
    
    // Calculate new acceleration due to external force
    a = -w0 * w0 * x;
    
    // Step 3: Update velocity (full-step)
    v = 0.5 * sqrt(c0) * a * dt + sqrt(c0) * v_half + vth * sqrt(1.0 - c0) * G2;
    
    // Store the results in the struct
    result.position = x;
    result.velocity = v;
    
    return result;
}

// Function to run the simulation
void run_simulation(double w0, double eta, double kB, double mass, double T, gsl_rng *k, double dt, int num_steps, const char *filename)
{
    // Initial conditions (arbitrary values)
    double x = 0.0;  // initial position
    double v = 0.0;  // initial velocity
    
    // Open the file to save the data
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file for writing.\n");
        return;
    }
    
    // Run the simulation for 'num_steps' steps
    for (int i = 0; i < num_steps; i++) {
        // Perform one step of the BD3 algorithm
        result_t result = BD3(x, v, w0, dt, eta, kB, mass, T, k);
        
        // Save the position and velocity data after the burn-in period
        if (i > num_steps / 10) {  // Discard the first 10% of the steps as burn-in
            fprintf(file, "%f\t%f\t%f\n", i*dt, result.position, result.velocity);
        }
        
        // Update x and v for the next step
        x = result.position;
        v = result.velocity;
    }
    
    // Close the file after saving the data
    fclose(file);
}

int main()
{
    // Set up the random number generator
    gsl_rng *k = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(k, 12345); // Set seed for reproducibility
    
    // Constants for the Low case
    int num_steps = 100000;        // Number of simulation steps
    double w0 = 19.477;              // 3.1 kHz to angular frequency (rad/ms)
    double kB = 1.380649e-8;         // Boltzmann constant μm^2μg/ms^2 K
    double mass = 3.0134e-5;         // Mass of the particle micrograms
    double T = 297.0;                // Temperature in K (room temperature)
    double dt_high = 0.05;           // Time step 0.05 ms 
    double dt_low = 0.001;           // Time step 0.001 ms
    double eta_low = 1.0 / 147.3e-3; // Relaxation time: 147.3 µs
    double eta_high = 1.0 / 48.5e-3; // Relaxation time: 48.5 µs
    
    // Run simulation for Low eta high dt
    run_simulation(w0, eta_low, kB, mass, T, k, dt_high, num_steps, "le_hdt.txt");
    // Run simulation for Low eta low dt
    run_simulation(w0, eta_low, kB, mass, T, k, dt_low, num_steps, "le_ldt.txt");
    
    // Run simulation for high eta high dt
    run_simulation(w0, eta_high, kB, mass, T, k, dt_high, num_steps, "he_hdt.txt");
    // Run simulation for high eta low dt
    run_simulation(w0, eta_high, kB, mass, T, k, dt_low, num_steps, "he_ldt.txt");
    

    gsl_rng_free(k);
    
    return 0;
}
