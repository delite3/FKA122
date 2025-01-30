#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "tools.h"
#include <math.h>

// Simulation Parameters
#define N       1000
#define TEQ     2000
#define TEMP    297.0
#define TAUH    147.3e-3 //microseconds
#define TAUL    48.5e-3 //microseconds

// Particle Constants
#define KB 1.380649e-8
#define MASS 3.0134e-5
#define W0 19.477 //rad per ms
#define PI 3.141592653589

// Struct to store the result
typedef struct {
    double position;
    double velocity; 
} result_t;

// BD3 function
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



void write_to_file(char *fname,double** array){
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "x,v,t\n");
    for (int i = 0; i < N; i++) {
       
        fprintf(fp, "%lf, %lf, %lf\n",array[i][0],array[i][1],array[i][2]);
    }
    fclose(fp);
}

int main(){
    const gsl_rng_type *T;
    gsl_rng *k;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    k = gsl_rng_alloc(T);
    gsl_rng_set(k, time(NULL));

    double** trajectory = create_2D_array(N,3);

    double etas[] = {1/TAUL, 1/TAUH};
    double dts[] = {0.001,0.005};

    for(int j = 0;j<2;j++){
        double eta = etas[j];
        for(int p = 0; p<2;p++){
            double dt = dts[p];
            double initial_position = 0.0;
            double initial_velocity = 0.0;
            //double a = -W0*W0*initial_position;
            result_t step = BD3(initial_position, initial_velocity, W0, dt, eta, KB, MASS, TEMP, k);

            for (int i = 0; i < TEQ; i++){
                step = BD3(initial_position, initial_velocity, W0, dt, eta, KB, MASS, TEMP, k);
            }
            //production run:
            for(int i = 0; i<N; i++){

                trajectory[i][0] = step.position;
                trajectory[i][1] = step.velocity;
                trajectory[i][2] = i*dt;
                step = BD3(initial_position, initial_velocity, W0, dt, eta, KB, MASS, TEMP, k);
            }   
            char filename[50];  
            sprintf(filename, "%d%d.csv", j, p);

            // Write to file with the unique filename
            write_to_file(filename, trajectory);
        }

    }
    
    return 0;
}