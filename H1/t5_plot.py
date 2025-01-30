import numpy as np
import matplotlib.pyplot as plt

# Constants
NSUPER = 256          # Number of supercells
MAL = 26.0/9649.0     # Simulation-Unit mass 
KB = 8.617333262e-5   # Boltzmann constant [eV/K]
KT = 1.0/760000.0     # AL isothermic compressibility as [Bulk Modulus]^-1
#SIMTIME = 10          # simulation time [ps]
#TMELT = 30            # Initial melt time [ps]
TSTABLE = 0.8          # percent for simulation to stabilize [ps]


filename = "Equilibrate_Liquid_dt0.0005.txt"

data = np.loadtxt(filename, delimiter=',')

time = data[:, 0]
potential_energy = data[:, 1]
kinetic_energy = data[:, 2]
total_energy = data[:, 3]
temperature = data[:, 4]
pressure = data[:, 5]
lattice_param = data[:, 6]


dt = time[1]- time[0]
simulation_time_tot = len(time)*dt # Simulation time in
averaging_time = int(len(time)*TSTABLE)  # time to average over simulation

kinetic_mean = np.mean(kinetic_energy[averaging_time:])
average_temperature = np.mean(temperature[averaging_time:])

kinetic_variance = np.mean((potential_energy[averaging_time:] - np.mean(potential_energy[averaging_time:]))**2)
#kinetic_variance = np.mean((kinetic_energy[averaging_time:] - np.mean(kinetic_energy[averaging_time:]))**2)
#kinetic_variance = 0.98
#average_kinetic_sq = (sum(kinetic_energy[sp:]) /len(time[sp:])) **2
#average_sq_kinetic = sum(kinetic_energy[sp:]**2) /len(time[sp:])
#kinetic_variance = average_sq_kinetic - average_kinetic_sq

term1 = (3.0 * NSUPER * KB)/2
term2 = (3.0 * NSUPER * KB * KB * average_temperature*average_temperature)

term3 = 1/ (1 - ((2*kinetic_variance)/term2))

heat_capacity = term1*term3
#print("kinetic_energy_sum_avg =     %.2f\n", average_kinetic_sq)
#print("kinetic_energy_sum_sq_avg =  %.2f\n", average_sq_kinetic)

print("##########")

print("kinetic variance: %f"% kinetic_variance)
print("stabilization period of: %f [ps]"% int(averaging_time*dt))
print("average temperature: %.2f [K]"% average_temperature)
print("heat capacity: %f [eV/K]"% heat_capacity)

#    double term1 = (3.0 * NSUPER * KB)/2;
    #double term2 = 1 /( 1 - ( (kinetic_variance*2.0) / (3.0 * NSUPER * KB * KB * Teq * Teq) ) );