import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

#TMELT = 0.1             # Initial melt time [ps]
TSTABLE = 0.25           # Time for simulation to stabilize [ps]

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


# Plot the results
plt.subplot(2, 2, 1)
plt.title("Total Energy (eV)", fontsize=14)
plt.plot(time, total_energy, color='green')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid()
plt.ticklabel_format(useOffset=False)
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))   

plt.subplot(2, 2, 2)
plt.plot(time, lattice_param, color='blue')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title("Lattice parameter [Å]", fontsize=14)
plt.grid()
plt.ticklabel_format(useOffset=False)

plt.subplot(2, 2, 3)
plt.plot(time, temperature, color='orange')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title("Temperature [K]", fontsize=14)
plt.grid()
plt.ticklabel_format(useOffset=False)


plt.subplot(2, 2, 4)
plt.plot(time, pressure, color='red')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title("Pressure [Bar]", fontsize=14)
plt.grid()
plt.ticklabel_format(useOffset=False)

plt.show()

cell_volume = 4*lattice_param
dt = time[1]- time[0]

dt = time[1]- time[0]
simulation_time_tot = len(time)*dt # Simulation time in
averaging_time = int(len(time)*TSTABLE)  # time to average over simulation


average_temperature = sum(temperature[averaging_time:])/len(time[averaging_time:])
average_pressure = sum(pressure[averaging_time:])/len(time[averaging_time:])
average_lattice = sum(lattice_param[averaging_time:])/len(time[averaging_time:])

print("stabilization period of : %ips"% simulation_time_tot)
print("average temperature: %.2f"% average_temperature)
print("averag pressure: %.2f"% average_pressure)
print("average lattice parmeter: %.4f"% average_lattice)

plt.subplot(1,2,1)
for i in range(3):
    plt.plot(time, at1[i])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel(f"Position [Å]", fontsize=16)
    plt.legend(['x', 'y', 'z'], fontsize=14, loc= 'center right')
    #plt.xlabel('time in ps')

plt.subplot(1,2,2)
for i in range(3):
    plt.plot(time, at2[i])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel(f"Position [Å]", fontsize=16,)
    plt.legend(['x', 'y', 'z'], fontsize=14, loc= 'center right')
    #plt.xlabel('time in ps')

plt.show()



MSD_at1 = (at1[0][:] - at1[0][0])**2 + (at1[1][:] - at1[1][0])**2 + (at1[2][:] - at1[2][0])**2
MSD_at2 = (at2[0][:] - at2[0][0])**2 + (at2[1][:] - at2[1][0])**2 + (at2[2][:] - at2[2][0])**2

MSD = (MSD_at1 + MSD_at2)/2

plt.plot(np.sqrt(MSD))
plt.show()