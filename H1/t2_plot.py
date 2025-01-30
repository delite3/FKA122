import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

dts = [0.0050, 0.0100, 0.0150]
dts = [0.0001, 0.0005, 0.0010, 0.0050, 0.0100, 0.0150]

n = len(dts)
for i in range(n):
    dt = dts[i]
    filename = "energy_time_" + str(format(dt, '.4f')) + ".txt"

    data = np.loadtxt(filename,delimiter=',')
    ttot =len(data)
    nn = int(3*ttot/4)
    nn = 0
    print(ttot)
    time = data[nn:, 0]
    potential_energy = data[nn:, 1]
    kinetic_energy = data[nn:, 2]
    total_energy = data[nn:, 3]
    temperature = data[nn:, 4]
    pressure = data[nn:, 5]
    lattice_constant = data[:, 6]

    # Plot the results
    plt.subplot(n, 2, (1+ i*2))
    if i == 0:
        plt.title("Total Energy [eV]", fontsize=14)
    plt.plot(time, total_energy, color='green')
    plt.ylabel(f"{dt:.4f}ps", rotation=0, labelpad=50, fontsize=14)
    plt.yticks(fontsize=14)
    #plt.xticks(fontsize=10)
    plt.grid(which='both')
    plt.ticklabel_format(useOffset=False)
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))   
    

    plt.subplot(n, 2, (2 +i*2))
    plt.plot(time, temperature, color='red')
    plt.yticks(fontsize=14)
    if i == 0:
        plt.title("Temperature [K]", fontsize=14)
    plt.grid()
    plt.ticklabel_format(useOffset=False)


plt.show()



