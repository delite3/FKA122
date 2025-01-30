import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt("energies.txt")
time = data[:, 0]
energies = data[:, 1:]

# Initialize variables for time-averaging
cumulative_energies = np.zeros(energies.shape[1])
time_averages = []

# Compute time-averaged energies for each time step
for i, t in enumerate(time):
    cumulative_energies += energies[i]  # Accumulate energies
    time_average = cumulative_energies / (i + 1)  # Time-average up to this point
    time_averages.append(time_average)

# Convert to NumPy array for easier slicing
time_averages = np.array(time_averages)

# Plot energies for modes 1 to 5
plt.figure(figsize=(12, 6))

# Linear scale plot
plt.subplot(1, 2, 1)
for k in range(31):
    plt.plot(time, time_averages[:, k], label=f"Mode {k+1}", linewidth=1)
plt.xlabel("Time (t)")
plt.ylabel("Time-Averaged Energy $\\langle E_k \\rangle_{time}(t)$")
plt.title("Time-Averaged Energies (Linear Scale)")
#plt.legend()
plt.grid()

# Log-log plot
plt.subplot(1, 2, 2)
for k in range(31):
    plt.loglog(time, time_averages[:, k], label=f"Mode {k+1}", linewidth=1)
plt.xlabel("Time (t)")
plt.ylabel("Time-Averaged Energy $\\langle E_k \\rangle_{time}(t)$")
plt.title("Time-Averaged Energies (Log-Log Scale)")
#plt.legend()
plt.grid()

# Show plots
plt.tight_layout()
plt.show()
