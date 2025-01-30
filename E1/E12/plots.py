import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt("energies.txt")
time = data[:, 0]
energies = data[:, 1:]

# Plot energies for modes 1 to 5
plt.figure(figsize=(10, 6))
for k in range(5):
    plt.subplot(1,2,1)
    plt.plot(time, energies[:, k], label=f"Mode {k+1}", linewidth = 1)
    plt.subplot(1,2,2)
    plt.loglog(time, energies[:, k], label=f"Mode {k+1}", linewidth = 1)

# Add labels and legend
plt.xlabel("Time (t)")
plt.ylabel("Energy $E_k(t)$")
plt.title("Energy of Normal Modes Over Time")
plt.legend()
plt.grid()
plt.show()
