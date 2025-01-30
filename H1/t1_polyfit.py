import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt("lattice_energy.txt")
lattice_constants = data[:, 0]
energies = data[:, 1]

# Perform quadratic fit
fit_coeffs = np.polyfit(lattice_constants, energies, 10)  # Fit a 2nd-degree polynomial
quadratic_fit = np.poly1d(fit_coeffs)

# Generate smooth curve for visualization
lattice_range = np.linspace(lattice_constants.min(), lattice_constants.max(), 500)
fitted_energies = quadratic_fit(lattice_range)

# Calculate equilibrium lattice constant
a_eq = min(fitted_energies)

# Plot the results
plt.figure(figsize=(10, 6))
plt.scatter(lattice_constants, energies, label="Data", color="red")
plt.plot(lattice_range, fitted_energies, label=f"Fitted curve", color="blue")
plt.plot(lattice_range[min(fitted_energies)==fitted_energies], a_eq, '*', color="gold", markersize = 10, label= f"min {a_eq:.3f} eV @  65.42 Å$^3$")
plt.xlabel("Lattice Constant (Å)")
plt.ylabel("Potential Energy (eV)")
plt.title("Lattice Constant vs. Potential Energy")
plt.legend()
plt.grid()
plt.show()

print(f"Theoretical equilibrium lattice constant: {a_eq:.3f} Å")
