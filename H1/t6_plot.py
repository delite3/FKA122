import numpy as np
import matplotlib.pyplot as plt

# Load g(r) data
data = np.loadtxt('gofr.txt',delimiter=',')
r = data[:, 0]
g = data[:, 1]

# Known parameters
# Replace this with your system's values:
rho = data[0,2]  # Example number density (units must match those used in gofr.txt)
# If you know delta_r:
delta_r = r[1] - r[0]

# Find the first peak (maximum) in g(r)
peak_index = np.argmax(g)
# Find the first minimum after the peak
min_index = peak_index + np.argmin(g[peak_index:])

r_m = r[min_index]

print(f"First minimum found at r_m = {r_m}")

# Compute the coordination number I(r_m) by approximating the integral using a simple Riemann sum
# Identify all bins up to r_m
bins_up_to_rm = np.where(r <= r_m)[0]
I_rm = 0.0
for i in bins_up_to_rm:
    r_bin = r[i]
    I_rm += g[i] * (r_bin**2) * delta_r

I_rm *= 4.0 * np.pi * rho

print(f"Coordination number I(r_m) = {I_rm}")

# Now we have r, g(r), r_m, and I(r_m).
plt.figure(figsize=(8,6))
plt.plot(r, g, label='g(r)')
plt.axvline(r_m, color='red', linestyle='--', label='r_m')
plt.xlabel('r')
plt.ylabel('g(r)')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title('Radial Distribution Function and First Minimum', fontsize=16)
plt.legend(fontsize=16)
plt.grid(True)
plt.show()
