import numpy as np
import matplotlib.pyplot as plt

# Constants (replace with your actual values)
k_B = 1.380649e-8       # Boltzmann constant μm^2μg/ms^2 K
T = 297.0               # Temperature (K)
mass = 3.0134e-5        # Mass of the particle micrograms
w0 = 19.477             # Natural frequency (rad/ms)
eta = 1.0 / 147.3e-3    # Damping coefficient 

# Frequency range for the analytical spectrum
omega_min = 0       # Minimum angular frequency (rad/s)
omega_max = 100000  # Maximum angular frequency (rad/s)
num_points = 100000   # Number of frequency points
omega = np.linspace(omega_min, omega_max, num_points)

# Analytical power spectrum P_x(ω) for a damped harmonic oscillator
P_x_analytical = (k_B * T / mass) * (2 * eta*omega**2) / ((omega**2 - w0**2)**2 + eta**2 * omega**2)

# Plot the analytical power spectrum
plt.figure(figsize=(10, 6))
plt.plot(omega/(2*np.pi), P_x_analytical, label='Analytical Power Spectrum', color='b')
plt.xlabel('Angular Frequency (rad/s)')
plt.ylabel('Power Spectrum $P_x(\omega)$')
plt.title('Analytical Power Spectrum for Damped Harmonic Oscillator')
plt.legend()
#plt.xlim(0, 20)  # Adjust frequency range as needed (to 10 kHz for example)
plt.show()