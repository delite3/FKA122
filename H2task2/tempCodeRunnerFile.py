import pandas as pd
import matplotlib.pyplot as plt

# Load data
data = pd.read_csv('results_mean_field.csv')

# Plot magnetization
plt.figure()
plt.plot(data['T'], data['m'], label='Magnetization m(T)')
plt.axvline(x=4.0, color='red', linestyle='--', label='Critical Temperature Tc')
plt.xlabel('Temperature (T)')
plt.ylabel('Magnetization m(T)')
plt.title('Magnetization vs Temperature')
plt.legend()
plt.grid()

# Plot energy
plt.figure()
plt.plot(data['T'], data['U'], label='Energy U(T)')
plt.xlabel('Temperature (T)')
plt.ylabel('Energy U(T)')
plt.title('Energy vs Temperature')
plt.legend()
plt.grid()

# Plot heat capacity
plt.figure()
plt.plot(data['T'], data['C'], label='Heat Capacity C(T)')
plt.axvline(x=4.0, color='red', linestyle='--', label='Critical Temperature Tc')
plt.xlabel('Temperature (T)')
plt.ylabel('Heat Capacity C(T)')
plt.title('Heat Capacity vs Temperature')
plt.legend()
plt.grid()

plt.show()
