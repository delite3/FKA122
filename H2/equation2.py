import pandas as pd
import matplotlib.pyplot as plt

# File to load
file_name = 'results_mean_field.csv'  # Update this if the file name changes

# Load data
try:
    data = pd.read_csv(file_name)
except FileNotFoundError:
    print(f"Error: File '{file_name}' not found.")
    exit()

# Verify column names in the CSV file
print("Columns in the file:", data.columns)

# Ensure the columns are as expected
required_columns = ['T', 'm', 'U', 'C']
for col in required_columns:
    if col not in data.columns:
        print(f"Error: Missing column '{col}' in the CSV file.")
        exit()

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
plt.plot(data['T'], data['U'], label='Energy U(T)', color='orange')
plt.xlabel('Temperature (T)')
plt.ylabel('Energy U(T)')
plt.title('Energy vs Temperature')
plt.legend()
plt.grid()

# Plot heat capacity
plt.figure()
plt.plot(data['T'], data['C'], label='Heat Capacity C(T)', color='green')
plt.axvline(x=4.0, color='red', linestyle='--', label='Critical Temperature Tc')
plt.xlabel('Temperature (T)')
plt.ylabel('Heat Capacity C(T)')
plt.title('Heat Capacity vs Temperature')
plt.legend()
plt.grid()

plt.show()
