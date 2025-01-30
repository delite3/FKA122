import pandas as pd
import matplotlib.pyplot as plt

# Read the simulation results CSV
try:
    data = pd.read_csv("simulation_results.csv")
except FileNotFoundError:
    print("Error: simulation_results.csv not found!")
    exit()

# Debugging: Print the raw data from the CSV file
print("Full CSV Data:")
print(data)

# Ensure required columns exist
required_columns = {"Start", "Temp", "InitialEnergy", "AverageEnergy"}
if not required_columns.issubset(data.columns):
    print(f"Error: Missing columns in CSV. Expected columns: {required_columns}")
    exit()

# Separate cold start data
cold_data = data[data["Start"] == "Cold"]

# Debugging: Print Cold Start dataset
print("\nCold Start Data:")
print(cold_data)

# Ensure no empty data
if cold_data.empty:
    print("Error: No data for Cold Start!")
    exit()

# Verify and cast columns as numeric
cold_data["Temp"] = pd.to_numeric(cold_data["Temp"])
cold_data["AverageEnergy"] = pd.to_numeric(cold_data["AverageEnergy"])

# Debugging: Confirm the exact data to be plotted
print("\nCold Start Temperatures:", cold_data["Temp"].tolist())
print("Cold Start Average Energies:", cold_data["AverageEnergy"].tolist())

# Plot Cold Start Average Energy
plt.figure(figsize=(10, 6))
plt.plot(
    cold_data["Temp"],
    cold_data["AverageEnergy"],
    'o-',
    label="Cold Start - Average Energy"
)

# Customize the plot
plt.title("Cold Start: Average Energy vs Temperature")
plt.xlabel("Temperature (K)")
plt.ylabel("Energy (meV)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show(block=True)
