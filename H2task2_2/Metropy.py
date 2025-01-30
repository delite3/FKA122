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

# Separate cold and warm start data
cold_data = data[data["Start"] == "Cold"]
warm_data = data[data["Start"] == "Warm"]

# Debugging: Print separated datasets
print("\nCold Start Data:")
print(cold_data)

print("\nWarm Start Data:")
print(warm_data)

# Ensure no empty dataframes
if cold_data.empty:
    print("Warning: No data for Cold Start!")
if warm_data.empty:
    print("Warning: No data for Warm Start!")

# Verify and cast columns as numeric
cold_data["Temp"] = pd.to_numeric(cold_data["Temp"])
cold_data["InitialEnergy"] = pd.to_numeric(cold_data["InitialEnergy"])
cold_data["AverageEnergy"] = pd.to_numeric(cold_data["AverageEnergy"])
warm_data["Temp"] = pd.to_numeric(warm_data["Temp"])
warm_data["InitialEnergy"] = pd.to_numeric(warm_data["InitialEnergy"])
warm_data["AverageEnergy"] = pd.to_numeric(warm_data["AverageEnergy"])

# Debugging: Confirm the exact data to be plotted
print("\nCold Start Temperatures:", cold_data["Temp"].tolist())
print("Cold Start Initial Energies:", cold_data["InitialEnergy"].tolist())
print("Cold Start Average Energies:", cold_data["AverageEnergy"].tolist())
print("\nWarm Start Temperatures:", warm_data["Temp"].tolist())
print("Warm Start Initial Energies:", warm_data["InitialEnergy"].tolist())
print("Warm Start Average Energies:", warm_data["AverageEnergy"].tolist())

# Plot initial and average energies for each temperature
plt.figure(figsize=(10, 6))

# Plot Cold Start Data
plt.plot(
    cold_data["Temp"],
    cold_data["InitialEnergy"],
    'o-',
    label="Cold Start - Initial Energy"
)
plt.plot(
    cold_data["Temp"],
    cold_data["AverageEnergy"],
    'o-',
    label="Cold Start - Average Energy"
)

# Plot Warm Start Data
plt.plot(
    warm_data["Temp"],
    warm_data["InitialEnergy"],
    's--',
    label="Warm Start - Initial Energy"
)
plt.plot(
    warm_data["Temp"],
    warm_data["AverageEnergy"],
    's--',
    label="Warm Start - Average Energy"
)

# Customize the plot
plt.title("Simulation Results: Energy vs Temperature")
plt.xlabel("Temperature (K)")
plt.ylabel("Energy (meV)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show(block=True)
