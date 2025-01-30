import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
# Load the data from the file
data = np.loadtxt('trajectory.dat')

# Extract columns
time = data[:, 0]
x_position = data[:, 1]
y_position = data[:, 2]

# Plot the trajectory
plt.figure(figsize=(10, 6))
plt.plot(x_position, y_position, label="Object Trajectory", color="b", markersize=3)

# Add labels and a title
plt.xlabel("X Position")
plt.ylabel("Y Position")
plt.title("Trajectory")
plt.legend()
plt.grid()

# Show the plot
plt.show()
