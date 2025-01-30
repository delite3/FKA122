import numpy as np
import matplotlib.pyplot as plt

# Constants
NSUPER = 256          # Number of atoms
MAL = 26.0/9649.0     # Simulation-Unit mass (not directly needed for g(r))
KB = 8.617333262e-5   # Boltzmann constant [eV/K]
KT = 1.0/760000.0     # AL isothermic compressibility (not directly needed for g(r))
TSTABLE = 0.8         # percent for simulation to stabilize [ps]

filename = "Positions_Solid_dt0.0005.txt"

data = np.loadtxt(filename, delimiter=',')
# data format: [time, CL, X, Y, Z]

time = data[:,0]
CL = data[:,1]
XValues = data[:,2]
YValues = data[:,3]
ZValues = data[:,4]

# Determine dt and number of steps:
dt = time[NSUPER] - time[0]  # time difference between two snapshots
steps = int((time[-1] - time[0])/dt) + 1

# Number of bins
k = 200
# We'll base delta_r on the box size at the first step:
L_initial = CL[0]
delta_r = L_initial / (2.0 * k)  # max r is L/2, so delta_r = (L/2)/k

bin_counts = np.zeros(k)

# Loop over time steps
for t in range(steps):
    # For each time step, atoms occupy rows [t*NSUPER : (t+1)*NSUPER]
    idx_offset = t * NSUPER
    L = CL[idx_offset]  # cell length for this snapshot

    # Compute pairs of atoms
    for atom1 in range(NSUPER):
        x1 = XValues[idx_offset + atom1]
        y1 = YValues[idx_offset + atom1]
        z1 = ZValues[idx_offset + atom1]
        for atom2 in range(atom1+1, NSUPER):
            x2 = XValues[idx_offset + atom2]
            y2 = YValues[idx_offset + atom2]
            z2 = ZValues[idx_offset + atom2]

            dx = x2 - x1
            dy = y2 - y1
            dz = z2 - z1

            # Minimum image convention
            if dx > L/2: dx -= L
            elif dx < -L/2: dx += L
            if dy > L/2: dy -= L
            elif dy < -L/2: dy += L
            if dz > L/2: dz -= L
            elif dz < -L/2: dz += L

            r = np.sqrt(dx*dx + dy*dy + dz*dz)

            bin_index = int(r / delta_r)
            if bin_index < k:
                # Each pair contributes to two atoms' local environment
                # but for g(r) normalization, counting 2 each time is a known convention
                bin_counts[bin_index] += 2.0

# Now we have bin_counts for all pairs over all steps
# Normalize to get g(r):

# Total number of snapshots is 'steps'
# The volume is L^3, but L might vary with time. For simplicity,
# if L doesn't vary much, you could use an average L or do a weighted average.
# Ideally, you should compute g(r) time-step by time-step and average.
# For demonstration, let's just use the average box size:

average_L = np.mean(CL[::NSUPER])  # Average L over time
volume = average_L**3
density = NSUPER / volume

# Save g(r) to file
with open("gofr.txt", "w") as f:
    for bin_index in range(k):
        r_bin = (bin_index + 0.5)*delta_r
        shell_volume = 4.0 * np.pi * (r_bin**2) * delta_r
        ideal_count = density * shell_volume * NSUPER * steps
        g_of_r = bin_counts[bin_index] / ideal_count if ideal_count != 0 else 0.0
        f.write(f"{r_bin} {g_of_r}\n")

# Optionally, plot g(r)
r_vals = (np.arange(k) + 0.5)*delta_r
# Compute g(r) again for plotting
g_vals = []
for bin_index in range(k):
    r_bin = (bin_index + 0.5)*delta_r
    shell_volume = 4.0 * np.pi * (r_bin**2) * delta_r
    ideal_count = density * shell_volume * NSUPER * steps
    g_val = bin_counts[bin_index]/ideal_count if ideal_count != 0 else 0.0
    g_vals.append(g_val)

plt.plot(r_vals, g_vals)
plt.xlabel("r [units]")
plt.ylabel("g(r)")
plt.title("Radial Distribution Function")
plt.show()
