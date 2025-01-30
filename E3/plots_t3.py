import numpy as np
import matplotlib.pyplot as plt

he_hdt = np.loadtxt('he_hdt.txt') #99.8kPa
he_ldt = np.loadtxt('he_ldt.txt') #99.8kPa
le_hdt = np.loadtxt('le_hdt.txt') #2.75kPa
le_ldt = np.loadtxt('le_ldt.txt') #2.75kPa

timehdt = he_hdt[:,0] - he_hdt[0,0]
timeldt = he_ldt[:,0] - he_ldt[0,0]
timehdt = timehdt[timehdt<2]
timeldt = timeldt[timeldt<2]

he_hdt = he_hdt[:len(timehdt)]
he_ldt = he_ldt[:len(timeldt)]
le_hdt = le_hdt[:len(timehdt)]
le_ldt = le_ldt[:len(timeldt)]

fig, ax = plt.subplots(2, 2)
fig.canvas.manager.set_window_title("Positions - Pressure and Time Step Comparisons")

ax[0,0].plot(timehdt, he_hdt[:,1], label="99.8 kPa, dt = 0.001")
ax[1,0].plot(timeldt, he_ldt[:,1], label="99.8 kPa, dt = 0.005")
ax[0,1].plot(timehdt, le_hdt[:,1], label="2.75 kPa, dt = 0.001")
ax[1,1].plot(timeldt, le_ldt[:,1], label="2.75 kPa, dt = 0.005")

for i in range(2):
    for j in range(2):
        ax[i,j].set_xlabel("Time (ms)")
        ax[i,j].set_ylabel("position (micrometers)")
        ax[i,j].grid()
        ax[i,j].legend()

plt.show()
