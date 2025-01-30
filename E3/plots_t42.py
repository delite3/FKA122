import numpy as np
import matplotlib.pyplot as plt

def compute_ps(v, sampling_rate, dt):
    v = v[::sampling_rate]
    dt = sampling_rate * dt  # New sampling interval
    fft_v = np.fft.fft(v)
    N = len(v)
    freqs = np.fft.fftfreq(N, dt)
    freqs = np.fft.fftshift(freqs)  # Shift frequency axis to center around 0
    P_v = np.abs(fft_v)**2 / N      # Normalize power spectrum
    P_v = np.fft.fftshift(P_v)      # Shift the power spectrum accordingly
    return freqs, P_v

def compute_segmented_pv(v, sampling_rate, M, dt=0.001e-3):
    # Downsample the data
    v_downsampled = v[::sampling_rate]
    dt = sampling_rate * dt  # New sampling interval
    
    segment_length = len(v_downsampled) // M
    P_v_avg = np.zeros(segment_length)
    
    for i in range(M):
        segment = v_downsampled[i*segment_length:(i+1)*segment_length]
        f, P_v = compute_ps(segment, dt)
        P_v_avg += P_v
    P_v_avg /= M  # Average the power spectra
    
    return f, P_v_avg


##### SIMULATION FFT SOLUTION ######
data_low = np.loadtxt('le_ldt.txt')  # Relaxation time: 147.3 µs
data_high = np.loadtxt('he_ldt.txt') # Relaxation time: 48.5 µs

v_low = data_low[:, 2]  
v_high = data_high[:, 2] 

dt = data_low[0,0] - data_low[1,0]


# Compute power spectrum for both time spacings
lf_50, lP_v_50 = compute_ps(v_low, 50, dt=dt)  
lf_25, lP_v_25 = compute_ps(v_low, 25, dt=dt)  
hf_50, hP_v_50 = compute_ps(v_high, 50, dt=dt) 
hf_25, hP_v_25 = compute_ps(v_high, 25, dt=dt) 


##### ANALYTICAL FFT SOLUTION ######
# Constants 
k_B = 1.380649e-8       # Boltzmann constant μm^2μg/ms^2 K
T = 297.0               # Temperature (K)
mass = 3.0134e-5        # Mass of the particle micrograms
w0 = 19.477             # Natural frequency (rad/ms)
eta_low = 1.0 / 147.3e-3    # Damping coefficient 
eta_high = 1.0 / 48.5e-3

# Frequency range for the analytical spectrum
omega_min = 0       # Minimum angular frequency (rad/s)
omega_max = 10000  # Maximum angular frequency (rad/s)
num_points = 100000   # Number of frequency points

omega = np.linspace(omega_min, omega_max, num_points)

# Analytical power spectrum P_x(ω) for a damped harmonic oscillator
lP_v_analytical = (k_B * T / mass) * (2 * eta_low*omega**2) / ((omega**2 - w0**2)**2 + (eta_low**2) * omega**2)
hP_v_analytical = (k_B * T / mass) * (2 * eta_high*omega**2) / ((omega**2 - w0**2)**2 + (eta_high**2) * omega**2)

# Plot for Low
plt.subplot(1,2,1)
plt.plot(lf_50, lP_v_50, label='0.05ms', color='blue')
plt.plot(lf_25, lP_v_25, label='0.025ms', color='green')
plt.plot(omega/(2*np.pi), lP_v_analytical*10, label='Analytical Power Spectrum', color='black')

plt.xlim(0, 20)
plt.xlabel('Frequency (Hz *$10^3$)')
plt.ylabel('Power Spectrum $P_v(f)$')
plt.title('2.75kPa')
plt.legend()

# Plot for High
plt.subplot(1,2,2)
plt.plot(hf_50, hP_v_50, label='0.05ms', color='blue')
plt.plot(hf_25, hP_v_25, label='0.025ms', color='green')
plt.plot(omega/(2*np.pi), hP_v_analytical*10, label='Analytical Power Spectrum', color='black')

plt.xlim(0, 20)
plt.xlabel('Frequency (Hz *$10^3$)')
plt.ylabel('Power Spectrum $P_v(f)$')
plt.title('99.8kPa')
plt.legend()


plt.show()
