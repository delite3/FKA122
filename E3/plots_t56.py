import numpy as np
import matplotlib.pyplot as plt

def compute_ps(v, dt):
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

def compute_autocorrelation(v):
    N = len(v)
    v_mean = np.mean(v)
    v_fluct = v - v_mean
    # Compute ACF via FFT for efficiency:
    # ACF = IFFT(FFT(v)*conjugate(FFT(v))) normalized appropriately
    fft_v = np.fft.fft(v_fluct, n=2*N)  # zero-padding to avoid wrap-around
    acf_full = np.fft.ifft(np.abs(fft_v)**2)
    acf_full = np.real(acf_full[:N])  # take first N points
    # Normalize by number of points and variance
    acf = acf_full / (N - np.arange(N))
    acf /= acf[0]  # normalize so ACF(0)=1
    return acf

def acf_to_psd(acf, dt):
    # FFT of the ACF gives PSD
    fft_acf = np.fft.fft(acf)
    psd = np.real(fft_acf) * dt
    N = len(acf)
    freqs = np.fft.fftfreq(N, dt)
    return freqs, psd

##### SIMULATION FFT SOLUTION ######
data_low = np.loadtxt('le_ldt.txt')  # Relaxation time: 147.3 µs
data_high = np.loadtxt('he_ldt.txt') # Relaxation time: 48.5 µs

v_low = data_low[:, 2]  
v_high = data_high[:, 2] 

dt = data_low[0,0] - data_low[1,0]

# Compute power spectrum for both time spacings
lf_50, lP_v_50 = compute_segmented_pv(v_low, 50, 100, dt=dt)  
lf_25, lP_v_25 = compute_segmented_pv(v_low, 25, 100, dt=dt)  
hf_50, hP_v_50 = compute_segmented_pv(v_high, 50, 100, dt=dt) 
hf_25, hP_v_25 = compute_segmented_pv(v_high, 25, 100, dt=dt) 

# ACF FFT  
acf_low = compute_autocorrelation(v_low)
acf_high = compute_autocorrelation(v_high)
# Compute PSD from ACF
freqs_acf_low, psd_acf_low = compute_segmented_pv(acf_low, 50, 100, dt=dt)
freqs_acf_high, psd_acf_high = compute_segmented_pv(acf_high, 50, 100, dt=dt)

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
plt.plot(freqs_acf_low, psd_acf_low, label='From ACF (Wiener-Khintchin)')

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
plt.plot(freqs_acf_high, psd_acf_high, label='From ACF (Wiener-Khintchin)')

plt.xlim(0, 20)
plt.xlabel('Frequency (Hz *$10^3$)')
plt.ylabel('Power Spectrum $P_v(f)$')
plt.title('99.8kPa')
plt.legend()


plt.show()
