import numpy as np
import matplotlib.pyplot as plt

def compute_ps(v, dt):
    # Compute FFT of the signal
    fft_v = np.fft.fft(v)
    N = len(v)
    # Frequency axis (shifted for visualization)
    freqs = np.fft.fftfreq(N, dt)
    freqs = np.fft.fftshift(freqs)  # Shift frequency axis to center around 0
    P_v = np.abs(fft_v)**2 / N  # Normalize power spectrum
    P_v = np.fft.fftshift(P_v)  # Shift the power spectrum accordingly
    return freqs, P_v

def compute_segmented_pv(v, dt, M):
    segment_length = len(v) // M
    P_v_avg = np.zeros(segment_length)
    for i in range(M):
        segment = v[i*segment_length:(i+1)*segment_length]
        f, P_v = compute_ps(segment, dt)
        P_v_avg += P_v
    P_v_avg /= M  # Average the power spectra
    return f, P_v_avg

# Load the velocity data (from file, for example)
data_low = np.loadtxt('le_ldt.txt')  # Relaxation time: 147.3 µs
data_high = np.loadtxt('he_ldt.txt') # Relaxation time: 48.5 µs

v_low = data_low[:, 2]  
v_high = data_high[:, 2] 

tau50 = 0.05 # 0.05ms
tau25 = 0.025 # 0.025ms

# Compute power spectrum for both time spacings
lf_50, lP_v_50 = compute_segmented_pv(v_low, tau50, 100)  
lf_25, lP_v_25 = compute_segmented_pv(v_low, tau25, 100)  
hf_50, hP_v_50 = compute_segmented_pv(v_high, tau50, 100)  
hf_25, hP_v_25 = compute_segmented_pv(v_high, tau25, 100)

# Plot for delta Tau = 50 * dt
plt.subplot(1,2,1)
#plt.plot(omega, P_x_analytical*10**5, label='Analytical', color='black')
plt.plot(lf_50, lP_v_50, label='0.05ms', color='blue')
plt.plot(lf_25, lP_v_25, label='0.025ms', color='green')

#plt.xlim(0, 5e3)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectrum $P_v(f)$')
plt.title('2.75kPa')
plt.legend()

# Plot for delta Tau = 25 * dt
plt.subplot(1,2,2)
plt.plot(hf_50, hP_v_50, label='0.05ms', color='blue')
plt.plot(hf_25, hP_v_25, label='0.025ms', color='green')

plt.xlim(0, 5e3)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectrum $P_v(f)$')
plt.title('99.8kPa')
plt.legend()

plt.show()

print("low case 25:")
print(lf_25[lP_v_25[:]>5])
print("low case 50:")
print(hf_50[lP_v_50[:]>5])



