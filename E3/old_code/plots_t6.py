import numpy as np
import matplotlib.pyplot as plt

def compute_time_correlation(v, max_lag):
    N = len(v)
    correlation = np.zeros(max_lag)  # Array to store the correlation values
    
    # Average over all possible starting times t_i
    for t in range(N - max_lag):
        for lag in range(max_lag):
            correlation[lag] += v[t + lag] * v[t]
    
    # Normalize by the number of averages and the initial velocity squared
    correlation /= (N - max_lag) * np.var(v)
    
    return correlation

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

# Example usage with velocity data (v)
data_high = np.loadtxt('he_ldt.txt')  # Low eta, low dt
data_low = np.loadtxt('le_ldt.txt')  # High eta, low dt
v_high = data_high[:, 2]  # velocity
v_low = data_low[:, 2]  # position

time_low = data_low[:, 0]
dt = time_low[1] - time_low[0]

# Set maximum lag time for correlation function (e.g., max_lag = 1000)
max_lag = 1000
#C_v_low = compute_time_correlation(v_low, max_lag)
C_v_low = np.loadtxt('C_v_low.txt')  # Low eta, low dt


# Compute the power spectrum from the time correlation function
freqs_from_correlation, P_v_from_correlation = compute_ps(C_v_low, 25*dt)

# Compute the power spectrum directly from the velocity signal for comparison
freqs_direct, P_v_direct = compute_ps(v_low, 25*dt)

freqs_segmented, P_v_segmented = compute_segmented_pv(v_low, 25, 100, dt=dt)


# Plot direct power spectrum from the velocity signal
#plt.scatter(freqs_direct, P_v_direct, marker='.', c='darkgreen', label='Direct Power Spectrum', linestyle='--')

# Plot segmented power spectrum from the velocity signal
#plt.plot(freqs_segmented, P_v_segmented, '-r', label='Segmented Power Spectrum')

# Plot power spectrum from time correlation
plt.plot(freqs_from_correlation, P_v_from_correlation, '-b',label='From Time Correlation')


plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectrum $P_v(f)$')
plt.title('Power Spectrum Comparison')
plt.legend()

#plt.xlim(0, 20)  # Limit frequency range to 0 < f < 10 kHz
plt.show()
