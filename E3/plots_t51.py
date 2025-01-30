import matplotlib.pyplot as plt
import numpy as np

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

# Example usage with velocity data (v)
data_high = np.loadtxt('he_ldt.txt') # Low eta, low dt
data_low = np.loadtxt('le_ldt.txt') # High eta, low dt
v_high = data_high[:, 2]  # velocity
v_low = data_low[:, 2]  # position

time_low = data_low[:, 0]
ldt = time_low[1] - time_low[0]

# Set maximum lag time for correlation function 
max_lag = 1000
C_v_low = compute_autocorrelation(v_low)
C_v_high = compute_autocorrelation(v_high)

plt.plot(np.arange(max_lag) * ldt, C_v_low[:max_lag], label='2.75kPa')
plt.plot(np.arange(max_lag) * ldt, C_v_high[:max_lag], label='99.8kPa')
plt.legend(loc="upper right")
plt.xlabel('Time (s)')
plt.ylabel('Time Correlation Function $C_v(t)$')
plt.title('Time Correlation Function for Low eta')
plt.show()
