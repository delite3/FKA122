import numpy as np
import matplotlib.pyplot as plt

def generate_signal(a, phi, f, dt, N):
    """
    Generate signal for task in C4b; h(t) = a * cos( 2 * pi * f * t + phi)
    :param a:   Amplitude
    :param phi: Phase
    :param f:   Frequency
    :param dt:  Time-step
    :param N:   Number of steps
    :return:    time-array, signal-array
    """
    t = np.arange(0, N*dt, dt)
    h = a*np.cos(2*np.pi*f*t + phi)
    return t, h


def my_dft(h):
    """
    Compute the Discrete Fourier Transform (DFT) of a signal.
    :param h: Real-valued input signal (1D array)
    :return: DFT of the signal (complex array)
    """
    N = len(h)
    H = np.zeros(N, dtype=complex)
    for n in range(N):
        for k in range(N):
            H[n] += h[k] * np.exp(-2j * np.pi * n * k / N)
    return H

def my_aps(h, dt):
    """
    Compute the asymmetric power spectrum of a signal using the DFT.
    Includes redundant positive frequencies for visualization.
    :param h: Real-valued input signal (1D array)
    :param dt: Time step between samples
    :return: Frequencies (Hz), Power spectrum
    """
    H = my_dft(h)  # Compute DFT
    N = len(h)  # Number of points
    P = (np.abs(H) ** 2) / N  # Compute power spectrum (asymmetric, unshifted)
    
    # Compute the actual frequency values
    freqs = np.arange(N) / (N * dt)
    
    return freqs, P


def my_sps(h, dt):
    """
    Compute the symmetric power spectrum of a signal using the DFT.
    Manually shift the spectrum for symmetric representation using a for-loop.
    :param h: Real-valued input signal (1D array)
    :param dt: Time step between samples
    :return: Shifted frequencies, Symmetric power spectrum
    """
    H = my_dft(h)  # Compute the DFT
    N = len(h)  # Number of points

    # Compute power spectrum
    P = (np.abs(H) ** 2) / N

    # Compute unshifted frequencies
    freqs = np.arange(N) / (N * dt)

    # Create arrays for symmetric spectrum and frequencies
    P_sym = np.zeros_like(P)
    freqs_sym = np.zeros_like(freqs)

    # Perform manual shifting
    mid = N // 2
    for i in range(N):
        if i < mid:
            # Negative frequencies  
            P_sym[i] = P[i + mid]
            freqs_sym[i] = (i - mid) / (N * dt)
        else:
            # Positive frequencies
            P_sym[i] = P[i - mid]
            freqs_sym[i] = (i - mid) / (N * dt)

    return freqs_sym, P_sym


def numpy_fft(h, dt):
    """
    Compute the Fourier transform and power spectrum of a signal using NumPy's FFT.
    
    :param h: Real-valued input signal (1D array)
    :param dt: Time step between samples
    :return: frequencies, power_spectrum (symmetric representation)
    """
    # Compute the FFT of the signal
    H = np.fft.fft(h)
    
    # Compute the number of points
    N = len(h)
    
    # Calculate the power spectrum
    P = (np.abs(H) ** 2) / N
    
    # Generate frequencies using FFT frequencies
    freqs = np.fft.fftfreq(N, dt)
    
    # Shift the frequencies and power spectrum for symmetric representation
    freqs_sym = np.fft.fftshift(freqs)
    P_sym = np.fft.fftshift(P)
    
    return freqs_sym, P_sym


# Generate a cosine signal
a=1
phi=0
dt=0.02
N=253
_, h1 = generate_signal(a=a, phi=phi, f=2, dt=dt, N=N)

_, h2 = generate_signal(a=a, phi=phi, f=6, dt=dt, N=N)

h = h1 

# Asymmetric spectrum
freqs, P_asym = numpy_fft(h, dt)
plt.figure(figsize=(10, 4))
plt.stem(freqs, P_asym, basefmt=" ")
plt.title("Asymmetric Power Spectrum")
plt.xlabel("Frequency index")
plt.ylabel("Power")
plt.grid()
plt.show()

# Symmetric spectrum
freqs_sym, P_sym = my_sps(h, dt)
plt.figure(figsize=(10, 4))
plt.stem(freqs_sym, P_sym, basefmt=" ")
plt.title("Symmetric Power Spectrum")
plt.xlabel("Frequency")
plt.ylabel("Power")
plt.grid()
plt.show()

print("")
