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


# Case 1: a=1, phi = 0, f=2, dt = 0.1, N = 250
t1, h1 = generate_signal(a=1, phi=0, f=2, dt=0.1, N=250)

# Case 1: a=1, phi = 0, f=2, dt = 0.1, N = 250
t2, h2 = generate_signal(a=1, phi=0, f=1, dt=0.1, N=250)

# Case 1: a=1, phi = pi/2, f=2, dt = 0.1, N = 250
t3, h3 = generate_signal(a=1, phi=np.pi/2, f=2, dt=0.1, N=250)

# Plot results
plt.figure(figsize=(12, 6))

plt.subplot(3, 1, 1)
plt.plot(t1, h1)
plt.title("Case 1: a=1, f=2, phi=0")
plt.xlabel("Time (t)")
plt.ylabel("h(t)")

plt.subplot(3, 1, 2)
plt.plot(t2, h2)
plt.title("Case 2: a=1, f=1, phi=0")
plt.xlabel("Time (t)")
plt.ylabel("h(t)")

plt.subplot(3, 1, 3)
plt.plot(t3, h3)
plt.title("Case 3: a=1, f=2, phi=π/2")
plt.xlabel("Time (t)")
plt.ylabel("h(t)")

plt.tight_layout()
plt.show()

#####


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


# Parameters for Case 1
a, f, phi, dt, N = 1, 2, 0, 0.1, 250
_, h = generate_signal(a, phi, f, dt, N)

import matplotlib.pyplot as plt

# Generate a cosine signal
N = 250
timestep = 0.1
t = np.arange(N) * timestep
f = 2  # Frequency in Hz
phi = 0  # Phase shift
h = np.cos(2 * np.pi * f * t + phi)

# Asymmetric spectrum
freqs, P_asym = my_aps(h,dt)
plt.figure(figsize=(10, 4))
plt.stem(freqs, P_asym, basefmt=" ")
plt.title("Asymmetric Power Spectrum")
plt.xlabel("Frequency index")
plt.ylabel("Power")
plt.grid()
plt.show()

# Symmetric spectrum
freqs_sym, P_sym = my_sps(h,dt)
plt.figure(figsize=(10, 4))
plt.stem(freqs_sym, P_sym, basefmt=" ", use_line_collection=True)
plt.title("Symmetric Power Spectrum")
plt.xlabel("Frequency")
plt.ylabel("Power")
plt.grid()
plt.show()

print("")
