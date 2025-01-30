import numpy as np
import matplotlib.pyplot as plt


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
    fs = 1/dt
    # Calculate the power spectrum
    P = (np.abs(H) ** 2) / N

    # Generate frequencies using FFT frequencies
    freqs = np.fft.fftfreq(N, dt)  
    
    # Shift the frequencies and power spectrum for symmetric representation
    freqs_sym = np.fft.fftshift(freqs)
    P_sym = np.fft.fftshift(P)
    
    return freqs_sym, P_sym



# Generate data
data = np.loadtxt('positions.txt')
time = data[:,0]
positions = data[:,1:]

dt = time[1] - time[0]

# Asymmetric spectrum
for i in range(3):
    freqs, P = numpy_fft(positions[:,i], dt)    
    plt.subplot(1,3,i+1)
    plt.stem(freqs, P, basefmt=" ")
    plt.title("Power Spectrum")
    plt.xlabel("Frequency index")
    plt.ylabel("Power")
    plt.grid()
    plt.xlim([-100,100])
    Ps = np.sort(abs(P))[-2]
    print("Frequencies: %f", freqs[P > Ps*0.95])
plt.show()
print("")

plt.subplot(1,2,1)
for i in range(3):
    plt.plot(time[:5000], positions[:5000,i])
    plt.title("Motion in time")
    plt.xlabel("time [ps]")
    plt.ylabel("motion [Å]")
    plt.grid()


data = np.loadtxt('energies.txt')
data = data[:,1:]
plt.subplot(1,2,2)
for i in range(3):
    plt.plot(time[:5000], data[:5000,i])
    plt.title("Motion in time")
    plt.xlabel("time [ps]")
    plt.ylabel("motion [Å]")
    plt.grid()
plt.show()