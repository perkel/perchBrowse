
import matplotlib.pyplot as plt
from scipy import signal

from numpy import load
from pylab import plot, show, title, xlabel, ylabel, subplot
from scipy import fft, arange
import matplotlib

## Makes spectrogram from perch torque measurements

def plotSpectrum(y,Fs):
    """
    Plots a Single-Sided Amplitude Spectrum of y(t)
    """
    power, freqs = matplotlib.mlab.psd(y, Fs=Fs, NFFT = int(2**8), pad_to=30000)
    print power[0:5], freqs[0:5]
    
    plot(freqs,power,'r') # plotting the spectrum
    xlabel('Freq (Hz)')
    ylabel('Power')
    show()

data = load("Perching_data_2016-07-01Trial_000.npy")
y = data[1]
print y[0:5]
print data[0][1]
Fs = 1.0/(data[0][1] - data[0][0])
plotSpectrum(y, Fs)



freqs, bins, im = signal.spectrogram(y, Fs)
print bins[0:5]
plt.pcolormesh(bins, freqs, im)
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.colorbar()
plt.show()

print "max im", im.max()