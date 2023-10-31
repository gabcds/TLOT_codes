import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
def anti_aliasing_filter(x,fs,cutoff_freq,filter_order):
    """Perform a anti_aliasing_filter

    Args:
        x (array_like): signal to be filtered
        fs (int): sampling frequency
        cutoff_freq (int): set your desired cutoff frequency
        filter_order (int): desired filter order

    Returns:
        filtered_signal (ndarray): the output of the digital filter

    """
    #t = np.linspace(0, 1, fs, endpoint=False)
    nyquist = fs/2
    cutoff = cutoff_freq/nyquist
    b, a = signal.butter(filter_order,cutoff, btype = "lowpass")
    filtered_signal = signal.lfilter(b,a,x)
    return filtered_signal
    

