import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, interpolate
Ac = .5
Am = 1
wm = 1500
wc = 4000
fs = 25000
initial_time = 0
final_time = 1
time = np.arange(initial_time, final_time, 1.0 / fs )
noise_amplitude = .5
noise = np.random.normal(0, noise_amplitude, len(time))
rpm = np.linspace(0, 100, int(len(time)/100))
rpm[:] = 4000 + np.random.normal(-5, 5, int(len(time)/100))
for i in range(0,len(rpm)):

    x = (Ac + Am*np.sin((wm/60)*time))*np.cos((rpm[i]/60)*time) +noise

