"""Plot the results of simulating a single argon atom on a field file"""
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from matplotlib.ticker import MaxNLocator
import load


def load_binary(filename, nrows, ncols, dtype=complex):
    """Load data from a binary file"""
    temporal = np.fromfile(filename, dtype=dtype)
    temporal.shape = (nrows, ncols)
    return temporal

def plot_wavefunction(r, l, wf, title=None, n=None, logthresh=1e-10):
    """Plot the wavefunction in radius and momentum: psi(r, l)"""
    fig, ax = plt.subplots()
    if n == None:
        I = abs(wf[::-1])**2
        m = ax.imshow(I, norm=colors.LogNorm(I.max()*logthresh, I.max()), aspect='auto',
                      extent=(r[0], r[-1], -0.5, l[-1]+0.5))
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_xlabel('radius')
        ax.set_ylabel('angular momentum')
        plt.colorbar(m, ax=ax)
        cmap = cm.get_cmap()
        cmap.set_bad(cmap(0))
        
    else:
        ax.plot(r, abs(wf[n])**2)
        ax.set_xlabel('radius')

    if title != None:
        ax.set_title(title)
    fig.tight_layout()
    return fig



# load configuration
config = load.load_config('input')
folder = config['folder']

# calculate space and momentum coordinates
Nr = config['argon/nr']
r = np.array([(i+0.5) * 0.25 for i in range(Nr)], dtype=float)
Nl = config['argon/nl']
l = np.array([l for l in range(Nl)], dtype=int)

# read in wavefunction
wf = load_binary(os.path.join(folder, 'test_argon_wavefunction.dat'), len(l), len(r))

# plot wavefunction
plot_wavefunction(r, l, wf, title='wavefunction after laser pulse')

# load on-axis field, ionization rate, electron density, and nonlinear dipole
time, onaxis = np.loadtxt(os.path.join(folder, 'test_argon_onaxis_field.dat'), unpack=True)
time, rate = np.loadtxt(os.path.join(folder, 'test_argon_ionization_rate.dat'), unpack=True)
time, density = np.loadtxt(os.path.join(folder, 'test_argon_electron_density.dat'), unpack=True)
time, dipole = np.loadtxt(os.path.join(folder, 'test_argon_nonlinear_dipole.dat'), unpack=True)

# time in femtoseconds
fs = time * 1e15

# plot nonlinear dipole
fig, ax = plt.subplots()
ax.plot(fs, dipole, label='nonlinear dipole')
ax.plot(fs, onaxis/onaxis.max() * dipole.max(), color='gray', label='field (scaled)')
ax.set_xlabel('time [fs]')
ax.set_ylabel('nonlinear dipole')
ax.legend()
fig.tight_layout()

# plot ionization rate
fig, ax = plt.subplots()
ax.plot(fs, rate, label='ionization rate')
ax.plot(fs, abs(onaxis/onaxis.max() * rate.max()), color='gray', label='field (scaled)')
ax.set_xlabel('time [fs]')
ax.set_ylabel('ionization rate [1/s]')
ax.legend()
fig.tight_layout()

# plot electron density
fig, ax = plt.subplots()
ax.plot(fs, density, label='electron_density')
ax.plot(fs, abs(onaxis/onaxis.max() * density.max()), color='gray', label='field (scaled)')
ax.set_xlabel('time [fs]')
ax.set_ylabel('electron density [1/m^3]')
ax.legend()
fig.tight_layout()

plt.show()
