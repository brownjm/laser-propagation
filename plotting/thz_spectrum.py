"""Plot the spectrum for a given input file and distance"""
import sys
from scipy.constants import pi
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200
plt.close('all')

import load


def plot(input_file, z):
    r = load.Results(input_file)
    nm = r.wavelength * 1e9
    thz = r.omega / (2*pi) / 1e12

    fig, ax = plt.subplots()
    z, spec = r.spectrum(z)
    spec /= spec.max()
    ax.plot(thz, spec)
        
    ax.set_xlim(0, thz.max())
    # ax.set_ylim(0, 1.1)
    ax.set_yscale('log')
    
    ax.set_xlabel('frequency [THz]')
    ax.set_ylabel('spectral intensity [scaled]')
    
    fig.tight_layout()
        
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Please provide an input file and distance")
        sys.exit()
    
    input_file = sys.argv[1]
    z = float(sys.argv[2])
    plot(input_file, z)
