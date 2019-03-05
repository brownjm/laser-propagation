"""Plot the electron density in space/time for a given input file and distance"""
import sys
import numpy
from scipy.constants import epsilon_0, c
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200
plt.close('all')

import load


def plot(input_file, z):
    r = load.Results(input_file)
    fs = r.time * 1e15
    mm = r.radius * 1e3
    mm = numpy.hstack((0, mm))
    fig, ax = plt.subplots()
    z, rho = r.electron_density(z)
    mesh = ax.pcolormesh(fs, mm, rho, cmap='magma', norm=LogNorm(vmin=rho.max()/1e6, vmax=rho.max()))
        
    ax.set_xlim(fs.min(), fs.max())
    ax.set_ylim(mm.min(), mm.max())
    
    ax.set_xlabel('time [fs]')
    ax.set_ylabel('radius [mm]')
    ax.set_title('Electron density [1/m^3], z = {:1.2f}m'.format(z))

    plt.colorbar(mesh, ax=ax)
    
    # fig.tight_layout()
        
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Please provide an input file and distance")
        sys.exit()
    
    input_file = sys.argv[1]
    z = float(sys.argv[2])
    plot(input_file, z)
