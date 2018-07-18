"""Plot the spectral field in frequency space for a given input file and distance"""
import sys
import numpy
from scipy.constants import pi, epsilon_0, c
from matplotlib import colors, cm
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200
plt.close('all')

import load


def plot(input_file, z):
    r = load.Results(input_file)

    fig, ax = plt.subplots()
    z, A = r.spectral_field(z)
    I = 0.5 * epsilon_0 * c * abs(A)**2
    I /= I.max()

    omega0 = 2*pi*c / r.config['laser/wavelength']
    order = r.omega / omega0
    
    if 'capillary/modes' in r.config:
        modes = r.config['capillary/modes']
        m = ax.imshow(I[::-1], norm=colors.LogNorm(1e-6, 1), aspect='auto',
                      extent=(order.min(), order.max(), -0.5, modes-0.5))
        ax.set_ylabel('modes')
        ax.set_yticks(numpy.arange(modes))

        
    else:
        kperp = numpy.hstack((0, r.kperp))
        m = ax.pcolormesh(order, kperp, I, norm=colors.LogNorm(1e-6, 1))
        ax.set_ylabel('kperp [1/m]')
        
    plt.colorbar(m, ax=ax)
    cmap = cm.get_cmap()
    cmap.set_bad(cmap(0))
    
    ax.set_xlabel('harmonic order')
    ax.set_title('Spectral intensity [scaled], z = {:1.2f}m'.format(z))

    
    fig.tight_layout()
        
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Please provide an input file and distance")
        sys.exit()
    
    input_file = sys.argv[1]
    z = float(sys.argv[2])
    plot(input_file, z)
