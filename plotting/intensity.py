"""Plot the intensity in space/time for a given input file and distance"""
import sys
from scipy.constants import epsilon_0, c
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200
plt.close('all')

import load


def plot(input_file, z):
    r = load.Results(input_file)
    fs = r.time * 1e15
    mm = r.radius * 1e3

    fig, ax = plt.subplots()
    z, E = r.electric_field(z)
    I = 0.5 * epsilon_0 * c * abs(E)**2
    mesh = ax.pcolormesh(fs, mm, I, cmap='inferno')
        
    ax.set_xlim(fs.min(), fs.max())
    ax.set_ylim(mm.min(), mm.max())
    
    ax.set_xlabel('time [fs]')
    ax.set_ylabel('radius [mm]')
    ax.set_title('Intensity [W/m^2], z = {:1.2f}m'.format(z))

    plt.colorbar(mesh, ax=ax)
    
    fig.tight_layout()
        
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Please provide an input file and distance")
        sys.exit()
    
    input_file = sys.argv[1]
    z = float(sys.argv[2])
    plot(input_file, z)
