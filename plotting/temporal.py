"""Plot the electric field in space/time for a given input file and distance"""
import sys
import numpy
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
    z, E = r.electric_field(z)
    v = max(E.real.max(), -E.real.min())
    mesh = ax.pcolormesh(fs, mm, E.real, vmin=-v, vmax=v, cmap='bwr')
        
    ax.set_xlim(fs.min(), fs.max())
    ax.set_ylim(mm.min(), mm.max())
    
    ax.set_xlabel('time [fs]')
    ax.set_ylabel('radius [mm]')
    ax.set_title('Electric field [V/m], z = {:1.2f}m'.format(z))

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
