"""Plot the electron density in space/time for a given input file and distance"""
import argparse
import numpy
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import load

def plot(input_file, z):
    r = load.Results(input_file)
    z, rho = r.electron_density(z)
    fs = r.time * 1e15
    mm = r.radius * 1e3
    cm = z * 100
    
    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(fs, mm, rho, cmap='magma',
                         norm=LogNorm(vmin=rho.max()/1e6, vmax=rho.max()))
        
    ax.set_xlim(fs.min(), fs.max())
    ax.set_ylim(mm.min(), mm.max())
    ax.set_xlabel('time [fs]')
    ax.set_ylabel('radius [mm]')
    ax.set_title('Electron density [1/m^3], z = {:1.2f}cm'.format(cm))
    plt.colorbar(mesh, ax=ax)
    fig.tight_layout()
        
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot the space-time electron density')
    parser.add_argument('input_file', help='input file from simulation')
    parser.add_argument('distance', type=float, help='distance along propagation')

    args = parser.parse_args()
    plot(args.input_file, args.distance)
