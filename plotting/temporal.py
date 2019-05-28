"""Plot the electric field in space/time for a given input file and distance"""
import argparse
import numpy
import matplotlib.pyplot as plt
import load

def plot(input_file, z):
    r = load.Results(input_file)
    z, E = r.electric_field(z)
    cm = z * 100
    fs = r.time * 1e15
    mm = r.radius * 1e3
    mm = numpy.hstack((0, mm))

    fig, ax = plt.subplots()
    v = max(E.real.max(), -E.real.min())
    mesh = ax.pcolormesh(fs, mm, E.real, vmin=-v, vmax=v, cmap='bwr')
        
    ax.set_xlim(fs.min(), fs.max())
    ax.set_ylim(mm.min(), mm.max())
    
    ax.set_xlabel('time [fs]')
    ax.set_ylabel('radius [mm]')
    ax.set_title('Electric field [V/m], z = {:1.2f}cm'.format(cm))

    plt.colorbar(mesh, ax=ax)
    fig.tight_layout()
        
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot the space-time electric field')
    parser.add_argument('input_file', help='input file from simulation')
    parser.add_argument('distance', type=float, help='distance along propagation')

    args = parser.parse_args()
    plot(args.input_file, args.distance)
