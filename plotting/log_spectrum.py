"""Plot the log spectrum for a given input file and distances"""
import argparse
import matplotlib.pyplot as plt
import load

def plot(input_file, z):
    r = load.Results(input_file)
    z, spec = r.spectrum(z)
    spec /= spec.max()
    nm = r.wavelength * 1e9
    cm = z * 100

    fig, ax = plt.subplots()
    ax.plot(nm, spec)
        
    ax.set_xlim(nm.min(), nm.max())
    ax.set_yscale('log')
    ax.set_ylim(1e-10, 2)
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel('spectral intensity [scaled]')
    ax.set_title('z = {:1.2f}cm'.format(cm))
    fig.tight_layout()
        
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot the spectrum of the laser in log scale')
    parser.add_argument('input_file', help='input file from simulation')
    parser.add_argument('distance', type=float, help='distance along propagation')

    args = parser.parse_args()
    plot(args.input_file, args.distance)
