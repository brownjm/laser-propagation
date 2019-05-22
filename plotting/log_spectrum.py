"""Plot the log spectrum for a given input file and distances"""
import sys
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200
plt.close('all')

import load


def plot(input_file, z):
    r = load.Results(input_file)
    nm = r.wavelength * 1e9

    fig, ax = plt.subplots()
    z, spec = r.spectrum(z)
    spec /= spec.max()
    ax.plot(nm, spec)
        
    ax.set_xlim(nm.min(), nm.max())
    ax.set_yscale('log')
    ax.set_ylim(1e-10, 2)
    
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel('spectral intensity [scaled]')
    ax.set_title('z = {:1.2f}m'.format(z))
    
    fig.tight_layout()
        
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Please provide an input file and distance")
        sys.exit()
    
    input_file = sys.argv[1]
    z = float(sys.argv[2])
    plot(input_file, z)
