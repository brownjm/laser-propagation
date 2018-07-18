"""Plot the maximum electron density along propagation"""
import sys
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200
plt.close('all')

import load


def plot(input_file):
    r = load.Results(input_file)
    z, density = r.max_electron_density()

    fig, ax = plt.subplots()
    ax.plot(z, density)

    ax.set_xlabel('distance [m]')
    ax.set_ylabel('electron density [1/m^3]')
    ax.set_yscale('log')
    fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Please provide an input file")
        sys.exit()

    input_file = sys.argv[1]
    plot(input_file)
