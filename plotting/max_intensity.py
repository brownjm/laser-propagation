"""Plot the maximum intensity along propagation"""
import sys
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200
plt.close('all')

import load


def plot(input_file):
    r = load.Results(input_file)
    z, intensity = r.max_intensity()

    fig, ax = plt.subplots()
    ax.plot(z, intensity)

    ax.set_xlabel('distance [m]')
    ax.set_ylabel('intensity [W/m^2]')
    ax.set_yscale('log')
    
    fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Please provide an input file")
        sys.exit()

    input_file = sys.argv[1]
    plot(input_file)
