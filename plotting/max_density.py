"""Plot the maximum electron density along propagation"""
import argparse
import matplotlib.pyplot as plt
import load

def plot(input_file):
    r = load.Results(input_file)
    z, density = r.max_electron_density()
    cm = z * 100

    fig, ax = plt.subplots()
    ax.plot(cm, density)

    ax.set_xlabel('distance [cm]')
    ax.set_ylabel('electron density [1/m^3]')
    ax.set_yscale('log')
    fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot the maximum electron density along propagation')
    parser.add_argument('input_file', help='input file from simulation')

    args = parser.parse_args()
    plot(args.input_file)
