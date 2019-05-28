"""Plot the laser pulse energy along propagation"""
import argparse
import matplotlib.pyplot as plt
import load

def plot(input_file):
    r = load.Results(input_file)
    z, energy = r.energy()
    cm = z * 100
    energy_mJ = energy * 1e3
    
    fig, ax = plt.subplots()
    ax.plot(cm, energy_mJ)
    
    ax.set_xlabel('distance [cm]')
    ax.set_ylabel('energy [mJ]')
    fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot the laser pulse energy along propagation')
    parser.add_argument('input_file', help='input file from simulation')

    args = parser.parse_args()
    plot(args.input_file)
