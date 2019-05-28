"""Plot the maximum intensity along propagation"""
import argparse
import matplotlib.pyplot as plt
import load

def plot(input_file):
    r = load.Results(input_file)
    z, intensity = r.max_intensity()
    cm = z * 100

    fig, ax = plt.subplots()
    ax.plot(cm, intensity)

    ax.set_xlabel('distance [cm]')
    ax.set_ylabel('intensity [W/m^2]')
    ax.set_yscale('log')
    
    fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot the maximum laser intensity along propagation')
    parser.add_argument('input_file', help='input file from simulation')

    args = parser.parse_args()
    plot(args.input_file)
