"""Plot the on-axis electric field for a given input file and distance"""
import argparse
import matplotlib.pyplot as plt
import load

def plot(input_file, z):
    r = load.Results(input_file)
    fs = r.time * 1e15
    z, E = r.electric_field(z)
    cm = z * 100

    fig, ax = plt.subplots()
    ax.plot(fs, E[0].real)
        
    ax.set_xlim(fs.min(), fs.max())
    ax.set_xlabel('time [fs]')
    ax.set_ylabel('electric field [V/m]')
    ax.set_title('z = {:1.2f}cm'.format(cm))
    fig.tight_layout()
        
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot the on-axis electric field')
    parser.add_argument('input_file', help='input file from simulation')
    parser.add_argument('distance', type=float, help='distance along propagation')

    args = parser.parse_args()
    plot(args.input_file, args.distance)
