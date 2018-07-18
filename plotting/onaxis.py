"""Plot the on-axis electric field for a given input file and distance"""
import sys
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200
plt.close('all')

import load


def plot(input_file, z):
    r = load.Results(input_file)
    fs = r.time * 1e15

    fig, ax = plt.subplots()
    z, E = r.electric_field(z)
    ax.plot(fs, E[0].real)
        
    ax.set_xlim(fs.min(), fs.max())
    
    ax.set_xlabel('time [fs]')
    ax.set_ylabel('electric field [V/m]')
    
    fig.tight_layout()
        
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Please provide an input file and distance")
        sys.exit()
    
    input_file = sys.argv[1]
    distances = float(sys.argv[2])
    plot(input_file, distances)
