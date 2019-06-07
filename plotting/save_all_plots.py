"""Save plots for all of the output data for a given input file"""
import sys
from os.path import join
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import load
import temporal
import spectral
import density
import intensity
import onaxis
import spectrum
import log_spectrum
import thz_spectrum
import energy
import max_intensity
import max_density


if len(sys.argv) < 2:
    print("Please provide an input file")
    sys.exit()


input_file = sys.argv[1]
r = load.Results(input_file)


results_files = ['temporal_field', 'spectral_field', 'electron_density']
results_files = [join('results', o) for o in results_files]
funcs = [temporal.plot, spectral.plot, density.plot]

for func, results in zip(funcs, results_files):
    for i, z in enumerate(r.distances):
        func(input_file, z)
        filename = r._make_filename(r.config[results], i).replace('dat', 'png')
        plt.savefig(filename)
        print('Wrote:', filename)
        plt.close('all')


for i, z in enumerate(r.distances):
    intensity.plot(input_file, z)
    filename = 'I{:03d}.png'.format(i)
    plt.savefig(filename)
    print('Wrote:', filename)
    plt.close('all')

    
for i, z in enumerate(r.distances):
    onaxis.plot(input_file, z)
    filename = 'onaxis{:03d}.png'.format(i)
    plt.savefig(filename)
    print('Wrote:', filename)
    plt.close('all')


for i, z in enumerate(r.distances):
    spectrum.plot(input_file, z)
    filename = 'spectrum{:03d}.png'.format(i)
    plt.savefig(filename)
    print('Wrote:', filename)
    plt.close('all')


for i, z in enumerate(r.distances):
    log_spectrum.plot(input_file, z)
    filename = 'log_spectrum{:03d}.png'.format(i)
    plt.savefig(filename)
    print('Wrote:', filename)
    plt.close('all')

for i, z in enumerate(r.distances):
    thz_spectrum.plot(input_file, z)
    filename = 'thz_spectrum{:03d}.png'.format(i)
    plt.savefig(filename)
    print('Wrote:', filename)
    plt.close('all')
    
energy.plot(input_file)
filename = 'energy.png'
plt.savefig(filename)
print('Wrote:', filename)
plt.close('all')

max_intensity.plot(input_file)
filename = 'max_intensity.png'
plt.savefig(filename)
print('Wrote:', filename)
plt.close('all')

max_density.plot(input_file)
filename = 'max_density.png'
plt.savefig(filename)
print('Wrote:', filename)
plt.close('all')
