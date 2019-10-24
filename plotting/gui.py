import sys
import os
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
import subprocess
import numpy as np

class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.grid()
        self.setup_menu()
        self.top = tk.Frame(self)
        self.top.grid(row=0, column=0)
        self.bottom = tk.Frame(self)
        self.bottom.grid(row=1, column=0)
        self.setup_distances()
        self.setup_input_file()
        self.create_widgets()

    def setup_menu(self):
        self.menubar = tk.Menu(self.master)
        self.menubar.add_command(label="Open", command=self.open_filename)
        self.master.config(menu=self.menubar)
        
    def setup_distances(self):
        self.distance_text = ttk.Label(self.top, text='Distance [m]: ')
        self.distance_text.grid(row=1, column=0)
            
        self.distances = ttk.Combobox(self.top)
        self.distances.grid(row=1, column=1, columnspan=2)
        
    def setup_input_file(self):
        self.input_file_text = ttk.Label(self.top, text='Input file: ')
        self.input_file_text.grid(row=0, column=0)

        self.input_file = tk.StringVar()
        self.input_file_entry = ttk.Entry(self.top, textvariable=self.input_file)
        self.input_file_entry.grid(row=0, column=1, columnspan=2)

        # self.input_file_button = ttk.Button(self.top, text='Refresh')
        # self.input_file_button['command'] = self.refresh_directory
        # self.input_file_button.grid(row=0, column=2)

    def open_filename(self):
        input_file = filedialog.askopenfilename(initialdir=self.input_file.get(),
                                                title='Please select an input file')
        if self.input_file:
            self.input_file.set(input_file)
            self.load_distances()
        
    def load_distances(self):
        filename = os.path.join(os.path.dirname(self.input_file.get()), 'distance.dat') 
        with open(filename, 'r') as f:
            distances = ['{:1.6f}'.format(float(z)) for z in f.readlines()]

        self.distances.config(values=distances)
        self.distances.current(0)

    def create_widgets(self):
        # first column
        self.temporal = ttk.Button(self.bottom, text='E(r,t)')
        self.temporal['command'] = lambda: self.launch_script('temporal.py',
                                                              self.input_file.get(),
                                                              self.distances.get())
        self.temporal.grid(row=3, column=0)
        
        self.onaxis = ttk.Button(self.bottom, text='E(r=0,t)')
        self.onaxis['command'] = lambda: self.launch_script('onaxis.py',
                                                            self.input_file.get(),
                                                            self.distances.get())
        self.onaxis.grid(row=4, column=0)

        self.intensity = ttk.Button(self.bottom, text='I(r,t)')
        self.intensity['command'] = lambda: self.launch_script('intensity.py',
                                                               self.input_file.get(),
                                                               self.distances.get())
        self.intensity.grid(row=5, column=0)
        
        self.max_intensity = ttk.Button(self.bottom, text='Max I(z)')
        self.max_intensity['command'] = lambda: self.launch_script('max_intensity.py',
                                                                   self.input_file.get())
        self.max_intensity.grid(row=6, column=0)

        self.energy = ttk.Button(self.bottom, text='Energy(z)')
        self.energy['command'] = lambda: self.launch_script('energy.py',
                                                               self.input_file.get())
        self.energy.grid(row=7, column=0)


        # second column
        self.spectral = ttk.Button(self.bottom, text='A(k,\u03c9)')
        self.spectral['command'] = lambda: self.launch_script('spectral.py',
                                                              self.input_file.get(),
                                                              self.distances.get())
        self.spectral.grid(row=3, column=1)
        
        self.spectrum = ttk.Button(self.bottom, text='S(\u03bb)')
        self.spectrum['command'] = lambda: self.launch_script('spectrum.py',
                                                              self.input_file.get(),
                                                              self.distances.get())
        self.spectrum.grid(row=4, column=1)

        self.log_spectrum = ttk.Button(self.bottom, text='log S(\u03bb)')
        self.log_spectrum['command'] = lambda: self.launch_script('log_spectrum.py',
                                                                  self.input_file.get(),
                                                                  self.distances.get())
        self.log_spectrum.grid(row=5, column=1)

        self.thz_spectrum = ttk.Button(self.bottom, text='log S(THz)')
        self.thz_spectrum['command'] = lambda: self.launch_script('thz_spectrum.py',
                                                                  self.input_file.get(),
                                                                  self.distances.get())
        self.thz_spectrum.grid(row=6, column=1)


        # third column
        self.density = ttk.Button(self.bottom, text='\u03c1(r,t)')
        self.density['command'] = lambda: self.launch_script('density.py',
                                                             self.input_file.get(),
                                                             self.distances.get())
        self.density.grid(row=3, column=2)

        self.max_density = ttk.Button(self.bottom, text='Max \u03c1(z)')
        self.max_density['command'] = lambda: self.launch_script('max_density.py',
                                                                 self.input_file.get())
        self.max_density.grid(row=4, column=2)


    def launch_script(self, name, *args):
        if not self.input_file.get():
            print("No input file found")
            return
        
        script_name = os.path.join(sys.path[0], name)
        subprocess.Popen(['python3', script_name, *args],
                         stdin=None, stdout=None, stderr=None, close_fds=True)
        

root = tk.Tk()
root.title('Plotting for laser-propagation')
#root.option_add("*Font", ("TkDefaultFont", 12))
app = Application(master=root)
app.mainloop()
