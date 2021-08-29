# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 08:17:12 2021

@author: felix
"""

#Radial plot

from t2listing import *
import matplotlib.pyplot as plt
from t2data import * #t2data is both a module and class.
import numpy as np
from mulgrids import *
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
import os
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import pandas as pd
from CoolProp.CoolProp import PropsSI
import matplotlib as mpl
import matplotlib.font_manager as font_manager
import matplotlib.pylab as pl
from matplotlib import rc

SMALL_SIZE = 12
BIGGER_SIZE = 13

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=11)    # legend fontsize

#Simulation parameters
model_number = '020'
years = 20
initial_temperature = 300
p_sat = 85.93 #(300, 85.93),(290, 74.45), (280, 64.17), (270, 55.03), (260, 46.92)
generator = 20 #kg/s
res_thickness = 500. #m

#Control radial dimension of the grid
first_block_radial = 0 # 10log
last_block_radial = 3 # 10log
number_innerblocks = 10
number_outerblocks = 100

#Dual porosity parameters
fracture_volume = 0.001 #[0.01, 0.09, 0.3, 0.6] 
fracture_permeability_minc = 7e-14
matrix_permeability_minc = 1e-18
matrix_permeability_sp = 3e-15
fracture_porosity = 0.9
matrix_porosity = 0.1
fracture_spacing = 50
fracture_planes = 1
fracture_extension = 50 #Fracture extension in number of blocks
rel_perm_minc = 3 # 3=Corey
rel_perm_sp = 1 # 1 = Xcurve
heat_conductivity_matrix = 2
heat_conductivity_fracture = 1.5

#Make changes to the dat files and run the files with TOUGH2 executable 
listing_files = []
    
dat = t2data('WELL1.dat')#Import dat file that will be modified    

dat.parameter['default_incons'][1] = initial_temperature #Change initial temperature. IT= 200C (single phase) 
dat.parameter['default_incons'][0] = p_sat*10**5 #Change initial pressure. IP= 42bar 
dat.parameter['tstart'] = 0.0           # assign start time
yr = 3600.*24.*365.25
dat.parameter['tstop'] = years*yr         # assign end time
dat.parameter['const_timestep'] = 24*3600.0*7.0     # assign initial time step to 1 day
dat.parameter['max_timestep'] = 24*3600.0*7.0     # assign max time step to 1 week
dat.parameter['print_interval'] = 10     # assign max time step to 1 week
dat.parameter['max_timesteps'] = 9999
dat.parameter['max_duration'] = 0
dat.parameter['option'][24] = 2
dat.output_times['time'] = [3*24*3600, 30*24*3600, 180*24*3600, 5*365*24*3600, 10*365*24*3600]
dat.short_output['frequency'] = 1
dat.short_output['generator'] = []
dat.short_output['block'] = []

# Modify radial thickness of grid
dr2_1 = [10**first_block_radial]*number_innerblocks
dr2_2 = np.logspace(first_block_radial, last_block_radial, number_outerblocks)
dr2_2 = dr2_2.tolist()
radial_thickness = dr2_1 + dr2_2
res_extension = sum(radial_thickness)
dat.grid = t2grid().radial(radial_thickness, [res_thickness])

#Remove and add new generator
dat.clear_generators()
gen = t2generator(name='wel 1', block='  a 1', type='MASS') #, gx=q*col.area
dat.add_generator(gen)
dat.generator['  a 1', 'wel 1'].gx = generator*-1
       
dat.relative_permeability['type'] = rel_perm_minc
dat.relative_permeability['parameters'] = [0, 0]

#Function that changes the names of the rocktypes after applying minc
def minc_rock(rockname, level):
    first = 'f' if level == 0 else 'm'
    return first + rockname[1:]  

indices = dat.grid.minc([fracture_volume, 1-fracture_volume], spacing = fracture_spacing,
num_fracture_planes = fracture_planes, minc_rockname = minc_rock)
dat.grid.rocktype['ffalt'].permeability= [fracture_permeability_minc]*3 #fracture
dat.grid.rocktype['ffalt'].porosity= fracture_porosity #fracture
dat.grid.rocktype['mfalt'].permeability= [matrix_permeability_minc]*3 #matrix
dat.grid.rocktype['mfalt'].porosity= matrix_porosity #matrix
dat.grid.rocktype['mfalt'].conductivity= heat_conductivity_matrix 
dat.grid.rocktype['ffalt'].conductivity= heat_conductivity_fracture 
i=1
dat.write('dual_porosity'+str(i)+'.dat')
dat.run(incon_filename= 'dual_porosity'+str(i)+'.dat', simulator='../tim-0.7.1-251-teaching/AUTOUGH2_5_teaching') 
listing_files.append('dual_porosity'+str(i))

# Create a figure with multiple plots organized in a grid.
fig = plt.figure(figsize=(9, 4))
gs1 = gridspec.GridSpec(1, 1)
ax1 = fig.add_subplot(gs1[0])
ax2 = ax1.twinx()

for i in listing_files:
    lst = t2listing(i+'.LISTING') #data+'.LISTING'
    time2 , enthalpy_gen = lst.history(('g', ('  a 1', 'wel 1'), 'Enthalpy' ))
    time2 , pressure_a1 = lst.history(('e', ('  a 1'), 'Pressure' ))
        
    if lst.times[-1] >= (years-years*0.05)*3600.*24.*365.25:

        ax1.plot(time2/(365*3600*24), pressure_a1*10**-5, 'r', linestyle= '-') #label = labels[z])
        ax2.plot(time2/(365*3600*24), enthalpy_gen*10**-3, 'b', linestyle= '-')# label = labels[z])
    else: 
        label_texts[k] = label_texts[k]+' N/A'
        
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Fluid pressure (bar)', color='r')
    ax2.set_ylabel('Flowing enthalpy (kJ/kg)', rotation=270, labelpad=16, color='b')
    ax1.tick_params(axis='y', colors='r')
    ax2.tick_params(axis='y', colors='b')

#
#
# SP model (using Xcurve instead of Corey curve)
#
#

#Make changes to the dat files and run the files with TOUGH2 executable 
listing_files = []
    
dat = t2data('WELL1.dat')#Import dat file that will be modified    

dat.parameter['default_incons'][1] = initial_temperature #Change initial temperature. IT= 200C (single phase) 
dat.parameter['default_incons'][0] = p_sat*10**5 #Change initial pressure. IP= 42bar 
dat.parameter['tstart'] = 0.0           # assign start time
yr = 3600.*24.*365.25
dat.parameter['tstop'] = years*yr         # assign end time
dat.parameter['const_timestep'] = 24*3600.0*7.0     # assign initial time step to 1 day
dat.parameter['max_timestep'] = 24*3600.0*7.0     # assign max time step to 1 week
dat.parameter['print_interval'] = 10     # assign max time step to 1 week
dat.parameter['max_timesteps'] = 9999
dat.parameter['max_duration'] = 0
dat.parameter['option'][24] = 2
dat.output_times['time'] = [3*24*3600, 30*24*3600, 180*24*3600, 5*365*24*3600, 10*365*24*3600]
dat.short_output['frequency'] = 1
dat.short_output['generator'] = []
dat.short_output['block'] = []

# Modify radial thickness of grid
dr2_1 = [10**first_block_radial]*number_innerblocks
dr2_2 = np.logspace(first_block_radial, last_block_radial, number_outerblocks)
dr2_2 = dr2_2.tolist()
radial_thickness = dr2_1 + dr2_2
res_extension = sum(radial_thickness)
dat.grid = t2grid().radial(radial_thickness, [res_thickness])  

#Remove and add new generator
dat.clear_generators()
gen = t2generator(name='wel 1', block='  a 1', type='MASS') #, gx=q*col.area
dat.add_generator(gen)
dat.generator['  a 1', 'wel 1'].gx = generator*-1      

#Adjust the relative permeability
dat.relative_permeability['type'] = rel_perm_sp
dat.relative_permeability['parameters'] = [0, 0, 1, 1]

i=1
dat.grid.rocktype['dfalt'].permeability= [matrix_permeability_sp]*3 #matrix
dat.grid.rocktype['dfalt'].porosity= matrix_porosity #matrix
dat.grid.rocktype['dfalt'].conductivity= heat_conductivity_matrix 
dat.write('single_porosity'+str(i)+'.dat')
dat.run(incon_filename= 'single_porosity'+str(i)+'.dat', simulator='../tim-0.7.1-251-teaching/AUTOUGH2_5_teaching') 
listing_files.append('single_porosity'+str(i))
label_texts.append('single_porosity'+str(i))

k = 0
for i in listing_files:
    lst = t2listing(i+'.LISTING') #data+'.LISTING'    
    time2 , enthalpy_gen = lst.history(('g', ('  a 1', 'wel 1'), 'Enthalpy' ))
    time2 , pressure_a1 = lst.history(('e', ('  a 1'), 'Pressure' ))
                
    if lst.times[-1] >= (years-years*0.05)*3600.*24.*365.25:
        ax1.plot(time2/(365*3600*24), pressure_a1*10**-5, 'r', linestyle= '--') #label = labels[z])
        ax2.plot(time2/(365*3600*24), enthalpy_gen*10**-3, 'b', linestyle= '--')# label = labels[z])
    
    else: 
        label_texts[k] = label_texts[k]+' N/A'
        
    labels = ['MINC model *','SP model']
    linestyle = ['-', '--']
    legend_elements = [Line2D([0], [0], color='k', linestyle= linestyle[0], label = labels[0]),
                    Line2D([0], [0], color='k', linestyle = linestyle[1], label = labels[1]),]
    ax2.legend(handles = legend_elements, loc=('center right'), framealpha=1)  
    plt.tight_layout()

    fig.savefig('MINC/SP_vs_MINC_timehistory.jpg', dpi = 150) 
    plt.show()
