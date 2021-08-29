# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 16:52:40 2021

@author: felix
"""

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
matrix_permeability_minc = 1e-18#3e-18
# matrix_permeability_sp = 1e-13#3e-18
# matrix_permeabilities_sp = [3e-15]#3e-18

fracture_porosity = 0.9
matrix_porosity = 0.1
fracture_spacing = 50
fracture_planes = 1
heat_conductivity_matrix = 2
heat_conductivity_fracture = 1.5
    
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
 
# Get center and connection coordinates for the blocks. Used for plotting.
cone_coord = []
cent_coord = []
j = 0    
for k in radial_thickness:    
    cent_coord.append(sum(radial_thickness[0:j])+k/2)
    cone_coord.append(sum(radial_thickness[0:j])+k)
    j +=1
cone_coord = cone_coord[:-1]    

#Remove and add new generator
dat.clear_generators()
gen = t2generator(name='wel 1', block='  a 1', type='MASS') #, gx=q*col.area
dat.add_generator(gen)
dat.generator['  a 1', 'wel 1'].gx = generator*-1     

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


#Simulate different relative permeability curves
i=1
rel_perm_function = 1 # Xcurve
dat.relative_permeability['type'] = rel_perm_function 
dat.relative_permeability['parameters'] = [0, 0, 1, 1]
dat.write('dual_porosity'+str(i)+'.dat')
dat.run(incon_filename= 'dual_porosity'+str(i)+'.dat', simulator='../tim-0.7.1-251-teaching/AUTOUGH2_5_teaching') 
listing_files = []
listing_files.append('dual_porosity'+str(i))
i+=1

rel_perm_function = 3 # 3=Corey      
dat.relative_permeability['type'] = rel_perm_function 
dat.relative_permeability['parameters'] = [0, 0]
dat.write('dual_porosity'+str(i)+'.dat')
dat.run(incon_filename= 'dual_porosity'+str(i)+'.dat', simulator='../tim-0.7.1-251-teaching/AUTOUGH2_5_teaching') 
listing_files.append('dual_porosity'+str(i))


# Create a figure with multiple plots organized in a grid.
fig = plt.figure(figsize=(9, 4))
gs1 = gridspec.GridSpec(1, 1)
ax1 = fig.add_subplot(gs1[0])
ax2 = ax1.twinx()

linestyle = ['-','--'] 

k = 0
for i in listing_files:
    # Open listing file. Make plots with the listing file
    lst = t2listing(i+'.LISTING') #data+'.LISTING'
    time2 , enthalpy_gen = lst.history(('g', ('  a 1', 'wel 1'), 'Enthalpy' ))
    time2 , pressure_a1 = lst.history(('e', ('  a 1'), 'Pressure' ))
                
    if lst.times[-1] >= (years-years*0.05)*3600.*24.*365.25:
        
        ax1.plot(time2/(365*3600*24), pressure_a1*10**-5, 'r', linestyle= linestyle[k]) 
        ax2.plot(time2/(365*3600*24), enthalpy_gen*10**-3, 'b', linestyle= linestyle[k])
    
    else: 
        label_texts[k] = label_texts[k]+' N/A'
        
    k+=1


    ax1.set_xlabel('Year')
    ax1.set_ylabel('Fluid pressure (bar)', color='r')
    ax2.set_ylabel('Flowing enthalpy (kJ/kg)', rotation=270, labelpad=16, color='b')
    ax1.tick_params(axis='y', colors='r')
    ax2.tick_params(axis='y', colors='b')      
    labels = ['X-curves', 'Corey curves *']
    
    legend_elements = [Line2D([0], [0], color='k', linestyle= linestyle[0], label = labels[0]),
                    Line2D([0], [0], color='k', linestyle = linestyle[1], label = labels[1])]
    
    ax2.legend(handles = legend_elements, loc=('lower right'), framealpha=1)  
    plt.tight_layout()
    
    fig.savefig('MINC/MINC_timehistory_XvsC.jpg', dpi = 150) 
    plt.show()
