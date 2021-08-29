# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 09:36:11 2021

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
fracture_permeability = 7e-14
matrix_permeability = 1e-18#3e-18
fracture_porosity = 0.9
matrix_porosity = 0.1
fracture_spacing = 50
fracture_planes = 1
fracture_extension = 50 #Fracture extension in number of blocks
rel_perm_function = 3 # 3=Corey
residual_liquid_saturation = 0 #[Srl, Srg]
residual_steam_saturation = 0
heat_conductivity_matrix = 2
heat_conductivity_fracture = 1.5

parameter_to_plot = 'matrix_permeability' # matrix_permeability, fracture_permeability, fracture_porosity, fracture spacing etc.
parameter_values = [1e-19, 5e-19, 1e-18, 5e-18]

# parameter_to_plot = 'fracture_spacing'
# parameter_values = [50, 100, 150, 200]


#Exit script if parameter name is misspelled 
if parameter_to_plot not in globals():
    sys.exit('Parameter name invalid')

# Get units for legend
if parameter_to_plot == 'matrix_permeability':
    parameter2 = ['$k_{m}$', '$m^{2}$']     
if parameter_to_plot == 'fracture_permeability':
    parameter2 = ['$k_{f}$', '$m^{2}$']  
if parameter_to_plot == 'fracture_porosity':
    parameter2 = [r'$\phi_{f}$', '']
if parameter_to_plot == 'matrix_porosity':
    parameter2 = [r'$\phi_{m}$', '']  
if parameter_to_plot == 'fracture_spacing':
    parameter2 = ['Fracture spacing', 'm']   
    

#Make changes to the dat files and run the files with TOUGH2 executable 
listing_files = []
label_texts = []
i=1   
for parameter in parameter_values:
    
    vars()[parameter_to_plot] = parameter #turn parameter_to_plot string into variable
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
    dat.output_times['time'] = [60*24*3600, 180*24*3600, 1*365*24*3600, 5*365*24*3600, 20*365*24*3600]
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
    
    #Adjust the relative permeability
    dat.relative_permeability['type'] = rel_perm_function 
    dat.relative_permeability['parameters'] = [residual_liquid_saturation, residual_steam_saturation]
    
    #Function that chanes the names of the rocktypes after applying minc
    def minc_rock(rockname, level):
        first = 'f' if level == 0 else 'm'
        return first + rockname[1:]  
    
    indices = dat.grid.minc([fracture_volume, 1-fracture_volume], spacing = fracture_spacing,
    num_fracture_planes = fracture_planes, minc_rockname = minc_rock)
    dat.grid.rocktype['ffalt'].permeability= [fracture_permeability]*3 #fracture
    dat.grid.rocktype['ffalt'].porosity= fracture_porosity #fracture
    dat.grid.rocktype['mfalt'].permeability= [matrix_permeability]*3 #matrix
    dat.grid.rocktype['mfalt'].porosity= matrix_porosity #matrix
    dat.grid.rocktype['mfalt'].conductivity= heat_conductivity_matrix 
    dat.grid.rocktype['ffalt'].conductivity= heat_conductivity_fracture 
    dat.write('dual_porosity'+str(i)+'.dat')
    dat.run(incon_filename= 'dual_porosity'+str(i)+'.dat', simulator='../tim-0.7.1-251-teaching/AUTOUGH2_5_teaching') 
    listing_files.append('dual_porosity'+str(i))
    label_texts.append(parameter2[0]+' = '+str(parameter)+parameter2[1])
    i+=1

# Create a figure with multiple plots organized in a grid.
fig = plt.figure(figsize=(9, 4))

gs1 = gridspec.GridSpec(1, 1)
ax1 = fig.add_subplot(gs1[0])
ax3 = ax1.twinx()


color = ['r','b','g','y']
linestyle = ['-', '--']
time = years*365*24*3600


total_mass_list = []
total_heat_list = []

linestyle = ['-', '--', '-.', ':']

k = 0
for i, j in zip(listing_files, color):
    # Open listing file. Make plots with the listing file
    lst = t2listing(i+'.LISTING') #data+'.LISTING'
    lst.time = time
    # [60*24*3600, 180*24*3600, 365*24*3600, 5*365*24*3600]
    
    pressure_cent = lst.element['Pressure'][0:len(cent_coord)]
    mass_flow = lst.connection['Mass flow'][0:len(cone_coord)] 
    temperature = lst.element['Temperature'][0:len(cent_coord)]
    vapor_saturation = lst.element['Vapour saturation'][0:len(cent_coord)]
    enthalpy_flowing = lst.connection['Enthalpy'][0:len(cone_coord)] 
    vapor_mass_fraction = lst.connection['Vapour mass flow'][0:len(cone_coord)] /(lst.connection['Vapour mass flow'][0:len(cone_coord)] +lst.connection['Liquid mass flow'][0:len(cone_coord)]) 
    liquid_mass_fraction = lst.connection['Liquid mass flow'][0:len(cone_coord)]/(lst.connection['Vapour mass flow'][0:len(cone_coord)] +lst.connection['Liquid mass flow'][0:len(cone_coord)]) 
     
    
    pressure_cone = np.interp(cone_coord, cent_coord, pressure_cent)
    sw_flowing_cone = []
    for x, p in zip(vapor_mass_fraction, pressure_cone):        
            #Calculate specific volumes
        try:
            v_w = PropsSI('D','P', p ,'Q', 0 ,'IF97::Water')**-1 #100% saturated liquid water
            v_s = PropsSI('D','P', p ,'Q', 1 ,'IF97::Water')**-1 #100% saturated steam
        
            steam_saturation_cone = x*v_s/(v_w+x*(v_s-v_w))
            sw_flowing_cone_i = 1-steam_saturation_cone
            sw_flowing_cone.append(sw_flowing_cone_i)  
        except ValueError:
            pass #If pressure drops to zero it will be denoted N/A in plot
        
        
    sw_static = []
    for x in vapor_saturation:
        sw_static.append(1-x)
     
    
    if lst.times[-1] >= (years-years*0.05)*3600.*24.*365.25:
        
        ax1.semilogx(cent_coord, pressure_cent*10**-5, 'r', linestyle= linestyle[k], label = j) 
        ax3.semilogx(cone_coord, enthalpy_flowing*10**-3, 'b', linestyle= linestyle[k], label = 'DP(fracture)')

    
    else: 
        # ax2.plot([], [], ' ', label =' N/A')
        label_texts[k] = label_texts[k]+' N/A'
    
    k+=1


parameter_to_plot = parameter_to_plot.replace('_', ' ')  
  

ax1.set_ylabel('Fluid pressure (bar)', color='r')
ax3.set_ylabel('Flowing enthalpy (kJ/kg)', rotation=270, labelpad=16, color='b')

ax1.set_xlabel('Radial distance (m)')
ax1.set_xlim([1, 8000])

ax1.tick_params(axis='y', colors='r')
ax3.tick_params(axis='y', colors='b')


legend_elements = [Line2D([0], [0], color='k', linestyle= linestyle[0], label = label_texts[0]),
                Line2D([0], [0], color='k', linestyle = linestyle[1], label = label_texts[1]),
                Line2D([0], [0], color='k', linestyle = linestyle[2], label = label_texts[2]),
                Line2D([0], [0], color='k', linestyle = linestyle[3], label = label_texts[3])]

ax3.legend(handles = legend_elements, loc=('center right'), framealpha=1, labelspacing=0.2)  

plt.tight_layout()    

# Save figure
fig.savefig('MINC/compare_matrix_permeability.jpg', dpi = 150) 
plt.show()

 
