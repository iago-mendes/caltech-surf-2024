# type: ignore
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import sys
from cycler import cycler
import h5py
import re
import os

colors = ['#4c72b0', '#55a868', '#c44e52', '#8172b3', '#937860', '#da8bc3', '#8c8c8c', '#ccb974', '#64b5cd']
markers = ['o', '^', 's', 'v', 'D', '<', '>', '*', '+']
plt.rcParams.update({
	'text.usetex' : True,
	'font.family' : 'serif',
	'font.serif' : ['Computer Modern Serif'],
	'font.size': 15,
	'axes.prop_cycle': cycler('color', colors) +
                     cycler('marker', markers) +
										 cycler('markerfacecolor', ['none'] * len(colors))
})

name = sys.argv[1]
dirs = [sys.argv[2], sys.argv[3]]

# If True, takes the difference between the current resolution and the highest
# resolution values. If False, takes the difference between the current
# resolution and the next lower resolution values.
ref_diff = True
# ref_diff = False

fig, axes = plt.subplots(3, 2, squeeze=True)
axes = np.transpose(axes)

data_min = [np.Infinity, np.Infinity, np.Infinity]
data_max = [-np.Infinity, -np.Infinity, -np.Infinity]

for i in range(2):
  dir = dirs[i]
  ax1, ax2, ax3 = axes[i]

  # List of directories
  directories = ["L0P02", "L0P03", "L1P09", "L1P10", "L2P05"]

  # Initialize a dictionary to store L and respective P values
  L_P_dict = {}

  # Regular expression to extract L and P values
  pattern = re.compile(r'L(\d+)P(\d+)')

  # Loop through the sub-directories
  for sub_dir in os.listdir(dir):
    match = pattern.match(sub_dir)
    if match:
      L_value = int(match.group(1))
      P_value = int(match.group(2))
        
      if L_value in L_P_dict:
        L_P_dict[L_value].append(P_value)
      else:
        L_P_dict[L_value] = [P_value]
  
  # Take L=1 data from "mixed convergece" (old data)
  old_data = False
  if len(L_P_dict) == 0 and os.path.exists(f'{dir}/R5P06'):
    old_data = True
    L_value = 1
    pattern = re.compile(r'R5P(\d+)')
    for sub_dir in os.listdir(dir):
      match = pattern.match(sub_dir)
      if match:
        P_value = int(match.group(1))
          
        if L_value in L_P_dict:
          L_P_dict[L_value].append(P_value)
        else:
          L_P_dict[L_value] = [P_value]
  
  print(L_P_dict)
  h_refinements = sorted(L_P_dict.keys())
  for L in h_refinements:
    adm_masses = []
    adm_linear_momenta = []
    centers_of_mass = []

    p_refinements = sorted(L_P_dict[L])
    used_P_values = []
    for P in p_refinements:
      try:
        reductions = h5py.File(f'{dir}/L{L}P{P:02}/BbhReductions.h5') if not old_data else h5py.File(f'{dir}/R5P{P:02}/BbhReductions.h5')

        adm_mass = reductions['AdmIntegrals.dat'][0,0]

        adm_linear_momentum_x = reductions['AdmIntegrals.dat'][0,1]
        adm_linear_momentum_y = reductions['AdmIntegrals.dat'][0,2]
        adm_linear_momentum_z = reductions['AdmIntegrals.dat'][0,3]
        adm_linear_momentum = np.sqrt(adm_linear_momentum_x**2 + adm_linear_momentum_y**2 + adm_linear_momentum_z**2)

        center_of_mass_x = reductions['AdmIntegrals.dat'][0,4]
        center_of_mass_y = reductions['AdmIntegrals.dat'][0,5]
        center_of_mass_z = reductions['AdmIntegrals.dat'][0,6]
        center_of_mass = np.sqrt(center_of_mass_x**2 + center_of_mass_y**2 + center_of_mass_z**2)

        adm_masses.append(adm_mass)
        adm_linear_momenta.append(adm_linear_momentum)
        centers_of_mass.append(center_of_mass)

        print(f'Using L = {L}, P = {P}')
        used_P_values.append(P)
      except Exception as e:
        print(f"Skipping L = {L}, P = {P}: {str(e)}")

    
    adm_masses = np.array(adm_masses)
    adm_linear_momenta = np.array(adm_linear_momenta)
    centers_of_mass = np.array(centers_of_mass)

    if ref_diff:
      # Take difference with the highest resolution value
      adm_masses = np.array(adm_masses)
      adm_linear_momenta = np.array(adm_linear_momenta)
      centers_of_mass = np.array(centers_of_mass)
      adm_masses -= adm_masses[-1]
      adm_linear_momenta -= adm_linear_momenta[-1]
      centers_of_mass -= centers_of_mass[-1]
      adm_masses = np.abs(adm_masses[:-1])
      adm_linear_momenta = np.abs(adm_linear_momenta[:-1])
      centers_of_mass = np.abs(centers_of_mass[:-1])
      P_values = used_P_values[:-1]
    else:
      # Take difference with the next lower resolution value
      adm_masses = np.abs(np.diff(adm_masses))
      adm_linear_momenta = np.abs(np.diff(adm_linear_momenta))
      centers_of_mass = np.abs(np.diff(centers_of_mass))
      P_values = used_P_values[1:]

    ax1.plot(P_values, adm_masses, label=f'$L = {L}$', color=colors[L], marker=markers[L])
    ax2.plot(P_values, adm_linear_momenta, label=f'$L = {L}$', color=colors[L], marker=markers[L])
    ax3.plot(P_values, centers_of_mass, label=f'$L = {L}$', color=colors[L], marker=markers[L])

    for i, row in enumerate([adm_masses, adm_linear_momenta, centers_of_mass]):
      data_min[i] = min(data_min[i], np.min(row))
      data_max[i] = max(data_max[i], np.max(row))

for ax in axes.flatten():
  ax.set_yscale('log')
  ax.xaxis.set_major_locator(MaxNLocator(integer=True))
  ax.grid('on', linestyle='--', alpha=0.5)

for i, row in enumerate(np.transpose(axes)):
  y_min = data_min[i]/2.
  # y_max = data_max[i]*2.
  # y_max = min(data_max[i]*2., 1.e1)
  y_max = min(data_max[i]*2., 1.e2)
  for ax in row:
    ax.set_ylim(y_min, y_max)

# Remove y tick labels
for ax in axes[1:,:].flatten():
  ax.set_yticklabels([])

# Remove x tick labels
for ax in axes[:,:-1].flatten():
  ax.set_xticklabels([])

# Add x axis labels
for bottom_ax in np.transpose(axes)[-1]:
  bottom_ax.set_xlabel(f'Polynomial order $P$')

# Add y axis labels
if ref_diff:
  axes[0,0].set_ylabel(r'$M_{ADM} - M_{ADM}|_{\max P}$')
  axes[0,1].set_ylabel(r'$P_{ADM} - P_{ADM}|_{\max P}$')
  axes[0,2].set_ylabel(r'$C_{CoM} - C_{CoM}|_{\max P}$')
else:
  axes[0,0].set_ylabel(r'$M_{ADM} - M_{ADM}|_{P-1}$')
  axes[0,1].set_ylabel(r'$P_{ADM} - P_{ADM}|_{P-1}$')
  axes[0,2].set_ylabel(r'$C_{CoM} - C_{CoM}|_{P-1}$')

for i in range(2):
  dir = dirs[i]
  top_ax = np.transpose(axes)[0, i]
  with open(f'{dir}/title.txt', 'r') as title:
    top_ax.set_title(f'{title.read()}', fontsize=15)

top_right_ax = axes[-1, 0]
top_right_ax.legend()

fig.set_size_inches(8, 8)
plt.tight_layout()
plt.subplots_adjust(wspace=0.02, hspace=0.03)
fig.savefig(f'bbh_resolution_convergence-{name}-{"ref_diff" if ref_diff else "next_diff"}.pdf', format='pdf', bbox_inches='tight')

plt.show()
