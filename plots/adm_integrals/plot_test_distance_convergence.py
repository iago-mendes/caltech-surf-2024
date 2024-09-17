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

# If True, takes the difference between the current distance and the highest
# distance values. If False, takes the difference with the reference value.
# maxR_diff = True
maxR_diff = False

fig, axes = plt.subplots(3, 2, squeeze=True)
axes = np.transpose(axes)

data_min = [np.Infinity, np.Infinity, np.Infinity]
data_max = [-np.Infinity, -np.Infinity, -np.Infinity]

for i in range(2):
  dir = dirs[i]
  ax1, ax2, ax3 = axes[i]

  adm_mass_entries = np.genfromtxt(f'{dir}/Test_AdmMass.output', delimiter=',')
  adm_linear_momentum_entries = np.genfromtxt(f'{dir}/Test_AdmLinearMomentum.output', delimiter=',')
  center_of_mass_entries = np.genfromtxt(f'{dir}/Test_CenterOfMass.output', delimiter=',')

  RL_dict = {}
  for row in range(len(adm_mass_entries)):
      R, L, P, adm_mass, expected_adm_mass = adm_mass_entries[row]
      _, _, _, adm_linear_momentum, expected_adm_linear_momentum = adm_linear_momentum_entries[row]
      _, _, _, center_of_mass, expected_center_of_mass = center_of_mass_entries[row]

      if R not in RL_dict:
        RL_dict[R] = {}
      if L not in RL_dict[R]:
        RL_dict[R][L] = {
          "P_values": [],
          "adm_masses": [],
          "expected_adm_masses": [],
          "adm_linear_momenta": [],
          "expected_adm_linear_momenta": [],
          "centers_of_mass": [],
          "expected_centers_of_mass": [],
        }
      
      RL_dict[R][L]["P_values"].append(P)
      RL_dict[R][L]["adm_masses"].append(adm_mass)
      RL_dict[R][L]["expected_adm_masses"].append(expected_adm_mass)
      RL_dict[R][L]["adm_linear_momenta"].append(adm_linear_momentum)
      RL_dict[R][L]["expected_adm_linear_momenta"].append(expected_adm_linear_momentum)
      RL_dict[R][L]["centers_of_mass"].append(center_of_mass)
      RL_dict[R][L]["expected_centers_of_mass"].append(expected_center_of_mass)

  LP_counter = 0
  for L, P in [(1,6), (2,6), (2,15)]:
    adm_mass_residuals = []
    adm_linear_momentum_residuals = []
    center_of_mass_residuals = []

    R_values = sorted(RL_dict.keys())
    maxR = R_values[-1]
    for R in RL_dict:
      index = RL_dict[R][L]["P_values"].index(P)
      if maxR_diff:
        adm_mass_residuals.append(abs(RL_dict[R][L]["adm_masses"][index] - RL_dict[maxR][L]["adm_masses"][index]))
        adm_linear_momentum_residuals.append(abs(RL_dict[R][L]["adm_linear_momenta"][index] - RL_dict[maxR][L]["adm_linear_momenta"][index]))
        center_of_mass_residuals.append(abs(RL_dict[R][L]["centers_of_mass"][index] - RL_dict[maxR][L]["centers_of_mass"][index]))
      else:
        adm_mass_residuals.append(abs(RL_dict[R][L]["adm_masses"][index] - RL_dict[R][L]["expected_adm_masses"][index]))
        adm_linear_momentum_residuals.append(abs(RL_dict[R][L]["adm_linear_momenta"][index] - RL_dict[R][L]["expected_adm_linear_momenta"][index]))
        center_of_mass_residuals.append(abs(RL_dict[R][L]["centers_of_mass"][index] - RL_dict[R][L]["expected_centers_of_mass"][index]))
    
    if maxR_diff:
      R_values = R_values[:-1]
      adm_mass_residuals = adm_mass_residuals[:-1]
      adm_linear_momentum_residuals = adm_linear_momentum_residuals[:-1]
      center_of_mass_residuals = center_of_mass_residuals[:-1]
    
    adm_mass_residuals = np.array(adm_mass_residuals)
    adm_linear_momentum_residuals = np.array(adm_linear_momentum_residuals)
    center_of_mass_residuals = np.array(center_of_mass_residuals)

    # tmp = np.array(adm_mass_residuals, adm_linear_momentum_residuals, center_of_mass_residuals).flatten()
    # print(np.min(tmp[np.nonzero(tmp)]))

    ax1.plot(R_values, adm_mass_residuals, label=f'$L = {L}, P = {P}$', color=colors[LP_counter], marker=markers[LP_counter])
    ax2.plot(R_values, adm_linear_momentum_residuals, label=f'$L = {L}, P = {P}$', color=colors[LP_counter], marker=markers[LP_counter])
    ax3.plot(R_values, center_of_mass_residuals, label=f'$L = {L}, P = {P}$', color=colors[LP_counter], marker=markers[LP_counter])
    LP_counter += 1

    
    for i, row in enumerate([adm_mass_residuals, adm_linear_momentum_residuals, center_of_mass_residuals]):
      if not len(np.nonzero(row)[0]) == 0:
        data_min[i] = min(data_min[i], np.min(row[np.nonzero(row)]))
      data_max[i] = max(data_max[i], np.max(row))

for ax in axes.flatten():
  ax.set_xscale('log')
  ax.grid('on', linestyle='--', alpha=0.5)

for i, row in enumerate(np.transpose(axes)):
  y_min = data_min[i]/2.
  y_max = data_max[i]*2.
  for ax in row:
    if y_min == np.Infinity:
      ax.set_yscale('linear')
    else:
      ax.set_yscale('log')
      ax.set_ylim(y_min, y_max)

# Remove y tick labels
for ax in axes[1:,:].flatten():
  ax.set_yticklabels([])

# Remove x tick labels
for ax in axes[:,:-1].flatten():
  ax.set_xticklabels([])

# Add x axis labels
for bottom_ax in np.transpose(axes)[-1]:
  bottom_ax.set_xlabel(f'Outer radii $R$')

# Add y axis labels
if maxR_diff:
  axes[0,0].set_ylabel(r'$\Bigg|M_{ADM} - M_{ADM}|_{\max R}\Bigg|$')
  axes[0,1].set_ylabel(r'$\Bigg|P_{ADM} - P_{ADM}|_{\max R}\Bigg|$')
  axes[0,2].set_ylabel(r'$\Bigg|C_{CoM} - C_{CoM}|_{\max R}\Bigg|$')
else:
  axes[0,0].set_ylabel(r'$\Bigg|M_{ADM} - M_{ADM}|_{ref}\Bigg|$')
  axes[0,1].set_ylabel(r'$\Bigg|P_{ADM} - P_{ADM}|_{ref}\Bigg|$')
  axes[0,2].set_ylabel(r'$\Bigg|C_{CoM} - C_{CoM}|_{ref}\Bigg|$')

for i in range(2):
  dir = dirs[i]
  top_ax = np.transpose(axes)[0, i]
  with open(f'{dir}/title.txt', 'r') as title:
    top_ax.set_title(f'{title.read()}', fontsize=15)

top_right_ax = axes[-1, 0]
top_right_ax.legend()

# fig.set_size_inches(8, 8)
fig.set_size_inches(10, 10)
plt.tight_layout()
plt.subplots_adjust(wspace=0.03, hspace=0.03)
# fig.savefig(f'test_distance_convergence-{name}-{"maxR_diff" if maxR_diff else "ref_diff"}.pdf', format='pdf', bbox_inches='tight')
fig.savefig(f'test_distance_convergence-{name}-{"maxR_diff" if maxR_diff else "ref_diff"}.png', format='png', bbox_inches='tight')

plt.show()
