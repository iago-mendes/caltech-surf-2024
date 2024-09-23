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
# ref_diff = True
# ref_diff = False

resolutions = [
  [1, 6],
  [1, 12],
  [2, 6],
]
distance_exps = [2, 3, 4, 5, 6, 7, 8]

fig, axes = plt.subplots(3, 2, squeeze=True)
axes = np.transpose(axes)

data_min = [np.Infinity, np.Infinity, np.Infinity]
data_max = [-np.Infinity, -np.Infinity, -np.Infinity]

for i in range(2):
  dir = dirs[i]
  ax1, ax2, ax3 = axes[i]

  print(dir)

  counter = 0
  for L, P in resolutions:
    distances = []
    adm_masses = []
    adm_linear_momenta = []
    centers_of_mass = []

    for R in distance_exps:
      subdir = f'R{R}L{L}P{P:02}'

      try:
        reductions = h5py.File(f'{dir}/{subdir}/BbhReductions.h5')

        adm_mass = reductions['AdmIntegrals.dat'][0,0]

        adm_linear_momentum_x = reductions['AdmIntegrals.dat'][0,1]
        adm_linear_momentum_y = reductions['AdmIntegrals.dat'][0,2]
        adm_linear_momentum_z = reductions['AdmIntegrals.dat'][0,3]
        adm_linear_momentum = np.sqrt(adm_linear_momentum_x**2 + adm_linear_momentum_y**2 + adm_linear_momentum_z**2)

        center_of_mass_x = reductions['AdmIntegrals.dat'][0,4]
        center_of_mass_y = reductions['AdmIntegrals.dat'][0,5]
        center_of_mass_z = reductions['AdmIntegrals.dat'][0,6]
        # center_of_mass = np.sqrt(center_of_mass_x**2 + center_of_mass_y**2 + center_of_mass_z**2)
        center_of_mass = np.abs(center_of_mass_z)

        distances.append(float(f'1.e{R}'))
        adm_masses.append(adm_mass)
        adm_linear_momenta.append(adm_linear_momentum)
        centers_of_mass.append(center_of_mass)

      except Exception as e:
        print(f'\tSkipped {subdir}')
    
    if len(distances) == 0:
      print(f'\tNo data for L={L}, P={P}')
      continue

    distances = np.array(distances)
    adm_masses = np.array(adm_masses)
    adm_linear_momenta = np.array(adm_linear_momenta)
    centers_of_mass = np.array(centers_of_mass)

    # adm_masses -= adm_masses[-1]
    # adm_linear_momenta -= adm_linear_momenta[-1]
    # with open(f'{dir}/ref_CoM.txt', 'r') as ref_CoM:
    #   centers_of_mass -= float(ref_CoM.read())

    # distances = distances[:-1]
    # adm_masses = np.abs(adm_masses[:-1])
    # adm_linear_momenta = np.abs(adm_linear_momenta[:-1])
    # centers_of_mass = np.abs(centers_of_mass[:-1])

    ax1.plot(distances, adm_masses, label=f'$L = {L}, P = {P}$', color=colors[counter], marker=markers[counter])
    ax2.plot(distances, adm_linear_momenta, label=f'$L = {L}, P = {P}$', color=colors[counter], marker=markers[counter])
    ax3.plot(distances, centers_of_mass, label=f'$L = {L}, P = {P}$', color=colors[counter], marker=markers[counter])

    for j, row in enumerate([adm_masses, adm_linear_momenta, centers_of_mass]):
      data_min[j] = min(data_min[j], np.min(row))
      data_max[j] = max(data_max[j], np.max(row))

    counter += 1


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
axes[0,0].set_ylabel(r'$\Bigg|M_{ADM} - M_{ADM}|_{\max R}\Bigg|$')
axes[0,1].set_ylabel(r'$\Bigg|P_{ADM} - P_{ADM}|_{\max R}\Bigg|$')
# axes[0,2].set_ylabel(r'$\Bigg|C_{CoM} - C_{CoM}|_{\max R}\Bigg|$')
axes[0,2].set_ylabel(r'$\Bigg|C_{CoM} - C_{CoM}|_{ref}\Bigg|$')

for i in range(2):
  dir = dirs[i]
  top_ax = np.transpose(axes)[0, i]
  with open(f'{dir}/title.txt', 'r') as title:
    top_ax.set_title(f'{title.read()}', fontsize=15)

top_right_ax = axes[-1, 0]
top_right_ax.legend()

fig.set_size_inches(10, 10)
plt.tight_layout()
plt.subplots_adjust(wspace=0.03, hspace=0.03)

plt.show()


# ###################################################################


# for ax in axes.flatten():
#   ax.set_yscale('log')
#   ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#   ax.grid('on', linestyle='--', alpha=0.5)

# for i, row in enumerate(np.transpose(axes)):
#   y_min = data_min[i]/2.
#   # y_max = data_max[i]*2.
#   # y_max = min(data_max[i]*2., 1.e1)
#   y_max = min(data_max[i]*2., 1.e2)
#   for ax in row:
#     ax.set_ylim(y_min, y_max)

# # Remove y tick labels
# for ax in axes[1:,:].flatten():
#   ax.set_yticklabels([])

# # Remove x tick labels
# for ax in axes[:,:-1].flatten():
#   ax.set_xticklabels([])

# # Add x axis labels
# for bottom_ax in np.transpose(axes)[-1]:
#   bottom_ax.set_xlabel(f'Polynomial order $P$')

# # Add y axis labels
# if ref_diff:
#   axes[0,0].set_ylabel(r'$M_{ADM} - M_{ADM}|_{\max P}$')
#   axes[0,1].set_ylabel(r'$P_{ADM} - P_{ADM}|_{\max P}$')
#   axes[0,2].set_ylabel(r'$C_{CoM} - C_{CoM}|_{\max P}$')
# else:
#   axes[0,0].set_ylabel(r'$M_{ADM} - M_{ADM}|_{P-1}$')
#   axes[0,1].set_ylabel(r'$P_{ADM} - P_{ADM}|_{P-1}$')
#   axes[0,2].set_ylabel(r'$C_{CoM} - C_{CoM}|_{P-1}$')

# for i in range(2):
#   dir = dirs[i]
#   top_ax = np.transpose(axes)[0, i]
#   with open(f'{dir}/title.txt', 'r') as title:
#     top_ax.set_title(f'{title.read()}', fontsize=15)

# top_right_ax = axes[-1, 0]
# top_right_ax.legend()

# fig.set_size_inches(8, 8)
# plt.tight_layout()
# plt.subplots_adjust(wspace=0.02, hspace=0.03)
# fig.savefig(f'bbh_resolution_convergence-{name}-{"ref_diff" if ref_diff else "next_diff"}.pdf', format='pdf', bbox_inches='tight')

# plt.show()
