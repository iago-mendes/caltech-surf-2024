# type: ignore
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import h5py
import os

colors = [
  '#4c72b0', '#55a868', '#c44e52', '#ccb974', '#937860', '#8172b3', '#8c8c8c',
  '#da8bc3', '#64b5cd'
]
markers = ['o', 's', 'D', '^', 'v', '<', '>', '*', '+']
plt.rcParams.update({
	'text.usetex' : True,
	'font.family' : 'serif',
	'font.serif' : ['Computer Modern Serif'],
	'font.size': 15,
	'axes.prop_cycle': cycler('color', colors) +
                     cycler('marker', markers) +
										 cycler('markerfacecolor', ['none'] * len(colors))
})

basedir = './data/resolution_convergence/horizon_finder'
dirs = [f'{basedir}/P09', f'{basedir}/P15']
# dirs = [f'{basedir}/P15', f'{basedir}/P15']
dir_titles = [
  r'L=0, P=9 (60000 points, $N^{1/3} \approx 39$)',
  r'L=0, P=15 (236544 points, $N^{1/3} \approx 62$)'
]

fig, axes = plt.subplots(1, 2, squeeze=True)

data_min = np.Infinity
data_max = -np.Infinity

for i, dir in enumerate(dirs):
  ax = axes[i]
  ax.set_title(dir_titles[i])

  L_values = []
  masses_irr_a = []
  masses_irr_b = []
  masses_a = []
  masses_b = []
  spins_a_x = []
  spins_a_y = []
  spins_a_z = []
  spins_b_x = []
  spins_b_y = []
  spins_b_z = []

  for L in range(5, 31+1):
    subdir = f'L{L:02}'
    try:
      AhA = np.loadtxt(f'{dir}/{subdir}/Horizons/AhA.dat', comments='#')
      AhB = np.loadtxt(f'{dir}/{subdir}/Horizons/AhB.dat', comments='#')

      mass_irr_a = AhA[1]
      mass_irr_b = AhB[1]

      mass_a = AhA[4]
      mass_b = AhB[4]

      spin_a_x = AhA[6]
      spin_a_y = AhA[7]
      spin_a_z = AhA[8]

      spin_b_x = AhB[6]
      spin_b_y = AhB[7]
      spin_b_z = AhB[8]

      L_values.append(L)
      masses_irr_a.append(mass_irr_a)
      masses_irr_b.append(mass_irr_b)
      masses_a.append(mass_a)
      masses_b.append(mass_b)
      spins_a_x.append(spin_a_x)
      spins_a_y.append(spin_a_y)
      spins_a_z.append(spin_a_z)
      spins_b_x.append(spin_b_x)
      spins_b_y.append(spin_b_y)
      spins_b_z.append(spin_b_z)
    except Exception as e:
      print(f'\tSkipped {subdir} due to error: {e}')

  L_values = np.array(L_values)
  masses_irr_a = np.array(masses_irr_a)
  masses_irr_b = np.array(masses_irr_b)
  masses_a = np.array(masses_a)
  masses_b = np.array(masses_b)
  spins_a_x = np.array(spins_a_x)
  spins_a_y = np.array(spins_a_y)
  spins_a_z = np.array(spins_a_z)
  spins_b_x = np.array(spins_b_x)
  spins_b_y = np.array(spins_b_y)
  spins_b_z = np.array(spins_b_z)

  L_values = L_values[:-1]
  D_masses_irr_a = np.abs(masses_irr_a - masses_irr_a[-1])[:-1]
  D_masses_irr_b = np.abs(masses_irr_b - masses_irr_b[-1])[:-1]
  D_masses_a = np.abs(masses_a - masses_a[-1])[:-1]
  D_masses_b = np.abs(masses_b - masses_b[-1])[:-1]
  D_spins_a_x = np.abs(spins_a_x - spins_a_x[-1])[:-1]
  D_spins_a_y = np.abs(spins_a_y - spins_a_y[-1])[:-1]
  D_spins_a_z = np.abs(spins_a_z - spins_a_z[-1])[:-1]
  D_spins_b_x = np.abs(spins_b_x - spins_b_x[-1])[:-1]
  D_spins_b_y = np.abs(spins_b_y - spins_b_y[-1])[:-1]
  D_spins_b_z = np.abs(spins_b_z - spins_b_z[-1])[:-1]

  D_spins_a = np.sqrt(D_spins_a_x**2 + D_spins_a_y**2 + D_spins_a_z**2)
  D_spins_b = np.sqrt(D_spins_b_x**2 + D_spins_b_y**2 + D_spins_b_z**2)

  ax.plot(L_values, D_masses_irr_a, label=r'$\Delta M^{irr}_A$', color=colors[0], marker=markers[0])
  ax.plot(L_values, D_masses_irr_b, label=r'$\Delta M^{irr}_B$', color=colors[0], marker=markers[0], linestyle='dashed')
  ax.plot(L_values, D_masses_a, label=r'$\Delta M_A$', color=colors[1], marker=markers[1])
  ax.plot(L_values, D_masses_b, label=r'$\Delta M_B$', color=colors[1], marker=markers[1], linestyle='dashed')
  ax.plot(L_values, D_spins_a, label=r'$\Delta \chi_A$', color=colors[2], marker=markers[2])
  ax.plot(L_values, D_spins_b, label=r'$\Delta \chi_B$', color=colors[2], marker=markers[2], linestyle='dashed')

  for values in [
    D_masses_irr_a,
    D_masses_irr_b,
    D_masses_a,
    D_masses_b,
    D_spins_a,
    D_spins_b
  ]:
    data_min = min(data_min, np.min(values))
    data_max = max(data_max, np.max(values))

for ax in axes:
  ax.set_xlabel(r'Horizon $L_{max}$')

  ax.set_xscale('linear')
  ax.set_yscale('log')
  ax.set_ylim(data_min/2., data_max*2.)

  ax.grid('on', linestyle='--', alpha=0.5)

axes[-1].legend(loc='center left', bbox_to_anchor=(0.99, 0.5))

# Remove y tick labels
for ax in axes[1:]:
  ax.set_yticklabels([])

fig.set_size_inches(12,7)
plt.tight_layout()
plt.subplots_adjust(wspace=0.02)

fig.savefig(f'horizon_finder_resolution.png', format='png', bbox_inches='tight')

plt.show()
