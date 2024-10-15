# type: ignore
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import h5py
import os

colors = [
  '#4c72b0', '#55a868', '#c44e52', '#8172b3', '#937860', '#da8bc3', '#8c8c8c',
  '#ccb974', '#64b5cd'
]
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

L_values = [0]
P_values = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12]

def get_data(dir):
  print(dir)

  N = []
  adm_masses = []
  adm_momenta = []
  centers_of_mass = []

  for L in L_values:
    for P in P_values:
      subdir = f'L{L}P{P:02}'

      try:
        adm_integrals = np.loadtxt(f'{dir}/{subdir}/BbhReductions/AdmIntegrals.dat', comments='#')

        number_of_points = adm_integrals[0]

        adm_mass = adm_integrals[1]

        adm_linear_momentum_x = adm_integrals[2]
        adm_linear_momentum_y = adm_integrals[3]
        adm_linear_momentum_z = adm_integrals[4]
        adm_linear_momentum = np.sqrt(adm_linear_momentum_x**2 + adm_linear_momentum_y**2 + adm_linear_momentum_z**2)

        center_of_mass_x = adm_integrals[5]
        center_of_mass_y = adm_integrals[6]
        center_of_mass_z = adm_integrals[7]
        center_of_mass = np.sqrt(center_of_mass_x**2 + center_of_mass_y**2 + center_of_mass_z**2)

        N.append(number_of_points)
        adm_masses.append(adm_mass)
        adm_momenta.append(adm_linear_momentum)
        centers_of_mass.append(center_of_mass)

      except Exception as e:
        print(f'\tSkipped {subdir} due to error: {e}')

  N = np.array(N)
  adm_masses = np.array(adm_masses)
  adm_momenta = np.array(adm_momenta)
  centers_of_mass = np.array(centers_of_mass)

  sort_indices = np.argsort(N)
  N = N[sort_indices]
  adm_masses = adm_masses[sort_indices]
  adm_momenta = adm_momenta[sort_indices]
  centers_of_mass = centers_of_mass[sort_indices]

  N3rt = N**(1./3.)

  N3rt = N3rt[:-1]
  D_adm_masses = np.abs(adm_masses - adm_masses[-1])[:-1]
  D_adm_momenta = np.abs(adm_momenta - adm_momenta[-1])[:-1]
  D_centers_of_mass = np.abs(centers_of_mass - centers_of_mass[-1])[:-1]

  return {
    'N3rt': N3rt,
    'D_adm_masses': D_adm_masses,
    'D_adm_momenta': D_adm_momenta,
    'D_centers_of_mass': D_centers_of_mass,
  }

fig, axes = plt.subplots(1, 2, squeeze=True)
col = 0

data_min = np.Infinity
data_max = -np.Infinity

def plot(title, data):
  global col, data_min, data_max
  ax = axes[col]
  col += 1

  ax.set_title(title)

  ax.plot(data['N3rt'], data['D_adm_masses'], label=r'$\Delta M_{ADM}$')
  ax.plot(data['N3rt'], data['D_adm_momenta'], label=r'$\Delta P_{ADM}$')
  ax.plot(data['N3rt'], data['D_centers_of_mass'], label=r'$\Delta C_{CoM}$')

  for values in [data['D_adm_masses'], data['D_adm_momenta'], data['D_centers_of_mass']]:
    data_min = min(data_min, np.min(values))
    data_max = max(data_max, np.max(values))

plot(
  'surface',
  get_data('./data/resolution_convergence/q3-surface')
)
plot(
  'volume',
  get_data('./data/resolution_convergence/q3-volume')
)

for ax in axes:
  ax.set_xlabel(r'$N^{1/3}$')

  ax.set_xscale('linear')
  ax.set_yscale('log')

  ax.set_ylim(data_min/2., data_max*2.)

  ax.grid('on', linestyle='--', alpha=0.5)

# Remove y tick labels
for ax in axes[1:]:
  ax.set_yticklabels([])

# Add legend to one axis
axes[0].legend()

fig.set_size_inches(10,6)
plt.tight_layout()
plt.subplots_adjust(wspace=0.02)

fig.savefig(f'resolution_convergence.png', format='png', bbox_inches='tight')

plt.show()
