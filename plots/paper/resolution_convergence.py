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

L_values = [0]
P_values = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

def get_data(dir):
  print(dir)

  N = []
  adm_masses = []
  masses_irr_a = []
  masses_irr_b = []
  masses_a = []
  masses_b = []
  adm_momenta_x = []
  adm_momenta_y = []
  adm_momenta_z = []
  centers_of_mass_x = []
  centers_of_mass_y = []
  centers_of_mass_z = []
  spins_a_x = []
  spins_a_y = []
  spins_a_z = []
  spins_b_x = []
  spins_b_y = []
  spins_b_z = []
  excision_areas = []

  for L in L_values:
    for P in P_values:
      subdir = f'L{L}P{P:02}'

      try:
        adm_integrals = np.loadtxt(f'{dir}/{subdir}/BbhReductions/AdmIntegrals.dat', comments='#')
        AhA = np.loadtxt(f'{dir}/{subdir}/Horizons/AhA.dat', comments='#')
        AhB = np.loadtxt(f'{dir}/{subdir}/Horizons/AhB.dat', comments='#')

        number_of_points = adm_integrals[0]

        adm_mass = adm_integrals[1]

        mass_irr_a = AhA[1]
        mass_irr_b = AhB[1]

        mass_a = AhA[4]
        mass_b = AhB[4]

        adm_linear_momentum_x = adm_integrals[2]
        adm_linear_momentum_y = adm_integrals[3]
        adm_linear_momentum_z = adm_integrals[4]

        center_of_mass_x = adm_integrals[5]
        center_of_mass_y = adm_integrals[6]
        center_of_mass_z = adm_integrals[7]

        spin_a_x = AhA[6]
        spin_a_y = AhA[7]
        spin_a_z = AhA[8]

        spin_b_x = AhB[6]
        spin_b_y = AhB[7]
        spin_b_z = AhB[8]

        excision_area = adm_integrals[8]

        N.append(number_of_points)
        adm_masses.append(adm_mass)
        masses_irr_a.append(mass_irr_a)
        masses_irr_b.append(mass_irr_b)
        masses_a.append(mass_a)
        masses_b.append(mass_b)
        adm_momenta_x.append(adm_linear_momentum_x)
        adm_momenta_y.append(adm_linear_momentum_y)
        adm_momenta_z.append(adm_linear_momentum_z)
        centers_of_mass_x.append(center_of_mass_x)
        centers_of_mass_y.append(center_of_mass_y)
        centers_of_mass_z.append(center_of_mass_z)
        spins_a_x.append(spin_a_x)
        spins_a_y.append(spin_a_y)
        spins_a_z.append(spin_a_z)
        spins_b_x.append(spin_b_x)
        spins_b_y.append(spin_b_y)
        spins_b_z.append(spin_b_z)
        excision_areas.append(excision_area)

      except Exception as e:
        print(f'\tSkipped {subdir} due to error: {e}')

  N = np.array(N)
  adm_masses = np.array(adm_masses)
  masses_irr_a = np.array(masses_irr_a)
  masses_irr_b = np.array(masses_irr_b)
  masses_a = np.array(masses_a)
  masses_b = np.array(masses_b)
  adm_momenta_x = np.array(adm_momenta_x)
  adm_momenta_y = np.array(adm_momenta_y)
  adm_momenta_z = np.array(adm_momenta_z)
  centers_of_mass_x = np.array(centers_of_mass_x)
  centers_of_mass_y = np.array(centers_of_mass_y)
  centers_of_mass_z = np.array(centers_of_mass_z)
  spins_a_x = np.array(spins_a_x)
  spins_a_y = np.array(spins_a_y)
  spins_a_z = np.array(spins_a_z)
  spins_b_x = np.array(spins_b_x)
  spins_b_y = np.array(spins_b_y)
  spins_b_z = np.array(spins_b_z)
  excision_areas = np.array(excision_areas)

  sort_indices = np.argsort(N)
  N = N[sort_indices]
  adm_masses = adm_masses[sort_indices]
  masses_irr_a = masses_irr_a[sort_indices]
  masses_irr_b = masses_irr_b[sort_indices]
  masses_a = masses_a[sort_indices]
  masses_b = masses_b[sort_indices]
  adm_momenta_x = adm_momenta_x[sort_indices]
  adm_momenta_y = adm_momenta_y[sort_indices]
  adm_momenta_z = adm_momenta_z[sort_indices]
  centers_of_mass_x = centers_of_mass_x[sort_indices]
  centers_of_mass_y = centers_of_mass_y[sort_indices]
  centers_of_mass_z = centers_of_mass_z[sort_indices]
  spins_a_x = spins_a_x[sort_indices]
  spins_a_y = spins_a_y[sort_indices]
  spins_a_z = spins_a_z[sort_indices]
  spins_b_x = spins_b_x[sort_indices]
  spins_b_y = spins_b_y[sort_indices]
  spins_b_z = spins_b_z[sort_indices]
  excision_areas = excision_areas[sort_indices]

  excision_masses_irr = np.sqrt(excision_areas / (16. * np.pi))

  N3rt = N**(1./3.)

  N3rt = N3rt[:-1]
  D_adm_masses = np.abs(adm_masses - adm_masses[-1])[:-1]
  D_masses_irr_a = np.abs(masses_irr_a - masses_irr_a[-1])[:-1]
  D_masses_irr_b = np.abs(masses_irr_b - masses_irr_b[-1])[:-1]
  D_masses_a = np.abs(masses_a - masses_a[-1])[:-1]
  D_masses_b = np.abs(masses_b - masses_b[-1])[:-1]
  D_adm_momenta_x = np.abs(adm_momenta_x - adm_momenta_x[-1])[:-1]
  D_adm_momenta_y = np.abs(adm_momenta_y - adm_momenta_y[-1])[:-1]
  D_adm_momenta_z = np.abs(adm_momenta_z - adm_momenta_z[-1])[:-1]
  D_centers_of_mass_x = np.abs(centers_of_mass_x - centers_of_mass_x[-1])[:-1]
  D_centers_of_mass_y = np.abs(centers_of_mass_y - centers_of_mass_y[-1])[:-1]
  D_centers_of_mass_z = np.abs(centers_of_mass_z - centers_of_mass_z[-1])[:-1]
  D_spins_a_x = np.abs(spins_a_x - spins_a_x[-1])[:-1]
  D_spins_a_y = np.abs(spins_a_y - spins_a_y[-1])[:-1]
  D_spins_a_z = np.abs(spins_a_z - spins_a_z[-1])[:-1]
  D_spins_b_x = np.abs(spins_b_x - spins_b_x[-1])[:-1]
  D_spins_b_y = np.abs(spins_b_y - spins_b_y[-1])[:-1]
  D_spins_b_z = np.abs(spins_b_z - spins_b_z[-1])[:-1]
  D_excision_masses_irr = np.abs(excision_masses_irr - excision_masses_irr[-1])[:-1]

  D_adm_momenta = np.sqrt(D_adm_momenta_x**2 + D_adm_momenta_y**2 + D_adm_momenta_z**2)
  D_centers_of_mass = np.sqrt(D_centers_of_mass_x**2 + D_centers_of_mass_y**2 + D_centers_of_mass_z**2)
  D_spins_a = np.sqrt(D_spins_a_x**2 + D_spins_a_y**2 + D_spins_a_z**2)
  D_spins_b = np.sqrt(D_spins_b_x**2 + D_spins_b_y**2 + D_spins_b_z**2)

  return {
    'N3rt': N3rt,
    'D_adm_masses': D_adm_masses,
    'D_masses_irr_a': D_masses_irr_a,
    'D_masses_irr_b': D_masses_irr_b,
    'D_masses_a': D_masses_a,
    'D_masses_b': D_masses_b,
    'D_adm_momenta': D_adm_momenta,
    'D_centers_of_mass': D_centers_of_mass,
    'D_spins_a': D_spins_a,
    'D_spins_b': D_spins_b,
    'D_excision_masses_irr': D_excision_masses_irr
  }

fig, axes = plt.subplots(1, 2, squeeze=True)
col = 0

data_min = np.Infinity
data_max = -np.Infinity

def plot(title, data, plot_horizon_quantities=False):
  global col, data_min, data_max
  ax = axes[col]
  col += 1

  ax.set_title(title)

  if plot_horizon_quantities:
    ax.plot(data['N3rt'], data['D_excision_masses_irr'], label=r'$\Delta \sqrt{A_{S_0}/16\pi}$', color=colors[0], marker=markers[0])
    ax.plot(data['N3rt'], data['D_masses_irr_a'], label=r'$\Delta M^{irr}_A$', color=colors[1], marker=markers[1])
    ax.plot(data['N3rt'], data['D_masses_irr_b'], label=r'$\Delta M^{irr}_B$', color=colors[1], marker=markers[1], linestyle='dashed')
    ax.plot(data['N3rt'], data['D_masses_a'], label=r'$\Delta M_A$', color=colors[2], marker=markers[2])
    ax.plot(data['N3rt'], data['D_masses_b'], label=r'$\Delta M_B$', color=colors[2], marker=markers[2], linestyle='dashed')
    ax.plot(data['N3rt'], data['D_spins_a'], label=r'$\Delta \chi_A$', color=colors[3], marker=markers[3])
    ax.plot(data['N3rt'], data['D_spins_b'], label=r'$\Delta \chi_B$', color=colors[3], marker=markers[3], linestyle='dashed')

    for values in [
      data['D_excision_masses_irr'],
      data['D_masses_irr_a'],
      data['D_masses_irr_b'],
      data['D_masses_a'],
      data['D_masses_b'],
      data['D_spins_a'],
      data['D_spins_b']
    ]:
      data_min = min(data_min, np.min(values))
      data_max = max(data_max, np.max(values))
  else:
    ax.plot(data['N3rt'], data['D_adm_masses'], label=r'$\Delta M_{ADM}$')
    ax.plot(data['N3rt'], data['D_adm_momenta'], label=r'$\Delta P_{ADM}$')
    ax.plot(data['N3rt'], data['D_centers_of_mass'], label=r'$\Delta C_{CoM}$')

    for values in [data['D_adm_masses'], data['D_adm_momenta'], data['D_centers_of_mass']]:
      data_min = min(data_min, np.min(values))
      data_max = max(data_max, np.max(values))

plot(
  'Horizon quantities',
  get_data('./data/resolution_convergence/q3-surface'),
  plot_horizon_quantities=True
)
plot(
  'Asymptotic surface integrals',
  get_data('./data/resolution_convergence/q3-surface')
)
# plot(
#   'Asymptotic volume integrals',
#   get_data('./data/resolution_convergence/q3-volume')
# )

for ax in axes:
  ax.set_xlabel(r'$N^{1/3}$')
  # ax.legend(loc='upper right')
  ax.legend()

  ax.set_xscale('linear')
  ax.set_yscale('log')

  ax.set_ylim(data_min/2., data_max*2.)

  ax.grid('on', linestyle='--', alpha=0.5)

# Remove y tick labels
for ax in axes[1:]:
  ax.set_yticklabels([])

fig.set_size_inches(10,7)
plt.tight_layout()
plt.subplots_adjust(wspace=0.02)

fig.savefig(f'resolution_convergence.png', format='png', bbox_inches='tight')

plt.show()
