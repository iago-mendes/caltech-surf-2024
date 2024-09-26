# type: ignore
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import h5py
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

def get_bbh_residuals(dir):
  print(dir)

  residuals = []

  for L in [0, 1, 2]:
    P_values = []
    adm_masses = []
    adm_momenta = []
    centers_of_mass = []

    for P in [2, 4, 6, 8, 10, 12, 14]:
      subdir = f'R5L{L}P{P:02}'

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
        center_of_mass = np.sqrt(center_of_mass_x**2 + center_of_mass_y**2 + center_of_mass_z**2)

        P_values.append(P)
        adm_masses.append(adm_mass)
        adm_momenta.append(adm_linear_momentum)
        centers_of_mass.append(center_of_mass)

      except Exception as e:
        print(f'\tSkipped {subdir} due to error: {e}')
    
    adm_masses = np.array(adm_masses)
    adm_momenta = np.array(adm_momenta)
    centers_of_mass = np.array(centers_of_mass)

    adm_masses -= adm_masses[-1]
    adm_momenta -= adm_momenta[-1]
    centers_of_mass -= centers_of_mass[-1]

    P_values = P_values[:-1]
    adm_masses = np.abs(adm_masses[:-1])
    adm_momenta = np.abs(adm_momenta[:-1])
    centers_of_mass = np.abs(centers_of_mass[:-1])
      
    residuals.append({
      'P_values': P_values,
      'adm_masses': adm_masses,
      'adm_momenta': adm_momenta,
      'centers_of_mass': centers_of_mass,
    })

  return residuals

def get_test_residuals(dir):
  print(dir)

  residuals = []

  for L in [0, 1, 2]:
    P_values = []
    adm_masses = []
    adm_momenta = []
    centers_of_mass = []

    adm_mass_entries = np.genfromtxt(f'{dir}/Test_AdmMass.output', delimiter=',')
    adm_linear_momentum_entries = np.genfromtxt(f'{dir}/Test_AdmLinearMomentum.output', delimiter=',')

    center_of_mass_fname = ''
    if os.path.exists(f'{dir}/Test_CenterOfMass-new.output'):
      center_of_mass_fname = f'{dir}/Test_CenterOfMass-new.output'
    else:
      center_of_mass_fname = f'{dir}/Test_CenterOfMass.output'
    center_of_mass_entries = np.genfromtxt(center_of_mass_fname, delimiter=',')

    for row in range(len(adm_mass_entries)):
      R, row_L, P, adm_mass, _ = adm_mass_entries[row]
      _, _, _, adm_linear_momentum, _ = adm_linear_momentum_entries[row]

      if R != 1.e5 or row_L != L:
        continue

      P_values.append(P)
      adm_masses.append(adm_mass)
      adm_momenta.append(adm_linear_momentum)
    
    for row in range(len(center_of_mass_entries)):
      R = 0.
      row_L = 0.
      center_of_mass = 0.
      if len(center_of_mass_entries[row]) == 7:
        R_, row_L_, _, center_of_mass_x, center_of_mass_y, center_of_mass_z, _ = center_of_mass_entries[row]
        R = R_
        row_L = row_L_
        center_of_mass = np.sqrt(center_of_mass_x**2 + center_of_mass_y**2 + center_of_mass_z**2)
      else:
        R_, row_L_, _, center_of_mass_, _ = center_of_mass_entries[row]
        R = R_
        row_L = row_L_
        center_of_mass = center_of_mass_

      if R != 1.e5 or row_L != L:
        continue
      
      centers_of_mass.append(center_of_mass)
    
    adm_masses = np.array(adm_masses)
    adm_momenta = np.array(adm_momenta)
    centers_of_mass = np.array(centers_of_mass)

    adm_masses -= adm_masses[-1]
    adm_momenta -= adm_momenta[-1]
    centers_of_mass -= centers_of_mass[-1]

    P_values = P_values[:-1]
    adm_masses = np.abs(adm_masses[:-1])
    adm_momenta = np.abs(adm_momenta[:-1])
    centers_of_mass = np.abs(centers_of_mass[:-1])
      
    residuals.append({
      'P_values': P_values,
      'adm_masses': adm_masses,
      'adm_momenta': adm_momenta,
      'centers_of_mass': centers_of_mass,
    })

  return residuals

fig, axes = plt.subplots(3, 5, squeeze=True)
axes = np.transpose(axes)
col = 0

def plot(title, residuals):
  global col
  ax1, ax2, ax3 = axes[col]
  col += 1

  ax1.set_title(title)

  for i, label in enumerate([r'$L = 0$', r'$L = 1$', r'$L = 2$']):
    P_values = residuals[i]['P_values']
    ax1.plot(P_values, residuals[i]['adm_masses'], label=label, color=colors[i], marker=markers[i])
    ax2.plot(P_values, residuals[i]['adm_momenta'], label=label, color=colors[i], marker=markers[i])
    ax3.plot(P_values, residuals[i]['centers_of_mass'], label=label, color=colors[i], marker=markers[i])

plot(
  'Isotropic (surface)',
  get_test_residuals('../adm_integrals/data/isotropic/surface')
)
plot(
  'Isotropic (volume)',
  get_test_residuals('../adm_integrals/data/isotropic/volume')
)
plot(
  'Boosted Isotro. (surf.)',
  get_test_residuals('../adm_integrals/data/isotropic_boosted/surface')
)
plot(
  'BBH (surface)',
  get_bbh_residuals('../adm_integrals/data/bbh/surface-new')
)
plot(
  'BBH (volume)',
  get_bbh_residuals('../adm_integrals/data/bbh/volume')
)

for ax in axes.flatten():
  ax.set_xscale('linear')
  ax.set_yscale('log')
  ax.grid('on', linestyle='--', alpha=0.5)

# Adjust specific axes
# Isotropic Padm surface
axes[0,1].plot([2, 14], [1.e-15, 1.e-15], alpha=0)
axes[0,1].set_yscale('symlog', linthresh=1e-17)
# Isotropic Padm volume
axes[1,1].plot([2, 14], [1.e-15, 1.e-15], alpha=0)
axes[1,1].set_yscale('symlog', linthresh=1e-17)

# Remove x tick labels
for ax in axes[:,:-1].flatten():
  ax.set_xticklabels([])

# Add x axis labels
for bottom_ax in np.transpose(axes)[-1]:
  bottom_ax.set_xlabel(r'$P$')

# Add y axis labels
axes[0,0].set_ylabel(r'$\Bigg| M_{ADM} - M_{ADM}|_{\max P} \Bigg|$')
axes[0,1].set_ylabel(r'$\Bigg| P_{ADM} - P_{ADM}|_{\max P} \Bigg|$')
axes[0,2].set_ylabel(r'$\Bigg| C_{CoM} - C_{CoM}|_{\max P} \Bigg|$')

# Add legend to one axis
axes[0, 0].legend()

fig.set_size_inches(14,10)
plt.tight_layout()
plt.subplots_adjust(wspace=0.31, hspace=0.03)

fig.savefig(f'resolution_convergence.pdf', format='pdf', bbox_inches='tight')

plt.show()
