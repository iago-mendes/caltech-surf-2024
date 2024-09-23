# type: ignore
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import h5py

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

distance_exps = [2, 3, 4, 5, 6, 7, 8]
L_fixed = 1
P_fixed = 12

def get_bbh_residuals(surface_dir, volume_dir):
  residuals = []

  for i, dir in enumerate([surface_dir, volume_dir]):
    print(dir) 

    distances = []
    adm_masses = []
    adm_momenta = []

    for R in distance_exps:
      subdir = f'R{R}L{L_fixed}P{P_fixed:02}'

      try:
        reductions = h5py.File(f'{dir}/{subdir}/BbhReductions.h5')

        adm_mass = reductions['AdmIntegrals.dat'][0,0]

        adm_linear_momentum_x = reductions['AdmIntegrals.dat'][0,1]
        adm_linear_momentum_y = reductions['AdmIntegrals.dat'][0,2]
        adm_linear_momentum_z = reductions['AdmIntegrals.dat'][0,3]
        adm_linear_momentum = np.sqrt(adm_linear_momentum_x**2 + adm_linear_momentum_y**2 + adm_linear_momentum_z**2)

        distances.append(float(f'1.e{R}'))
        adm_masses.append(adm_mass)
        adm_momenta.append(adm_linear_momentum)

      except Exception as e:
        print(f'\tSkipped {subdir} due to error: {e}')
    
    adm_masses = np.array(adm_masses)
    adm_momenta = np.array(adm_momenta)

    adm_masses = np.abs(adm_masses - adm_masses[-1])
    adm_momenta = np.abs(adm_momenta - 0.0)
      
    residuals.append({
      'adm_mass_distances': distances,
      'adm_masses': adm_masses,
      'adm_momentum_distances': distances,
      'adm_momenta': adm_momenta,
    })

  return residuals

def get_test_residuals(surface_dir, volume_dir):
  residuals = []

  for i, dir in enumerate([surface_dir, volume_dir]):
    print(dir)

    adm_mass_distances = []
    adm_masses = []
    adm_momentum_distances = []
    adm_momenta = []

    adm_mass_entries = np.genfromtxt(f'{dir}/Test_AdmMass.output', delimiter=',')
    adm_linear_momentum_entries = np.genfromtxt(f'{dir}/Test_AdmLinearMomentum.output', delimiter=',')

    for row in range(len(adm_mass_entries)):
      distance, L, P, adm_mass, expected_adm_mass = adm_mass_entries[row]
      _, _, _, adm_linear_momentum, expected_adm_linear_momentum = adm_linear_momentum_entries[row]

      if L != L_fixed or P != P_fixed:
        continue
      
      # if adm_mass - expected_adm_mass == 0.0:
      #   print(f'\tSkipped {distance} (Madm) due to 0.0 value')
      # else:
      adm_mass_distances.append(distance)
      adm_masses.append(np.abs(adm_mass - expected_adm_mass))

      # if adm_linear_momentum - expected_adm_linear_momentum == 0.0:
      #   print(f'\tSkipped {distance} (Padm) due to 0.0 value')
      # else:
      adm_momentum_distances.append(distance)
      adm_momenta.append(np.abs(adm_linear_momentum - expected_adm_linear_momentum))
    
    residuals.append({
      'adm_mass_distances': adm_mass_distances,
      'adm_masses': adm_masses,
      'adm_momentum_distances': adm_momentum_distances,
      'adm_momenta': adm_momenta,
    })

  return residuals

fig, axes = plt.subplots(2, 5, squeeze=True)
axes = np.transpose(axes)
col = 0

def plot(title, residuals):
  global col
  ax1, ax2 = axes[col]
  col += 1

  ax1.set_title(title)

  for i, label in enumerate([r'surface', r'volume']):
    ax1.plot(residuals[i]['adm_mass_distances'], residuals[i]['adm_masses'], label=label, color=colors[i], marker=markers[i])
    ax2.plot(residuals[i]['adm_momentum_distances'], residuals[i]['adm_momenta'], label=label, color=colors[i], marker=markers[i])

plot(
  'Isotropic',
  get_test_residuals('../adm_integrals/data/isotropic/surface', '../adm_integrals/data/isotropic/volume')
)
plot(
  'Kerr-Schild',
  get_test_residuals('../adm_integrals/data/KerrSchild/surface', '../adm_integrals/data/KerrSchild/volume')
)
plot(
  'Boosted Isotropic',
  get_test_residuals('../adm_integrals/data/isotropic_boosted/surface', '../adm_integrals/data/isotropic_boosted/volume')
)
plot(
  'Boosted Kerr-Schild',
  get_test_residuals('../adm_integrals/data/KerrSchild_boosted/surface', '../adm_integrals/data/KerrSchild_boosted/volume')
)
plot(
  'BBH',
  get_bbh_residuals('../adm_integrals/data/bbh/surface', '../adm_integrals/data/bbh/volume')
)

for ax in axes.flatten():
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.grid('on', linestyle='--', alpha=0.5)

# Adjust specific axes
# Isotropic Madm
axes[0,0].set_ylim(2.e-17, 8.e-14)
# Isotropic Padm
axes[0,1].plot([1.e2, 1.e8], [1.e-15, 1.e-15], alpha=0)
axes[0,1].set_yscale('symlog', linthresh=1e-17)
# Kerr-Schild Padm
axes[1,1].set_ylim(1.2e-18, 1.5e-15)
# Boosted Isotropic Padm
axes[2,1].set_ylim(2.e-7, 2.e+1)
axes[2,1].set_yticks([1.e+1, 1.e-1, 1.e-3, 1.e-5])
# BBH Madm
axes[4,0].set_yscale('symlog', linthresh=1e-18)
axes[4,0].set_yticks([1.e-4, 1.e-8, 1.e-12, 1.e-16, 0])
# BBH Padm
axes[4,1].set_ylim(1.2e-13, 1.4e-10)

# Remove x tick labels
for ax in axes[:,:-1].flatten():
  ax.set_xticklabels([])

# Add x axis labels
for bottom_ax in np.transpose(axes)[-1]:
  bottom_ax.set_xlabel(r'$R$')

# Add y axis labels
axes[0,0].set_ylabel(r'$\Bigg|M_{ADM} - M_{ref}\Bigg|$')
axes[0,1].set_ylabel(r'$\Bigg|P_{ADM} - P_{ref}\Bigg|$')

# Add legend to one axis
axes[0, 0].legend()

fig.set_size_inches(14,7)
plt.tight_layout()
plt.subplots_adjust(wspace=0.31, hspace=0.03)

plt.show()

fig.savefig(f'distance_convergence_Madm_Padm.pdf', format='pdf', bbox_inches='tight')
