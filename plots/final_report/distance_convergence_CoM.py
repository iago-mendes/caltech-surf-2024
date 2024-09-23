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

def get_bbh_residuals(unshifted_dir, shifted_dir):
  residuals = []

  for i, dir in enumerate([unshifted_dir, shifted_dir]):
    print(dir)
    z_shift = [0.0, 0.1][i]

    distances = []
    centers_of_mass_x = []
    centers_of_mass_y = []
    centers_of_mass_z = []

    for R in distance_exps:
      subdir = f'R{R}L{L_fixed}P{P_fixed:02}'

      try:
        reductions = h5py.File(f'{dir}/{subdir}/BbhReductions.h5')

        center_of_mass_x = reductions['AdmIntegrals.dat'][0,4]
        center_of_mass_y = reductions['AdmIntegrals.dat'][0,5]
        center_of_mass_z = reductions['AdmIntegrals.dat'][0,6]

        if center_of_mass_x == 0.0 or center_of_mass_y == 0.0 or center_of_mass_z - z_shift == 0.0:
          print(f'\tSkipped {subdir} due to 0.0 value')
          continue

        distances.append(float(f'1.e{R}'))
        centers_of_mass_x.append(np.abs(center_of_mass_x))
        centers_of_mass_y.append(np.abs(center_of_mass_y))
        centers_of_mass_z.append(np.abs(center_of_mass_z - z_shift))

      except Exception as e:
        print(f'\tSkipped {subdir} due to error: {e}')
      
    residuals.append({
      'distances': distances,
      'centers_of_mass_x': centers_of_mass_x,
      'centers_of_mass_y': centers_of_mass_y,
      'centers_of_mass_z': centers_of_mass_z,
    })

  return residuals

def get_test_residuals(unshifted_dir, shifted_dir):
  residuals = []

  for i, dir in enumerate([unshifted_dir, shifted_dir]):
    print(dir)
    z_shift = [0.0, 0.1][i]

    distances = []
    centers_of_mass_x = []
    centers_of_mass_y = []
    centers_of_mass_z = []

    center_of_mass_entries = np.genfromtxt(f'{dir}/Test_CenterOfMass.output', delimiter=',')

    for row in range(len(center_of_mass_entries)):
      distance, L, P, center_of_mass_x, center_of_mass_y, center_of_mass_z, _ = center_of_mass_entries[row]

      if L != L_fixed or P != P_fixed:
        continue
      
      if center_of_mass_x == 0.0 or center_of_mass_y == 0.0 or center_of_mass_z - z_shift == 0.0:
        print(f'\tSkipped {distance} due to 0.0 value')
        continue

      distances.append(distance)
      centers_of_mass_x.append(np.abs(center_of_mass_x))
      centers_of_mass_y.append(np.abs(center_of_mass_y))
      centers_of_mass_z.append(np.abs(center_of_mass_z - z_shift))
    
    residuals.append({
      'distances': distances,
      'centers_of_mass_x': centers_of_mass_x,
      'centers_of_mass_y': centers_of_mass_y,
      'centers_of_mass_z': centers_of_mass_z,
    })

  return residuals

fig, axes = plt.subplots(3, 4, squeeze=True)
axes = np.transpose(axes)
col = 0

def plot(title, residuals):
  global col
  ax1, ax2, ax3 = axes[col]
  col += 1

  ax1.set_title(title)

  for i, label in enumerate([r'$C_0^z=0$', r'$C_0^z=0.1$']):
    distances = residuals[i]['distances']
    ax1.plot(distances, residuals[i]['centers_of_mass_x'], label=label, color=colors[i], marker=markers[i])
    ax2.plot(distances, residuals[i]['centers_of_mass_y'], label=label, color=colors[i], marker=markers[i])
    ax3.plot(distances, residuals[i]['centers_of_mass_z'], label=label, color=colors[i], marker=markers[i])

plot(
  'Isotropic (surface)',
  get_test_residuals('../adm_integrals/data/isotropic/surface', '../adm_integrals/data/isotropic_shifted/surface')
)
plot(
  'Isotropic (volume)',
  get_test_residuals('../adm_integrals/data/isotropic/volume', '../adm_integrals/data/isotropic_shifted/volume')
)
plot(
  'BBH (surface)',
  get_bbh_residuals('../adm_integrals/data/bbh/surface', '../adm_integrals/data/bbh_shifted/surface')
)
plot(
  'BBH (volume)',
  get_bbh_residuals('../adm_integrals/data/bbh/volume', '../adm_integrals/data/bbh_shifted/volume')
)

for ax in axes.flatten():
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.grid('on', linestyle='--', alpha=0.5)

# Remove x tick labels
for ax in axes[:,:-1].flatten():
  ax.set_xticklabels([])

# Add x axis labels
for bottom_ax in np.transpose(axes)[-1]:
  bottom_ax.set_xlabel(r'$R$')

# Add y axis labels
axes[0,0].set_ylabel(r'$\Bigg| C_{CoM}^x \Bigg|$')
axes[0,1].set_ylabel(r'$\Bigg| C_{CoM}^y \Bigg|$')
axes[0,2].set_ylabel(r'$\Bigg| C_{CoM}^z - C_0^z \Bigg|$')

# Add legend to one axis
axes[0, 0].legend()

fig.set_size_inches(12,10)
plt.tight_layout()
plt.subplots_adjust(wspace=0.275, hspace=0.025)

plt.show()

fig.savefig(f'distance_convergence_CoM.pdf', format='pdf', bbox_inches='tight')
