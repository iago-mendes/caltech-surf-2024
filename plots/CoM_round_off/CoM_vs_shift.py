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

distance_fixed = 100000
L_fixed = 1
P_fixed = 12

def get_test_residuals(dir):
  residuals = []

  for i, dir in enumerate([dir]):
    print(dir)

    z_shifts = []
    centers_of_mass_x = []
    centers_of_mass_y = []
    centers_of_mass_z = []

    center_of_mass_entries = np.genfromtxt(f'{dir}/Test_CenterOfMass.output', delimiter=',')

    for row in range(len(center_of_mass_entries)):
      distance, L, P, center_of_mass_x, center_of_mass_y, center_of_mass_z, z_shift = center_of_mass_entries[row]

      if distance != distance_fixed or L != L_fixed or P != P_fixed:
        continue
      
      # if center_of_mass_x == 0.0 or center_of_mass_y == 0.0 or center_of_mass_z == 0.0:
      #   print(f'\tSkipped {distance} due to 0.0 value')
      #   continue

      z_shifts.append(z_shift)
      centers_of_mass_x.append(np.abs(center_of_mass_x))
      centers_of_mass_y.append(np.abs(center_of_mass_y))
      centers_of_mass_z.append(np.abs(center_of_mass_z))
      # centers_of_mass_z.append(center_of_mass_z)
    
    # centers_of_mass_z = np.array(centers_of_mass_z)
    # centers_of_mass_z -= centers_of_mass_z[-1]
    # centers_of_mass_z = np.abs(centers_of_mass_z[:-1])
    
    residuals.append({
      'z_shifts': z_shifts,
      'centers_of_mass_x': centers_of_mass_x,
      'centers_of_mass_y': centers_of_mass_y,
      'centers_of_mass_z': centers_of_mass_z,
    })

  return residuals

fig, axes = plt.subplots(3, 1, squeeze=True)
# axes = np.transpose(axes)
# col = 0

def plot(title, residuals):
  global col
  # ax1, ax2, ax3 = axes[col]
  ax1, ax2, ax3 = axes
  # col += 1

  ax1.set_title(title)

  for i, label in enumerate([r'$C_0^z=0.1$']):
    distances = residuals[i]['z_shifts']
    ax1.plot(distances[:-1], residuals[i]['centers_of_mass_x'][:-1], label=label, color=colors[i], marker=markers[i])
    ax2.plot(distances[:-1], residuals[i]['centers_of_mass_y'][:-1], label=label, color=colors[i], marker=markers[i])
    ax3.plot(distances[:-1], residuals[i]['centers_of_mass_z'][:-1], label=label, color=colors[i], marker=markers[i])

plot(
  r'Isotropic [$(\psi^4 - \langle \psi^4 \rangle) n^i$]',
  get_test_residuals('./data/isotropic_varying_shift')
)

for ax in axes.flatten():
  # ax.set_xscale('log')
  # ax.set_yscale('log')
  ax.grid('on', linestyle='--', alpha=0.5)

# # Remove x tick labels
# for ax in axes[:,:-1].flatten():
#   ax.set_xticklabels([])

# # Add x axis labels
# for bottom_ax in np.transpose(axes)[-1]:
#   bottom_ax.set_xlabel(r'$R$')

# # Add y axis labels
# axes[0,0].set_ylabel(r'$\Bigg| C_{CoM}^x \Bigg|$')
# axes[0,1].set_ylabel(r'$\Bigg| C_{CoM}^y \Bigg|$')
# axes[0,2].set_ylabel(r'$\Bigg| C_{CoM}^z \Bigg|$')

# # Add legend to one axis
# axes[0, 0].legend()

fig.set_size_inches(12,10)
plt.tight_layout()
plt.subplots_adjust(wspace=0.275, hspace=0.025)

fig.savefig(f'CoM_vs_shift.pdf', format='pdf', bbox_inches='tight')

plt.show()
