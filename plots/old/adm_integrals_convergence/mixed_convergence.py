# type: ignore
import numpy as np
import matplotlib.pyplot as plt
import h5py
import os

plt.rcParams.update({
	'text.usetex' : True,
	'font.family' : 'serif',
	'font.serif' : ['Computer Modern Serif'],
	'font.size': 15
})

# root_dir = './L1'
# fname = 'mixed_convergence_L1'
# dir_distances = [2, 3, 5]
# dir_resolutions = [2, 3, 4, 5, 6, 7, 8, 9, 10]

# root_dir = './L2'
# fname = 'mixed_convergence_L2'
# dir_distances = [2, 3, 5]
# dir_resolutions = [2, 3, 4, 5, 6, 7, 8]

root_dir = './lower_tolerance'
fname = 'mixed_convergence_lower_tolerance'
dir_distances = [2, 3, 5]
dir_resolutions = [2, 3, 4, 5, 6, 7, 8, 9, 10]

# root_dir = './volume'
# fname = 'mixed_convergence_volume'
# dir_distances = [2, 3, 5]
# dir_resolutions = [2, 3, 4, 5, 6, 7, 8, 9, 10]

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, squeeze=True)

# dir_distances = [5]
# dir_resolutions = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]

# dir_distances = [3, 5]
# dir_resolutions = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

# dir_distances = [3, 5, 7]
# dir_resolutions = [2, 3, 4, 5, 6, 7, 8, 9, 10]

# dir_distances = [3, 5, 7, 9]
# dir_resolutions = [4, 5, 6, 7, 8]


markers = ['o','^','s', 'D']
# ['o', 's', 'D', '^', 'v', 'p', '*', 'h', 'H', '+']

for i, R in enumerate(dir_distances):
  resolutions = []
  adm_masses = []
  adm_linear_momenta = []
  centers_of_mass = []

  for P in dir_resolutions:
    folder = ''
    # if (R == 5):
    #   folder = f'P-{P:02}'
    # elif (P == 6):
    #   folder = f'R-{R:02}'
    # else:
      # folder = f'R-{R:02}-P-{P:02}'
    folder = f'R{R:02}-P{P:02}'
    # if not os.path.isdir(folder):
    #   continue

    reductions = h5py.File(f'{root_dir}/{folder}/BbhReductions.h5')
    # if 'AdmIntegrals.dat' not in reductions.keys():
    #   continue

    print(f'Using {folder}')

    adm_mass = reductions['AdmIntegrals.dat'][0,0]

    adm_linear_momentum_x = reductions['AdmIntegrals.dat'][0,1]
    adm_linear_momentum_y = reductions['AdmIntegrals.dat'][0,2]
    adm_linear_momentum_z = reductions['AdmIntegrals.dat'][0,3]
    adm_linear_momentum = np.sqrt(adm_linear_momentum_x**2 + adm_linear_momentum_y**2 + adm_linear_momentum_z**2)

    center_of_mass_x = reductions['AdmIntegrals.dat'][0,4]
    center_of_mass_y = reductions['AdmIntegrals.dat'][0,5]
    center_of_mass_z = reductions['AdmIntegrals.dat'][0,6]
    center_of_mass = np.sqrt(center_of_mass_x**2 + center_of_mass_y**2 + center_of_mass_z**2)

    resolutions.append(P)
    adm_masses.append(adm_mass)
    adm_linear_momenta.append(adm_linear_momentum)
    centers_of_mass.append(center_of_mass)

  resolutions = np.array(resolutions)
  adm_masses = np.array(adm_masses)
  adm_linear_momenta = np.array(adm_linear_momenta)
  centers_of_mass = np.array(centers_of_mass)

  adm_masses -= adm_masses[-1]
  adm_linear_momenta -= adm_linear_momenta[-1]
  centers_of_mass -= centers_of_mass[-1]

  resolutions = np.abs(resolutions[:-1])
  adm_masses = np.abs(adm_masses[:-1])
  adm_linear_momenta = np.abs(adm_linear_momenta[:-1])
  centers_of_mass = np.abs(centers_of_mass[:-1])

  # resolutions = resolutions[1:]
  # adm_masses = np.abs(np.diff(adm_masses))
  # adm_linear_momenta = np.abs(np.diff(adm_linear_momenta))
  # centers_of_mass = np.abs(np.diff(centers_of_mass))

  ax1.plot(resolutions, adm_masses, label=f'$R = 10^{R}$', marker=markers[i])
  ax2.plot(resolutions, adm_linear_momenta, label=f'$R = 10^{R}$', marker=markers[i])
  ax3.plot(resolutions, centers_of_mass, label=f'$R = 10^{R}$', marker=markers[i])

# plt.suptitle(r'Convergence test of ADM integrals')
# ax1.set_title(r'$\Delta M_{ADM}$')
# ax2.set_title(r'$\Delta P_{ADM}$')
# ax3.set_title(r'$\Delta CoM$')
ax1.set_title(r'$M_{ADM} - M_{ADM}|_{ref}$')
ax2.set_title(r'$P_{ADM} - P_{ADM}|_{ref}$')
ax3.set_title(r'$C_{CoM} - C_{CoM}|_{ref}$')

# ax.set_xscale('log')
# ax.set_yscale('log')

for ax in [ax1, ax2, ax3]:
  ax.set_yscale('log')
  ax.grid('on', linestyle='--', alpha=0.3)
  ax.set_xlabel(r'$P$')

ax1.legend()

fig.set_size_inches(16, 6)
plt.tight_layout()
plt.subplots_adjust(wspace=0.15)
fig.savefig(f'{fname}.pdf', format='pdf', bbox_inches='tight')

plt.show()
