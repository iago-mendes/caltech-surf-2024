# type: ignore
import numpy as np
import matplotlib.pyplot as plt
import h5py

plt.rcParams.update({
	'text.usetex' : True,
	'font.family' : 'serif',
	'font.serif' : ['Computer Modern Serif'],
	'font.size': 15
})

all_resolutions = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

resolutions = []
adm_masses = []
adm_linear_momenta = []
centers_of_mass = []

for P in all_resolutions:
  reductions = h5py.File(f'P-{P:02}/BbhReductions.h5')
  if 'AdmIntegrals.dat' in reductions.keys():
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

fig, ax = plt.subplots(1, 1, squeeze=True)

resolutions = np.array(resolutions)
adm_masses = np.array(adm_masses)
adm_linear_momenta = np.array(adm_linear_momenta)
centers_of_mass = np.array(centers_of_mass)

# adm_masses -= adm_masses[-1]
# adm_linear_momenta -= adm_linear_momenta[-1]
# centers_of_mass -= centers_of_mass[-1]

# resolutions = np.abs(resolutions[:-1])
# adm_masses = np.abs(adm_masses[:-1])
# adm_linear_momenta = np.abs(adm_linear_momenta[:-1])
# centers_of_mass = np.abs(centers_of_mass[:-1])

resolutions = resolutions[:-1]
adm_masses = np.abs(np.diff(adm_masses))
adm_linear_momenta = np.abs(np.diff(adm_linear_momenta))
centers_of_mass = np.abs(np.diff(centers_of_mass))

ax.plot(resolutions, adm_masses, label=r'$\Delta M_{ADM}$', marker='o')
ax.plot(resolutions, adm_linear_momenta, label=r'$\Delta P_{ADM}$', marker='o')
ax.plot(resolutions, centers_of_mass, label=r'$\Delta CoM$', marker='o')

ax.set_title(r'Convergence test of ADM integrals with resolution')
ax.set_xlabel(r'Polynomial order $P$')
ax.set_ylabel(r'Difference with highest resolution')

ax.set_yscale('log')

ax.grid('on', linestyle='--', alpha=0.3)
ax.legend()

plt.tight_layout()
fig.set_size_inches(10, 6)

plt.show()
