# type: ignore
import numpy as np
import matplotlib.pyplot as plt
import sys

plt.rcParams.update({
	'text.usetex' : True,
	'font.family' : 'serif',
	'font.serif' : ['Computer Modern Serif'],
	'font.size': 15
})

data_file = sys.argv[1]
data = np.genfromtxt(data_file, delimiter=',')

distances = np.unique(data[:,0])
markers = ['o', 's', 'D', '^', 'v', 'p', '*', 'h', 'H', '+']
# resolutions = data[:,1]
# adm_masses = data[:,2]
# adm_linear_momenta = data[:,2]
# centers_of_mass = data[:,3]

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, squeeze=True)

for i, R in enumerate(distances):
  resolutions = []
  adm_masses = []
  adm_linear_momenta = []
  centers_of_mass = []

  for row in range(len(data)):
    if data[row,0] != R:
      continue
    resolutions.append(data[row,1])
    adm_masses.append(data[row,2])
    adm_linear_momenta.append(data[row,3])
    centers_of_mass.append(data[row,4])
  
  adm_masses = np.array(adm_masses)
  adm_linear_momenta = np.array(adm_linear_momenta)
  centers_of_mass = np.array(centers_of_mass)

  adm_masses -= adm_masses[-1]
  adm_linear_momenta -= adm_linear_momenta[-1]
  centers_of_mass -= centers_of_mass[-1]
  resolutions = resolutions[:-1]
  adm_masses = np.abs(adm_masses[:-1])
  adm_linear_momenta = np.abs(adm_linear_momenta[:-1])
  centers_of_mass = np.abs(centers_of_mass[:-1])

  # resolutions = resolutions[:-1]
  # adm_masses = np.abs(np.diff(adm_masses))
  # adm_linear_momenta = np.abs(np.diff(adm_linear_momenta))
  # centers_of_mass = np.abs(np.diff(centers_of_mass))

  ax1.plot(resolutions, adm_masses, label=f'$R = 10^{np.log10(R):.0f}$', marker=markers[i])
  ax2.plot(resolutions, adm_linear_momenta, label=f'$R = 10^{np.log10(R):.0f}$', marker=markers[i])
  ax3.plot(resolutions, centers_of_mass, label=f'$R = 10^{np.log10(R):.0f}$', marker=markers[i])

# plt.suptitle(r'Convergence test of ADM integrals (boosted isotropic Schwarzschild)')
ax1.set_title(r'$M_{ADM} - M_{ADM}|_{ref}$')
ax2.set_title(r'$P_{ADM} - P_{ADM}|_{ref}$')
ax3.set_title(r'$C_{CoM} - C_{CoM}|_{ref}$')

for ax in [ax1, ax2, ax3]:
  ax.set_yscale('log')
  ax.grid('on', linestyle='--', alpha=0.3)
  ax.set_xlabel(r'$P$')

ax1.legend()

fig.set_size_inches(16, 6)
plt.tight_layout()
plt.subplots_adjust(wspace=0.15)
fig.savefig(f'{data_file[:-7]}.pdf', format='pdf', bbox_inches='tight')

plt.show()
