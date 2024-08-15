# type: ignore
import numpy as np
import matplotlib.pyplot as plt
import sxs

plt.rcParams.update({
	'text.usetex' : True,
	'font.family' : 'serif',
	'font.serif' : ['Computer Modern Serif'],
	'font.size': 15
})


# simulation = "SXS:BBH:0303" # high mass ratio
simulation = "SXS:BBH:1146" # high spin
metadata = sxs.load(f"{simulation}/Lev/metadata.json")
horizons = sxs.load(f"{simulation}/Lev/Horizons.h5")

fig, ax = plt.subplots(1, 1, squeeze=True)

# ax.plot(horizons.A.time, horizons.A.christodoulou_mass)
# ax.set_ylabel(r'Christodoulou mass of AhA')

ax.plot(horizons.A.time, horizons.A.chi_inertial_mag)
ax.set_ylabel(r'Dimensionless spin of AhA')

ax.axvline(metadata.relaxation_time, ymin=0, ymax=1, linestyle='dashed', color='black')

# ax.set_yscale('log')
ax.grid('on', linestyle='--', alpha=0.3)

ax.set_xlabel(r'Evolution time')

plt.tight_layout()
fig.set_size_inches(10, 5)
# fig.savefig(f'analyze_junk.pdf', format='pdf', bbox_inches='tight')

plt.show()
