# type: ignore
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
	'text.usetex' : True,
	'font.family' : 'serif',
	'font.serif' : ['Computer Modern Serif'],
	'font.size': 15
})


fname = 'mass_residuals.output'
entries = np.genfromtxt(fname, delimiter=',')

Ls = np.array(entries[:,0])
Ps = np.array(entries[:,1])
distances = np.array(entries[:,2])
residuals = np.array(entries[:,3])

fig, ax = plt.subplots(1, 1, squeeze=True)

for L in np.unique(Ls):
  for d in np.unique(distances):
    group_Ps = []
    group_residuals = []

    for i in range(len(residuals)):
      if Ls[i] == L and distances[i] == d:
        group_Ps.append(Ps[i])
        group_residuals.append(residuals[i])
  
    ax.plot(group_Ps, np.abs(group_residuals), label=f'L = {L:.0f}, R = {d:.0e}', marker='o')

ax.plot(np.unique(Ps), 200*10**(-np.unique(Ps)), linestyle='dashed', color='black', label='Test')

# ax.legend()
ax.legend(loc='center left', bbox_to_anchor=(1., 0.5))
ax.set_yscale('log')
ax.grid('on', linestyle='--', alpha=0.3)

ax.set_title('ADM Mass for Schwarzschild in Isotropic Coordinates')
ax.set_ylabel(r'$|M_{ADM} - M|$')
ax.set_xlabel('Polynomial order')

plt.tight_layout()
fig.set_size_inches(15, 5)
fig.savefig(f'mass_residuals.pdf', format='pdf', bbox_inches='tight')
