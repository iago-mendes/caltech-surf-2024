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


fname = 'residuals_vs_distance.output'
entries = np.genfromtxt(fname, delimiter=',')

R = np.array(entries[:,0])
residuals = np.array(entries[:,1])

fig, ax = plt.subplots(1, 1, squeeze=True)

# ax.plot(L_Ps, L_residuals, label=f'Refinement level = {L}', marker='o')
ax.plot(np.log10(R), np.log10(residuals), marker='o')

a, b = np.polyfit(np.log10(R), np.log10(residuals), 1)
ax.plot(np.log10(R), a*np.log10(R)+b, linestyle='dashed', label=f'linear fit: y = {a:.2f}x + {b:.2f}', color='black')

ax.legend()
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.grid('on', linestyle='--', alpha=0.3)
ax.set_ylabel(r'$\log_{10} |P_{ADM}^z - \gamma M v|$')
ax.set_xlabel(r'$\log_{10} R_{outer}$')

plt.tight_layout()
# plt.subplots_adjust(hspace=0.03)
fig.set_size_inches(10, 4)
fig.savefig(f'adm_linear_momentum_residual_vs_distance.pdf', format='pdf', bbox_inches='tight')

plt.show()
