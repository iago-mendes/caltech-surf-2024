import numpy as np
import matplotlib.pyplot as plt # type: ignore
import sys

plt.rcParams.update({
	'text.usetex' : True,
	'font.family' : 'serif',
	'font.serif' : ['Computer Modern Serif'],
	'font.size': 15
})


test = sys.argv[1]
fname = f'{test}/ControlParamsData.txt'
entries = np.genfromtxt(fname, delimiter=',')
if len(entries.shape) < 2:
	print('ERROR: not enough iterations to plot.')
	exit()


iterations = np.array(entries[:,0])
residual_M_A = np.array(entries[:,9])
residual_M_B = np.array(entries[:,10])
residual_chi_A_x = np.array(entries[:,11])
residual_chi_A_y = np.array(entries[:,12])
residual_chi_A_z = np.array(entries[:,13])
residual_chi_B_x = np.array(entries[:,14])
residual_chi_B_y = np.array(entries[:,15])
residual_chi_B_z = np.array(entries[:,16])


def plotAndFit(data, ax, label, linewidth=1):
	ax.plot(data, label=label, linewidth=linewidth)
	a, b = np.polyfit(iterations, data, 1)
	ax.plot(iterations, a*iterations+b, linestyle='dashed', label=f'linear fit: y = {a:.2f}x + {b:.2f}', color='black')


fig, (ax1, ax2) = plt.subplots(2, 1, squeeze=True)

plotAndFit(np.log10(np.abs(residual_M_A)), ax1, label=r'$\log_{10}|M_A - M^*_A|$', linewidth=8)
plotAndFit(np.log10(np.abs(residual_M_B)), ax1, label=r'$\log_{10}|M_B - M^*_B|$', linewidth=3)
# ax1.plot(np.abs(residual_M_A), label=r'$|M_A - M^*_A|$', linewidth=8)
# ax1.plot(np.abs(residual_M_B), label=r'$|M_B - M^*_B|$', linewidth=3)
# ax1.set_yscale('log')

ax2.plot(np.abs(residual_chi_A_x), label=r'$|\chi^x_A - \chi^{*x}_A|$', linewidth=8)
ax2.plot(np.abs(residual_chi_A_y), label=r'$|\chi^y_A - \chi^{*y}_A|$', linewidth=8)
ax2.plot(np.abs(residual_chi_A_z), label=r'$|\chi^z_A - \chi^{*z}_A|$', linewidth=8)

ax2.plot(np.abs(residual_chi_B_x), label=r'$|\chi^x_B - \chi^{*x}_B|$', linewidth=3)
ax2.plot(np.abs(residual_chi_B_y), label=r'$|\chi^y_B - \chi^{*y}_B|$', linewidth=3)
ax2.plot(np.abs(residual_chi_B_z), label=r'$|\chi^z_B - \chi^{*z}_B|$', linewidth=3)

ax1.legend()
ax2.legend()
ax2.set_yscale('log')
ax1.grid('on', linestyle='--', alpha=0.3)
ax2.grid('on', linestyle='--', alpha=0.3)
ax2.set_xlabel('Iteration')

plt.tight_layout()
plt.subplots_adjust(hspace=0.03)
fig.set_size_inches(10, 8)
fig.savefig(f'{test}/{test}-Residuals.pdf', format='pdf', bbox_inches='tight')

# plt.show()
