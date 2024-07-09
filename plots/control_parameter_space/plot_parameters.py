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


iterations = np.array(entries[:,0])
M_A = np.array(entries[:,1])
M_B = np.array(entries[:,2])
chi_A_x = np.array(entries[:,3])
chi_A_y = np.array(entries[:,4])
chi_A_z = np.array(entries[:,5])
chi_B_x = np.array(entries[:,6])
chi_B_y = np.array(entries[:,7])
chi_B_z = np.array(entries[:,8])


if len(iterations) < 2:
	print('Not enough iterations to plot!')
	exit()


fig, (ax1, ax2) = plt.subplots(2, 1, squeeze=True)

ax1.plot(M_A, label=r'$M_A$', linewidth=8)
ax1.plot(M_B, label=r'$M_B$', linewidth=3)

ax2.plot(np.abs(chi_A_x), label=r'$\chi^x_A$', linewidth=8)
ax2.plot(np.abs(chi_A_y), label=r'$\chi^y_A$', linewidth=8)
ax2.plot(np.abs(chi_A_z), label=r'$\chi^z_A$', linewidth=8)

ax2.plot(np.abs(chi_B_x), label=r'$\chi^x_B$', linewidth=3)
ax2.plot(np.abs(chi_B_y), label=r'$\chi^y_B$', linewidth=3)
ax2.plot(np.abs(chi_B_z), label=r'$\chi^z_B$', linewidth=3)

ax1.legend()
ax2.legend()
# ax1.set_yscale('log')
# ax2.set_yscale('log')
ax1.grid('on', linestyle='--', alpha=0.3)
ax2.grid('on', linestyle='--', alpha=0.3)
ax2.set_xlabel('Iteration')

plt.tight_layout()
plt.subplots_adjust(hspace=0.03)
fig.set_size_inches(10, 8)
fig.savefig(f'{test}/{test}-Parameters.pdf', format='pdf', bbox_inches='tight')

# plt.show()
