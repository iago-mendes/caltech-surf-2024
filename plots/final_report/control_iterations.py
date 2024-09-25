# type: ignore
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from cycler import cycler

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

fig, axes = plt.subplots(3, 3, squeeze=True)
axes = np.transpose(axes)

data_min = [np.Infinity, np.Infinity, np.Infinity]
data_max = [-np.Infinity, -np.Infinity, -np.Infinity]

col = 0
def plot(title, dir):
  global col
  ax1, ax2, ax3 = axes[col]
  col += 1

  ax1.set_title(title)

  fname = f'{dir}/ControlParamsData.txt'
  entries = np.genfromtxt(fname, delimiter=',')

  iterations = np.array(entries[:,0])
  residual_M_A = np.abs(np.array(entries[:,1 if len(entries[0]) == 15 else 9]))
  residual_M_B = np.abs(np.array(entries[:,2 if len(entries[0]) == 15 else 10]))
  residual_chi_A_x = np.array(entries[:,3 if len(entries[0]) == 15 else 11])
  residual_chi_A_y = np.array(entries[:,4 if len(entries[0]) == 15 else 12])
  residual_chi_A_z = np.array(entries[:,5 if len(entries[0]) == 15 else 13])
  residual_chi_B_x = np.array(entries[:,6 if len(entries[0]) == 15 else 14])
  residual_chi_B_y = np.array(entries[:,7 if len(entries[0]) == 15 else 15])
  residual_chi_B_z = np.array(entries[:,8 if len(entries[0]) == 15 else 16])
  residual_CoM_x = np.array(entries[:,9 if len(entries[0]) == 15 else 17])
  residual_CoM_y = np.array(entries[:,10 if len(entries[0]) == 15 else 18])
  residual_CoM_z = np.array(entries[:,11 if len(entries[0]) == 15 else 19])
  residual_Padm_x = np.array(entries[:,12 if len(entries[0]) == 15 else 20])
  residual_Padm_y = np.array(entries[:,13 if len(entries[0]) == 15 else 21])
  residual_Padm_z = np.array(entries[:,14 if len(entries[0]) == 15 else 22])
  
  residual_chi_A = np.sqrt(residual_chi_A_x**2 + residual_chi_A_y**2 + residual_chi_A_z**2)
  residual_chi_B = np.sqrt(residual_chi_B_x**2 + residual_chi_B_y**2 + residual_chi_B_z**2)
  residual_CoM = np.sqrt(residual_CoM_x**2 + residual_CoM_y**2 + residual_CoM_z**2)
  residual_Padm = np.sqrt(residual_Padm_x**2 + residual_Padm_y**2 + residual_Padm_z**2)

  ax1.plot(iterations, residual_M_A, label=r'$|M_A - M^*_A|$', color=colors[0], marker=markers[0])
  ax1.plot(iterations, residual_M_B, label=r'$|M_B - M^*_B|$', color=colors[1], marker=markers[1])
  ax2.plot(iterations, residual_chi_A, label=r'$|\chi_A - \chi^*_A|$', color=colors[2], marker=markers[2])
  ax2.plot(iterations, residual_chi_B, label=r'$|\chi_B - \chi^*_B|$', color=colors[3], marker=markers[3])
  ax3.plot(iterations, residual_CoM, label=r'$|C_{CoM}|$', color=colors[4], marker=markers[4])
  ax3.plot(iterations, residual_Padm, label=r'$|P_{ADM}|$', color=colors[5], marker=markers[5])

  data_rows = [
    np.concatenate([residual_M_A, residual_M_B]),
    np.concatenate([residual_chi_A, residual_chi_B]),
    np.concatenate([residual_CoM, residual_Padm]),
  ]
  for i, row in enumerate(data_rows):
    data_min[i] = min(data_min[i], np.min(row))
    data_max[i] = max(data_max[i], np.max(row))

plot(r'\texttt{q1}', '../control_iterations/data/q1')
plot(r'\texttt{q3}', '../control_iterations/data/q3')
plot(r'\texttt{q3-spin}', '../control_iterations/data/q3-spin')

for ax in axes.flatten():
  ax.set_yscale('log')
  ax.xaxis.set_major_locator(MaxNLocator(integer=True))
  ax.grid('on', linestyle='--', alpha=0.3)

for i, row in enumerate(np.transpose(axes)):
  for ax in row:
    ax.set_ylim(data_min[i]/2., data_max[i]*2.)

for ax in axes[1:].flatten():
  ax.set_yticklabels([])

for right_ax in axes[-1]:
  right_ax.legend(loc='center left', bbox_to_anchor=(0.99, 0.5))

for bottom_ax in np.transpose(axes)[-1]:
  bottom_ax.set_xlabel(r'$k$')

fig.set_size_inches(12, 8)
plt.tight_layout()
plt.subplots_adjust(wspace=0.03, hspace=0.04)

fig.savefig(f'control_iterations.pdf', format='pdf', bbox_inches='tight')

plt.show()
