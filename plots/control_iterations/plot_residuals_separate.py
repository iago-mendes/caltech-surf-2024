# type: ignore
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import sys
from cycler import cycler

colors = ['#4c72b0', '#55a868', '#c44e52', '#8172b3', '#937860', '#da8bc3', '#8c8c8c', '#ccb974', '#64b5cd']
markers = ['o', 's', 'D', '^', 'v', '<', '>', '*', '+']
plt.rcParams.update({
	'text.usetex' : True,
	'font.family' : 'serif',
	'font.serif' : ['Computer Modern Serif'],
	'font.size': 15,
	'axes.prop_cycle': cycler('color', colors) +
                     cycler('marker', markers) +
										 cycler('markerfacecolor', ['none'] * len(colors))
})

dir = sys.argv[1]
fname = f'{dir}/ControlParamsData.txt'
entries = np.genfromtxt(fname, delimiter=',')
if len(entries.shape) < 2:
	print('ERROR: not enough iterations to plot.')
	exit()

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

plot_counter = 0
def plot(data1, data2, label1, label2):
  global plot_counter
  fig, ax = plt.subplots(1, 1, squeeze=True)
  ax.plot(iterations, data1, label=label1, color=colors[plot_counter], marker=markers[plot_counter])
  ax.plot(iterations, data2, label=label2, color=colors[plot_counter+1], marker=markers[plot_counter+1])
  plot_counter += 2
  
  ax.set_xlabel('Control iteration')

  ax.legend()
  ax.set_yscale('log')
  ax.xaxis.set_major_locator(MaxNLocator(integer=True))
  ax.grid('on', linestyle='--', alpha=0.3)

  fig.set_size_inches(7, 4)
  plt.tight_layout()
  fig.savefig(f'{dir}-residuals{plot_counter//2}.pdf', format='pdf', bbox_inches='tight')

  plt.show()


plot(np.abs(residual_M_A), np.abs(residual_M_B), r'$|M_A - M^*_A|$', r'$|M_B - M^*_B|$')
plot(residual_chi_A, residual_chi_B, r'$|\chi_A - \chi^*_A|$', r'$|\chi_B - \chi^*_B|$')
plot(residual_CoM, residual_Padm, r'$|C_{CoM}|$', r'$|P_{ADM}|$')
