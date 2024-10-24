# type: ignore
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from cycler import cycler

colors = [
  '#4c72b0', '#55a868', '#c44e52', '#8172b3', '#937860', '#da8bc3', '#8c8c8c',
  '#ccb974', '#64b5cd'
]
markers = ['o', '^', 's', 'v', 'D', '<', '*', '>', '+']
markersizes = [7] * len(markers)
markersizes[6] = 10
plt.rcParams.update({
	'text.usetex' : True,
	'font.family' : 'serif',
	'font.serif' : ['Computer Modern Serif'],
	'font.size': 15,
	'axes.prop_cycle': cycler('color', colors) +
                     cycler('marker', markers) +
										 cycler('markerfacecolor', ['none'] * len(colors)) +
										 cycler('markersize', markersizes) 
})

# title = r'L=0, P=9 (60000 points, $N^{1/3}\approx39$)'
# fname_ext = 'P09'
# dir = './data/root_finding/damping-P09'

# title = r'L=0, P=9 (60000 points, $N^{1/3}\approx39$)'
# fname_ext = 'P9-no1'
# dir = './data/root_finding/damping-mid_resolution-no1'

title = r'L=0, P=12 (128778 points, $N^{1/3}\approx51$)'
fname_ext = 'P12'
dir = './data/root_finding/damping-P12'

# title = r'L=0, P=15 (236544 points, $N^{1/3}\approx62$)'
# fname_ext = 'P15'
# dir = './data/root_finding/damping-lower_tolerance'

# title = r'L=1, P=6'
# fname_ext = 'L1P6'
# dir = './data/root_finding/diagonal'

# title = r'delay'
# fname_ext = 'delay'
# dir = './data/root_finding/delay_first_2_its'

# title = r'L=1, P=8'
# fname_ext = 'L1P8'
# dir = './data/root_finding/delay-L1P8'

control_data = np.loadtxt(f'{dir}/ControlParamsData.txt', delimiter=',', dtype=float)

iterations = control_data[:,0].astype(int)
residuals = np.max(np.abs(control_data[:,1:]), axis=1)
residuals_mass_a = np.abs(control_data[:,1])
residuals_mass_b = np.abs(control_data[:,2])
residuals_chi_a_x = np.abs(control_data[:,3])
residuals_chi_a_y = np.abs(control_data[:,4])
residuals_chi_a_z = np.abs(control_data[:,5])
residuals_chi_b_x = np.abs(control_data[:,6])
residuals_chi_b_y = np.abs(control_data[:,7])
residuals_chi_b_z = np.abs(control_data[:,8])
residuals_CoM_x = np.abs(control_data[:,9])
residuals_CoM_y = np.abs(control_data[:,10])
residuals_CoM_z = np.abs(control_data[:,11])
residuals_Padm_x = np.abs(control_data[:,12])
residuals_Padm_y = np.abs(control_data[:,13])
residuals_Padm_z = np.abs(control_data[:,14])

data_min = np.min(np.abs(control_data[:,1:].flatten()))
data_max = np.max(np.abs(control_data[:,1:].flatten()))

fig, axes = plt.subplots(2, 3, squeeze=True)

axes[0,0].plot(iterations, residuals, label=r'Control residual')

axes[0,1].plot(iterations, residuals_mass_a, label=r'$|M_A-M_A^*|$')
axes[0,1].plot(iterations, residuals_mass_b, label=r'$|M_B-M_B^*|$')

axes[0,2].plot(iterations, residuals_chi_a_x, label=r'$|\chi_A^x-\chi_A^{x,*}|$')
axes[0,2].plot(iterations, residuals_chi_a_y, label=r'$|\chi_A^y-\chi_A^{y,*}|$')
axes[0,2].plot(iterations, residuals_chi_a_z, label=r'$|\chi_A^z-\chi_A^{z,*}|$')

axes[1,0].plot(iterations, residuals_chi_b_x, label=r'$|\chi_B^x-\chi_B^{x,*}|$')
axes[1,0].plot(iterations, residuals_chi_b_y, label=r'$|\chi_B^y-\chi_B^{y,*}|$')
axes[1,0].plot(iterations, residuals_chi_b_z, label=r'$|\chi_B^z-\chi_B^{z,*}|$')

axes[1,1].plot(iterations, residuals_Padm_x, label=r'$|P_{ADM}^x|$')
axes[1,1].plot(iterations, residuals_Padm_y, label=r'$|P_{ADM}^y|$')
axes[1,1].plot(iterations, residuals_Padm_z, label=r'$|P_{ADM}^z|$')

axes[1,2].plot(iterations, residuals_CoM_x, label=r'$|C_{CoM}^x|$')
axes[1,2].plot(iterations, residuals_CoM_y, label=r'$|C_{CoM}^y|$')
axes[1,2].plot(iterations, residuals_CoM_z, label=r'$|C_{CoM}^z|$')

plt.suptitle(title)

for ax in axes.flatten():
  ax.set_xlabel(r'$k$')
  ax.legend(loc='upper right')

  ax.set_xscale('linear')
  ax.set_yscale('log')

  ax.set_ylim(data_min/2., data_max*2.)

  ax.xaxis.set_major_locator(MaxNLocator(integer=True))
  ax.grid('on', linestyle='--', alpha=0.5)

# Remove x tick labels
for ax in axes[:-1,:].flatten():
  ax.set_xticklabels([])

# Remove y tick labels
for ax in axes[:,1:].flatten():
  ax.set_yticklabels([])

fig.set_size_inches(16,8)
plt.tight_layout()
plt.subplots_adjust(wspace=0.01, hspace=0.02)

fig.savefig(f'root_finding-{fname_ext}.png', format='png', bbox_inches='tight')

plt.show()
