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

def get_data(dir):
  print(dir)

  control_data = np.loadtxt(f'{dir}/ControlParamsData.txt', delimiter=',', dtype=float)

  iterations = control_data[:,0].astype(int)
  residuals = []

  for i in iterations:
    residuals.append(np.max(np.abs(control_data[i,1:])))

  return {
    'iterations': iterations,
    'residuals': residuals,
  }

fig, ax = plt.subplots(1, 1, squeeze=True)

def plot(title, args):
  ax.set_title(title)

  for label, data in args:
    ax.plot(data['iterations'], data['residuals'], label=label)

plot("Modified root-finding methods", [
  [r'Delay first 2 iterations of $P_{ADM}$ and $C_{CoM}$', get_data('./data/root_finding/2024-11-11-delay')],
  [r'Diagonal Jacobian', get_data('./data/root_finding/2024-11-11-Jacobian-diagonal')],
  [r'Damp $\Delta u$ with $1-e^{-k}$', get_data('./data/root_finding//2024-11-11-FreeData-damping-expk')],
  [r'Damp $\Delta J$ with $1-e^{-|\Delta F|^2/|\Delta u|^2}$', get_data('./data/root_finding/2024-11-11-Jacobian-damping-expDFDu')],
  [r'Damp $\Delta J$ with $1-e^{-k}$', get_data('./data/root_finding/2024-11-15-Jacobian-damping-expk')],
  [r'Damp $\Delta J$ with $1-e^{-k/2}$', get_data('./data/root_finding/2024-11-15-Jacobian-damping-expko2')],
  [r'Damp $\Delta J$ with $1-e^{-k/4}$', get_data('./data/root_finding/2024-11-15-Jacobian-damping-expko4')],
])

ax.set_xlabel(r'$k$')
ax.legend()

ax.set_xscale('linear')
ax.set_yscale('log')

ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.grid('on', linestyle='--', alpha=0.5)

ax.set_ylabel(r'Control residual')


fig.set_size_inches(12,8)
plt.tight_layout()

fig.savefig(f'root_finding_methods.png', format='png', bbox_inches='tight')

plt.show()
