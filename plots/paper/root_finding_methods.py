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

fig, axes = plt.subplots(1, 3, squeeze=True)

col = 0
data_min = np.Infinity
data_max = -np.Infinity

def plot(title, args):
  global col, data_min, data_max
  ax = axes[col]
  col += 1

  ax.set_title(title)

  for label, data in args:
    ax.plot(data['iterations'], data['residuals'], label=label)

    data_min = min(data_min, np.min(data['residuals']))
    data_max = max(data_max, np.max(data['residuals']))

plot("Delayed/damped Broyden's method", [
  [r'Delay first 2 iterations', get_data('./data/root_finding/delay_first_2_its')],
  [r'Damp with $(1-e^{-k})$', get_data('./data/root_finding/damp_initial_iterations')],
  [r'Diagonal', get_data('./data/root_finding/diagonal')],
  [r'{scipy(broyden1)}: alpha = 1', get_data('./data/root_finding/scipy_noline')],
  [r'{scipy(broyden1)}: alpha = 0.25', get_data('./data/root_finding/scipy_noline_damped')],
])
plot(r"\texttt{scipy(broyden1)} with line search", [
  ['line_search = armijo, alpha = 1', get_data('./data/root_finding/scipy')],
  ['line_search = armijo, alpha = 0.25', get_data('./data/root_finding/scipy_damped')],
  ['line_search = wolfe, alpha = 1', get_data('./data/root_finding/scipy_wolfe_alpha1')],
  ['line_search = wolfe, alpha = 0.25', get_data('./data/root_finding/scipy_wolfe')],
])
plot(r"\texttt{scipy(diagbroyden)}", [
  ['line_search = None, alpha = 0.25', get_data('./data/root_finding/scipy_diag_noline_damped')],
  ['line_search = armijo, alpha = 1', get_data('./data/root_finding/scipy_diag')],
  ['line_search = armijo, alpha = 0.25', get_data('./data/root_finding/scipy_diag_damped')],
  ['line_search = wolfe, alpha = 0.25', get_data('./data/root_finding/scipy_diag_damped_wolfe')],
])

for ax in axes:
  ax.set_xlabel(r'$k$')
  ax.legend()

  ax.set_xscale('linear')
  ax.set_yscale('log')

  ax.set_ylim(data_min/1.5, data_max*1.5)

  ax.xaxis.set_major_locator(MaxNLocator(integer=True))
  ax.grid('on', linestyle='--', alpha=0.5)

axes[0].set_ylabel(r'Control residual')

axes[1].set_xlim(-0.5, 20.5)
axes[2].set_xlim(-0.5, 20.5)


fig.set_size_inches(18,8)
plt.tight_layout()
plt.subplots_adjust(wspace=0.02)

# Remove y tick labels
for ax in axes[1:]:
  ax.set_yticklabels([])

fig.savefig(f'root_finding_methods.png', format='png', bbox_inches='tight')

plt.show()
