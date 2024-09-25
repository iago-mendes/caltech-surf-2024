# type: ignore
import numpy as np
import matplotlib.pyplot as plt
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

L_fixed = 1
P_fixed = 12

def get_test_values(dir):
  print(dir)

  distances = []
  avgs_conformal_factor_4 = []

  entries = np.genfromtxt(f'{dir}/Test_CenterOfMass-AvgConformalFactor4.output', delimiter=',')

  for row in range(len(entries)):
    distance, L, P, avg_conformal_factor_4 = entries[row]

    if L != L_fixed or P != P_fixed:
      continue
    
    distances.append(distance)
    avgs_conformal_factor_4.append((avg_conformal_factor_4 - 1.)/2.)
    
  return {
    'distances': distances,
    'avgs_conformal_factor_4': avgs_conformal_factor_4,
  }

fig, ax = plt.subplots(1, 1, squeeze=True)

def plot(label, values):
  ax.plot(values['distances'], values['avgs_conformal_factor_4'], label=label)

plot(r'$C_0^z=0$', get_test_values('./data/isotropic'))
plot(r'$C_0^z=0.1$', get_test_values('./data/isotropic_shifted'))

ax.set_title('Isotropic')
ax.set_xlabel(r'$R$')
ax.set_ylabel(r'$2(\langle \psi^4 \rangle - 1)$')

ax.set_xscale('log')
ax.set_yscale('log')
ax.grid('on', linestyle='--', alpha=0.5)
ax.legend()

fig.set_size_inches(10,6)
plt.tight_layout()
plt.subplots_adjust(wspace=0.275, hspace=0.025)

fig.savefig(f'avg_conformal_factor_4.pdf', format='pdf', bbox_inches='tight')

plt.show()
