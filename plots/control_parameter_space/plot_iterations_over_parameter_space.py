# type: ignore
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.markers import MarkerStyle
import sys

plt.rcParams.update({
	'text.usetex' : True,
	'font.family' : 'serif',
	'font.serif' : ['Computer Modern Serif'],
	'font.size': 15
})


# test = sys.argv[1]
# fname = f'{test}/ControlParamsData.txt'
# entries = np.genfromtxt(fname, delimiter=',')


# iterations = np.array(entries[:,0])
# M_A = np.array(entries[:,1])
# M_B = np.array(entries[:,2])
# chi_A_x = np.array(entries[:,3])
# chi_A_y = np.array(entries[:,4])
# chi_A_z = np.array(entries[:,5])
# chi_B_x = np.array(entries[:,6])
# chi_B_y = np.array(entries[:,7])
# chi_B_z = np.array(entries[:,8])


# if len(iterations) < 2:
# 	print('Not enough iterations to plot!')
# 	exit()


fig, ax = plt.subplots(1, 1, squeeze=True)

# ax1.plot(M_A, label=r'$M_A$', linewidth=8)
# ax1.plot(M_B, label=r'$M_B$', linewidth=3)

# ax2.plot(np.abs(chi_A_x), label=r'$\chi_{A,x}$', linewidth=8)
# ax2.plot(np.abs(chi_A_y), label=r'$\chi_{A,y}$', linewidth=8)
# ax2.plot(np.abs(chi_A_z), label=r'$\chi_{A,z}$', linewidth=8)

# ax2.plot(np.abs(chi_B_x), label=r'$\chi_{B,x}$', linewidth=3)
# ax2.plot(np.abs(chi_B_y), label=r'$\chi_{B,y}$', linewidth=3)
# ax2.plot(np.abs(chi_B_z), label=r'$\chi_{B,z}$', linewidth=3)

# ax1.legend()
# ax2.legend()
# # ax1.set_yscale('log')
# # ax2.set_yscale('log')
# ax1.grid('on', linestyle='--', alpha=0.3)
# ax2.grid('on', linestyle='--', alpha=0.3)
# ax2.set_xlabel('Iteration')

# plt.tight_layout()
# plt.subplots_adjust(hspace=0.03)
# fig.set_size_inches(10, 8)
# fig.savefig(f'{test}/{test}-Parameters.pdf', format='pdf', bbox_inches='tight')

# plt.show()

# Create sample data
x = np.array([1, 2, 3])
y = np.array([2, 4, 1])
values = np.array([5, 8, 3])

q_values = [1, 2, 3, 4, 5]
chi_values = [0.0, 0.2, 0.5, 0.7, 0.9]

points_success_q = []
points_success_chi = []
points_success_num_of_iterations = []

points_fail_q = []
points_fail_chi = []

for q in q_values:
  for chi in chi_values:
    fname = f'q{q}chi{chi:.1f}/ControlParamsData.txt'
    entries = np.genfromtxt(fname, delimiter=',')

    if len(np.shape(entries)) != 2:
      points_fail_q.append(q)
      points_fail_chi.append(chi)
    else:
      points_success_q.append(q)
      points_success_chi.append(chi)
      iterations = np.array(entries[:,0])
      points_success_num_of_iterations.append(iterations[-1])

# q_mesh, chi_mesh = np.meshgrid(q_values, chi_values)

# Create a scatter plot
ax.scatter(points_success_q,
           points_success_chi,
           color='green',
           s=300,
           marker=MarkerStyle('o', fillstyle='none'),
           label=r'Number of control iterations needed to reach a residual $< 10^{-6}$',
)
ax.scatter(points_fail_q, points_fail_chi, color='red', s=300, marker='x')


# Annotate each point with its corresponding value
for i, num_of_iterations in enumerate(points_success_num_of_iterations):
  ax.annotate(f'{num_of_iterations:.0f}',
               (points_success_q[i], points_success_chi[i]),
               (points_success_q[i]-0.03, points_success_chi[i]-0.02),
  )


# ax.set_title(r'Number of control iterations needed to reach a residual of $10^{-6}$')
ax.set_xlabel(r'Mass ratio $q^*=M^*_A/M^*_B$')
ax.set_ylabel(r'Dimensionless spin $\chi^* = \chi^{*z}_A = \chi^{*z}_B$')

ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15))
ax.grid('on', linestyle='--', alpha=0.3)
ax.set_xticks(q_values)
ax.set_yticks(chi_values)

plt.tight_layout()
fig.set_size_inches(10, 5)
fig.savefig(f'control_iterations_over_parameter_space.pdf', format='pdf', bbox_inches='tight')

plt.show()
