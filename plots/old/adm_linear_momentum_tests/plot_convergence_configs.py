# type: ignore
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
	'text.usetex' : True,
	'font.family' : 'serif',
	'font.serif' : ['Computer Modern Serif'],
	'font.size': 15
})


distances = np.array([
  1000,
  10000,
  100000,
  1000000,
  10000000,
  100000000,
  1000000000,
  10000000000,
])
mass_residuals = {
  'lin': np.array([
    -0.002051201106249456,
    -0.001368369125293922,
    -0.001299893473913905,
    -0.001293057283517918,
    -0.001292370503239271,
    -0.001292258668839974,
    -0.001300861227510053,
    -0.001327328499668212,
  ]),
  'lin-p': np.array([
    -0.002051201106249456,
    -0.000595121473521143,
    -0.0005242796785178694,
    -0.0005171915924093717,
    -0.0005165490058056665,
    -0.0005165398618094397,
    -0.0005253400388858243,
    -0.0005660293022529128,
  ]),
  'log': np.array([
    -0.006119469240669639,
    -0.07481381607408499,
    -0.5913498580421278,
    -3.517940640643124,
    -17.61249484687846,
    -78.52856938416946,
    -322.445185900725,
    -1246.099665617596,
  ]),
  'log-p': np.array([
    -0.006119469240669639,
    -8.569082967602171e-05,
    -1.472099881594247e-05,
    -7.62222384009803e-06,
    -6.92236942168023e-06,
    -7.032048803479185e-06,
    -9.526098139867045e-06,
    8.608495553730222e-06,
  ]),
  'inv': np.array([
    -0.0007875791285105738, 
    -7.963689814460828e-05,
    -8.619163775147598e-06,
    -1.515147296959896e-06,
    -8.047120692022958e-07,
    -7.335467866464995e-07,
    -7.31091952754781e-07,
    -6.832487067232051e-07,
  ]),
  'inv-p': np.array([
    -0.0007875791285105738, 
    -7.891241452329112e-05,
    -7.893445575390956e-06,
    -7.884640711441904e-07,
    -7.91535339494942e-08,
    -3.642880419540973e-08,
    1.589405620450179e-07,
    -2.35270721238301e-06,
  ]),
}
momentum_residuals = {
  'lin': np.array([
    -0.001078374089097278,
    -0.000537663210089967,
    -0.0004833464960790623,
    -0.0004779135052742634,
    -0.0004773749179605158,
    -0.0004773763625428584,
    -0.0004771632129443315,
    -0.0004813759829204178,
  ]),
  'lin-p': np.array([
    -0.001078374089097278,
    -0.0002529241062603216,
    -0.0001969578038379938,
    -0.0001913602133136738,
    -0.0001908052504625557,
    -0.0001907754611085721,
    -0.0001905418010318405,
    -0.0001917313040438184,
  ]),
  'log': np.array([
    0.000814444167387518,
    -0.01727862514734324,
    -0.2003392887269249,
    -1.274708634853444,
    -6.475698696692864,
    -28.97201444935374,
    -119.0615281385565,
    -460.2154342818607,
  ]),
  'log-p': np.array([
    0.000814444167387518,
    -6.513206984504194e-05,
    -8.787458908643409e-06,
    -3.150329750289949e-06,
    -2.587223030925401e-06,
    -2.563572120806157e-06,
    -2.574492598772515e-06,
    -4.099644883348574e-06,
  ]),
  'inv': np.array([
    -0.0007012245182380639,
    -0.0001406544054193271,
    -8.43002696014894e-05,
    -7.866187172178574e-05,
    -7.80979616252786e-05,
    -7.80419612939065e-05,
    -7.803997659283635e-05,
    -7.80241183936381e-05,
  ]),
  'inv-p': np.array([
    -0.0007012245182380639,
    -6.2791228222836e-05,
    -6.282033481408433e-06,
    -6.280989670592874e-07,
    -6.247417472238226e-08,
    -1.128869664412946e-08,
    2.316862046658485e-08,
    -6.846613345157238e-07,
  ]),
}

fig, (ax1, ax2) = plt.subplots(2, 1, squeeze=True)

for config in mass_residuals.keys():
  ax1.plot(distances, np.abs(mass_residuals[config]), label=config, marker='o')
  ax2.plot(distances, np.abs(momentum_residuals[config]), label=config, marker='o')

ax1.plot(distances, 10. / distances, label='test', linestyle='dashed', color='black')
ax2.plot(distances, 10. / distances, label='test', linestyle='dashed', color='black')

ax1.set_xscale('log')
ax2.set_xscale('log')
ax1.set_yscale('log')
ax2.set_yscale('log')

ax1.grid('on', linestyle='--', alpha=0.3)
ax2.grid('on', linestyle='--', alpha=0.3)
ax1.legend(loc='center left', bbox_to_anchor=(1., -0.02))

ax1.set_title('Convergence domain options (with fixed L and P)', pad=10)
ax1.set_ylabel(r'$M_{ADM}$ residual')
ax2.set_ylabel(r'$P_{ADM}^z$ residual')
ax2.set_xlabel('Outer radius')

plt.tight_layout()
plt.subplots_adjust(hspace=0.02)
fig.set_size_inches(10, 8)
fig.savefig(f'convergence_configs.pdf', format='pdf', bbox_inches='tight')

# plt.show()
