# type: ignore
import numpy as np
import matplotlib.pyplot as plt

data = {
  'L':            np.array([                     0,                      0,                      0,                      0,                      1,                      1,                      1,]),
  'P':            np.array([                     3,                      6,                      9,                     12,                      3,                      6,                      9,]),
  'rest_mass':    np.array([2.0772534027862068e+00, 1.9851106019853038e+00, 1.9849708645663506e+00, 1.9849685179498329e+00, 2.9798804409747697e+00, 2.9774506653351689e+00, 2.9774683948080769e+00,]),
  'boosted_mass': np.array([2.3887172550703495e+00, 2.2639204074634485e+00, 2.2591981494020787e+00, 2.2586624542170357e+00, 3.3983479114678383e+00, 3.3883058064285896e+00, 3.3879244204496421e+00,]),
  'momentum':     np.array([1.0924732866020506e+00, 1.0665994424872418e+00, 1.0645233642262937e+00, 1.0642814533022198e+00, 1.5997675706213774e+00, 1.5965556602257960e+00, 1.5963894555312503e+00,]),
}

expected_speed = 0.5
expected_lorentz_factor = 1. / np.sqrt(1 - expected_speed**2)

lorentz_factors_from_mass = data['boosted_mass'] / data['rest_mass']
mass_residuals = np.abs(lorentz_factors_from_mass - expected_lorentz_factor)
print(mass_residuals)
plt.plot(data['P'][-3:], mass_residuals[-3:])

speeds_from_momentum = data['momentum'] / data['boosted_mass']
momentum_residuals = np.abs(speeds_from_momentum - expected_speed)
print(momentum_residuals)
plt.plot(data['P'][-3:], momentum_residuals[-3:])

plt.show()
