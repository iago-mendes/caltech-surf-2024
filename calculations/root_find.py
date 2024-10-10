from scipy import optimize

def F(x):
  Fx = [x[0]**2, x[1]**2]
  print(f'F({x}) = {Fx}')
  return Fx

result = optimize.root(F, [0.1, 0.1], method='broyden1',
                       options={
                          'fatol': 1.e-4,
                          'maxiter': 15,
                          'jac_options': {'alpha': .25}
                       })

print(result)
