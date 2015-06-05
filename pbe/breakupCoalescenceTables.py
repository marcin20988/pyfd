import fluid
import numpy as np
import matplotlib.pyplot as plt

path = 'validationData/tables/1/'
f = open(path + 'parameters', 'w')
f_br = open(path + 'breakup', 'w')

F = fluid.fluid('simmonsAzzopardi', 0, 0)
F.C1 = 1.52e-03
F.C2 = 6.78e-02
F.C3 = 1.06e-12
F.C4 = 5.13e13

f.write('C1 ' + repr(F.C1) + '\n')
f.write('C2 ' + repr(F.C2) + '\n')
f.write('C3 ' + repr(F.C3) + '\n')
f.write('C4 ' + repr(F.C4) + '\n')
f.write('rhoc ' + repr(F.rhoc) + '\n')
f.write('rhod ' + repr(F.rhod) + '\n')
f.write('muc ' + repr(F.muc) + '\n')
f.write('epsilon ' + repr(F.epsilon) + '\n')
f.write('sigma ' + repr(F.sigma) + '\n')
f.write('alpha ' + repr(F.alpha) + '\n')
f.write('V ' + repr(F.V) + '\n')
f.close()

dMin = 0.1 * F.d0
dMax = 10.0 * F.d0
dd = 0.1 * F.d0
d = np.arange(dMin, dMax, dd)
gamma = F.gamma(d)
Q = F.Q(d, d)
result_breakup = zip(d, gamma)
np.savetxt(f_br, result_breakup)

fig = plt.figure()

plt.plot(d, Q)
plt.show()
