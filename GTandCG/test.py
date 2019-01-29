#import cplex
#import pyomo
#from pyomo.opt import SolverFactory
#from pyomo.core import Var
#from pyomo.environ import *
import math
import numpy as np
import matplotlib.pyplot as plt

a = [i/1000 for i in range(10,100)]
exp_a = [np.exp(-i) for i in a]
approx_exp_a = [1-i for i in a]
ratio = [(exp_a[i]-approx_exp_a[i])/exp_a[i] for i in range(len(a))]
#print(exp_a)
#print(approx_exp_a)
#plt.plot(range(100), a, label = 'a')
#plt.plot(a, exp_a, label = 'exp_a')
#plt.plot(a, approx_exp_a, label = 'approx_exp_a')
#plt.plot(a, ratio, label = 'ratio')
V_LO2 = [100, 400, 700, 1000, 1500]
dec_LO2 = 48
steps = [np.exp(-V_LO2[i]/dec_LO2) for i in range(len(V_LO2))]
print(steps)
plt.plot(range(len(V_LO2)), steps)
plt.legend()
plt.show()