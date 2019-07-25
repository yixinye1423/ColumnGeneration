import cplex
import pyomo
#from pyomo.opt import SolverFactory
#from pyomo.core import Var
from pyomo.environ import *
from construct import *
import itertools
import numpy
import math
import os
import random
from scipy import sparse
from scipy.sparse import diags
from research_supportFull import *
import time
import json
import pickle

with open('data_forpaper_3433_5.p', 'rb') as fp:
    mstDat = pickle.load(fp)

#create model objects
Master = MP(AbstractModel())
instance = Master.create_instance({None:mstDat})
instance.dual = Suffix(direction=Suffix.IMPORT)
opt = SolverFactory('cplex', executable="C:/Program Files/IBM/ILOG/CPLEX_Studio128/cplex/bin/x64_win64/cplex")
#opt = SolverFactory('cplex', executable="/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex")

opt.solve(instance,tee = True).write()	

for kh in instance.KH:
	print(kh, instance.z[kh].value)
'''
print('-----------------------------------------------------')
print(max(instance.z_bar[h_bar].value for h_bar in instance.H_bar))
'''
print('-----------------------------------------------------')
for n in instance.N:
	print(n, instance.x_LO2[n].value)
for n in instance.N:
	print(n, instance.x_LN2[n].value)
print('-----------------------------------------------------')
print(sum(instance.rfinv_LO2[n,h_bar].value for (n,h_bar) in instance.N*instance.H_bar))
print(sum(instance.rfinv_LN2[n,h_bar].value for (n,h_bar) in instance.N*instance.H_bar))
for n in instance.N:
	for h_bar in instance.H_bar:
		print((n,h_bar),instance.rfinv_LO2[n,h_bar].value)

for n in instance.N:
	for h_bar in instance.H_bar:
		print((n,h_bar),instance.rfinv_LN2[n,h_bar].value)
#instance.dual.display()
#instance.pprint()

'''
for index in instance.logic4: 
	if instance.dual[instance.logic4[index]] != 0:
		print(index, instance.dual[instance.logic4[index]])

for c in instance.component_objects(Constraint, active=True):
	duals = set()
	for index in c:
		duals.add(instance.dual[c[index]])
	print ("   Constraint",c,duals)
'''