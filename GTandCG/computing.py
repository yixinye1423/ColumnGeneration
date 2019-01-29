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

def run(dataFile, model):
	with open(dataFile, 'rb') as fp:
	    mstDat = pickle.load(fp)
	print(mstDat)
	#create model objects
	Master = model(AbstractModel())
	instance = Master.create_instance({None:mstDat})
	instance.dual = Suffix(direction=Suffix.IMPORT)
	opt.solve(instance,tee = True).write()	
	instance.pprint()
	return instance

def saveResult(resultFile, instance):
	solution = dict()
	for khhr in instance.KHHR:
		solution[khhr]= instance.z[khhr].value
		print(solution[khhr])
	print(solution)
	with open(resultFile, 'wb') as fp:
		pickle.dump(solution, fp, protocol=pickle.HIGHEST_PROTOCOL)

def procedure():
	instance = run('data_2323_full_sep.p',MP)
	saveResult('result_2323_full_sep.p',instance)
	#run('data_2323_full.p',MP)
	return

opt = SolverFactory('cplex', executable="C:/Program Files/IBM/ILOG/CPLEX_Studio128/cplex/bin/x64_win64/cplex")
#opt = SolverFactory('cplex', executable="/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex")
procedure()


#for kh in instance.KH:
#	print(kh, instance.z[kh].value)
'''
print('-----------------------------------------------------')
print(max(instance.z_bar[h_bar].value for h_bar in instance.H_bar))
'''
'''
print('-----------------------------------------------------')
for n in instance.N:
	print(n, instance.x_LO2[n].value)
for n in instance.N:
	print(n, instance.x_LN2[n].value)
print('-----------------------------------------------------')
print(sum(instance.rfinv_LO2[n,k,h].value for (n,k,h) in instance.N*instance.KH))
print(sum(instance.rfinv_LN2[n,k,h].value for (n,k,h) in instance.N*instance.KH))
for n in instance.N:
	for (k,h) in instance.KH:
		print((n,k,h),instance.rfinv_LO2[n,k,h].value)
for n in instance.N:
	for (k,h) in instance.KH:
		print((n,k,h),instance.rfinv_LN2[n,k,h].value)
'''
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