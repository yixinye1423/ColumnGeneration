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
from init import *
from perturb import *
from deviation import *



def run(dataFile, model):
	with open(dataFile, 'rb') as fp:
	    mstDat = pickle.load(fp)
	#create model object
	instance = model(AbstractModel()).create_instance({None:mstDat})
	instance.dual = Suffix(direction=Suffix.IMPORT)
	opt.solve(instance,tee = True)	
	#instance.pprint()
	return instance

def saveResult_stagewise(stageFile, instance):
	solution = dict()
	for (k,h) in instance.KH:
		solution[(k,h)]= instance.z[k,h].value
	print(solution)
	#for index in instance.invLO21:
	#	if instance.dual[instance.invLO21[index]]!=0:
	#		print(instance.dual[instance.invLO21[index]])
	with open(stageFile, 'rb') as fp:
		stageData = pickle.load(fp)[None]
	for k in range(len(stageData)):
		for h in range(len(stageData[k])):
			stageData[k][h]['selected'] = False
			if abs(solution[(k,h)] - 1) < 10**(-5):
				stageData[k][h]['selected'] = True
	print('unit cost',sum(instance.c_hat[k,h]*instance.z[k,h].value for (k,h) in instance.KH))
	print('LO2 tank cost',sum(instance.x_LO2[n].value*instance.c_LO2[n] for n in instance.N))
	print('LN2 tank cost',sum(instance.x_LN2[n].value*instance.c_LN2[n] for n in instance.N))	
	print('LO2 penalty',sum(instance.rfinv_LO2[n,k,h].value for (n,k,h) in instance.N*instance.KH))
	print('LN2 penalty',sum(instance.rfinv_LN2[n,k,h].value for (n,k,h) in instance.N*instance.KH))
	with open(stageFile, 'wb') as fp:
		pickle.dump({None:stageData}, fp, protocol=pickle.HIGHEST_PROTOCOL)

def saveResult_system(stageFile, candilogFile, instance):#to be changed
	with open(stageFile, 'rb') as fp:
		stageData = pickle.load(fp)[None]
	with open(candilogFile, 'rb') as fp:
		candilog = pickle.load(fp)[None]
	h_bar_star = 0
	for h_bar in instance.H_bar:
		if abs(instance.z_bar[h_bar].value-1)<10**(-5):
			print(h_bar)
	for n in instance.N:
		if abs(instance.x_LO2[n].value-1)<10**(-5):
			print(n)
		if abs(instance.x_LN2[n].value-1)<10**(-5):
			print(n)
	for h_bar in instance.H_bar:
		if abs(instance.z_bar[h_bar].value-1)<10**(-5):
			h_bar_star = h_bar
			break
			#solution = candilog[h_bar]
			#break
	#print(solution)
	print(candilog[h_bar_star])
	for k in range(len(stageData)):
		for h in range(len(stageData[k])):
			stageData[k][h]['selected'] = False
			if candilog[h_bar_star][k] == h:
				stageData[k][h]['selected'] = True

	print('unit cost',sum(instance.c_hat[h_bar]*instance.z_bar[h_bar].value for h_bar in instance.H_bar))
	print('LO2 tank cost',sum(instance.x_LO2[n].value*instance.c_LO2[n] for n in instance.N))
	print('LN2 tank cost',sum(instance.x_LN2[n].value*instance.c_LN2[n] for n in instance.N))	
	print('LO2 penalty',sum(instance.rfinv_LO2[n,h_bar].value for (n,h_bar) in instance.N*instance.H_bar))
	print('LN2 penalty',sum(instance.rfinv_LN2[n,h_bar].value for (n,h_bar) in instance.N*instance.H_bar))

	with open(stageFile, 'wb') as fp: 
		pickle.dump({None:stageData}, fp, protocol=pickle.HIGHEST_PROTOCOL)

def isConverged(stageFile, instance):
	with open(stageFile, 'rb') as fp:
		stageData = pickle.load(fp)[None]
	for k in range(len(stageData)):
		for h in range(len(stageData[k])):
			if abs(instance.z[k,h].value - 1) < 10**(-5):
				print(k,h) 
	print('unit cost',sum(instance.c_hat[k,h]*instance.z[k,h].value for (k,h) in instance.KH))
	print('LO2 tank cost',sum(instance.x_LO2[n].value*instance.c_LO2[n] for n in instance.N))
	print('LN2 tank cost',sum(instance.x_LN2[n].value*instance.c_LN2[n] for n in instance.N))	
	print('LO2 penalty',sum(instance.rfinv_LO2[n,k,h].value for (n,k,h) in instance.N*instance.KH))
	print('LN2 penalty',sum(instance.rfinv_LN2[n,k,h].value for (n,k,h) in instance.N*instance.KH))
	for k in range(len(stageData)):
		for h in range(len(stageData[k])):
			if abs(instance.z[k,h].value - 1) < 10**(-5):
				#print(k,h)
				if stageData[k][h]['selected'] == False:
					return False
	return True
'''				
def peel(rawDataFile):
	with open(rawDataFile, 'rb') as fp:
	    rawData = pickle.load(fp)
	    print(rawData)
	print(tuple([rawData[key] for key in rawData]))
	return tuple([rawData[key] for key in rawData])
'''
def procedure(rawDataFile):
	open('stageData.p','w')
	open('data_sep.p','w')
	open('data_sys.p','w')
	open('candilog.p','w')
	with open(rawDataFile, 'rb') as fp:
	    rawData = pickle.load(fp)
	parameters = rawData['parameters']
	V_LO2, V_LN2 = rawData['V_LO2'], rawData['V_LN2']
	dec_LO2, dec_LN2 = rawData['dec_LO2'], rawData['dec_LN2']
	pn_LO2, pn_LN2 = rawData['pn_LO2'], rawData['pn_LN2']
	c_LO2, c_LN2 = rawData['c_LO2'], rawData['c_LN2']
	cap = rawData['cap']
	Init('stageData.p','data_sep.p', parameters,V_LO2, V_LN2, dec_LO2, dec_LN2, pn_LO2, pn_LN2, c_LO2, c_LN2, cap)
	instance = run('data_sep.p',MP_Indep)
	saveResult_stagewise('stageData.p',instance)
	count = 1
	while True:
		pair('stageData.p','data_sep.p',V_LO2, V_LN2, dec_LO2, dec_LN2)#calculate complementary parameters
		ECinstance = run('data_sep.p',MP_Indep)#equilibrium checking optimization
		#saveResult_stagewise('stageData.p' %count, ECinstance)
		if isConverged('stageData.p',ECinstance):#check if the result indicates equilibrium
			print('converged')
			break
		perturb('stageData.p','data_sep.p', 'data_sys.p', 'candilog.p',V_LO2, V_LN2, dec_LO2, dec_LN2)#generate system level data
		POinstance = run('data_sys.p', MP_extend)#partial optimization
		saveResult_system('stageData.p', 'candilog.p', POinstance)
		count += 1

 

#opt = SolverFactory('cplex', executable="C:/Program Files/IBM/ILOG/CPLEX_Studio128/cplex/bin/x64_win64/cplex")
opt = SolverFactory('cplex', executable="/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex")
opt.options['mipgap']=0
procedure('rawData_3433_5.p')

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