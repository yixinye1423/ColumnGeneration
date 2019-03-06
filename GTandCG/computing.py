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

#basic data
options = ['active','standby','being repaired']
failureModes = [['rotor','bearing','Gearbox','LubeOil','motorBearing','motor'],#MAC
['generalFailure'],#PPF
['rotor','bearing','Gearbox','LubeOil','motorBearing','motor'],#BAC
['generalFailure']]#LO2 PUMP

parameters = [{'lambdas':[[0.00018265,0.00027397,0.00010959,0.00054795,0.00018265,0.00010969]
,[0.00019265,0.00028397,0.00011959,0.00055795,0.00019265,0.0001196]],
'mus':[[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222]]},#MAC 
{'lambdas':[[0.00018265],[0.00019265],[0.00020265]],'mus':[[0.2],[0.2],[0.2]]},#PPF
{'lambdas':[[0.00018265,0.00027397,0.00010959,0.00054795,0.00018265,0.00010969]
,[0.00019265,0.00028397,0.00011959,0.00055795,0.00019265,0.0001196]],
'mus':[[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222]]},#BAC
{'lambdas':[[0.00054795],[0.00055795],[0.00056795]],'mus':[[2.4],[2.4],[2.4]]}]#LO2 PUMP
unitNum = [2,3,2,3]
cap = [[1250,1200],[520,500,480],[1000,950],[150,145,140]]


parameters = [{'lambdas':[[0.00018265,0.00027397,0.00010959,0.00054795,0.00018265,0.00010969]
,[0.00019265,0.00028397,0.00011959,0.00055795,0.00019265,0.00011969],
[0.00020265,0.00029397,0.00012959,0.00056795,0.00020265,0.00012969]],
'mus':[[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222]]},#MAC 
{'lambdas':[[0.00018265],[0.00019265],[0.00020265],[0.00021265]],'mus':[[0.2],[0.2],[0.2],[0.2]]},#PPF
{'lambdas':[[0.00018265,0.00027397,0.00010959,0.00054795,0.00018265,0.00010969]
,[0.00019265,0.00028397,0.00011959,0.00055795,0.00019265,0.0001196],
[0.00020265,0.00029397,0.00012959,0.00056795,0.00020265,0.00012969]],
'mus':[[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222]]},#BAC
{'lambdas':[[0.00054795],[0.00055795],[0.00056795]],'mus':[[2.4],[2.4],[2.4]]}]#LO2 PUMP
unitNum = [3,4,3,3]
#cap = [[1250,1200,1150],[520,500,480,460],[1000,950,900],[150,145,140]]
cap = [[1250,1240,1230],[520,515,510,505],[1000,990,980],[150,145,140]]


parameters = [{'lambdas':[[0.00018265,0.00027397,0.00010959,0.00054795,0.00018265,0.00010969],
[0.00018275,0.00027407,0.00010969,0.00054805,0.00018275,0.00010979],
[0.00018285,0.00027417,0.00010979,0.00054815,0.00018285,0.00010989]],
'mus':[[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222]]},#MAC 
{'lambdas':[[0.00018265],[0.00018275],[0.00018285],[0.00018295]],'mus':[[0.2],[0.2],[0.2],[0.2]]},#PPF
{'lambdas':[[0.00018265,0.00027397,0.00010959,0.00054795,0.00018265,0.00010969],
[0.00018275,0.00027407,0.00010969,0.00054805,0.00018275,0.00010979],
[0.00018285,0.00027417,0.00010979,0.00054815,0.00018285,0.00010989]],
'mus':[[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222]]},#BAC
#{'lambdas':[[0.00054795],[0.00054805],[0.00054815]],'mus':[[2.4],[2.4],[2.4]]}]#LO2 PUMP
#{'lambdas':[[0.001644],[0.001654],[0.001664]],'mus':[[1.2],[1.2],[1.2]]}]#LO2 PUMP
{'lambdas':[[0.0054795],[0.0054805],[0.0054815]],'mus':[[2.4],[2.4],[2.4]]}]#LO2 PUMP
unitNum = [3,4,3,3]
#cap = [[1250,1200,1150],[520,500,480,460],[1000,950,900],[150,145,140]]
#cap = [[1250,1249,1248],[520,519,518,517],[1000,999,998],[150,149,148]]
cap = [[1250,1250,1250],[520,520,520,520],[1000,1000,1000],[150,150,150]]



V_LO2 = [100, 400, 700, 1000, 1500]
V_LN2 = [100, 400, 700, 1000, 1500]
c_LO2 = [55, 237, 427, 621, 951]
c_LN2 = [50, 215, 388, 565, 864]

pn_LO2 = 2000
pn_LN2 = 2000
dec_LO2 = 48
dec_LN2 = 60

def run(dataFile, model):
	with open(dataFile, 'rb') as fp:
	    mstDat = pickle.load(fp)
	#create model object
	instance = model(AbstractModel()).create_instance({None:mstDat})
	instance.dual = Suffix(direction=Suffix.IMPORT)
	opt.solve(instance,tee = False).write()	
	#instance.pprint()
	return instance

def saveResult_stagewise(stageFile, instance):
	solution = dict()
	for (k,h) in instance.KH:
		solution[(k,h)]= instance.z[k,h].value
	print(solution)
	for index in instance.invLO21:
		if instance.dual[instance.invLO21[index]]!=0:
			print(instance.dual[instance.invLO21[index]])
	with open(stageFile, 'rb') as fp:
		stageData = pickle.load(fp)[None]
	for k in range(len(stageData)):
		for h in range(len(stageData[k])):
			stageData[k][h]['selected'] = False
			if abs(solution[(k,h)] - 1) < 10**(-5):
				stageData[k][h]['selected'] = True
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
				


def procedure():
	open('stageData.p','w')
	open('data_sep.p','w')
	open('data_sys.p','w')
	open('candilog.p','w')
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