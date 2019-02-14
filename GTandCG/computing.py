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
repa = [[5, 2.4, 3, 3.5, 2, 20],#MAC 
[5],#PPF
[5, 2.4, 3, 3.5, 2, 20],#BAC
[5]]#LO2 PUMP

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
cap = [[1250,1200,1150],[520,500,480,460],[1000,950,900],[150,145,140]]

parameters = [{'lambdas':[[0.00018265,0.00027397,0.00010959,0.00054795,0.00018265,0.00010969],
[0.00018365,0.00027497,0.00011059,0.00054895,0.00018365,0.00011069],
[0.00018465,0.00027597,0.00011159,0.00054995,0.00018465,0.00011169]],
'mus':[[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222]]},#MAC 
{'lambdas':[[0.00018265],[0.00018365],[0.00018465],[0.00018565]],'mus':[[0.2],[0.2],[0.2],[0.2]]},#PPF
{'lambdas':[[0.00018265,0.00027397,0.00010959,0.00054795,0.00018265,0.00010969],
[0.00018365,0.00027497,0.00011059,0.00054895,0.00018365,0.00011069],
[0.00018465,0.00027597,0.00011159,0.00054995,0.00018465,0.00011169]],
'mus':[[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222]]},#BAC
{'lambdas':[[0.00054795],[0.00054895],[0.00054995]],'mus':[[2.4],[2.4],[2.4]]}]#LO2 PUMP
unitNum = [3,4,3,3]
cap = [[1250,1200,1150],[520,500,480,460],[1000,950,900],[150,145,140]]

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
	print(mstDat)
	#create model object
	instance = model(AbstractModel()).create_instance({None:mstDat})
	instance.dual = Suffix(direction=Suffix.IMPORT)
	opt.solve(instance,tee = True).write()	
	instance.pprint()
	return instance

def saveResult_stagewise(resultFile, instance):
	solution = dict()
	for (k,h) in instance.KH:
		solution[(k,h)]= instance.z[k,h].value
		print(solution[k,h])
	print(solution)
	with open(resultFile, 'wb') as fp: 
		pickle.dump(solution, fp, protocol=pickle.HIGHEST_PROTOCOL)

def saveResult_system(resultFile, instance):
	solution = dict()
	for h_bar in instance.H_bar:
		solution[h_bar]= instance.z_bar[h_bar].value
		print(solution[h_bar])
	print(solution)
	with open(resultFile, 'wb') as fp: 
		pickle.dump(solution, fp, protocol=pickle.HIGHEST_PROTOCOL)

def procedure():
	Init('stageData.p','data_3433_sep_0.p')
	instance = run('data_3433_sep_0.p',MP_Indep)
	saveResult_stagewise('result_3433_sep_0.p',instance)
	count = 1
	while True:
		pair('stageData.p','result_3433_sys_%s.p' %(count-1),'data_3433_sep_%s.p' %count)#calculate complementary parameters
		ECinstance = run('data_3433_sep_%s.p' %count,MP_Indep)#equilibrium checking optimization
		saveResult_stagewise('result_3433_sep_%s.p' %count, ECinstance)
		if isConverged('result_3433_sep%s.p'):#check if the result indicates equilibrium
			break
		perturb('stageData.p','result_3433_sys_%s.p' %(count-1),'data_3433_sys_%s.p' %count)#generate system level data
		POinstance = run('data_3433_sys_%s.p' %count, MP_extend)#partial optimization
		saveResult_system('result_3433_sys_%s.p' %count,POinstance)
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