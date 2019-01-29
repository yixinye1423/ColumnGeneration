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
from SUBfunc import *
from MPfuncs import *

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
V_LO2 = [100, 400, 700, 1000, 1500]
V_LN2 = [100, 400, 700, 1000, 1500]
c_LO2 = [55, 237, 427, 621, 951]
c_LN2 = [50, 215, 388, 565, 864]
pn_LO2 = 1000
pn_LN2 = 1000
dec_LO2 = 48
dec_LN2 = 60

#create model objects
Master = MP(AbstractModel())

Subs = list()
for k in range(len(unitNum)):
	Subs.append(PP(AbstractModel()))
#opt = SolverFactory('cplex', executable="C:/Program Files/IBM/ILOG/CPLEX_Studio128/cplex/bin/x64_win64/cplex")
opt = SolverFactory('cplex', executable="/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex")

def main():	
	def output():
		for kh in MPinstance.KH:
			print(kh, MPinstance.z[kh].value)
		print('-----------------------------------------------------')
		print(max(MPinstance.z_bar[h_bar].value for h_bar in MPinstance.H_bar))
		print('-----------------------------------------------------')
		for n in MPinstance.N:
			print(n, MPinstance.x_LO2[n].value)
		for n in MPinstance.N:
			print(n, MPinstance.x_LN2[n].value)	

	#Instantiate sub problems and master problems
	for k in range(len(unitNum)):
		subinstance_k = Subs[k].create_instance(subdataInit(k))
	MPinstance = Master.create_instance(MPdataInit())
	#instance.dual = Suffix(direction=Suffix.IMPORT)
	opt.solve(MPinstance).write(tee=True)

	#loop for adding new extreme points
	for count in range(50):
		stageDataNew = list()
		#loop for each stage
		for k in range(len(unitNum)):
			subinstance_k = Subs[k].create_instance(subdataUpdate(MPinstance, k))#update subproblem data from master solution
			opt.solve(subinstance_k).write(tee=True)#solve subproblem
			stageDataNew.append(extractSltn(subinstance_k))
		MPdataUpdate(stageDataNew)#make changes to the old stageData, design01, and MPinstance data
		opt.solve(MPinstance).write(tee=True)
		output()
#4 functions to write:
#subdataInit --working
#subdataUpdate
#extractSltn
#MPdataUpdate

main()
#instance.dual.display()