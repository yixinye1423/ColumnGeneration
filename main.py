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

def powerSet(n):
	sets = list()
	for m in range(1,n+1):
		sets += itertools.combinations(range(n),m)
	return sets
def listCost(cap):
	Cap = dict()
	for k in range(len(cap)):
		designs = powerSet(len(cap[k]))
		if k==1:
			designs = [Set for Set in designs if len(Set)>=2]
		for h in range(len(designs)):
			cost = sum(cap[k][j] for j in designs[h])
			Cap[(k,h)]=cost
	return Cap
def numberAndRavel(matrix_k, states_k):
	for h in range(len(states_k)):
		for state in states_k[h]:
			state['label'] = h
	states_k = list(itertools.chain(*states_k))
	return (matrix_k, states_k)
def pseudoMat(k):
	matrix_k = list()
	states_k = list()
	for level in range(1,unitNum[k]+1):
		if (k==1 and level==1): continue
		matrix_k_level, states_k_level = generatePseudoTMatrix(k, level, parameters[k])#list of dictionaries
		matrix_k.extend(matrix_k_level)
		states_k.extend(states_k_level)
	#matrix_k, states_k = numberAndRavel(matrix_k, states_k)
	return matrix_k, states_k 
def getData(k):
	matrix_k, states_k = pseudoMat(k)#pseudo matrix of stage l
	data = list()
	for h in range(len(matrix_k)):
		mat = numpy.asarray(matrix_k[h])
		diag = mat.sum(axis = 1)
		numpy.fill_diagonal(mat, -diag)
		Qt = numpy.concatenate((numpy.transpose(mat), numpy.ones((1,mat.shape[0]))))
		b = numpy.concatenate((numpy.zeros((mat.shape[0],)), numpy.ones((1,))))
		pi = numpy.linalg.lstsq(Qt, b)
		isFailure = [state['isFailure'] for state in states_k[h]]
		data.append({'pi': pi[0], 'diag': diag, 'isFailure':isFailure})
	return data

def combine(dist):#calculate fInv_LO2, fInv_LN2
	def multiply(vectors):
		hat = vectors[0]
		for k in range(1, len(vectors)):
			hat = hat.multiply(vectors[k])
		return hat
	def add(vectors):
		hat = vectors[0]
		for k in range(1, len(vectors)):
			hat = hat + vectors[k]
		return hat
	def extend(vectors, func):
		leng = [vector.shape[0] for vector in vectors]
		pi_krons = [None]*len(vectors)
		for k in range(len(vectors)):
			A = sparse.csr_matrix(numpy.ones(int(numpy.prod(leng[k+1:]))))
			B = sparse.csr_matrix(vectors[k])
			C = sparse.csr_matrix(numpy.ones(int(numpy.prod(leng[:k]))))
			pi_krons[k] = sparse.kron(sparse.kron(A,B),C)
		return func(pi_krons)

	pis = numpy.asarray([stage['pi'] for stage in dist])
	diags = numpy.asarray([stage['diag'] for stage in dist])
	fails = [stage['isFailure'] for stage in dist]
	failureInd = sparse.csr_matrix([numpy.prod(tup) for tup in itertools.product(*fails)])
	pi_hat = extend(pis, multiply)
	pi_hat /= pi_hat.sum()
	diag_hat = extend(diags, add)#full length vector

	fInv_LO2 = [None]*5
	fInv_LN2 = [None]*5
	N = len(V_LO2)	
	for n in range(N):
		dh_dense = diag_hat.todense().tolist()[0]
		diag_exp_LO2 = [exp(-V_LO2[n]/dec_LO2*dh_dense[s]) for s in range(len(dh_dense))]
		diag_exp_LN2 = [exp(-V_LO2[n]/dec_LN2*dh_dense[s]) for s in range(len(dh_dense))]
		fInv_LO2[n] = failureInd.multiply(pi_hat.multiply(diag_hat.multiply(diag_exp_LO2))).sum()
		fInv_LN2[n] = failureInd.multiply(pi_hat.multiply(diag_hat.multiply(diag_exp_LN2))).sum()
	return (fInv_LO2, fInv_LN2)		

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
'''
Subs = list()
for k in range(len(unitNum)):
	Subs.append(PP(AbstractModel()))
'''
mstDat = dict()
mstDat['K'] = {None:list(range(len(unitNum)))}
hlist = [len(powerSet(j)) for j in unitNum]
print(hlist)
mstDat['H'] = {None:list(range(max(hlist)))}
mstDat['c_hat'] = listCost(cap)

mstDat['N'] = {None:list(range(len(V_LO2)))}
mstDat['c_LO2'] = {n: c_LO2[n] for n in range(len(c_LO2))}
mstDat['c_LN2'] = {n: c_LN2[n] for n in range(len(c_LN2))}
stageData = list()
for k in range(len(unitNum)):
	data = getData(k)
	stageData.append(data)

designs = [None]*len(unitNum)
for k in range(len(unitNum)):
	designs[k] = powerSet(unitNum[k])
	if k==1:
		designs[k] = [Set for Set in designs[k] if len(Set)>=2]
hs = [list(range(len(designs[k]))) for k in range(len(designs))]
list2d = [[(k,h) for h in hs[k]] for k in range(len(hs))]
mstDat['KH'] = {None:[val for sublist in list2d for val in sublist]}
hbars = list(itertools.product(*hs))
mstDat['H_bar'] = {None:list(range(len(hbars)))}
list2d = [[(k,hbars[i][k],i) for k in range(len(hbars[i]))] for i in range(len(hbars))]
mstDat['D'] = {None:[val for sublist in list2d for val in sublist]}

fInv = list()
for comb in itertools.product(*stageData):
	fInv.append(combine(comb))
fInv_LO2 = dict()
fInv_LN2 = dict()
for h_bar in range(len(fInv)):
	for n in range(len(fInv[h_bar][0])):
		fInv_LO2[(n,h_bar)] = fInv[h_bar][0][n]
		fInv_LN2[(n,h_bar)] = fInv[h_bar][1][n]
mstDat['finv_LO2'] = fInv_LO2
mstDat['finv_LN2'] = fInv_LN2
instance = Master.create_instance({None:mstDat})
instance.dual = Suffix(direction=Suffix.IMPORT)
opt = SolverFactory('cplex', executable="C:/Program Files/IBM/ILOG/CPLEX_Studio128/cplex/bin/x64_win64/cplex")
#opt = SolverFactory('cplex', executable="/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex")

opt.solve(instance).write()	

for kh in instance.KH:
	print(kh, instance.z[kh].value)
print('-----------------------------------------------------')
print(max(instance.z_bar[h_bar].value for h_bar in instance.H_bar))
print('-----------------------------------------------------')
for n in instance.N:
	print(n, instance.x_LO2[n].value)
for n in instance.N:
	print(n, instance.x_LN2[n].value)

#instance.dual.display()