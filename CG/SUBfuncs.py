import cplex
import pyomo
#from pyomo.opt import SolverFactory
#from pyomo.core import Var
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

def listCost_k(cap, design, k):
	Cap = dict()
	for h in range(len(design)):
		cost = sum(cap[k][j] for j in design[h])
		Cap[h] = cost
	return Cap

def numberAndRavel(states_k):
	for h in range(len(states_k)):
		for state in states_k[h]:
			state['label'] = h
	states_k_ravel = list(itertools.chain(*states_k))
	return states_k_ravel

def pseudoMat(k, design_k, options, failureModes, parameters):
	matrix_k = list()
	states_k = list()
	for design_k_h in design_k:
		matrix_k_h, states_k_h = generateTMatrix(k, design_k_h, parameters[k])#list of dictionaries
		matrix_k.append(matrix_k_h)
		states_k.append(states_k_h)
	#matrix_k, states_k = numberAndRavel(matrix_k, states_k)
	return matrix_k, states_k 

def getData(k, design_k, options, failureModes, parameters):
	#print(design_k)
	matrix_k, states_k = pseudoMat(k, design_k, options, failureModes, parameters)#pseudo matrix of stage k
	#print(matrix_k)
	data = list()
	for h in range(len(matrix_k)):
		mat = numpy.asarray(matrix_k[h])
		diag = mat.sum(axis = 1)
		numpy.fill_diagonal(mat, -diag)
		#Qt = numpy.concatenate((numpy.transpose(mat), numpy.ones((1,mat.shape[0]))))
		#b = numpy.concatenate((numpy.zeros((mat.shape[0],)), numpy.ones((1,))))
		#pi = numpy.linalg.lstsq(Qt, b)#solve the linear system in advance
		isFailure = [state['isFailure'] for state in states_k[h]]
		#data.append({'pi': pi[0], 'diag': diag, 'isFailure':isFailure})
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
	failureInd = sparse.csr_matrix([1-numpy.prod(1-numpy.asarray(tup)) for tup in itertools.product(*fails)])
	pi_hat = extend(pis, multiply)#do the kronecker product
	pi_hat /= pi_hat.sum()#normalization
	diag_hat = extend(diags, add)#do the kronecker sum

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


def subdataInit(k, unitNum, options, failureModes, parameters, cap):
'''
s, hs, sf, q
'''
	subDat = dict()
	design = powerSet(unitNum[k])
	if k==1:
		design = [Set for Set in design if len(Set)>=2]
	design01 = copy.deepcopy(design)
	print(design01)
	for h in range(len(design)):
		design01[h] = [0]*unitNum[k]
		for j in range(unitNum[k]):
			design01[h][j] = 1 if j in design[h] else 0
	subDat['H'] = {None:list(range(len(design)))}#number of designs for each stage
	subDat['inst'] = listCost_k(cap,design,k)#cap cost of each design
	matrix_k, states_k = pseudoMat(k, design01, options, failureModes, parameters)
	subDat['Q'] = {None: putogether(matrix_k)}
	states_k_ravel = numberAndRavel(states_k)
	subDat['S'] = {None: len(states_k_ravel)}
	isFailure = [states_k_ravel.index(state) for state in states_k_ravel if state['isFailure']==1]
	subDat['U'] = 
	return {None:subDat}	

def subdataUpdate():
	pass

#instance.dual.display()
'''
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
'''