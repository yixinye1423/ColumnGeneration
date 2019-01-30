import cplex
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
import matplotlib.pyplot as plt

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
'''
parameters = [{'lambdas':[[0.00018265,0.00027397,0.00010959,0.00054795,0.00018265,0.00010969]
,[0.00019265,0.00028397,0.00011959,0.00055795,0.00019265,0.00011969],
[0.00020265,0.00029397,0.00012959,0.00056795,0.00020265,0.00012969]],
'mus':[[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222]]},#MAC 
{'lambdas':[[0.00018265],[0.00019265],[0.00020265]],'mus':[[0.2],[0.2],[0.2]]},#PPF
{'lambdas':[[0.00018265,0.00027397,0.00010959,0.00054795,0.00018265,0.00010969]
,[0.00019265,0.00028397,0.00011959,0.00055795,0.00019265,0.0001196],
[0.00020265,0.00029397,0.00012959,0.00056795,0.00020265,0.00012969]],
'mus':[[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222]]},#BAC
{'lambdas':[[0.00054795],[0.00055795],[0.00056795]],'mus':[[2.4],[2.4],[2.4]]}]#LO2 PUMP
unitNum = [3,3,3,3]
cap = [[1250,1200,1150],[520,500,480],[1000,950,900],[150,145,140]]
'''
V_LO2 = [100, 400, 700, 1000, 1500]
V_LN2 = [100, 400, 700, 1000, 1500]
c_LO2 = [55, 237, 427, 621, 951]
c_LN2 = [50, 215, 388, 565, 864]

pn_LO2 = 2000
pn_LN2 = 2000
dec_LO2 = 48
dec_LN2 = 60
def combine(dist):#calculate fInv_LO2, fInv_LN2
	def multiply(vectors):
		spsVecs = copy.deepcopy(vectors)
		for k in range(len(vectors)):
			spsVecs[k] = sparse.csr_matrix(vectors[k])
		for k in range(1, len(vectors)):
			spsVecs[k] = sparse.kron(spsVecs[k], spsVecs[k-1])
		return spsVecs[-1]
	def add(vectors):
		leng = [vector.shape[0] for vector in vectors]
		fullVec = [None]*len(vectors)
		for k in range(len(vectors)):
			A = sparse.csr_matrix(numpy.ones(int(numpy.prod(leng[k+1:]))))
			B = sparse.csr_matrix(vectors[k])
			C = sparse.csr_matrix(numpy.ones(int(numpy.prod(leng[:k]))))
			fullVec[k] = sparse.kron(sparse.kron(A,B),C)
		hat = fullVec[0]
		for k in range(1, len(fullVec)):
			hat = hat + fullVec[k]
		return hat	
	
	pis = numpy.asarray([stage['pi'] for stage in dist])
	diags = numpy.asarray([stage['diag'] for stage in dist])
	fails = [stage['isFailure'] for stage in dist]

	start = time.time()
	pi_hat = multiply(pis)
	end = time.time()
	print('prod', end-start)

	pi_hat /= pi_hat.sum()

	start = time.time()
	diag_hat = add(diags)#full length vector
	end = time.time()
	print('add', end-start)

	#failureInd = sparse.csr_matrix([1-numpy.prod(1-numpy.asarray(tup)) for tup in itertools.product(*fails[::-1])])#all failures
	failureInd = sparse.csr_matrix([1 if sum(tup)==1 else 0 for tup in itertools.product(*fails[::-1]) ])#no multiple failures

	N = len(V_LO2)
	fInv_LO2 = [None]*N
	fInv_LN2 = [None]*N
	#starts = time.time()
	for n in range(N):
		dh_dense = diag_hat.todense().tolist()[0]
		diag_exp_LO2 = sparse.csr_matrix([exp(-V_LO2[n]/dec_LO2*dh_dense[s]) for s in range(len(dh_dense))])
		diag_exp_LN2 = sparse.csr_matrix([exp(-V_LN2[n]/dec_LN2*dh_dense[s]) for s in range(len(dh_dense))])
		fInv_LO2[n] = failureInd.multiply(pi_hat.multiply(diag_hat.multiply(diag_exp_LO2))).sum()
		fInv_LN2[n] = failureInd.multiply(pi_hat.multiply(diag_hat.multiply(diag_exp_LN2))).sum()
	#ends = time.time()
	#print(ends - starts)
	return (fInv_LO2, fInv_LN2)

def updateData():
	def isConverged():
		for comb in selectCombs:
			k = comb[0]
			h = comb[1]
			hr = comb[2]
			if stageData[k][h]['complementary'][hr] != tuple([selectDesign[l] for l in list(range(len(stageData))).remove(k)]):
				return False
		return True
	def addNewPoints():#add new subdesigns that are part of the current solution
		for k in range(len(stageData)):
			comb = [stageData[k][selectDesign[l]] for l in list(range(len(stageData))).remove(k)]
			f,g = combine(comb)


	#get basic data
	with open('stageData.p', 'rb') as fp:
		stageData = pickle.load(fp)['None']
	with open('result_2323_full_sep.p', 'rb') as fp:
		zkhhr = pickle.load(fp)
	print(zkhhr)

	selectCombs = {k for k,v in zkhhr.items() if v==1}#remove unselected comb's
	print(selectCombs)
	selectDesign = {comb[0]:comb[1] for comb in selectCombs}

	if isConverged():
		print('Oh yeah')
		return
	else:
		addNewPoints()

