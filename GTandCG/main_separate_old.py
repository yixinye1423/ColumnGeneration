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

def powerSet(n):
	sets = list()
	for m in range(1,n+1):
		sets += itertools.combinations(range(n),m)
	return sets
def listCost(cap, hs):
	Cap = dict()
	for k in range(len(cap)):
		designs = powerSet(len(cap[k]))
		if k==1:
			designs = [Set for Set in designs if len(Set)>=2]
		for h in hs[k]:
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
		pi = numpy.linalg.lstsq(Qt, b)[0]
		isFailure = [state['isFailure'] for state in states_k[h]]

		diag_exp_LO2 = [numpy.asarray([numpy.exp(-V_LO2[n]/dec_LO2*diag[s]) for s in range(len(diag))]) for n in range(len(V_LO2))]
		diag_exp_LN2 = [numpy.asarray([numpy.exp(-V_LN2[n]/dec_LN2*diag[s]) for s in range(len(diag))]) for n in range(len(V_LN2))]
		approx_diag_exp_LO2 = [numpy.asarray([1-V_LO2[n]/dec_LO2*diag[s] for s in range(len(diag))]) for n in range(len(V_LO2))]

		singlepn_LO2 = [3650*pn_LO2*sum(numpy.multiply(numpy.asarray(isFailure), numpy.multiply(numpy.asarray(pi), numpy.multiply(numpy.asarray(diag), diag_exp_LO2[n])))) for n in range(len(V_LO2))]
		singlepn_LN2 = [3650*pn_LN2*sum(numpy.multiply(numpy.asarray(isFailure), numpy.multiply(numpy.asarray(pi), numpy.multiply(numpy.asarray(diag), diag_exp_LN2[n])))) for n in range(len(V_LN2))]
		print(singlepn_LO2, singlepn_LN2)
		data.append({'pi': pi, 'diag': diag, 'isFailure':isFailure, 'singlepn':{'LO2':singlepn_LO2, 'LN2': singlepn_LN2}})
	return data


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

#First level indices
mstDat = dict()
mstDat['K'] = {None:list(range(len(unitNum)))}
hlist = [len(powerSet(j)) for j in unitNum]
mstDat['H'] = {None:list(range(max(hlist)))}


mstDat['N'] = {None:list(range(len(V_LO2)))}
mstDat['c_LO2'] = {n: c_LO2[n] for n in range(len(c_LO2))}
mstDat['c_LN2'] = {n: c_LN2[n] for n in range(len(c_LN2))}

#matrix
stageData = list()
for k in range(len(unitNum)):
	data = getData(k)
	stageData.append(data)

hs = [list(range(len(stageData[k]))) for k in range(len(stageData))]#full

print(hs)
list2d = [[(k,h) for h in hs[k]] for k in range(len(hs))]

mstDat['KH'] = {None:[val for sublist in list2d for val in sublist]}
mstDat['c_hat'] = listCost(cap, hs)

stageData = [[stageData[k][h] for h in hs[k]] for k in range(len(hs))]

fInv_LO2 = dict()
fInv_LN2 = dict()
for k in range(len(hs)):
	for h in range(len(hs[k])):
		for n in range(len(V_LO2)):
			fInv_LO2[(n,k,h)] = stageData[k][h]['singlepn']['LO2'][n]
			fInv_LN2[(n,k,h)] = stageData[k][h]['singlepn']['LN2'][n]
mstDat['finv_LO2'] = fInv_LO2
mstDat['finv_LN2'] = fInv_LN2

with open('data_2323_full_sep.p', 'wb') as fp:
    pickle.dump(mstDat, fp, protocol=pickle.HIGHEST_PROTOCOL)

with open('stageData.p', 'wb') as fp:
	pickle.dump({'None':stageData}, fp, protocol=pickle.HIGHEST_PROTOCOL)

