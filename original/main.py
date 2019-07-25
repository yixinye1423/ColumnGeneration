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
		data.append({'pi': pi, 'diag': diag, 'isFailure':isFailure})
	#plt.plot(range(len(matrix_k)), justk)
			#print(k,h, n, 3650*pn_LO2*sum(numpy.multiply(numpy.asarray(isFailure), numpy.multiply(numpy.asarray(pi), numpy.multiply(numpy.asarray(diag), diag_exp)))))
	return data

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
	
	dist = dist[::-1]
	pis = numpy.asarray([stage['pi'] for stage in dist])
	diags = numpy.asarray([stage['diag'] for stage in dist])
	fails = [stage['isFailure'] for stage in dist]

	start = time.time()
	pi_hat = multiply(pis)
	end = time.time()
	#print('prod', end-start)

	#pi_hat /= pi_hat.sum()

	start = time.time()
	diag_hat = add(diags)#full length vector
	end = time.time()
	#print('add', end-start)

	#failureInd = sparse.csr_matrix([1-numpy.prod(1-numpy.asarray(tup)) for tup in itertools.product(*fails[::-1])])#all failures
	failureInd = sparse.csr_matrix([1 if sum(tup)>=1 else 0 for tup in itertools.product(*fails[::-1]) ])#no multiple failures
	#print(pi_hat.toarray())
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
'''
#parameters,V_LO2, V_LN2, dec_LO2, dec_LN2, pn_LO2, pn_LN2, c_LO2, c_LN2, cap = peel('C:/Users/yye1/Documents/GitHub/ColumnGeneration/GTandCG/rawData.p')
#First level indices
parameters = [{'lambdas':[[0.00018265,0.00027397,0.00010959,0.00054795,0.00018265,0.00010969],
[0.00018275,0.00027407,0.00010969,0.00054805,0.00018275,0.00010979]],
'mus':[[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222]]},#MAC 
{'lambdas':[[0.00018265],[0.00018275],[0.00018285]],'mus':[[0.2],[0.2],[0.2]]},#PPF
{'lambdas':[[0.00018265,0.00027397,0.00010959,0.00054795,0.00018265,0.00010969],
[0.00018275,0.00027407,0.00010969,0.00054805,0.00018275,0.00010979]],
'mus':[[0.1,0.2,0.143,3,0.5,0.0222],[0.1,0.2,0.143,3,0.5,0.0222]]},#BAC
#{'lambdas':[[0.00054795],[0.00054805],[0.00054815]],'mus':[[2.4],[2.4],[2.4]]}]#LO2 PUMP
#{'lambdas':[[0.001644],[0.001654],[0.001664]],'mus':[[1.2],[1.2],[1.2]]}]#LO2 PUMP
{'lambdas':[[0.0054795],[0.0054805],[0.0054815]],'mus':[[2.4],[2.4],[2.4]]}]#LO2 PUMP
unitNum = [2,3,2,3]
#cap = [[1250,1200,1150],[520,500,480,460],[1000,950,900],[150,145,140]]
#cap = [[1250,1249,1248],[520,519,518,517],[1000,999,998],[150,149,148]]
cap = [[1250,1250],[520,520,520],[1000,1000],[150,150,150]]

V_LO2 = [100, 400, 700, 1000, 1500]
V_LN2 = [100, 400, 700, 1000, 1500]
c_LO2 = [4*v**0.6 for v in V_LO2]
c_LN2 = [3.5*v**0.6 for v in V_LN2]
c_LO2 = [55, 237, 427, 621, 951]
c_LN2 = [50, 215, 388, 565, 864]

pn_LO2 = 2000
pn_LN2 = 2000
dec_LO2 = 48
dec_LN2 = 60
'''

def peel(rawDataFile):
	with open(rawDataFile, 'rb') as fp:
	    rawData = pickle.load(fp)
	return tuple([rawData[key] for key in rawData])

parameters,V_LO2, V_LN2, dec_LO2, dec_LN2, pn_LO2, pn_LN2, c_LO2, c_LN2, cap = peel('rawData_4544_1.p')
#print(parameters)
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
#Second level indices
hs = [list(range(len(stageData[k]))) for k in range(len(stageData))]#full

print(hs)
list2d = [[(k,h) for h in hs[k]] for k in range(len(hs))]

mstDat['KH'] = {None:[val for sublist in list2d for val in sublist]}
mstDat['c_hat'] = listCost(cap, hs)

hbars = list(itertools.product(*hs[::-1]))
#for hbar in range(len(hbars)):
#	print(hbar, hbars[hbar][::-1])
mstDat['H_bar'] = {None:list(range(len(hbars)))}
list2d = [[(k,hbars[i][len(unitNum)-1-k],i) for k in range(len(hbars[i]))] for i in range(len(hbars))]
mstDat['D'] = {None:[val for sublist in list2d for val in sublist]}

stageData = [[stageData[k][i] for i in hs[k]] for k in range(len(stageData))]

#Calculate frequencies
fInv = list()
count = 0
for comb in itertools.product(*stageData[::-1]):
	fInv.append(combine(comb))
	count += 1
	print(count)
print(len(fInv))

fInv_LO2 = dict()
fInv_LN2 = dict()
#forPlot = list()

for h_bar in range(len(fInv)):
	for n in range(len(fInv[h_bar][0])):
		fInv_LO2[(n,h_bar)] = 3650*pn_LO2*fInv[h_bar][0][n]
		fInv_LN2[(n,h_bar)] = 3650*pn_LN2*fInv[h_bar][1][n]
	#forPlot.append(fInv_LO2[(2,h_bar)])
#plt.plot(list(range(len(forPlot))), forPlot)
#plt.show()
#print(forPlot)
mstDat['finv_LO2'] = fInv_LO2
mstDat['finv_LN2'] = fInv_LN2



with open('data_forpaper_4544_1.p', 'wb') as fp:
    pickle.dump(mstDat, fp, protocol=pickle.HIGHEST_PROTOCOL)

