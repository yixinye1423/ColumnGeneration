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
	pi_hat = multiply(pis)
	pi_hat /= pi_hat.sum()
	diag_hat = add(diags)#full length vector
	dh_dense = diag_hat.todense().tolist()[0]
	#failureInd = sparse.csr_matrix([1-numpy.prod(1-numpy.asarray(tup)) for tup in itertools.product(*fails[::-1])])#all failures
	#failureInd = sparse.csr_matrix([1 if sum(tup)==1 else 0 for tup in itertools.product(*fails[::-1]) ])#no multiple failures
	avaiInd = sparse.csr_matrix([1 if sum(tup)==0 else 0 for tup in itertools.product(*fails[::-1]) ])
	N = len(V_LO2)
	f_LO2 = [None]*N
	g_LO2 = [None]*N
	f_LN2 = [None]*N
	g_LN2 = [None]*N
	#starts = time.time()
	for n in range(N):
		#g:pi*diag_exp; f:pi*diag_exp*diag
		diag_exp_LO2 = sparse.csr_matrix([numpy.exp(-V_LO2[n]/dec_LO2*dh_dense[s]) for s in range(len(dh_dense))])
		diag_exp_LN2 = sparse.csr_matrix([numpy.exp(-V_LN2[n]/dec_LN2*dh_dense[s]) for s in range(len(dh_dense))])
		f_LO2[n] = avaiInd.multiply(pi_hat.multiply(diag_hat.multiply(diag_exp_LO2))).sum()
		g_LO2[n] = avaiInd.multiply(pi_hat.multiply(diag_exp_LO2)).sum()
		f_LN2[n] = avaiInd.multiply(pi_hat.multiply(diag_hat.multiply(diag_exp_LN2))).sum()
		g_LN2[n] = avaiInd.multiply(pi_hat.multiply(diag_exp_LN2)).sum()
	print(f_LO2, g_LO2, f_LN2, g_LN2)
	#ends = time.time()
	#print(ends - starts)
	return f_LO2, g_LO2, f_LN2, g_LN2

def updateData():
	def isConverged():
		for comb in selectCombs:
			k = comb[0]
			h = comb[1]
			hr = comb[2]
			complem = [l for l in list(range(len(stageData))) if l != k]
			if stageComplem[k][hr] != tuple([selectDesign[l] for l in complem]):
				return False
		return True
	def addNewPoints():#add new subdesigns that are part of the current solution
		for k in range(len(stageData)):
			stageComplem[k] = {k:v for k,v in stageComplem[k].items() if v!=None}
			complem = [l for l in list(range(len(stageData))) if l != k]
			newSubdesign = tuple([selectDesign[l] for l in complem])
			if newSubdesign in [v for k,v in stageComplem[k].items()]:
				continue
			else:
				currSize = len(stageComplem[k])
				stageComplem[k][currSize] = newSubdesign#add to record
				comb = [stageData[l][selectDesign[l]] for l in complem]
				f_LO2, g_LO2, f_LN2, g_LN2 = combine(comb)
				for h in range(len(stageData[k])):
					pi = stageData[k][h]['pi']
					diag = stageData[k][h]['diag']
					isFailure = stageData[k][h]['isFailure']
					diag_exp_LO2 = [numpy.asarray([numpy.exp(-V_LO2[n]/dec_LO2*diag[s]) for s in range(len(diag))]) for n in range(len(V_LO2))]
					diag_exp_LN2 = [numpy.asarray([numpy.exp(-V_LN2[n]/dec_LN2*diag[s]) for s in range(len(diag))]) for n in range(len(V_LN2))]
					newPN_LO2 = [3650*pn_LO2*(sum(numpy.multiply(numpy.asarray(isFailure), numpy.multiply(numpy.asarray(pi), numpy.multiply(numpy.asarray(diag), diag_exp_LO2[n]))))*g_LO2[n]
						+sum(numpy.multiply(numpy.asarray(isFailure), numpy.multiply(numpy.asarray(pi), diag_exp_LO2[n])))*f_LO2[n]) for n in range(len(V_LO2))]
					newPN_LN2 = [3650*pn_LN2*(sum(numpy.multiply(numpy.asarray(isFailure), numpy.multiply(numpy.asarray(pi), numpy.multiply(numpy.asarray(diag), diag_exp_LN2[n]))))*g_LN2[n]
						+sum(numpy.multiply(numpy.asarray(isFailure), numpy.multiply(numpy.asarray(pi), diag_exp_LN2[n])))*f_LN2[n]) for n in range(len(V_LN2))]
					stageData[k][h]['singlepn']['LO2'][currSize] = newPN_LO2
					stageData[k][h]['singlepn']['LN2'][currSize] = newPN_LN2
	def changeModelData():
		khr = list()
		for k in range(len(stageComplem)):
			for hr in range(len(stageComplem[k])):
				khr.append((k, hr))
		mstDat['KHR'] = {None:khr}
		mstDat['HR'] = {None:list(range(max([len(stageComplem[k]) for k in range(len(stageComplem[k]))])))}
		prods = itertools.product(mstDat['KH'][None], mstDat['KHR'][None])
		KHHRlist = list()
		for prod in prods:
			if prod[0][0]==prod[1][0]:
				KHHRlist.append((prod[0][0],prod[0][1], prod[1][1]))
		mstDat['KHHR'] = {None:list(set(KHHRlist))}
		fInv_LO2 = dict()
		fInv_LN2 = dict()
		for k in range(len(stageData)):
			for h in range(len(stageData[k])):
				for hr in range(len(stageComplem[k])):
					for n in range(len(V_LO2)):
						fInv_LO2[(n,k,h,hr)] = stageData[k][h]['singlepn']['LO2'][hr][n]
						fInv_LN2[(n,k,h,hr)] = stageData[k][h]['singlepn']['LN2'][hr][n]
		mstDat['finv_LO2'] = fInv_LO2
		mstDat['finv_LN2'] = fInv_LN2


	#get basic data
	with open('stageData.p', 'rb') as fp:
		stageData = pickle.load(fp)['None']
	with open('stageComplem.p', 'rb') as fp:
		stageComplem = pickle.load(fp)['None']	
	with open('data_3433_full_sep.p', 'rb') as fp:
	    mstDat = pickle.load(fp)
	with open('result_3433_full_sep.p', 'rb') as fp:
		zkhhr = pickle.load(fp)
	print(zkhhr)

	selectCombs = {k for k,v in zkhhr.items() if abs(v-1)<0.000001}#remove unselected comb's
	print(selectCombs)
	selectDesign = {comb[0]:comb[1] for comb in selectCombs}
	print(selectDesign)
	if isConverged():
		print('Oh yeah')
		return
	else:
		addNewPoints()
		print(stageComplem)
		changeModelData()

	with open('stageData.p', 'wb') as fp:
		pickle.dump({'None':stageData}, fp, protocol=pickle.HIGHEST_PROTOCOL)

	with open('stageComplem.p', 'wb') as fp:
		pickle.dump({'None':stageComplem}, fp, protocol=pickle.HIGHEST_PROTOCOL)

	with open('data_3433_full_sep.p', 'wb') as fp:
	    pickle.dump(mstDat, fp, protocol=pickle.HIGHEST_PROTOCOL)
updateData()