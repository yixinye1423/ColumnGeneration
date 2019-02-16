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

def perturb(stageFile, sepDataFile,sysDataFile, candilogFile):
	with open(stageFile, 'rb') as fp:
		stageData = pickle.load(fp)[None]
	with open(dataFile, 'rb') as fp:
	    mstDat = pickle.load(fp)#add h_bar, change fInv, remove K, H, KH
	zkh = dict()
	for k in range(len(stageData)):
		for h in range(len(stageData[k])):
			if stageData[k][h]['selected'] == True:
				zkh[k] = h
	fInv_LO2 = list()
	fInv_LN2 = list()
	candilog = list()
	N = len(mstDat['N'][None])
	for k in range(len(stageData)):
		for h in range(len(stageData[k])): 
			firstPiece_LO2 = [sepDataFile['fInv_LO2'][(n,k,h)] for n in range(N)]
			firstPiece_LN2 = [sepDataFile['fInv_LN2'][(n,k,h)] for n in range(N)]
			candilog.append({k:h})
			for l in range(len(stageData)):
				if l == k:
					continue
				else:
					candilog[-1][l]=zkh[l]
					comb = list()
					for j in range(len(stageData)):
						if j == l:
							continue
						elif j == k:
							comb.append(stageData[j][h])
						else:
							comb.append(stageData[j][zkh[j]])
				Phi_LO2, Theta_LO2, Phi_LN2, Theta_LN2 = combine(comb)
				for n in range(N):
					firstPiece_LO2[n] += stageData[l][h]['singlepn']['LO2'][n]*Phi_LO2[n] + stageData[l][h]['stagePhi']['LO2'][n]*Theta_LO2[n]
					firstPiece_LN2[n] += stageData[l][h]['singlepn']['LN2'][n]*Phi_LN2[n] + stageData[l][h]['stagePhi']['LN2'][n]*Theta_LN2[n]
			fInv_LO2.append(firstPiece_LO2)
			fInv_LN2.append(firstPiece_LN2)
	finv_LO2 = dict()
	finv_LN2 = dict()
	for h_bar in range(len(candilog)):
		for n in range(N):
			finv_LO2[(n,h_bar)] = fInv_LO2[h_bar][n]
			finv_LN2[(n,h_bar)] = fInv_LN2[h_bar][n]
	mstDat['H_bar'] = {None:list(len(candilog))}
	mstDat['finv_LO2'] = finv_LO2
	mstDat['finv_LN2'] = fnv_LN2
	mstDat.pop('K', None)
	mstDat.pop('H', None)
	mstDat.pop('KH', None)
	with open(sysDataFile, 'wb') as fp:
		pickle.dump(mstDat, fp, protocol=pickle.HIGHEST_PROTOCOL)
	with open(candilogFile, 'wb') as fp:
		pickle.dump({None:candilog}, fp, protocol=pickle.HIGHEST_PROTOCOL)






