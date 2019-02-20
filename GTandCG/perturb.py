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

def combine(dist,V_LO2, V_LN2, dec_LO2, dec_LN2):#calculate fInv_LO2, fInv_LN2
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
	#print(f_LO2, g_LO2, f_LN2, g_LN2)
	#ends = time.time()
	#print(ends - starts)
	return f_LO2, g_LO2, f_LN2, g_LN2

def perturb(stageFile, sepDataFile,sysDataFile, candilogFile, V_LO2, V_LN2, dec_LO2, dec_LN2):
	with open(stageFile, 'rb') as fp:
		stageData = pickle.load(fp)[None]
	if os.path.getsize(sysDataFile) == 0:
		with open(sepDataFile, 'rb') as fp:
		    mstDat = pickle.load(fp)
	else:
		with open(sysDataFile, 'rb') as fp:
			mstDat = pickle.load(fp)
		print(a)
	zkh = dict()
	for k in range(len(stageData)):
		for h in range(len(stageData[k])):
			if stageData[k][h]['selected'] == True:
				zkh[k] = h
	print("perturbing based on: ", zkh)
	c_hat = mstDat['c_hat']
	fInv_LO2 = list()
	fInv_LN2 = list()
	candilog = list()
	cost = list()
	N = len(mstDat['N'][None])
	with open(candilogFile, 'rb') as fp:
		oldCandilog = pickle.load(fp)[None]

	for k in range(len(stageData)):
		for h in range(len(stageData[k])): 
			candi = {k:h}
			for l in range(len(stageData)):
				if l != k:
					candi[l]=zkh[l]
			if candi in oldCandilog:
				continue	
			candilog.append(candi)
			firstPiece_LO2 = [mstDat['finv_LO2'][(n,k,h)] for n in range(N)]
			firstPiece_LN2 = [mstDat['finv_LN2'][(n,k,h)] for n in range(N)]
			cost.append(c_hat[(k,h)])
			for l in range(len(stageData)):
				if l == k:
					continue
				else:
					cost[-1]+=c_hat[(l,zkh[l])]
					comb = list()
					for j in range(len(stageData)):
						if j == l:
							continue
						elif j == k:
							comb.append(stageData[j][h])
						else:
							comb.append(stageData[j][zkh[j]])
				Phi_LO2, Theta_LO2, Phi_LN2, Theta_LN2 = combine(comb,V_LO2, V_LN2, dec_LO2, dec_LN2)
				for n in range(N):
					firstPiece_LO2[n] += (stageData[l][zkh[l]]['singlepn']['LO2'][n]*Phi_LO2[n] + stageData[l][zkh[l]]['stagePhi']['LO2'][n]*Theta_LO2[n])
					firstPiece_LN2[n] += (stageData[l][zkh[l]]['singlepn']['LN2'][n]*Phi_LN2[n] + stageData[l][zkh[l]]['stagePhi']['LN2'][n]*Theta_LN2[n])
			fInv_LO2.append(firstPiece_LO2)
			fInv_LN2.append(firstPiece_LN2)
	candilog = oldCandilog + candilog
	finv_LO2 = dict()
	finv_LN2 = dict()
	for h_bar in range(len(oldCandilog), len(candilog)):
		for n in range(N):
			finv_LO2[(n,h_bar)] = fInv_LO2[h_bar-len(oldCandilog)][n]
			finv_LN2[(n,h_bar)] = fInv_LN2[h_bar-len(oldCandilog)][n]	
	mstDat['H_bar'] = {None:list(range(len(candilog)))}
	mstDat['finv_LO2'].update(finv_LO2)
	mstDat['finv_LN2'].update(finv_LN2)
	mstDat['c_hat'].update({h_bar+len(oldCandilog): cost[h_bar] for h_bar in range(len(cost))})
	mstDat.pop('K', None)
	mstDat.pop('H', None)
	mstDat.pop('KH', None)
	#print(mstDat)
	with open(sysDataFile, 'wb') as fp:
		pickle.dump(mstDat, fp, protocol=pickle.HIGHEST_PROTOCOL)
	with open(candilogFile, 'wb') as fp:
		pickle.dump({None:candilog}, fp, protocol=pickle.HIGHEST_PROTOCOL)






