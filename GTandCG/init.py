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
#Independent design selection
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
		#approx_diag_exp_LO2 = [numpy.asarray([1-V_LO2[n]/dec_LO2*diag[s] for s in range(len(diag))]) for n in range(len(V_LO2))]

		singlepn_LO2 = [3650*pn_LO2*sum(numpy.multiply(numpy.asarray(isFailure), numpy.multiply(numpy.asarray(pi), numpy.multiply(numpy.asarray(diag), diag_exp_LO2[n])))) for n in range(len(V_LO2))]
		singlepn_LN2 = [3650*pn_LN2*sum(numpy.multiply(numpy.asarray(isFailure), numpy.multiply(numpy.asarray(pi), numpy.multiply(numpy.asarray(diag), diag_exp_LN2[n])))) for n in range(len(V_LN2))]
		stagePhi_LO2 = [3650*pn_LO2*sum(numpy.multiply(numpy.asarray(isFailure), numpy.multiply(numpy.asarray(pi), diag_exp_LO2[n]))) for n in range(len(V_LO2))]
		stagePhi_LN2 = [3650*pn_LN2*sum(numpy.multiply(numpy.asarray(isFailure), numpy.multiply(numpy.asarray(pi), diag_exp_LN2[n]))) for n in range(len(V_LN2))]

		print(singlepn_LO2, singlepn_LN2)
		data.append({'pi': pi, 'diag': diag, 'isFailure':isFailure, 'singlepn':{'LO2':singlepn_LO2, 'LN2': singlepn_LN2},'stagePhi':{'LO2':stagePhi_LO2, 'LN2':stagePhi_LN2}, 'selected':False})
		#'complementary' list of the format:{1:(1,2,4), ...}
	return data

def Init(stageFile, dataFile):
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
	stageData = [[stageData[k][h] for h in hs[k]] for k in range(len(hs))]

	print(hs)
	list2d = [[(k,h) for h in hs[k]] for k in range(len(hs))]

	mstDat['KH'] = {None:[val for sublist in list2d for val in sublist]}
	mstDat['c_hat'] = listCost(cap, hs)

	fInv_LO2 = dict()
	fInv_LN2 = dict()
	for k in range(len(hs)):
		for h in range(len(hs[k])):
			for n in range(len(V_LO2)):
				fInv_LO2[(n,k,h)] = stageData[k][h]['singlepn']['LO2'][n]
				fInv_LN2[(n,k,h)] = stageData[k][h]['singlepn']['LN2'][n]
	mstDat['finv_LO2'] = fInv_LO2
	mstDat['finv_LN2'] = fInv_LN2

	with open(stageFile, 'wb') as fp:
		pickle.dump({'None':stageData}, fp, protocol=pickle.HIGHEST_PROTOCOL)

	with open(dataFile, 'wb') as fp:
	    pickle.dump(mstDat, fp, protocol=pickle.HIGHEST_PROTOCOL)