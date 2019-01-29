import itertools
import numpy
import math
import os
import random
from scipy import sparse
from scipy.sparse import diags
from research_supportFull import *

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
	matrix_k, states_k = numberAndRavel(matrix_k, states_k)
	return matrix_k, states_k   

def listStateSubsets(states_k):
	hq = list()
	currh = -1
	for q in range(len(states_k)):
		h = states_k[q].get('label',0)
		if h == currh:
			hq[currh].append(q)
		else:
			hq.append([q])
			currh += 1
	return hq

def listQQHH(khq, origInd, fileQQHH):
	def belongTogether(design,state,khq):
		K = len(khq)
		for k in list(range(len(design))):#k 0,1,2,3->stage 4,3,2,1
			if state[k] not in khq[K-k-1][design[k]]:
				return False
		return True
	listOfhForProduct = [list(range(len(khq[k]))) for k in range(len(khq))]
	listOfqForProduct = [list(itertools.chain(*khq[k])) for k in range(len(khq))]
	sysDesigns = list(itertools.product(*listOfhForProduct[::-1]))
	sysStates = list(itertools.product(*listOfqForProduct[::-1]))
	#sysDesigns and sysStates are lists of tuples. Each tuple contain a design or state number in stage 4,3,2,1
	qqhh = list()
	for hh in range(len(sysDesigns)):
		for qq in range(len(origInd)):
			design = sysDesigns[hh]
			state = sysStates[origInd[qq]]
			if belongTogether(design,state,khq):
				qqhh.append((qq,hh))
	print('calculated')
	output = str()
	for (qq,hh) in qqhh:
		output += 'qqhh("%s","%s")= yes; \n' %(hh+1,qq+1)
	f1 = open(fileQQHH,'w')
	f1.write(output)
	f1.close()

def listqqf(states, origInd, fileqqf):
	kq = list()
	for k in range(len(states)):
		kq.append(list(range(len(states[k]))))
	sysStates = list(itertools.product(*kq[::-1]))
	output = str()
	K = len(kq)
	for qq in range(len(origInd)):
		sysState = sysStates[origInd[qq]]
		for k in range(len(sysState)):#k 0,1,2,3->stage 4,3,2,1
			q_k = sysState[k]
			realk = K-k-1
			if states[realk][q_k]['isFailure'] == 1:
				output += 'qqf("%s") = yes; \n' %(qq+1)
				break
	f = open(fileqqf,'w')
	f.write(output)
	f.close()

def sparseKron(dimG,dimL,tMatrix):
	A = sparse.csr_matrix(numpy.identity(dimG))
	B = sparse.csr_matrix(tMatrix)
	C = sparse.csr_matrix(numpy.identity(dimL))
	return sparse.kron(sparse.kron(A,B),C)

def sparseCalQMatrix(tMatrixList,tMatrixDimList):
	dimQQ = numpy.prod(tMatrixDimList)
	preQMatrix = list()
	indexTuple = list()
	for k in range(len(tMatrixList)):
		(dimG,dimL) = calEyeDim(k, tMatrixDimList)
		preQMatrix.append(sparseKron(dimG,dimL,tMatrixList[k]))
	qmatrix = ((sum(preQMatrix)*10000000000).rint()/10000000000).tolil()
	return qmatrix

def output(SM):
	output = str()
	cooSM = SM.tocoo()
	for i,j,v in zip(cooSM.row, cooSM.col, cooSM.data):
		output += 'WM("%d","%d")= %.10f; \n' %(i+1,j+1,v)
	return output

def reduceScen(qmatrix, states):#remove multiple failures from different stages
	kq = list()
	for k in range(len(states)):
		kq.append(list(range(len(states[k]))))
	sysStates = list(itertools.product(*kq[::-1]))
	K = len(kq)
	toStay = list()
	for qq in range(len(sysStates)):
		sysState = sysStates[qq]
		failure = 0
		for k in range(len(sysState)):#k 0,1,2,3->stage 4,3,2,1
			q_k = sysState[k]
			realk = K-k-1
			if states[realk][q_k]['isFailure'] == 1:
				failure += 1
		if failure <= 1:#remove multiple failures from different stages
			toStay.append(qq)
	mask = sparse.csr_matrix(sparse.eye(len(sysStates)))
	mask = mask[toStay]
	Lmat = mask.copy()
	Rmat = mask.transpose()
	print(Lmat.shape)
	reducedQ = (Lmat.dot(qmatrix)).dot(Rmat)
	reducedQ = reducedQ.tolil()
	sumOverCol = sparse.csr_matrix(-reducedQ.sum(axis = 1))
	reducedQ.setdiag(sumOverCol.data,k=0)
	return reducedQ, toStay


	
def listCost(l,skip, design, fileCost):
	output = str()
	fixedCost = 0
	for k in skip:
		for j in range(len(design[k])):
			fixedCost += cap[k][j]*design[k][j]
	hh = 1
	cap_l = cap[l]
	for level in range(1,unitNum[l]+1):
		if (l==1 and level==1):
			continue
		alters = itertools.combinations(range(1,unitNum[l]+1), level)
		for alter in alters:
			varyCost = 0
			for j in alter:
				varyCost += cap_l[j-1]
			output += 'inst("%s") = %d; \n' %(hh, fixedCost+varyCost)
			hh += 1
	f = open(fileCost,'w')
	f.write(output)
	f.close()
'''
def listRepa(l,skip,states, origInd, fileRepa): 
	kOrder = [l]+skip
	kq = list()#correspondance between stage k and state q
	for k in range(len(states)):
		kq.append(list(range(len(states[k]))))
	sysStates = list(itertools.product(*kq[::-1]))
	output = str()
	K = len(kq)
	for qq in range(len(origInd)):
		repaCost = 0
		sysState = sysStates[origInd[qq]]
		for k in range(len(sysState)):#k 0,1,2,3->stage 4,3,2,1
			q_k = sysState[k]
			realk = K-k-1
			for unit in range(len(states[realk][q_k]['name'])):
				unitState = states[realk][q_k]['name'][unit]
				if unitState in failureModes[kOrder[realk]]:
					f = failureModes[kOrder[realk]].index(unitState)
					repaCost += repa[kOrder[realk]][f]
		output += 'c_repa("%s")=%d; \n' %(qq+1, repaCost)
	f = open(fileRepa,'w')
	f.write(output)
	f.close()
'''
def listDesignCandidates(l, fileDesignCand):
	numOfUnit = unitNum[l]#total number of potential units
	output = str()
	for level in range(1,numOfUnit+1):
		if (l==1 and level==1): continue
		alters = itertools.combinations(range(numOfUnit), level)#all combinations of this level
		for alter in alters:
			subdesign = [0] * numOfUnit
			for element in alter:
				subdesign[element] = 1
			designm = copy.deepcopy(design)
			designm[l] = subdesign
			output += str(designm)+'\n'
	f = open(fileDesignCand,'w')
	f.write(output)
	f.close()
#get individual stage design
def recoverDesign(base):
	base = [stage.split(', ') for stage in base[2:-2].split('], [')]
	print(base)
	for k in range(len(base)):
		for j in range(len(base[k])):
			base[k][j] = int(base[k][j])
	print(base)
	return base

def findPricingBasis():
	designFile = path+'pool.txt'
	f=open(designFile, 'r')
	designs=f.read().split('\n')[:-1]
	f.close()
	resultFile = path+'hh.txt'
	f=open(resultFile,'r')
	zbr=f.read().split('\n')[:-1]
	f.close()
	print(designs)
	print(zbr)
	for hh in range(len(zbr)):
		if float(zbr[hh])>0.5:
			return recoverDesign(designs[hh])

def preparePricing(l, fileQQHH, fileqqf, fileQMatrix, fileCost, fileRepa, fileDesignCand):
	matrices = list()
	states = list()
	matrix_l, states_l = pseudoMat(l)#pseudo matrix of stage l
	matrices.append(putogether(matrix_l))
	states.append(states_l)
	#exact matrices of the other stages
	skip = list(range(len(unitNum)))
	skip.remove(l)
	for k in skip:
		print('fixed',design[k])
		matrix_k, states_k = generateTMatrix(k, design[k], parameters[k])#list of dictionaries
		print(matrix_k)
		matrices.append(matrix_k)
		states.append(states_k)
	matrixDims = [len(matrices[k]) for k in range(len(matrices))]
	qmatrix = sparseCalQMatrix(matrices,matrixDims)
	qmatrix, origInd = reduceScen(qmatrix, states)
	matrix = output(qmatrix)
	file = open(fileQMatrix,'w')
	file.write(matrix)
	file.close()
	khq = list()#correspondance between stage k, design h and state q
	for k in range(len(states)):
		khq.append(listStateSubsets(states[k]))
	print(khq)
	listQQHH(khq, origInd, fileQQHH)
	listqqf(states, origInd, fileqqf)
	listCost(l,skip, design, fileCost)
	#listRepa(l,skip,states, origInd, fileRepa)
	listDesignCandidates(l, fileDesignCand)

path = os.getcwd() +'/'
fpool = path+'pool.txt'
fhh = path+'hh.txt'
design = findPricingBasis()
print(design)
for l in range(len(unitNum)):
	pathl = path + str(l+1)
	fileQQHH = pathl + 'qqhh.txt'
	fileqqf = pathl + 'qqf.txt'
	fileQMatrix = pathl + 'qMatrix.txt'
	fileCost = pathl + 'inst.txt'
	fileRepa = pathl + 'repa.txt'
	fileDesignCand = pathl + 'desCand.txt'
	preparePricing(l, fileQQHH, fileqqf, fileQMatrix, fileCost, fileRepa,fileDesignCand)


vLO2 = path+'vLO2.txt'
vLN2 = path+'vLN2.txt'
def translateMargin(file, var):
	f = open(file,'r')
	numberList = f.read().split('\n')[:-1]
	f.close()
	output = str()
	for n in range(len(numberList)):
		numberList[n] = var+'("'+str(n+1)+'")='+numberList[n]+';\n'
	output = ''.join(numberList)
	f = open(file,'w')
	f.write(output)
	f.close()
translateMargin(vLO2, 'vLO2')
translateMargin(vLN2, 'vLN2')