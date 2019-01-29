import os
path = os.getcwd() +'/'

def newDesign(stageNum):
	fileDesignCand = path+stageNum + 'desCand.txt'
	f = open(fileDesignCand)
	designs = f.read().split('\n')[:-1]
	f.close()
	resultFile = path+'subhh_'+stageNum+'.txt'
	f=open(resultFile,'r')
	zbr=f.read().split('\n')[:-1]
	f.close()
	for hh in range(len(zbr)):
		if int(float(zbr[hh]))==1:
			return designs[hh]

def cummulateEP(currHH, stageNum):#collect storage needs and costs of selected designs
	def convertResultToParameter(file,filePrime,var):
	    f1 = open(file,'r')
	    String = f1.read()#remove the last '\n'
	    f1.close()
	    List = String.split('\n')[:-1]
	    output = str()
	    if len(List)==1:
	        item = List[0]
	        value = float(item)
	        output += var + '("%s") = %s; \n' %(currHH+1, value)    	
	    else:
	        for i in range(len(List)):
	            List[i] = List[i].replace('"', '').split(',')#each item becomes a list looking like k,j,value
	            index = int(List[i][0])
	            value = float(List[i][1])
	            output += var + '("%s","%s") = %.10f; \n' %(index, currHH+1, value)
	    f = open(filePrime, 'a')
	    f.write(output)
	    f.close()
	convertResultToParameter(LO2+stageNum+'.txt', LO2_cum, 'freq_down_LO2')
	convertResultToParameter(LN2+stageNum+'.txt', LN2_cum, 'freq_down_LN2')
	convertResultToParameter(repa+stageNum+'.txt', repa_cum, 'repa')
	convertResultToParameter(cap+stageNum+'.txt', cap_cum, 'cap')
	f = open(fpool,'a')
	f.write(newDesign(stageNum)+'\n')#add new designs to list
	f.close()
	return True

def findObj(stageNum):
	f = open(obj+stageNum+'.txt','r')
	objective = float(f.read())
	f.close()
	return objective

def addEP():
	f = open(repa_cum,'r')
	repair = f.read()
	if repair!= '':
		currHH = len(repair.split('\n'))-1
	else:
		currHH = 0
	f.close()
	objs = list()
	for k in range(4):
		objs.append(findObj(str(k+1)))
	print(objs)
	if min(objs)>0:
		print("Time to branch")
		return
	#findNextPricingBasis(objs)
	for k in range(4):
		if objs[k]<0:
			cummulateEP(currHH, str(k+1))
			currHH += 1

LO2 = path + 'freq_down_LO2_'
LN2 = path + 'freq_down_LN2_'
repa = path + 'repa_'
cap = path + 'cap_'
obj = path + 'obj_'
LO2_cum = path + 'freq_down_LO2.txt'
LN2_cum = path + 'freq_down_LN2.txt'
repa_cum = path + 'repa_cum.txt'
cap_cum = path + 'cap_cum.txt'
fpool = path + 'pool.txt'
addEP()
