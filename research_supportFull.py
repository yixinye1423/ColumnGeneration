from __future__ import print_function
import itertools
import numpy
import copy
options = ['active','standby','being repaired']
unitNum = [2,3,2,3]
failureModes = [['rotor','bearing','Gearbox','LubeOil','motorBearing','motor'],#MAC
['generalFailure'],#PPF
['rotor','bearing','Gearbox','LubeOil','motorBearing','motor'],#BAC
['generalFailure']]#LO2 PUMP
#general basic tool (public)
def printSquareMatrix(mat):
    output = str()      
    for i in range(len(mat)):
        output += '\t'+str(i+1)
    #output += '\n'+'-'*len(mat)*8
    output += '\n'

    for i in range(len(mat)):
        output += str(i+1)+'\t'
        for ele in mat[i]:
            output += str(ele)+'\t'
        output += '\n'
    return output

def almostIs(a,b):
    if abs(a-b) < 0.000001:
        return True

def calEyeDim(k,tMatrixDimList):
    dimG = 1
    dimL = 1
    for l in range(k+1,len(tMatrixDimList)):
        dimG *= tMatrixDimList[l]
    for l in range(k):
        dimL *= tMatrixDimList[l]
    return (dimG, dimL)

def putogether(matrixList):
    length = 0
    for mat in matrixList:
        length += len(mat)
    diagMatrix = [[0]*length for i in range(length)]
    marker = 0
    for mat in matrixList:
        for row in range(len(mat)):
            for col in range(len(mat)):
                diagMatrix[marker+row][marker+col] = mat[row][col]
        marker += len(mat)
    return diagMatrix

def listMatrix(mat):
    output = str()
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            if mat[i][j] == 0:
                continue
            else:
                output += 'WM("%s","%s") = %.5f; \n' %(str(i+1), str(j+1), mat[i][j])
    return output

#matrix value assign (public)
def transit(k, stateFrom, stateTo, parameters):#locate failure mode
    def isFailure():#The transition is caused by a failure
        if not (stateFrom[:origAct] == stateTo[:origAct] 
                and stateFrom[newAct+1:] == stateTo[newAct+1:]
                and stateFrom[newAct]=='standby'
                and stateTo[origAct] in failureModes[k]):
            return False
        for i in range(origAct+1,newAct):
            if not stateFrom[i] == stateTo[i]:
                return False
        return True

    def standbyIsRepaired():#The transition is caused by a failed unit being repaired and turn to a standby
        fromChange = []
        toChange = []
        repairedIndex = None
        for i in range(len(stateFrom)): 
            if stateFrom[i] != stateTo[i]:
                fromChange.append(stateFrom[i])
                toChange.append(stateTo[i])
                repairedIndex = i
        if not (len(fromChange)==1 and fromChange[0] in failureModes[k] and toChange == ['standby']):
            return None
        else:
            return repairedIndex


    def activeIsRepaired():#The transition is caused by a failed unit being repaired and turn to active, 
        if not (stateFrom[newAct] in failureModes[k] and stateTo[origAct] == 'standby'):
            return False
        stateFrom_other = stateFrom[:newAct]+stateFrom[newAct+1:origAct]+stateFrom[origAct+1:]
        stateTo_other = stateTo[:newAct]+stateTo[newAct+1:origAct]+stateTo[origAct+1:]
        if not (stateFrom_other == stateTo_other):
            return False
        return True
    #---------------------------------------------------------------------       

    if ('active' not in stateFrom) and ('active' not in stateTo):#Both are all-down states, not communicating
        return None
    if 'active' not in stateTo:#change to all-down state:last active went down (have to check if it turns to repair or infeasible transition)
        if 'standby' in stateFrom:
            return None
        else:
            origAct = stateFrom.index('active')
            if ((stateTo[0:origAct] == stateFrom[0:origAct]) 
                and (stateTo[origAct+1:] == stateFrom[origAct+1:])):
                if stateTo[origAct] in failureModes[k]:
                    failInd = failureModes[k].index(stateTo[origAct])
                    return parameters['lambdas'][origAct][failInd]
            else:
                return None
    if 'active' not in stateFrom:#change from all-down state: have to check if it's repaired or infeasible transition
        if 'standby' in stateTo:
            return None
        else:
            newAct = stateTo.index('active')
            if ((stateFrom[0:newAct] == stateTo[0:newAct]) 
                and (stateTo[newAct+1:] == stateFrom[newAct+1:])):
                if stateFrom[newAct] in failureModes[k]:
                    failInd = failureModes[k].index(stateFrom[newAct])
                    return parameters['mus'][newAct][failInd]
            else:
                return None
    elif ('active' in stateFrom and 'active' in stateTo):#other normal states
        origAct = stateFrom.index('active');newAct = stateTo.index('active')
        if origAct < newAct:#might be breakdown happens
            if isFailure():
                failInd = failureModes[k].index(stateTo[origAct])
                value = parameters['lambdas'][origAct][failInd]
                #print(stateFrom, stateTo, value)
                return value
            else:
                return None
        elif origAct == newAct:#might be a unit repaired to become standby
            repairedIndex = standbyIsRepaired()
            if repairedIndex != None:
                failInd = failureModes[k].index(stateFrom[repairedIndex])
                value = parameters['mus'][repairedIndex][failInd]
                #print(stateFrom, stateTo, value)
                return value
            else:
                return None
        elif origAct > newAct:#might be a unit repaired to become active
            if activeIsRepaired():
                failInd = failureModes[k].index(stateFrom[newAct])
                value = parameters['mus'][newAct][failInd]
                #print(stateFrom, stateTo, value)
                return value
            else:
                return None
def transit2(k, stateFrom, stateTo, parameters):
    def isFailure():#The transition is caused by a failure
        if not (stateFrom[:origAct] == stateTo[:origAct] 
                and stateFrom[newAct+1:] == stateTo[newAct+1:]
                and stateFrom[newAct]=='standby'
                and stateTo[origAct] in failureModes[k]):
            return False
        for i in range(origAct+1,newAct):
            if not stateFrom[i] == stateTo[i]:
                return False
        return True

    def standbyIsRepaired():#The transition is caused by a failed unit being repaired and turn to a standby
        fromChange = []
        toChange = []
        repairedIndex = None
        for i in range(len(stateFrom)): 
            if stateFrom[i] != stateTo[i]:
                fromChange.append(stateFrom[i])
                toChange.append(stateTo[i])
                repairedIndex = i
        if not (fromChange in failureModes[k] and toChange == ['standby']):
            return None
        else:
            return repairedIndex


    def activeIsRepaired():#The transition is caused by a failed unit being repaired and turn to active
        if not (stateFrom[newAct] in failureModes[k] and stateTo[origAct] == 'standby'):
            return False
        for i in range(origAct+1,len(stateFrom)):
            if not (stateFrom[i]==stateTo[i]):
                return False
        return True
    #---------------------------------------------------------------------       

    if (stateFrom.count('active')<2) and (stateTo.count('active')<2):#Both are fail states, not communicating
        return None
    if stateTo.count('active')<2:#jump to fail state
        if 'standby' in stateFrom:
            return None
        else:
            compare = [stateTo[i] != stateFrom[i] for i in range(len(stateTo))]
            if sum(compare) == 1:
                failUnit = compare.index(True)
                failInd = failureModes[k].index(stateTo[failUnit])
                return parameters['lambdas'][failUnit][failInd]
            else:
                return None
    if stateFrom.count('active')<2:#jump from fail state
        if 'standby' in stateTo:
            return None
        else:
            compare = [stateTo[i] != stateFrom[i] for i in range(len(stateTo))]
            if sum(compare) == 1:
                repairedUnit = compare.index(True)
                failInd = failureModes[k].index(stateFrom[repairedUnit])
                return parameters['mus'][repairedUnit][failInd]
            else:
                return None
    elif (stateFrom.count('active')>=2) and (stateTo.count('active')>=2):#other normal states
        if stateFrom.count('standby') > stateTo.count('standby'):#failure
            origAct = [index for index, value in enumerate(stateFrom) if value == 'active']
            newAct = [index for index, value in enumerate(stateTo) if value == 'active']
            compareRepair = [stateTo[i]=='generalFailure' for i in origAct]
            compareStandby = [stateFrom[i]=='standby' for i in newAct]
            print(compareRepair)
            if (sum(compareRepair)==1 and sum(compareStandby)==1):
                return parameters['lambdas'][origAct[compareRepair.index(True)]][0]
        else:
            origAct = [index for index, value in enumerate(stateFrom) if value == 'active']
            newAct = [index for index, value in enumerate(stateTo) if value == 'active']
            compareStandby = [stateTo[i]=='standby' for i in origAct]
            compareRepair = [stateFrom[i]=='generalFailure' for i in newAct]
            if (sum(compareRepair)==1 and sum(compareStandby)==1):
                return parameters['mus'][newAct[compareRepair.index(True)]][0]

def isValidState(state):#decide whether a state is valid
    if (('standby' in state) and ('active' in state) and (state.index('standby')<state.index('active'))):
        return False
    if state.count('active') > 1:
        return False
    if ('active' not in state and 'standby' in state):
        return False
    #if (len(state)>=3 and state[-1] == 'being repaired'):#remove low probability scenarios in single failure mode stages
    #    return False
    return True

def isValidState2(state):
    if len(state) < 2:#at least 2 units
        return False
    if state.count('active') > 2:#at most 2 active units
        return False
    if state.count('being repaired') > len(state)-1:#at most 2 active units
        return False
    if (('standby' in state) and (state.count('active')<=1)):#if there's standby, at least 2 units are good
        return False
    if (('standby' in state) and (state.count('active')>=2)):
        standbyPos = state.index('standby')
        actPos = [index for index, value in enumerate(state) if value == 'active']
        for ele in actPos:
            if standbyPos < ele:
                return False
    #if (state.count('being repaired') >= len(state)):#scenario reduction
    #    return False    
    #if (state.count('being repaired') >= 2 and state.index('active') != 0):
    #    return False    

    return True
def generateTMatrix(k, design, parameter):#for selected design of stage k
    states = list()
    statesToEliminate = [state for state in 
    list(itertools.product(failureModes[0], repeat = 2))+[('active','motorBearing'),('active','Gearbox'),
    ('active','LubeOil'),('active','rotor'),('active','bearing')] if state not in [('motor','motor')]]
    if k == 1:
        for state in itertools.product(options, repeat = sum(design)):
            if isValidState2(state):
                stateSuite = dict()
                if 'being repaired' not in state:
                    stateSuite['name'] = state
                    stateSuite['isFailure'] = 0
                    #stateSuite['design'] = design
                    states.append(stateSuite)
                else:#replace 'being repaired' with failure modes
                    failNum = state.count('being repaired')
                    replacements = itertools.product(failureModes[k], repeat = failNum)
                    for replacement in replacements:
                        replacement = list(replacement)
                        newState = list()
                        for s in state:
                            if s == 'being repaired':
                                newState.append(replacement.pop())
                            else:
                                newState.append(s)
                        newState = tuple(newState)
                        if newState not in statesToEliminate:    
                            stateSuite['name'] = newState
                            stateSuite['isFailure'] = 1 if failNum==sum(design)-1 else 0
                            #stateSuite['design'] = design
                            states.append(stateSuite)
    else:
        for state in itertools.product(options, repeat = sum(design)):
            if isValidState(state):
                if 'being repaired' not in state:
                    stateSuite = dict()
                    stateSuite['name'] = state
                    stateSuite['isFailure'] = 0
                    #stateSuite['design'] = design
                    states.append(stateSuite)
                else:#replace 'being repaired' with failure modes
                    failNum = state.count('being repaired')
                    replacements = itertools.product(failureModes[k], repeat = failNum)
                    for replacement in replacements:
                        replacement = list(replacement)
                        newState = list()
                        for s in state:
                            if s == 'being repaired':
                                newState.append(replacement.pop())
                            else:
                                newState.append(s)
                        newState = tuple(newState)
                        if newState not in statesToEliminate:  
                            stateSuite = dict()  
                            stateSuite['name'] = newState
                            stateSuite['isFailure'] = 1 if failNum==sum(design) else 0
                            #stateSuite['design'] = design
                            states.append(stateSuite)

    selectParameters = dict()
    for rate in parameter:#parameter like {'lambdas':[[],[],..],'mus':[[],[],..]}
        selectParameters[rate] = list()
        for j in range(len(design)):
            if design[j] == 1:
                selectParameters[rate].append(parameter[rate][j])
        selectParameters[rate] = tuple(selectParameters[rate])

    matrix = [[0]*len(states) for i in range(len(states))]
    if k == 1:
        for i in range(len(states)):
            stateFrom = states[i]['name']
            for j in range(len(states)):
                stateTo = states[j]['name']
                value = transit2(k,stateFrom, stateTo, selectParameters)
                matrix[i][j] = 0 if value is None else value
    else:
        for i in range(len(states)):
            stateFrom = states[i]['name']
            for j in range(len(states)):
                stateTo = states[j]['name']
                value = transit(k,stateFrom, stateTo, selectParameters)
                matrix[i][j] = 0 if value is None else value
    return (matrix,states)

def generatePseudoTMatrix(k, level, parameter):#for stage k, total number of units = level
    numOfUnit = unitNum[k]#total number of potential units
    alters = itertools.combinations(range(numOfUnit), level)#all combinations of this level
    matrixList = list()
    stateList = list()
    for alter in alters:
        design = [0] * numOfUnit
        for element in alter:
            design[element] = 1
        print(design)
        matrix, states = generateTMatrix(k,design,parameter)
        for stage in states:
            print(stage)
        print('---------------------------------------------------')
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                print(matrix[i][j],end='  ')
            print(';\n')
        matrixList.append(matrix)
        stateList.append(states)
    return matrixList, stateList

def calTmatrices(parameters):#include all stage
    def calTmatrix(k, parameter):#for only stage k
        largeDiagMatrix = list()
        for i in range(unitNum[k]):
            largeDiagMatrix.append(generatePseudoTMatrix(k, i+1 ,parameter))
        tmatrix = putogether(largeDiagMatrix)
        return tmatrix

    tMatrixList = list()
    tMatrixDimList = list()
    for k in range(len(parameters)):
        tMatrixList.append(calTmatrix(k, parameters[k]))
        tMatrixDimList.append(len(tMatrixList[k]))
    return (tMatrixList,tMatrixDimList)


def readDesigns(fileDesign):
    designs = [[]*1 for i in range(len(unitNum))]
    f = open(fileDesign,'r')
    designString = f.read()[:-1]
    f.close()
    designList = designString.split('\n')
    for i in range(len(designList)):
        designList[i] = designList[i].replace('"', '').split(',')
        stageNum = int(designList[i][0])
        y = int(float(designList[i][2]))
        designs[stageNum-1].append(y)
    return designs


###########################################################################################
  
###########################################################################################


def savePastResult(filePrime, filePrimeOld, fileCounting, var):
    f1 = open(filePrime,'r')
    toSave = f1.read()
    f1.close()
    toSaveList = toSave.split(var+'(')
    fc = open(fileCounting,'r')
    num = fc.read()
    fc.close()
    newHead = var+'("'+num+'",'
    toSave = newHead.join(toSaveList)
    f2 = open(filePrimeOld,'a')
    f2.write(toSave)
    f2.close()


def convertResultToParameter(file,filePrime,var):
    f1 = open(file,'r')
    String = f1.read()[:-1]#remove the last '\n'
    f1.close()
    List = String.split('\n')
    output = str()
    for i in range(len(List)):
        List[i] = List[i].replace('"', '').split(',')#each item becomes a list looking like k,j,value
        index1 = int(List[i][0])
        index2 = int(List[i][1])
        value = int(float(List[i][2]))
        output += var + '("%s","%s") = %s; \n' %(index1, index2, value)
    f = open(filePrime, 'w')
    f.write(output)
    f.close()



'''
def calQMatrix(tMatrixList,tMatrixDimList):
    preQMatrix = list()
    for k in range(len(tMatrixList)):
        (dimG,dimL) = calEyeDim(k, tMatrixDimList)
        preQMatrix.append(numpy.kron(numpy.identity(dimG),numpy.kron(tMatrixList[k], numpy.identity(dimL))))#matrix is too large
    for k in range(1,len(preQMatrix)):
        for row in range(len(preQMatrix[0])):
            for col in range(len(preQMatrix[0][0])):
                preQMatrix[0][row][col] += preQMatrix[k][row][col]
    qmatrix = copy.deepcopy(preQMatrix[0])
    for row in range(len(qmatrix)):
        for col in range(len(qmatrix[row])):
            qmatrix[row][col] = round(qmatrix[row][col],5)
        qmatrix[row][row] = -sum(qmatrix[row])
    return qmatrix
'''




'''
def locateTMatrix(design, options, parameter):#put the Tmatrix to its designated place in the pseudoTmatrix 
    originalMatrix, originalStates = generateTMatrix(design, options, parameter)
    numOfUnit = len(parameter['lambdas'])
    alters = list()    
    for level in range(1,numOfUnit+1):
        alters += list(itertools.combinations(range(numOfUnit), level))
    designs = list()
    statesBefore = list()
    for alter in alters:
        currDesign = [0] * numOfUnit
        for element in alter:
            currDesign[element] = 1
        if currDesign != design:#currDesign must have fewer or same number of units as design 
            for state in itertools.product(options, repeat = sum(currDesign)):
                if isValidState(state):
                    statesBefore.append(state)
        else:#once the iteration reaches the design, stop
            break
    states = [None]*len(statesBefore)+originalStates
    print('original')
    print(originalMatrix)
    for i in range(len(originalMatrix)):
        originalMatrix[i] = [0]*len(statesBefore) + originalMatrix[i]
    matrix = [[0]*(len(statesBefore)+len(originalStates)) for i in range(len(statesBefore))] + originalMatrix
    return (matrix, states)    
'''


'''def generatePseudoTMatrix(level, options, parameters):#include all combinations in one stage with same number of units
    #generate state products for the design of certain number of units
    print(parameters)
    states = list()
    for state in itertools.product(options, repeat = level):
        if isValidState(state):
            states.append(state)
    #apply the state products for all combinations of this number of units
    numOfUnit = len(parameters['lambdas'])#total number of potential units
    alters = itertools.combinations(range(numOfUnit), level)#all combinations of this level (number of units in the design)
    matrixList = list()#prepare a matrix array
    for alter in alters:
        print(alter)
        #generate parameters for different combinations
        selectParameters = dict()
        for index in parameters:
            selectParameters[index] = []
            for i in alter:
                selectParameters[index].append(parameters[index][i])
            selectParameters[index] = tuple(selectParameters[index])

        matrix = [[0]*len(states) for i in range(len(states))]
        for i in range(len(states)):
            stateFrom = states[i]
            for j in range(len(states)):
                stateTo = states[j]
                value = transit(stateFrom, stateTo, selectParameters)
                matrix[i][j] = 0 if value is None else value
        matrixList.append(matrix)
    #print(matrixList)
    return putogether(matrixList)
'''
