import numpy as np
import math
import os
import sys
import time
import random
import matplotlib.pyplot as plt

#Setting parameters
def getTEParams(workDir,dMode,pMode) :
    TEP = TEParams()
    
    TEP.progMode = pMode
    TEP.workPath = workDir
    TEP.timingFileName = workDir + "/" + "timing.txt"
    TEP.execPath = "../src/DTE"
    TEP.fileTEOp = workDir+"/"+"TE.txt"
    TEP.saveDataFile = workDir+"/"+"savedData"
    TEP.numNodes = 2
    TEP.delayMode = dMode
    TEP.calcMode = 1
    TEP.Tmax = 10000
    TEP.rate = 4.1
    TEP.delayMin = 2
    TEP.delayMax = 4
    TEP.delDelay = 0.05
    TEP.binWidth = 0.05

    return TEP

#Forward delay model parameters

def getDelayModelParams(TEP) :

    #Various delay models
    #1. Constant delay
    dPConst = dParams()
    dPConst.name = 'Constant'
    dPConst.numParams=1
    dPConst.paramList.append(3.0)

    dPUnif = dParams()
    dPUnif.name = 'Uniform'
    dPUnif.numParams=2
    dPUnif.paramList.append(2.5)
    dPUnif.paramList.append(3.5)

    dPGauss = dParams()
    dPGauss.name = 'Gauss'
    dPGauss.numParams=2
    dPGauss.paramList.append(3.0)
    dPGauss.paramList.append(0.2)
    
    dPRet = dParams()
    if ((TEP.delayMode == "Constant") or (TEP.delayMode == "constant")) :
        dPRet = dPConst
    elif ((TEP.delayMode == "Uniform") or (TEP.delayMode == "uniform")):
        dPRet = dPUnif
    elif ((TEP.delayMode == "Gauss") or (TEP.delayMode == "gauss")) :
        dPRet = dPGauss

    return dPRet

#Parent Poisson process
def getParentProcess(TEPar):
    
    X=list()
    t=0

    while (t < TEPar.Tmax):
        X.append(t)
        delT = -math.log(np.random.rand())/TEPar.rate;
        delT = math.ceil(delT*10000)/10000 #4 decimal places
        t=t+delT
    return X

def getDelayedProcess(X,dP) :

    Y=list()
   
    if ((dP.name == 'Constant') or (dP.name == 'constant')):
        constDelay = dP.paramList[0]
        for t in X :
            Y.append(t+constDelay)

    if ((dP.name == 'Uniform') or (dP.name == 'uniform')):
        a = dP.paramList[0]
        b = dP.paramList[1]

        for t in X :
            Y.append(t+np.random.uniform(a,b))

    if ((dP.name == 'Gauss') or (dP.name == 'gauss')):
        mu = dP.paramList[0]
        std = dP.paramList[1]

        for t in X :
            Y.append(t+np.random.normal(mu,std))
    
    Y.sort()

    return Y

def printProcessesToFile(pX,pY,TEPar) :
    

    fName = TEPar.timingFileName
    f = open(fName, 'w')
    strX=''
    strY=''
    
    #Write the first process
    count = 0
    for t in pX :
        if (count == 0) :
            strX = str(t)
        else :
            strX = strX + ',' + str(t)
        count = count + 1
    
    #Write the second process
    count = 0
    for t in pY :
        if (count == 0) :
            strY = str(t)
        else :
            strY = strY + ',' + str(t)
        count = count + 1

    f.write(strX)
    f.write('\n')
    f.write(strY)
    f.write('\n')

    f.close()

def plotTransferEntropyProfile(TEPar) :

    if ((TEPar.progMode == "Save") or (TEPar.progMode == "save")) :
        numDelays = int(math.ceil((TEPar.delayMax-TEPar.delayMin)/TEPar.delDelay)+1)
        sweepDelay = np.linspace(TEPar.delayMin,TEPar.delayMax,numDelays)
        sweptTE = np.zeros((TEPar.numNodes,TEPar.numNodes,numDelays))
        
        for idx in range(0,numDelays) :
            execCmd = TEPar.execPath+" "+TEP.workPath+" "+TEPar.timingFileName+" "+str(TEPar.binWidth)+" "+str(sweepDelay[idx])+" "+ str(TEPar.calcMode)
            print execCmd
            #Run the TE code
            os.system(execCmd)
            #Gather the results and store
            gatherResults(TEPar,sweptTE,idx)
        
        # Adding the delDelay parameter finally while plotting to account for the next index in Transfer Entropy definition
        sweepDelay = sweepDelay+TEPar.binWidth
        invokePlotter(sweptTE,sweepDelay)
        
        #Save the data
        np.savez(TEPar.saveDataFile,sweptTE=sweptTE,sweepDelay=sweepDelay)

    elif ((TEPar.progMode == "Load") or (TEPar.progMode == "load")) :
        
        fileName = str(TEPar.saveDataFile) + str(".npz")
        if (os.path.isfile(fileName)) :
            savedFile = np.load(fileName)
            sweptTE = savedFile['sweptTE']
            sweepDelay = savedFile['sweepDelay']
            invokePlotter(sweptTE,sweepDelay)
        else :
            print "Unable to read saved data from " + fileName

def gatherResults(TEPar,sweptTE,idx) :
    
    #TE Data for this run
    f = open(TEPar.fileTEOp,'r')
    dataTE = f.readlines() # will append in the list out
    #Remove the comment (first line)
    dataTE.pop(0)
    f.close()

    #Parse lines and store results
    for line in dataTE :
        words = line.split()
        if words :
            #print "sTE("+str(int(words[0])) + "," + str(int(words[1])) + "," + str(idx) + ") = " + str(float(words[2]))
            sweptTE[int(words[0]),int(words[1]),idx] = float(words[2])

def invokePlotter(sTE,sDelay) :

    #Plotter function currently plots only the TE values for node pairs [0,1] and [1,0]

    x = sDelay
    nDelays = sDelay.size

    te01 = np.zeros((nDelays))
    te10 = np.zeros((nDelays))
    
    #print "\n From plotter : \n"
    #print "Size of sTE = " + str(sTE.shape)
    
    for i in range(0,nDelays) :
        
        #print "Index : " + str(i) + "  out of " + str(nDelays)
        te01[i] = sTE[0,1,i]
        te10[i] = sTE[1,0,i]

    #Plotting
    plt.figure(1)
    plt.plot(x,te01,label = 'DTE (0 -> 1)')
    plt.plot(x,te10,color = 'r',label = 'DTE (1 -> 0)')
    plt.legend(loc=1, ncol=3, shadow=True)
    plt.xlabel('Delay in time units :')
    plt.ylabel('Transfer Entropy in bits : ')
    plt.title('Swept Delayed Transfer Entropy')
    
    plt.show()

#Class definitions (essentially working as structs in this script)
class dParams:
    def __init__(self) :
        self.name = ''
        self.paramList = []
        self.numParams = 0

class TEParams:
    def __init__(self) :

        self.progMode = ''
        self.workPath = ''
        self.timingFileName = ''
        self.execPath = ''
        self.fileTEOp = ''
        self.numNodes = 0
        self.delayMode = ''
        self.Tmax = 0
        self.rate = 0
        self.delayMin = 0
        self.delayMax = 0
        self.delDelay = 0
        self.binWidth = 0



#Main function

totalArgs = len(sys.argv)

if ((totalArgs == 3) or (totalArgs == 4)) :

    dirName= (str)(sys.argv[1])
    progMode = (str)(sys.argv[2])
    delayMode = ""

    if (totalArgs == 4) :
        delayMode = (str)(sys.argv[3])

    workPath = "../workingDir/" + dirName

    if ((progMode == "Save") or (progMode == "save")) :
        
        if (delayMode == "") :
            print "Please enter a delay mode and rerun. Current choices : Constant or Uniform or Gauss"
        elif(os.path.exists(workPath)):   # Check for path
            print "Warning: " + dirName + " folder already exists. Please choose a new name to avoid overwriting.\n"
        else:
            os.makedirs(workPath)
        
            TEP = getTEParams(workPath,delayMode,progMode)
            dP = getDelayModelParams(TEP)
            processX = getParentProcess(TEP)
            processY = getDelayedProcess(processX,dP)
            printProcessesToFile(processX,processY,TEP)
            plotTransferEntropyProfile(TEP)

    elif ((progMode == "Load") or (progMode == "load")) :

        if(os.path.exists(workPath)):   # Check for path
            TEP = getTEParams(workPath,delayMode,progMode)
            plotTransferEntropyProfile(TEP)
        else:
            print "Warning: " + dirName + " folder does not exist. Please choose a new graph name.\n"

    else :
        print "Nothing to do."

else:
    print "Usage :"
    print "python DTEDriver.py directoryName programMode delayMode"
    print "delayMode is needed only if programMode is save and not needed when programMode is load"

