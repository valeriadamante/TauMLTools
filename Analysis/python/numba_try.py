from numba import jit, njit
import numpy as np
import time
import uproot4
import awkward as ak
import awkward.numba
import math
start_0 = time.time()
fileName = 'DataSetTrainingWeight.root'
treeName = 'L2TauTrainTuple'
absolute_path='/home/Valeria/Desktop/Dottorato/L2SkimmedTuples/DataSetTraining/'
filePath = absolute_path+fileName
file = uproot4.open(filePath)


tree = file[treeName]
start_2 = time.time()
print (("time to get file {}").format(start_2-start_0))

awkArray = tree.arrays(how="zip")
start_3 = time.time()
print (("time to get awkarray {}").format(start_3-start_2))


'''
vecVars = {
    "caloRecHit_e" : ['energy', 'DeltaEta','DeltaPhi','chi2'],
    "caloRecHit_had" : ['energy', 'DeltaEta','DeltaPhi','chi2'],
    "patatrack" : ['pt','eta','phi','DeltaEta','DeltaPhi','DeltaR','chi2','ndof','charge','dxy','dz','hasVertex','vert_z','vert_ptv2','vert_chi2','vert_ndof']
}
'''
flatVars = [ "nVertices", "l1Tau_pt", "l1Tau_eta", "l1Tau_phi","l1Tau_hwIso" ]
#caloVars = [ ["caloRecHit_e",  'energy'], ["caloRecHit_e",  'DeltaEta'], ["caloRecHit_e", 'DeltaPhi'], ["caloRecHit_e", 'chi2'], ["caloRecHit_had",  'energy'], ["caloRecHit_had",  'DeltaEta'], ["caloRecHit_had", 'DeltaPhi'], ["caloRecHit_had", 'chi2']]
caloVars = [ ["caloRecHit_e",  'energy']]#, ["caloRecHit_e",  'DeltaEta'], ["caloRecHit_e", 'DeltaPhi'], ["caloRecHit_e", 'chi2'], ["caloRecHit_had",  'energy'], ["caloRecHit_had",  'DeltaEta'], ["caloRecHit_had", 'DeltaPhi'], ["caloRecHit_had", 'chi2']]
pataVars= [ ["patatrack", 'pt'], ["patatrack" , 'eta'], ["patatrack" , 'phi'], ["patatrack" , 'DeltaEta'], ["patatrack" , 'DeltaPhi'], ["patatrack" , 'DeltaR'], ["patatrack" , 'chi2'], ["patatrack" , 'ndof'], ["patatrack" , 'charge'], ["patatrack" , 'dxy'], ["patatrack" , 'dz'], ["patatrack" , 'hasVertex'], ["patatrack" , 'vert_z'], ["patatrack" , 'vert_ptv2'], ["patatrack" , 'vert_chi2'], ["patatrack" , 'vert_ndof'] ]

n_cellsX = 5
n_cellsY = 5
#nEvents = len(getattr(awkArray, flatVars[0]))
#n_vars = len(flatVars)
#CellGrid = np.zeros((nEvents, n_cellsX, n_cellsY, n_vars))

#1. fill flat vars
@jit
def getCellGridMatrix(awkArray, flatVars, caloVars, pataVars, n_cellsX, n_cellsY, verbose=False):
    print("writing feature matrix")
    nVars = 32
    #nEvents = len(awkArray[flatVars[0]])
    nTaus = 2
    CellGrid = np.zeros((nTaus, n_cellsX, n_cellsY, nVars))
    dPhi_min = -0.5
    dPhi_max = 0.5
    dEta_min = -0.5
    dEta_max = 0.5
    dPhi_width = (dPhi_max-dPhi_min)/(n_cellsX-1)
    dEta_width = (dEta_max-dEta_min)/(n_cellsX-1)
    varIndex = 0
    # 1. fill flat vars
    for var in flatVars:
        for i in range(0,nTaus):
            k=awkArray[var][i]
            for j in range(0,n_cellsX):
                CellGrid[i][j][j][varIndex] = k
        varIndex += 1
    for var in caloVars:
        a = var[0]
        b = var[1]
        for i in range(0, nTaus):
            for j in range(0,len(awkArray[a][b][i])):
                varToFill = awkArray[a][b][i][j]
                dphi = awkArray[a]['DeltaPhi'][i][j]
                deta = awkArray[a]['DeltaEta'][i][j]
                energy = awkArray[a]['energy'][i][j]
                if(b == "chi2" and varToFill <0 ):
                    continue
                phi_idx = int(math.floor((dphi + dPhi_max) / dPhi_width + 0.5))
                eta_idx = int(math.floor((deta + dEta_max) / dEta_width + 0.5))
                CellGrid[i][phi_idx][eta_idx][varIndex]  = varToFill * energy
            #weigh in energy
        varIndex += 1
    return CellGrid

CellGrid = getCellGridMatrix(awkArray,flatVars, caloVars, pataVars,  n_cellsX, n_cellsY, True)
assert not getCellGridMatrix.nopython_signatures # this was not compiled in nopython mode
start_4 = time.time()
print (("time to define grid {}").format(start_4-start_3))
'''
# 2. fill energy


@jit
def fillFlatVars(CellGrid, awkArray, var, varIndex, n_cellsX, n_cellsY, verbose=False): # Function is compiled to machine code when called the first time
    for i in range(0,n_cellsX):
        dEta_current = dEta_min + i * dEta_width
        for j in range(0, n_cellsY):
            dPhi_current = dPhi_min + j * dPhi_width
            if(verbose):
                print(dEta_current, dPhi_current)
            # 1. copy global variables : nVertex - l1Tau_pt - l1Tau_eta - l1Tau_phi - l1Tau_hwIso
            CellGrid[:, i, j, varIndex] = getattr(awkArray, var)


# DO NOT REPORT THIS... COMPILATION TIME IS INCLUDED IN THE EXECUTION TIME!
start = time.time()
varIndex = 0
for var in flatVars:
    fillFlatVars(CellGrid, awkArray, var, varIndex, n_cellsX, n_cellsY)
    varIndex +=1
end = time.time()
print("Elapsed (with compilation) = %s" % (end - start))

# NOW THE FUNCTION IS COMPILED, RE-TIME IT EXECUTING FROM CACHE
start = time.time()
varIndex = 0
for var in flatVars:
    fillFlatVars(CellGrid, awkArray, var, varIndex, n_cellsX, n_cellsY)
    varIndex +=1
end = time.time()
print("Elapsed (after compilation) = %s" % (end - start))
#print( CellGrid)
'''
