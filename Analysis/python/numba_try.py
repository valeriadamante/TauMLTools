from numba import jit, njit
from numba.typed import List
import numpy as np
import time
import uproot
import awkward as ak
#import awkward.numba
import math
import argparse
import array
parser = argparse.ArgumentParser()
parser.add_argument('--time', required=False, type=bool, default=False , help='if set to true, display time information')
parser.add_argument('--n_max_events', required=False, type=int, default=-1, help='max number of events to be processed')
parser.add_argument('--input_file', required=False, type=str, default='DataSetTrainingWeight.root', help='input file name')
parser.add_argument('--input_tuple', required=False, type=str, default='L2TauTrainTuple', help='input tree name')
parser.add_argument('--n_cellsX', required=False, type=int, default=5, help='number of cells along X dir')
parser.add_argument('--n_cellsY', required=False, type=int, default=5, help='number of cells along Y dir')
#parser.add_argument('--output_tuple', required=False, type=str, default='L2TauTrainTuple.root', help='output file name')
#parser.add_argument('-p', required=False, type=bool, default=True, help=' True = process , False = just print')
#parser.add_argument('--eff',required=False, type=bool, default=False, help=' True = evaluate efficiency/rate ')
#parser.add_argument('--eff_algo', required=False, type=bool, default=False, help=' True = evaluate algo efficiency ')
#parser.add_argument('--verbose', required=False, type=bool, default=0)
args = parser.parse_args()

start_0 = time.time()
absolute_path='/Users/valeriadamante/Desktop/Dottorato/L2SkimmedTuples/DataSetTraining/'
#absolute_path='/home/valeria/'
fileName = args.input_file
treeName = args.input_tuple
filePath = absolute_path+fileName
file = uproot.open(filePath)
tree = file[treeName]
start_1 = time.time()
if(args.time):
    print (("time to get file {}").format(start_1-start_0))
#awkArray = tree.arrays(how="zip")
if(args.n_max_events!=-1):
    awkArray =tree.arrays(entry_stop=args.n_max_events)
else:
    awkArray =tree.arrays()
start_2 = time.time()
if(args.time):
    print (("time to get awkarray {}").format(start_2-start_1))


flatVars = [ "nVertices", "l1Tau_pt", "l1Tau_eta", "l1Tau_phi","l1Tau_hwIso" ]
#caloVars = [ ["caloRecHit_e",  'energy'], ["caloRecHit_e",  'DeltaEta'], ["caloRecHit_e", 'DeltaPhi'], ["caloRecHit_e", 'chi2'], ["caloRecHit_had",  'energy'], ["caloRecHit_had",  'DeltaEta'], ["caloRecHit_had", 'DeltaPhi'], ["caloRecHit_had", 'chi2']]
caloVars = [ "caloRecHit_e_energy", "caloRecHit_e_DeltaEta", "caloRecHit_e_DeltaPhi", "caloRecHit_e_chi2", "caloRecHit_had_energy", "caloRecHit_had_DeltaEta", "caloRecHit_had_DeltaPhi", "caloRecHit_had_chi2"]
pataVars= [ "patatrack_pt", "patatrack_eta", "patatrack_phi", "patatrack_DeltaEta", "patatrack_DeltaPhi", "patatrack_DeltaR", "patatrack_chi2", "patatrack_ndof", "patatrack_charge", "patatrack_dxy", "patatrack_dz", "patatrack_hasVertex", "patatrack_vert_z", "patatrack_vert_ptv2", "patatrack_vert_chi2", "patatrack_vert_ndof" ]

n_cellsX = args.n_cellsX
n_cellsY = args.n_cellsY
nEvents = len(getattr(awkArray, flatVars[0]))
# define flat Vars --> they will be repeated
flatMatrix = np.zeros((len(flatVars), nEvents))
for i in range(0, len(flatVars)):
    flatMatrix[i, :] = getattr(awkArray, flatVars[i])

#patatrackPtList = List()
#for i in range(0, nEvents):
#    patatrackPtArray = np.array(ak.to_numpy(awkArray.patatrack_pt[i]))
#    patatrackPtList.append(patatrackPtArray)

#PtEnMatrixList.append(hcalEnergyList)
#PtEnMatrixList.append(patatrackPtList)

# define ECAL Vars
ecalVarsList = List()

ecalEnergyList = List()
for i in range(0, nEvents):
    ecalEnergyArray = awkArray.caloRecHit_e_energy[i]
    ecalEnergyList.append(ecalEnergyArray)
ecalVarsList.append(ecalEnergyList)

ecalDeltaEtaList = List()
for i in range(0, nEvents):
    ecalDeltaEtaArray = awkArray.caloRecHit_e_DeltaEta[i]
    ecalDeltaEtaList.append(ecalDeltaEtaArray)
ecalVarsList.append(ecalDeltaEtaList)

ecalDeltaPhiList = List()
for i in range(0, nEvents):
    ecalDeltaPhiArray = awkArray.caloRecHit_e_DeltaPhi[i]
    ecalDeltaPhiList.append(ecalDeltaPhiArray)
ecalVarsList.append(ecalDeltaPhiList)

ecalChi2List = List()
for i in range(0, nEvents):
    ecalChi2Array = awkArray.caloRecHit_e_chi2[i][awkArray.caloRecHit_e_chi2[i]>0]
    ecalChi2List.append(ecalChi2Array)
ecalVarsList.append(ecalChi2List)

# define HCAL Vars
hcalVarsList = List()

hcalEnergyList = List()
for i in range(0, nEvents):
    hcalEnergyArray = np.array(ak.to_numpy(awkArray.caloRecHit_had_energy[i]))
    hcalEnergyList.append(hcalEnergyArray)
hcalVarsList.append(hcalEnergyList)

hcalDeltaEtaList = List()
for i in range(0, nEvents):
    hcalDeltaEtaArray = np.array(ak.to_numpy(awkArray.caloRecHit_had_DeltaEta[i]))
    hcalDeltaEtaList.append(hcalDeltaEtaArray)
hcalVarsList.append(hcalDeltaEtaList)

hcalDeltaPhiList = List()
for i in range(0, nEvents):
    hcalDeltaPhiArray = np.array(ak.to_numpy(awkArray.caloRecHit_had_DeltaPhi[i]))
    hcalDeltaPhiList.append(hcalDeltaPhiArray)
hcalVarsList.append(hcalDeltaPhiList)

hcalChi2List = List()
for i in range(0, nEvents):
    hcalChi2Array = np.array(ak.to_numpy(awkArray.caloRecHit_had_chi2[i][awkArray.caloRecHit_had_chi2[i]>0]))
    hcalChi2List.append(hcalDeltaPhiArray)
hcalVarsList.append(hcalChi2List)


@njit()
def addCellGridSumVar(CellGrid, nTaus, varList, varPos, varIndex, dR_min, dR_max, dPhi_width, dEta_width):
    for tau in range(0,nTaus):
        #print(len(varList[varPos][tau]))
        for item in range(0, len(varList[varPos][tau])):
            deta = varList[1][tau][item]
            dphi = varList[2][tau][item]
            phi_idx = int(math.floor((dphi + dR_max) / dPhi_width + 0.5))
            eta_idx = int(math.floor((deta + dR_max) / dEta_width + 0.5))
            CellGrid[tau][varIndex][phi_idx][eta_idx]+= varList[varPos][tau][item]
    varIndex += 1
    return varIndex
@njit()
def addCellGridWeightedSumVar(CellGrid, nTaus, varList, varPos, varIndex,dR_min, dR_max, dPhi_width, dEta_width):
    for tau in range(0,nTaus):
        EnSum = 0
        #print(len(varList[varPos][tau]))
        for item in range(0, len(varList[varPos][tau])):
            energy = varList[0][tau][item]
            deta = varList[1][tau][item]
            dphi = varList[2][tau][item]
            phi_idx = int(math.floor((dphi + dR_max) / dPhi_width + 0.5))
            eta_idx = int(math.floor((deta + dR_max) / dEta_width + 0.5))
            CellGrid[tau][varIndex][phi_idx][eta_idx]+= varList[varPos][tau][item] * energy
            #print("has added")
            EnSum += energy
        CellGrid[tau][varIndex][phi_idx][eta_idx]/=EnSum
    varIndex += 1
    return varIndex
@njit()
def addCellGridDevStdVar(CellGrid, nTaus, varList, varIndex,dR_min, dR_max, dPhi_width, dEta_width):
    for tau in range(0,nTaus):
        EnMean = 0
        for item in range(0, len(varList[0][tau])):
            energy = varList[0][tau][item]
            EnMean += energy/len(varList[0][tau])
        for item in range(0, len(varList[0][tau])):
            energy = varList[0][tau][item]
            dev_std = math.sqrt(math.pow((energy - EnMean),2)/(len(varList[0][tau])-1))
            deta = varList[1][tau][item]
            dphi = varList[2][tau][item]
            phi_idx = int(math.floor((dphi + dR_max) / dPhi_width + 0.5))
            eta_idx = int(math.floor((deta + dR_max) / dEta_width + 0.5))
            CellGrid[tau][varIndex][phi_idx][eta_idx]+= dev_std
    varIndex += 1
    return varIndex
@njit()
def addCellGridSizeVar(CellGrid, nTaus, varList, varIndex,dR_min, dR_max, dPhi_width, dEta_width):
    for tau in range(0,nTaus):
        for item in range(0, len(varList[0][tau])):
            energy = varList[0][tau][item]
            deta = varList[1][tau][item]
            dphi = varList[2][tau][item]
            phi_idx = int(math.floor((dphi + dR_max) / dPhi_width + 0.5))
            eta_idx = int(math.floor((deta + dR_max) / dEta_width + 0.5))
            CellGrid[tau][varIndex][phi_idx][eta_idx]+= 1
    varIndex += 1
    return varIndex


#@njit()
def getCellGridMatrix(nVars, flatVarsMatrix, ecalVarsList, hcalVarsList, n_cellsX, n_cellsY, verbose=False):
    #print("writing feature matrix")
    nTaus = len(flatVarsMatrix[0])
    nFlatVars = len(flatVarsMatrix)
    CellGrid = np.zeros((nTaus, nVars, n_cellsX, n_cellsY))
    dR_min = -0.5
    dR_max = 0.5
    dPhi_width = (dR_max-dR_min)/(n_cellsX-1)
    dEta_width = (dR_max-dR_min)/(n_cellsY-1)
    varIndex = 0
    # 1. fill flat vars
    print("filling flat vars")
    for var in range(0, nFlatVars):
        for tau in range(0,nTaus):
            for phi_idx in range(0,n_cellsX):
                for eta_idx in range(0,n_cellsY):
                    #print(i,j,k,varIndex)
                    CellGrid[tau][varIndex][phi_idx][eta_idx] = flatVarsMatrix[var][tau]
                    #print(flatVarsMatrix[var][i])
        varIndex += 1
    # 2. ECAL vars

    print(varIndex)
    print("filling ECAL vars - Sum")
    varIndex = addCellGridSumVar(CellGrid, nTaus, ecalVarsList, 0, varIndex,dR_min, dR_max, dPhi_width, dEta_width) # sum of EcalEnergy
    assert addCellGridSumVar.nopython_signatures # this was compiled in nopython mode

    print(varIndex)
    print("filling ECAL vars - WSum")
    for i in range(1,4): # delta eta, delta phi, chi2
        varIndex = addCellGridWeightedSumVar(CellGrid, nTaus, ecalVarsList, i, varIndex,dR_min, dR_max, dPhi_width, dEta_width)
        assert addCellGridWeightedSumVar.nopython_signatures # this was compiled in nopython mode

    print(varIndex)
    print("filling ECAL vars - DevSTD")
    varIndex = addCellGridDevStdVar(CellGrid, nTaus, ecalVarsList, varIndex,dR_min, dR_max, dPhi_width, dEta_width)
    assert addCellGridDevStdVar.nopython_signatures # this was compiled in nopython mode

    print(varIndex)
    print("filling ECAL vars - Size")
    varIndex = addCellGridSizeVar(CellGrid, nTaus, ecalVarsList, varIndex,dR_min, dR_max, dPhi_width, dEta_width)
    assert addCellGridSizeVar.nopython_signatures # this was compiled in nopython mode

    # 3. HCAL vars
    '''
    print(varIndex)
    print("filling HCAL vars - Sum")
    varIndex = addCellGridSumVar(CellGrid, nTaus, hcalVarsList, 0, varIndex,dR_min, dR_max, dPhi_width, dEta_width) # sum of EcalEnergy
    assert addCellGridSumVar.nopython_signatures # this was compiled in nopython mode

    print(varIndex)
    print("filling HCAL vars - WSum")
    for i in [1,2,3]: # delta eta, delta phi, chi2
        varIndex = addCellGridWeightedSumVar(CellGrid, nTaus, hcalVarsList, i, varIndex,dR_min, dR_max, dPhi_width, dEta_width)
        assert addCellGridWeightedSumVar.nopython_signatures # this was compiled in nopython mode

    print(varIndex)
    print("filling HCAL vars - DevSTD")
    varIndex = addCellGridDevStdVar(CellGrid, nTaus, hcalVarsList, varIndex,dR_min, dR_max, dPhi_width, dEta_width)
    assert addCellGridDevStdVar.nopython_signatures # this was compiled in nopython mode

    print(varIndex)
    print("filling HCAL vars - Size")
    varIndex = addCellGridSizeVar(CellGrid, nTaus, hcalVarsList, varIndex,dR_min, dR_max, dPhi_width, dEta_width)
    assert addCellGridSizeVar.nopython_signatures # this was compiled in nopython mode

    print(varIndex)
    '''

    return CellGrid


#variablesToAdd =[ "nVertices", "l1Tau_pt", "l1Tau_eta", "l1Tau_phi", "l1Tau_hwIso", "caloRecHit_e_DeltaEta_weighted", "caloRecHit_e_DeltaPhi_weighted", "caloRecHit_e_chi2_weighted", "caloRecHit_had_DeltaEta_weighted", "caloRecHit_had_DeltaPhi_weighted","caloRecHit_had_chi2_weighted",  "patatrack_DeltaEta_weighted", "patatrack_DeltaPhi_weighted", "(patatrack_chi2_Over_patatrack_ndof)_weighted", "(patatrack_vert_z_Minus_PV_z)_weighted","(patatrack_vert_chi2/_Over_patatrack_vert_ndof)_weighted","patatrack_vert_ptv2_weighted","caloRecHit_e_energy_sum", "caloRecHit_e_energy_dev_std_sum",
variablesToAdd =[ "nVertices", "l1Tau_pt", "l1Tau_eta", "l1Tau_phi", "l1Tau_hwIso", "caloRecHit_e_energy_sum", "caloRecHit_e_DeltaEta_weighted", "caloRecHit_e_DeltaPhi_weighted", "caloRecHit_e_chi2_weighted","caloRecHit_e_energy_size"]#, "caloRecHit_had_energy_sum", "caloRecHit_had_DeltaEta_weighted", "caloRecHit_had_DeltaPhi_weighted", "caloRecHit_had_chi2_weighted","caloRecHit_had_energy_size"]#,  "patatrack_DeltaEta_weighted", "caloRecHit_had_energy_sum", #"caloRecHit_e_energy_dev_std_sum",  "caloRecHit_had_energy_dev_std _sum", "patatrack_pt_sum", "patatrack_pt_dev_std_sum", "patatrack_charge_sum", , "caloRecHit_had_energy_size", "patatrack_size_withVertex", "patatrack_size_withoutVertex"]

# (where PV is the vertex with the highest ptv2 amoung ones available for this tau)

nVars = len(variablesToAdd)
CellGrid = getCellGridMatrix(nVars, flatMatrix, ecalVarsList, hcalVarsList,  n_cellsX, n_cellsY, True)
np.save("CellGridECAL.npy",CellGrid)
#print(CellGrid[1][8][:][:])
#print(CellGrid.shape)
#assert getCellGridMatrix.nopython_signatures # this was compiled in nopython mode
#print( CellGrid)
start_3 = time.time()
if(args.time):
    print (("time to define grid {}").format(start_3-start_2))
