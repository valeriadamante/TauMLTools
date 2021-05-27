from numba import jit, njit
import numpy as np
import time
import uproot
import awkward as ak
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import enum
import math
import argparse
import os
parser = argparse.ArgumentParser()
parser.add_argument('--time', required=False, type=bool, default=False , help='if set to true, display time information')
parser.add_argument('--n_max_events', required=False, type=int, default=-1, help='max number of events to be processed')
parser.add_argument('--input_file', required=False, type=str, default='DataSetTrainingWeight.root', help='input file name')
parser.add_argument('--output_file', required=False, type=str, default='CellGrid.npy', help='output file name')
parser.add_argument('--outputNorm_file', required=False, type=str, default='CellGridNorm.npy', help='output file name')
parser.add_argument('--input_tuple', required=False, type=str, default='L2TauTrainTuple', help='input tree name')
parser.add_argument('--n_cellsX', required=False, type=int, default=5, help='number of cells along X dir')
parser.add_argument('--n_cellsY', required=False, type=int, default=5, help='number of cells along Y dir')
parser.add_argument('--saveNormGrid', required=False, type=bool, default=False , help='if set to true, save normalized cell grid')
parser.add_argument('--verbose', required=False, type=int, default=0)
args = parser.parse_args()

class NNInputs (enum.IntEnum):
    nVertices = 0
    l1Tau_pt = 1
    l1Tau_eta = 2
    l1Tau_hwIso = 3
    EcalEnergySum = 4
    EcalSize = 5
    EcalEnergyStdDev = 6
    EcalDeltaEta = 7
    EcalDeltaPhi = 8
    EcalChi2 = 9
    EcalEnergySumForPositiveChi2 = 10
    EcalSizeForPositiveChi2 = 11
    HcalEnergySum = 12
    HcalSize = 13
    HcalEnergyStdDev = 14
    HcalDeltaEta = 15
    HcalDeltaPhi = 16
    HcalChi2 = 17
    HcalEnergySumForPositiveChi2 = 18
    HcalSizeForPositiveChi2 = 19
    PatatrackPtSum = 20
    PatatrackSize = 21
    PatatrackSizeWithVertex = 22
    PatatrackPtSumWithVertex = 23
    PatatrackChargeSum = 24
    PatatrackDeltaEta = 25
    PatatrackDeltaPhi = 26
    PatatrackChi2OverNdof = 27
    PatatrackNdof = 28
    PatatrackDxy = 29
    PatatrackDz = 30
    #PatavertPtv2 = 32
    #PatavertChi2OverNdof = 33
    #PatavertNdof = 34
    #PatavertZ = 35


dict = {}

def StandardizeSingleVar(CellGrid, varPos, varName, min=None, max=None, mean=None, std=None, mask=np.array(None), **kwargs):
    start= time.time()
    outDir = kwargs['outDir']
    if not os.path.exists(outDir) :
        os.mkdir(outDir)
        if(kwargs['verbose']>1):
            print(("created {} directory").format(outDir))
    else:
        if(kwargs['verbose']>1):
            print(("{} directory already exists").format(outDir))
    if mask.any() == None:
        mask = True
    if(kwargs['verbose']>1):
        print(("cell grid of first elements \n {}").format(CellGrid[1,:,:,varPos]))
    maArray = np.ma.MaskedArray(CellGrid[:,:,:,varPos],mask=1-mask)
    if(kwargs['verbose']>1):
        print(("masked array of first elements \n {}").format(maArray[1,:,:]))
    if mean ==None :
        #mean = np.ma.MaskedArray.round((maArray).mean(),4)
        mean = maArray.mean()
        if(kwargs['verbose']>1):
            print(("mean before put 4 float = {}").format(mean))
    if std == None :
        std = maArray.std()
        if(kwargs['verbose']>1):
            print(("std before put 4 float = {}").format(std))
    maArrayStd = np.ma.MaskedArray(((CellGrid[:,:,:,varPos]-mean)/std),mask=1-mask)
    if(kwargs['verbose']>1):
        print(("masked std array of first elements \n {}").format(maArrayStd[1,:,:]))
    if min==None :
        min = np.ma.MaskedArray.min(maArrayStd)
        if(kwargs['verbose']>1):
            print(("min before put 4 float = {}").format(min))
    if max == None :
        max = np.ma.MaskedArray.max(maArrayStd)
        if(kwargs['verbose']>1):
            print(("max before put 4 float = {}").format(max))
    mean = np.round(mean, 4)
    std = np.round(std, 4)
    min = np.round(min, 4)
    max = np.round(max, 4)
    if(kwargs['verbose']>0):
        print(("mean = {}, std = {}, min = {}, max = {}").format(mean, std, min, max))
    minVar = np.ma.MaskedArray.min(maArray)
    maxVar = np.ma.MaskedArray.max(maArray)
    minVar = np.round(minVar, 4)
    maxVar = np.round(maxVar, 4)
    if(kwargs['verbose']>0):
        print(("var {}, minVar {}, maxVar {}").format(varName, minVar,maxVar))
    if(kwargs['plot']):
        nbins = kwargs["nBins"]
        bins=np.linspace(minVar,maxVar,nbins)
        if kwargs["logX"]==True:
            plt.xscale('log')
        if kwargs["logY"]==True:
            plt.yscale('log')
        fig1 = plt.figure()
        plt.title(varName)
        if kwargs["logX"]==True:
            plt.xscale('log')
        if kwargs["logY"]==True:
            plt.yscale('log')
        low_high = (minVar, maxVar)
        plt.hist(CellGrid[:,:,:,varPos].flatten(), color='b', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True)
        plt.title(('{} before normalizing').format(varName))
        plt.xlabel(varName)
        plt.ylabel("events")
        text = [(("mean = {}\nstd = {}\nminVar = {}\nmaxVar = {}").format(mean, std, minVar, maxVar))]
        plt.legend(text, loc="upper right")
        fig1.savefig(('{}/{}_beforeNormalizing.pdf').format(outDir, varName))
    CellGrid[:,:,:,varPos]= np.clip(np.where(mask, ((CellGrid[:,:,:,varPos]-mean)/std), 0 ), min, max)
    dict[varName]={"mean ": mean, "std ": std,  "min ": min,  "max ": max }
    if(kwargs['verbose']>0):
        print(dict)
    if(kwargs['plot']):
        fig2 = plt.figure()
        plt.title(varName)
        if kwargs["logX"]==True:
            plt.xscale('log')
        if kwargs["logY"]==True:
            plt.yscale('log')
        bins=np.linspace(min,max,nbins)
        low_high = (min, max)
        plt.hist(CellGrid[:,:,:,varPos].flatten(), color='b', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True)
        plt.title(('{} after normalizing').format(varName))
        plt.xlabel(varName)
        plt.ylabel("events")
        text = [(("mean = {}\nstd = {}\nmin = {}\nmax = {}").format(mean, std, min, max))]
        #plt.text(maxVar*0.8,, text)
        plt.legend(text, loc="upper right")
        #plt.savefig(('{}/{}_beforeNormalizing.pdf').format(outDir, varName))
        fig2.savefig(('{}/{}_afterNormalizing.pdf').format(outDir, varName))
    final_time = time.time()-start
    if(kwargs["time"]==True):
        print(("time to fill {} = {} ").format(varName, final_time))
    plt.close('all')



def StandardizeVars(CellGrid, verbose=0, timeInfo=True):
    kwArgs = {'outDir': 'outDistributions',  'logX':False, 'logY':False, 'nBins':50, 'verbose':verbose, 'time':timeInfo, 'plot':True}
    # nVertices
    mean = int(round(CellGrid[:,:,:,np.intp(NNInputs.nVertices)].mean(),0))
    std = int(round(CellGrid[:,:,:,np.intp(NNInputs.nVertices)].std(),0))
    StandardizeSingleVar(CellGrid, NNInputs.nVertices, "nVertices", -5, 5, mean, std,  **kwArgs)
    # l1Pt  --> rescaling to l1Pt max
    l1Ptmax = 256
    StandardizeSingleVar(CellGrid, NNInputs.l1Tau_pt, "l1Tau_pt", 0, 1, 0, l1Ptmax,  **kwArgs)
    # l1 Eta  --> rescaling from [-a, a] to [-1,1]
    minEta = (CellGrid[:,:,:,np.intp(NNInputs.l1Tau_eta)]).min()
    maxEta = (CellGrid[:,:,:,np.intp(NNInputs.l1Tau_eta)]).max()
    std = (maxEta-minEta)/2
    mean = 0
    StandardizeSingleVar(CellGrid, NNInputs.l1Tau_eta, "l1Tau_eta", -1, 1, mean, std, **kwArgs)
    # l1Tau hwIso --> no rescaling
    StandardizeSingleVar(CellGrid, NNInputs.l1Tau_hwIso, "l1Tau_hwIso", 0, 1, 0, 1, **kwArgs)
    # Ecal EnergySum --> rescaling to l1Pt max
    enMask =  CellGrid[:,:,:,np.intp(NNInputs.EcalEnergySum)]>0
    kwArgs['logY'] = True
    StandardizeSingleVar(CellGrid, NNInputs.EcalEnergySum, "EcalEnergySum", 0, 5, 0, l1Ptmax, enMask, **kwArgs)
    # Ecal EnergyStdDev --> rescaling to l1Pt max
    StandardizeSingleVar(CellGrid, NNInputs.EcalEnergyStdDev, "EcalEnergyStdDev", -5, 5, None, None, enMask, **kwArgs)
    # Ecal Size
    kwArgs['logY'] = False
    mean = round(CellGrid[:,:,:,NNInputs.EcalSize].mean(),0)
    std = round(CellGrid[:,:,:,NNInputs.EcalSize].std(),0)
    StandardizeSingleVar(CellGrid, NNInputs.EcalSize, "EcalSize", 0, 10, 0, std, enMask, **kwArgs)
    # Ecal Delta Eta --> take in range from [-0.5, 0.5] to [-1,1]
    minEta = (CellGrid[:,:,:,np.intp(NNInputs.EcalDeltaEta)]).min()
    maxEta = (CellGrid[:,:,:,np.intp(NNInputs.EcalDeltaEta)]).max()
    std = (maxEta-minEta)/2
    mean = 0
    StandardizeSingleVar(CellGrid, NNInputs.EcalDeltaEta, "EcalDeltaEta",  -1,1, mean, std, enMask, **kwArgs)
    # Ecal Delta Phi --> take in range from [-0.5, 0.5] to [-1,1]
    minPhi = (CellGrid[:,:,:,np.intp(NNInputs.EcalDeltaPhi)]).min()
    maxPhi = (CellGrid[:,:,:,np.intp(NNInputs.EcalDeltaPhi)]).max()
    std = (maxPhi-minPhi)/2
    mean = 0
    StandardizeSingleVar(CellGrid, NNInputs.EcalDeltaPhi, "EcalDeltaPhi",  -1,1, mean, std, enMask, **kwArgs)
    # Ecal Chi2 --> NEW MASK
    chi2Mask =  CellGrid[:,:,:,np.intp(NNInputs.EcalChi2)]>0
    StandardizeSingleVar(CellGrid, NNInputs.EcalChi2, "EcalChi2",  -5, 5, None, None, chi2Mask, **kwArgs)
    # Ecal EcalEnergySumForPositiveChi2  --> rescaling to l1Pt max with Chi2 mask
    kwArgs['logY'] = True
    StandardizeSingleVar(CellGrid, NNInputs.EcalEnergySumForPositiveChi2, "EcalEnergySumForPositiveChi2",  0, 5, 0, l1Ptmax, chi2Mask, **kwArgs)
    kwArgs['logY'] = False
    # ecal energy size for positive chi2 --> with Chi2 mask
    mean = round(CellGrid[:,:,:,NNInputs.EcalSizeForPositiveChi2].mean(),0)
    std = round(CellGrid[:,:,:,NNInputs.EcalSizeForPositiveChi2].std(),0)
    StandardizeSingleVar(CellGrid, NNInputs.EcalSizeForPositiveChi2, "EcalSizeForPositiveChi2",  0, 10, 0, std, chi2Mask, **kwArgs)
    # Hcal EnergySum --> rescaling to l1Pt max
    enMask =  CellGrid[:,:,:,np.intp(NNInputs.HcalEnergySum)]>0
    kwArgs['logY'] = True
    StandardizeSingleVar(CellGrid, NNInputs.HcalEnergySum, "HcalEnergySum", 0, 5, 0, l1Ptmax, enMask, **kwArgs)
    # Hcal EnergyStdDev --> rescaling to l1Pt max
    StandardizeSingleVar(CellGrid, NNInputs.HcalEnergyStdDev, "HcalEnergyStdDev", -5, 5, None, None, enMask, **kwArgs)
    # Hcal Size
    kwArgs['logY'] = False
    mean = round(CellGrid[:,:,:,NNInputs.HcalSize].mean(),0)
    std = round(CellGrid[:,:,:,NNInputs.HcalSize].std(),0)
    StandardizeSingleVar(CellGrid, NNInputs.HcalSize, "HcalSize", 0, 10, 0, std, enMask, **kwArgs)
    # Hcal Delta Eta --> take in range from [-0.5, 0.5] to [-1,1]
    minEta = (CellGrid[:,:,:,np.intp(NNInputs.HcalDeltaEta)]).min()
    maxEta = (CellGrid[:,:,:,np.intp(NNInputs.HcalDeltaEta)]).max()
    std = (maxEta-minEta)/2
    mean = 0
    StandardizeSingleVar(CellGrid, NNInputs.HcalDeltaEta, "HcalDeltaEta",  -1,1, mean, std, enMask, **kwArgs)
    # Hcal Delta Phi --> take in range from [-0.5, 0.5] to [-1,1]
    minPhi = (CellGrid[:,:,:,np.intp(NNInputs.HcalDeltaPhi)]).min()
    maxPhi = (CellGrid[:,:,:,np.intp(NNInputs.HcalDeltaPhi)]).max()
    std = (maxPhi-minPhi)/2
    mean = 0
    StandardizeSingleVar(CellGrid, NNInputs.HcalDeltaPhi, "HcalDeltaPhi",  -1,1, mean, std, enMask, **kwArgs)
    # Hcal Chi2 --> NEW MASK
    chi2Mask =  CellGrid[:,:,:,np.intp(NNInputs.HcalChi2)]>0
    StandardizeSingleVar(CellGrid, NNInputs.HcalChi2, "HcalChi2",  -5, 5, None, None, chi2Mask, **kwArgs)
    # Hcal HcalEnergySumForPositiveChi2  --> rescaling to l1Pt max with Chi2 mask
    kwArgs['logY'] = True
    StandardizeSingleVar(CellGrid, NNInputs.HcalEnergySumForPositiveChi2, "HcalEnergySumForPositiveChi2",  0, 5, 0, l1Ptmax, chi2Mask, **kwArgs)
    kwArgs['logY'] = False
    # hcal energy size for positive chi2 --> with Chi2 mask
    mean = round(CellGrid[:,:,:,NNInputs.HcalSizeForPositiveChi2].mean(),0)
    std = round(CellGrid[:,:,:,NNInputs.HcalSizeForPositiveChi2].std(),0)
    StandardizeSingleVar(CellGrid, NNInputs.HcalSizeForPositiveChi2, "HcalSizeForPositiveChi2",  0, 10, 0, std, chi2Mask, **kwArgs)
    # Patatrack PtSum  --> rescaling to l1Pt max
    enMask =  CellGrid[:,:,:,np.intp(NNInputs.PatatrackPtSum)]>0
    kwArgs['logY'] = True
    StandardizeSingleVar(CellGrid, NNInputs.PatatrackPtSum, "PatatrackPtSum", 0, 5, 0, l1Ptmax, enMask, **kwArgs)
    # Patatrack PtSum WithVertex --> rescaling to l1Pt max
    StandardizeSingleVar(CellGrid, NNInputs.PatatrackPtSumWithVertex, "PatatrackPtSumWithVertex", 0, 5, 0, l1Ptmax, enMask, **kwArgs)
    kwArgs['logY'] = False
    # Patatrack Size
    mean = round(CellGrid[:,:,:,NNInputs.PatatrackSize].mean(),0)
    std = round(CellGrid[:,:,:,NNInputs.PatatrackSize].std(),0)
    StandardizeSingleVar(CellGrid, NNInputs.PatatrackSize, "PatatrackSize", 0, 10, 0, std, enMask, **kwArgs)
    # Patatrack Charge Sum
    StandardizeSingleVar(CellGrid, NNInputs.PatatrackChargeSum, "PatatrackChargeSum", -5, 5, None, None, enMask, **kwArgs)
    # Patatrack Size With Vertex
    mean = round(CellGrid[:,:,:,NNInputs.PatatrackSizeWithVertex].mean(),0)
    std = round(CellGrid[:,:,:,NNInputs.PatatrackSizeWithVertex].std(),0)
    StandardizeSingleVar(CellGrid, NNInputs.PatatrackSizeWithVertex, "PatatrackSizeWithVertex", 0, 10, 0, std, enMask, **kwArgs)
    # Patatrack Ndof
    StandardizeSingleVar(CellGrid, NNInputs.PatatrackNdof, "PatatrackNdof", -5, 5, None, None, enMask, **kwArgs)
    # Patatrack Delta Eta --> take in range from [-0.5, 0.5] to [-1,1]
    minEta = (CellGrid[:,:,:,np.intp(NNInputs.PatatrackDeltaEta)]).min()
    maxEta = (CellGrid[:,:,:,np.intp(NNInputs.PatatrackDeltaEta)]).max()
    std = (maxEta-minEta)/2
    mean = 0
    StandardizeSingleVar(CellGrid, NNInputs.PatatrackDeltaEta, "PatatrackDeltaEta",  -1,1, mean, std, enMask, **kwArgs)
    # Patatrack Delta Phi --> take in range from [-0.5, 0.5] to [-1,1]
    minPhi = (CellGrid[:,:,:,np.intp(NNInputs.PatatrackDeltaPhi)]).min()
    maxPhi = (CellGrid[:,:,:,np.intp(NNInputs.PatatrackDeltaPhi)]).max()
    std = (maxPhi-minPhi)/2
    mean = 0
    StandardizeSingleVar(CellGrid, NNInputs.PatatrackDeltaPhi, "PatatrackDeltaPhi",  -1,1, mean, std, enMask, **kwArgs)
    # Patatrack Chi2 --> new mask
    chi2Mask =  CellGrid[:,:,:,np.intp(NNInputs.PatatrackChi2OverNdof)]>0
    StandardizeSingleVar(CellGrid, NNInputs.PatatrackChi2OverNdof, "PatatrackChi2OverNdof",  -5, 5, None, None, chi2Mask, **kwArgs)
    # patatrack Dxy --> new mask for size > 0
    sizeMask =  CellGrid[:,:,:,np.intp(NNInputs.PatatrackSize)]>0
    StandardizeSingleVar(CellGrid, NNInputs.PatatrackDxy, "PatatrackDxy",  -5, 5, None, None, sizeMask, **kwArgs)
    # patatrack Dz
    sizeMask =  CellGrid[:,:,:,np.intp(NNInputs.PatatrackSize)]>0
    StandardizeSingleVar(CellGrid, NNInputs.PatatrackDz, "PatatrackDz",  -5, 5, None, None, sizeMask, **kwArgs)



@jit(nopython=True)
def getCellGridMatrix(nVars, n_cellsX, n_cellsY, nVertices,
                        l1Tau_pt, l1Tau_eta, l1Tau_hwIso,
                        caloRecHit_e_energy, caloRecHit_e_DeltaEta, caloRecHit_e_DeltaPhi, caloRecHit_e_chi2,
                        caloRecHit_had_energy, caloRecHit_had_DeltaEta, caloRecHit_had_DeltaPhi, caloRecHit_had_chi2,
                        patatrack_pt, patatrack_DeltaEta, patatrack_DeltaPhi, patatrack_chi2, patatrack_ndof, patatrack_charge, patatrack_dxy, patatrack_dz, patatrack_hasVertex, patatrack_vert_z, patatrack_vert_ptv2, patatrack_vert_chi2, patatrack_vert_ndof):
    nTaus = len(nVertices)
    CellGrid = np.zeros((nTaus, n_cellsX, n_cellsY, nVars))
    dR_min = -0.5
    dR_max = 0.5
    dPhi_width = (dR_max-dR_min)/(n_cellsX)
    dEta_width = (dR_max-dR_min)/(n_cellsY)
    for tau in range(nTaus):
        # 1. fill flat vars
        CellGrid[tau,:,:,np.intp(NNInputs.nVertices)] = nVertices[tau]
        CellGrid[tau,:,:,np.intp(NNInputs.l1Tau_pt)] = l1Tau_pt[tau]
        CellGrid[tau,:,:,np.intp(NNInputs.l1Tau_eta)] = l1Tau_eta[tau]
        CellGrid[tau,:,:,np.intp(NNInputs.l1Tau_hwIso)] = l1Tau_hwIso[tau]
        # 2. filling ecal
        nEcal = len(caloRecHit_e_energy[tau])
        for item in range(nEcal):
            deta = caloRecHit_e_DeltaEta[tau][item]
            dphi = caloRecHit_e_DeltaPhi[tau][item]
            phi_idx = int(math.floor((dphi + dR_max) / dPhi_width + 0.5))
            eta_idx = int(math.floor((deta + dR_max) / dEta_width + 0.5))
            CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.EcalEnergySum)]+= caloRecHit_e_energy[tau][item]
            CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.EcalSize)]+= 1
            CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.EcalEnergyStdDev)]+= (caloRecHit_e_energy[tau][item])**2
            CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.EcalDeltaEta)] += caloRecHit_e_DeltaEta[tau][item] * caloRecHit_e_energy[tau][item]
            CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.EcalDeltaPhi)] += caloRecHit_e_DeltaPhi[tau][item] * caloRecHit_e_energy[tau][item]
            if(caloRecHit_e_chi2[tau][item]>=0):
                CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.EcalChi2)]+=  caloRecHit_e_chi2[tau][item]*caloRecHit_e_energy[tau][item]
                CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.EcalEnergySumForPositiveChi2)]+=  caloRecHit_e_energy[tau][item]
                CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.EcalSizeForPositiveChi2)]+= 1
        # 2.1 calculate energy std dev
        CellGrid[tau,:,:,np.intp(NNInputs.EcalEnergyStdDev)] = np.where( CellGrid[tau,:,:,np.intp(NNInputs.EcalSize)]>1 , ( CellGrid[tau,:,:,np.intp(NNInputs.EcalEnergyStdDev)] - ((CellGrid[tau,:,:,np.intp(NNInputs.EcalEnergySum)])**2/CellGrid[tau,:,:,np.intp(NNInputs.EcalSize)]) ) /(CellGrid[tau,:,:,np.intp(NNInputs.EcalSize)]-1) , 0)

        # 2.2 weight variables to EcalEnergySum
        CellGrid[tau,:,:,np.intp(NNInputs.EcalDeltaEta)] = np.where(CellGrid[tau,:,:,np.intp(NNInputs.EcalEnergySum)]>0, CellGrid[tau,:,:,np.intp(NNInputs.EcalDeltaEta)]/CellGrid[tau,:,:,np.intp(NNInputs.EcalEnergySum)], CellGrid[tau,:,:,np.intp(NNInputs.EcalDeltaEta)])

        CellGrid[tau,:,:,np.intp(NNInputs.EcalDeltaPhi)] = np.where(CellGrid[tau,:,:,np.intp(NNInputs.EcalEnergySum)]>0, CellGrid[tau,:,:,np.intp(NNInputs.EcalDeltaPhi)]/CellGrid[tau,:,:,np.intp(NNInputs.EcalEnergySum)], CellGrid[tau,:,:,np.intp(NNInputs.EcalDeltaPhi)])

        CellGrid[tau,:,:,np.intp(NNInputs.EcalChi2)] = np.where(CellGrid[tau,:,:,np.intp(NNInputs.EcalEnergySumForPositiveChi2)]>0,  CellGrid[tau,:,:,np.intp(NNInputs.EcalChi2)]/CellGrid[tau,:,:,np.intp(NNInputs.EcalEnergySumForPositiveChi2)], CellGrid[tau,:,:,np.intp(NNInputs.EcalChi2)])

        # 3. filling Hcal
        nHcal = len(caloRecHit_had_energy[tau])
        for item in range(nHcal):
          deta = caloRecHit_had_DeltaEta[tau][item]
          dphi = caloRecHit_had_DeltaPhi[tau][item]
          phi_idx = int(math.floor((dphi + dR_max) / dPhi_width + 0.5))
          eta_idx = int(math.floor((deta + dR_max) / dEta_width + 0.5))
          CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.HcalEnergySum)]+= caloRecHit_had_energy[tau][item]
          CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.HcalSize)]+= 1
          CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.HcalEnergyStdDev)]+= (caloRecHit_had_energy[tau][item])**2
          CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.HcalDeltaEta)] += caloRecHit_had_DeltaEta[tau][item] * caloRecHit_had_energy[tau][item]
          CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.HcalDeltaPhi)] += caloRecHit_had_DeltaPhi[tau][item] * caloRecHit_had_energy[tau][item]
          if(caloRecHit_had_chi2[tau][item]>=0):
              CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.HcalChi2)]+=  caloRecHit_had_chi2[tau][item]*caloRecHit_had_energy[tau][item]
              CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.HcalEnergySumForPositiveChi2)]+=  caloRecHit_had_energy[tau][item]
              CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.HcalSizeForPositiveChi2)]+= 1
        # 3.1 calculate energy std dev
        CellGrid[tau,:,:,np.intp(NNInputs.HcalEnergyStdDev)] = np.where( CellGrid[tau,:,:,np.intp(NNInputs.HcalSize)]>1 , ( CellGrid[tau,:,:,np.intp(NNInputs.HcalEnergyStdDev)] - ((CellGrid[tau,:,:,np.intp(NNInputs.HcalEnergySum)])**2/CellGrid[tau,:,:,np.intp(NNInputs.HcalSize)]) ) /(CellGrid[tau,:,:,np.intp(NNInputs.HcalSize)]-1) , 0)

        # 3.2 weight variables to HcalEnergySum
        CellGrid[tau,:,:,np.intp(NNInputs.HcalDeltaEta)] = np.where(CellGrid[tau,:,:,np.intp(NNInputs.HcalEnergySum)]>0, CellGrid[tau,:,:,np.intp(NNInputs.HcalDeltaEta)]/CellGrid[tau,:,:,np.intp(NNInputs.HcalEnergySum)], CellGrid[tau,:,:,np.intp(NNInputs.HcalDeltaEta)])

        CellGrid[tau,:,:,np.intp(NNInputs.HcalDeltaPhi)] = np.where(CellGrid[tau,:,:,np.intp(NNInputs.HcalEnergySum)]>0, CellGrid[tau,:,:,np.intp(NNInputs.HcalDeltaPhi)]/CellGrid[tau,:,:,np.intp(NNInputs.HcalEnergySum)], CellGrid[tau,:,:,np.intp(NNInputs.HcalDeltaPhi)])

        CellGrid[tau,:,:,np.intp(NNInputs.HcalChi2)] = np.where(CellGrid[tau,:,:,np.intp(NNInputs.HcalEnergySumForPositiveChi2)]>0,  CellGrid[tau,:,:,np.intp(NNInputs.HcalChi2)]/CellGrid[tau,:,:,np.intp(NNInputs.HcalEnergySumForPositiveChi2)], CellGrid[tau,:,:,np.intp(NNInputs.HcalChi2)])


        # 4. filling patatrack
        nPatatrack = len(patatrack_pt[tau])
        for item in range(nPatatrack):
          deta = patatrack_DeltaEta[tau][item]
          dphi = patatrack_DeltaPhi[tau][item]
          phi_idx = int(math.floor((dphi + dR_max) / dPhi_width + 0.5))
          eta_idx = int(math.floor((deta + dR_max) / dEta_width + 0.5))
          CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.PatatrackPtSum)]+= patatrack_pt[tau][item]
          CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.PatatrackSize)]+= 1
          CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.PatatrackChargeSum)]+= patatrack_charge[tau][item]
          CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.PatatrackDeltaEta)] += patatrack_DeltaEta[tau][item] * patatrack_pt[tau][item]
          CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.PatatrackDeltaPhi)] += patatrack_DeltaPhi[tau][item] * patatrack_pt[tau][item]
          CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.PatatrackDz)] += patatrack_dz[tau][item] * patatrack_pt[tau][item]
          CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.PatatrackDxy)] += patatrack_dxy[tau][item] * patatrack_pt[tau][item]
          CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.PatatrackNdof)] += patatrack_ndof[tau][item] * patatrack_pt[tau][item]
          if(patatrack_hasVertex[tau][item]>0):
              CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.PatatrackSizeWithVertex)]+= 1
              CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.PatatrackPtSumWithVertex)]+= patatrack_pt[tau][item]
          if(patatrack_ndof[tau][item]>0):
              CellGrid[tau,phi_idx,eta_idx,np.intp(NNInputs.PatatrackChi2OverNdof)]+= patatrack_chi2[tau][item] * patatrack_pt[tau][item] / patatrack_ndof[tau][item]

        # normalize variables
        CellGrid[tau,:,:,np.intp(NNInputs.PatatrackDeltaEta)] = np.where(CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)]>0, CellGrid[tau,:,:,np.intp(NNInputs.PatatrackDeltaEta)]/CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)], CellGrid[tau,:,:,np.intp(NNInputs.PatatrackDeltaEta)])

        CellGrid[tau,:,:,np.intp(NNInputs.PatatrackDeltaPhi)] = np.where(CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)]>0, CellGrid[tau,:,:,np.intp(NNInputs.PatatrackDeltaPhi)]/CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)], CellGrid[tau,:,:,np.intp(NNInputs.PatatrackDeltaPhi)])

        CellGrid[tau,:,:,np.intp(NNInputs.PatatrackDxy)] = np.where(CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)]>0, CellGrid[tau,:,:,np.intp(NNInputs.PatatrackDxy)]/CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)], CellGrid[tau,:,:,np.intp(NNInputs.PatatrackDxy)])

        CellGrid[tau,:,:,np.intp(NNInputs.PatatrackDz)] = np.where(CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)]>0, CellGrid[tau,:,:,np.intp(NNInputs.PatatrackDz)]/CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)], CellGrid[tau,:,:,np.intp(NNInputs.PatatrackDz)])

        CellGrid[tau,:,:,np.intp(NNInputs.PatatrackNdof)] = np.where(CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)]>0, CellGrid[tau,:,:,np.intp(NNInputs.PatatrackNdof)]/CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)], CellGrid[tau,:,:,np.intp(NNInputs.PatatrackNdof)])

        #CellGrid[tau,:,:,np.intp(NNInputs.PatavertPtv2)] = np.where(CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)]>0, CellGrid[tau,:,:,np.intp(NNInputs.PatavertPtv2)]/CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)], CellGrid[tau,:,:,np.intp(NNInputs.PatavertPtv2)])

        CellGrid[tau,:,:,np.intp(NNInputs.PatatrackChi2OverNdof)] = np.where(CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)]>0,  CellGrid[tau,:,:,np.intp(NNInputs.PatatrackChi2OverNdof)]/CellGrid[tau,:,:,np.intp(NNInputs.PatatrackPtSum)], CellGrid[tau,:,:,np.intp(NNInputs.PatatrackChi2OverNdof)])

    return CellGrid



start_0 = time.time()
absolute_path='/Users/valeriadamante/Desktop/Dottorato/L2SkimmedTuples/DataSetTraining/'
#absolute_path='/home/valeria/'
inFile_name = args.input_file
inFile = absolute_path+inFile_name
treeName = args.input_tuple
outFile_name =  args.output_file
outFileNorm_name =  args.outputNorm_file
outFile = absolute_path+outFile_name
outFileNorm = absolute_path+outFileNorm_name
n_cellsX = args.n_cellsX
n_cellsY = args.n_cellsY
nVars = len(NNInputs)
if(not os.path.exists(outFile)):
    tree = uproot.open(inFile)[treeName]
    start_1 = time.time()
    if(args.time):
        print (("time to get file {}").format(start_1-start_0))
    if(args.n_max_events!=-1):
        awkArray =tree.arrays(entry_stop=args.n_max_events)
    else:
        awkArray =tree.arrays()
    start_2 = time.time()
    if(args.time):
        print (("time to get awkarray {}").format(start_2-start_1))
    CellGrid = getCellGridMatrix(nVars, n_cellsX, n_cellsY, awkArray.nVertices,
                            awkArray.l1Tau_pt, awkArray.l1Tau_eta, awkArray.l1Tau_hwIso, awkArray.caloRecHit_e_energy, awkArray.caloRecHit_e_DeltaEta, awkArray.caloRecHit_e_DeltaPhi, awkArray.caloRecHit_e_chi2, awkArray.caloRecHit_had_energy, awkArray.caloRecHit_had_DeltaEta, awkArray.caloRecHit_had_DeltaPhi, awkArray.caloRecHit_had_chi2, awkArray.patatrack_pt, awkArray.patatrack_DeltaEta, awkArray.patatrack_DeltaPhi, awkArray.patatrack_chi2, awkArray.patatrack_ndof, awkArray.patatrack_charge, awkArray.patatrack_dxy, awkArray.patatrack_dz, awkArray.patatrack_hasVertex, awkArray.patatrack_vert_z, awkArray.patatrack_vert_ptv2, awkArray.patatrack_vert_chi2, awkArray.patatrack_vert_ndof)

    np.save(outFile,CellGrid)
    print(CellGrid.shape)
    start_3 = time.time()
    if(args.time):
        print (("time to define grid {}").format(start_3-start_2))
    StandardizeVars(CellGrid, args.verbose, args.time)
else:
    CellGrid = np.load(outFile)
    start_1 = time.time()
    if(args.time):
        print (("time to get file {}").format(start_1-start_0))
    #print(CellGrid.shape)
    StandardizeVars(CellGrid, args.verbose, args.time)
    if(args.saveNormGrid == True):
        np.save(outFileNorm,CellGrid)
