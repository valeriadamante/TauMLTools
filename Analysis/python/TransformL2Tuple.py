import os
import argparse

parent_directory =  "/gpfs/ddn/srm/cms/store/user/vdamante/L2TauTuple/"
parent_directory_2 = "/gpfs/ddn/srm/cms/store/user/vdamante/L2TausTuple/" # EphemeralHLTPhysics3 and ZprimeToTauTau
files_directories = {
"EphemeralHLTPhysics1": "crab_EphemeralHLTPhysics1",
"EphemeralHLTPhysics2": "crab_EphemeralHLTPhysics2",
"EphemeralHLTPhysics3": "crab_EphemeralHLTPhysics3",
"EphemeralHLTPhysics4": "crab_EphemeralHLTPhysics4",
"EphemeralHLTPhysics5": "crab_EphemeralHLTPhysics5",
"EphemeralHLTPhysics6": "crab_EphemeralHLTPhysics6",
"EphemeralHLTPhysics7": "crab_EphemeralHLTPhysics7",
"EphemeralHLTPhysics8": "crab_EphemeralHLTPhysics8",
"DYToLL_M-50_TuneCP5_14TeV-pythia8":"crab_DYToLL_M-50",
"QCD_Pt_120to170_TuneCP5_14TeV_pythia8":"crab_QCD_Pt_120to170",
"QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8":"crab_QCD_Pt-15to3000_Flat",
"QCD_Pt-15to7000_TuneCP5_Flat_14TeV_pythia8":"crab_QCD_Pt-15to7000_Flat",
"QCD_Pt-15to7000_TuneCP5_Flat_14TeV_pythia8":"crab_QCD_Pt-15to7000_Flat_ext1",
"QCD_Pt_170to300_TuneCP5_14TeV_pythia8":"crab_QCD_Pt_170to300",
"QCD_Pt_300to470_TuneCP5_14TeV_pythia8":"crab_QCD_Pt_300to470",
"QCD_Pt_30to50_TuneCP5_14TeV_pythia8":"crab_QCD_Pt_30to50",
"QCD_Pt_470to600_TuneCP5_14TeV_pythia8":"crab_QCD_Pt_470to600",
"QCD_Pt_50to80_TuneCP5_14TeV_pythia8":"crab_QCD_Pt_50to80",
"QCD_Pt_600oInf_TuneCP5_14TeV_pythia8":"crab_QCD_Pt_600oInf",
"QCD_Pt_80to120_TuneCP5_14TeV_pythia8":"crab_QCD_Pt_80to120",
"TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8":"crab_TTToSemiLeptonic",
"TT_TuneCP5_14TeV-powheg-pythia8":"crab_TT",
"TT_TuneCP5_14TeV-powheg-pythia8":"crab_TT_ext1",
"VBFHToTauTau_M125_TuneCUETP8M1_14TeV_powheg_pythia8":"crab_VBFHToTauTau",
"WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8":"crab_WJetsToLNu",
"ZprimeToTauTau_M-4000_TuneCP5_14TeV-pythia8-tauola":"crab_ZprimeToTauTau_M-4000",
}



parser = argparse.ArgumentParser()
#parser.add_argument('machine', type=str, default="local", choices=["lxplus", "local"])
parser.add_argument('--sample', required=True, type=str, help= "input reference file", choices=['DY', 'QCD', 'TT', 'VBF', 'WJets', 'ZPrime'])
parser.add_argument('--n-threads', required=False, type=int, default=1, help='Number of threads')
parser.add_argument('--start-entry', required=False, type=int, default=0 , help='start entry')
parser.add_argument('--end-entry', required=False, type=int, default=-1, help='end entry')

args = parser.parse_args()

cmd = 'L2TrainingTupleProducer --input "VBFHToTauTau_M125_TuneCUETP8M1_14TeV_powheg_pythia8_output1.root" --output "VBFProvaoutput_2.root" --isQCD false'
inputs =[]
for k in files_directories:
    if(args.sample in k):
        parent = parent_directory
        if(k == "ZprimeToTauTau_M-4000_TuneCP5_14TeV-pythia8-tauola" or k == "EphemeralHLTPhysics3" ):
            parent = parent_directory_2
        path = parent + k + "/"+files_directories[k]+"/"
        for file in os.listdir( path ):
            print path+file
            #inputs.append()
