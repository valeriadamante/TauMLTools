#count average Taus
import ROOT
import numpy as np
import pandas as pd
import root_pandas
#pd.set_option('display.max_columns', 500)
import sys
import argparse
import os
from TauMLTools.Analysis.getPathNames import *
#tt = root_pandas.read_root(inpath + tt_file, 'events_3j1t')#, where='Jet_score_best>0.7 || Jet_score_secondbest>0.7')
#wjets = root_pandas.read_root(inpath + wj_file, 'events_3j1t')#, where='Jet_score_best>0.7 || Jet_score_secondbest>0.7')
#stbb = root_pandas.read_root(inpath + stbb_file, 'events_3j1t' )#,where='Jet_score_best>0.7 || Jet_score_secondbest>0.7')
#stbq = root_pandas.read_root(inpath+stbq_file, 'events_3j1t', where='Jet_score_best>0.7 || Jet_score_secondbest>0.7')
parent_directory =  "/gpfs/ddn/srm/cms/store/user/vdamante/L2TauTuple/"
parent_directory_2 = "/gpfs/ddn/srm/cms/store/user/vdamante/L2TausTuple/"
files_directories = {
    "Data" : {
        "EphemeralHLTPhysics1":[ "crab_EphemeralHLTPhysics1"],
        "EphemeralHLTPhysics2":[ "crab_EphemeralHLTPhysics2"],
        "EphemeralHLTPhysics3":[ "crab_EphemeralHLTPhysics3"],
        "EphemeralHLTPhysics4":[ "crab_EphemeralHLTPhysics4"],
        "EphemeralHLTPhysics5":[ "crab_EphemeralHLTPhysics5"],
        "EphemeralHLTPhysics6":[ "crab_EphemeralHLTPhysics6"],
        "EphemeralHLTPhysics7":[ "crab_EphemeralHLTPhysics7"],
        "EphemeralHLTPhysics8":[ "crab_EphemeralHLTPhysics8"]
    },
    "VBF": {
        "VBFHToTauTau_M125_TuneCUETP8M1_14TeV_powheg_pythia8":["crab_VBFHToTauTau"],
    }
}


ROOT.gInterpreter.Declare(
    """
    using Vec_t = const ROOT::VecOps::RVec<float>;
    float l1pt_threshold = 32. ;
    ROOT::VecOps::RVec<int> FindL1Taus(Vec_t& l1Tau_pt,Vec_t& l1Tau_hwIso) {
      ROOT::VecOps::RVec<int> n;
      for(auto k = 0; k < l1Tau_pt.size(); k++){
          if(l1Tau_pt.at(k) >= l1pt_threshold  && (l1Tau_hwIso.at(k) > 0 || l1Tau_pt.at(k) >= 70))
              n.push_back(k);
      }

      return n;
    }
    """
)
parser = argparse.ArgumentParser()
parser.add_argument('--sample', required=True, type=str, help= "input reference file", choices=['Data', 'VBF'])
parser.add_argument('--range', required=False, type=int, default=0 , help='number of event to process')
parser.add_argument('--n_max_files', required=False, type=int, default=100000, help='max number of files to be processed')
parser.add_argument('--last_module_name', required=False, type=str, default="", help='max number of files to be processed')
parser.add_argument('-p', required=False, type=bool, default=True, help=' True = process , False = just print')
parser.add_argument('--eff',required=False, type=bool, default=False, help=' True = evaluate efficiency/rate ')
parser.add_argument('--eff_algo', required=False, type=bool, default=False, help=' True = evaluate algo efficiency ')
parser.add_argument('--verbose', required=False, type=bool, default=0)
args = parser.parse_args()

def search_files(inputs):
    n_files = 1
    for k in files_directories[args.sample]:
            parent = parent_directory
            if(k == "ZprimeToTauTau_M-4000_TuneCP5_14TeV-pythia8-tauola" or k == "EphemeralHLTPhysics3" ):
                parent = parent_directory_2
            for i in files_directories[args.sample][k]:
                path = parent + k + "/"+ i + "/"
                #print path
                for file in os.listdir(path):
                    subdir_1 = os.listdir(path+file)
                    for l in subdir_1:
                        sub_dir2 = os.listdir(path+file+"/"+l)
                        complete_path = path+file+"/"+l+"/*.root"#+j
                        inputs.append(complete_path)
                        print complete_path
                        '''
                        #for j in sub_dir2:
                        #    complete_path = path+file+"/"+l+"/#+j
                            inputs.append(complete_path)
                            #print complete_path
                            if(args.n_max_files>0 and n_files == args.n_max_files):
                                return inputs
                            n_files += 1
                            '''

    return inputs
inputs = []
search_files(inputs)

n_zerotaus=0
n_onetau=0
all_events = 0
for i in range(0, len(inputs)):
    df = ROOT.RDataFrame('taus',inputs[i])
    pandataframe = LoadIdFrames(inputs[i].replace("*.root", "output_1.root"), ["module"])
    bigOrIndex = GetDataIdHash(pandataframe["module"], "hltL1sDoubleTauBigOR")
    df = df.Filter(("defaultDiTauPath_lastModuleIndex >={}").format(bigOrIndex["hltL1sDoubleTauBigOR"]+1)).Define("n_l1Taus", "FindL1Taus(l1Tau_pt,l1Tau_hwIso)")
    n_zerotaus+=df.Filter("n_l1Taus.size()==0").Count().GetValue()
    n_onetau+=df.Filter("n_l1Taus.size()<2").Count().GetValue()
print n_zerotaus
print n_onetau
