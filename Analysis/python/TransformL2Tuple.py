import os
import argparse

parent_directory =  "/gpfs/ddn/srm/cms/store/user/vdamante/L2TauTuple/"
parent_directory_2 = "/gpfs/ddn/srm/cms/store/user/vdamante/L2TausTuple/" # EphemeralHLTPhysics3 and ZprimeToTauTau
files_directories = {
"EphemeralHLTPhysics1":[ "crab_EphemeralHLTPhysics1"],
"EphemeralHLTPhysics2":[ "crab_EphemeralHLTPhysics2"],
"EphemeralHLTPhysics3":[ "crab_EphemeralHLTPhysics3"],
"EphemeralHLTPhysics4":[ "crab_EphemeralHLTPhysics4"],
"EphemeralHLTPhysics5":[ "crab_EphemeralHLTPhysics5"],
"EphemeralHLTPhysics6":[ "crab_EphemeralHLTPhysics6"],
"EphemeralHLTPhysics7":[ "crab_EphemeralHLTPhysics7"],
"EphemeralHLTPhysics8":[ "crab_EphemeralHLTPhysics8"],
"DYToLL_M-50_TuneCP5_14TeV-pythia8":["crab_DYToLL_M-50"],
"QCD_Pt_120to170_TuneCP5_14TeV_pythia8":["crab_QCD_Pt_120to170"],
"QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8":["crab_QCD_Pt-15to3000_Flat"],
"QCD_Pt-15to7000_TuneCP5_Flat_14TeV_pythia8":["crab_QCD_Pt-15to7000_Flat","crab_QCD_Pt-15to7000_Flat_ext1"],
"QCD_Pt_170to300_TuneCP5_14TeV_pythia8":["crab_QCD_Pt_170to300"],
"QCD_Pt_300to470_TuneCP5_14TeV_pythia8":["crab_QCD_Pt_300to470"],
"QCD_Pt_30to50_TuneCP5_14TeV_pythia8":["crab_QCD_Pt_30to50"],
"QCD_Pt_470to600_TuneCP5_14TeV_pythia8":["crab_QCD_Pt_470to600"],
"QCD_Pt_50to80_TuneCP5_14TeV_pythia8":["crab_QCD_Pt_50to80"],
"QCD_Pt_600oInf_TuneCP5_14TeV_pythia8":["crab_QCD_Pt_600oInf"],
"QCD_Pt_80to120_TuneCP5_14TeV_pythia8":["crab_QCD_Pt_80to120"],
"TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8":["crab_TTToSemiLeptonic"],
"TT_TuneCP5_14TeV-powheg-pythia8":["crab_TT","crab_TT_ext1"],
"VBFHToTauTau_M125_TuneCUETP8M1_14TeV_powheg_pythia8":["crab_VBFHToTauTau"],
"WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8":["crab_WJetsToLNu"],
"ZprimeToTauTau_M-4000_TuneCP5_14TeV-pythia8-tauola":["crab_ZprimeToTauTau_M-4000"],
}



parser = argparse.ArgumentParser()
#parser.add_argument('machine', type=str, default="local", choices=["lxplus", "local"])
parser.add_argument('--sample', required=True, type=str, help= "input reference file", choices=['DY', 'QCD_Pt-15to3000', 'QCD_Pt-15to7000', 'QCD_Pt_170to300', 'QCD_Pt_120to170', 'QCD_Pt_300to470', 'QCD_Pt_30to50', 'QCD_Pt_470to600', 'QCD_Pt_50to80', 'QCD_Pt_600oInf', 'QCD_Pt_80to120', 'TT', 'VBF', 'WJets', 'Zprime', 'EphemeralHLTPhysics1', 'EphemeralHLTPhysics2', 'EphemeralHLTPhysics3', 'EphemeralHLTPhysics4', 'EphemeralHLTPhysics5', 'EphemeralHLTPhysics6', 'EphemeralHLTPhysics7', 'EphemeralHLTPhysics8'])
parser.add_argument('--n_threads', required=False, type=int, default=1, help='Number of threads')
parser.add_argument('--start_entry', required=False, type=int, default=0 , help='start entry')
parser.add_argument('--end_entry', required=False, type=str, default="std::numeric_limits<Long64_t>::max()", help='end entry')
parser.add_argument('--n_max_files', required=False, type=int, default=100, help='max number of files to be processed')
parser.add_argument('-p', required=False, type=bool, help=' True = process , False = just print')

args = parser.parse_args()


cmd_in = 'L2TrainingTupleProducer '
inputs =[]
directory = "/home/users/damante/L2SkimmedTuples/"
def search_files(inputs):
    for k in files_directories:
        if(args.sample in k):
            parent = parent_directory
            if(k == "ZprimeToTauTau_M-4000_TuneCP5_14TeV-pythia8-tauola" or k == "EphemeralHLTPhysics3" ):
                parent = parent_directory_2
            for i in files_directories[k]:
                path = parent + k + "/"+ i + "/"
                #print "subdir = "+ path
                for file in os.listdir(path):
                    subdir_1 = os.listdir(path+file)
                    for i in subdir_1:
                        sub_dir2 = os.listdir(path+file+"/"+i)
                        for j in sub_dir2:
                            inputs.append(path+file+"/"+i+"/"+j)
    return inputs


search_files(inputs)
print directory
#"isQCDDataVBF", "0 = QCD; 1 = TT,DY,ZPrime; 2 = VBF; 3=Data"

isQCD_dict = {'DY':'1',
                'QCD_Pt-15to3000':'0',
                'QCD_Pt-15to7000':'0',
                'QCD_Pt_170to300':'0',
                'QCD_Pt_120to170':'0',
                'QCD_Pt_300to470':'0',
                'QCD_Pt_30to50':'0',
                'QCD_Pt_470to600':'0',
                'QCD_Pt_50to80':'0',
                'QCD_Pt_600oInf':'0',
                'QCD_Pt_80to120':'0',
                'TT':'1',
                'VBF':'2',
                'WJets':'1',
                'Zprime':'2',
                'EphemeralHLTPhysics1':'3',
                'EphemeralHLTPhysics2':'3',
                'EphemeralHLTPhysics3':'3',
                'EphemeralHLTPhysics4':'3',
                'EphemeralHLTPhysics5':'3',
                'EphemeralHLTPhysics6':'3',
                'EphemeralHLTPhysics7':'3',
                'EphemeralHLTPhysics8':'3'}
isQCD = isQCD_dict[args.sample]

#print directory + "\n"
k=1
for i in inputs:
    if k > args.n_max_files:
        pass
    else:
        cmd = cmd_in + ' --input \"' + i + '\"  --output \"' + directory
        filename  = ""
        for j in i.split("/"):
            if("crab_" in j):
                cmd =cmd + j.replace("crab_", "") + '/'
                if not os.path.exists(directory+j.replace("crab_", "")):
                    os.makedirs(directory+j.replace("crab_", ""))
            elif(args.sample in j):
                filename = j
            elif("output" in j):
                #print j
                cmd = cmd + filename +"_"+j + '\"'
        cmd = cmd + ' --isQCD ' + isQCD
        if(args.n_threads>1):
            cmd = cmd +  ' --n-threads ' + str(args.n_threads)
        if(args.start_entry>0):
            cmd = cmd +  ' --start-entry ' + str(args.start_entry)
        if(args.end_entry!="std::numeric_limits<Long64_t>::max()"):
            cmd = cmd +  ' --end-entry ' + args.end_entry
        print cmd+ "\n"
        if(args.p):
            print args.p
            os.system(cmd)
        k+=1
