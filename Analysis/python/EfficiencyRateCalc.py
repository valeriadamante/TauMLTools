import ROOT
import os
import argparse
import numpy as np
from array import array

dictionary_Pt1_filters = {
"df_num_30_20" : "genLepton_vis_pt[n_tauHad[0]] > 0 && genLepton_vis_pt[n_tauHad[0]] <= 30 ",
"df_num_40_30" : "genLepton_vis_pt[n_tauHad[0]] > 30 && genLepton_vis_pt[n_tauHad[0]] <= 40 ",
"df_num_50_40" : "genLepton_vis_pt[n_tauHad[0]] > 40 && genLepton_vis_pt[n_tauHad[0]] <= 50 ",
"df_num_60_50" : "genLepton_vis_pt[n_tauHad[0]] > 50 && genLepton_vis_pt[n_tauHad[0]] <= 60 ",
"df_num_70_60" : "genLepton_vis_pt[n_tauHad[0]] > 60 && genLepton_vis_pt[n_tauHad[0]] <= 70 ",
"df_num_80_70" : "genLepton_vis_pt[n_tauHad[0]] > 70 && genLepton_vis_pt[n_tauHad[0]] <= 80 ",
"df_num_90_80" : "genLepton_vis_pt[n_tauHad[0]] > 80 && genLepton_vis_pt[n_tauHad[0]] <= 90 ",
"df_num_100_90" : "genLepton_vis_pt[n_tauHad[0]] > 90 && genLepton_vis_pt[n_tauHad[0]] <= 100 ",
"df_num_120_100" : "genLepton_vis_pt[n_tauHad[0]] > 100 && genLepton_vis_pt[n_tauHad[0]] <= 120 ",
"df_num_140_120" : "genLepton_vis_pt[n_tauHad[0]] > 120 && genLepton_vis_pt[n_tauHad[0]] <= 140 ",
"df_num_160_100" : "genLepton_vis_pt[n_tauHad[0]] > 140 && genLepton_vis_pt[n_tauHad[0]] <= 160 ",
"df_num_180_100" : "genLepton_vis_pt[n_tauHad[0]] > 160 && genLepton_vis_pt[n_tauHad[0]] <= 180 ",
"df_num_200_100" : "genLepton_vis_pt[n_tauHad[0]] > 180 && genLepton_vis_pt[n_tauHad[0]] <= 200 ",
"df_num_250_100" : "genLepton_vis_pt[n_tauHad[0]] > 200 && genLepton_vis_pt[n_tauHad[0]] <= 250 ",
"df_num_300_100" : "genLepton_vis_pt[n_tauHad[0]] > 250 && genLepton_vis_pt[n_tauHad[0]] <= 300 ",
"df_num_350_100" : "genLepton_vis_pt[n_tauHad[0]] > 300 && genLepton_vis_pt[n_tauHad[0]] <= 350 "
}
dictionary_Pt2_filters = {
"df_num_30_20" : "genLepton_vis_pt[n_tauHad[1]] > 0 && genLepton_vis_pt[n_tauHad[1]] <= 30 ",
"df_num_40_30" : "genLepton_vis_pt[n_tauHad[1]] > 30 && genLepton_vis_pt[n_tauHad[1]] <= 40 ",
"df_num_50_40" : "genLepton_vis_pt[n_tauHad[1]] > 40 && genLepton_vis_pt[n_tauHad[1]] <= 50 ",
"df_num_60_50" : "genLepton_vis_pt[n_tauHad[1]] > 50 && genLepton_vis_pt[n_tauHad[1]] <= 60 ",
"df_num_70_60" : "genLepton_vis_pt[n_tauHad[1]] > 60 && genLepton_vis_pt[n_tauHad[1]] <= 70 ",
"df_num_80_70" : "genLepton_vis_pt[n_tauHad[1]] > 70 && genLepton_vis_pt[n_tauHad[1]] <= 80 ",
"df_num_90_80" : "genLepton_vis_pt[n_tauHad[1]] > 80 && genLepton_vis_pt[n_tauHad[1]] <= 90 ",
"df_num_100_90" : "genLepton_vis_pt[n_tauHad[1]] > 90 && genLepton_vis_pt[n_tauHad[1]] <= 100 ",
"df_num_120_100" : "genLepton_vis_pt[n_tauHad[1]] > 100 && genLepton_vis_pt[n_tauHad[1]] <= 120 ",
"df_num_140_120" : "genLepton_vis_pt[n_tauHad[1]] > 120 && genLepton_vis_pt[n_tauHad[1]] <= 140 ",
"df_num_160_100" : "genLepton_vis_pt[n_tauHad[1]] > 140 && genLepton_vis_pt[n_tauHad[1]] <= 160 ",
"df_num_180_100" : "genLepton_vis_pt[n_tauHad[1]] > 160 && genLepton_vis_pt[n_tauHad[1]] <= 180 ",
"df_num_200_100" : "genLepton_vis_pt[n_tauHad[1]] > 180 && genLepton_vis_pt[n_tauHad[1]] <= 200 ",
"df_num_250_100" : "genLepton_vis_pt[n_tauHad[1]] > 200 && genLepton_vis_pt[n_tauHad[1]] <= 250 ",
"df_num_300_100" : "genLepton_vis_pt[n_tauHad[1]] > 250 && genLepton_vis_pt[n_tauHad[1]] <= 300 ",
"df_num_350_100" : "genLepton_vis_pt[n_tauHad[1]] > 300 && genLepton_vis_pt[n_tauHad[1]] <= 350 "
}
#df_num_40_30 = df_num.Filter(true_tau1_pt > 40 && true_tau1_pt <= 50 && true_tau2_pt > 30 && true_tau2_pt <= 40)
#df_den_40_30 = df_den.Filter(....)

parent_directory =  "/gpfs/ddn/srm/cms/store/user/vdamante/L2TauTuple/"
parent_directory_2 = "/gpfs/ddn/srm/cms/store/user/vdamante/L2TausTuple/" # EphemeralHLTPhysics3 and ZprimeToTauTau
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
class EffRate:
    opt_stat = 1111
    files_open = []
    numerator = 0
    denominator = 0
    Pt1Bins= [0, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300, 350]
    Pt2Bins= [0, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300, 350]
    all_numerators= np.zeros((len(Pt1Bins)-1)*(len(Pt2Bins)-1))
    all_denominators= np.zeros((len(Pt1Bins)-1)*(len(Pt2Bins)-1))
    all_efficiencies= np.zeros((len(Pt1Bins)-1)*(len(Pt2Bins)-1))

    def __init__(self, input_files, tree, range):
        self.input_files = input_files
        self.tree = tree
        self.range = range

    def rate_calculator(self):
        for file in self.input_files:
            print "\nfile = " , file,
            df = ROOT.RDataFrame(self.tree, file)
            if(self.range>0):
                df = df.Range(int(self.range))
            df_num = df.Filter("defaultDiTauPath_lastModuleIndex >=39")
            numerator = df_num.Count().GetValue()
            self.numerator+=(numerator)
            print "partial numerator = ", numerator
            df_den = df.Filter("defaultDiTauPath_lastModuleIndex < 39")
            denominator = df_den.Count().GetValue()
            self.denominator+=(denominator)
            print "partial denominator = ", denominator
            print "partial rate = ", float(numerator)/float(denominator)

    def efficiency_calculator(self):
        ROOT.gInterpreter.Declare(
            """
            using Vec_t = const ROOT::VecOps::RVec<float>;
            float l1pt_threshold= 20.0;
            float genlep_eta_threshold = 2.1 ;
            float genlep_VisPt_threshold = 15.0 ;
            float delta_R_threshold = 0.5;
            ROOT::VecOps::RVec<int> FindTwoHadTaus(Vec_t& kind,Vec_t& vis_pt,Vec_t& vis_eta, float pt_threshold) {
              ROOT::VecOps::RVec<int> n;
              for(auto i = 0; i < kind.size(); i++){
                  if(kind.at(i) == 5 && vis_pt.at(i)>genlep_VisPt_threshold && std::abs(vis_eta.at(i))<genlep_eta_threshold)
                      n.push_back(i);
              }
              return n;
            }
            """
        )
        last_file_index = 0

        for file in self.input_files:
          ROOT.gStyle.SetOptStat(0)
          print "file = " , file, "\n"
          pt_threshold = 15.0
          df = ROOT.RDataFrame(self.tree, file)
          if(self.range>0):
              df = df.Range(int(self.range))
          df = df.Define("n_tauHad", "FindTwoHadTaus(genLepton_kind,genLepton_vis_pt,genLepton_vis_eta, {})".format(pt_threshold))
          df = df.Filter("n_tauHad.size() == 2")
          denominator = df.Count().GetValue()
          df_num = df.Filter("defaultDiTauPath_lastModuleIndex >=40")
          #df_num = df_num.Filter("for(auto i=0; i<n_tauHad.size(); i++) {return genLepton_vis_pt[n_tauHad.at(i)]>50;};")
          k =0
          for j in range(0, len(self.Pt1Bins)-1):
              for i in range(0, len(self.Pt2Bins)-1):
                  #print "genLepton_vis_pt[n_tauHad[0]] > {} && genLepton_vis_pt[n_tauHad[0]] <= {} ".format(self.Pt1Bins[j], self.Pt1Bins[j+1])
                  df_denominator_ptbin = df.Filter("genLepton_vis_pt[n_tauHad[0]] > {} && genLepton_vis_pt[n_tauHad[0]] <= {} ".format(self.Pt1Bins[j], self.Pt1Bins[j+1]))
                  df_ptbin = df_num.Filter("genLepton_vis_pt[n_tauHad[0]] > {} && genLepton_vis_pt[n_tauHad[0]] <= {} ".format(self.Pt1Bins[j], self.Pt1Bins[j+1]))

                  #print "genLepton_vis_pt[n_tauHad[1]] > {} && genLepton_vis_pt[n_tauHad[1]] <= {} ".format(self.Pt1Bins[i], self.Pt1Bins[i+1])
                  df_denominator_ptbin = df_denominator_ptbin.Filter("genLepton_vis_pt[n_tauHad[0]] > {} && genLepton_vis_pt[n_tauHad[0]] <= {} ".format(self.Pt1Bins[j], self.Pt1Bins[j+1]))
                  df_ptbin = df_ptbin.Filter("genLepton_vis_pt[n_tauHad[1]] > {} && genLepton_vis_pt[n_tauHad[1]] <= {} ".format(self.Pt1Bins[i], self.Pt1Bins[i+1]))

                  #print float(df_ptbin.Count().GetValue())/float(denominator)
                  self.all_numerators[k]+=(float(df_ptbin.Count().GetValue()))
                  self.all_denominators[k]+=float(df_denominator_ptbin.Count().GetValue())
                  k+=1

          numerator = df_num.Count().GetValue()
          self.numerator+=(numerator)
          print "partial numerator = ", numerator
          print "partial denominator = ", denominator
          self.denominator+=(denominator)
          efficiency = float(numerator)/float(denominator)
          print "partial efficiency = ", efficiency
          last_file_index += 1
          if(last_file_index == len(self.input_files)):
              for i in range(0, len(self.all_efficiencies)):
                  self.all_efficiencies[i]=float(self.all_numerators[i])/float(self.all_denominators[i])
              self.fill_histogram(self.all_efficiencies)



    def fill_histogram(self, efficiencies):
        #uniformTau = ROOT.TGraph(len(self.Pt1Bins)-1,array('d',self.Pt1Bins),array('d',efficiencies))
        canvas=ROOT.TCanvas()
        canvas.cd()
        uniformTau = ROOT.TH2D("efficiency","",len(self.Pt1Bins)-1,array('d',self.Pt1Bins),len(self.Pt2Bins)-1,array('d',self.Pt2Bins))
        uniformTau.GetXaxis().SetTitle("#tau 1 P_{T_} (GeV)")
        uniformTau.GetYaxis().SetTitle("#tau 2 P_{T_} (GeV)")
        for i in range(0,len(efficiencies)):
            uniformTau.SetBinContent(i+1, efficiencies[i])
        uniformTau.Draw("COLZ")
        canvas.Update()
        canvas.Print("efficiency.png", "png")
        raw_input()




parser = argparse.ArgumentParser()
parser.add_argument('--sample', required=True, type=str, help= "input reference file", choices=['Data', 'VBF'])
parser.add_argument('--range', required=False, type=int, default=0 , help='number of event to process')
parser.add_argument('--n_max_files', required=False, type=int, default=100, help='max number of files to be processed')
parser.add_argument('-p', required=False, type=bool, help=' True = process , False = just print')
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
                        for j in sub_dir2:
                            complete_path = path+file+"/"+l+"/"+j
                            inputs.append(complete_path)
                            #print complete_path
                            if(args.n_max_files>0 and n_files == args.n_max_files):
                                return inputs
                            n_files += 1

    return inputs
inputs = []
search_files(inputs)
ER_Calc= EffRate(inputs, "taus", args.range)
if(args.sample == "VBF"):
    ER_Calc.efficiency_calculator()
    print "\n total numerator = ", ER_Calc.numerator
    print "total denominator = ", ER_Calc.denominator
    print "total efficiency = ", float(ER_Calc.numerator)/float(ER_Calc.denominator)

elif(args.sample == "Data"):
    ER_Calc.rate_calculator()
    print "\ntotal numerator = ", ER_Calc.numerator
    print "total denominator = ", ER_Calc.denominator
    print "total rate = ", float(ER_Calc.numerator)/float(ER_Calc.denominator)
