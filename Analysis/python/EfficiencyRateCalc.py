import ROOT
import os
import argparse
import numpy as np
from array import array
from TauMLTools.Analysis.getPathNames import *

dictionary_Pt1_filters = {
"df_num_25_20" : "genLepton_vis_pt[n_tauHad[0]] > 0 && genLepton_vis_pt[n_tauHad[0]] <= 30 ",
#"df_num_30_25" : "genLepton_vis_pt[n_tauHad[0]] > 25 && genLepton_vis_pt[n_tauHad[0]] <= 30 ",
"df_num_35_30" : "genLepton_vis_pt[n_tauHad[0]] > 30 && genLepton_vis_pt[n_tauHad[0]] <= 35 ",
"df_num_40_35" : "genLepton_vis_pt[n_tauHad[0]] > 35 && genLepton_vis_pt[n_tauHad[0]] <= 40 ",
"df_num_50_40" : "genLepton_vis_pt[n_tauHad[0]] > 40 && genLepton_vis_pt[n_tauHad[0]] <= 50 ",
"df_num_60_50" : "genLepton_vis_pt[n_tauHad[0]] > 50 && genLepton_vis_pt[n_tauHad[0]] <= 60 ",
"df_num_70_60" : "genLepton_vis_pt[n_tauHad[0]] > 60 && genLepton_vis_pt[n_tauHad[0]] <= 70 ",
"df_num_80_70" : "genLepton_vis_pt[n_tauHad[0]] > 70 && genLepton_vis_pt[n_tauHad[0]] <= 80 ",
"df_num_90_80" : "genLepton_vis_pt[n_tauHad[0]] > 80 && genLepton_vis_pt[n_tauHad[0]] <= 90 ",
"df_num_90_80" : "genLepton_vis_pt[n_tauHad[1]] > 80 && genLepton_vis_pt[n_tauHad[1]] <= 90 ",
"df_num_130_90" : "genLepton_vis_pt[n_tauHad[0]] > 90 && genLepton_vis_pt[n_tauHad[0]] <= 130 ",
#"df_num_170_130" : "genLepton_vis_pt[n_tauHad[0]] > 130 && genLepton_vis_pt[n_tauHad[0]] <= 170 ",
#"df_num_210_170" : "genLepton_vis_pt[n_tauHad[0]] > 170 && genLepton_vis_pt[n_tauHad[0]] <= 210 ",
#"df_num_250_210" : "genLepton_vis_pt[n_tauHad[0]] > 210 && genLepton_vis_pt[n_tauHad[0]] <= 250 ",
"df_num_350_170" : "genLepton_vis_pt[n_tauHad[0]] > 130 && genLepton_vis_pt[n_tauHad[0]] <= 250 ",
#"df_num_350_300" : "genLepton_vis_pt[n_tauHad[0]] > 300 && genLepton_vis_pt[n_tauHad[0]] <= 350 "
}
dictionary_Pt2_filters = {
"df_num_25_20" : "genLepton_vis_pt[n_tauHad[1]] > 0 && genLepton_vis_pt[n_tauHad[1]] <= 30 ",
#"df_num_30_25" : "genLepton_vis_pt[n_tauHad[1]] > 25 && genLepton_vis_pt[n_tauHad[1]] <= 30 ",
"df_num_35_30" : "genLepton_vis_pt[n_tauHad[1]] > 30 && genLepton_vis_pt[n_tauHad[1]] <= 35 ",
"df_num_40_35" : "genLepton_vis_pt[n_tauHad[1]] > 35 && genLepton_vis_pt[n_tauHad[1]] <= 40 ",
"df_num_40_35" : "genLepton_vis_pt[n_tauHad[1]] > 40 && genLepton_vis_pt[n_tauHad[1]] <= 45 ",
"df_num_50_40" : "genLepton_vis_pt[n_tauHad[1]] > 45 && genLepton_vis_pt[n_tauHad[1]] <= 50 ",
"df_num_60_50" : "genLepton_vis_pt[n_tauHad[1]] > 50 && genLepton_vis_pt[n_tauHad[1]] <= 60 ",
"df_num_70_60" : "genLepton_vis_pt[n_tauHad[1]] > 60 && genLepton_vis_pt[n_tauHad[1]] <= 70 ",
"df_num_80_70" : "genLepton_vis_pt[n_tauHad[1]] > 70 && genLepton_vis_pt[n_tauHad[1]] <= 80 ",
"df_num_90_80" : "genLepton_vis_pt[n_tauHad[1]] > 80 && genLepton_vis_pt[n_tauHad[1]] <= 90 ",
"df_num_130_90" : "genLepton_vis_pt[n_tauHad[1]] > 90 && genLepton_vis_pt[n_tauHad[1]] <= 350 ",
#"df_num_170_130" : "genLepton_vis_pt[n_tauHad[1]] > 130 && genLepton_vis_pt[n_tauHad[1]] <= 170 ",
#"df_num_210_170" : "genLepton_vis_pt[n_tauHad[1]] > 170 && genLepton_vis_pt[n_tauHad[1]] <= 210 ",
#"df_num_250_210" : "genLepton_vis_pt[n_tauHad[1]] > 210 && genLepton_vis_pt[n_tauHad[1]] <= 250 ",
#"df_num_350_170" : "genLepton_vis_pt[n_tauHad[1]] > 130 && genLepton_vis_pt[n_tauHad[1]] <= 250 ",
#"df_num_350_300" : "genLepton_vis_pt[n_tauHad[1]] > 300 && genLepton_vis_pt[n_tauHad[1]] <= 350 "
}

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
    float genlep_eta_threshold = 2.1 ;
    ROOT::VecOps::RVec<int> FindTwoHadTaus(ROOT::VecOps::RVec<int>& kind,Vec_t& vis_pt,Vec_t& vis_eta, float pt_threshold) {
      ROOT::VecOps::RVec<int> n;
      for(auto i = 0; i < kind.size(); i++){
          if(kind.at(i) == 5 && vis_pt.at(i)>=pt_threshold && std::abs(vis_eta.at(i))<=genlep_eta_threshold)
              n.push_back(i);
      }
      return n;
    }
    """
)

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
ROOT.gInterpreter.Declare(
    """
    ROOT::VecOps::RVec<float> ReorderHadTaus(ROOT::VecOps::RVec<float>& vis_pt, ROOT::VecOps::RVec<float>& branch_to_order, ROOT::VecOps::RVec<int> NTauHad ) {

        ROOT::VecOps::RVec<float> not_ordered_pt;
        ROOT::VecOps::RVec<float> not_ordered_branch;
        ROOT::VecOps::RVec<float> reordered_branch;
        ROOT::VecOps::RVec<int> already_taken_index;

        int k=0;
        int pt_max_index =0;
        for(auto& i : NTauHad){
            not_ordered_branch.push_back(branch_to_order.at(i));
            not_ordered_pt.push_back(vis_pt.at(i));
        }

        int not_ordered_branch_size = not_ordered_branch.size();

        while(k<not_ordered_branch_size){
              float pt_max = 0.;
              for(int i = 0; i < not_ordered_branch.size(); i++){
                  if (std::find(already_taken_index.begin(),already_taken_index.end(), i) != already_taken_index.end()) {
                      continue;
                  }
                  if(not_ordered_pt.at(i)>=pt_max){
                      pt_max_index = i;
                      pt_max = not_ordered_pt.at(i);
                  }
              }
              reordered_branch.push_back(not_ordered_branch.at(pt_max_index));
              already_taken_index.push_back(pt_max_index);
              k++;
        }
        return reordered_branch;
    }
    """

)

ROOT.gInterpreter.Declare(
    """
    auto pi = TMath::Pi();

    float DeltaPhi(Float_t phi1, Float_t phi2)
    {
        float dphi = phi1 - phi2;
        if(dphi > pi)
            dphi -= 2*pi;
        else if(dphi <= -pi)
            dphi += 2*pi;
        return dphi;
    }
    float deltaR(Float_t phi1,Float_t eta1,Float_t phi2,Float_t eta2) {
        auto dphi = DeltaPhi(phi1, phi2);
        auto deta = eta1-eta2;
        return (std::sqrt(deta * deta + dphi * dphi));

    }

    ROOT::VecOps::RVec<float> FindCorrespondenceL1GenTau(ROOT::VecOps::RVec<float>& branch_to_find_corr, ROOT::VecOps::RVec<float>& gen_tau_vis_pt, ROOT::VecOps::RVec<float>& gen_tau_vis_eta, ROOT::VecOps::RVec<float>& gen_tau_vis_phi, ROOT::VecOps::RVec<float>& l1Tau_pt, ROOT::VecOps::RVec<float>& l1Tau_eta, ROOT::VecOps::RVec<float>& l1Tau_phi, ROOT::VecOps::RVec<int>& l1Tau_hwIso) {

        ROOT::VecOps::RVec<float> branch_with_correspondence ;
        ROOT::VecOps::RVec<int> taus_indices;

        for (int i = 0 ; i < gen_tau_vis_pt.size(); i++){
            int l1tau_previous_index = -100;
            float previous_deltaR = 100.;

            for (int j = 0 ; j<l1Tau_pt.size(); j++){

                if (std::find(taus_indices.begin(),taus_indices.end(), j) != taus_indices.end()) {
                    continue;
                } // end of if find indices


                if (l1Tau_pt.at(j) >= 32.  && (l1Tau_hwIso.at(j) > 0 || l1Tau_pt.at(j) >= 70)) {

                    float current_deltaR = deltaR(l1Tau_phi.at(j), l1Tau_eta.at(j), gen_tau_vis_phi.at(i), gen_tau_vis_eta.at(i));

                    if (current_deltaR < 0.3 && current_deltaR < previous_deltaR ){
                        previous_deltaR = current_deltaR;
                        l1tau_previous_index = j;

                    } // end of if (dR)

                } // end of if on L1 taus

            } // end of loop over L1 taus


            if(l1tau_previous_index!=-100) taus_indices.push_back(l1tau_previous_index);
            if(l1tau_previous_index!=-100) branch_with_correspondence.push_back(branch_to_find_corr.at(i));

        } // end of loop over gentaus


        return branch_with_correspondence;
    } // end of function


    """
)

class EffRate:
    opt_stat = 0
    files_open = []
    numerator = 0
    denominator = 0
    Pt1Bins= [20,30,35,40,50,60,70,80,90,350]#, 600]
    Pt2Bins= [20,30,35,40,50,60,70,80,90,350]#, 600]
    all_numerators= np.zeros([(len(Pt1Bins)-1),(len(Pt2Bins)-1)])
    all_denominators= np.ones([(len(Pt1Bins)-1),(len(Pt2Bins)-1)])
    all_efficiencies= np.zeros([(len(Pt1Bins)-1),(len(Pt2Bins)-1)])
    average_taus=  []
    last_module_name = "hltDoubleL2IsoTau26eta2p2"
    pt_threshold = 15.0
    n_events_1l1tau =0

    def __init__(self, input_files, tree, range, verbose):
        self.input_files = input_files
        self.tree = tree
        self.range = range
        self.verbose = verbose

    def average_l1tau_calculator(self, last_module_name):
        self.average_taus= []
        print "ciao"
        self.n_events_1l1tau=0
        print self.verbose
        for file in self.input_files:
            df = ROOT.RDataFrame(self.tree, file)
            pandataframe = LoadIdFrames(file.replace("*.root", "output_1.root"), ["module"])
            last_module_index = GetDataIdHash(pandataframe["module"], last_module_name)
            bigOrIndex = GetDataIdHash(pandataframe["module"], "hltL1sDoubleTauBigOR")
            if(self.range>0):
                df = df.Range(int(self.range))
            df = df.Filter(("defaultDiTauPath_lastModuleIndex >={}").format(bigOrIndex["hltL1sDoubleTauBigOR"]+1)).Define("n_l1Taus", "FindL1Taus(l1Tau_pt,l1Tau_hwIso)")
            self.average_taus.append(df.Define("number_L1Taus","n_l1Taus.size()").Mean("number_L1Taus").GetValue())
            if(self.verbose>1):
                print(df.Display("n_l1Taus", 100).Print())
                print(df.Filter("n_l1Taus.size()<2").Count().GetValue())
                print(df.Filter("n_l1Taus.size()<2").Count().GetValue())
            self.n_events_1l1tau+=df.Filter("n_l1Taus.size()<2").Count().GetValue()


    def rate_calculator(self, last_module_name):
        self.numerator=0.
        self.demominator=0.
        for file in self.input_files:
            df = ROOT.RDataFrame(self.tree, file)
            #print(file.replace("*.root", "output_1.root"))
            pandataframe = LoadIdFrames(file.replace("*.root", "output_1.root"), ["module"])
            last_module_index = GetDataIdHash(pandataframe["module"], last_module_name)
            if(self.range>0):
                df = df.Range(int(self.range))
                if(self.verbose>1):
                    print(df.Display("defaultDiTauPath_lastModuleIndex", 100).Print())
            df_num = df.Filter(("defaultDiTauPath_lastModuleIndex >={}").format(last_module_index[last_module_name]+1))
            numerator = df_num.Count().GetValue()
            self.numerator+=(numerator)
            #print "partial numerator = ", numerator
            df_den = df
            denominator = df_den.Count().GetValue()
            self.denominator+=(denominator)
            #print "partial denominator = ", denominator
            #print "partial rate = ", float(numerator)/float(denominator)

    def algo_rate_calculator(self, last_module_name):
        self.numerator=0.
        self.demominator=0.
        for file in self.input_files:
            df = ROOT.RDataFrame(self.tree, file)
            #print(file.replace("*.root", "output_1.root"))
            pandataframe = LoadIdFrames(file.replace("*.root", "output_1.root"), ["module"])
            last_module_index = GetDataIdHash(pandataframe["module"], last_module_name)
            bigOrIndex = GetDataIdHash(pandataframe["module"], "hltL1sDoubleTauBigOR")
            if(self.range>0):
                df = df.Range(int(self.range))
                if(self.verbose>1):
                    print(df.Display("defaultDiTauPath_lastModuleIndex", 100).Print())
            df_num = df.Filter(("defaultDiTauPath_lastModuleIndex >={}").format(last_module_index[last_module_name]+1))
            numerator = df_num.Count().GetValue()
            self.numerator+=(numerator)
            #print "partial numerator = ", numerator
            df_den = df.Filter(("defaultDiTauPath_lastModuleIndex >={}").format(bigOrIndex["hltL1sDoubleTauBigOR"]+1))
            denominator = df_den.Count().GetValue()
            self.denominator+=(denominator)
            #print "partial denominator = ", denominator
            #print "partial rate = ", float(numerator)/float(denominator)


    def full_path_rate_calculator(self):
        rate_calculator("hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg")

    def efficiency_calculator(self, last_module_name):
        if(self.verbose>0):
            print("calculating efficiency ... ")
        self.numerator=0.
        self.demominator=0.
        last_file_index = 0
        for file in self.input_files:
          if(self.verbose>0):
            print("analyzing file ... ", file)
          ROOT.gStyle.SetOptStat(self.opt_stat)
          df = ROOT.RDataFrame(self.tree, file)
          pandataframe = LoadIdFrames(file.replace("*.root", "output_1.root").replace("*.root", "output_1.root"), ["module"])
          last_module_index = GetDataIdHash(pandataframe["module"], last_module_name)
          #bigOrIndex = GetDataIdHash(pandataframe["module"], "hltL1sDoubleTauBigOR")
          #df_den = df.Filter(("defaultDiTauPath_lastModuleIndex >={}").format(bigOrIndex["hltL1sDoubleTauBigOR"]+1))
          if(self.verbose>1):
              print("load tree ... ", self.tree)
          if(self.range>0):
              df = df.Range(int(self.range))
              if(self.verbose>1):
                  print(df.Display("defaultDiTauPath_lastModuleIndex", 100).Print())

          print("inizialmente c'erano ... ", df.Count().GetValue())

          df = df.Define("n_tauHad", "FindTwoHadTaus(genLepton_kind,genLepton_vis_pt,genLepton_vis_eta, {})".format(self.pt_threshold))
          df = df.Filter("n_tauHad.size() == 2")
          print("DOPO due tau had ci sono ... ", df.Count().GetValue())
          print("saving denominator with two had taus ")



          df = df.Define("reordered_vis_pt", "ReorderHadTaus(genLepton_vis_pt,genLepton_vis_pt, n_tauHad)")
          df = df.Define("reordered_vis_phi", "ReorderHadTaus(genLepton_vis_pt,genLepton_vis_phi, n_tauHad)")
          df = df.Define("reordered_vis_eta", "ReorderHadTaus(genLepton_vis_pt,genLepton_vis_eta, n_tauHad)")

          df_den = df

          df = df.Define("n_l1Taus", "FindL1Taus(l1Tau_pt,l1Tau_hwIso)")
          df = df.Filter("n_l1Taus.size()>=1 ")
          print("DOPO almeno un l1tau presente ci sono ... ", df.Count().GetValue())
          df = df.Define("vis_pt_with_corr", "FindCorrespondenceL1GenTau(reordered_vis_pt, reordered_vis_pt, reordered_vis_eta, reordered_vis_phi, l1Tau_pt, l1Tau_eta, l1Tau_phi, l1Tau_hwIso)")

          df = df.Filter("vis_pt_with_corr.size()>=1 ")
          print("DOPO almeno una corrispondenza presente ci sono ... ", df.Count().GetValue())

          df = df.Filter("vis_pt_with_corr.size()==2 ")
          print("DOPO due corrispondenze presenti ci sono ... ", df.Count().GetValue())

          df = df.Filter(("defaultDiTauPath_lastModuleIndex >={}").format(last_module_index[last_module_name]+1))
          print("dopo last module index ci sono... ", df.Count().GetValue())
          print("saving numerator with 2 had taus with match and last module index ")
          df_num = df

          if(self.verbose>1):
              print(df.Display({"n_tauHad", "genLepton_vis_pt"}, 100).Print())
              print(df.Filter("reordered_vis_pt[1]>100").Display({"n_tauHad", "genLepton_vis_pt", "reordered_vis_pt"}).Print())
              print(df.Display("n_l1Taus", 100).Print())

          denominator = df_den.Count().GetValue()
          numerator = df_num.Count().GetValue()
          #print("numerator = ")
          if(self.verbose>1):
              print("computed  numerator & denominator... ", file)
          self.numerator+=(numerator)
          self.denominator+=(denominator)
          efficiency = float(numerator)/float(denominator)
          k = 0
          for pt1_index in range(0, len(self.Pt1Bins)-1):
               for pt2_index in range(0, len(self.Pt2Bins)-1):
                   if(pt1_index>=pt2_index):
                       print(("filling bin .. [ {}, {} ]").format(pt1_index, pt2_index))
                       # quello piu energetico
                       df_denominator_ptbin = df_den.Filter("reordered_vis_pt[0] > {} && reordered_vis_pt[0] <= {} ".format(self.Pt1Bins[pt1_index], self.Pt1Bins[pt1_index+1]))
                       df_numerator_ptbin = df_num.Filter("vis_pt_with_corr[0] > {} && vis_pt_with_corr[0] <= {} ".format(self.Pt1Bins[pt1_index], self.Pt1Bins[pt1_index+1]))
                       df_denominator_ptbin = df_denominator_ptbin.Filter("reordered_vis_pt[1] > {} && reordered_vis_pt[1] <= {} ".format(self.Pt2Bins[pt2_index], self.Pt2Bins[pt2_index+1]))
                       df_numerator_ptbin = df_numerator_ptbin.Filter("vis_pt_with_corr[1] > {} && vis_pt_with_corr[1] <= {} ".format(self.Pt2Bins[pt2_index], self.Pt2Bins[pt2_index+1]))
                       numerator = float(df_numerator_ptbin.Count().GetValue())
                       #print numerator
                       denominator = float(df_denominator_ptbin.Count().GetValue())
                       #print denominator
                       if(k==0):
                           self.all_numerators[pt1_index][pt2_index]=(numerator)
                           if(denominator!=0):
                                self.all_denominators[pt1_index][pt2_index]=float(denominator)
                       else:
                           self.all_numerators[pt1_index][pt2_index]+=(numerator)
                           self.all_denominators[pt1_index][pt2_index]+=float(denominator)
                   k+=1
          last_file_index += 1
          if(self.verbose>0):
              print("end of calculation, now filling histogram")
          if(last_file_index == len(self.input_files)):
              self.fill_histogram("normal", last_module_name, True)
              self.save_numbers("normal", last_module_name)

    def algo_efficiency_calculator(self, last_module_name):
            if(self.verbose>0):
                print("calculating efficiency ... ")
            self.numerator=0.
            self.demominator=0.
            last_file_index = 0
            for file in self.input_files:
              if(self.verbose>0):
                print("analyzing file ... ", file)
              ROOT.gStyle.SetOptStat(self.opt_stat)
              df = ROOT.RDataFrame(self.tree, file)
              pandataframe = LoadIdFrames(file.replace("*.root", "output_1.root").replace("*.root", "output_1.root"), ["module"])
              last_module_index = GetDataIdHash(pandataframe["module"], last_module_name)
              bigOrIndex = GetDataIdHash(pandataframe["module"], "hltL1sDoubleTauBigOR")

              if(self.verbose>1):
                  print("load tree ... ", self.tree)
              if(self.range>0):
                  df = df.Range(int(self.range))
                  if(self.verbose>1):
                      print(df.Display("defaultDiTauPath_lastModuleIndex", 100).Print())

              print("inizialmente c'erano ... ", df.Count().GetValue())

              df = df.Define("n_tauHad", "FindTwoHadTaus(genLepton_kind,genLepton_vis_pt,genLepton_vis_eta, {})".format(self.pt_threshold))
              df = df.Filter("n_tauHad.size() == 2")
              print("DOPO due tau had ci sono ... ", df.Count().GetValue())
              print("saving denominator with two had taus ")


              df = df.Define("n_l1Taus", "FindL1Taus(l1Tau_pt,l1Tau_hwIso)")
              df = df.Filter("n_l1Taus.size()>=1 ")
              print("DOPO almeno un l1tau presente ci sono ... ", df.Count().GetValue())

              df = df.Define("reordered_vis_pt", "ReorderHadTaus(genLepton_vis_pt,genLepton_vis_pt, n_tauHad)")
              df = df.Define("reordered_vis_phi", "ReorderHadTaus(genLepton_vis_pt,genLepton_vis_phi, n_tauHad)")
              df = df.Define("reordered_vis_eta", "ReorderHadTaus(genLepton_vis_pt,genLepton_vis_eta, n_tauHad)")

              df = df.Define("vis_pt_with_corr", "FindCorrespondenceL1GenTau(reordered_vis_pt, reordered_vis_pt, reordered_vis_eta, reordered_vis_phi, l1Tau_pt, l1Tau_eta, l1Tau_phi, l1Tau_hwIso)")

              df = df.Filter("vis_pt_with_corr.size()>=1 ")
              print("DOPO almeno una corrispondenza presente ci sono ... ", df.Count().GetValue())

              df = df.Filter("vis_pt_with_corr.size()==2 ")
              print("DOPO due corrispondenze presenti ci sono ... ", df.Count().GetValue())
              df_den = df.Filter(("defaultDiTauPath_lastModuleIndex >={}").format(bigOrIndex["hltL1sDoubleTauBigOR"]+1))
              print("dopo big or index ci sono... ", df_den.Count().GetValue())

              print("saving denominator with: 2 had taus, match with l1 taus, and last module index Big or")
              df_num = df.Filter(("defaultDiTauPath_lastModuleIndex >={}").format(last_module_index[last_module_name]+1))
              print("dopo last module index ci sono... ", df_num.Count().GetValue())


              if(self.verbose>1):
                  print(df.Display({"n_tauHad", "genLepton_vis_pt"}, 100).Print())
                  print(df.Filter("reordered_vis_pt[1]>100").Display({"n_tauHad", "genLepton_vis_pt", "reordered_vis_pt"}).Print())
                  print(df.Display("n_l1Taus", 100).Print())

              denominator = df_den.Count().GetValue()
              numerator = df_num.Count().GetValue()
              #print("numerator = ")
              if(self.verbose>1):
                  print("computed  numerator & denominator... ", file)
              self.numerator+=(numerator)
              self.denominator+=(denominator)
              efficiency = float(numerator)/float(denominator)
              k = 0
              for pt1_index in range(0, len(self.Pt1Bins)-1):
                   for pt2_index in range(0, len(self.Pt2Bins)-1):
                       if(pt1_index>=pt2_index):
                           # quello piu energetico
                           print(("filling bin .. [ {}, {} ]").format(pt1_index, pt2_index))
                           df_denominator_ptbin = df_den.Filter("vis_pt_with_corr[0] > {} && vis_pt_with_corr[0] <= {} ".format(self.Pt1Bins[pt1_index], self.Pt1Bins[pt1_index+1]))
                           df_numerator_ptbin = df_num.Filter("vis_pt_with_corr[0] > {} && vis_pt_with_corr[0] <= {} ".format(self.Pt1Bins[pt1_index], self.Pt1Bins[pt1_index+1]))
                           df_denominator_ptbin = df_denominator_ptbin.Filter("vis_pt_with_corr[1] > {} && vis_pt_with_corr[1] <= {} ".format(self.Pt2Bins[pt2_index], self.Pt2Bins[pt2_index+1]))
                           df_numerator_ptbin = df_numerator_ptbin.Filter("vis_pt_with_corr[1] > {} && vis_pt_with_corr[1] <= {} ".format(self.Pt2Bins[pt2_index], self.Pt2Bins[pt2_index+1]))
                           numerator = float(df_numerator_ptbin.Count().GetValue())
                           #print numerator
                           denominator = float(df_denominator_ptbin.Count().GetValue())
                           #print denominator
                           if(k==0):
                               self.all_numerators[pt1_index][pt2_index]=(numerator)
                               if(denominator!=0):
                                    self.all_denominators[pt1_index][pt2_index]=float(denominator)
                           else:
                               self.all_numerators[pt1_index][pt2_index]+=(numerator)
                               self.all_denominators[pt1_index][pt2_index]+=float(denominator)
                       k+=1
              last_file_index += 1
              if(self.verbose>0):
                  print("end of calculation, now filling histogram")
              if(last_file_index == len(self.input_files)):
                  self.fill_histogram("algo", last_module_name, True)
                  self.save_numbers("algo", last_module_name)




    def save_numbers(self, last_module_name, type):
        import json
        # save to npy file
        full_name = "efficiency_NumAndDen_"+type+"_"+last_module_name+".json"
        dict = {}
        for pt1_index in range(0, len(self.Pt1Bins)-1):
            dict[self.Pt1Bins[pt1_index]]={}
            for pt2_index in range(0, len(self.Pt2Bins)-1):
                dict[self.Pt1Bins[pt1_index]][self.Pt2Bins[pt2_index]]=[self.all_numerators[pt1_index][pt2_index], self.all_denominators[pt1_index][pt2_index]]
        jsonString = json.dumps(dict)
        jsonFile = open(full_name, "w")
        jsonFile.write(jsonString)
        jsonFile.close()


    def fill_histogram(self, type, last_module_name,add1Dplot=False):
        ROOT.gStyle.SetPaintTextFormat(".2f")
        if(self.verbose>0):
            print("opening canvas, histogram")
        #uniformTau = ROOT.TGraph(len(self.Pt1Bins)-1,array('d',self.Pt1Bins),array('d',efficiencies))
        canvas=ROOT.TCanvas()
        canvas.cd()
        uniformTau = ROOT.TH2D("efficiency_"+type,"efficiency_"+type,len(self.Pt1Bins)-1,array('d',self.Pt1Bins),len(self.Pt2Bins)-1,array('d',self.Pt2Bins))
        hPassed = ROOT.TH1D("passed_"+type,"passed_"+type,len(self.Pt1Bins)-1,array('d',self.Pt1Bins))
        hTotal = ROOT.TH1D("total_"+type,"total_"+type,len(self.Pt1Bins)-1,array('d',self.Pt1Bins))
        uniformTau.GetXaxis().SetTitle("#tau 1 P_{T_} (GeV)")
        uniformTau.GetYaxis().SetTitle("#tau 2 P_{T_} (GeV)")
        x, y = array( 'd' ), array( 'd' )
        for pt1_index in range(0, len(self.Pt1Bins)-1):
             for pt2_index in range(0, len(self.Pt2Bins)-1):
                 efficiency = float(self.all_numerators[pt1_index][pt2_index])/float(self.all_denominators[pt1_index][pt2_index])
                 if(self.verbose>1):
                     print "bin X [", self.Pt1Bins[pt1_index], self.Pt1Bins[pt1_index+1], "] bin Y [", self.Pt2Bins[pt2_index], self.Pt2Bins[pt2_index+1], "] Valore " , efficiency
                 uniformTau.SetBinContent((pt1_index+1),(pt2_index+1), efficiency)
                 if(pt1_index==pt2_index and add1Dplot):
                     x.append(self.Pt1Bins[pt1_index])
                     y.append(np.sqrt(efficiency))
                     if(self.verbose>0):
                         print "x = ", x[pt1_index], "\t num = ", np.sqrt(self.all_numerators[pt1_index][pt2_index]), "\t den = ",  np.sqrt(self.all_denominators[pt1_index][pt2_index])
                         print "bin X [", self.Pt1Bins[pt1_index], self.Pt1Bins[pt1_index+1], "] bin Y [", self.Pt2Bins[pt2_index], self.Pt2Bins[pt2_index+1], "] num= " , self.all_numerators[pt1_index][pt2_index], " den = ", self.all_denominators[pt1_index][pt2_index]
                     hPassed.SetBinContent((pt1_index+1), np.sqrt(self.all_numerators[pt1_index][pt2_index]))
                     hTotal.SetBinContent((pt1_index+1), np.sqrt(self.all_denominators[pt1_index][pt2_index]))
        n = len(x)
        uniformTau.Draw("TEXT2 COLZ")
        canvas.SetLogx()
        canvas.SetLogy()
        canvas.Update()
        full_name = "efficiency_"+type+"_"+last_module_name
        canvas.Print(full_name+".png", "png")
        #save histograms
        myfile = ROOT.TFile( full_name+".root", 'RECREATE' )
        uniformTau.Write()
        if(add1Dplot):
            canvas3=ROOT.TCanvas()
            canvas3.cd()
            EffGraph = ROOT.TEfficiency(hPassed,hTotal)
            hPassed.Write()
            hTotal.Write()
            EffGraph.Write()
            EffGraph.SetMarkerColor( 4 )
            EffGraph.SetMarkerStyle( 21 )
            EffGraph.SetTitle("Efficiency;#tau P_{T_} (GeV);#epsilon")
            EffGraph.SetMarkerColor(4)
            EffGraph.SetMarkerStyle(20)
            EffGraph.Draw("AP")
            #canvas3.SetLogx()
            canvas3.Update()
            canvas3.Print(full_name+"_Histo1D.png", "png")
        myfile.Close()
        #raw_input()


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
                        #complete_path = path+file+"/"+l+"/output_1.root"#+j
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

ER_Calc= EffRate(inputs, "taus", args.range, args.verbose)
if(args.last_module_name != ""):
    ER_Calc.last_module_name==args.last_module_name

if(args.p == True):
    if(args.sample == "VBF"):
        if(args.eff==True):
            '''
            ER_Calc.efficiency_calculator("hltL1sDoubleTauBigOR")
            print "\nafter hltL1sDoubleTauBigOR"
            print "total numerator = ", ER_Calc.numerator
            print "total denominator = ", ER_Calc.denominator
            print "total efficiency = ", float(ER_Calc.numerator)/float(ER_Calc.denominator)
            '''
            ER_Calc.efficiency_calculator("hltDoubleL2IsoTau26eta2p2")
            print "\nafter hltDoubleL2IsoTau26eta2p2"
            print "total numerator = ", ER_Calc.numerator
            print "total denominator = ", ER_Calc.denominator
            print "total efficiency = ", float(ER_Calc.numerator)/float(ER_Calc.denominator)

        if(args.eff_algo==True):
            ER_Calc.algo_efficiency_calculator("hltDoubleL2IsoTau26eta2p2")
            print "\n after hltDoubleL2IsoTau26eta2p2 - algo efficiency"
            print "total numerator = ", ER_Calc.numerator
            print "total denominator = ", ER_Calc.denominator
            print "total efficiency = ", float(ER_Calc.numerator)/float(ER_Calc.denominator)

    elif(args.sample == "Data"):
        if(args.eff==True):
            '''
            ER_Calc.rate_calculator("hltDoubleL2IsoTau26eta2p2")
            print "\nafter hltDoubleL2IsoTau26eta2p2"
            print "total numerator = ", ER_Calc.numerator
            print "total denominator = ", ER_Calc.denominator
            print "total rate = ", float(ER_Calc.numerator)*75817.94/float(ER_Calc.denominator)
            print "\nafter hltDoubleL2IsoTau26eta2p2"
            '''
            ER_Calc.average_l1tau_calculator("hltDoubleL2IsoTau26eta2p2")
            sum =0.
            for i in ER_Calc.average_taus:
                #print "average taus = ", i
                sum+=i
            mean=sum/float(len(ER_Calc.average_taus))
            print "total mean = ", mean
            print "events with 1 tau =", ER_Calc.n_events_1l1tau

        if(args.eff_algo==True):
            '''
            ER_Calc.full_path_rate_calculator()
            print "\n after hltBoolEnd"
            print "total numerator = ", ER_Calc.numerator
            print "total denominator = ", ER_Calc.denominator
            print "total rate = ", float(ER_Calc.numerator)*75817.94/float(ER_Calc.denominator)

            ER_Calc.rate_calculator("hltL1sDoubleTauBigOR")
            print "\nafter  hltL1sDoubleTauBigOR"
            print "total numerator = ", ER_Calc.numerator
            print "total denominator = ", ER_Calc.denominator
            print "total rate = ", float(ER_Calc.numerator)*75817.94/float(ER_Calc.denominator)
            '''
            ER_Calc.algo_rate_calculator("hltDoubleL2IsoTau26eta2p2")
            print "\nafter  hltDoubleL2IsoTau26eta2p2"
            print "total numerator = ", ER_Calc.numerator
            print "total denominator = ", ER_Calc.denominator
            print "total rate = ", float(ER_Calc.numerator)*75817.94/float(ER_Calc.denominator)
