#include <iostream>
#include<cmath>
#include <vector>
#include <iterator>
#include <algorithm> // for std::transform
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TH1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "ROOT/RDataFrame.hxx"
#include "TRandom3.h"
#include "TLegend.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include <Math/Vector4Dfwd.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include "TStyle.h"
#include <string>


std::string absolute_path = "/Users/valeriadamante/Desktop/Dottorato/gridui/L2SkimmedTuples/DataSetTraining/";
std::string outDir_path = "/Users/valeriadamante/Desktop/Dottorato/gridui/CMSSW_11_2_1_Patatrack/src/outVars/";
//std::string absolute_path = "/home/users/damante/L2SkimmedTuples/DataSetTraining/";
int offset_value =10;

template <typename T>
void plot(T sig, T bkgFiltered, T bkg, const std::string &VarName)
{
   // Canvas and general style options
   gStyle->SetOptStat(0);
   gStyle->SetTextFont(42);
   auto c = new TCanvas("c", "", 800, 700);
   //TLegend *legend = new TLegend(0.421053, 0.82087 , 0.552632,0.902609 );
   //TLegend *legend = new TLegend( 0.62, 0.70, 0.82, 0.88 );
   TLegend *legend = new TLegend(0.74812, 0.820741, 0.901003, 0.902222 );
   c->SetLeftMargin(0.15);
   // Get signal and background histograms and stack them to show Higgs signal
   // on top of the background process
   auto hsig = *sig;
   auto hbkgFiltered = *bkgFiltered;
   auto hbkg = *bkg;
   hsig.SetTitle(VarName.c_str());
   hsig.SetStats(0);
   hsig.GetXaxis()->SetTitle(VarName.c_str());
   hsig.GetYaxis()->SetTitle("A.U.");
   hsig.SetLineWidth(2);
   hsig.SetLineColor(kBlue);
   //hsig.SetFillColor(kBlue);
   //hsig.SetFillStyle(3001);
   //hsig.SetMarkerStyle(20);
   //hsig.Rebin(2);
   hsig.SetMarkerColor(kBlue);
   hsig.Scale(100/hsig.GetEntries());
   legend->AddEntry(&hsig, "signal", "l");

   hbkgFiltered.SetTitle(VarName.c_str());
   hbkgFiltered.SetStats(0);
   hbkgFiltered.GetXaxis()->SetTitle(VarName.c_str());
   hbkgFiltered.GetYaxis()->SetTitle("A.U.");
   hbkgFiltered.SetLineColor(kRed);
   //hbkgFiltered.SetFillColor(kRed);
   //hbkgFiltered.SetFillStyle(3001);
   //hbkgFiltered.SetMarkerStyle(20);
   //hbkgFiltered.Rebin(2);
   hbkgFiltered.SetMarkerColor(kRed);
   hbkgFiltered.Scale(100/hbkgFiltered.GetEntries());
   legend->AddEntry(&hbkgFiltered, "bkgFiltered", "l");

   hbkg.SetTitle(VarName.c_str());
   hbkg.SetStats(0);
   hbkg.GetXaxis()->SetTitle(VarName.c_str());
   hbkg.GetYaxis()->SetTitle("A.U.");
   hbkg.SetLineColor(kGreen);
   //hbkg.SetFillColor(kGreen);
   //hbkg.SetFillStyle(3001);
   //hbkg.SetMarkerStyle(20);
   //hbkg.Rebin(2);
   hbkg.SetMarkerColor(kRed);
   hbkg.Scale(100/hbkg.GetEntries());
   legend->AddEntry(&hbkg, "bkg", "l");

   if(hbkg.GetBinContent(hbkg.GetMaximumBin()) > hsig.GetBinContent(hsig.GetMaximumBin()) ){
     hbkg.DrawClone("HIST");
     hsig.DrawClone("HISTSAME");
     hbkgFiltered.DrawClone("HISTSAME");
   }
   else{
     hsig.DrawClone("HIST");
     hbkg.DrawClone("HISTSAME");
     hbkgFiltered.DrawClone("HISTSAME");
   }

    legend->Draw("SAME");

   // Save plot
   c->SaveAs((outDir_path+VarName+".pdf").c_str());
}

std::vector<int> CountEvents(const std::vector<std::string> signals){
  std::vector<int>  total_events;
  for (auto& i: signals){
    ROOT::RDataFrame df("L2TauTrainTuple", i);
    total_events.push_back(df.Count().GetValue());
  }
  return total_events;
}
using vec_f = ROOT::VecOps::RVec<float>;
using vec_i = ROOT::VecOps::RVec<int>;

std::map<int, std::string> all_branches =  { {3, "genEventWeight"},{8, "genLepton_vis_pt"},{9, "genLepton_vis_eta"},{10, "genLepton_vis_phi"},{11, "genLepton_vis_mass"},{13, "l1Tau_pt" },{14, "l1Tau_eta"},{15, "l1Tau_phi"},{16, "l1Tau_mass"} , {25, "caloRecHit_ee_rho"},  {26, "caloRecHit_eb_rho"},  {27, "caloRecHit_ee_eta"},  {28, "caloRecHit_eb_eta"},  {29, "caloRecHit_ee_phi"},  {30, "caloRecHit_eb_phi"},  {31, "caloRecHit_ee_energy"},  {32, "caloRecHit_eb_energy"},  {33, "caloRecHit_ee_time"},  {34, "caloRecHit_eb_time"},  {37, "caloRecHit_ee_chi2"},  {38, "caloRecHit_eb_chi2"},  {39, "caloRecHit_ee_energyError"},  {40, "caloRecHit_eb_energyError"},  {41, "caloRecHit_ee_timeError"},  {42, "caloRecHit_eb_timeError"},  {51, "caloRecHit_hbhe_rho"},  {52, "caloRecHit_hbhe_eta"},  {53, "caloRecHit_hbhe_phi"},  {54, "caloRecHit_hbhe_energy"},  {55, "caloRecHit_hbhe_time"},  {57, "caloRecHit_hbhe_chi2"},  {59, "caloRecHit_ho_rho"},  {60, "caloRecHit_ho_eta"},  {61, "caloRecHit_ho_phi"},  {62, "caloRecHit_ho_energy"},  {63, "caloRecHit_ho_time"},  {67, "caloRecHit_hf_rho"},  {68, "caloRecHit_hf_eta"},  {69, "caloRecHit_hf_phi"},  {70, "caloRecHit_hf_energy"},  {71, "caloRecHit_hf_time"},  {74, "caloRecHit_hf_timeFalling"},  {77, "patatrack_pt"},  {78, "patatrack_eta"},  {79, "patatrack_phi"},  {80, "patatrack_chi2"},  {84, "patatrack_dxy"},  {85, "patatrack_dz"},  {87, "patavert_z"},  {88, "patavert_weight"},  {89, "patavert_ptv2"},  {90, "patavert_chi2"} , {35, "caloRecHit_ee_detId"},{36, "caloRecHit_eb_detId"},{43, "caloRecHit_ee_flagsBits"},{44, "caloRecHit_eb_flagsBits"},{45, "caloRecHit_ee_isRecovered"},{46, "caloRecHit_eb_isRecovered"},{47, "caloRecHit_ee_isTimeValid"},{48, "caloRecHit_eb_isTimeValid"},{49, "caloRecHit_ee_isTimeErrorValid"},{50, "caloRecHit_eb_isTimeErrorValid"},{56, "caloRecHit_hbhe_detId"},{58, "caloRecHit_hbhe_flags"},{64, "caloRecHit_ho_detId"},{65, "caloRecHit_ho_aux"},{66, "caloRecHit_ho_flags"},{72, "caloRecHit_hf_detId"},{73, "caloRecHit_hf_flags"},{75, "caloRecHit_hf_auxHF"},{76, "caloRecHit_hf_aux"},{81, "patatrack_ndof"},{82, "patatrack_charge"},{83, "patatrack_quality"},{86, "patatrack_vertex_id"},{91, "patavert_ndof"}, {0, "run"},{1, "lumi"},{2, "evt"},{4, "sampleType"},{5, "genLepton_index"},{6, "genLepton_kind"},{7, "genLepton_charge"},{12, "l1Tau_index"},{17, "l1Tau_hwIso"},{18, "l1Tau_hwQual"},{19, "l1Tau_towerIEta"},{20, "l1Tau_towerIPhi"},{21, "l1Tau_rawEt"},{22, "l1Tau_isoEt"},{23, "l1Tau_hasEM"},{24, "l1Tau_isMerged"} };


std::map<std::string, std::vector<float>> all_branches_binning = {{"genEventWeight", {50, -0.5,1.5}},{"genLepton_vis_pt",{100, 0, 1000}},{"genLepton_vis_eta", {100, -4.5, 5.5}},{"genLepton_vis_phi", {70, -3.5, 3.5}},{"genLepton_vis_mass", {100, -0.1, 9.9}},{"l1Tau_pt",{50, 0, 250}},{"l1Tau_eta", {100, -4.5, 5.5}},{"l1Tau_phi", {70, -3.5, 3.5}},{"l1Tau_mass",  {100, -0,1, 9.9}},{"caloRecHit_ee_rho", {100, 0, 300}},{"caloRecHit_eb_rho", {100, 0, 300}},{"caloRecHit_ee_eta",{100, -4.5, 5.5}},{"caloRecHit_eb_eta",{100, -4.5, 5.5}},{"caloRecHit_ee_phi",{70, -3.5, 3.5}},{"caloRecHit_eb_phi",{70, -3.5, 3.5}},{"caloRecHit_ee_energy",{60, 0, 6}},{"caloRecHit_eb_energy",{100, 0, 80}},{"caloRecHit_ee_time",{50, -0.5, 9.5}},{"caloRecHit_eb_time",{50, -0.5, 9.5}},{"caloRecHit_ee_chi2",{100, -0.5, 60.5}},{"caloRecHit_eb_chi2",{100, -0.5, 60.5}},{"caloRecHit_ee_energyError",{50, -0.5, 9.5}},{"caloRecHit_eb_energyError",{50, -0.5, 9.5}},{"caloRecHit_ee_timeError",{50, -0.5, 9.5}},{"caloRecHit_eb_timeError",{50, -0.5, 9.5}},{"caloRecHit_hbhe_rho", {100, 0, 300}},{"caloRecHit_hbhe_eta",{100, -4.5, 5.5}},{"caloRecHit_hbhe_phi",{70, -3.5, 3.5}},{"caloRecHit_hbhe_energy",{100, 0, 80}},{"caloRecHit_hbhe_time",{50, -0.5, 9.5}},{"caloRecHit_hbhe_chi2",{100, -0.5, 60.5}},{"caloRecHit_ho_rho", {100, 0, 300}},{"caloRecHit_ho_eta",{100, -4.5, 5.5}},{"caloRecHit_ho_phi",{70, -3.5, 3.5}},{"caloRecHit_ho_energy",{100, 0, 80}},{"caloRecHit_ho_time",{50, -0.5, 9.5}},{"caloRecHit_hf_rho", {100, 0, 300}}, {"caloRecHit_hf_eta",{100, -4.5, 5.5}},{"caloRecHit_hf_phi",{70, -3.5, 3.5}},{"caloRecHit_hf_energy",{100, 0, 80}},{"caloRecHit_hf_time",{50, -0.5, 9.5}},{"caloRecHit_hf_timeFalling",{50, -0.5, 9.5}},{"patatrack_pt",{50, -0.5, 20.5}},{"patatrack_eta",{100, -4.5, 5.5}},{"patatrack_phi",{70, -3.5, 3.5}},{"patatrack_chi2",{100, -0.5, 60.5}},{"patatrack_dxy", {60, -0.3, 0.3}},{"patatrack_dz", {300, -15.5, 15.5}},{"patavert_z",{300, -15.5, 15.5}},{"patavert_weight", {100, 0, 2000}},{"patavert_ptv2", {100}},{"patavert_chi2",{100, -0.5, 1000.5}},{"caloRecHit_ee_detId", {1000}},{"caloRecHit_eb_flagsBits", {100}},{"caloRecHit_ee_isRecovered", {100}},{"caloRecHit_eb_isRecovered", {100}},{"caloRecHit_ee_isTimeValid", {100}},{"caloRecHit_eb_isTimeValid", {100}},{"caloRecHit_ee_isTimeErrorValid", {100}},{"caloRecHit_eb_isTimeErrorValid", {100}},{"caloRecHit_hbhe_detId", {100}},{"caloRecHit_hbhe_flags", {100}},{"caloRecHit_ho_detId", {100}},{"caloRecHit_ho_aux", {100}},{"caloRecHit_hf_flags", {100}},{"caloRecHit_hf_auxHF", {100}}, {"patatrack_charge", {50, -0.5, 1.5}},{"patatrack_quality", {50, -0.5, 1.5}},{"patatrack_vertex_id", {50, 0, 100}},{"patavert_ndof", {100, 0, 200}},{"run", {10, -0.5,1.5}},{"lumi", {1000}},{"evt", {1000}},{"sampleType", {10, -0.5,1.5}},{"genLepton_index", {60, -0.1, 5,9}},{"genLepton_kind", {100, -0.5, 9.5}},{"genLepton_charge", {10, -0.5, 1.5}}, {"caloRecHit_eb_detId", {1000}}, {"caloRecHit_ee_flagsBits", {100}}, {"caloRecHit_hf_detId", {100}}, {"caloRecHit_hf_aux", {100}},{"patatrack_ndof", {10, 0, 10}}, {"l1Tau_hwIso", {10, -0.5,1.5}}, {"l1Tau_hwQual", {10, -0.5,1.5}}, {"l1Tau_towerIEta",{10, -0.5,1.5}}, {"l1Tau_towerIPhi",{10, -0.5,1.5}}, {"l1Tau_rawEt",{10, -0.5,1.5}}, {"l1Tau_isoEt",{10, -0.5,1.5}}, {"l1Tau_hasEM",{10, -0.5,1.5}}, {"l1Tau_isMerged" ,{10, -0.5,1.5}}, {"l1Tau_index", {50, 0, 10}}, {"caloRecHit_ho_flags", {100}} };

void GetHistogramsAndNumbers(int CountOrPlot, int n_var=0, bool use_binning = true){

  /* CountOrPlot 0 == count / CountOrPlot 1 == plot */

  /* open DataFrames */

  const std::vector<std::string> signals = {absolute_path+"all_DY.root",absolute_path+"all_TT.root", absolute_path+"all_WJetsToLNu.root"};
  const std::vector<std::string> backgroundFiltered = {absolute_path+"QCDFiltered.root"};
  const std::vector<std::string> backgrounds = {absolute_path+"all_QCD.root"};

  ROOT::RDataFrame dfSignal("L2TauTrainTuple", signals);
  ROOT::RDataFrame dfBackground("L2TauTrainTuple", backgrounds);
  ROOT::RDataFrame dfQCDFiltered("L2TauTrainTuple", backgroundFiltered);

  /* count events */
  if(CountOrPlot == 0){
    int previous_bckg_events = dfBackground.Count().GetValue();
    int bckg_current_events = dfQCDFiltered.Count().GetValue();
    std::cout<< "before "<< dfBackground.Count().GetValue()<<" after "<<bckg_current_events << "\n signal events = " << dfSignal.Count().GetValue() <<std::endl;

    std::vector<int> all_events = CountEvents(signals);
    all_events.push_back(bckg_current_events);
    int sum = 0;
    for (auto&i : all_events){
      std::cout <<i << std::endl;
      sum += i;
    }
    // sum = 3303918;
    std::cout<<"sum = " << sum<<std::endl;
    std::vector<float> all_events_norm(all_events.size());
    float rescale_to_duecento = 200./float(sum);
    std::transform(all_events.begin(), all_events.end(), all_events_norm.begin(),
                   [&](int i) { return i * rescale_to_duecento; });
    for (auto&i : all_events_norm){
       std::cout<<i<<std::endl;
    }
  }
  else {
    if (all_branches.find(n_var) != all_branches.end()){
    std::string variabile = all_branches.at(n_var);
    dfSignal.Display(variabile)->Print();
      std::cout << "************************************************************" << std::endl ;
      std::cout << "plotted observable  = " << variabile << std::endl;
      std::cout << "************************************************************" << std::endl ;
      auto maximum_histogram_value = std::max(dfQCDFiltered.Max(variabile).GetValue(), dfSignal.Max(variabile).GetValue())*(1.5);
      std::cout << "previously maximum_histogram_value was " << maximum_histogram_value <<std::endl;
      auto minimum_histogram_value = std::min(dfQCDFiltered.Min(variabile).GetValue(), dfSignal.Min(variabile).GetValue());
      std::cout << "previously minimum_histogram_value was " << minimum_histogram_value <<std::endl;
      auto bins = maximum_histogram_value-minimum_histogram_value;
      int number_of_bin = static_cast<int>(bins);
      std::cout << "previously number of bins was " << number_of_bin <<std::endl;
      std::cout << std::endl ;
      if(minimum_histogram_value<0){
        minimum_histogram_value = minimum_histogram_value*1.5;
        std::cout << "after rescaling maximum_histogram_value is " << maximum_histogram_value <<std::endl;
      }
      else{
        minimum_histogram_value = minimum_histogram_value*0.75;
        std::cout << "after rescaling minimum_histogram_value is " << minimum_histogram_value <<std::endl;
      }
      if(minimum_histogram_value == maximum_histogram_value){
        maximum_histogram_value = minimum_histogram_value+offset_value;
        std::cout << "min = max, so we changed them into \nminimum =  " << minimum_histogram_value << "\nmaximum = " << maximum_histogram_value <<std::endl;
      }
      if (number_of_bin == 0){
        number_of_bin = 10;
      }
      else if (number_of_bin> 0 && number_of_bin <= 10){
        number_of_bin = number_of_bin*10;
      }
      else{
        number_of_bin = 100;
      }
      std::cout << "after rescaling number of bin was " << number_of_bin <<std::endl;
      std::cout << std::endl ;

      if(use_binning && all_branches_binning.at(variabile).size()>1){
        number_of_bin = static_cast<int>(all_branches_binning.at(variabile).at(0));
        minimum_histogram_value = all_branches_binning.at(variabile).at(1);
        maximum_histogram_value = all_branches_binning.at(variabile).at(2);
        std::cout << "number of bins = " << number_of_bin << std::endl;
        std::cout << "minimum_histogram_value = " << minimum_histogram_value << std::endl;
        std::cout << "maximum_histogram_value = " << maximum_histogram_value << std::endl;
      }

      auto signalHist = dfSignal.Histo1D({"signal", "signal", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, variabile);
      auto QCDFilteredHist = dfQCDFiltered.Histo1D({"backgroundFiltered", "backgroundFiltered", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, variabile);
      auto backgroundHist = dfBackground.Histo1D({"background", "background", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, variabile);
      plot(signalHist, QCDFilteredHist, backgroundHist, variabile);


    }
  }

}
