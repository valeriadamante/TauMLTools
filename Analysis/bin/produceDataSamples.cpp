#include <iostream>
#include<cmath>
#include <vector>
#include <iterator>
#include <algorithm> // for std::transform
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "ROOT/RDataFrame.hxx"
#include "TRandom3.h"
std::string absolute_path = "/home/users/damante/L2SkimmedTuples/";
int produceDataSamples(){

  /* open DataFrames */

  const std::vector<std::string> WJSignals = {absolute_path+"WJetsToLNu/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8_output*.root"};
  const std::vector<std::string> TTSignals = {absolute_path+"TT/TT_TuneCP5_14TeV-powheg-pythia8_output*.root", absolute_path+"TT_ext1/TT_TuneCP5_14TeV-powheg-pythia8_output*.root", absolute_path+"TTToSemiLeptonic/TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8_output*.root" };
  const std::vector<std::string> VBFSignals = {absolute_path+"VBFHToTauTau/VBFHToTauTau_M125_TuneCUETP8M1_14TeV_powheg_pythia8_output_*.root"};
  const std::vector<std::string> DYSignals = {absolute_path+"DYToLL_M-50/DYToLL_M-50_TuneCP5_14TeV-pythia8_output*.root"};
  const std::vector<std::string> ZPrimeSignals = {absolute_path+"ZprimeToTauTau_M-4000/ZprimeToTauTau_M-4000_TuneCP5_14TeV-pythia8-tauola_output*.root"};
  const std::vector<std::string> BackgroundsExt = {absolute_path+"QCD_Pt-15to3000_Flat/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8_output_*.root"};
  const std::vector<std::string> Background = { absolute_path+"QCD_Pt-15to3000_Flat/QCDFiltered.root"};

  ROOT::RDataFrame dfWJSignals("L2TauTrainTuple",WJSignals);
  ROOT::RDataFrame dfTTSignals("L2TauTrainTuple",TTSignals);
  ROOT::RDataFrame dfVBFSignals("L2TauTrainTuple",VBFSignals);
  ROOT::RDataFrame dfDYSignals("L2TauTrainTuple",DYSignals);
  ROOT::RDataFrame dfZPrimeSignals("L2TauTrainTuple",ZPrimeSignals);
  ROOT::RDataFrame dfBackgroundExt("L2TauTrainTuple", BackgroundsExt);
  ROOT::RDataFrame dfBackground("L2TauTrainTuple", Background);
  //dfBackground.Display({"l1Tau_pt"})->Print();

  // WJets, TT + TText1  + TTToSemiLep, DY, ZPrime, QCD
  std::vector<int> all_events = {26362, 139019, 182810, 198525, 792712};
  std::vector <int> all_events_norm = {4, 21, 27, 30, 118};
  dfWJSignals.Snapshot("L2TauTrainTuple", absolute_path+"DataSetTraining/all_WJetsToLNu.root");
  dfTTSignals.Snapshot("L2TauTrainTuple", absolute_path+"DataSetTraining/all_TT.root" );
  dfVBFSignals.Snapshot("L2TauTrainTuple", absolute_path+"DataSetTraining/all_VBF.root");
  dfDYSignals.Snapshot("L2TauTrainTuple", absolute_path+"DataSetTraining/all_DY.root");
  dfZPrimeSignals.Snapshot("L2TauTrainTuple", absolute_path+"DataSetTraining/all_ZPrime.root");
  dfBackgroundExt.Snapshot("L2TauTrainTuple", absolute_path+"DataSetTraining/all_QCD.root");
  dfBackground.Snapshot("L2TauTrainTuple", absolute_path+"DataSetTraining/all_QCD_filtered.root");

  return 0;

}
