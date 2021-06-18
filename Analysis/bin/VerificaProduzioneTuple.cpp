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
std::string absolute_path="/home/users/damante/L2SkimmedTuples/DataSetTraining/";

int VerificaProduzioneTuple(){

  /* open DataFrames */

  const std::vector<std::string> signals = {absolute_path+"WJetsToLNu/*.root", absolute_path+"TT/*.root", absolute_path+"TT_ext1/*.root", absolute_path+"DYToLL_M-50/*.root", absolute_path+"TTToSemiLeptonic/*.root" };
  const std::vector<std::string> backgrounds = {absolute_path+"QCD_Pt-15to3000_Flat/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8_output_*.root"};
  ROOT::RDataFrame df_DY("L2TauTrainTuple",absolute_path+"all_DY.root");
  ROOT::RDataFrame df_TT("L2TauTrainTuple",absolute_path+"all_TT.root");
  ROOT::RDataFrame df_WJetsToLNu("L2TauTrainTuple",absolute_path+"all_WJetsToLNu.root");
  ROOT::RDataFrame dfQCDFiltered("L2TauTrainTuple",absolute_path+"QCDFiltered.root");
  ROOT::RDataFrame df_DataSetTraining("L2TauTrainTuple",absolute_path+"DataSetTraining.root");


  //df_DY.Range(49);
  std::cout << " df_DY " << std::endl ;
  df_DY.Display({"l1Tau_pt"},49)->Print();

  //df_TT.Range(37);
  std::cout << " df_TT " << std::endl ;
  df_TT.Display({"l1Tau_pt"},37)->Print();

  //df_WJetsToLNu.Range(7);
  std::cout << " df_WJetsToLNu " << std::endl ;
  df_WJetsToLNu.Display({"l1Tau_pt"},7)->Print();

  //dfQCDFiltered.Range(107);
  std::cout << " dfQCDFiltered " << std::endl ;
  dfQCDFiltered.Display({"l1Tau_pt"},107)->Print();

  //df_DataSetTraining.Range();
  std::cout << " df_DataSetTraining " << std::endl ;
  df_DataSetTraining.Display({"l1Tau_pt"}, 200)->Print();


  return 0;

}
