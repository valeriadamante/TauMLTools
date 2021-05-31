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
using vec_f = ROOT::VecOps::RVec<float>;
using vec_i = ROOT::VecOps::RVec<int>;
std::string absolute_path = "/Users/valeriadamante/Desktop/Dottorato/gridui/L2SkimmedTuples/DataSetTraining/";
std::string outDir_path = "/Users/valeriadamante/Desktop/Dottorato/gridui/CMSSW_11_2_1_Patatrack/src/outVars/";
//std::string absolute_path = "/home/users/damante/L2SkimmedTuples/DataSetTraining/";
float minimum_histogram_value = 10.;
float maximum_histogram_value = 275.;
int number_of_bin = 100;
float bin_width = (maximum_histogram_value-minimum_histogram_value)/static_cast<float>(number_of_bin) ;

std::map<float, double> FillMap(TH1D DataSetHist){
  std::map<float, double> integral_map;
  for (int i=0; i<= DataSetHist.GetNbinsX(); ++i){
    integral_map.insert(std::pair<float, double>(DataSetHist.GetBinLowEdge(i),  DataSetHist.Integral(i, i+1) ));
    //std::cout << "bin = " << i ;
    //std::cout << "\t bin width = " << bin_width;
    //std::cout << "\t low edge = " << DataSetHist.GetBinLowEdge(i);
    //std::cout << "\t up edge = " << DataSetHist.GetBinLowEdge(i)+ bin_width;
    //std::cout << "\t bin content = " << DataSetHist.GetBinContent(i);
    //std::cout << "\t bin content * width = " << DataSetHist.GetBinContent(i)*bin_width/2;
    //std::cout << "\t bin integral = " << DataSetHist.Integral(i, i+1)<< std::endl;
  }
  return integral_map;
}

float GetWeightFromMap(float& l1Tau_pt, int& genLepton_kind, std::map<float, double>  SignalMap, std::map<float, double>  QCDMap){
  float weight = 1. ;
  for(auto& pair : SignalMap){
    if(l1Tau_pt>= pair.first &&  l1Tau_pt<pair.first+bin_width){
      if(genLepton_kind == 5){
        if(pair.second!=0){
          weight = 1/pair.second;
        }
        else{
          weight = 1.;
        }
      }
      else if(genLepton_kind==6){
        if(pair.second!=0){
          weight = QCDMap.at(pair.first)/pair.second;
        }
        else{
          weight = 1.;
        }
      }
    //std::cout<<weight << std::endl;
    }
  }
  return weight;
}

void GetWeight(){
  /* open DataFrames */
  const std::vector<std::string> DataSetPath = {absolute_path+"DataSetTraining.root"};
  ROOT::RDataFrame dfDataSet("L2TauTrainTuple", DataSetPath);
  auto dfSignals = dfDataSet.Filter("genLepton_kind==5");
  auto SignalHist = *dfSignals.Histo1D({"signal", "signal", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, "l1Tau_pt");
  std::map<float, double> SignalMap = FillMap(SignalHist);
  /*for(std::map<float, double>::const_iterator it = SignalMap.begin();
    it != SignalMap.end(); ++it)
  {
      std::cout << "bin low edge = " << it->first << "\t value = " << it->second  << "\n";
  }*/
  auto dfQCD = dfDataSet.Filter("genLepton_kind == 6 ");
  auto QCDHist = *dfQCD.Histo1D({"signal", "signal", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, "l1Tau_pt");
  std::map<float, double> QCDMap = FillMap(QCDHist);
  /*for(std::map<float, double>::const_iterator it = QCDMap.begin();
    it != QCDMap.end(); ++it)
  {
      std::cout << " qcd map bin low edge = " << it->first << "\t value = " << it->second  << "\n";
  }*/
  auto GetWeightFromMaps= [&SignalMap, &QCDMap](float &L1Pt, int &genlepKind){
      return GetWeightFromMap( L1Pt,genlepKind, SignalMap, QCDMap);
  };

  auto dfDataSetWeight = dfDataSet.Define("weight", GetWeightFromMaps, {"l1Tau_pt", "genLepton_kind"});
  //dfDataSetWeight.Display("weight", 201)->Print();
  dfDataSetWeight.Snapshot("L2TauTrainTuple", absolute_path+"DataSetTrainingWeight.root");

}
