#include <memory>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "ROOT/RDataFrame.hxx"
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
#include "TRandom3.h"
#include "TLegend.h"
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


//using vec_f = ROOT::VecOps::RVec<float>;
//using vec_i = ROOT::VecOps::RVec<int>;

#ifndef NEWTUPLE_H
#define NEWTUPLE_H

class NewTuple{
public:
    using vec_f = ROOT::VecOps::RVec<float>;
    using vec_i = ROOT::VecOps::RVec<int>;
    NewTuple();
    NewTuple(const std::vector<std::string> _newTupleFile, const std::vector<std::string> _oldTupleFile);
    NewTuple(const std::vector<std::string> _newTupleFile, const std::vector<std::string> _oldTupleFile, const std::string _absolute_path);
    void SetAbsolutePath(std::string _absolute_path) {absolute_path = _absolute_path;}
    void SetPlotDir(std::string _plotDir) {plotDir = _plotDir;}
    void SetNewTupleFile(std::vector<std::string> _newTupleFile){newTupleFile=_newTupleFile;}
    void SetOldTupleFile(std::vector<std::string> _oldTupleFile){oldTupleFile=_oldTupleFile;}
    std::string GetAbsolutePath() {return absolute_path;}
    std::string GetPlotDir() {return plotDir;}
    std::vector<std::string> GetNewTupleFile(){return newTupleFile;}
    std::vector<std::string> GetOldTupleFile(){return oldTupleFile;}
    void ValidateFiles();

    void DisplayColumnNames();
    void GetWeight(std::string dataFileName, std::string dataTupleName);
    void DrawOnlyHistos(std::string dataFileName, std::string dataTupleName);
    void LookAtDataHistogram();

private:
  float GetWeightFromHisto(float& l1Tau_pt, float& genLepton_vis_pt, bool& genLepton_isTau, TH1D SignalHist,TH1D QCDHist,TH1D DataHist);
  template <typename T>
  void LookAtHistogram(T Histogram);
private:
  std::string absolute_path, plotDir;
  std::vector<std::string> newTupleFile, oldTupleFile;
  float maximum_histogram_value, minimum_histogram_value;
  int number_of_bin;
  float bin_width;
};
#endif
