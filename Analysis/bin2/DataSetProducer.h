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

#ifndef DATASETPRODUCER_H
#define DATASETPRODUCER_H

class DataSetProducer{
public:
    using vec_f = ROOT::VecOps::RVec<float>;
    using vec_i = ROOT::VecOps::RVec<int>;
    DataSetProducer();
    DataSetProducer(const std::vector<std::string> _SignalFiles, const std::vector<std::string> _BackgroundFiles);
    DataSetProducer(const std::vector<std::string> _SignalFiles, const std::vector<std::string> _BackgroundFiles, const std::vector<std::string> _QCDFilteredFile);

    /* Setters */

    void SetAbsolutePath(std::string _absolute_path) {absolute_path = _absolute_path;}
    void SetSignalFiles(std::vector<std::string> _SignalFiles){ SignalFiles=_SignalFiles;}
    void SetWJFile(std::vector<std::string> _WJFile){WJFile=_WJFile;}
    void SetDYFile(std::vector<std::string> _DYFile){DYFile=_DYFile;}
    void SetTTFile(std::vector<std::string> _TTFile){TTFile=_TTFile;}
    void SetVBFFile(std::vector<std::string> _VBFFile){VBFFile=_VBFFile;}
    void SetZPFile(std::vector<std::string> _ZPFile){ZPFile=_ZPFile;}
    void SetEph1File(std::vector<std::string> _Eph1File){Eph1File=_Eph1File;}
    void SetEph2File(std::vector<std::string> _Eph2File){Eph2File=_Eph2File;}
    void SetEph3File(std::vector<std::string> _Eph3File){Eph3File=_Eph3File;}
    void SetEph4File(std::vector<std::string> _Eph4File){Eph4File=_Eph4File;}
    void SetEph5File(std::vector<std::string> _Eph5File){Eph5File=_Eph5File;}
    void SetEph6File(std::vector<std::string> _Eph6File){Eph6File=_Eph6File;}
    void SetEph7File(std::vector<std::string> _Eph7File){Eph7File=_Eph7File;}
    void SetEph8File(std::vector<std::string> _Eph8File){Eph8File=_Eph8File;}
    void SetDataFile(std::vector<std::string> _DataFile){DataFile=_DataFile;}
    void SetQCDFile(std::vector<std::string> _QCDFile){QCDFile=_QCDFile;}
    void SetQCDFilteredFile(std::vector<std::string> _QCDFilteredFile){QCDFilteredFile=_QCDFilteredFile;}
    void SetDataSetTrainingFile(std::vector<std::string> _DataSetTrainingFile){DataSetTrainingFile=_DataSetTrainingFile;}

    /* Getters */
    std::string GetAbsolutePath() {return absolute_path;}
    std::vector<std::string> GetSignalFiles(){return  SignalFiles;}
    std::vector<std::string> GetWJFile(){return WJFile;}
    std::vector<std::string> GetDYFile(){return DYFile;}
    std::vector<std::string> GetTTFile(){return TTFile;}
    std::vector<std::string> GetVBFFile(){return VBFFile;}
    std::vector<std::string> GetZPFile(){return ZPFile;}
    std::vector<std::string> GetEph1File(){return Eph1File;}
    std::vector<std::string> GetEph2File(){return Eph2File;}
    std::vector<std::string> GetEph3File(){return Eph3File;}
    std::vector<std::string> GetEph4File(){return Eph4File;}
    std::vector<std::string> GetEph5File(){return Eph5File;}
    std::vector<std::string> GetEph6File(){return Eph6File;}
    std::vector<std::string> GetEph7File(){return Eph7File;}
    std::vector<std::string> GetEph8File(){return Eph8File;}
    std::vector<std::string> GetDataFile(){return DataFile;}
    std::vector<std::string> GetQCDFile(){return QCDFile;}
    std::vector<std::string> GetQCDFilteredFile(){return QCDFilteredFile;}
    std::vector<std::string> GetDataSetTrainingFile(){return DataSetTrainingFile;}

    /* Other public functions */
    void CountEventInSignalAndBackground();
    void GetHistogramsSignalQCD(int n_var, bool use_binning);
    void GetHistogramsDataValidation(int n_var, bool use_binning);
    void QCDFilter();
    void ProduceDataSample(std::vector<std::string> inputFile, std::string outFile_suffix);

private:
    template <typename T>
    void plot(T sig, T bkgFiltered, T bkg, const std::string &VarName);
    template <typename T>
    void plot1(T h1,   const std::string &VarName, const std::string& outDir_path);
    template <typename T>
    void plot2(T h1, T h2,  const std::string &VarName, const std::string& outDir_path);
    template <typename T>
    void plot3(T h1, T h2, T h3, const std::string &VarName, const std::string& outDir_path) ;
    bool FilterEventsFromRatio(std::vector<float>& ratios, std::vector<float>& low_edges, float L1Pt);
    std::vector<int> CountEvents();
private:
  std::string absolute_path ;
  std::vector<std::string> DYFile,  TTFile,  WJFile,   SignalFiles;
  std::vector<std::string> QCDFile, QCDFilteredFile;
  std::vector<std::string> Eph1File, Eph2File, Eph3File, Eph4File, Eph5File, Eph6File, Eph7File, Eph8File, DataFile;
  std::vector<std::string> VBFFile,ZPFile;
  std::map<int, std::string> all_branches;
  std::map<std::string, std::vector<float>> all_branches_binning;
  std::vector<std::string> DataSetTrainingFile;
};
#endif
