#include "DataSetProducer.h"
#include "DataSetProducer.C"


void CountEvents(){
  DataSetProducer dataset;
  dataset.SetAbsolutePath("/Users/valeriadamante/Desktop/Dottorato/L2SkimmedTuples/");
  //dataset.SetAbsolutePath("/home/users/damante/L2SkimmedTuples/");
  // Signal files
  std::vector<std::string> SignalFiles = {dataset.GetAbsolutePath()+"DataSetTraining/all_TT.root",dataset.GetAbsolutePath()+"DataSetTraining/all_DY.root",dataset.GetAbsolutePath()+"DataSetTraining/all_WJetsToLNu.root"};
  dataset.SetSignalFiles(SignalFiles);
  // QCD files
  std::vector<std::string> QCDFile = {dataset.GetAbsolutePath()+"DataSetTraining/all_QCD.root"};
  dataset.SetQCDFile(QCDFile);
  std::vector<std::string> QCDFilteredFile = {dataset.GetAbsolutePath()+"DataSetTraining/QCDFiltered.root"};
  dataset.SetQCDFilteredFile(QCDFilteredFile);
  // Data files
  std::vector<std::string> DataFile = {dataset.GetAbsolutePath()+"DataSetTraining/all_Data_forReweighting.root"};
  dataset.SetDataFile(DataFile);
  // VBF files
  std::vector<std::string> VBFFile = {dataset.GetAbsolutePath()+"DataSetTraining/all_VBF.root"};
  dataset.SetVBFFile(VBFFile);
  // ZP files
  std::vector<std::string> ZPFile = {dataset.GetAbsolutePath()+"DataSetTraining/all_ZPrime.root"};
  dataset.SetZPFile(ZPFile);
  //dataset.CountEventInSignalAndBackground();
  dataset.FromDFToTTree();

}
