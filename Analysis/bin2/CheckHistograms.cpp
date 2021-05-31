#include "DataSetProducer.h"
#include "DataSetProducer.C"


void CheckHistograms(){
  DataSetProducer dataset;
  dataset.SetAbsolutePath("/Users/valeriadamante/Desktop/Dottorato/gridui/L2SkimmedTuples/");
  /*
  //ßß Signal files
  std::vector<std::string> SignalFiles = {dataset.GetAbsolutePath()+"DataSetTraining/all_TT.root",dataset.GetAbsolutePath()+"DataSetTraining/all_DY.root",dataset.GetAbsolutePath()+"DataSetTraining/all_WJetsToLNu.root"};
  dataset.SetSignalFiles(SignalFiles);

  //ßß QCD files
  std::vector<std::string> QCDFile = {dataset.GetAbsolutePath()+"DataSetTraining/all_QCD.root"};
  dataset.SetQCDFile(QCDFile);
  std::vector<std::string> QCDFilteredFile = {dataset.GetAbsolutePath()+"DataSetTraining/QCDFiltered.root"};
  dataset.SetQCDFilteredFile(QCDFilteredFile);

  auto colNames = ROOT::RDataFrame("L2TauTrainTuple", QCDFile).GetColumnNames();
  int i =0;
  for (auto &&colName : colNames){
    std::cout << "Getting histogram for  " << colName << std::endl;
    dataset.GetHistogramsSignalQCD(i, true);
    i++;
  }*/
  // Data files
  std::vector<std::string> DataFile = {dataset.GetAbsolutePath()+"DataSetTraining/all_Data.root"};
  dataset.SetDataFile(DataFile);

  // VBF files
  std::vector<std::string> VBFFile = {dataset.GetAbsolutePath()+"DataSetTraining/all_VBF.root"};
  dataset.SetVBFFile(VBFFile);

  // ZP files
  std::vector<std::string> ZPFile = {dataset.GetAbsolutePath()+"DataSetTraining/all_ZPrime.root"};
  dataset.SetZPFile(ZPFile);

  auto colNames = ROOT::RDataFrame("L2TauTrainTuple", VBFFile).GetColumnNames();
  int i =0;
  for (auto &&colName : colNames){
    std::cout << "Getting histogram for  " << colName << std::endl;
    dataset.GetHistogramsDataValidation(i, true);
    i++;
  } 

}
