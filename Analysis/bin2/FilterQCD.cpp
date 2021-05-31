#include "DataSetProducer.h"
#include "DataSetProducer.C"


void FilterQCD(){
  DataSetProducer dataset;
  dataset.SetAbsolutePath("/home/users/damante/L2SkimmedTuples/");

  /* Signal file */
  std::vector<std::string> SignalFiles = {dataset.GetAbsolutePath()+"DataSetTraining/all_TT.root",dataset.GetAbsolutePath()+"DataSetTraining/all_DY.root",dataset.GetAbsolutePath()+"DataSetTraining/all_WJetsToLNu.root"};
  dataset.SetSignalFiles(SignalFiles);

  /* QCD files */
  std::vector<std::string> QCDFile = {dataset.GetAbsolutePath()+"DataSetTraining/all_QCD.root"};
  dataset.SetQCDFile(QCDFile);

  dataset.QCDFilter();


}
