#include "TauMLTools/Analysis/Studies/DataSetProducer.h"
#include "TauMLTools/Analysis/Studies/DataSetProducer.C"


void MergeFiles(){
  DataSetProducer dataset;
  //dataset.SetAbsolutePath("/home/users/damante/L2SkimmedTuples/");
  dataset.SetAbsolutePath("/gpfs/ddn/cms/user/damante/L2SkimmedTuples/");
  //std::vector<std::string> WJetsToLNuFile = {dataset.GetAbsolutePath()+"WJetsToLNu/*.root"};
  //dataset.SetWJFile(WJetsToLNuFile);
  //dataset.ProduceDataSample(dataset.GetWJFile(), "WJetsToLNu");
  //std::cout << "WJetsToLNu Finito " << std::endl;

  //std::vector<std::string> TTFile = {dataset.GetAbsolutePath()+"TT/*.root", dataset.GetAbsolutePath()+"TT_ext1/*.root", dataset.GetAbsolutePath()+"TTToSemiLeptonic/*.root"};
  //dataset.SetTTFile(TTFile);
  //dataset.ProduceDataSample(dataset.GetTTFile(), "TT");
  //std::cout << "TT Finito " << std::endl;

  //std::vector<std::string> DYToLL_MFile = {dataset.GetAbsolutePath()+"DYToLL_M-50/*.root"};
  //dataset.SetDYFile(DYToLL_MFile);
  //dataset.ProduceDataSample(dataset.GetDYFile(), "DY");
  //std::cout << "DY Finito " << std::endl;

  //std::vector<std::string> VBFHToTauTauFile = {dataset.GetAbsolutePath()+"VBFHToTauTau/*.root"};
  //dataset.SetVBFFile(VBFHToTauTauFile);
  //dataset.ProduceDataSample(dataset.GetVBFFile(), "VBF");
  //std::cout << "VBF Finito " << std::endl;

  //std::vector<std::string> ZprimeToTauTau_MFile = {dataset.GetAbsolutePath()+"ZprimeToTauTau_M-4000/*.root"};
  //dataset.SetZPFile(ZprimeToTauTau_MFile);
  //dataset.ProduceDataSample(dataset.GetZPFile(), "ZPrime");
  //std::cout << "ZPrime Finito " << std::endl;

  //std::vector<std::string> QCD_PtFile = {dataset.GetAbsolutePath()+"QCD_Pt-15to3000_Flat/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8*.root"};
  //dataset.SetQCDFile(QCD_PtFile);
  //dataset.ProduceDataSample(dataset.GetQCDFile(), "QCD");
  //std::cout << "QCD Finito " << std::endl;

  std::vector<std::string> EphemeralHLTPhysics1File = {dataset.GetAbsolutePath()+"EphemeralHLTPhysics1/*.root"};
  dataset.SetEph1File(EphemeralHLTPhysics1File);
  dataset.ProduceDataSample(dataset.GetEph1File(), "Eph1");
  std::cout << "Eph1 Finito " << std::endl;

  //std::vector<std::string> EphemeralHLTPhysics2File = {dataset.GetAbsolutePath()+"EphemeralHLTPhysics2/*.root"};
  //dataset.SetEph2File(EphemeralHLTPhysics2File);
  //dataset.ProduceDataSample(dataset.GetEph2File(), "Eph2");
  //std::cout << "Eph2 Finito " << std::endl;

  //std::vector<std::string> EphemeralHLTPhysics3File = {dataset.GetAbsolutePath()+"EphemeralHLTPhysics3/*.root"};
  //dataset.SetEph3File(EphemeralHLTPhysics3File);
  //dataset.ProduceDataSample(dataset.GetEph3File(), "Eph3");
  //std::cout << "Eph3 Finito " << std::endl;

  //std::vector<std::string> EphemeralHLTPhysics4File = {dataset.GetAbsolutePath()+"EphemeralHLTPhysics4/*.root"};
  //dataset.SetEph4File(EphemeralHLTPhysics4File);
  //dataset.ProduceDataSample(dataset.GetEph4File(), "Eph4");
  //std::cout << "Eph4 Finito " << std::endl;

  //std::vector<std::string> EphemeralHLTPhysics5File = {dataset.GetAbsolutePath()+"EphemeralHLTPhysics5/*.root"};
  //dataset.SetEph5File(EphemeralHLTPhysics5File);
  //dataset.ProduceDataSample(dataset.GetEph5File(), "Eph5");
  //std::cout << "Eph5 Finito " << std::endl;

  //std::vector<std::string> EphemeralHLTPhysics6File = {dataset.GetAbsolutePath()+"EphemeralHLTPhysics6/*.root"};
  //dataset.SetEph6File(EphemeralHLTPhysics6File);
  //dataset.ProduceDataSample(dataset.GetEph6File(), "Eph6");
  //std::cout << "Eph6 Finito " << std::endl;

  //std::vector<std::string> EphemeralHLTPhysics7File = {dataset.GetAbsolutePath()+"EphemeralHLTPhysics7/*.root"};
  //dataset.SetEph7File(EphemeralHLTPhysics7File);
  //dataset.ProduceDataSample(dataset.GetEph7File(), "Eph7");
  //std::cout << "Eph7 Finito " << std::endl;

  //std::vector<std::string> EphemeralHLTPhysics8File = {dataset.GetAbsolutePath()+"EphemeralHLTPhysics8/*.root"};
  //dataset.SetEph8File(EphemeralHLTPhysics8File);
  //dataset.ProduceDataSample(dataset.GetEph8File(), "Eph8");
  //std::cout << "Eph8 Finito " << std::endl;

  //std::vector<std::string> DataFile = {dataset.GetAbsolutePath()+"EphemeralHLTPhysics1/*.root",dataset.GetAbsolutePath()+"EphemeralHLTPhysics2/*.root",dataset.GetAbsolutePath()+"EphemeralHLTPhysics3///*.root",dataset.GetAbsolutePath()+"EphemeralHLTPhysics4/*.root",dataset.GetAbsolutePath()+"EphemeralHLTPhysics5/*.root",dataset.GetAbsolutePath()+"EphemeralHLTPhysics6///*.root",dataset.GetAbsolutePath()+"EphemeralHLTPhysics7/*.root",dataset.GetAbsolutePath()+"EphemeralHLTPhysics8/*.root"};
  //dataset.SetDataFile(DataFile);
  //dataset.ProduceDataSample(dataset.GetDataFile(), "Data");
  //std::cout << "Data Finito " << std::endl;


}
