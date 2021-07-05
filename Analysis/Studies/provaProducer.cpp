#include "DataSetProducer.h"
#include "DataSetProducer.C"

void MergeAllHistos(DataSetProducer& dataset){
  const std::vector<std::string> WJSignals = {absolute_path+"WJetsToLNu/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8_output*.root"};
  const std::vector<std::string> TTSignals = {absolute_path+"TT/TT_TuneCP5_14TeV-powheg-pythia8_output*.root", absolute_path+"TT_ext1/TT_TuneCP5_14TeV-powheg-pythia8_output*.root", absolute_path+"TTToSemiLeptonic/TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8_output*.root" };
  const std::vector<std::string> VBFSignals = {absolute_path+"VBFHToTauTau/VBFHToTauTau_M125_TuneCUETP8M1_14TeV_powheg_pythia8_output_*.root"};
  const std::vector<std::string> DYSignals = {absolute_path+"DYToLL_M-50/DYToLL_M-50_TuneCP5_14TeV-pythia8_output*.root"};
  const std::vector<std::string> QCDSample = {absolute_path+"QCD_Pt-15to3000_Flat/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8_output_*.root"};

  dataset.ProduceDataSamples();

}

void ProvaProducer(){
  DataSetProducer dataset(DYSignals,VBFSignals, TTSignals, WJSignals, QCDSample);
  MergeAllHistos(dataset);


}
