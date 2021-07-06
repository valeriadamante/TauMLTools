#include "NewTuple.h"
/* initializer with anything */
NewTuple::NewTuple(){}
/* constructor that takes newTuple file and oldTuple file*/
NewTuple::NewTuple(const std::vector<std::string> _newTupleFile, const std::vector<std::string> _oldTupleFile){
  newTupleFile = _newTupleFile;
  oldTupleFile = _oldTupleFile;
}
/* constructor that takes newTuple file oldTuple and absolute path file*/
NewTuple::NewTuple(const std::vector<std::string> _newTupleFile, const std::vector<std::string> _oldTupleFile, const std::string _absolute_path){
  absolute_path = _absolute_path;
  newTupleFile = _newTupleFile;
  oldTupleFile = _oldTupleFile;
}

/* check column names for old and new tuple */
void NewTuple::DisplayColumnNames(){
  ROOT::RDataFrame dfNew("L2TauTrainTuple", newTupleFile);
  ROOT::RDataFrame dfOld("L2TauTrainTuple", oldTupleFile);
  std::cout << "old column Names \n[" ;
  for (auto& col : dfOld.GetColumnNames()){
    std::cout << col << ",\t";
  }
  std::cout << "]\nnew column Names \n[" ;
  for (auto& col : dfNew.GetColumnNames()){
    std::cout << col << ",\t";
  }
  std::cout << "]";
}

/* check if minituple is correctly filled  */
void NewTuple::ValidateFiles(){
  ROOT::RDataFrame dfNew("L2TauTrainTuple", newTupleFile);
  ROOT::RDataFrame dfOld("L2TauTrainTuple", oldTupleFile);
  //"caloRecHit_ee_energy", "caloRecHit_eb_energy",
  std::cout << "old calorechits ecal energy" << std::endl;
  std::cout << dfOld.Define("calorechit_ee_size", [](ROOT::VecOps::RVec<float> vector){ROOT::VecOps::RVec<float> vector2; for(auto& i: vector){if(i!=0){vector2.push_back(i);}}return vector2.size();}, {"caloRecHit_ee_energy"}).Define("calorechit_eb_size", [](ROOT::VecOps::RVec<float> vector){ROOT::VecOps::RVec<float> vector2; for(auto&i: vector){if(i!=0){vector2.push_back(i);}}return vector2.size();}, {"caloRecHit_eb_energy"}).Display({"calorechit_ee_size", "calorechit_eb_size" }, 10)->AsString()<< std::endl;
  std::cout << dfOld.Define("calorechit_ee_noZeroes", [](ROOT::VecOps::RVec<float> vector){ROOT::VecOps::RVec<float> vector2; for(auto& i: vector){if(i!=0){vector2.push_back(i);}}return vector2;}, {"caloRecHit_ee_energy"}).Define("calorechit_eb_noZeroes", [](ROOT::VecOps::RVec<float> vector){ROOT::VecOps::RVec<float> vector2; for(auto& i: vector){if(i!=0){vector2.push_back(i);}}return vector2;}, {"caloRecHit_eb_energy"}).Display({"calorechit_ee_noZeroes","calorechit_eb_noZeroes"},1)->AsString()<<std::endl;
  //"caloRecHit_e_energy",
  std::cout << "new calorechits ecal energy" << std::endl;
  std::cout << dfNew.Define("calorechit_e_size", [](ROOT::VecOps::RVec<float> vector){return vector.size();}, {"caloRecHit_e_energy"}).Display({ "calorechit_e_size"}, 10)->AsString()<< std::endl;
  std::cout << dfNew.Display({"caloRecHit_e_energy"},1)->AsString() << std::endl;
}

/* Get info for histogram */
template <typename T>
void NewTuple::LookAtHistogram(T Histogram){
  auto histogram = *Histogram;
  int j = 0;
  for (int i=0; i<= histogram.GetNbinsX()+1; ++i){
    std::cout << "bin = " << i ;
    std::cout << "\t bin width = " << histogram.GetBinWidth(i);
    std::cout << "\t low edge = " << histogram.GetBinLowEdge(i);
    std::cout << "\t up edge = " << histogram.GetBinLowEdge(i)+ histogram.GetBinWidth(i);
    std::cout << "\t bin content = " << histogram.GetBinContent(i);
    std::cout<< std::endl;
    j+= histogram.GetBinContent(i);
  }
  std::cout << "nentries tot = " << j << std::endl;
}

float NewTuple::GetWeightFromHisto(float& l1Tau_pt, float& genLepton_vis_pt, bool& genLepton_isTau, TH1D SignalHist,TH1D QCDHist,TH1D DataHist){
  float weight =1.;
  if(genLepton_isTau==0){
    for (int i=1; i<= QCDHist.GetNbinsX(); ++i){
      if(QCDHist.GetBinContent(i)==0 ){
        std::cout << "QCDHist bin bin on 0 content: min =  " <<  QCDHist.GetBinLowEdge(i) << " max = " << QCDHist.GetBinLowEdge(i+1) << std::endl;
      }
      if(DataHist.GetBinContent(i)==0 ){
        std::cout << "DataHist bin con 0 content: min =  " <<  DataHist.GetBinLowEdge(i) << " max = " << DataHist.GetBinLowEdge(i+1) << std::endl;
      }
      if(l1Tau_pt >= QCDHist.GetBinLowEdge(i) && l1Tau_pt<QCDHist.GetBinLowEdge(i+1)){
          //weight = (1./QCDHist.GetBinContent(i))*(n_tot_entries/(2.*QCDHist.GetNbinsX())); // to get an uniform
          weight = (DataHist.GetBinContent(i)/QCDHist.GetBinContent(i))*(QCDHist.Integral(1,QCDHist.GetNbinsX())/DataHist.Integral(1,DataHist.GetNbinsX())) ;
        }
      }
    }
  else{
    for (int i=1; i<= SignalHist.GetNbinsX(); ++i){
      if(SignalHist.GetBinContent(i)==0 ){
        std::cout << "bin con 0 content: min =  " <<  SignalHist.GetBinLowEdge(i) << " max = " << SignalHist.GetBinLowEdge(i+1) << std::endl;
      }
      if(genLepton_vis_pt >= SignalHist.GetBinLowEdge(i) && genLepton_vis_pt<SignalHist.GetBinLowEdge(i+1) ){
        weight = (1./SignalHist.GetBinContent(i))*(SignalHist.Integral(1,SignalHist.GetNbinsX())/(SignalHist.GetNbinsX()));
      }
    }
  }
  return weight;
}


void NewTuple::GetWeight(std::string dataFileName, std::string dataTupleName="L2TauTrainTuple"){
  gStyle->SetOptStat(111111111);

  /* open DataFrames */
  ROOT::RDataFrame dfNew("L2TauTrainTuple", newTupleFile); // MC
  ROOT::RDataFrame dfData(dataTupleName, absolute_path+dataFileName); // Data

  auto dfSignals = dfNew.Filter("genLepton_isTau==1");
  auto dfQCD = dfNew.Filter("genLepton_isTau==0");

  /* static binning, no more needed but I keep it for MC histogram */
  minimum_histogram_value = 15.;
  maximum_histogram_value = 275.;
  number_of_bin = 100;
  bin_width = (maximum_histogram_value-minimum_histogram_value)/static_cast<float>(number_of_bin) ;

  /* need to count number of total entries, it's the same for each flat variable */
  auto allDataSetHist = *dfNew.Histo1D({"allDataSetHist", "allDataSetHist", number_of_bin, minimum_histogram_value,maximum_histogram_value },"l1Tau_pt");
  int n_tot_entries = allDataSetHist.GetEntries();

  /* signal histogram*/
  //Float_t SignalBins[] = {15.002, 25.962, 36.922, 47.882, 58.842, 69.802, 80.762, 91.722, 102.682, 113.642, 124.602,  157.482, 1406.93};
  Float_t SignalBins[] = {15.002, 20.,27., 36., 47., 58., 69., 80.762, 91.722, 102.682, 113.642, 124.602,  157.482, 1044.38};
  // this binning was found by looking the signal histogram with static binning and unifying the empty/low stat bins
  Int_t  binnum = sizeof(SignalBins)/sizeof(Float_t) - 1;
  auto SignalToLookAt = dfSignals.Histo1D({"", "", binnum, SignalBins},"genLepton_vis_pt");
  LookAtHistogram(SignalToLookAt);
  auto SignalHist = *SignalToLookAt;

  /*QCD histogram*/
  Float_t BckgBins[] = {20.2, 24.5, 30. ,38.22, 47.882, 58.842, 69.802, 80.762, 91.722, 102.682, 113.642, 124.602,  157.482, 275.};
  // this binning was found by looking the bckg histogram with static binning and unifying the empty/low stat bins
  Int_t  binnum_b = sizeof(BckgBins)/sizeof(Float_t) - 1;
  auto QCDHistToLookAt = dfQCD.Histo1D({"", "",binnum_b, BckgBins}, "l1Tau_pt");
  std::cout << "Looking at QCD Histogram " << std::endl;
  LookAtHistogram(QCDHistToLookAt);
  auto QCDHist = *QCDHistToLookAt;

  /* Data histogram*/
  // this must have the same bin of QCD histogram
  auto DataHistToLookAt = dfData.Histo1D({"", "", binnum_b, BckgBins}, "l1Tau_pt");
  std::cout << "Looking at Data Histogram " << std::endl;
  LookAtHistogram(DataHistToLookAt);
  auto DataHist = *DataHistToLookAt;


  /* define lambda function to get weights from histograms */
  auto GetWeightFromHistos = [&](float &L1Pt, float& genLepton_vis_pt, bool &genLepton_isTau){
    return GetWeightFromHisto(L1Pt, genLepton_vis_pt, genLepton_isTau, SignalHist, QCDHist, DataHist);
  };


  /* define new dataset with weights */
  auto dfDataSetWeight = dfNew.Define("weight", GetWeightFromHistos, {"l1Tau_pt","genLepton_vis_pt","genLepton_isTau"});
  //dfDataSetWeight.Display("weight", 201)->Print();
  /* save dataset with only weights and with variables+weights*/
  dfDataSetWeight.Snapshot("L2TauTrainTuple", absolute_path+"DataSetTrainingOnlyWeight.root", {"weight"});
  dfDataSetWeight.Snapshot("L2TauTrainTuple", absolute_path+"DataSetTrainingWeight.root");

  /* let's print all stuff */
  std::cout << "all entries = " << n_tot_entries << std::endl;
  std::cout << "all entries/2 = " << static_cast<float>(n_tot_entries)/2. << std::endl;
  std::cout << "all weights' sum " << *dfDataSetWeight.Sum("weight")<<std::endl;
  std::cout << "sig weights' sum " << *dfDataSetWeight.Filter("genLepton_isTau==1").Sum("weight")<<std::endl;
  std::cout << "qcd weights' sum " << *dfDataSetWeight.Filter("genLepton_isTau==0").Sum("weight")<<std::endl;

  std::cout << "all weights " << std::endl;
  std::cout << dfDataSetWeight.Max("weight").GetValue() << std::endl;
  std::cout << dfDataSetWeight.Min("weight").GetValue() << std::endl;
  std::cout << "ratio all weights " << std::endl;
  std:: cout << dfDataSetWeight.Max("weight").GetValue() /dfDataSetWeight.Min("weight").GetValue() << std::endl;

  std::cout << "signal weights " << std::endl;
  std::cout << dfDataSetWeight.Filter("genLepton_isTau==1").Max("weight").GetValue() << std::endl;
  std::cout << dfDataSetWeight.Filter("genLepton_isTau==1").Min("weight").GetValue() << std::endl;
  std::cout << "ratio signal weights " << std::endl;
  std::cout << dfDataSetWeight.Filter("genLepton_isTau==1").Max("weight").GetValue()/ dfDataSetWeight.Filter("genLepton_isTau==1").Min("weight").GetValue() << std::endl;

  std::cout << "QCD weights " << std::endl;
  std::cout << dfDataSetWeight.Filter("genLepton_isTau==0").Max("weight").GetValue() << std::endl;
  std::cout << dfDataSetWeight.Filter("genLepton_isTau==0").Min("weight").GetValue() << std::endl;
  std::cout << "ratio qcd weights " << std::endl;
  std:: cout << dfDataSetWeight.Filter("genLepton_isTau==0").Max("weight").GetValue()/ dfDataSetWeight.Filter("genLepton_isTau==0").Min("weight").GetValue() << std::endl;
}


void NewTuple::DrawOnlyHistos(std::string dataFileName, std::string dataTupleName="L2TauTrainTuple"){
  gStyle->SetOptStat(111111111);

  /* open DataFrames */
  ROOT::RDataFrame dfNew("L2TauTrainTuple", newTupleFile); // MC
  ROOT::RDataFrame dfData(dataTupleName, absolute_path+dataFileName); // Data
  ROOT::RDataFrame dfDataSetWeight("L2TauTrainTuple", absolute_path+"DataSetTrainingWeight.root");


  auto dfSignals = dfNew.Filter("genLepton_isTau==1");
  auto dfQCD = dfNew.Filter("genLepton_isTau==0");

  /* static binning, no more needed but I keep it for MC histogram */
  minimum_histogram_value = 15.;
  maximum_histogram_value = 275.;
  number_of_bin = 100;
  bin_width = (maximum_histogram_value-minimum_histogram_value)/static_cast<float>(number_of_bin) ;

  /* need to count number of total entries, it's the same for each flat variable */
  auto allDataSetHist = *dfNew.Histo1D({"allDataSetHist", "allDataSetHist", number_of_bin, minimum_histogram_value,maximum_histogram_value },"l1Tau_pt");
  int n_tot_entries = allDataSetHist.GetEntries();

  /* signal histogram*/
  //Float_t SignalBins[] = {15.002, 25.962, 36.922, 47.882, 58.842, 69.802, 80.762, 91.722, 102.682, 113.642, 124.602,  157.482, 1406.93};
  Float_t SignalBins[] = {15.002, 20.,27., 36., 47., 58., 69., 80.762, 91.722, 102.682, 113.642, 124.602,  157.482, 1044.38};
  // this binning was found by looking the signal histogram with static binning and unifying the empty/low stat bins
  Int_t  binnum = sizeof(SignalBins)/sizeof(Float_t) - 1;
  auto SignalToLookAt = dfSignals.Histo1D({"", "", binnum, SignalBins},"genLepton_vis_pt");
  LookAtHistogram(SignalToLookAt);
  auto SignalHist = *SignalToLookAt;

  /*QCD histogram*/
  Float_t BckgBins[] = {20.2, 24.5, 30. ,38.22, 47.882, 58.842, 69.802, 80.762, 91.722, 102.682, 113.642, 124.602,  157.482, 275.};
  // this binning was found by looking the bckg histogram with static binning and unifying the empty/low stat bins
  Int_t  binnum_b = sizeof(BckgBins)/sizeof(Float_t) - 1;
  auto QCDHistToLookAt = dfQCD.Histo1D({"", "",binnum_b, BckgBins}, "l1Tau_pt");
  std::cout << "Looking at QCD Histogram " << std::endl;
  LookAtHistogram(QCDHistToLookAt);
  auto QCDHist = *QCDHistToLookAt;

  /* Data histogram*/
  // this must have the same bin of QCD histogram
  auto DataHistToLookAt = dfData.Histo1D({"", "", binnum_b, BckgBins}, "l1Tau_pt");
  std::cout << "Looking at Data Histogram " << std::endl;
  LookAtHistogram(DataHistToLookAt);
  auto DataHist = *DataHistToLookAt;

  /* Draw all histograms before reweighting */
  // 1. signal histogram
  TCanvas c1("c1", "c1", 10000,10000);
  SignalHist.GetXaxis()->SetTitle(("gen #tau p_{T} (GeV)"));
  SignalHist.GetYaxis()->SetTitle("A.U.");
  SignalHist.SetStats(0);
  //SignalHist.SetFillColorAlpha(kBlue, 0.35);
  SignalHist.SetLineColor(kBlue);
  SignalHist.Draw("HIST");
  c1.SaveAs((plotDir+"Reweighing/SignalHistBefore.pdf").c_str());

  // 2. qcd histogam
  TCanvas c2("c2", "c2", 10000,10000);
  QCDHist.GetXaxis()->SetTitle(("L1 #tau p_{T} (GeV)"));
  QCDHist.GetYaxis()->SetTitle("A.U.");
  QCDHist.SetStats(0);
  QCDHist.Scale(1/QCDHist.GetEntries());
  //QCDHist.SetFillColorAlpha(kBlue, 0.35);
  QCDHist.SetLineColor(kBlue);
  QCDHist.Draw("HIST");
  c2.SaveAs((plotDir+"Reweighing/QCDHistBefore.pdf").c_str());
  // 3. data histogram
  TCanvas c3("c3", "c3", 10000,10000);
  DataHist.GetXaxis()->SetTitle(("L1 #tau p_{T}"));
  DataHist.GetYaxis()->SetTitle("A.U.");
  DataHist.SetStats(0);
  DataHist.Scale(1/DataHist.GetEntries());
  //DataHist.SetFillColorAlpha(kRed, 0.35);
  DataHist.SetLineColor(kRed);
  DataHist.Draw("HIST");
  c3.SaveAs((plotDir+"Reweighing/DataHistBefore.pdf").c_str());
  /* draw superimposed data/qcd hists before reweighting*/
  TCanvas c3_1("c3_1", "c3_1", 10000,10000);
  TLegend *legend = new TLegend(0.74812, 0.820741, 0.901, 0.902 );
  QCDHist.SetStats(0);
  DataHist.SetStats(0);
  legend->AddEntry(&QCDHist, "QCD", "l");
  legend->AddEntry(&DataHist, "Data", "l");
  if(DataHist.GetBinContent(DataHist.GetMaximumBin()) > QCDHist.GetBinContent(QCDHist.GetMaximumBin()) ){
    DataHist.Draw("HIST");
    QCDHist.Draw("HIST SAME");
  }
  else{
    QCDHist.Draw("HIST");
    DataHist.Draw("HIST SAME");
  }
  legend->Draw("SAME");
  c3_1.SaveAs((plotDir+"Reweighing/QCDDataHistBefore.pdf").c_str());


  /* draw histograms after reweighting*/
  auto WeightedSignalHist = dfDataSetWeight.Filter("genLepton_isTau==1").Histo1D({"","",binnum, SignalBins}, "genLepton_vis_pt", "weight");
  auto WeightedQCDHist = dfDataSetWeight.Filter("genLepton_isTau==0").Histo1D({"","",binnum_b, BckgBins},"l1Tau_pt", "weight");


  // 1. signal
  auto hWeightedSignalHist = *WeightedSignalHist;
  auto hWeightedQCDHist = *WeightedQCDHist;
  TCanvas c4("c4", "c4", 10000,10000);
  hWeightedSignalHist.GetXaxis()->SetTitle(("gen #tau p_{T} (GeV)"));
  hWeightedSignalHist.GetYaxis()->SetTitle("A.U.");
  hWeightedSignalHist.SetStats(0);
  //hWeightedSignalHist.SetFillColorAlpha(kBlue, 0.35);
  hWeightedSignalHist.SetLineColor(kBlue);
  hWeightedSignalHist.Draw("EHIST");
  c4.SaveAs((plotDir+"Reweighing/WeightedSignalHist.pdf").c_str());

  // 2. background
  TCanvas c5("c5", "c5", 10000,10000);
  hWeightedQCDHist.GetXaxis()->SetTitle(("L1 #tau p_{T}"));
  hWeightedQCDHist.GetYaxis()->SetTitle("A.U.");
  hWeightedQCDHist.SetStats(0);
  //hWeightedQCDHist.SetFillColorAlpha(kBlue, 0.35);
  hWeightedQCDHist.Scale(1/hWeightedQCDHist.GetEntries());
  hWeightedQCDHist.SetLineColor(kBlue);
  hWeightedQCDHist.Draw("HIST");
  c5.SaveAs((plotDir+"Reweighing/WeightedQCDHist.pdf").c_str());

  //3. data
  TCanvas c6("c6", "c6", 10000,10000);
  TLegend *legend2 = new TLegend(0.74812, 0.820741, 0.901003, 0.902222 );
  hWeightedQCDHist.SetStats(0);
  hWeightedQCDHist.Scale(1/hWeightedQCDHist.GetEntries());
  DataHist.SetStats(0);
  legend2->AddEntry(&hWeightedQCDHist, "QCD", "l");
  legend2->AddEntry(&DataHist, "Data", "l");
  if(DataHist.GetBinContent(DataHist.GetMaximumBin()) > hWeightedQCDHist.GetBinContent(hWeightedQCDHist.GetMaximumBin()) ){
    DataHist.Draw("HIST");
    hWeightedQCDHist.Draw("HIST SAME");
  }
  else{
    hWeightedQCDHist.Draw("HIST");
    DataHist.Draw("HIST SAME");
  }
  legend2->Draw("SAME");

  c6.SaveAs((plotDir+"Reweighing/WeightedQCDDataHist.pdf").c_str());
}
