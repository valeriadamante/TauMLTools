#include "DataSetProducer.h"

DataSetProducer::DataSetProducer(){
}
DataSetProducer::DataSetProducer(const std::vector<std::string> _SignalFiles, const std::vector<std::string> _BackgroundFiles){
  absolute_path = "/Users/valeriadamante/Desktop/Dottorato/gridui/L2SkimmedTuples";
  SignalFiles=_SignalFiles;
  QCDFile=_BackgroundFiles;
}
DataSetProducer::DataSetProducer(const std::vector<std::string> _SignalFiles, const std::vector<std::string> _BackgroundFiles, const std::vector<std::string> _QCDFilteredFile){
  absolute_path = "/Users/valeriadamante/Desktop/Dottorato/gridui/L2SkimmedTuples";
  SignalFiles=_SignalFiles;
  QCDFile=_BackgroundFiles;
  QCDFilteredFile=_QCDFilteredFile;
}

void DataSetProducer::ProduceDataSample(std::vector<std::string> inputFile, std::string outFile_suffix){
  ROOT::RDataFrame dfToSnapshot("L2TauTrainTuple", inputFile);
  dfToSnapshot.Snapshot("L2TauTrainTuple", absolute_path+"DataSetTraining/all_"+outFile_suffix+".root");
}

bool DataSetProducer::FilterEventsFromRatio(std::vector<float>& ratios, std::vector<float>& low_edges, float L1Pt){
    /* 1. search L1pt interval */
    int index =0;
    for(std::vector<float, std::allocator<float> >::size_type i=0; i<low_edges.size(); i++){
      if(L1Pt>low_edges.at(i) && L1Pt<low_edges.at(i+1)){
        index = i ;
      }
    }
    /* 2. generate random number */
    if(gRandom->Uniform(0.,1.)<ratios.at(index)){
      return true;
    }
    return false;
}

std::vector<int> DataSetProducer::CountEvents(){
  std::vector<int>  total_events;
  for (auto& i: SignalFiles){
    ROOT::RDataFrame df("L2TauTrainTuple", i);
    total_events.push_back(df.Count().GetValue());
  }
  return total_events;
}


void DataSetProducer::QCDFilter(){
    ROOT::RDataFrame dfSignal("L2TauTrainTuple", SignalFiles);
    ROOT::RDataFrame dfBackground("L2TauTrainTuple", QCDFile);
    const float maximum_histogram_value = std::max(dfBackground.Max<float>("l1Tau_pt").GetValue(), dfSignal.Max<float>("l1Tau_pt").GetValue()) + 10.;
    const float minimum_histogram_value = 10.;
    const float number_of_bin = (maximum_histogram_value-minimum_histogram_value)/10.;
    auto signalHist = dfSignal.Histo1D<float>({"signal", "signal", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, "l1Tau_pt");
    auto backgroundHist = dfBackground.Histo1D<float>({"background", "background", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, "l1Tau_pt");

    std::vector<float> ratios;
    std::vector<float> low_edges;
    for(int i =0; i<number_of_bin+1 ; ++i){
        ratios.push_back(2*((*signalHist).GetBinContent(i)/(*backgroundHist).GetBinContent(i)));
        low_edges.push_back((*signalHist).GetBinLowEdge(i));
    }
    //std::cout << "Lower edge \t\t Ratio " << std::endl;
    for (std::vector<float>::size_type i =0; i<low_edges.size() ; i++ ){
      std::cout << "  " << low_edges.at(i) << "\t\t" << ratios.at(i)<< std::endl;
    }
    float ratio_min = 10.;
    for (auto&j : ratios){
      if(j<=ratio_min){
        ratio_min = j;
      }
    }
    auto filterEvents= [&](float L1Pt){
        return FilterEventsFromRatio(ratios, low_edges, L1Pt);
    };
    auto dfQCDFiltered = dfBackground.Define("randomNumber", " gRandom->Uniform(0.,1.) ")
                     .Filter(filterEvents,{"l1Tau_pt"});

    dfQCDFiltered.Snapshot("L2TauTrainTuple", absolute_path+"QCDFiltered.root");
    // auto displ=dfQCDFiltered.Display({"randomNumber","Filter","l1Tau_pt"}, 1000);
    // displ->Print();

}
void DataSetProducer::CountEventInSignalAndBackground(){
    ROOT::RDataFrame dfBackground("L2TauTrainTuple", QCDFile);
    ROOT::RDataFrame dfQCDFiltered("L2TauTrainTuple", QCDFilteredFile);
    int previous_bckg_events = dfBackground.Count().GetValue();
    int bckg_current_events = dfQCDFiltered.Count().GetValue();
    std::cout<< "QCD events before filtering "<< dfBackground.Count().GetValue()<<"\nQCD events after filtering "<<bckg_current_events <<std::endl;
    //std::cout <<"signal events for single sample: \t";
    std::vector<int> all_events = DataSetProducer::CountEvents();
    std::vector<std::string> all_files = SignalFiles;
    all_files.insert( all_files.end(), QCDFilteredFile.begin(), QCDFilteredFile.end() );
    all_events.push_back(bckg_current_events);
    int sum = 0;
    for (int i =0; i< all_events.size(); i++){
        std::cout << "sample = " << all_files[i] << "\t \t \t events "<< all_events[i]<<std::endl;
        sum += all_events[i];
    }

    // sum = 3303918;
    int sigEvts = sum-bckg_current_events;
    std::cout<<"Total signal events  = " << sigEvts <<std::endl;
    std::cout<<"Total events  = " << sum<<std::endl;
    std::cout<<"difference between signal & background  = " << std::abs(sigEvts - bckg_current_events)<<std::endl;
    std::cout<<"difference between signal & background relative to sig = " << float(float(std::abs(sigEvts - bckg_current_events))/float(sigEvts))<<std::endl;
    std::cout<<"difference between signal & background relative to bckg = " << float(float(std::abs(sigEvts - bckg_current_events))/float(bckg_current_events))<<std::endl;
    std::cout<<"difference between signal & background relative to all = " << float(float(std::abs(sigEvts - bckg_current_events))/float(sum))<<std::endl;
    std::vector<float> all_events_norm(all_events.size());
    float rescale_to_duecento = 200./float(sum);
    std::transform(all_events.begin(), all_events.end(), all_events_norm.begin(),
                   [&](int i) { return i * rescale_to_duecento; });
    for (int i =0; i< all_events_norm.size(); i++){
        std::cout << "sample = " << all_files[i] << "\t \t \tevents scaled to 200 "<< all_events_norm[i]<<std::endl;
    }

}
template <typename T>
void DataSetProducer::plot1(T h1, const std::string &VarName, const std::string& outDir_path, const std::string sigBckg)
{
   //gStyle->SetOptStat(0);
   gStyle->SetTextFont(42);
   auto c = new TCanvas("c", "", 800, 700);
   //TLegend *legend = new TLegend(0.421053, 0.82087 , 0.552632,0.902609 );
   //TLegend *legend = new TLegend( 0.62, 0.70, 0.82, 0.88 );
   TLegend *legend = new TLegend(0.74812, 0.820741, 0.901003, 0.902222 );
   c->SetLeftMargin(0.15);
   auto hist1 = *h1;
   hist1.SetStats(0);
   hist1.GetXaxis()->SetTitle(VarName.c_str());
   hist1.GetYaxis()->SetTitle("A.U.");
   hist1.SetLineWidth(2);
   hist1.SetLineColor(kBlue);
   //hist1.SetFillColor(kBlue);
   //hist1.SetFillStyle(3001);
   //hist1.SetMarkerStyle(20);
   //hist1.Rebin(2);
   hist1.SetMarkerColor(kBlue);
   hist1.Scale(100/hist1.GetEntries());
   legend->AddEntry(&hist1, hist1.GetTitle(), "l");
   hist1.DrawClone("HIST");
   legend->Draw("SAME");
   // Save plot
   c->SaveAs((outDir_path+VarName+"_"+sigBckg+".pdf").c_str());
}
template <typename T>
void DataSetProducer::plot2(T h1, T h2, const std::string &VarName, const std::string& outDir_path)
{
   // Canvas and general style options
   gStyle->SetOptStat(0);
   gStyle->SetTextFont(42);
   auto c = new TCanvas("c", "", 800, 700);
   //TLegend *legend = new TLegend(0.421053, 0.82087 , 0.552632,0.902609 );
   //TLegend *legend = new TLegend( 0.62, 0.70, 0.82, 0.88 );
   TLegend *legend = new TLegend(0.74812, 0.820741, 0.901003, 0.902222 );
   c->SetLeftMargin(0.15);

   auto hist1 = *h1;
   hist1.SetStats(0);
   hist1.GetXaxis()->SetTitle(VarName.c_str());
   hist1.GetYaxis()->SetTitle("A.U.");
   hist1.SetLineWidth(2);
   hist1.SetLineColor(kBlue);
   //hist1.SetFillColor(kBlue);
   //hist1.SetFillStyle(3001);
   //hist1.SetMarkerStyle(20);
   //hist1.Rebin(2);
   hist1.SetMarkerColor(kBlue);
   hist1.Scale(100/hist1.GetEntries());
   legend->AddEntry(&hist1, hist1.GetTitle(), "l");

   auto hist2 = *h2;
   hist2.SetStats(0);
   hist2.GetXaxis()->SetTitle(VarName.c_str());
   hist2.GetYaxis()->SetTitle("A.U.");
   hist2.SetLineWidth(2);
   hist2.SetLineColor(kRed);
   //hist2.SetFillColor(kRed);
   //hist2.SetFillStyle(3001);
   //hist2.SetMarkerStyle(20);
   //hist2.Rebin(2);
   hist2.SetMarkerColor(kRed);
   hist2.Scale(100/hist2.GetEntries());
   legend->AddEntry(&hist2, hist2.GetTitle(), "l");

   if(hist1.GetBinContent(hist1.GetMaximumBin()) > hist2.GetBinContent(hist2.GetMaximumBin()) ){
     hist1.DrawClone("HIST");
     hist2.DrawClone("HISTSAME");
   }
   else{
     hist2.DrawClone("HIST");
     hist1.DrawClone("HISTSAME");
   }

    legend->Draw("SAME");
   // Save plot
   c->SaveAs((outDir_path+VarName+".pdf").c_str());
}


template <typename T>
void DataSetProducer::plot3(T h1, T h2, T h3, const std::string &VarName, const std::string& outDir_path)
{
  // Canvas and general style options
  gStyle->SetOptStat(0);
  gStyle->SetTextFont(42);
  auto c = new TCanvas("c", "", 800, 700);
  //TLegend *legend = new TLegend(0.421053, 0.82087 , 0.552632,0.902609 );
  //TLegend *legend = new TLegend( 0.62, 0.70, 0.82, 0.88 );
  TLegend *legend = new TLegend(0.74812, 0.820741, 0.901003, 0.902222 );
  c->SetLeftMargin(0.15);

  auto hist1 = *h1;
  hist1.SetStats(0);
  hist1.GetXaxis()->SetTitle(VarName.c_str());
  hist1.GetYaxis()->SetTitle("A.U.");
  hist1.SetLineWidth(2);
  hist1.SetLineColor(kBlue);
  //hist1.SetFillColor(kBlue);
  //hist1.SetFillStyle(3001);
  //hist1.SetMarkerStyle(20);
  //hist1.Rebin(2);
  hist1.SetMarkerColor(kBlue);
  hist1.Scale(100/hist1.GetEntries());
  legend->AddEntry(&hist1, hist1.GetTitle(), "l");

  auto hist2 = *h2;
  hist2.SetStats(0);
  hist2.GetXaxis()->SetTitle(VarName.c_str());
  hist2.GetYaxis()->SetTitle("A.U.");
  hist2.SetLineWidth(2);
  hist2.SetLineColor(kRed);
  //hist2.SetFillColor(kRed);
  //hist2.SetFillStyle(3001);
  //hist2.SetMarkerStyle(20);
  //hist2.Rebin(2);
  hist2.SetMarkerColor(kRed);
  hist2.Scale(100/hist2.GetEntries());
  legend->AddEntry(&hist2, hist2.GetTitle(), "l");

  auto hist3 = *h3;
  hist3.SetStats(0);
  hist3.GetXaxis()->SetTitle(VarName.c_str());
  hist3.GetYaxis()->SetTitle("A.U.");
  hist3.SetLineWidth(2);
  hist3.SetLineColor(kGreen);
  //hist3.SetFillColor(kGreen);
  //hist3.SetFillStyle(3001);
  //hist3.SetMarkerStyle(20);
  //hist3.Rebin(2);
  hist3.SetMarkerColor(kGreen);
  hist3.Scale(100/hist3.GetEntries());
  legend->AddEntry(&hist3, hist3.GetTitle(), "l");

  if(hist1.GetBinContent(hist1.GetMaximumBin()) >= hist2.GetBinContent(hist2.GetMaximumBin()) && hist1.GetBinContent(hist1.GetMaximumBin()) >= hist3.GetBinContent(hist3.GetMaximumBin()) ){
    hist1.DrawClone("HIST");
    hist2.DrawClone("HISTSAME");
    hist3.DrawClone("HISTSAME");
  }
  else if (hist2.GetBinContent(hist2.GetMaximumBin()) > hist1.GetBinContent(hist1.GetMaximumBin()) && hist2.GetBinContent(hist2.GetMaximumBin()) > hist3.GetBinContent(hist3.GetMaximumBin()) ){
    hist2.DrawClone("HIST");
    hist1.DrawClone("HISTSAME");
    hist3.DrawClone("HISTSAME");
  }
  else{
    hist3.DrawClone("HIST");
    hist1.DrawClone("HISTSAME");
    hist2.DrawClone("HISTSAME");
  }

   legend->Draw("SAME");
  // Save plot
  c->SaveAs((outDir_path+VarName+".pdf").c_str());
}

void DataSetProducer::GetHistogramsSignalQCD(int n_var, bool use_binning){
  ROOT::RDataFrame dfSignal("L2TauTrainTuple", SignalFiles);
  ROOT::RDataFrame dfQCD("L2TauTrainTuple", QCDFile);
  ROOT::RDataFrame dfQCDFiltered("L2TauTrainTuple", QCDFilteredFile);
  all_branches_binning = {{"genEventWeight",{50, -0.5,1.5}},{"genLepton_vis_pt",{100, 0, 1000}},{"genLepton_vis_eta",{100, -4.5, 5.5}},{"genLepton_vis_phi",{70, -3.5, 3.5}},{"genLepton_vis_mass",{100, -0.1, 9.9}},{"l1Tau_pt",{50, 0, 250}},{"l1Tau_eta",{100, -4.5, 5.5}},{"l1Tau_phi",{70, -3.5, 3.5}},{"l1Tau_mass",  {100, -0,1, 9.9}},{"caloRecHit_ee_rho",{100, 0, 300}},{"caloRecHit_eb_rho",{100, 0, 300}},{"caloRecHit_ee_eta",{100, -4.5, 5.5}},{"caloRecHit_eb_eta",{100, -4.5, 5.5}},{"caloRecHit_ee_phi",{70, -3.5, 3.5}},{"caloRecHit_eb_phi",{70, -3.5, 3.5}},{"caloRecHit_ee_energy",{60, 0, 6}},{"caloRecHit_eb_energy",{100, 0, 80}},{"caloRecHit_ee_time",{50, -0.5, 9.5}},{"caloRecHit_eb_time",{50, -0.5, 9.5}},{"caloRecHit_ee_chi2",{100, -0.5, 60.5}},{"caloRecHit_eb_chi2",{100, -0.5, 60.5}},{"caloRecHit_ee_energyError",{50, -0.5, 9.5}},{"caloRecHit_eb_energyError",{50, -0.5, 9.5}},{"caloRecHit_ee_timeError",{50, -0.5, 9.5}},{"caloRecHit_eb_timeError",{50, -0.5, 9.5}},{"caloRecHit_hbhe_rho",{100, 0, 300}},{"caloRecHit_hbhe_eta",{100, -4.5, 5.5}},{"caloRecHit_hbhe_phi",{70, -3.5, 3.5}},{"caloRecHit_hbhe_energy",{100, 0, 80}},{"caloRecHit_hbhe_time",{50, -0.5, 9.5}},{"caloRecHit_hbhe_chi2",{100, -0.5, 60.5}},{"caloRecHit_ho_rho",{100, 0, 300}},{"caloRecHit_ho_eta",{100, -4.5, 5.5}},{"caloRecHit_ho_phi",{70, -3.5, 3.5}},{"caloRecHit_ho_energy",{100, 0, 80}},{"caloRecHit_ho_time",{50, -0.5, 9.5}},{"caloRecHit_hf_rho",{100, 0, 300}},{"caloRecHit_hf_eta",{100, -4.5, 5.5}},{"caloRecHit_hf_phi",{70, -3.5, 3.5}},{"caloRecHit_hf_energy",{100, 0, 80}},{"caloRecHit_hf_time",{50, -0.5, 9.5}},{"caloRecHit_hf_timeFalling",{50, -0.5, 9.5}},{"patatrack_pt",{50, -0.5, 20.5}},{"patatrack_eta",{100, -4.5, 5.5}},{"patatrack_phi",{70, -3.5, 3.5}},{"patatrack_chi2",{100, -0.5, 60.5}},{"patatrack_dxy",{60, -0.3, 0.3}},{"patatrack_dz",{300, -15.5, 15.5}},{"patavert_z",{300, -15.5, 15.5}},{"patavert_weight",{100, 0, 2000}},{"patavert_ptv2",{100}},{"patavert_chi2",{100, -0.5, 1000.5}},{"caloRecHit_ee_detId",{1000}},{"caloRecHit_eb_flagsBits",{100}},{"caloRecHit_ee_isRecovered",{100}},{"caloRecHit_eb_isRecovered",{100}},{"caloRecHit_ee_isTimeValid",{100}},{"caloRecHit_eb_isTimeValid",{100}},{"caloRecHit_ee_isTimeErrorValid",{100}},{"caloRecHit_eb_isTimeErrorValid",{100}},{"caloRecHit_hbhe_detId",{100}},{"caloRecHit_hbhe_flags",{100}},{"caloRecHit_ho_detId",{100}},{"caloRecHit_ho_aux",{100}},{"caloRecHit_hf_flags",{100}},{"caloRecHit_hf_auxHF",{100}},{"patatrack_charge",{50, -0.5, 1.5}},{"patatrack_quality",{50, -0.5, 1.5}},{"patatrack_vertex_id",{50, 0, 100}},{"patavert_ndof",{100, 0, 200}},{"run",{10, -0.5,1.5}},{"lumi",{1000}}, {"defaultDiTauPath_result",{10}}, {"defaultDiTauPath_lastModuleIndex",{1000}},{"evt",{1000}},{"sampleType",{10, -0.5,1.5}},{"genLepton_index",{60, -0.1, 5,9}},{"genLepton_kind",{100, -0.5, 9.5}},{"genLepton_charge",{10, -0.5, 1.5}},{"caloRecHit_eb_detId",{1000}},{"caloRecHit_ee_flagsBits",{100}},{"caloRecHit_hf_detId",{100}},{"caloRecHit_hf_aux",{100}},{"patatrack_ndof",{10, 0, 10}},{"l1Tau_hwIso",{10, -0.5,1.5}},{"l1Tau_hwQual",{10, -0.5,1.5}},{"l1Tau_towerIEta",{10, -0.5,1.5}},{"l1Tau_towerIPhi",{10, -0.5,1.5}},{"l1Tau_rawEt",{10, -0.5,1.5}},{"l1Tau_isoEt",{10, -0.5,1.5}},{"l1Tau_hasEM",{10, -0.5,1.5}},{"l1Tau_isMerged" ,{10, -0.5,1.5}},{"l1Tau_index",{50, 0, 10}},{"caloRecHit_ho_flags",{100}} };
  bool found = false;
  if (all_branches.find(n_var) != all_branches.end()){
      found=true;
  }
  if(found==false) {return;}
  std::string variabile = all_branches.at(n_var);
  std::cout << "************************************************************" << std::endl ;
  std::cout << "    plotted observable  = " << variabile << std::endl;
  std::cout << "************************************************************" << std::endl ;
  auto maximum_histogram_value = std::max(dfQCD.Max(variabile).GetValue(), dfSignal.Max(variabile).GetValue())*(1.5);
  auto minimum_histogram_value = std::min(dfQCD.Min(variabile).GetValue(), dfSignal.Min(variabile).GetValue());
  auto bins = maximum_histogram_value-minimum_histogram_value;
  int number_of_bin = static_cast<int>(bins);
  std::cout << "previously minimum_histogram_value was " << minimum_histogram_value <<"\npreviously maximum_histogram_value was " << maximum_histogram_value << "\npreviously number of bins was " << number_of_bin <<std::endl;
  std::cout << std::endl ;
  if(minimum_histogram_value<0){
    minimum_histogram_value = minimum_histogram_value*1.5;
    std::cout << "after rescaling maximum_histogram_value is " << maximum_histogram_value <<std::endl;
  }
  else{
    minimum_histogram_value = minimum_histogram_value*0.75;
    std::cout << "after rescaling minimum_histogram_value is " << minimum_histogram_value <<std::endl;
  }
  if(minimum_histogram_value == maximum_histogram_value){
    maximum_histogram_value = minimum_histogram_value+10.;
    std::cout << "min = max, so we changed them into \nminimum =  " << minimum_histogram_value << "\nmaximum = " << maximum_histogram_value <<std::endl;
  }
  if (number_of_bin == 0){
    number_of_bin = 10;
  }
  else if (number_of_bin> 0 && number_of_bin <= 10){
    number_of_bin = number_of_bin*10;
  }
  else{
    number_of_bin = 100;
  }
  auto signalHistBef = dfSignal.Histo1D({"signal bef", "signal bef", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, variabile);
  auto QCDFilteredHistBef = dfQCDFiltered.Histo1D({"QCDFiltered bef", "QCDFiltered bef", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, variabile);
  auto QCDHistBef = dfQCD.Histo1D({"QCD bef", "QCD bef", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, variabile);
  std::string outDir_path1 = "/Users/valeriadamante/Desktop/Dottorato/plots/outVars/QCDSignalHists/Histo1D_bef/";
  DataSetProducer::plot1(signalHistBef, variabile,outDir_path1, "sig");
  DataSetProducer::plot1(QCDHistBef, variabile,outDir_path1, "QCD");
  DataSetProducer::plot1(QCDFilteredHistBef, variabile,outDir_path1, "QCDFiltered");
  std::cout << "after rescaling number of bin was " << number_of_bin <<std::endl;
  std::cout << std::endl ;
  bool foundBinningInMap = false;
  if (all_branches.find(n_var) != all_branches.end()){
      foundBinningInMap=true;
  }
  if(use_binning && foundBinningInMap&& all_branches_binning.at(variabile).size()>1){
    number_of_bin = static_cast<int>(all_branches_binning.at(variabile).at(0));
    minimum_histogram_value = all_branches_binning.at(variabile).at(1);
    maximum_histogram_value = all_branches_binning.at(variabile).at(2);
    std::cout << "number of bins = " << number_of_bin << std::endl;
    std::cout << "minimum_histogram_value = " << minimum_histogram_value << std::endl;
    std::cout << "maximum_histogram_value = " << maximum_histogram_value << std::endl;
  }


  auto signalHist = dfSignal.Histo1D({"signal", "signal", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, variabile);
  auto QCDFilteredHist = dfQCDFiltered.Histo1D({"QCDFiltered", "QCDFiltered", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, variabile);
  auto QCDHist = dfQCD.Histo1D({"QCD", "QCD", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, variabile);
  std::string outDir_path11 = "/Users/valeriadamante/Desktop/Dottorato/plots/outVars/QCDSignalHists/Histo1D/";
  DataSetProducer::plot1(signalHist, variabile,outDir_path11, "sig");
  DataSetProducer::plot1(QCDHist, variabile,outDir_path11, "QCD");
  DataSetProducer::plot1(QCDFilteredHist, variabile,outDir_path11,"QCDFiltered");

  std::string outDir_path2 = "/Users/valeriadamante/Desktop/Dottorato/plots/outVars/QCDSignalHists/SigVsQCD/";
  DataSetProducer::plot2(signalHist, QCDHist, variabile,outDir_path2);
  std::string outDir_path3 = "/Users/valeriadamante/Desktop/Dottorato/plots/outVars/QCDSignalHists/SigVsQCDVsQCDFiltered/";
  DataSetProducer::plot3(signalHist, QCDHist,  QCDFilteredHist, variabile,outDir_path3);

}
void DataSetProducer::GetHistogramsDataValidation(int n_var, bool use_binning){
  ROOT::RDataFrame dfVBF("L2TauTrainTuple", VBFFile);
  ROOT::RDataFrame dfData("L2TauTrainTuple", DataFile);
  ROOT::RDataFrame dfZP("L2TauTrainTuple", ZPFile);
  all_branches_binning = {{"genEventWeight",{50, -0.5,1.5}},{"genLepton_vis_pt",{100, 0, 1000}},{"genLepton_vis_eta",{100, -4.5, 5.5}},{"genLepton_vis_phi",{70, -3.5, 3.5}},{"genLepton_vis_mass",{100, -0.1, 9.9}},{"l1Tau_pt",{50, 0, 250}},{"l1Tau_eta",{100, -4.5, 5.5}},{"l1Tau_phi",{70, -3.5, 3.5}},{"l1Tau_mass",  {100, -0,1, 9.9}},{"caloRecHit_ee_rho",{100, 0, 300}},{"caloRecHit_eb_rho",{100, 0, 300}},{"caloRecHit_ee_eta",{100, -4.5, 5.5}},{"caloRecHit_eb_eta",{100, -4.5, 5.5}},{"caloRecHit_ee_phi",{70, -3.5, 3.5}},{"caloRecHit_eb_phi",{70, -3.5, 3.5}},{"caloRecHit_ee_energy",{60, 0, 6}},{"caloRecHit_eb_energy",{100, 0, 80}},{"caloRecHit_ee_time",{50, -0.5, 9.5}},{"caloRecHit_eb_time",{50, -0.5, 9.5}},{"caloRecHit_ee_chi2",{100, -0.5, 60.5}},{"caloRecHit_eb_chi2",{100, -0.5, 60.5}},{"caloRecHit_ee_energyError",{50, -0.5, 9.5}},{"caloRecHit_eb_energyError",{50, -0.5, 9.5}},{"caloRecHit_ee_timeError",{50, -0.5, 9.5}},{"caloRecHit_eb_timeError",{50, -0.5, 9.5}},{"caloRecHit_hbhe_rho",{100, 0, 300}},{"caloRecHit_hbhe_eta",{100, -4.5, 5.5}},{"caloRecHit_hbhe_phi",{70, -3.5, 3.5}},{"caloRecHit_hbhe_energy",{100, 0, 80}},{"caloRecHit_hbhe_time",{50, -0.5, 9.5}},{"caloRecHit_hbhe_chi2",{100, -0.5, 60.5}},{"caloRecHit_ho_rho",{100, 0, 300}},{"caloRecHit_ho_eta",{100, -4.5, 5.5}},{"caloRecHit_ho_phi",{70, -3.5, 3.5}},{"caloRecHit_ho_energy",{100, 0, 80}},{"caloRecHit_ho_time",{50, -0.5, 9.5}},{"caloRecHit_hf_rho",{100, 0, 300}},{"caloRecHit_hf_eta",{100, -4.5, 5.5}},{"caloRecHit_hf_phi",{70, -3.5, 3.5}},{"caloRecHit_hf_energy",{100, 0, 80}},{"caloRecHit_hf_time",{50, -0.5, 9.5}},{"caloRecHit_hf_timeFalling",{50, -0.5, 9.5}},{"patatrack_pt",{50, -0.5, 20.5}},{"patatrack_eta",{100, -4.5, 5.5}},{"patatrack_phi",{70, -3.5, 3.5}},{"patatrack_chi2",{100, -0.5, 60.5}},{"patatrack_dxy",{60, -0.3, 0.3}},{"patatrack_dz",{300, -15.5, 15.5}},{"patavert_z",{300, -15.5, 15.5}},{"patavert_weight",{100, 0, 2000}},{"patavert_ptv2",{100}},{"patavert_chi2",{100, -0.5, 1000.5}},{"caloRecHit_ee_detId",{1000}},{"caloRecHit_eb_flagsBits",{100}},{"caloRecHit_ee_isRecovered",{100}},{"caloRecHit_eb_isRecovered",{100}},{"caloRecHit_ee_isTimeValid",{100}},{"caloRecHit_eb_isTimeValid",{100}},{"caloRecHit_ee_isTimeErrorValid",{100}},{"caloRecHit_eb_isTimeErrorValid",{100}},{"caloRecHit_hbhe_detId",{100}},{"caloRecHit_hbhe_flags",{100}},{"caloRecHit_ho_detId",{100}},{"caloRecHit_ho_aux",{100}},{"caloRecHit_hf_flags",{100}},{"caloRecHit_hf_auxHF",{100}},{"patatrack_charge",{50, -0.5, 1.5}},{"patatrack_quality",{50, -0.5, 1.5}},{"patatrack_vertex_id",{50, 0, 100}},{"patavert_ndof",{100, 0, 200}},{"run",{10, -0.5,1.5}},{"lumi",{1000}}, {"defaultDiTauPath_result",{10}}, {"defaultDiTauPath_lastModuleIndex",{1000}},{"evt",{1000}},{"sampleType",{10, -0.5,1.5}},{"genLepton_index",{60, -0.1, 5,9}},{"genLepton_kind",{100, -0.5, 9.5}},{"genLepton_charge",{10, -0.5, 1.5}},{"caloRecHit_eb_detId",{1000}},{"caloRecHit_ee_flagsBits",{100}},{"caloRecHit_hf_detId",{100}},{"caloRecHit_hf_aux",{100}},{"patatrack_ndof",{10, 0, 10}},{"l1Tau_hwIso",{10, -0.5,1.5}},{"l1Tau_hwQual",{10, -0.5,1.5}},{"l1Tau_towerIEta",{10, -0.5,1.5}},{"l1Tau_towerIPhi",{10, -0.5,1.5}},{"l1Tau_rawEt",{10, -0.5,1.5}},{"l1Tau_isoEt",{10, -0.5,1.5}},{"l1Tau_hasEM",{10, -0.5,1.5}},{"l1Tau_isMerged" ,{10, -0.5,1.5}},{"l1Tau_index",{50, 0, 10}},{"caloRecHit_ho_flags",{100}} };
  bool found = false;
  if (all_branches.find(n_var) != all_branches.end()){
      found=true;
  }
  if(found==false) return;
  std::string variabile = all_branches.at(n_var);
  std::cout << "************************************************************" << std::endl ;
  std::cout << "     plotted observable  = " << variabile << std::endl;
  std::cout << "************************************************************" << std::endl ;
  auto maximum_histogram_value = std::max(dfZP.Max(variabile).GetValue(), dfVBF.Max(variabile).GetValue())*(1.5);
  auto minimum_histogram_value = std::min(dfZP.Min(variabile).GetValue(), dfVBF.Min(variabile).GetValue());
  auto bins = maximum_histogram_value-minimum_histogram_value;
  int number_of_bin = static_cast<int>(bins);
  std::cout << "previously minimum_histogram_value was " << minimum_histogram_value <<"\npreviously maximum_histogram_value was " << maximum_histogram_value << "\npreviously number of bins was " << number_of_bin <<std::endl;
  std::cout << std::endl ;
  if(minimum_histogram_value<0){
    minimum_histogram_value = minimum_histogram_value*1.5;
    std::cout << "after rescaling maximum_histogram_value is " << maximum_histogram_value <<std::endl;
  }
  else{
    minimum_histogram_value = minimum_histogram_value*0.75;
    std::cout << "after rescaling minimum_histogram_value is " << minimum_histogram_value <<std::endl;
  }
  if(minimum_histogram_value == maximum_histogram_value){
    maximum_histogram_value = minimum_histogram_value+10.;
    std::cout << "min = max, so we changed them into \nminimum =  " << minimum_histogram_value << "\nmaximum = " << maximum_histogram_value <<std::endl;
  }
  if (number_of_bin == 0){
    number_of_bin = 10;
  }
  else if (number_of_bin> 0 && number_of_bin <= 10){
    number_of_bin = number_of_bin*10;
  }
  else{
    number_of_bin = 100;
  }
  std::cout << "after rescaling number of bin was " << number_of_bin <<std::endl;
  std::cout << std::endl ;
  if(use_binning && all_branches_binning.at(variabile).size()>1){
    number_of_bin = static_cast<int>(all_branches_binning.at(variabile).at(0));
    minimum_histogram_value = all_branches_binning.at(variabile).at(1);
    maximum_histogram_value = all_branches_binning.at(variabile).at(2);
    std::cout << "number of bins = " << number_of_bin << std::endl;
    std::cout << "minimum_histogram_value = " << minimum_histogram_value << std::endl;
    std::cout << "maximum_histogram_value = " << maximum_histogram_value << std::endl;
  }
  auto VBFHist = dfVBF.Histo1D({"VBF", "VBF", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, variabile);
  auto ZPHist = dfZP.Histo1D({"ZPrime", "ZPrime", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, variabile);
  auto DataHist = dfData.Histo1D({"Data", "Data", static_cast<int>(number_of_bin), minimum_histogram_value, maximum_histogram_value}, variabile);
  std::string outDir_path1 = "/Users/valeriadamante/Desktop/Dottorato/plots/outVars/VBFZPDataHists/VBFZPDataHists1D/";
  DataSetProducer::plot1(DataHist, variabile,outDir_path1, "data");
  DataSetProducer::plot1(VBFHist, variabile,outDir_path1, "VBF");
  DataSetProducer::plot1(ZPHist, variabile,outDir_path1, "ZP");
  std::string outDir_path2 = "/Users/valeriadamante/Desktop/Dottorato/plots/outVars/VBFZPDataHists/VBFZPDataHists2D/";
  DataSetProducer::plot2(VBFHist, ZPHist, variabile,outDir_path2);
}
