#include <iostream>
#include<cmath>
#include <vector>
#include <iterator>
#include <algorithm> // for std::transform
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "ROOT/RDataFrame.hxx"
#include "TRandom3.h"
#include "TLegend.h"

//std::string absolute_path = "/home/users/damante/L2SkimmedTuples/";
std::string absolute_path = "/Users/valeriadamante/Desktop/Dottorato/gridui/L2SkimmedTuples/";
std::string mid_path_ZP = "ZprimeToTauTau_M-4000/Histo_ZprimeToTauTau_M-4000_TuneCP5_14TeV-pythia8-tauola_output_";
std::string mid_path_VBF = "VBFHToTauTau/Histo_VBFHToTauTau_M125_TuneCUETP8M1_14TeV_powheg_pythia8_output_";
std::string end_path = ".root";
TH1F * Histo_Merger(std::vector<std::string>& files){
  TH1F *hsum = nullptr;
  for(std::vector<std::string>::size_type i =0 ; i< files.size(); i++){
    TFile *f = new TFile(files.at(i).c_str(), "READ");
    TH1F *histogram = (TH1F*) f->Get("Epi");
    if(i==0){
      hsum =(TH1F*)histogram->Clone("hsum");
    }
    else {
      hsum->Add(histogram);
    }
  }
  return hsum;
}


int Merge_Histograms(){
    const std::vector<std::string> file_numbers_ZP  = {"1-1","1-2","1-3","1-5","1-6","1-7","1","3","6"};
    const std::vector<std::string> file_numbers_VBF = {"10","11","12","13","14","15","16","17","18","1","2","3","4","5","6","7","8","9"};
    std::vector<std::string> files_ZP, files_VBF;
    for(auto& i: file_numbers_ZP){
      files_ZP.push_back(absolute_path + mid_path_ZP + i + end_path);
    }
    for(auto& i: file_numbers_VBF){
      files_VBF.push_back(absolute_path + mid_path_VBF + i + end_path);
    }
    TH1F *hsum_ZP = Histo_Merger(files_ZP);
    TH1F *hsum_VBF = Histo_Merger(files_VBF);
    TCanvas *c = new TCanvas ("","", 800,600);
    TLegend *legend = new TLegend(0.421053, 0.82087 , 0.552632,0.902609 );

    //auto legend = new TLegend(0.1,0.7,0.48,0.9);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    hsum_ZP->SetTitle("hSumZP");
    hsum_ZP->SetStats(0);
    hsum_ZP->Rebin(2);
    hsum_ZP->GetXaxis()->SetTitle("#frac{E_{#pi_{+}}-E_{#pi_{0}}}{E_{#pi_{+}}+E_{#pi_{0}}}");
    hsum_ZP->GetYaxis()->SetTitle("A.U.");
    hsum_ZP->SetLineColor(kBlue);
    hsum_ZP->SetFillColor(kBlue);
    //hsum_ZP->SetMarkerStyle(20);
    hsum_ZP->SetMarkerColor(kBlue);
    hsum_ZP->Scale(100/hsum_ZP->GetEntries());
    legend->AddEntry(hsum_ZP, "ZPrime", "lep");

    hsum_VBF->SetTitle("hSumVBF");
    hsum_VBF->SetStats(0);
    hsum_VBF->GetXaxis()->SetTitle("(E_{#pi^{+}}-E_{#pi^{0}})/(E_{#pi^{+}}+E_{#pi^{0}})");
    hsum_VBF->GetYaxis()->SetTitle("A.U.");
    hsum_VBF->SetFillColor(kRed);
    hsum_VBF->SetLineColor(kRed);
    hsum_VBF->Rebin(2);
    //hsum_VBF->SetMarkerStyle(20);
    hsum_VBF->SetMarkerColor(kRed);
    hsum_VBF->Scale(100/hsum_VBF->GetEntries());
    if(hsum_VBF->GetBinContent(hsum_VBF->GetMaximumBin())>hsum_ZP->GetBinContent(hsum_ZP->GetMaximumBin())){
      hsum_VBF->Draw("E");
      hsum_ZP->Draw("ESAME");
    }
    else{

      hsum_ZP->Draw("E");
      hsum_VBF->Draw("ESAME");
    }

    legend->AddEntry(hsum_VBF, "VBF","lep");
    legend->Draw("SAME");
    c->Update();
    c->Print("SumHistogram.pdf", "pdf");
    c->Print("SumHistogram.png", "png");
    return 0;
}
