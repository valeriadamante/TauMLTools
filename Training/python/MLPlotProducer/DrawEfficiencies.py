# algo efficiency for test dataset to compare with sqrt(diagonal) of event cutbased efficiency
# how to build this algorithm:

# 1. goal = compare algo efficiency per taus (NOT per events) on testing dataset with sqrt(diagonal) of cut based (algo?) efficiency (on VBF) per event (NOT per taus) in pt bins
import statsmodels.stats.proportion as ssp
import numpy as np
from array import array
import os
import ROOT
import math
full_name_evtTuple = "/Users/valeriadamante/Desktop/Dottorato/cmssimphase2/outputs_EvtTuples/correct_efficiencies/efficiency_algo_hltDoubleL2IsoTau26eta2p2.root"
full_name_tauTuple = "/Users/valeriadamante/Desktop/Dottorato/cmssimphase2/outputs_TauTuples/Efficiencies.root"
full_name_evtTuple_BOR = "/Users/valeriadamante/Desktop/Dottorato/cmssimphase2/outputs_EvtTuples/correct_efficiencies/efficiency_normal_hltL1sDoubleTauBigOR.root"

ROOT.gStyle.SetPaintTextFormat(".2f")

def beautify_plot_canvas(histNum_cb, histDen_cb, histNum_dnn, histDen_dnn, names, save=True):
    H_ref = 1000;
    W_ref = 1000;
    W = W_ref
    H  = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref
    L = 0.12*W_ref
    R = 0.1*W_ref
    pt_bins = [20,30,35,40,50,60,70,80,90,350]
    k=0
    n_bins = histDen_cb.GetNbinsX()
    gr_cb= ROOT.TGraphAsymmErrors(n_bins)
    gr_dnn= ROOT.TGraphAsymmErrors(n_bins)

    for i in range(1,n_bins+1):
        num_cb = histNum_cb.GetBinContent(i)
        den_cb = histDen_cb.GetBinContent(i)
        num_dnn = histNum_dnn.GetBinContent(i)
        den_dnn = histDen_dnn.GetBinContent(i)
        eff_cb = math.sqrt(num_cb/den_cb)
        eff_dnn = math.sqrt(num_dnn/den_dnn)
        x_p = pt_bins[i] - ((pt_bins[i]-pt_bins[k])/2)
        print(k)
        gr_cb.SetPointX(k,x_p)
        gr_dnn.SetPointX(k,x_p)
        gr_cb.SetPointY(k,eff_cb)
        gr_dnn.SetPointY(k,eff_dnn)
        print(x_p)
        print(eff_cb)
        err_x_low = (pt_bins[i]-pt_bins[k])/2
        err_x_up= (pt_bins[i]-pt_bins[k])/2
        c_low_cb, c_up_cb = ssp.proportion_confint(num_cb, den_cb, alpha=1-0.68, method='beta')
        c_low_dnn, c_up_dnn = ssp.proportion_confint(num_dnn, den_dnn, alpha=1-0.68, method='beta')
        err_y_low_cb = eff_cb-math.sqrt(c_low_cb)
        err_y_up_cb = math.sqrt(c_up_cb) - eff_cb
        err_y_low_dnn = eff_dnn-math.sqrt(c_low_dnn)
        err_y_up_dnn = math.sqrt(c_up_dnn) - eff_dnn
        gr_cb.SetPointError(k, err_x_low, err_x_up, err_y_low_cb, err_y_up_cb)
        gr_dnn.SetPointError(k, err_x_low, err_x_up, err_y_low_dnn, err_y_up_dnn)
        k+=1

    canvas=ROOT.TCanvas()
    canvas.cd()
    legend = ROOT.TLegend(0.75,0.6,0.88,0.7)
    gr_cb.SetMarkerColor( ROOT.kRed )
    gr_cb.SetLineColor( ROOT.kRed )
    gr_cb.SetMarkerStyle( 20 )
    gr_cb.SetTitle(("{};{};{}").format(names["histTitle"],names["xAxisTitle"],names["yAxisTitle"]))
    legend.AddEntry(gr_cb, "cut based", "lep")
    #legend.SetBorderSize(2)
    gr_dnn.SetMarkerColor( ROOT.kBlue )
    gr_dnn.SetLineColor( ROOT.kBlue )
    gr_dnn.SetMarkerStyle( 20 )
    gr_dnn.SetTitle(("{};{};{}").format(names["histTitle"],names["xAxisTitle"],names["yAxisTitle"]))
    legend.AddEntry(gr_dnn, "dnn based", "lep")
    gr_dnn.Draw("AP")
    gr_cb.Draw("PSAME")
    legend.Draw("SAME")
    canvas.Update()
    if(save):
        canvas.Update()
        #canvas.Print(("{}.png").format(names["outFile"]),"png")
        canvas.Print(("{}.pdf").format(names["outFile"]),"pdf")
    #input()


def beautify_plot_canvas_1hist(histNum, histDen, names, save=True):
    H_ref = 1000;
    W_ref = 1000;
    W = W_ref
    H  = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref
    L = 0.12*W_ref
    R = 0.1*W_ref
    pt_bins = [20,30,35,40,50,60,70,80,90,350]
    k=0
    gr= ROOT.TGraphAsymmErrors(len(pt_bins))
    n_bins = histDen.GetNbinsX()

    for i in range(1,n_bins+1):
        num = histNum.GetBinContent(i)
        den = histDen.GetBinContent(i)
        eff = math.sqrt(num/den)
        gr.SetPointY(k,eff)
        x_p = pt_bins[i] - ((pt_bins[i]-pt_bins[k])/2)
        gr.SetPointX(k,x_p)
        err_x_low = (pt_bins[i]-pt_bins[k])/2
        err_x_up= (pt_bins[i]-pt_bins[k])/2
        c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
        err_y_low = eff-math.sqrt(c_low)
        err_y_up = math.sqrt(c_up) - eff
        gr.SetPointError(k, err_x_low, err_x_up, err_y_low, err_y_up)
        k+=1
    canvas=ROOT.TCanvas()
    canvas.cd()
    gr.SetMarkerColor( ROOT.kRed )
    gr.SetLineColor( ROOT.kRed )
    gr.SetMarkerStyle( 20 )
    gr.SetTitle(("{};{};{}").format(names["histTitle"],names["xAxisTitle"],names["yAxisTitle"]))
    gr.Draw("AP")
    canvas.Update()
    if(save):
        canvas.Update()
        canvas.Print(("{}.png").format(names["outFile"]),"png")
        canvas.Print(("{}.pdf").format(names["outFile"]),"pdf")


absolute_path ="/Users/valeriadamante/Desktop/Dottorato/cmssimphase2/"
evtTuple_dir= absolute_path+"outputs_EvtTuples/correct_efficiencies/"
tauTuple_dir= absolute_path+"outputs_TauTuples/"
all_files=[evtTuple_dir+"efficiency_algo_hltDoubleL2IsoTau26eta2p2.root",evtTuple_dir+"efficiency_normal_hltDoubleL2IsoTau26eta2p2.root",evtTuple_dir+"efficiency_normal_hltL1sDoubleTauBigOR.root", tauTuple_dir+"Efficiencies.root"]
numHist_names = ["passed_algo", "passed_normal", "passed_normal", "passed_d"]
denHist_names = ["total_algo", "total_normal", "total_normal", "total_d"]
x_y_titles = ["gen #tau_{1} p_{T} (GeV)", "gen #tau_{2} p_{T} (GeV)"]

outFile_dir = absolute_path+"final_plot_for_real/"

names = {
    "outFile":" ",
    "histTitle":"Efficiency",
    "xAxisTitle":" gen #tau p_{T} (GeV)",
    "yAxisTitle":"algorithmic efficiency"
}

# algo eff cutbased - algo eff dnn
eventTuple_file = ROOT.TFile((evtTuple_dir+"efficiency_algo_hltDoubleL2IsoTau26eta2p2.root"), "READ")
tauTuple_file = ROOT.TFile((tauTuple_dir+"Efficiencies.root"), "READ")

histNum_cb = ROOT.TH1D(eventTuple_file.Get( "passed_algo" ) )
histDen_cb = ROOT.TH1D(eventTuple_file.Get( "total_algo" ) )

histNum_dnn = ROOT.TH1D(tauTuple_file.Get( "passed_d" ) )
histDen_dnn = ROOT.TH1D(tauTuple_file.Get( "total_d") )
names["histTitle"] = "Algorithmic Efficiency"
names["outFile"] = outFile_dir+"/Algorithmic_Efficiency_diagonal"
beautify_plot_canvas(histNum_cb, histDen_cb, histNum_dnn, histDen_dnn, names)

# eff cutbased - eff DNN
eventTuple_file = ROOT.TFile((evtTuple_dir+"efficiency_normal_hltDoubleL2IsoTau26eta2p2.root"), "READ")
tauTuple_file = ROOT.TFile((tauTuple_dir+"Efficiencies.root"), "READ")

histNum_cb = ROOT.TH1D(eventTuple_file.Get( "passed_normal" ) )
histDen_cb = ROOT.TH1D(eventTuple_file.Get( "total_normal" ) )

histNum_dnn = ROOT.TH1D(tauTuple_file.Get( "passed_d" ) )
histDen_dnn = ROOT.TH1D(eventTuple_file.Get( "total_normal") )
names["histTitle"] = "Absolute L2 Efficiency"
names["outFile"] = outFile_dir+"/L2_Efficiency_diagonal"
beautify_plot_canvas(histNum_cb, histDen_cb, histNum_dnn, histDen_dnn, names)


# BOR eff cutbased - BOR  eff  dnn (it should be the same) - 2 different plots
eventTuple_file = ROOT.TFile((evtTuple_dir+"efficiency_normal_hltL1sDoubleTauBigOR.root"), "READ")
tauTuple_file = ROOT.TFile((tauTuple_dir+"Efficiencies.root"), "READ")

histNum_cb = ROOT.TH1D(eventTuple_file.Get( "passed_normal" ) )
histDen_cb = ROOT.TH1D(eventTuple_file.Get( "total_normal" ) )
names["histTitle"] = "Absolute L1 Efficiency"
names["outFile"] = outFile_dir+"/L1eff_abs_diag_cb"
beautify_plot_canvas_1hist(histNum_cb, histDen_cb, names)

histNum_dnn = ROOT.TH1D(tauTuple_file.Get( "passed_d" ) )
histDen_dnn = ROOT.TH1D(eventTuple_file.Get( "total_normal") )
names["histTitle"] = "Absolute L1 Efficiency"
names["outFile"] = outFile_dir+"/L1_Efficiency_diagonal"
beautify_plot_canvas_1hist(histNum_dnn, histDen_dnn, names)
