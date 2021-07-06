# algo efficiency for test dataset to compare with sqrt(diagonal) of event cutbased efficiency
# how to build this algorithm:

# 1. goal = compare algo efficiency per taus (NOT per events) on testing dataset with sqrt(diagonal) of cut based (algo?) efficiency (on VBF) per event (NOT per taus) in pt bins

from  TauMLTools.Training.python.produceGridDatasets import *
import argparse
import scipy
import statsmodels.stats.proportion as ssp
import ROOT

ROOT.gStyle.SetPaintTextFormat(".2f")
parser = argparse.ArgumentParser()
parser.add_argument('--machine', required=False, type=str, default="local", choices=["local", "cmssimphase2"])
parser.add_argument('--corr', required=False, type=str, default="noCorr", choices = ["noCorr","Corr"])
parser.add_argument('--sample', required=False, type=str, default="VBF", choices = ["Zprime","VBF"])
parser.add_argument('--verbose', required=False, type=int, default=0)
args = parser.parse_args()
params = {
    'num_dense_layers':3,
    'num_CNN1x1_layers':4,
    'num_CNN2x2_layers':4,
    'activation_dense':'relu',
    'activation_CNN1x1':'relu',
    'activation_CNN2x2':'relu',
    'nFilters_dense_0':40,
    'nFilters_dense_1':40,
    'nFilters_dense_2':20,
    'nFilters_CNN1x1_0':80,
    'nFilters_CNN1x1_1':60,
    'nFilters_CNN1x1_2':40,
    'nFilters_CNN1x1_3':20,
    'nFilters_CNN2x2_0':20,
    'nFilters_CNN2x2_1':20,
    'nFilters_CNN2x2_2':20,
    'nFilters_CNN2x2_3':40,
    'dropout_dense_layers':0,#0.2,
    'dropout_CNN1x1_layers':0,#0.2,
    'dropout_CNN2x2_layers':0,#0.2,
    #'n_units' : len(CellGridMatrix[0,0,0]),
    'batch_size':200,
    'train_fraction': 0.6102414433536747,
    'validation_fraction': 0.19474661713982488,
    #'learning_rate':0.002,
    'learning_rate':0.001,
    'monitor_es':'val_loss',
    'patience_es':10,
    'epochs':100000,
    'bigOrRate': 13603.37,

    'opt_threshold_3': 0.180858813224404,
    'opt_threshold_4': 0.12267940863785043,
    'opt_threshold_5': 0.08411243185219064,
}
evtTuplePath = "/Users/valeriadamante/Desktop/Dottorato/cmssimphase2/"
evtfileDir = ("{}/{}_{}").format(evtTuplePath, args.sample, args.corr)
CNNfileDir = ("{}/{}_CNN/").format(evtTuplePath, args.sample)
outDir= GetOutPath(args.machine)
plotDir = GetPlotDir(args.machine)

def beautify_plot_canvas(histNum_cb, histDen_cb, histNum_cnn3, histDen_cnn3, histNum_cnn4,histDen_cnn4,histNum_cnn5, histDen_cnn5, names):
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
    gr_cnn3= ROOT.TGraphAsymmErrors(n_bins)
    gr_cnn4= ROOT.TGraphAsymmErrors(n_bins)
    gr_cnn5= ROOT.TGraphAsymmErrors(n_bins)

    for i in range(1,n_bins+1):
        num_cb = histNum_cb.GetBinContent(i)
        den_cb = histDen_cb.GetBinContent(i)
        num_cnn3 = histNum_cnn3.GetBinContent(i)
        den_cnn3 = histDen_cnn3.GetBinContent(i)
        num_cnn4 = histNum_cnn4.GetBinContent(i)
        den_cnn4 = histDen_cnn4.GetBinContent(i)
        num_cnn5 = histNum_cnn5.GetBinContent(i)
        den_cnn5 = histDen_cnn5.GetBinContent(i)

        eff_cb = math.sqrt(num_cb/den_cb)
        eff_cnn3 = math.sqrt(num_cnn3/den_cnn3)
        eff_cnn4 = math.sqrt(num_cnn4/den_cnn4)
        eff_cnn5 = math.sqrt(num_cnn5/den_cnn5)
        if(names["verbose"]>0):
            print(("num cb = {}, den cb = {}, eff = {}").format(num_cb, den_cb, eff_cb))
            print(("num cnn3 = {}, den cnn3 = {}, eff = {}").format(num_cnn3, den_cnn3, eff_cnn3))
            print(("num cnn4 = {}, den cnn4 = {}, eff = {}").format(num_cnn4, den_cnn4, eff_cnn4))
            print(("num cnn5 = {}, den cnn5 = {}, eff = {}").format(num_cnn5, den_cnn5, eff_cnn5))
            print("\n")
        x_p = pt_bins[i] - ((pt_bins[i]-pt_bins[k])/2)
        gr_cb.SetPointX(k,x_p)
        gr_cnn3.SetPointX(k,x_p)
        gr_cnn4.SetPointX(k,x_p)
        gr_cnn5.SetPointX(k,x_p)
        gr_cb.SetPointY(k,eff_cb)
        gr_cnn3.SetPointY(k,eff_cnn3)
        gr_cnn4.SetPointY(k,eff_cnn4)
        gr_cnn5.SetPointY(k,eff_cnn5)
        err_x_low = (pt_bins[i]-pt_bins[k])/2
        err_x_up= (pt_bins[i]-pt_bins[k])/2
        c_low_cb, c_up_cb = ssp.proportion_confint(num_cb, den_cb, alpha=1-0.68, method='beta')
        c_low_cnn3, c_up_cnn3 = ssp.proportion_confint(num_cnn3, den_cnn3, alpha=1-0.68, method='beta')
        c_low_cnn4, c_up_cnn4 = ssp.proportion_confint(num_cnn4, den_cnn4, alpha=1-0.68, method='beta')
        c_low_cnn5, c_up_cnn5 = ssp.proportion_confint(num_cnn5, den_cnn5, alpha=1-0.68, method='beta')
        err_y_low_cb = eff_cb-math.sqrt(c_low_cb)
        err_y_up_cb = math.sqrt(c_up_cb) - eff_cb
        err_y_low_cnn3 = eff_cnn3-math.sqrt(c_low_cnn3)
        err_y_up_cnn3 = math.sqrt(c_up_cnn3) - eff_cnn3
        err_y_low_cnn4 = eff_cnn4-math.sqrt(c_low_cnn4)
        err_y_up_cnn4 = math.sqrt(c_up_cnn4) - eff_cnn4
        err_y_low_cnn5 = eff_cnn5-math.sqrt(c_low_cnn5)
        err_y_up_cnn5 = math.sqrt(c_up_cnn5) - eff_cnn5
        gr_cb.SetPointError(k, err_x_low, err_x_up, err_y_low_cb, err_y_up_cb)
        gr_cnn3.SetPointError(k, err_x_low, err_x_up, err_y_low_cnn3, err_y_up_cnn3)
        gr_cnn4.SetPointError(k, err_x_low, err_x_up, err_y_low_cnn4, err_y_up_cnn4)
        gr_cnn5.SetPointError(k, err_x_low, err_x_up, err_y_low_cnn5, err_y_up_cnn5)
        k+=1
    canvas=ROOT.TCanvas()
    canvas.cd()
    legend = ROOT.TLegend(0.75,0.6,0.88,0.7)
    gr_cb.SetMarkerColor( ROOT.kRed )
    gr_cb.SetLineColor( ROOT.kRed )
    gr_cb.SetMarkerStyle( 20 )
    gr_cb.SetTitle(("{};{};{}").format(names["histTitle"],names["xAxisTitle"],names["yAxisTitle"]))
    legend.AddEntry(gr_cb, "cut based", "lep")
    gr_cnn3.SetMarkerColor( ROOT.kGreen )
    gr_cnn3.SetLineColor( ROOT.kGreen )
    gr_cnn3.SetMarkerStyle( 20 )
    gr_cnn3.SetTitle(("{};{};{}").format(names["histTitle"],names["xAxisTitle"],names["yAxisTitle"]))
    legend.AddEntry(gr_cnn3, "CNN rate 3kHz", "lep")
    gr_cnn4.SetMarkerColor( ROOT.kBlue )
    gr_cnn4.SetLineColor( ROOT.kBlue )
    gr_cnn4.SetMarkerStyle( 20 )
    gr_cnn4.SetTitle(("{};{};{}").format(names["histTitle"],names["xAxisTitle"],names["yAxisTitle"]))
    legend.AddEntry(gr_cnn4, "CNN rate 4kHz", "lep")
    gr_cnn5.SetMarkerColor( ROOT.kOrange )
    gr_cnn5.SetLineColor( ROOT.kOrange )
    gr_cnn5.SetMarkerStyle( 20 )
    gr_cnn5.SetTitle(("{};{};{}").format(names["histTitle"],names["xAxisTitle"],names["yAxisTitle"]))
    legend.AddEntry(gr_cnn5, "CNN rate 5kHz", "lep")
    x_min = names["x_min"]
    x_max = names["x_max"]
    y_min = names["y_min"]
    y_max = names["y_max"]
    gr_cnn3.GetXaxis().SetRangeUser(x_min,x_max)
    gr_cnn3.GetYaxis().SetRangeUser(y_min,y_max)
    gr_cnn3.Draw("AP")
    gr_cnn4.Draw("PSAME")
    gr_cnn5.Draw("PSAME")
    gr_cb.Draw("PSAME")
    legend.Draw("SAME")
    canvas.Update()
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
        if(names["verbose"]>0):
            print(("num = {}, den = {}, eff = {}").format(num, den, eff))
            print("\n")
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
    x_min = names["x_min"]
    x_max = names["x_max"]
    y_min = names["y_min"]
    y_max = names["y_max"]
    gr.GetXaxis().SetRangeUser(x_min,x_max)
    gr.GetYaxis().SetRangeUser(y_min,y_max)
    gr.SetTitle(("{};{};{}").format(names["histTitle"],names["xAxisTitle"],names["yAxisTitle"]))
    gr.Draw("AP")
    canvas.Update()
    if(save):
        canvas.Update()
        canvas.Print(("{}.pdf").format(names["outFile"]),"pdf")

names = {
    "outFile":" ",
    "histTitle":"Efficiency",
    "xAxisTitle":" gen #tau p_{T} (GeV)",
    "yAxisTitle":"algorithmic efficiency",
    "x_min":20,
    "x_max":400,
    "y_min":0.,
    "y_max":1.1,
    "verbose":args.verbose,
}

# *** Draw superimposed algo efficiencies *****

Eff_algo_cb = ("{}/efficiency_algo_hltDoubleL2IsoTau26eta2p2_{}.root").format(evtfileDir, args.sample)
eventTuple_file_algo = ROOT.TFile(Eff_algo_cb, "READ")
histNum_algo_cb = ROOT.TH1D(eventTuple_file_algo.Get( "passed_algo" ) )
histDen_algo_cb = ROOT.TH1D(eventTuple_file_algo.Get( "total_algo" ) )

rate = 3
eff_CNN_3kHz = ("{}/EfficienciesCNN_forRate{}kHz.root").format(GetRateDir(args.machine,rate), rate)
tauTuple_file3 = ROOT.TFile(eff_CNN_3kHz, "READ")
histNum_cnn3 = ROOT.TH1D(tauTuple_file3.Get( "passed_d" ) )
histDen_cnn3 = ROOT.TH1D(tauTuple_file3.Get( "total_d") )

rate = 4
eff_CNN_4kHz = ("{}/EfficienciesCNN_forRate{}kHz.root").format(GetRateDir(args.machine,rate), rate)
tauTuple_file4 = ROOT.TFile(eff_CNN_4kHz, "READ")
histNum_cnn4 = ROOT.TH1D(tauTuple_file4.Get( "passed_d" ) )
histDen_cnn4 = ROOT.TH1D(tauTuple_file4.Get( "total_d") )

rate = 5
eff_CNN_5kHz = ("{}/EfficienciesCNN_forRate{}kHz.root").format(GetRateDir(args.machine,rate), rate)
tauTuple_file5 = ROOT.TFile(eff_CNN_5kHz, "READ")
histNum_cnn5 = ROOT.TH1D(tauTuple_file5.Get( "passed_d" ) )
histDen_cnn5 = ROOT.TH1D(tauTuple_file5.Get( "total_d") )

names["histTitle"] = "Algorithmic Efficiency"
names["yAxisTitle"] = "algo efficiency"
names["y_min"] = 0.55
names["y_max"] = 1.009
names["outFile"] =("{}/Algorithmic_Efficiency_diagonal_{}_{}").format(plotDir, args.sample, args.corr)
beautify_plot_canvas(histNum_algo_cb, histDen_algo_cb, histNum_cnn3, histDen_cnn3, histNum_cnn4,histDen_cnn4,histNum_cnn5, histDen_cnn5, names)

#  eff cutbased - eff DNN

Eff_abs_cb =("{}/efficiency_normal_hltDoubleL2IsoTau26eta2p2_{}.root").format(evtfileDir, args.sample)
eventTuple_file_abs = ROOT.TFile(Eff_abs_cb, "READ")
histNum_abs_cb = ROOT.TH1D(eventTuple_file_abs.Get( "passed_normal" ) )
histDen_abs_cb = ROOT.TH1D(eventTuple_file_abs.Get( "total_normal" ) )

names["histTitle"] = "Absolute Efficiency"
names["yAxisTitle"] = "efficiency"
names["outFile"] =("{}/Absolute_Efficiency_diagonal_{}_{}").format(plotDir, args.sample, args.corr)
beautify_plot_canvas(histNum_abs_cb, histDen_abs_cb, histNum_cnn3, histDen_abs_cb, histNum_cnn4,histDen_abs_cb,histNum_cnn5, histDen_abs_cb, names)

# 1d single plots - quite unuseful !!
names["histTitle"] = "Absolute Efficiency for CB"
names["outFile"] = plotDir+"/L2eff_abs_diag_cb"
names["outFile"] =("{}/AbsEff_CNN3kHz_{}_{}").format(plotDir, args.sample, args.corr)
beautify_plot_canvas_1hist(histNum_abs_cb, histDen_abs_cb, names)

names["histTitle"] = "Absolute L1 Efficiency"
names["outFile"] = ('{}/Rate_3kHz/L1_Efficiency_diagonal_3kHz').format(plotDir)
beautify_plot_canvas_1hist(histNum_cnn3, histDen_abs_cb, names)

names["histTitle"] = "Absolute L1 Efficiency"
names["outFile"] = ('{}/Rate_4kHz/L1_Efficiency_diagonal_4kHz').format(plotDir)
beautify_plot_canvas_1hist(histNum_cnn4, histDen_abs_cb, names)

names["histTitle"] = "Absolute L1 Efficiency"
names["outFile"] = ('{}/Rate_5kHz/L1_Efficiency_diagonal_5kHz').format(plotDir)
beautify_plot_canvas_1hist(histNum_cnn5, histDen_abs_cb, names)


# BOR eff cutbased - BOR  eff  dnn (it should be the same) - 2 different plots

Eff_BOR_cb =("{}/efficiency_normal_hltL1sDoubleTauBigOR_{}.root").format(evtfileDir, args.sample)
eventTuple_BOR_file = ROOT.TFile(Eff_BOR_cb, "READ")
histNum_BOR_cb = ROOT.TH1D(eventTuple_BOR_file.Get( "passed_normal" ) )
histDen_BOR_forAll = ROOT.TH1D(eventTuple_BOR_file.Get( "total_normal" ) )
names["y_min"]=0.
names["y_max"]=1.01
names["histTitle"] = "Absolute L1 Efficiency"
names["yAxisTitle"] = "efficiency"
names["outFile"] =("{}/Absolute_L1Efficiency_diagonal_{}_{}").format(plotDir, args.sample, args.corr)
#beautify_plot_canvas_1hist(histNum_BOR_cb, histDen_BOR_forAll, names)


#names["histTitle"] = "Absolute L1 Efficiency for CB"
#names["outFile"] = plotDir+"/L1eff_abs_diag_cb"
#beautify_plot_canvas_1hist(histNum_BOR_cb, histDen_BOR_forAll, names)

#names["histTitle"] = "Absolute L1 Efficiency"
#names["outFile"] = ('{}/Rate_3kHz/L1_Efficiency_diagonal_3kHz').format(plotDir)
#beautify_plot_canvas_1hist(histNum_cnn3, histDen_BOR_forAll, names)

#names["histTitle"] = "Absolute L1 Efficiency"
#names["outFile"] = ('{}/Rate_4kHz/L1_Efficiency_diagonal_4kHz').format(plotDir)
#beautify_plot_canvas_1hist(histNum_cnn4, histDen_BOR_forAll, names)

#names["histTitle"] = "Absolute L1 Efficiency"
#names["outFile"] = ('{}/Rate_5kHz/L1_Efficiency_diagonal_5kHz').format(plotDir)
#beautify_plot_canvas_1hist(histNum_cnn5, histDen_BOR_forAll, names)
