# algo efficiency for test dataset to compare with sqrt(diagonal) of event cutbased efficiency
# how to build this algorithm:

# 1. goal = compare algo efficiency per taus (NOT per events) on testing dataset with sqrt(diagonal) of cut based (algo?) efficiency (on VBF) per event (NOT per taus) in pt bins

from  TauMLTools.Training.python.produceGridDatasets import *
import ROOT
from array import array
import argparse

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

def beautify_plot_canvas(names, save=True):
    H_ref = 1000;
    W_ref = 1000;
    W = W_ref
    H  = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref
    L = 0.12*W_ref
    R = 0.1*W_ref
    file = ROOT.TFile( names["inFile"], "READ" )
    if(len(names["histogram"]))<2 :
        effHisto2D = ROOT.TH2D(file.Get(names["histogram"][0]))
    else:
        effHisto2D = ROOT.TH2D(file.Get(names["histogram"][0]))
        effHisto2D.Divide( ROOT.TH2D(file.Get(names["histogram"][1])))
        effHisto2D.SetStats(0)
    canvas=ROOT.TCanvas("", "",H_ref,W_ref)
    canvas.cd()
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    canvas.SetTopMargin( T/H )
    canvas.SetBottomMargin( B/H )
    canvas.SetTickx(0)
    canvas.SetTicky(0)
    if(names["logX"]):
        canvas.SetLogx()
        effHisto2D.GetXaxis().SetMoreLogLabels()
        effHisto2D.GetXaxis().SetTitleOffset(1.7)
    if(names["logY"]):
        canvas.SetLogy()
        effHisto2D.GetYaxis().SetMoreLogLabels()
        effHisto2D.GetYaxis().SetTitleOffset(1.7)

    effHisto2D.SetTitle(names["histogramTitle"])
    effHisto2D.GetXaxis().SetTitle(names["xAxis"]) #("gen #tau_{1} p_{T} (GeV)")
    effHisto2D.GetYaxis().SetTitle(names["yAxis"]) #("gen #tau_{2} p_{T} (GeV)")
    print(names["minYRange"])
    effHisto2D.GetZaxis().SetRangeUser(names["minYRange"], names["maxYRange"]) #("gen #tau_{2} p_{T} (GeV)")
    effHisto2D.Draw("TEXT2 COLZ")
    canvas.Update()
    canvas.Print(("{}.png").format(names["outFile"]),"png")
    canvas.Print(("{}.pdf").format(names["outFile"]),"pdf")
    #input()

def draw_AbsoluteCNNEff(names, save=True):
    H_ref = 1000;
    W_ref = 1000;
    W = W_ref
    H  = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref
    L = 0.12*W_ref
    R = 0.1*W_ref
    bin_dict = {"35": {"35": [124.0, 487.0], "70": [0.0, 1.0], "40": [0.0, 1.0], "80": [0.0, 1.0], "50": [0.0, 1.0], "20": [54.0, 1794.0], "90": [0.0, 1.0], "60": [0.0, 1.0], "30": [138.0, 965.0]}, "70": {"35": [299.0, 744.0], "70": [74.0, 96.0], "40": [585.0, 1032.0], "80": [0.0, 1.0], "50": [363.0, 549.0], "20": [101.0, 1888.0], "90": [0.0, 1.0], "60": [294.0, 397.0], "30": [220.0, 874.0]}, "40": {"35": [535.0, 1734.0], "70": [0.0, 1.0], "40": [639.0, 1531.0], "80": [0.0, 1.0], "50": [0.0, 1.0], "20": [151.0, 3388.0], "90": [0.0, 1.0], "60": [0.0, 1.0], "30": [330.0, 1733.0]}, "80": {"35": [201.0, 504.0], "70": [149.0, 187.0], "40": [358.0, 606.0], "80": [47.0, 54.0], "50": [263.0, 371.0], "20": [69.0, 1419.0], "90": [0.0, 1.0], "60": [190.0, 247.0], "30": [164.0, 637.0]}, "50": {"35": [488.0, 1367.0], "70": [0.0, 1.0], "40": [1153.0, 2345.0], "80": [0.0, 1.0], "50": [498.0, 822.0], "20": [176.0, 2982.0], "90": [0.0, 1.0], "60": [0.0, 1.0], "30": [347.0, 1471.0]}, "20": {"35": [0.0, 1.0], "70": [0.0, 1.0], "40": [0.0, 1.0], "80": [0.0, 1.0], "50": [0.0, 1.0], "20": [13.0, 1801.0], "90": [0.0, 1.0], "60": [0.0, 1.0], "30": [0.0, 1.0]}, "90": {"35": [448.0, 1028.0], "70": [375.0, 435.0], "40": [868.0, 1363.0], "80": [287.0, 313.0], "50": [662.0, 883.0], "20": [231.0, 3623.0], "90": [544.0, 586.0], "60": [494.0, 606.0], "30": [341.0, 1296.0]}, "60": {"35": [453.0, 1126.0], "70": [0.0, 1.0], "40": [921.0, 1685.0], "80": [0.0, 1.0], "50": [640.0, 971.0], "20": [131.0, 2430.0], "90": [0.0, 1.0], "60": [185.0, 267.0], "30": [267.0, 1149.0]}, "30": {"35": [0.0, 1.0], "70": [0.0, 1.0], "40": [0.0, 1.0], "80": [0.0, 1.0], "50": [0.0, 1.0], "20": [35.0, 1841.0], "90": [0.0, 1.0], "60": [0.0, 1.0], "30": [32.0, 451.0]}}
    pt_bins = [20,30,35,40,50,60,70,80,90,350]
    #{"pt_1_bin":{pt_2_bin:[num, den]}}
    den2DHist = ROOT.TH2D("den","den",len(pt_bins)-1,array('d',pt_bins),len(pt_bins)-1,array('d',pt_bins))
    for j in range(0,len(pt_bins)-1):
        for k in range(0,len(pt_bins)-1):
            if(j>=k):
                den = bin_dict[str(pt_bins[j])][str(pt_bins[k])][1]
                den2DHist.SetBinContent(j+1,k+1,den)
    c = ROOT.TCanvas()
    c.cd()
    den2DHist.Draw("TEXT2 COLZ")
    c.Update()
    #input()
    file = ROOT.TFile( names["inFile"], "READ" )
    effHisto2D = ROOT.TH2D(file.Get(names["histogram"][0]))
    effHisto2D.Divide(den2DHist)
    effHisto2D.SetStats(0)
    canvas=ROOT.TCanvas("", "",H_ref,W_ref)
    canvas.cd()
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    canvas.SetTopMargin( T/H )
    canvas.SetBottomMargin( B/H )
    canvas.SetTickx(0)
    canvas.SetTicky(0)
    if(names["logX"]):
        canvas.SetLogx()
        effHisto2D.GetXaxis().SetMoreLogLabels()
        effHisto2D.GetXaxis().SetTitleOffset(1.7)
    if(names["logY"]):
        canvas.SetLogy()
        effHisto2D.GetYaxis().SetMoreLogLabels()
        effHisto2D.GetYaxis().SetTitleOffset(1.7)

    effHisto2D.SetTitle(names["histogramTitle"])
    effHisto2D.GetXaxis().SetTitle(names["xAxis"]) #("gen #tau_{1} p_{T} (GeV)")
    effHisto2D.GetYaxis().SetTitle(names["yAxis"]) #("gen #tau_{2} p_{T} (GeV)")
    print(names["minYRange"])
    effHisto2D.GetZaxis().SetRangeUser(names["minYRange"], names["maxYRange"]) #("gen #tau_{2} p_{T} (GeV)")
    effHisto2D.Draw("TEXT2 COLZ")
    canvas.Update()
    #canvas.Print(("{}.png").format(names["outFile"]),"png")
    canvas.Print(("{}.pdf").format(names["outFile"]),"pdf")
    #input()
names = {
    "inFile":"",
    "outFile":"",
    "histogram":[],
    "histogramTitle":"",
    "xAxis":"",
    "yAxis":"",
    "logX":True,
    "logY":True,
    "minYRange":0,
    "maxYRange":0.99,
}

# *** Define all files *****

Eff_BOR_cb =("{}/efficiency_normal_hltL1sDoubleTauBigOR_{}.root").format(evtfileDir, args.sample)
eventTuple_BOR_file = ROOT.TFile(Eff_BOR_cb, "READ")

Eff_algo_cb = ("{}/efficiency_algo_hltDoubleL2IsoTau26eta2p2_{}.root").format(evtfileDir, args.sample)
eventTuple_file_algo = ROOT.TFile(Eff_algo_cb, "READ")

Eff_abs_cb =("{}/efficiency_normal_hltDoubleL2IsoTau26eta2p2_{}.root").format(evtfileDir, args.sample)
eventTuple_file_abs = ROOT.TFile(Eff_abs_cb, "READ")


rate = 3
eff_CNN_3kHz = ("{}/EfficienciesCNN_forRate{}kHz.root").format(GetRateDir(args.machine,rate), rate)
tauTuple_file3 = ROOT.TFile(eff_CNN_3kHz, "READ")

rate = 4
eff_CNN_4kHz = ("{}/EfficienciesCNN_forRate{}kHz.root").format(GetRateDir(args.machine,rate), rate)
tauTuple_file4 = ROOT.TFile(eff_CNN_4kHz, "READ")

rate = 5
eff_CNN_5kHz = ("{}/EfficienciesCNN_forRate{}kHz.root").format(GetRateDir(args.machine,rate), rate)
tauTuple_file5 = ROOT.TFile(eff_CNN_5kHz, "READ")

plotDir = GetPlotDir(args.machine)+"/"

#histNum_algo_cb = ROOT.TH1D(eventTuple_file_algo.Get( "passed_algo" ) )
#histDen_algo_cb = ROOT.TH1D(eventTuple_file_algo.Get( "total_algo" ) )

all_files=[Eff_algo_cb, Eff_abs_cb, Eff_BOR_cb, eff_CNN_3kHz,eff_CNN_4kHz,eff_CNN_5kHz]
x_y_titles = ["gen #tau_{1} p_{T} (GeV)", "gen #tau_{2} p_{T} (GeV)"]


#histogram_names = [["efficiency_algo"], ["efficiency_normal"], ["efficiency_normal"], ]
#histogram_titles = ["Algorithmic L2 Efficiency", "Absolute L2 Efficiency ", "Absolute L1 Efficiency"]
#output_files = ["L2_algorithmic_efficiency","L2_absolute_efficiency", "L1_absolute_efficiency" ]
histogram_names = [["efficiency_algo"], ["efficiency_normal"], ["efficiency_normal"], ["passed", "total"],["passed", "total"],["passed", "total"]]
histogram_titles = ["Algorithmic L2 Efficiency", "Absolute L2 Efficiency ", "Absolute L1 Efficiency", "Algorithmic CNN Efficiency","Algorithmic CNN Efficiency","Algorithmic CNN Efficiency"]
output_files = ["L2_algorithmic_efficiency","L2_absolute_efficiency", "L1_absolute_efficiency", "CNN3_algorithmic_efficiency","CNN3_absolute_efficiency",  "CNN4_algorithmic_efficiency","CNN4_absolute_efficiency",  "CNN4_algorithmic_efficiency","CNN4_absolute_efficiency" ]

for i in range(0, len(all_files)):
    names["inFile"]= all_files[i]
    names["outFile"]= plotDir+output_files[i]
    names["histogram"]=histogram_names[i]

    names["histogramTitle"]=histogram_titles[i]
    names["xAxis"]=x_y_titles[0]
    names["yAxis"]=x_y_titles[1]
    names["logX"]= True
    names["logY"]=True
    if(output_files[i]=="L2_algorithmic_efficiency"):
        names["minYRange"]=0.5
        names["maxYRange"]=0.91
    else:
        names["minYRange"]=0
    beautify_plot_canvas(names)




all_CNN_files=[ eff_CNN_3kHz,eff_CNN_4kHz,eff_CNN_5kHz]

histogram_names = [["passed", "total"],["passed", "total"],["passed", "total"]]
rates = [3,4,5]
output_files = ["CNN3_algorithmic_efficiency","CNN4_algorithmic_efficiency", "CNN5_algorithmic_efficiency" ]
histogram_titles = ["Algorithmic CNN Efficiency","Algorithmic CNN Efficiency","Algorithmic CNN Efficiency"]


for i in range(0, len(all_files)):
    names["inFile"]= all_CNN_files[i]
    names["outFile"] = ('{}/Rate_{}kHz/{}').format(plotDir, rates[i], output_files[i])

    #names["outFile"]= plotDir+output_files[i]
    names["histogram"]=histogram_names[i]
    names["histogramTitle"]=histogram_titles[i]
    names["xAxis"]=x_y_titles[0]
    names["yAxis"]=x_y_titles[1]
    names["logX"]= True
    names["logY"]=True
    names["minYRange"]=0.8
    beautify_plot_canvas(names)

x_y_titles = ["gen #tau_{1} p_{fT} (GeV)", "gen #tau_{2} p_{T} (GeV)"]
names["minYRange"]=0.
names["histogram"]= ["passed"]
names["histogramTitle"]= "Absolute CNN Efficiency"
names["xAxis"]=x_y_titles[0]
names["yAxis"]=x_y_titles[1]
names["logX"]= True
names["logY"]=True
names["inFile"]= eff_CNN_3kHz
names["outFile"] = ('{}/Rate_3kHz/CNN_absolute_efficiency_3kHz').format(plotDir)
draw_AbsoluteCNNEff(names)
names["inFile"]= eff_CNN_4kHz
names["outFile"] = ('{}/Rate_4kHz/CNN_absolute_efficiency_4kHz').format(plotDir)
draw_AbsoluteCNNEff(names)
names["inFile"]= eff_CNN_5kHz
names["outFile"] = ('{}/Rate_5kHz/CNN_absolute_efficiency_5kHz').format(plotDir)
draw_AbsoluteCNNEff(names)
