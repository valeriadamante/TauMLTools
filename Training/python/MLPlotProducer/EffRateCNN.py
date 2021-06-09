from  TauMLTools.Training.python.produceGridDatasets import *
import tensorflow as tf
import root_pandas
import scipy
import statsmodels.stats.proportion as ssp
import argparse
import ROOT
from array import array
ROOT.gStyle.SetPaintTextFormat(".2f")

parser = argparse.ArgumentParser()
parser.add_argument('--machine', required=False, type=str, default="local", choices=["local", "cmssimphase2"]) #aggiungi pccms65
parser.add_argument('--effRate', required= False, type=str, default = 'test', choices=['test', 'eff','rate', 'Zprime'])
parser.add_argument('--rateValue', required= False, type=int, default = 4, choices=[3,4,5])
parser.add_argument('--n_max_events', required=False, type=int, default=-1, help='max number of events to be processed')
parser.add_argument('--McTruth_file', required=False, type=str, default='MCTruth.npy', help='output file name')
parser.add_argument('--Weights_file', required=False, type=str, default='weights.npy', help='output file name')
parser.add_argument('--n_cellsX', required=False, type=int, default=5, help='number of cells along X dir')
parser.add_argument('--n_cellsY', required=False, type=int, default=5, help='number of cells along Y dir')
parser.add_argument('--verbose', required=False, type=int, default=1)
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
    'opt_threshold_3':0.1277466341834952,
    'opt_threshold_4': 0.08153110369858041,
    'opt_threshold_5': 0.051069704813926364,
}
# ***** Get cell grid Matrix *****

kwArgs = {'n_max_events':args.n_max_events, 'n_cellsX':args.n_cellsX, 'n_cellsY':args.n_cellsY, 'timeInfo' : False, 'verbose' : args.verbose}
CellGridMatrix = GetCellGridNormMatrix(args.machine, args.effRate, **kwArgs)
# ***** Get MC and weights Matrices *****
if(args.effRate== 'test'):
    MCTruth, weights = GetMCTruthWeights(args.machine, args.McTruth_file, args.Weights_file, **kwArgs)
    # ****** Get train - test - validation samples ********
    if(args.verbose)>0:
        print("preparing train/test/val samples")
    number_of_batches = len(MCTruth)/params['batch_size']
    x_train, y_train, w_train, x_test, y_test, w_test, x_val, y_val, w_val = GetTrainTestFraction(MCTruth, CellGridMatrix, weights, params['train_fraction'],  params['validation_fraction'], args.verbose)
print(CellGridMatrix.shape)

# ****** Get DataFrames ******
#y_predict_test = model.predict(x_test).reshape(y_test.shape)
variables = ['evt','l1Tau_pt', 'l1Tau_hwIso']
if(args.effRate != 'rate'):
    variables = ['evt','l1Tau_pt', 'l1Tau_hwIso', 'genLepton_vis_pt','genLepton_vis_eta','genLepton_isTau']
inFile = GetRootPath(args.machine, args.effRate)
treeName = 'L2TauTrainTuple'

absolute_path=GetDataSetPath(args.machine)
dataFrameWithPredName = ("{}/dfWithPredictions_{}.root").format(absolute_path,GetNameForEffRate(args.effRate))

# ***** save dataframe with predictions *****
if(not os.path.exists(dataFrameWithPredName)):
    # ***** load model ******
    model = tf.keras.models.load_model(GetModelPath(args.machine, params))
    print("model successfully loaded")
    df = root_pandas.read_root(inFile, treeName, variables)
    df = df[(df.l1Tau_pt>=32) & ((df.l1Tau_hwIso>0) | (df.l1Tau_pt>=70))]
    df['y_predict'] = model.predict(CellGridMatrix).reshape(df['evt'].shape)
    root_pandas.to_root(df, dataFrameWithPredName, treeName)

df = root_pandas.read_root(dataFrameWithPredName,treeName)
#print(df.head())

def get_n_evt(df_orig):
    df = df_orig.groupby("evt").count()
    df = df[df.l1Tau_pt >= 2]
    return df.shape[0]

if(args.effRate=='rate'):
    print(("evaluating for desired rate of {}").format(args.rateValue))
    df_for_num = df[(df.y_predict>params[('opt_threshold_{}').format(args.rateValue)])]
    den = get_n_evt(df)
    num = get_n_evt(df_for_num)
    eff_rate = num*params['bigOrRate']/den
    eff = num/den
    c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
    c_low_rate = c_low*params['bigOrRate']
    c_up_rate = c_up*params['bigOrRate']

    print(("total numerator = \t {} \ntotal denominator = \t {} \nefficiency = \t {} \nvar_up = \t {} \nvar_down = \t {} \nunc_up = \t {} \nunc_down = \t {}").format(num, den, eff, c_up, c_low, c_up-eff, eff-c_low))
    print(("total Rate = \t {} \nvar_up = \t {} \nvar_down = \t {} \nunc_up \t {} \nunc_down \t {}").format(eff_rate,  c_up_rate, c_low_rate, c_up_rate-eff_rate, eff_rate-c_low_rate))


pt_bins = [20,30,35,40,50,60,70,80,90,350]

def save_dataframe(df, path_name, tuple_name, rateValue, num=False):
    if(num):
        df = df[(df.y_predict>params[('opt_threshold_{}').format(rateValue)])]
    df = df.sort_values(by=['evt'])
    df_count = df.groupby('evt').size().reset_index(name='count')
    values = []
    prev_index = 0
    for i in range(0,df.shape[0]):
        for j in range(prev_index,df_count.shape[0]):
            temp_i = int(df.iloc[i]['evt'])
            temp_j = int(df_count.iloc[j]['evt'])
            if(temp_i == temp_j):
                values.append(df_count.iloc[j]['count'])
                prev_index = j
                break
    df['hasTwoHadTaus']=values
    root_pandas.to_root(df, path_name, tuple_name)

if(args.effRate=='eff'):
    dataframeDenName = ("{}/df_for_den_VBF.root").format(absolute_path)
    dataframeNumName = ('{}/df_for_num_VBF_forRate{}kHz.root').format(absolute_path, args.rateValue)
    if(not os.path.exists(dataframeDenName)):
        save_dataframe(df, dataframeDenName, treeName, args.rateValue, False)
    df_for_den = root_pandas.read_root(dataframeDenName, treeName)

    if(not os.path.exists(dataframeNumName) ):
        save_dataframe(df, dataframeNumName, treeName, args.rateValue, True)
    df_for_num = root_pandas.read_root(dataframeNumName, treeName)

    df_for_num_2 = df_for_num[df_for_num.hasTwoHadTaus==2]
    df_for_den_2 = df_for_den[df_for_den.hasTwoHadTaus==2]
    df_for_den_count = df_for_den_2.groupby('evt').count()
    df_for_num_count = df_for_num_2.groupby('evt').count()

    den = df_for_den_count[(df_for_den_count.genLepton_vis_pt==2)].shape[0]
    num = df_for_num_count[(df_for_num_count.genLepton_vis_pt==2)].shape[0]
    print(num)
    print(den)

    eff = num/den
    c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
    print(("total numerator = \t {} \ntotal denominator = \t {} \nefficiency = \t {} \nvar_up = \t {} \nvar_down = \t {} \nunc_up = \t {} \nunc_down = \t {}").format(num, den, eff, c_up, c_low, c_up-eff, eff-c_low))
    #print(("total Rate = \t {} \nvar_up = \t {} \nvar_down = \t {} \nunc_up \t {} \nunc_down \t {}").format(eff_rate,  c_up_rate, c_low_rate, c_up_rate-eff_rate, eff_rate-c_low_rate))

    print(df_for_num_2.shape[0])

    df_dict = {'num': df_for_num_2, 'den': df_for_den_2}

    hist_num_d = ROOT.TH1D("passed_d" ,"passed_d" ,len(pt_bins)-1,array('d',pt_bins))
    hist_den_d = ROOT.TH1D("total_d" ,"total_d" ,len(pt_bins)-1,array('d',pt_bins))
    hist_num_d_sqrt = ROOT.TH1D("passed_d_sqrt" ,"passed_d_sqrt" ,len(pt_bins)-1,array('d',pt_bins))
    hist_den_d_sqrt = ROOT.TH1D("total_d_sqrt" ,"total_d_sqrt" ,len(pt_bins)-1,array('d',pt_bins))
    hist_num = ROOT.TH2D("passed" , "passed" , len(pt_bins)-1, array('d',pt_bins), len(pt_bins)-1,array('d',pt_bins) )
    hist_den = ROOT.TH2D("total" , "total" , len(pt_bins)-1, array('d',pt_bins), len(pt_bins)-1,array('d',pt_bins) )
    hists = { 'num': hist_num, 'den': hist_den  }
    hists_d = { 'num': hist_num_d, 'den': hist_den_d  }
    hists_d_sqrt = { 'num': hist_num_d_sqrt, 'den': hist_den_d_sqrt  }

    for name in ['num', 'den']:
        for ptx_index in range(0, len(pt_bins)-1):
            for pty_index in range(0, len(pt_bins)-1):
                if(pty_index > ptx_index):
                    continue
                x_min, x_max = pt_bins[ptx_index], pt_bins[ptx_index+1]
                y_min, y_max = pt_bins[pty_index], pt_bins[pty_index+1]
                def inside_bin(values, index):
                    x_idx = -1
                    y_idx = -1
                    for n in range(len(values)):
                        pass_x = values[n]>= x_min and values[n]<x_max
                        pass_y = values[n]>= y_min and values[n]<y_max
                        if(pass_x and (x_idx<0 or x_idx==y_idx)):
                            x_idx= n
                        if(pass_y and (y_idx<0 or y_idx==x_idx)):
                            y_idx= n
                    flag = 1 if x_idx>=0 and y_idx >=0  and x_idx != y_idx else 0
                    return flag
                #print(sum(df_dict[name].groupby('evt')['genLepton_vis_pt'].agg(inside_bin, engine="numba")))
                labels = df_dict[name].groupby('evt')['genLepton_vis_pt'].agg(inside_bin, engine="numba")
                n_evt = np.sum(labels)
                hists[name].SetBinContent(ptx_index + 1, pty_index + 1, n_evt)
                hists[name].SetBinError(ptx_index + 1, pty_index + 1, math.sqrt(n_evt))
                if(ptx_index == pty_index):
                    if(name == 'den'):
                        print(("bin [{},{}], evt = {}").format(ptx_index+1, pty_index+1, n_evt))
                    print(name)
                    print(math.sqrt(n_evt))
                    hists_d[name].SetBinContent(ptx_index + 1, n_evt)
                    hists_d_sqrt[name].SetBinContent(ptx_index + 1, math.sqrt(n_evt))
                    #hists_d[name].SetBinError(ptx_index + 1, math.sqrt(math.sqrt(n_evt)))
                #print(n_evt)
    for name in ['num','den']:
        hists[name].GetXaxis().SetTitle("#tau 1 P_{T_} (GeV)")
        hists[name].GetYaxis().SetTitle("#tau 2 P_{T_} (GeV)")
        hists_d[name].GetXaxis().SetTitle("#tau P_{T_} (GeV)")
    print( ("{}/Rate_{}kHz/EfficienciesCNN_forRate{}kHz.root").format(GetPlotDir(args.machine), args.rateValue, args.rateValue))
    myfile = ROOT.TFile( ("{}/Rate_{}kHz/EfficienciesCNN_forRate{}kHz.root").format(GetPlotDir(args.machine), args.rateValue, args.rateValue), 'RECREATE' )
    hist_num.Write()
    hist_den.Write()
    hist_num_d.Write()
    hist_den_d.Write()
    hist_num_d_sqrt.Write()
    hist_den_d_sqrt.Write()
    canvas1=ROOT.TCanvas()
    canvas1.cd()
    canvas1.SetLogx()
    canvas1.SetLogy()
    EffGraph = ROOT.TEfficiency(hists['num'],hists['den'])
    EffGraph.SetTitle("Algorithmic Efficiency; #tau 1 P_{T} (GeV); #tau 2 P_{T} (GeV)")
    EffGraph.Draw("TEXT2 COLZ")
    canvas1.Update()
    canvas1.Print(('{}/Rate_{}kHz/Efficiency_VBF_forRate{}kHz.png').format(GetPlotDir(args.machine),args.rateValue, args.rateValue), "png")
    canvas1.Print(('{}/Rate_{}kHz/Efficiency_VBF_forRate{}kHz.pdf').format(GetPlotDir(args.machine),args.rateValue, args.rateValue), "pdf")

    canvas2=ROOT.TCanvas()
    canvas2.cd()
    canvas2.SetLogx()
    EffGraph_d = ROOT.TEfficiency(hists_d['num'],hists_d['den'])
    EffGraph_d.SetMarkerColor( 4 )
    EffGraph_d.SetMarkerStyle( 21 )
    EffGraph_d.SetTitle("Efficiency;#tau P_{T_} (GeV);#epsilon")
    EffGraph_d.SetMarkerColor(4)
    EffGraph_d.SetMarkerStyle(20)
    EffGraph_d.Draw()
    canvas2.Update()
    canvas2.Print(('{}/Rate_{}kHz/Efficiency_VBF_Diag_forRate{}kHz.png').format(GetPlotDir(args.machine),args.rateValue, args.rateValue), "png")
    canvas2.Print(('{}/Rate_{}kHz/Efficiency_VBF_Diag_forRate{}kHz.pdf').format(GetPlotDir(args.machine),args.rateValue, args.rateValue), "pdf")

    canvas3=ROOT.TCanvas()
    canvas3.cd()
    canvas3.SetLogx()
    EffGraph_d_sqrt = ROOT.TEfficiency(hists_d_sqrt['num'],hists_d_sqrt['den'])
    EffGraph_d_sqrt.SetMarkerColor( 4 )
    EffGraph_d_sqrt.SetMarkerStyle( 21 )
    EffGraph_d_sqrt.SetTitle("Efficiency;#tau P_{T_} (GeV);#epsilon")
    EffGraph_d_sqrt.SetMarkerColor(4)
    EffGraph_d_sqrt.SetMarkerStyle(20)
    EffGraph_d_sqrt.Draw()
    canvas3.Update()
    canvas3.Print(('{}/Rate_{}kHz/Efficiency_VBF_sqrt_forRate{}kHz.png').format(GetPlotDir(args.machine),args.rateValue, args.rateValue), "png")
    canvas3.Print(('{}/Rate_{}kHz/Efficiency_VBF_sqrt_forRate{}kHz.pdf').format(GetPlotDir(args.machine),args.rateValue, args.rateValue), "pdf")


    #save histograms

    EffGraph.Write()
    EffGraph_d.Write()
    EffGraph_d_sqrt.Write()
    myfile.Close()



if(args.effRate=='Zprime'):
    dataframeDenName = ("{}/df_for_den_Zprime.root").format(absolute_path)
    dataframeNumName = ('{}/df_for_num_Zprime_forRate{}kHz.root').format(absolute_path, args.rateValue)
    if(not os.path.exists(dataframeDenName)):
        save_dataframe(df, dataframeDenName, treeName, args.rateValue, False)
    df_for_den = root_pandas.read_root(dataframeDenName, treeName)

    if(not os.path.exists(dataframeNumName) ):
        save_dataframe(df, dataframeNumName, treeName, args.rateValue, True)
    df_for_num = root_pandas.read_root(dataframeNumName, treeName)

    df_for_num_2 = df_for_num[df_for_num.hasTwoHadTaus==2]
    df_for_den_2 = df_for_den[df_for_den.hasTwoHadTaus==2]
    df_for_den_count = df_for_den_2.groupby('evt').count()
    df_for_num_count = df_for_num_2.groupby('evt').count()

    den = df_for_den_count[(df_for_den_count.genLepton_vis_pt==2)].shape[0]
    num = df_for_num_count[(df_for_num_count.genLepton_vis_pt==2)].shape[0]
    print(num)
    print(den)

    eff = num/den
    c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
    print(("total numerator = \t {} \ntotal denominator = \t {} \nefficiency = \t {} \nvar_up = \t {} \nvar_down = \t {} \nunc_up = \t {} \nunc_down = \t {}").format(num, den, eff, c_up, c_low, c_up-eff, eff-c_low))
    #print(("total Rate = \t {} \nvar_up = \t {} \nvar_down = \t {} \nunc_up \t {} \nunc_down \t {}").format(eff_rate,  c_up_rate, c_low_rate, c_up_rate-eff_rate, eff_rate-c_low_rate))

    print(df_for_num_2.shape[0])

    df_dict = {'num': df_for_num_2, 'den': df_for_den_2}

    hist_num_d = ROOT.TH1D("passed_d" ,"passed_d" ,len(pt_bins)-1,array('d',pt_bins))
    hist_den_d = ROOT.TH1D("total_d" ,"total_d" ,len(pt_bins)-1,array('d',pt_bins))
    hist_num_d_sqrt = ROOT.TH1D("passed_d_sqrt" ,"passed_d_sqrt" ,len(pt_bins)-1,array('d',pt_bins))
    hist_den_d_sqrt = ROOT.TH1D("total_d_sqrt" ,"total_d_sqrt" ,len(pt_bins)-1,array('d',pt_bins))
    hist_num = ROOT.TH2D("passed" , "passed" , len(pt_bins)-1, array('d',pt_bins), len(pt_bins)-1,array('d',pt_bins) )
    hist_den = ROOT.TH2D("total" , "total" , len(pt_bins)-1, array('d',pt_bins), len(pt_bins)-1,array('d',pt_bins) )
    hists = { 'num': hist_num, 'den': hist_den  }
    hists_d = { 'num': hist_num_d, 'den': hist_den_d  }
    hists_d_sqrt = { 'num': hist_num_d_sqrt, 'den': hist_den_d_sqrt  }

    for name in ['num', 'den']:
        for ptx_index in range(0, len(pt_bins)-1):
            for pty_index in range(0, len(pt_bins)-1):
                if(pty_index > ptx_index):
                    continue
                x_min, x_max = pt_bins[ptx_index], pt_bins[ptx_index+1]
                y_min, y_max = pt_bins[pty_index], pt_bins[pty_index+1]
                def inside_bin(values, index):
                    x_idx = -1
                    y_idx = -1
                    for n in range(len(values)):
                        pass_x = values[n]>= x_min and values[n]<x_max
                        pass_y = values[n]>= y_min and values[n]<y_max
                        if(pass_x and (x_idx<0 or x_idx==y_idx)):
                            x_idx= n
                        if(pass_y and (y_idx<0 or y_idx==x_idx)):
                            y_idx= n
                    flag = 1 if x_idx>=0 and y_idx >=0  and x_idx != y_idx else 0
                    return flag
                #print(sum(df_dict[name].groupby('evt')['genLepton_vis_pt'].agg(inside_bin, engine="numba")))
                labels = df_dict[name].groupby('evt')['genLepton_vis_pt'].agg(inside_bin, engine="numba")
                n_evt = np.sum(labels)
                hists[name].SetBinContent(ptx_index + 1, pty_index + 1, n_evt)
                hists[name].SetBinError(ptx_index + 1, pty_index + 1, math.sqrt(n_evt))
                if(ptx_index == pty_index):
                    if(name == 'num'):
                        print(("bin [{},{}], evt = {}").format(ptx_index+1, pty_index+1, n_evt))
                    print(name)
                    print(math.sqrt(n_evt))
                    hists_d[name].SetBinContent(ptx_index + 1, n_evt)
                    hists_d_sqrt[name].SetBinContent(ptx_index + 1, math.sqrt(n_evt))
                    #hists_d[name].SetBinError(ptx_index + 1, math.sqrt(math.sqrt(n_evt)))
                #print(n_evt)
    for name in ['num','den']:
        hists[name].GetXaxis().SetTitle("#tau 1 P_{T_} (GeV)")
        hists[name].GetYaxis().SetTitle("#tau 2 P_{T_} (GeV)")
        hists_d[name].GetXaxis().SetTitle("#tau P_{T_} (GeV)")
    print( ("{}/Rate_{}kHz/EfficienciesCNN_ZPrime_forRate{}kHz.root").format(GetPlotDir(args.machine), args.rateValue, args.rateValue))
    myfile = ROOT.TFile( ("{}/Rate_{}kHz/EfficienciesCNN_ZPrime_forRate{}kHz.root").format(GetPlotDir(args.machine), args.rateValue, args.rateValue), 'RECREATE' )
    hist_num.Write()
    hist_den.Write()
    hist_num_d.Write()
    hist_den_d.Write()
    hist_num_d_sqrt.Write()
    hist_den_d_sqrt.Write()
    canvas1=ROOT.TCanvas()
    canvas1.cd()
    canvas1.SetLogx()
    canvas1.SetLogy()
    EffGraph = ROOT.TEfficiency(hists['num'],hists['den'])
    EffGraph.SetTitle("Algorithmic Efficiency; #tau 1 P_{T} (GeV); #tau 2 P_{T} (GeV)")
    EffGraph.Draw("TEXT2 COLZ")
    canvas1.Update()
    canvas1.Print(('{}/Rate_{}kHz/Efficiency_Zprime_forRate{}kHz.png').format(GetPlotDir(args.machine),args.rateValue, args.rateValue), "png")
    canvas1.Print(('{}/Rate_{}kHz/Efficiency_Zprime_forRate{}kHz.pdf').format(GetPlotDir(args.machine),args.rateValue, args.rateValue), "pdf")

    canvas2=ROOT.TCanvas()
    canvas2.cd()
    canvas2.SetLogx()
    EffGraph_d = ROOT.TEfficiency(hists_d['num'],hists_d['den'])
    EffGraph_d.SetMarkerColor( 4 )
    EffGraph_d.SetMarkerStyle( 21 )
    EffGraph_d.SetTitle("Efficiency;#tau P_{T_} (GeV);#epsilon")
    EffGraph_d.SetMarkerColor(4)
    EffGraph_d.SetMarkerStyle(20)
    EffGraph_d.Draw()
    canvas2.Update()
    canvas2.Print(('{}/Rate_{}kHz/Efficiency_Zprime_Diag_forRate{}kHz.png').format(GetPlotDir(args.machine),args.rateValue, args.rateValue), "png")
    canvas2.Print(('{}/Rate_{}kHz/Efficiency_Zprime_Diag_forRate{}kHz.pdf').format(GetPlotDir(args.machine),args.rateValue, args.rateValue), "pdf")

    canvas3=ROOT.TCanvas()
    canvas3.cd()
    canvas3.SetLogx()
    EffGraph_d_sqrt = ROOT.TEfficiency(hists_d_sqrt['num'],hists_d_sqrt['den'])
    EffGraph_d_sqrt.SetMarkerColor( 4 )
    EffGraph_d_sqrt.SetMarkerStyle( 21 )
    EffGraph_d_sqrt.SetTitle("Efficiency;#tau P_{T_} (GeV);#epsilon")
    EffGraph_d_sqrt.SetMarkerColor(4)
    EffGraph_d_sqrt.SetMarkerStyle(20)
    EffGraph_d_sqrt.Draw()
    canvas3.Update()
    canvas3.Print(('{}/Rate_{}kHz/Efficiency_Zprime_sqrt_forRate{}kHz.png').format(GetPlotDir(args.machine),args.rateValue, args.rateValue), "png")
    canvas3.Print(('{}/Rate_{}kHz/Efficiency_Zprime_sqrt_forRate{}kHz.pdf').format(GetPlotDir(args.machine),args.rateValue, args.rateValue), "pdf")


    #save histograms

    EffGraph.Write()
    EffGraph_d.Write()
    EffGraph_d_sqrt.Write()
    myfile.Close()
