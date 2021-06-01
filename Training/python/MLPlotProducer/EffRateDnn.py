#count average Taus
import ROOT
import numpy as np
import pandas as pd
import root_pandas
#pd.set_option('display.max_columns', 500)
#import matplotlib.pyplot as plt
import root_pandas
import sys
from array import array
from TauMLTools.Training.python.DataLoaderL2 import *
#from DataLoaderL2 import *
import tensorflow as tf
import awkward as ak
import numpy as np
from sklearn.metrics import classification_report, roc_auc_score,roc_curve, auc, accuracy_score, f1_score, confusion_matrix
import os
import argparse
import ROOT
import math
import scipy
import statsmodels.stats.proportion as ssp
ROOT.gStyle.SetPaintTextFormat(".2f")
parser = argparse.ArgumentParser(description='')
parser.add_argument('--machine', required=True, type=str, choices =['local', 'gridui'], default = 'local', help="input directory")
parser.add_argument('--effRate', required=True, type=str, choices =['eff', 'rate'], default = 'rate', help="input directory")
parser.add_argument('--verbose', required=False, type=int, default = 0)
parser.add_argument('--dropout', required=False, type=float, default = 0.2)
parser.add_argument('--num_den_layers', required=False, type=int, default = 5)
parser.add_argument('--learning_rate', required=False, type=float, default =0.002)
#parser.add_argument('--n-threads', required=False, type=int, default=1, help="number of threads")
args = parser.parse_args()

fileName = ""
outfileName = ""
datasetFile_name = ""
dataset_path = ""
absolute_path = ""
tupleName = 'L2TauTrainTuple'
trainDir_path = '/TauMLTools/Training/python/output/'
datasetDir_path= "/L2SkimmedTuples/DataSetTraining/"
flatVars = ['nVertices','l1Tau_pt', 'l1Tau_eta', 'l1Tau_hwIso']
vecVars = {
        'caloRecHit_e':'energy',
        'caloRecHit_had':'energy',
        'patatrack':'pt'
}


if(args.machine =='local'):
    absolute_path = "/Users/valeriadamante/Desktop/Dottorato/"
    dnn_path = absolute_path+'/cmssimphase2/'+trainDir_path
elif (args.machine == 'gridui'):
    absolute_path = "/home/users/damante/"
    dnn_path = absolute_path+'/CMSSW_11_2_1_Patatrack/src/'+trainDir_path

dataset_path = absolute_path+datasetDir_path
variables = []

if(args.effRate=='rate'):
    fileName = 'miniTuple_Data.root'
    outfileName = 'miniEvtTuple_Data.root'
    datasetFile_name = 'datasetData.npy'
    variables = {'evt','l1Tau_pt', 'l1Tau_hwIso'}
elif(args.effRate=='eff'):
    fileName = 'miniTuple_VBF.root'
    outfileName = 'miniEvtTuple_VBF.root'
    datasetFile_name = 'datasetVBF.npy'
    variables = {'evt','l1Tau_pt', 'l1Tau_hwIso', 'genLepton_vis_pt','genLepton_vis_eta','genLepton_isTau'}

#check if npy file with dataset exists

if not os.path.exists(("{}/{}").format(dataset_path, datasetFile_name)):
    awkArray = GetAwkArrayData(dataset_path+fileName, tupleName)
    featMatrix = GetFeatureMatrix(awkArray, flatVars, vecVars, False)
    np.save(("{}/{}").format(dataset_path, datasetFile_name), featMatrix)
else:
    if(args.verbose>0):
        print(("file {}/{} exists and taking data from there").format(dataset_path, datasetFile_name))
    featMatrix = np.load(('{}/{}').format(dataset_path, datasetFile_name))

params = {
    'activation_dense':'relu',
    'l1Pt_position':1,
    'batch_size':200,
    'train_fraction': 0.6102414433536747,
    'validation_fraction': 0.19474661713982488,
    'opt_threshold' : 0.2141452357154776,
    'bigOrRate': 13603.37
}

params['dropout_rate_den_layers']= args.dropout
params['learning_rate'] = args.learning_rate
params['num_den_layers'] = args.num_den_layers
file_suffix = "model_{}Layers_{:.2f}Dropout_{:.3f}LearningRate".format(params['num_den_layers'],params['dropout_rate_den_layers'],params['learning_rate']).replace(".","p")
# load model

model_path = dnn_path+ file_suffix +"/"+file_suffix
model = keras.models.load_model(model_path)
if(args.verbose>1):
    print(("model successfully loaded from {}").format(model_path))

# get Dataframe
df = root_pandas.read_root(dataset_path+fileName, tupleName, variables)
df ['y_predict'] = model.predict(featMatrix)

def get_n_evt(df_orig):
    df = df_orig.groupby("evt").count()
    df = df[df.l1Tau_pt >= 2]
    return df.shape[0]


if(args.effRate=='rate'):
    df_for_rate = df[(df.l1Tau_pt>=32) & ((df.l1Tau_hwIso>0) | (df.l1Tau_pt>=70))]
    df_for_den = df_for_rate
    df_for_num = df_for_rate[(df_for_rate.y_predict>params['opt_threshold'])]
    den = get_n_evt(df_for_den)
    num = get_n_evt(df_for_num)
    eff_rate = num*params['bigOrRate']/den
    eff = num/den

    c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
    c_low_rate = c_low*params['bigOrRate']
    c_up_rate = c_up*params['bigOrRate']

    print(("num \t {} \nden \t {} \neff \t {} \nunc_up \t {} \t unc_down \t {} \neff* hltL1sDoubleTauBigOR rate \t {} \nunc_up \t {} \t unc_down \t {}").format(num, den, eff, c_up, c_low, eff_rate, c_up_rate, c_low_rate))
    print(("eff-c_low \t {} \t c_up-eff \t {} \n(eff-c_low)*params['bigOrRate'] \t {} \t (c_up-eff)*params['bigOrRate'] \t {}").format(eff-c_low, c_up-eff, eff_rate-c_low_rate, c_up_rate-eff_rate))

    '''
    den = df.groupby('evt').count().shape[0]
    print(df[df.l1Tau_pt<2].shape[0])
    df_for_num = df[(df.l1Tau_pt>=32) & ((df.l1Tau_hwIso>0) | (df.l1Tau_pt>=70)) & (df.y_predict>params['opt_threshold'])].groupby('evt').count()
    print(df_for_num.head())
    num = df_for_num.shape[0]
    c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
    eff = num/den
    print(("num \t {} \nden \t {} \neff \t {} \nunc_up \t {} \t unc_down \t {} \neff* hltL1sDoubleTauBigOR rate \t {} \nunc_up \t {} \t unc_down \t {}").format(num, den, eff, c_low, c_up, eff*params['bigOrRate'], (c_low)*params['bigOrRate'],(c_up)*params['bigOrRate']))
    print(("eff-c_low \t {} \t c_up-eff \t {} \n(eff-c_low)*params['bigOrRate'] \t {} \t (c_up-eff)*params['bigOrRate'] \t {}").format(eff-c_low,c_up-eff, (eff-c_low)*params['bigOrRate'],(c_up-eff)*params['bigOrRate']))
    '''
pt_bins = [20,30,35,40,50,60,70,80,90,350]

def save_dataframe(df_orig, path_name, tuple_name):
    df_2 = df_orig.sort_values(by=['evt'])
    df_2_count = df_2.groupby('evt').size().reset_index(name='count')
    values = []
    prev_index = 0
    for i in range(0,df_2.shape[0]):
        for j in range(prev_index,df_2_count.shape[0]):
            temp_i = int(df_2.iloc[i]['evt'])
            temp_j = int(df_2_count.iloc[j]['evt'])
            if(temp_i == temp_j):
                values.append(df_2_count.iloc[j]['count'])
                prev_index = j
                break
    df_2['hasTwoHadTaus']=values
    root_pandas.to_root(df_2, path_name, tuple_name)

if(args.effRate=='eff'):
    if(not os.path.exists(dataset_path+'/df_for_den.root')):
        save_dataframe(df, dataset_path+'/df_for_den.root', tupleName)
    df_for_den = root_pandas.read_root(dataset_path+'/df_for_den.root', tupleName)

    if(not os.path.exists(dataset_path+'/df_for_num.root')):
        save_dataframe(df, dataset_path+'/df_for_num.root', tupleName)
    df_for_num = root_pandas.read_root(dataset_path+'/df_for_num.root', tupleName)

    df_for_den_count = df_for_den.groupby('evt').count()
    df_for_num_count = df_for_num.groupby('evt').count()
    den = df_for_den_count[(df_for_den_count.genLepton_vis_pt==2)].shape[0]
    num = df_for_num_count[(df_for_num_count.genLepton_vis_pt==2)].shape[0]
    print(num)
    print(den)

    eff = num/den
    c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
    print(("num \t {} \nden \t {} \neff \t {} \nunc_up \t {} \t unc_down \t {}").format(num, den, eff, c_up, c_low))
    print(("eff-c_low \t {} \t c_up-eff \t {} ").format(eff-c_low, c_up-eff))
    df_for_num_2 = df_for_num[df_for_num.hasTwoHadTaus==2]
    df_for_den_2 = df_for_den[df_for_den.hasTwoHadTaus==2]
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
    myfile = ROOT.TFile( "Efficiencies.root", 'RECREATE' )
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
    EffGraph.SetTitle("Efficiency; #tau 1 P_{T} (GeV); #tau 2 P_{T} (GeV)")
    EffGraph.Draw("TEXT2 COLZ")
    canvas1.Update()
    canvas1.Print("Efficiency_VBF.png", "png")

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
    canvas2.Print("Efficiency_VBF_Diag.png", "png")

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
    canvas3.Print("Efficiency_VBF_Diag_sqrt.png", "png")

    #save histograms

    EffGraph.Write()
    EffGraph_d.Write()
    EffGraph_d_sqrt.Write()
    myfile.Close()
