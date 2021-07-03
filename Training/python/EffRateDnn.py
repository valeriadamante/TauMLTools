#count average Taus
import ROOT
import numpy as np
import pandas as pd
import root_pandas
#pd.set_option('display.max_columns', 500)
#import matplotlib.pyplot as plt
import root_pandas
import sys
#from TauMLTools.Training.DataLoaderL2 import *
from DataLoaderL2 import *
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

parser = argparse.ArgumentParser(description='')
parser.add_argument('--machine', required=False, type=str, choices =['local', 'gridui'], default = 'local', help="input directory")
parser.add_argument('--verbose', required=False, type=int, default = 0)
parser.add_argument('--dropout', required=False, type=float, default = 0.2)
parser.add_argument('--num_den_layers', required=False, type=int, default = 5)
parser.add_argument('--learning_rate', required=False, type=float, default =0.002)
#parser.add_argument('--n-threads', required=False, type=int, default=1, help="number of threads")
args = parser.parse_args()

absolute_path = ""
if(args.machine =='local'):
    absolute_path = "/Users/valeriadamante/Desktop/Dottorato/gridui/L2SkimmedTuples/DataSetTraining/"
elif (args.machine == 'gridui'):
    absolute_path = "/home/users/damante/L2SkimmedTuples/DataSetTraining/"
fileName = 'miniTuple_Data.root'
tupleName = 'L2TauTrainTuple'
trainDir_path = 'TauMLTools/Training/python/output/'
cmssimpath = '/Users/valeriadamante/Desktop/Dottorato/cmssimphase2/'
#cmssimpath = '/home/users/damante/L2SkimmedTuples/'
flatVars = ['nVertices','l1Tau_pt', 'l1Tau_eta', 'l1Tau_hwIso']
vecVars = {
        'caloRecHit_e':'energy',
        'caloRecHit_had':'energy',
        'patatrack':'pt'
}

if not os.path.exists(("{}/datasetData.npy").format(cmssimpath)):
    awkArray = GetAwkArrayData(absolute_path+fileName, tupleName)
    featMatrix = GetFeatureMatrix(awkArray, flatVars, vecVars, False)
    np.save(("{}/datasetData.npy").format(cmssimpath), featMatrix)
else:
    if(args.verbose>0):
        print(("file {}/datasetData.npy exists and taking data from there").format(cmssimpath))
    featMatrix = np.load(('{}/datasetData.npy').format(cmssimpath))

params = {
    'activation_dense':'relu',
    'num_units_den_layers':int(2*len(featMatrix[0])/0.8),
    'l1Pt_position':1,
    'batch_size':200,
    'train_fraction': 0.6102414433536747,
    'validation_fraction': 0.19474661713982488,
    'opt_threshold' : 0.10631134733375802,
    'bigOrRate': 13603.37
}
params['dropout_rate_den_layers']= args.dropout
params['learning_rate'] = args.learning_rate
params['num_den_layers'] = args.num_den_layers
file_suffix = "model_{}Layers_{:.2f}Dropout_{:.3f}LearningRate".format(params['num_den_layers'],params['dropout_rate_den_layers'],params['learning_rate']).replace(".","p")
model_path = cmssimpath+trainDir_path+ file_suffix +"/"+file_suffix

model = keras.models.load_model(model_path)
if(args.verbose>1):
    print(("model successfully loaded from {}").format(model_path))
df = root_pandas.read_root(absolute_path+fileName, tupleName, {'evt','l1Tau_pt', 'l1Tau_hwIso'})

df ['y_predict_data'] = model.predict(featMatrix)
den = df.groupby('evt').count().shape[0]
#print(df.head())
df_for_num = df[(df.l1Tau_pt>=32) & ((df.l1Tau_hwIso>0) | (df.l1Tau_pt>70)) & (df.y_predict_data>params['opt_threshold'])].groupby('evt').count()
num = df_for_num[df_for_num.l1Tau_pt<2].shape[0]
c_low, c_up = ssp.proportion_confint(num, den, alpha=1-0.68, method='beta')
eff = num/den
print(("num \t {} \nden \t {} \neff \t {} \nunc_up \t {} \t unc_down \t {} \neff* hltL1sDoubleTauBigOR rate \t {} \nunc_up \t {} \t unc_down \t {}").format(num, den, eff, c_low, c_up, eff*params['bigOrRate'], (c_low)*params['bigOrRate'],(c_up)*params['bigOrRate']))
print(("eff-c_low \t {} \t c_up-eff \t {} \n(eff-c_low)*params['bigOrRate'] \t {} \t (c_up-eff)*params['bigOrRate'] \t {}").format(eff-c_low,c_up-eff, (eff-c_low)*params['bigOrRate'],(c_up-eff)*params['bigOrRate']))
