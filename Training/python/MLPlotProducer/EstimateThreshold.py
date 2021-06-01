#Estimate DNN thresholds
import ROOT
import numpy as np
import pandas as pd
import root_pandas
#pd.set_option('display.max_columns', 500)
#import matplotlib.pyplot as plt
import root_pandas
import sys
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

parser = argparse.ArgumentParser(description='')
parser.add_argument('--machine', required=True, type=str, choices =['local', 'gridui'], default = 'local', help="input directory")
parser.add_argument('--verbose', required=False, type=int, default = 0)
parser.add_argument('--dropout', required=False, type=float, default = 0.2)
parser.add_argument('--num_den_layers', required=False, type=int, default = 5)
parser.add_argument('--learning_rate', required=False, type=float, default =0.002)
#parser.add_argument('--n-threads', required=False, type=int, default=1, help="number of threads")
args = parser.parse_args()

dnn_path = ""

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


variables = {'evt','l1Tau_pt', 'l1Tau_hwIso'}

fileName_data = 'miniTuple_Data.root'
outfileName_data = 'miniEvtTuple_Data.root'
datasetFile_name_data = 'datasetData.npy'

fileName_test = 'DataSetTrainingWeight.root'
outfileName_test = 'miniEvtTuple.root'
datasetFile_name_test ='dataset.npy'

#check if npy file with dataset exists
# load dataset for data
if not os.path.exists(("{}/{}").format(dataset_path, datasetFile_name_data)):
    awkArray_data = GetAwkArrayData(dataset_path+fileName_data, tupleName)
    featMatrix_data = GetFeatureMatrix(awkArray_data, flatVars, vecVars, False)
    np.save(("{}/{}").format(dataset_path, datasetFile_name_data), featMatrix_data)
else:
    if(args.verbose>0):
        print(("file {}/{} exists and taking data from there").format(dataset_path, datasetFile_name_data))
    featMatrix_data = np.load(('{}/{}').format(dataset_path, datasetFile_name_data))

# load dataset for test
# dataset - x
if not os.path.exists(("{}/{}").format(dataset_path, datasetFile_name_test)):
    awkArray_test = GetAwkArrayData(dataset_path+fileName_test, tupleName)
    featMatrix_test = GetFeatureMatrix(awkArray_test, flatVars, vecVars, False)
    np.save(("{}/{}").format(dataset_path, datasetFile_name_test), featMatrix_test)
else:
    if(args.verbose>0):
        print(("file {}/{} exists and taking data from there").format(dataset_path, datasetFile_name_test))
    featMatrix_test = np.load(('{}/{}').format(dataset_path, datasetFile_name_test))

# dataset - MCtruth
if not os.path.exists(('{}/MCTruth.npy').format(dataset_path)):
    awkArray_test = GetAwkArray(dataset_path+fileName_test, tupleName)
    MCTruth = GetMCTruth(awkArray_test, 'genLepton_isTau')
else:
    if(args.verbose>0):
        print(("file {}/MCTruth.npy exists and taking data from there").format(dataset_path))
    MCTruth = np.load(('{}/MCTruth.npy').format(dataset_path))



params = {
    'activation_dense':'relu',
    'num_units_den_layers':int(2*len(featMatrix_test[0])/0.8),
    'l1Pt_position':1,
    'batch_size':200,
    'train_fraction': 0.6102414433536747,
    'validation_fraction': 0.19474661713982488,
    'opt_threshold' : 0.30631134733375802,
    'bigOrRate': 13603.37
}

params['dropout_rate_den_layers']= args.dropout
params['learning_rate'] = args.learning_rate
params['num_den_layers'] = args.num_den_layers
file_suffix = "model_{}Layers_{:.2f}Dropout_{:.3f}LearningRate".format(params['num_den_layers'],params['dropout_rate_den_layers'],params['learning_rate']).replace(".","p")
# dataset - weight
if not os.path.exists(("{}/weights_test.npy").format(dataset_path)):
    awkArray_test = GetAwkArrayData(dataset_path+fileName_test, tupleName)
    w_train, w_test, w_val = GetTrainTestWeight(awkArray_test, params['train_fraction'],params['validation_fraction'])
    np.save(("{}/weights_test.npy").format(dataset_path), w_test)
else:
    if(args.verbose>0):
        print(("file {}/weights_test.npy exists and taking data from there").format(dataset_path))
    w_test = np.load(("{}/weights_test.npy").format(dataset_path))

x_train, y_train, x_test, y_test, x_val, y_val = GetTrainTestFraction(MCTruth, featMatrix_test, params['train_fraction'],
                                                                      params['validation_fraction'], False)

# load model

model_path = dnn_path+ file_suffix +"/"+file_suffix
model = keras.models.load_model(model_path)
if(args.verbose>1):
    print(("model successfully loaded from {}").format(model_path))

# get DataFrame for test
#df_test = root_pandas.read_root(dataset_path+fileName_test, tupleName, variables)
#df_test['y_predict'] = model.predict(featMatrix_test)
y_predict_test = model.predict(featMatrix_test)

# get Dataframe for data
df_data = root_pandas.read_root(dataset_path+fileName_data, tupleName, variables)
df_data['y_predict'] = model.predict(featMatrix_data)

# estimate the rate


def get_n_evt(df_orig):
    df = df_orig.groupby("evt").count()
    df = df[df.l1Tau_pt >= 2]
    return df.shape[0]

df_for_rate = df_data[(df_data.l1Tau_pt>=32) & ((df_data.l1Tau_hwIso>0) | (df_data.l1Tau_pt>=70))]
def rate_calculator(thr):
    df_for_den = df_for_rate
    df_for_num = df_for_rate[(df_for_rate.y_predict>thr)]
    den = get_n_evt(df_for_den)
    num = get_n_evt(df_for_num)
    eff = num*params['bigOrRate']/den
    return eff

target_rate = 5000
n_per_evt = math.floor(pd.DataFrame.mean(df_data['l1Tau_pt']))


def get_delta_rate(dnn_thr):
    return rate_calculator(dnn_thr) - target_rate
opt_dnn_thr = scipy.optimize.bisect(get_delta_rate, 0, 1)
print(opt_dnn_thr)
print(rate_calculator(opt_dnn_thr))
