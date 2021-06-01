# algo efficiency for test dataset to compare with sqrt(diagonal) of event cutbased efficiency
# how to build this algorithm:

# 1. goal = compare algo efficiency per taus (NOT per events) on testing dataset with sqrt(diagonal) of cut based (algo?) efficiency (on VBF) per event (NOT per taus) in pt bins

from TauMLTools.Training.python.DataLoaderL2 import *
import tensorflow as tf
import awkward as ak
import numpy as np
from sklearn.metrics import classification_report, roc_auc_score,roc_curve, auc, accuracy_score, f1_score, confusion_matrix
import os
import argparse
import ROOT
import math
import scipy

parser = argparse.ArgumentParser(description='')

parser = argparse.ArgumentParser(description='')
parser.add_argument('--machine', required=True, type=str, choices =['local', 'gridui'], default = 'local', help="input directory")
parser.add_argument('--verbose', required=False, type=int, default = 0)
parser.add_argument('--dropout', required=False, type=float, default = 0.2)
parser.add_argument('--num_den_layers', required=False, type=int, default = 5)
parser.add_argument('--learning_rate', required=False, type=float, default =0.002)
#parser.add_argument('--n-threads', required=False, type=int, default=1, help="number of threads")
args = parser.parse_args()



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
    'l1Pt_position':1,
    'batch_size':200,
    'train_fraction': 0.6102414433536747,
    'validation_fraction': 0.19474661713982488
}
params['dropout_rate_den_layers']= args.dropout
params['learning_rate'] = args.learning_rate
params['num_den_layers'] = args.num_den_layers

file_suffix = "model_{}Layers_{:.2f}Dropout_{:.3f}LearningRate".format(params['num_den_layers'],params['dropout_rate_den_layers'],params['learning_rate']).replace(".","p")

model_path = absolute_path+trainDir_path+ file_suffix +"/"+file_suffix

pts = [20,30,35,40,50,60,70,80,90,350]

numerator = 0.
denominator = 0.
numerator_fpr = 0.
denominator_fpr = 0.
algo_eff = 0.
th = 0.
eff = 0.
my_fpr = 0.
fpr_atTrheshold =0.
zero_order=0.
#  create all vars array, mctruth array and feature matrix
if(args.verbose>1):
    print("obtained featMatrix and MCTruth, now computing bin subsamples ")
x_train, y_train, x_test, y_test, x_val, y_val = GetTrainTestFraction(MCTruth, featMatrix_test, params['train_fraction'],
                                                                      params['validation_fraction'], False)

if not os.path.exists(("{}/weights_test.npy").format(absolute_path)):
    awkArray = GetAwkArray(absolute_path+fileName, tupleName)
    w_train, w_test, w_val = GetTrainTestWeight(awkArray, params['train_fraction'],params['validation_fraction'])
    np.save("weights_test.npy", w_test)
else:
    if(args.verbose>0):
        print(("file {}/weights_test.npy exists and taking data from there").format(absolute_path))
    w_test = np.load(('{}/weights_test.npy').format(absolute_path))

if(args.verbose>1):
    print(("loading model from {}").format(model_path))

if not os.path.exists(model_path):
    print(("{} does not exist").format(model_path))

model = keras.models.load_model(model_path)
if(args.verbose>1):
    print(("model successfully loaded from {}").format(model_path))

y_predict_test = model.predict(x_test)
y_predict_test.resize(len(y_test))
# 2.3 evaluate algo efficiency

#numerator = 0
#denominator = 0
#th = 0.

fpr, tpr, thresholds = roc_curve(y_test, y_predict_test,sample_weight=w_test)


rate = 5000/13603.37
print(math.floor(2.0686249073387692))
n_per_evt = math.floor(2.0686249073387692)
def f(fpr):
    return 1 - scipy.stats.binom.cdf(1, n_per_evt, fpr) - rate
opt_fpr = scipy.optimize.bisect(f, 0.01, 0.99)
weights_fake = w_test[y_test==0]
dnn_thr_fake = y_predict_test[y_test==0]
den = np.sum(weights_fake)
def f2(dnn_thr):
    fpr = np.sum(weights_fake[dnn_thr_fake>dnn_thr])/den
    return fpr-opt_fpr
opt_dnn_thr= scipy.optimize.bisect(f2, 0., 1.)
print(opt_dnn_thr)
#print(scipy.optimize.bisect(f, 0.01, 0.99))

'''
prev_fpr = 0.
for i in range(0, len(thresholds)):
    tpr_check = round(tpr[i], 2)
    fpr_check = round(fpr[i], 8)
    opt_fpr_check = round(opt_fpr,8)
    #if(i%10000==0):
    #    print(("i = {}, tpr_check = {}, f(fpr_check) = {}, fpr_check = {}").format(i, tpr[i], opt_fpr_check, fpr[i]))
    if(fpr[i]>=opt_fpr and opt_fpr > prev_fpr):
        print(("sono simili in questa finestra : {}, e valgono {} {} ").format(i,fpr[i], opt_fpr))
        print(( "sono uguali e valgono {} {}").format(fpr_check,opt_fpr_check))

    prev_fpr=fpr[i]
    #if(tpr_check>=0.98 ):
    #    print(("tpr = {} \nthreshold = {} \nfpr = {}").format(tpr_check, thresholds[i], fpr[i]))
    #    th = thresholds[i]
    #    fpr_atTrheshold = fpr[i]
    #    zero_order=f(fpr_check)
    #    break
#print(("f(tpr) for tpr {} = {}").format(0.98, zero_order))
for i in range(0, len(y_predict_test)):
    if y_test[i] == 0:
        denominator_fpr += w_test[i]
        if y_predict_test[i]>th:
            numerator_fpr +=w_test[i]
my_fpr = numerator_fpr/denominator_fpr
print(("fpr with weights = {}").format(my_fpr))

for i in range(0, len(y_predict_test)):
    if y_test[i] == 1:
        denominator += 1
        if y_predict_test[i]>th:
            numerator += 1
algo_eff = numerator / denominator
print(("algo efficiency is {}").format(algo_eff))
'''
