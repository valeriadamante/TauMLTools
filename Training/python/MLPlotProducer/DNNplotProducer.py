import tensorflow as tf
import awkward as ak
from TauMLTools.Training.python.DataLoaderL2 import *
import numpy as np
from sklearn.metrics import classification_report, roc_auc_score,roc_curve, auc, accuracy_score, f1_score, confusion_matrix
import os
import matplotlib.pyplot as plt
import argparse

import pandas as pd
import ROOT
import math
import scipy
import statsmodels.stats.proportion as ssp
ROOT.gStyle.SetPaintTextFormat(".2f")
parser = argparse.ArgumentParser(description='')
parser.add_argument('--machine', required=True, type=str, choices =['local', 'gridui'], default = 'local', help="input directory")
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

params = {
    'activation_dense':'relu',
    'l1Pt_position':1,
    'batch_size':200,
    'train_fraction': 0.6102414433536747,
    'validation_fraction': 0.19474661713982488,
    'opt_threshold' : 0.2141452357154776,
    'bigOrRate': 13603.37
}

if(args.machine =='local'):
    absolute_path = "/Users/valeriadamante/Desktop/Dottorato/"
    dnn_path = absolute_path+'/cmssimphase2/'+trainDir_path
elif (args.machine == 'gridui'):
    absolute_path = "/home/users/damante/"
    dnn_path = absolute_path+'/CMSSW_11_2_1_Patatrack/src/'+trainDir_path

dataset_path = absolute_path+datasetDir_path
variables = []

fileName = 'DataSetTrainingWeight.root'
outfileName = 'miniEvtTuple.root'
datasetFile_name ='dataset.npy'

print(dataset_path)

# load dataset for test
# dataset - x
if not os.path.exists(("{}/{}").format(dataset_path, datasetFile_name)):
    awkArray = GetAwkArrayData(dataset_path+fileName, tupleName)
    featMatrix = GetFeatureMatrix(awkArray, flatVars, vecVars, False)
    np.save(("{}/{}").format(dataset_path, datasetFile_name), featMatrix)
else:
    if(args.verbose>0):
        print(("file {}/{} exists and taking data from there").format(dataset_path, datasetFile_name))
    featMatrix = np.load(('{}/{}').format(dataset_path, datasetFile_name))

# dataset - MCtruth
if not os.path.exists(('{}/MCTruth.npy').format(dataset_path)):
    awkArray = GetAwkArray(dataset_path+fileName, tupleName)
    MCTruth = GetMCTruth(awkArray, 'genLepton_isTau')
else:
    if(args.verbose>0):
        print(("file {}/MCTruth.npy exists and taking data from there").format(dataset_path))
    MCTruth = np.load(('{}/MCTruth.npy').format(dataset_path))
x_train, y_train, x_test, y_test, x_val, y_val = GetTrainTestFraction(MCTruth, featMatrix, params['train_fraction'],
                                                                      params['validation_fraction'], False)
if not os.path.exists(("{}/weights_test1.npy").format(dataset_path)):
    awkArray = GetAwkArray(dataset_path+fileName, tupleName)
    w_train, w_test, w_val = GetTrainTestWeight(awkArray, params['train_fraction'],params['validation_fraction'])
    np.save(("{}/weights_test1.npy").format(dataset_path), w_test)
else:
    if(args.verbose>0):
        print(("file {}/weights_test.npy exists and taking data from there").format(dataset_path))
    w_test = np.load(('{}/weights_test1.npy').format(dataset_path))

if not os.path.exists(("{}/weights_val1.npy").format(dataset_path)):
    awkArray = GetAwkArray(dataset_path+fileName, tupleName)
    w_train, w_test, w_val = GetTrainTestWeight(awkArray, params['train_fraction'],params['validation_fraction'])
    np.save(("{}/weights_val1.npy").format(dataset_path), w_val)
else:
    if(args.verbose>0):
        print(("file {}/weights_val.npy exists and taking data from there").format(dataset_path))
    w_val = np.load(('{}/weights_val1.npy').format(dataset_path))
print(len(w_test))

params['dropout_rate_den_layers']= args.dropout
params['learning_rate'] = args.learning_rate
params['num_den_layers'] = args.num_den_layers
file_suffix = "model_{}Layers_{:.2f}Dropout_{:.3f}LearningRate".format(params['num_den_layers'],params['dropout_rate_den_layers'],params['learning_rate']).replace(".","p")
# load model

model_path = dnn_path+ file_suffix +"/"+file_suffix
model = keras.models.load_model(model_path)
if(args.verbose>1):
    print(("model successfully loaded from {}").format(model_path))

y_pred = model.predict(x_test)
y_pred.resize(len(y_test))
predictions = [round(value) for value in y_pred]
# evaluate predictions
accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))
df = pd.read_csv(dnn_path+ file_suffix +"/"+file_suffix+".csv")
#get_ipython().run_line_magic('matplotlib', 'inline')

# plot loss vs epoch for validation and training
fig = plt.figure(figsize=(5,5))
#plt.figure()
plt.plot(df['loss'], label='loss')
plt.plot(df['val_loss'], label='val_loss')
plt.legend(loc="upper right")
plt.title('loss VS epoch')
plt.xlabel('epoch')
plt.ylabel('loss')

fig.savefig('loss_plot.pdf')
fig.savefig('loss_plot.png')

# plot accuracy vs epoch for validation and training
fig = plt.figure(figsize=(5,5))
plt.plot(df['accuracy'], label='acc')
plt.plot(df['val_accuracy'], label='val_acc')
plt.legend(loc="lower left")
plt.title('accuracy VS epoch')
plt.xlabel('epoch')
plt.ylabel('acc')
fig.savefig('accuracy_plot.pdf')
fig.savefig('accuracy_plot.png')
# Plot ROC
# generate a no skill prediction (majority class)
ns_probs = [0 for _ in range(len(y_test))]
y_predict = model.predict(x_test)
print(len(y_predict))
from sklearn.metrics import roc_curve, auc
fpr, tpr, thresholds = roc_curve(y_test, y_predict,sample_weight=w_test)
#print(type(fpr))
#print(fpr.shape)
#print((fpr))
#print(type(tpr))
#print(tpr.shape)
#print((tpr))
#print(type(thresholds))
#print(thresholds.shape)
#print((thresholds))
roc_auc = auc(fpr, tpr)
fig = plt.figure(figsize=(5,5))
plt.plot(fpr, tpr, lw=2, ms=2.1, color='cyan', label='auc = %.3f' % (roc_auc))
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, ms=2.1, color='k', label='random chance')
plt.xlim([0, 1.0])
plt.ylim([0, 1.0])
plt.xlabel('false positive rate')
plt.ylabel('true positive rate')
plt.title('receiver operating curve')
plt.legend(loc="lower right")
fig.savefig('roc_onlyTest_plot.pdf')
fig.savefig('roc_onlyTest_plot.png')

# plot distribution signal and background for train and test
decisions, low_high=compare_train_test(model, x_train, y_train, x_test, y_test)

fig = plt.figure(figsize=(5,5))

bins=30
plt.hist(decisions[0], color='b', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True,
  label='S (train)')
plt.hist(decisions[1], color='r', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True,
  label='B (train)')
hist, bins = np.histogram(decisions[2], bins=bins, range=low_high, density=True)
scale = len(decisions[2]) / sum(hist)
err = np.sqrt(hist * scale) / scale
width = (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.errorbar(center, hist, yerr=err, fmt='o', c='b', label='S (test)')

hist, bins = np.histogram(decisions[3], bins=bins, range=low_high, density=True)
scale = len(decisions[3]) / sum(hist)
err = np.sqrt(hist * scale) / scale

plt.errorbar(center, hist, yerr=err, fmt='o', c='r', label='B (test)')
plt.title('S,B distribution')
plt.xlabel("Output")
plt.ylabel("Arbitrary units")
plt.legend(loc='best')
#plt.show()
for par in params:
  print(par , "= ", params[par])

fig.savefig('compare_train_test_plot.pdf')
fig.savefig('compare_train_test_plot.png')

fig = plt.figure()
from sklearn.metrics import roc_curve, auc
# calculate scores
ns_probs_test = [0 for _ in range(len(y_test))]
y_predict_test = model.predict(x_test)

ns_probs_val = [0 for _ in range(len(y_val))]
y_predict_val = model.predict(x_val)

ns_auc = roc_auc_score(y_val, ns_probs_val, sample_weight=w_val)
ts_auc = roc_auc_score(y_test, y_predict_test, sample_weight=w_test)
lr_auc = roc_auc_score(y_val, y_predict_val, sample_weight=w_val)

# summarize scores
print('No Skill: ROC AUC=%.3f' % (ns_auc))
print('L2Model test: ROC AUC=%.3f' % (ts_auc))
print('L2Model val: ROC AUC=%.3f' % (lr_auc))
# calculate roc curves
ns_fpr, ns_tpr, ns_thresholds = roc_curve(y_val, ns_probs_val, sample_weight=w_val)
lr_fpr, lr_tpr, lr_thresholds = roc_curve(y_val, y_predict_val, sample_weight=w_val)
ts_fpr, ts_tpr, ts_thresholds = roc_curve(y_test, y_predict_test, sample_weight=w_test)

# plot the roc curve for the model
plt.plot(ts_fpr, ts_tpr, color = 'blue', label=('Test auc = {:.3f}').format(ts_auc))
plt.plot(ns_fpr, ns_tpr, linestyle='--', label=('No Skill  auc = {:.3f}').format(ns_auc))
plt.plot(lr_fpr, lr_tpr, color='cyan',  label=('Val auc = {:.3f}').format(lr_auc))

# axis labels
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
# show the legend
plt.legend()
# show the plot
fig.savefig('ROC_val_'+file_suffix+'.pdf')
fig.savefig('ROC_val_'+file_suffix+'.png')
#model.save("output/model_"+file_suffix)\
