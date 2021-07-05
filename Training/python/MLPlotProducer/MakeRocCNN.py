import tensorflow as tf
import uproot
import awkward as ak
import pickle
import numpy as np
import json
import time
from  TauMLTools.Training.python.produceGridDatasets import *
import matplotlib
matplotlib.use('agg')
import pylab as plt
from sklearn.metrics import classification_report, roc_auc_score, accuracy_score, f1_score, confusion_matrix
import os
import argparse
import shutil

def plot_ROC(outDir, ns_fpr, ns_tpr, lr_fpr, lr_tpr, ts_fpr, ts_tpr, file_suffix):
    fig = plt.figure()
    # calculate roc curves
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
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    plot_directory =outDir+'/plots'
    if not os.path.exists(plot_directory):
        os.mkdir(plot_directory)
    print(("saving plot into {}").format(plot_directory))
    fig.savefig(plot_directory+'/ROC_val_'+file_suffix+'.pdf')
    fig.savefig(plot_directory+'/ROC_val_'+file_suffix+'.png')
    #model.save("output/model_"+file_suffix)

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser()
parser.add_argument('--machine', required=False, type=str, default="local", choices=["local", "cmssimphase2"]) #aggiungi pccms65
parser.add_argument('--effRate', required= False, type=str, default = 'test', choices=['test', 'eff','rate', 'Zprime'])
parser.add_argument('--n_max_events', required=False, type=int, default=-1, help='max number of events to be processed')
parser.add_argument('--McTruth_file', required=False, type=str, default='MCTruth.npy', help='output file name')
parser.add_argument('--Weights_file', required=False, type=str, default='weights.npy', help='output file name')
parser.add_argument('--n_cellsX', required=False, type=int, default=5, help='number of cells along X dir')
parser.add_argument('--n_cellsY', required=False, type=int, default=5, help='number of cells along Y dir')
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
    'train_fraction': 0.60003292,
    'validation_fraction': 0.19998354,
    #'learning_rate':0.002,
    'learning_rate':0.001,
    'monitor_es':'val_loss',
    'patience_es':10,
    'epochs':100000,
    'bigOrRate': 13603.37
}
# ***** Get cell grid Matrix *****
kwArgs = {'n_max_events':args.n_max_events, 'n_cellsX':args.n_cellsX, 'n_cellsY':args.n_cellsY, 'timeInfo' : True, 'verbose' : args.verbose}
CellGridMatrix = GetCellGridNormMatrix(args.machine, args.effRate, **kwArgs)
# ***** Get MC and weights Matrices *****
if(args.effRate== 'test'):
    MCTruth, weights = GetMCTruthWeights(args.machine, args.McTruth_file, args.Weights_file, **kwArgs)
    # ****** Get train - test - validation samples ********
    if(args.verbose)>0:
        print("preparing train/test/val samples")
    number_of_batches = len(MCTruth)/params['batch_size']
    x_train, y_train, w_train, x_test, y_test, w_test, x_val, y_val, w_val = GetTrainTestFraction(MCTruth, CellGridMatrix, weights, params['train_fraction'],  params['validation_fraction'], args.verbose)



# ******* parameter setting ******
# do =
# lr =
# nlayersdense =
# nlayerscnn2=
# nlayerscnn1=
# params['dropout_rate_den_layers']= do
# params['learning_rate'] = lr
# params['num_den_layers'] = ly


# ***** load model ******
model = tf.keras.models.load_model(GetModelPath(args.machine, params))
print("model successfully loaded")

best_auc_val = 0.
best_auc_test = 0.
best_acc_test = 0.
best_acc_val = 0.
best_auc_val_file = ""
best_auc_test_file = ""

from sklearn.metrics import roc_curve, auc
# calculate scores
ns_probs_test = [0 for _ in range(len(y_test))]
y_predict_test_nonshaped = model.predict(x_test)
print("prediction is done on test")
ns_probs_val = [0 for _ in range(len(y_val))]
y_predict_val_nonshaped = model.predict(x_val)
print("prediction is done on val")
y_predict_test = y_predict_test_nonshaped.reshape((y_test.shape))
y_predict_val  = y_predict_val_nonshaped.reshape((y_val.shape))
#y_pred_test_reshaped = y_predict_test[:,0,0,0]

ns_auc = roc_auc_score(y_val, ns_probs_val, sample_weight=w_val)
ts_auc = roc_auc_score(y_test, y_predict_test, sample_weight=w_test)
lr_auc = roc_auc_score(y_val, y_predict_val, sample_weight=w_val)

# summarize scores
print('No Skill: ROC AUC=%.3f' % (ns_auc))
print('L2Model test: ROC AUC=%.3f' % (ts_auc))
print('L2Model val: ROC AUC=%.3f' % (lr_auc))
y_predict_test.resize(len(y_test))
predictions_test = [round(value) for value in y_predict_test]
accuracy_test = accuracy_score(y_test, predictions_test)
print("Accuracy test: %.2f%%" % (accuracy_test * 100.0))

y_predict_val.resize(len(y_val))
predictions_val = [round(value) for value in y_predict_val]
accuracy_val = accuracy_score(y_val, predictions_val)
print("Accuracy val: %.2f%%" % (accuracy_val * 100.0))

ns_fpr, ns_tpr, ns_thresholds = roc_curve(y_val, ns_probs_val, sample_weight=w_val)
lr_fpr, lr_tpr, lr_thresholds = roc_curve(y_val, y_predict_val, sample_weight=w_val)
ts_fpr, ts_tpr, ts_thresholds = roc_curve(y_test, y_predict_test, sample_weight=w_test)
if(lr_auc > best_auc_val):
    best_auc_val = lr_auc
    best_auc_val_file = GetModelPath(args.machine, params)
    best_acc_val = accuracy_val
if(ts_auc > best_auc_test):
    best_auc_test =  ts_auc
    best_auc_test_file = GetModelPath(args.machine, params)
    best_acc_test = accuracy_test
outDir=GetOutPath(args.machine)
plot_ROC(outDir,ns_fpr, ns_tpr, lr_fpr, lr_tpr, ts_fpr, ts_tpr,GetFileSuffix(params))


print(("best AUC validation is : {:.3f} from the file {}").format(best_auc_val, best_auc_val_file) )
print(("best AUC test is : {:.3f} from the file {}").format(best_auc_test, best_auc_test_file) )
print(("best accuracy for validation is : {}").format(best_acc_val))
print(("best accuracy for test is : {}").format(best_acc_test))
'''
from tensorflow import keras
for ly in range(5,11):
    for lr in np.arange(0.001, 0.02, 0.001):
        for do in np.arange(0.2, 0.5, 0.1):
            #do = 0.3
            #ly = 5
            #lr = 0.001
            params['dropout_rate_den_layers']= do
            params['learning_rate'] = lr
            params['num_den_layers'] = ly
            file_suffix ="model_{}Layers_{:.2f}Dropout_{:.3f}LearningRate".format(params['num_den_layers'],params['dropout_rate_den_layers'],params['learning_rate']).replace(".","p")
            model_path = absolute_path+trainDir_path+ file_suffix +"/"+file_suffix
            if not os.path.exists(model_path):
                print(("{} does not exist").format(model_path))
                continue
            model = keras.models.load_model(model_path)
            print("model successfully loaded")

            from sklearn.metrics import roc_curve, auc
            # calculate scores
            ns_probs_test = [0 for _ in range(len(y_test))]
            y_predict_test = model.predict(x_test)
            print("prediction is done on test")
            ns_probs_val = [0 for _ in range(len(y_val))]
            y_predict_val = model.predict(x_val)
            print("prediction is done on val")
            ns_auc = roc_auc_score(y_val, ns_probs_val, sample_weight=w_val)
            ts_auc = roc_auc_score(y_test, y_predict_test, sample_weight=w_test)
            lr_auc = roc_auc_score(y_val, y_predict_val, sample_weight=w_val)

            # summarize scores
            print('No Skill: ROC AUC=%.3f' % (ns_auc))
            print('L2Model test: ROC AUC=%.3f' % (ts_auc))
            print('L2Model val: ROC AUC=%.3f' % (lr_auc))
            y_predict_test.resize(len(y_test))
            predictions_test = [round(value) for value in y_predict_test]
            accuracy_test = accuracy_score(y_test, predictions_test)
            print("Accuracy test: %.2f%%" % (accuracy_test * 100.0))

            y_predict_val.resize(len(y_val))
            predictions_val = [round(value) for value in y_predict_val]
            accuracy_val = accuracy_score(y_val, predictions_val)
            print("Accuracy val: %.2f%%" % (accuracy_val * 100.0))

            ns_fpr, ns_tpr, ns_thresholds = roc_curve(y_val, ns_probs_val, sample_weight=w_val)
            lr_fpr, lr_tpr, lr_thresholds = roc_curve(y_val, y_predict_val, sample_weight=w_val)
            ts_fpr, ts_tpr, ts_thresholds = roc_curve(y_test, y_predict_test, sample_weight=w_test)
            if(lr_auc > best_auc_val):
                best_auc_val = lr_auc
                best_auc_val_file = model_path
                best_acc_val = accuracy_val
            if(ts_auc > best_auc_test):
                best_auc_test =  ts_auc
                best_auc_test_file = model_path
                best_acc_test = accuracy_test
            plot_ROC(absolute_path+trainDir_path,ns_fpr, ns_tpr, lr_fpr, lr_tpr, ts_fpr, ts_tpr)


print(("best AUC validation is : {:.3f} from the file {}").format(best_auc_val, best_auc_val_file) )
print(("best AUC test is : {:.3f} from the file {}").format(best_auc_test, best_auc_test_file) )
print(("best accuracy for validation is : {}").format(best_acc_val))
print(("best accuracy for test is : {}").format(best_acc_test))

for ly in range(5,11):
    for lr in np.arange(0.001, 0.04, 0.01):
        for do in np.arange(0.2, 0.5, 0.1):
            #do = 0.3
            #ly = 5
            #lr = 0.001
            params['dropout_rate_den_layers']= do
            params['learning_rate'] = lr
            params['num_den_layers'] = ly
            file_suffix ="model_{}Layers_{:.2f}Dropout_{:.3f}LearningRate".format(params['num_den_layers'],params['dropout_rate_den_layers'],params['learning_rate']).replace(".","p")
            model_path = absolute_path+trainDir_path+ file_suffix +"/"+file_suffix
            if not os.path.exists(model_path):
                print(("{} does not exist").format(model_path))
                continue
            shutil.rmtree(model_path, ignore_errors=True)
'''
