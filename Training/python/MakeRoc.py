import tensorflow as tf
import uproot
import awkward as ak
import pickle
import numpy as np
import json
from DataLoaderL2 import *
import matplotlib
matplotlib.use('agg')
import pylab as plt
from sklearn.metrics import classification_report, roc_auc_score, accuracy_score, f1_score, confusion_matrix
import os
import argparse
import shutil
parser = argparse.ArgumentParser(description='Choose best ROC.')
parser.add_argument('--machine', required=False, type=str, choices =['local', 'cmssimphase2'], default = 'local', help="input directory")
#parser.add_argument('--output', required=True, type=str, help="output directory")
#parser.add_argument('--n-threads', required=False, type=int, default=1, help="number of threads")
args = parser.parse_args()


fileName = 'DataSetTrainingWeight.root'
tupleName = 'L2TauTrainTuple'
trainDir_path = 'TauMLTools/Training/python/output/'
absolute_path = '/Users/valeriadamante/Desktop/Dottorato/cmssimphase2/'
outDir_name = '/output4'

if(args.machine == 'cmssimphase2'):
    absolute_path = '/home/valeria/'
    outDir_name = '/output3'
#trainDir_path = 'TauMLTools/Training/python/output/'

def plot_ROC(dir_path, ns_fpr, ns_tpr, lr_fpr, lr_tpr, ts_fpr, ts_tpr):
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
    outdir = dir_path+"/"+outDir_name
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    plot_directory =outdir+'/plots'
    if not os.path.exists(plot_directory):
        os.mkdir(plot_directory)
    print(("saving plot into {}").format(plot_directory))
    fig.savefig(plot_directory+'/ROC_val_'+file_suffix+'.pdf')
    fig.savefig(plot_directory+'/ROC_val_'+file_suffix+'.png')
    #model.save("output/model_"+file_suffix)

# create all vars array, mctruth array and feature matrix
awkArray = GetAwkArray(absolute_path+fileName, tupleName)
if not os.path.exists(absolute_path+'MCTruth.npy'):
    print("file "+absolute_path+"MCTruth.npy exists and taking data from there")
    MCTruth = GetMCTruth(awkArray, 'genLepton_isTau')
else:
    MCTruth = np.load(absolute_path+'MCTruth.npy')
if not os.path.exists(absolute_path+'FeatMatrix.npy'):
    featMatrix = GetFeatureMatrix(awkArray, flatVars, vecVars, False)
else:
    print("file "+absolute_path+"FeatMatrix.npy exists and taking data from there")
    featMatrix = np.load(absolute_path+'FeatMatrix.npy')

params = {
    'activation_dense':'relu',
    #'num_units_den_layers':int(2*len(featMatrix[0])/0.8),
    #'dropout_rate_den_layers':0.2,
    'batch_size':200,
    'train_fraction': 0.6102414433536747,
    'validation_fraction': 0.19474661713982488,
    #'learning_rate':0.002,
    'monitor_es':'val_loss',
    'patience_es':10,
    'epochs':100000,
}
# #print the array
# print(featMatrix)
# print(MCTruth)

# define train, test and validation samples
print("preparing train/test/val samples")
number_of_batches = len(MCTruth)/params['batch_size']
x_train, y_train, x_test, y_test, x_val, y_val = GetTrainTestFraction(MCTruth, featMatrix, params['train_fraction'],
                                                                      params['validation_fraction'], False)
w_train, w_test, w_val = GetTrainTestWeight(awkArray, params['train_fraction'],params['validation_fraction'])


best_auc_val = 0.
best_auc_test = 0.
best_acc_test = 0.
best_acc_val = 0.
best_auc_val_file = ""
best_auc_test_file = ""

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
'''
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
