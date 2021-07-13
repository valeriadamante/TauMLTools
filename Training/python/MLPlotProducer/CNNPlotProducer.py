from sklearn.metrics import classification_report, roc_auc_score,roc_curve, auc, accuracy_score, f1_score, confusion_matrix
from  TauMLTools.Training.python.produceGridDatasets import *
import tensorflow as tf
import root_pandas
import pandas as pd
import numpy as np
import scipy
import statsmodels.stats.proportion as ssp
import argparse
import ROOT
from array import array
ROOT.gStyle.SetPaintTextFormat(".2f")
parser = argparse.ArgumentParser()
parser.add_argument('--machine', required=False, type=str, default="local", choices=["local", "cmssimphase2"]) #aggiungi pccms65
parser.add_argument('--rateValue', required= False, type=int, default = 4, choices=[3,4,5])
parser.add_argument('--effRate', required= False, type=str, default = 'test', choices=['test', 'eff','rate'])
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
    'bigOrRate': 13603.37,


    'opt_threshold_3': 0.180858813224404,
    'opt_threshold_4': 0.12267940863785043,
    'opt_threshold_5': 0.08411243185219064,
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


# ***** create model directory *****

model = tf.keras.models.load_model(GetModelPath(args.machine, params))
print("model successfully loaded")


# ****** loss vs epoch for validation and training
def round_for_th(th, value):
    return 1 if value>=th else 0
from sklearn.metrics import roc_curve, auc
y_pred = model.predict(x_test)
y_pred.resize(len(y_test))
predictions = [round_for_th(params[('opt_threshold_{}').format(args.rateValue)], value) for value in y_pred]
#predictions = y_pred

# evaluate predictions
accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))
df = pd.read_csv(GetModelPath(args.machine, params)+".csv")
print(df.head())
# plot loss vs epoch for validation and training
fig = plt.figure(figsize=(5,5))
#plt.figure()
plt.plot(df['loss'], label='loss')
plt.plot(df['val_loss'], label='val_loss')
plt.legend(loc="upper right")
plt.title('loss VS epoch')
plt.xlabel('epoch')
plt.ylabel('loss')

fig.savefig(('{}/loss_plot_{}.pdf').format(GetRateDir(args.machine,args.rateValue),GetFileSuffix(params)))
fig.savefig(('{}/loss_plot_{}.png').format(GetRateDir(args.machine,args.rateValue),GetFileSuffix(params)))


# plot accuracy vs epoch for validation and training
fig = plt.figure(figsize=(5,5))
plt.plot(df['accuracy'], label='acc')
plt.plot(df['val_accuracy'], label='val_acc')
plt.legend(loc="lower left")
plt.title('accuracy VS epoch')
plt.xlabel('epoch')
plt.ylabel('acc')
fig.savefig(('{}/accuracy_plot_{}.pdf').format(GetRateDir(args.machine,args.rateValue),GetFileSuffix(params)))
fig.savefig(('{}/accuracy_plot_{}.png').format(GetRateDir(args.machine,args.rateValue),GetFileSuffix(params)))
# Plot ROC
# generate a no skill prediction (majority class)
ns_probs = [0 for _ in range(len(y_test))]
y_predict = model.predict(x_test)
#print(len(y_predict))
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

fig.savefig(('{}/ROC_onlyTest_{}.pdf').format(GetRateDir(args.machine,args.rateValue),GetFileSuffix(params)))
fig.savefig(('{}/ROC_onlyTest_{}.png').format(GetRateDir(args.machine,args.rateValue),GetFileSuffix(params)))
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
#for par in params:
#  print(par , "= ", params[par])

fig.savefig(('{}/compare_train_test_plot_{}.pdf').format(GetRateDir(args.machine,args.rateValue),GetFileSuffix(params)))
fig.savefig(('{}/compare_train_test_plot_{}.png').format(GetRateDir(args.machine,args.rateValue),GetFileSuffix(params)))


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

fig.savefig(('{}/ROC_TestVal_{}.pdf').format(GetRateDir(args.machine,args.rateValue),GetFileSuffix(params)))
fig.savefig(('{}/ROC_TestVal_{}.png').format(GetRateDir(args.machine,args.rateValue),GetFileSuffix(params)))
#model.save("output/model_"+GetFileSuffix(params))\
