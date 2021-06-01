
from  TauMLTools.Training.python.produceGridDatasets import *
import tensorflow as tf
import ROOT
import pandas as pd
import root_pandas
import argparse
import scipy
import statsmodels.stats.proportion as ssp
parser = argparse.ArgumentParser()
parser.add_argument('--machine', required=False, type=str, default="local", choices=["local", "cmssimphase2"]) #aggiungi pccms65
parser.add_argument('--effRate', required= False, type=str, default = 'rate', choices=['test', 'eff','rate'])
parser.add_argument('--rateValue', required= False, type=int, default = 4, choices=[3,4,5])
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
    'train_fraction': 0.6102414433536747,
    'validation_fraction': 0.19474661713982488,
    #'learning_rate':0.002,
    'learning_rate':0.001,
    'monitor_es':'val_loss',
    'patience_es':10,
    'epochs':100000,
    'bigOrRate': 13603.37
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

# ****** Get DataFrames ******
#y_predict_test = model.predict(x_test).reshape(y_test.shape)
variables = ['evt','l1Tau_pt', 'l1Tau_hwIso']
if(args.effRate != 'rate'):
    variables = ['evt','l1Tau_pt', 'l1Tau_hwIso', 'genLepton_vis_pt','genLepton_vis_eta','genLepton_isTau']

inFile = GetRootPath(args.machine, args.effRate)
treeName = 'L2TauTrainTuple'
df = root_pandas.read_root(inFile, treeName, variables)
df['y_predict'] = model.predict(CellGridMatrix).reshape(df['evt'].shape)
if(args.verbose>0):
    print("prediction done, now evaluating the threshold")
# ***** Estimate the rate ******
def get_n_evt(df_orig):
    df = df_orig.groupby("evt").count()
    df = df[df.l1Tau_pt >= 2]
    return df.shape[0]

df_for_rate = df[(df.l1Tau_pt>=32) & ((df.l1Tau_hwIso>0) | (df.l1Tau_pt>=70))]
def rate_calculator(thr):
    df_for_den = df_for_rate
    df_for_num = df_for_rate[(df_for_rate.y_predict>thr)]
    den = get_n_evt(df_for_den)
    num = get_n_evt(df_for_num)
    eff = num*params['bigOrRate']/den
    return eff

target_rate = args.rateValue * 1000
n_per_evt = math.floor(pd.DataFrame.mean(df_for_rate['l1Tau_pt']))


def get_delta_rate(dnn_thr):
    return rate_calculator(dnn_thr) - target_rate
opt_dnn_thr = scipy.optimize.bisect(get_delta_rate, 0, 1)

print(("Treshold = {}").format(opt_dnn_thr))
print(("value of rate at that threshold = {}").format(rate_calculator(opt_dnn_thr)))
