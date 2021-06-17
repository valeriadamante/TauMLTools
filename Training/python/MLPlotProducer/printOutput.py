from  TauMLTools.Training.python.produceGridDatasets import *
import tensorflow as tf
import root_pandas
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import scipy
import statsmodels.stats.proportion as ssp
import argparse
import ROOT
from array import array
ROOT.gStyle.SetPaintTextFormat(".2f")

parser = argparse.ArgumentParser()
parser.add_argument('--machine', required=False, type=str, default="local", choices=["local", "cmssimphase2"]) #aggiungi pccms65
parser.add_argument('--effRate', required= False, type=str, default = 'eff', choices=['test', 'eff','rate', 'Zprime'])
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
varDict = {
    0:"nVertices",
    1:"l1Tau_pt",
    2:"l1Tau_eta",
    3:"l1Tau_hwIso",
    4:"EcalEnergySum",
    5:"EcalSize",
    6:"EcalEnergyStdDev",
    7:"EcalDeltaEta",
    8:"EcalDeltaPhi",
    9:"EcalChi2",
    10:"EcalEnergySumForPositiveChi2",
    11:"EcalSizeForPositiveChi2",
    12:"HcalEnergySum",
    13:"HcalSize",
    14:"HcalEnergyStdDev",
    15:"HcalDeltaEta",
    16:"HcalDeltaPhi",
    17:"HcalChi2",
    18:"HcalEnergySumForPositiveChi2",
    19:"HcalSizeForPositiveChi2",
    20:"PatatrackPtSum",
    21:"PatatrackSize",
    22:"PatatrackSizeWithVertex",
    23:"PatatrackPtSumWithVertex",
    24:"PatatrackChargeSum",
    25:"PatatrackDeltaEta",
    26:"PatatrackDeltaPhi",
    27:"PatatrackChi2OverNdof",
    28:"PatatrackNdof",
    29:"PatatrackDxy",
    30:"PatatrackDz",
    }

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
#print(CellGridMatrix.shape)

#for i in range(len(CellGridMatrix)):
#    if(CellGridMatrix[i][0][0][1]==35.5/256.0):
#        print(i)
#        print(CellGridMatrix[i][0][0][1])
        #for j in range(len(CellGridMatrix[i][0][0])):
            #print(("var = {}").format(varDict[j]))
            #print(CellGridMatrix[i][0][0][:])
i = 26201
#print(CellGridMatrix[i][0][0][1]*256.)
'''
for phi_idx in range(len(CellGridMatrix[i])):
    for eta_idx in range(len(CellGridMatrix[i][phi_idx])):
        for j in range(len(CellGridMatrix[i][phi_idx][eta_idx])):
            long_value = CellGridMatrix[i][phi_idx][eta_idx][j]
            CellGridMatrix[i][phi_idx][eta_idx][j] = round(long_value,6)
            print(("var name = {} \t tau_idx = {} \t eta_idx = {} \t phi_idx = {} \tvalue = {}").format(varDict[j], i, eta_idx, phi_idx,CellGridMatrix[i][phi_idx][eta_idx][j]))
            #print(("{}, ").format(CellGridMatrix[i][phi_idx][eta_idx][j]))

'''
#print(CellGridMatrix)
# ****** Get DataFrames ******
#y_predict_test = model.predict(x_test).reshape(y_test.shape)
variables = ['evt','run','lumi','l1Tau_pt', 'l1Tau_hwIso']
if(args.effRate != 'rate'):
    variables = ['evt','run','lumi','l1Tau_pt', 'l1Tau_hwIso', 'genLepton_vis_pt','genLepton_vis_eta','genLepton_isTau']
inFile = GetRootPath(args.machine, args.effRate)
treeName = 'L2TauTrainTuple'

absolute_path=GetDataSetPath(args.machine)
dataFrameWithPredName = ("{}/dfWithPredictions_{}_3.root").format(absolute_path,GetNameForEffRate(args.effRate))

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

print(df[(df.lumi == 136)].head(10))
