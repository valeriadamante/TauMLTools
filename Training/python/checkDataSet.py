from  TauMLTools.Training.python.produceGridDatasets import *
import tensorflow as tf
import root_pandas
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--machine', required=False, type=str, default="local", choices=["local", "cmssimphase2"]) #aggiungi pccms65
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
print(CellGridMatrix[0:1,:,:,6])
inFile, outFile, outFileNorm, dictFileToRead, dictFileToSave = GetDictFile(args.machine, args.effRate, **kwArgs)

print(("max = {}").format(np.amax(CellGridMatrix)))
print(("nan = {}").format(np.count_nonzero(np.isnan(CellGridMatrix))))
print(CellGridMatrix.shape)
for n in range(CellGridMatrix.shape[3]):
  if np.count_nonzero(np.isnan(CellGridMatrix[:,:,:,n])) > 0:
    print(n)
dictFileToLoad = open(dictFileToRead)
print(("for {}").format(inFile))
dict = json.load(dictFileToLoad)
#for i in dict:
    #print(dict[i]['std'])
print(("for test"))
dictFileToLoad = open(dictFileToSave)
dict = json.load(dictFileToLoad)
#for i in dict:
    #print(dict[i]['std'])
