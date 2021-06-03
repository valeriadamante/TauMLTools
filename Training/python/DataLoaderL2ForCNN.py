
import tensorflow as tf
import uproot
import awkward as ak
import numpy as np
#from  CellGridProducer import *
from  TauMLTools.Training.python.CellGridProducer import *
import json
from tensorflow import keras

def GetAwkArray(filePath, treeName):
    file = uproot.open(filePath)
    tree = file[treeName]
    entryStop = 200*3769
    return tree.arrays(entry_stop = entryStop, how="zip")

def GetAwkArrayData(filePath, treeName):
    file = uproot.open(filePath)
    tree = file[treeName]
    return tree.arrays(how="zip")

def GetMCTruth(awkArray, truthVar):
    nEvents = len(getattr(awkArray, truthVar))
    MCTruth = np.zeros(nEvents)
    MCTruth[:] = getattr(awkArray, truthVar)
    return MCTruth

def getNormCellGridMatrix( n_cellsX, n_cellsY, inFile, treeName, outFile, outFileNorm, dictFileToSave, dictFileToRead, distribFile, timeInfo, n_max_events, verbose, isTrainingDataSet, plot=True):
    start_0 = time.time()
    nVars = len(NNInputs)
    if(not os.path.exists(outFile) and not os.path.exists(outFileNorm)):
        tree = uproot.open(inFile)[treeName]
        start_1 = time.time()
        if(timeInfo):
            print (("time to get file {}").format(start_1-start_0))
        if(n_max_events!=-1):
            awkArray =tree.arrays(entry_stop=n_max_events)
        else:
            awkArray =tree.arrays()
        start_2 = time.time()
        if(timeInfo):
            print (("time to get awkarray {}").format(start_2-start_1))
        CellGrid = getCellGridMatrix(nVars, n_cellsX, n_cellsY, awkArray.nVertices,
                                awkArray.l1Tau_pt, awkArray.l1Tau_eta, awkArray.l1Tau_hwIso, awkArray.caloRecHit_e_energy, awkArray.caloRecHit_e_DeltaEta, awkArray.caloRecHit_e_DeltaPhi, awkArray.caloRecHit_e_chi2, awkArray.caloRecHit_had_energy, awkArray.caloRecHit_had_DeltaEta, awkArray.caloRecHit_had_DeltaPhi, awkArray.caloRecHit_had_chi2, awkArray.patatrack_pt, awkArray.patatrack_DeltaEta, awkArray.patatrack_DeltaPhi, awkArray.patatrack_chi2, awkArray.patatrack_ndof, awkArray.patatrack_charge, awkArray.patatrack_dxy, awkArray.patatrack_dz, awkArray.patatrack_hasVertex, awkArray.patatrack_vert_z, awkArray.patatrack_vert_ptv2, awkArray.patatrack_vert_chi2, awkArray.patatrack_vert_ndof)
        np.save(outFile,CellGrid)
        if(verbose>0):
            print(("cellGrid shape = {}").format(CellGrid.shape))
        start_3 = time.time()
        if(timeInfo):
            print (("time to define grid {}").format(start_3-start_2))
        dict={}
        if(isTrainingDataSet):
            StandardizeVars(CellGrid, dict, distribFile, verbose, timeInfo, plot)
        else:
            StandardizeVarsForDifferentDataSet(CellGrid, dict, distribFile, dictFileToRead, verbose, timeInfo, plot)
        with open(dictFileToSave, 'w') as fp:
            json.dump(dict, fp)
        np.save(outFileNorm, CellGrid)
    elif(os.path.exists(outFile) and not os.path.exists(outFileNorm)):
        CellGrid = np.load(outFile)
        start_1 = time.time()
        if(timeInfo):
            print (("time to get file {}").format(start_1-start_0))
        if(verbose>0):
            print(("cellGrid shape = {}").format(CellGrid.shape))
        dict={}
        if(isTrainingDataSet):
            StandardizeVars(CellGrid, dict, verbose, timeInfo)
        else:
            StandardizeVarsForDifferentDataSet(CellGrid, dict, distribFile, dictFileToRead, verbose=0, timeInfo=True, plot = True)
        with open(dictFileToSave, 'w') as fp:
            json.dump(dict, fp)
        np.save(outFileNorm, CellGrid)
    elif(os.path.exists(outFile) and os.path.exists(outFileNorm)):
        if(verbose)>0:
            print(("file {} exists and taking data from there").format(outFileNorm))
        CellGrid = np.load(outFileNorm)
    return CellGrid

def GetTrainTestFraction(MCTruth, CellMatrix, weights, trainFraction, valFraction, verbose):
    testFraction = 1 - (trainFraction+valFraction)
    nEvents = len(MCTruth)
    nTrain = int(round(nEvents * trainFraction, 0))
    nVal = nTrain+int(round(nEvents * valFraction, 0))
    nTest = nVal+int(round(nEvents * testFraction, 0))
    if verbose>0:
        print("train fraction = ", trainFraction, "\ntest fraction = ", testFraction)
        print("total events = ", nEvents, "\ntrain events = ", int(round(nEvents * trainFraction, 0)), "\ntest events = ",
                  int(round(nEvents * testFraction, 0)), "\nvalidation events = ", int(round(nEvents * valFraction, 0)))
        print("train+test+validation events = ",  nTest)
    x_train,y_train,w_train = CellMatrix[0:nTrain, :,:,:],MCTruth[0:nTrain], weights[0:nTrain]
    x_val, y_val,w_val  = CellMatrix[nTrain:nVal,:,:,:],MCTruth[nTrain:nVal], weights[nTrain:nVal]
    x_test,y_test,w_test = CellMatrix[nVal:nEvents,:,:,:],MCTruth[nVal:nEvents], weights[nVal:nEvents]
    return x_train,y_train,w_train,x_test,y_test,w_test,x_val,y_val,w_val


def GetTestSubSample(x_test,y_test,w_test, low_threshold, high_threshold, var_pos):
    n_events = int(len(y_test))
    n_vars = int(len(x_test[0]))
    print (("n events {} - n vars {}").format(n_events, n_vars))
    y_test_new = np.zeros( (n_events) )
    w_test_new = np.zeros( (n_events) )
    x_test_new = np.zeros( (n_events, n_vars)  )
    max_size = 0
    for i in range(0, len(y_test)):
        if(x_test[i][var_pos]>=low_threshold and x_test[i][var_pos]<=high_threshold):
            x_test_new[max_size][:]=x_test[i][:]
            y_test_new[max_size]=y_test[i]
            w_test_new[max_size]=w_test[i]
            max_size +=1
    print(max_size)
    y_test_bin = y_test_new[0:max_size]
    w_test_bin = w_test_new[0:max_size]
    x_test_bin = x_test_new[0:max_size][:]
    return x_test_bin, y_test_bin, w_test_bin

    testFraction = 1 - (trainFraction+valFraction)
    nEvents = len(MCTruth)
    nTrain = int(round(nEvents * trainFraction, 0))
    nVal = nTrain+int(round(nEvents * valFraction, 0))
    nTest = nVal+int(round(nEvents * testFraction, 0))
    if verbose:
        print("train fraction = ", trainFraction, "\ntest fraction = ", testFraction)
        print("total events = ", nEvents, "\ntrain events = ", nTrain, "\ntest events = ",
                  nTest, "\nvalidation events = ", nVal)
        print("train+test+validation events = ",  nTest)
    x_train,y_train = featMatrix[0:nTrain],MCTruth[0:nTrain]
    x_val,y_val  = featMatrix[nTrain:nVal],MCTruth[nTrain:nVal]
    x_test,y_test = featMatrix[nVal:nEvents],MCTruth[nVal:nEvents]
    return x_train,y_train,x_test,y_test,x_val,y_val




from tensorflow.keras.models import Model
class CNNModel(Model):
    def __init__(self, params):
        super(CNNModel, self).__init__()
        # Dense Block
        self.cnn_1x1 = []
        self.activation_cnn_1x1 = []
        self.dropout_cnn_1x1 = []
        self.batch_norm_cnn_1x1 = []
        self.cnn_2x2 = []
        self.activation_cnn_2x2 = []
        self.dropout_cnn_2x2 = []
        self.batch_norm_cnn_2x2 = []
        self.dense = []
        self.dropout_dense = []
        self.batch_norm_dense = []
        for i in range(params['num_CNN1x1_layers']):
            n_units_cnn1x1 = params[('nFilters_CNN1x1_{}').format(i)]
            # conv2d layer
            self.cnn_1x1.append(tf.keras.layers.Conv2D(filters=n_units_cnn1x1,
                                    kernel_size = (1,1), name='CNN1x1_{}'.format(i)))
            # batch norm
            self.batch_norm_cnn_1x1.append(tf.keras.layers.BatchNormalization(name='batch_normalization_CNN1x1_{}'.format(i)))
            # activation
            self.activation_cnn_1x1.append(tf.keras.layers.Activation(params['activation_CNN1x1']))
            # dropout
            if params['dropout_CNN1x1_layers'] > 0:
              self.dropout_cnn_1x1.append(tf.keras.layers.Dropout(params['dropout_CNN1x1_layers'],
                                                    name='dropout_CNN1x1_{}'.format(i)))
        for i in range(params['num_CNN2x2_layers']):
            n_units_cnn2x2 = params[('nFilters_CNN2x2_{}').format(i)]
            # conv2d layer
            self.cnn_2x2.append(tf.keras.layers.Conv2D(filters=n_units_cnn2x2,
                                    kernel_size = (2,2), name='CNN2x2_{}'.format(i)))
            # batch norm
            self.batch_norm_cnn_2x2.append(tf.keras.layers.BatchNormalization(name='batch_normalization_CNN2x2_{}'.format(i)))
            # activation
            self.activation_cnn_2x2.append(tf.keras.layers.Activation(params['activation_CNN2x2']))
            # dropout
            if params['dropout_CNN2x2_layers'] > 0:
              self.dropout_cnn_2x2.append(tf.keras.layers.Dropout(params['dropout_CNN2x2_layers'],
                                                    name='dropout_CNN2x2_{}'.format(i)))
        for i in range(params['num_dense_layers']):
            n_units_dense = params[('nFilters_dense_{}').format(i)]
            self.dense.append(tf.keras.layers.Dense(n_units_dense,
                                activation=params['activation_dense'],
                                name='dense_{}'.format(i)))
            self.batch_norm_dense.append(tf.keras.layers.BatchNormalization(name='batch_normalization{}'.format(i)))
            if params['dropout_dense_layers'] > 0:
                self.dropout_dense.append(tf.keras.layers.Dropout(params['dropout_dense_layers'],name='dropout_dense_{}'.format(i)))
        self.final_dense = tf.keras.layers.Dense(1, activation="sigmoid", name='output')
    @tf.function
    def call(self, x):
        for i in range(len(self.cnn_1x1)):
            x = self.cnn_1x1[i](x)
            if len(self.batch_norm_cnn_1x1) > i:
                x = self.batch_norm_cnn_1x1[i](x)
            if len(self.activation_cnn_1x1) > i:
                x = self.activation_cnn_1x1[i](x)
            if len(self.dropout_cnn_1x1) > i:
                x = self.dropout_cnn_1x1[i](x)
        for i in range(len(self.cnn_2x2)):
            x = self.cnn_2x2[i](x)
            if len(self.batch_norm_cnn_2x2) > i:
                x = self.batch_norm_cnn_2x2[i](x)
            if len(self.activation_cnn_2x2) > i:
                x = self.activation_cnn_2x2[i](x)
            if len(self.dropout_cnn_2x2) > i:
                x = self.dropout_cnn_2x2[i](x)
        shape_of_x = tf.shape(x)
        x = tf.reshape(x, (shape_of_x[0],shape_of_x[1]*shape_of_x[2]*shape_of_x[3]))
        for i in range(len(self.dense)):
            x = self.dense[i](x)
            if len(self.batch_norm_dense) > i:
                x = self.batch_norm_dense[i](x)
            if len(self.dropout_dense) > i:
                x = self.dropout_dense[i](x)
        x = self.final_dense(x)
        return x

from tensorflow.keras.callbacks import Callback
class SaveBestModel(Callback):
    def __init__(self, output, verbose=0, metric='val_loss'):
        self.output = output
        self.verbose= verbose
        self.metric = metric
        self.best = 1000000000.

    def on_epoch_end(self, epoch, logs=None):
        current = logs.get(self.metric)
        #print(current)
        if np.less(current, self.best):
                 if self.verbose > 0:
                   print(('\nEpoch {}: {} improved from {:.5f} to {:.5f}, saving model to {}').format(epoch + 1, self.metric, self.best, current, self.output))
                 self.model.save(self.output)
                 self.best = current


def compare_train_test(clf, X_train, y_train, X_test, y_test):
    decisions = []
    for X,y in ((X_train, y_train), (X_test, y_test)):
        d1 = clf.predict(X[y>0.5]).ravel()
        d2 = clf.predict(X[y<0.5]).ravel()
        decisions += [d1, d2]
    low = min(np.min(d) for d in decisions)
    high = max(np.max(d) for d in decisions)
    low_high = (low,high)
    return decisions,low_high

def computeArrayInVarBins(variable_to_cut, low_threshold, high_threshold, awkArray, verbose):
    if(verbose>0):
        print(('pt min = {}, pt_max = {}').format(pt_low_threshold, pt_high_threshold))
    var_cond0 = getattr(awkArray, variable_to_cut)
    # apply first condition
    if(verbose>0):
        print('applying first condition')
    condition1 = (var_cond0)>=low_threshold
    var_cond1 = var_cond0[condition1]
    awkArray_cond1 = awkArray[condition1]
    if(verbose>0):
        print('applying first condition')
    condition2 = (var_cond1)<=high_threshold
    var_cond2 = var_cond0[condition2]
    awkArray_cond2 = awkArray[condition2]
    print("after no cond")
    print(len(getattr(awkArray, variable_to_cut)))
    print("after 1st cond")
    print(len(getattr(awkArray, variable_to_cut)))
    print("after 2nd cond")
    print(len(getattr(awkArray, variable_to_cut)))
    return awkArray_cond2
