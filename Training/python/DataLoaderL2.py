
import tensorflow as tf
import uproot
import awkward as ak
import numpy as np
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
    np.save('/home/valeria/MCTruth.npy', MCTruth)
    return MCTruth

vecVars = {
    # class : [variable to get size and sum, variables to make weightedSum]
    "Patatracks " : []
}
def getCellGridMatrix(awkArray, flatVars, vecVars, n_cellsX, n_cellsY, verbose=False):
    print("writing feature matrix")
    nVars = 27
    nEvents = len(getattr(awkArray, flatVars[0]))
    dPhi_min = -0.5
    dPhi_max = 0.5
    dEta_min = -0.5
    dEta_max = 0.5
    dPhi_witdh = (dPhi_max-dPhi_min)/n_cellsX
    dEta_witdh = (dEta_max-dEta_min)/n_cellsY
    nEvents = len(variableToAdd)
    CellGrid = np.zeros((nEvents, n_cellsX, n_cellsY, n_vars))
    for i in range(0,n_cellsX):
        dEta_current = dEta_min + i * dEta_witdh
        for j in range(0, n_cellsY):
            dPhi_current = dPhi_min + j * dPhi_witdh
            if(verbose):
                print(('dEta  = {}, dPhi = {}').format(dEta_current, dPhi_current))
            # define variable index
            varIndex = 0
            for var in flatVars:
                # 1. copy global variables : nVertex - l1Tau_pt - l1Tau_eta - l1Tau_phi - l1Tau_hwIso
                CellGrid[:, i, j, varIndex] = getattr(awkArray, var)
            varIndex = varIndex + len(flatVars)
            for var in vecVars:
                dPhi_cond0 = getattr(getattr(awkArray, var),'DeltaPhi')
                dEta_cond0 = getattr(getattr(awkArray, var),'DeltaEta')
                # 2. sum of objects within cell
                VarToSumSize = getattr(getattr(awkArray, var),var[0])

                if(verbose):
                   print('applying first condition')
                cond1_dPhi = (dPhi_cond0)>=dPhi_current
                dPhi_cond1 = dPhi_cond0[cond1_dPhi]
                VarToSumSize_condPhi1=VarToSumSize[cond1_dPhi]

                cond1_dEta = (dEta_cond0)>=dEta_current
                dEta_cond1 = dEta_cond0[cond1_dEta]
                VarToSumSize_condEta1=VarToSumSize_condPhi1[cond1_dEta]

                if(verbose):
                   print('applying second condition')
                cond2_dPhi = (dPhi_cond0)<dPhi_current+dPhi_width
                dPhi_cond2 = dPhi_cond0[cond2_dPhi]
                VarToSumSize_condPhi2=VarToSumSize_condEta1[cond2_dPhi]

                cond2_dEta = (dEta_cond0)<dEta_current+dEta_width
                dEta_cond2 = dEta_cond0[cond2_dEta]
                VarToSumSize_condEta2=VarToSumSize_condPhi2[cond2_dEta]

                # get the sum and the size of each vector
                if(verbose):
                    print('defining sum and size arrays')
                arr_sum = np.zeros(nEvents)
                arr_size = np.zeros(nEvents)
                if(verbose):
                    print('filling sum and size arrays')

                for v in range(0, len(VarToSumSize_condEta2)-1):
                    sumen = ak.sum(VarToSumSize_condEta2[v])
                    arr_sum[v] = sumen
                    sizen = ak.size(VarToSumSize_condEta2[v])
                    arr_size[v] = sizen
                CellGrid[:, i, j, varIndex] = arr_sum
                varIndex +=1
                CellGrid[:, i, j, varIndex] = arr_size


def getArraySumdRbins(allVars, allVarsIndex, variableToAdd, deltaR, verbose =False):
    dR_bin = [0.,0.1,0.2,0.3,0.4,0.5]
    nEvents = len(variableToAdd)
    for i in range(0,len(dR_bin)-1):
        if(verbose):
            print('dR bin = ', dR_bin[i], ' - ' ,dR_bin[i+1])
        # define arrays of DeltaR and energy
        if(verbose):
            print('defining deltaR and var')
        dR_arr0 = deltaR
        en_arr0 = variableToAdd
        # apply first condition
        if(verbose):
            print('applying first condition')
        condition1 = (dR_arr0)>=dR_bin[i]
        dR_cond1 = dR_arr0[condition1]
        en_cond1 = en_arr0[condition1]
        # apply second condition
        if(verbose):
            print('applying second condition')
        condition2 = (dR_cond1)<dR_bin[i+1]
        en_cond2 = en_cond1[condition2]
        dR_cond2 = dR_cond1[condition2]
        # get the sum and the size of each vector
        if(verbose):
            print('defining sum and size arrays')
        arr_sum = np.zeros(nEvents)
        arr_size = np.zeros(nEvents)
        if(verbose):
            print('filling sum and size arrays')
        for j in range(0, len(en_cond2)-1):
            sumen = ak.sum(en_cond2[j])
            arr_sum[j] = sumen
            sizen = ak.size(en_cond2[j])
            arr_size[j] = sizen
        if(verbose):
            print('filled sum and size arrays, now filling allVars')
        allVars[:,allVarsIndex]= arr_sum
        allVarsIndex=allVarsIndex+1
        allVars[:,allVarsIndex]= arr_size
        allVarsIndex=allVarsIndex+1
        if(verbose):
            print ('current var index = ',allVarsIndex, ' at end of bin ', dR_bin[i], ' - ' ,dR_bin[i+1], '\n')

def GetFeatureMatrix(awkArray, flatVars, vecVars, verbose=False):
    print("writing feature matrix")
    nVars = len(flatVars)+10*len(vecVars)
    nEvents = len(getattr(awkArray, flatVars[0]))
    allVars = np.zeros((nEvents,nVars))
    if(verbose):
        print("nVars = ", nVars, "\nnEvents = ", nEvents)
    j = 0
    for i in flatVars:
        allVars[:,j]=getattr(awkArray,i)
        j=j+1
    if(verbose):
        print('first ', j, ' vars added to allVars: ', flatVars)
    for i in vecVars:
        print('adding', i)
        getArraySumdRbins(allVars, j, getattr(getattr(awkArray, i), vecVars[i]), getattr(getattr(awkArray, i), 'DeltaR'), verbose)
        j=j+10
        if(verbose):
            print(j, ' vars added to allVars: ', vecVars)
    #np.save('/home/valeria/FeatMatrix.npy', allVars)
    return allVars

def GetTrainTestFraction(MCTruth, featMatrix, trainFraction, valFraction, verbose = False):
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

def GetTrainTestWeight(awkArray, trainFraction, valFraction,weightColumnName='weight'):
    weights= getattr(awkArray,weightColumnName)
    weights = np.array(ak.to_list(awkArray.weight ))
    ak.to_numpy(weights)
    nEvents = len(weights)
    testFraction = 1 - (trainFraction+valFraction)
    nTrain = int(round(nEvents * trainFraction, 0))
    nVal = nTrain+int(round(nEvents * valFraction, 0))
    trainWeights=  weights[0:nTrain]
    valWeights= weights[nTrain:nVal]
    testWeights= weights[nVal:nEvents]
    return trainWeights, testWeights, valWeights

def ReshapeVariables(x, y, w, batch_size):
    n_events = len(y)
    n_vars = len(x[0])
    n_batches = int(n_events/batch_size)
    x_reshaped = x.reshape(batch_size,n_batches,n_vars)
    y_reshaped = y.reshape(batch_size,n_batches)
    w_reshaped = w.reshape(batch_size,n_batches)
    return x_reshaped, y_reshaped, w_reshaped




from tensorflow.keras.models import Model
class BasicModel(Model):
    def __init__(self, params):
        super(BasicModel, self).__init__()

        # Dense Block
        self.dense = []
        self.dropout_dense = []
        self.batch_norm_dense = []
        for i in range(params['num_den_layers']):
            self.dense.append(tf.keras.layers.Dense(params['num_units_den_layers'],
                                                        activation=params['activation_dense'],
                                                        name='dense_{}'.format(i)))
            self.batch_norm_dense.append(tf.keras.layers.BatchNormalization(name='batch_normalization{}'.format(i)))
            if params['dropout_rate_den_layers'] > 0:
                self.dropout_dense.append(tf.keras.layers.Dropout(params['dropout_rate_den_layers'],
                                                      name='dropout_dense_{}'.format(i)))


        self.final_dense = tf.keras.layers.Dense(1, activation="sigmoid", name='output')
    @tf.function
    def call(self, x):
        for i in range(len(self.dense)):
            x = self.dense[i](x)
            if len(self.batch_norm_dense) > i:
                x = self.batch_norm_dense[i](x)
            if len(self.dropout_dense) > i:
                x = self.dropout_dense[i](x)

        last_pre = x

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
