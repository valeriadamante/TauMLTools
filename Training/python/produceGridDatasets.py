# produce datasets so that I never have to reload them
from TauMLTools.Training.python.DataLoaderL2ForCNN import *
#from DataLoaderL2ForCNN import *
import numpy as np
# ****** Define directories *******
dir_dict = {} # add also pccms65
dir_dict["cmssimphase2"]={}
dir_dict["cmssimphase2"]["data"]="/data/tau-l2/DataSetTraining/"
dir_dict["cmssimphase2"]["model"]="/home/valeria/model/"
dir_dict["cmssimphase2"]["output"]="/home/valeria/output/"
dir_dict["local"]={}
dir_dict["local"]["data"]="/Users/valeriadamante/Desktop/Dottorato/L2SkimmedTuples/DataSetTraining/"
dir_dict["local"]["model"]="/Users/valeriadamante/Desktop/Dottorato/cmssimphase2/model/"
dir_dict["local"]["output"]="/Users/valeriadamante/Desktop/Dottorato/cmssimphase2/output/"
effRate_dict ={"eff":"VBF", "rate":"data", "test":"DataSetTraining","Zprime":"Zprime"}
def GetDataSetPath(machine):
    return dir_dict[machine]["data"]
def GetOutPath(machine):
    return dir_dict[machine]["output"]
def GetCellGridNormMatrix(machine, effRate, **kwargs):
    n_max_events= kwargs['n_max_events']
    n_cellsX=kwargs['n_cellsX']
    n_cellsY=kwargs['n_cellsY']
    timeInfo=kwargs['timeInfo']
    verbose=kwargs['verbose']
    # ******** Define file & tuple names ********
    absolute_path =  dir_dict[machine]["data"]
    outDir = dir_dict[machine]["output"]
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    plotDir = outDir+"/plots"
    if not os.path.exists(plotDir):
        os.mkdir(plotDir)
    plotDir = plotDir+("/StdDistrib_{}").format(effRate_dict[effRate])
    if not os.path.exists(plotDir):
        os.mkdir(plotDir)
    plotDir+="/"
    treeName = 'L2TauTrainTuple'
    inFile_name = 'DataSetTrainingWeight.root'
    outFile_name =  'CellGrid.npy'
    outFileNorm_name='CellGridNorm.npy'
    dictName = 'NormalizationDict.json'
    if(n_max_events>0):
        outFile_name =  ('CellGrid_{}Evts.npy').format(n_max_events)
        outFileNorm_name=('CellGridNorm_{}Evts.npy').format(n_max_events)
    dictFileToRead = absolute_path+dictName
    isTrainingDataSet = True
    if(effRate=='eff'):
        inFile_name = 'miniTuple_VBF.root'
        outFile_name =  'CellGridVBF.npy'
        outFileNorm_name='CellGridVBFNorm.npy'
        dictName = 'NormalizationDictVBF.json'
        if(n_max_events>0):
            outFile_name =  ('CellGridVBF_{}Evts.npy').format(n_max_events)
            outFileNorm_name=('CellGridNormVBF_{}Evts.npy').format(n_max_events)
        isTrainingDataSet = False
    if(effRate=='Zprime'):
        inFile_name = 'miniTuple_Zprime.root'
        outFile_name =  'CellGridZprime.npy'
        outFileNorm_name='CellGridZprimeNorm.npy'
        dictName = 'NormalizationDictZprime.json'
        if(n_max_events>0):
            outFile_name =  ('CellGridZprime_{}Evts.npy').format(n_max_events)
            outFileNorm_name=('CellGridNormZprime_{}Evts.npy').format(n_max_events)
        isTrainingDataSet = False
    elif(effRate=='rate'):
        inFile_name = 'miniTuple_Data.root'
        outFile_name =  'CellGridData.npy'
        outFileNorm_name='CellGridDataNorm.npy'
        dictName = 'NormalizationDictData.json'
        if(n_max_events>0):
            outFile_name =  ('CellGridData_{}Evts.npy').format(n_max_events)
            outFileNorm_name=('CellGridNormData_{}Evts.npy').format(n_max_events)
        isTrainingDataSet = False

    inFile = absolute_path+inFile_name
    outFile = absolute_path+outFile_name
    outFileNorm = absolute_path+outFileNorm_name
    dictFileToSave = absolute_path+dictName

    return getNormCellGridMatrix(n_cellsX, n_cellsY, inFile, treeName, outFile, outFileNorm, dictFileToSave, dictFileToRead, plotDir, timeInfo, n_max_events, verbose, isTrainingDataSet, plot=True)

def GetMCTruthWeights(machine, McTruth_file,Weights_file, **kwargs):
    n_max_events= kwargs['n_max_events']
    n_cellsX=kwargs['n_cellsX']
    n_cellsY=kwargs['n_cellsY']
    verbose=kwargs['verbose']
    # ******** Define file & tuple names ********
    absolute_path =  dir_dict[machine]["data"]
    treeName = 'L2TauTrainTuple'
    inFile_name = 'DataSetTrainingWeight.root'
    inFile = absolute_path+inFile_name
    MCTruthFile = absolute_path+McTruth_file
    WeightsFile = absolute_path+Weights_file
    # ***** Load MCTruth Matrix *******
    if(verbose)>0:
        print("Loading MCTruthFile")
    if not os.path.exists(MCTruthFile):
        awkArray = GetAwkArray(inFile, treeName)
        MCTruth = GetMCTruth(awkArray, 'genLepton_isTau')
        np.save(MCTruthFile, MCTruth)
    else:
        if(verbose)>0:
            print(("file {} exists and taking data from there").format(MCTruthFile))
        MCTruth = np.load(MCTruthFile)
    # ****** Load Weights matrix *******
    if not os.path.exists(WeightsFile):
        if(verbose)>0:
            print("Loading WeightsFile")
        awkArray = GetAwkArray(inFile, treeName)
        weights = GetMCTruth(awkArray, 'weight')
        np.save(WeightsFile, weights)
    else:
        if(verbose)>0:
            print(("file {} exists and taking data from there").format(WeightsFile))
        weights = np.load(WeightsFile)

    return MCTruth, weights

def GetModelDir(machine):
    return dir_dict[machine]["model"]

def GetModelPath(machine, params):
    file_suffix=("model_{}D{}CNN1{}CNN2_{:.2f}Dropout_{:.4f}LearningRate").format(params['num_dense_layers'],params['num_CNN1x1_layers'],params['num_CNN2x2_layers'],params['dropout_dense_layers'],params['learning_rate'])
    print(file_suffix)
    file_suffix = file_suffix.replace(".","p")
    model_directory=dir_dict[machine]["model"]+"/"+file_suffix
    if not os.path.exists(model_directory):
        print(("model directory {} does not exist!").format(model_directory))
        return ''
    model_path = model_directory+"/"+file_suffix
    if os.path.exists(model_path):
        print(("{} exists").format(model_path))
    else:
        print(("{} does not exist").format(model_path))
    return model_path

def GetPlotDir(machine):
    plot_directory = dir_dict[machine]["output"]+"/plots"
    if not os.path.exists(plot_directory):
        print(("{} does not exist, creating it").format(plot_directory))
        os.mkdir(plot_directory)
    return plot_directory

def GetRootPath(machine, effRate):
    absolute_path =  dir_dict[machine]["data"]
    inFile_name = 'DataSetTrainingWeight.root'
    if(effRate=='eff'):
        inFile_name = 'miniTuple_VBF.root'
    elif(effRate=='rate'):
        inFile_name = 'miniTuple_Data.root'
    elif(effRate=='Zprime'):
        inFile_name = 'miniTuple_Zprime.root'
    inFile = absolute_path+inFile_name
    return inFile

def GetDictFile(machine, effRate, **kwargs):
    n_max_events= kwargs['n_max_events']
    n_cellsX=kwargs['n_cellsX']
    n_cellsY=kwargs['n_cellsY']
    timeInfo=kwargs['timeInfo']
    verbose=kwargs['verbose']
    # ******** Define file & tuple names ********
    absolute_path =  dir_dict[machine]["data"]
    treeName = 'L2TauTrainTuple'
    inFile_name = 'DataSetTrainingWeight.root'
    outFile_name =  'CellGrid.npy'
    outFileNorm_name='CellGridNorm.npy'
    dictName = 'NormalizationDict.json'
    if(n_max_events>0):
        outFile_name =  ('CellGrid_{}Evts.npy').format(n_max_events)
        outFileNorm_name=('CellGridNorm_{}Evts.npy').format(n_max_events)
    dictFileToRead = absolute_path+dictName
    isTrainingDataSet = True
    if(effRate=='eff'):
        inFile_name = 'miniTuple_VBF.root'
        outFile_name =  'CellGridVBF.npy'
        outFileNorm_name='CellGridVBFNorm.npy'
        dictName = 'NormalizationDictVBF.json'
        if(n_max_events>0):
            outFile_name =  ('CellGridVBF_{}Evts.npy').format(n_max_events)
            outFileNorm_name=('CellGridNormVBF_{}Evts.npy').format(n_max_events)
        isTrainingDataSet = False
    elif(effRate=='rate'):
        inFile_name = 'miniTuple_Data.root'
        outFile_name =  'CellGridData.npy'
        outFileNorm_name='CellGridDataNorm.npy'
        dictName = 'NormalizationDictData.json'
        if(n_max_events>0):
            outFile_name =  ('CellGridData_{}Evts.npy').format(n_max_events)
            outFileNorm_name=('CellGridNormData_{}Evts.npy').format(n_max_events)
        isTrainingDataSet = False

    inFile = absolute_path+inFile_name
    outFile = absolute_path+outFile_name
    outFileNorm = absolute_path+outFileNorm_name
    dictFileToSave = absolute_path+dictName
    return inFile, outFile, outFileNorm, dictFileToRead, dictFileToSave

def GetNameForEffRate(effRate):
    return effRate_dict[effRate]

def GetFileSuffix(params):
    file_suffix=("model_{}D{}CNN1{}CNN2_{:.2f}Dropout_{:.4f}LearningRate").format(params['num_dense_layers'],params['num_CNN1x1_layers'],params['num_CNN2x2_layers'],params['dropout_dense_layers'],params['learning_rate'])
    file_suffix = file_suffix.replace(".","p")
    return file_suffix
