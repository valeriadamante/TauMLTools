# algo efficiency for test dataset to compare with sqrt(diagonal) of event cutbased efficiency
# how to build this algorithm:

# 1. goal = compare algo efficiency per taus (NOT per events) on testing dataset with sqrt(diagonal) of cut based (algo?) efficiency (on VBF) per event (NOT per taus) in pt bins

from TauMLTools.Training.python.DataLoaderL2 import *
import tensorflow as tf
import awkward as ak
import numpy as np
from sklearn.metrics import classification_report, roc_auc_score,roc_curve, auc, accuracy_score, f1_score, confusion_matrix
import os
import argparse
import ROOT

parser = argparse.ArgumentParser(description='')
parser.add_argument('--machine', required=True, type=str, choices =['local', 'gridui'], default = 'local', help="input directory")
parser.add_argument('--verbose', required=False, type=int, default = 0)
parser.add_argument('--dropout', required=False, type=float, default = 0.2)
parser.add_argument('--num_den_layers', required=False, type=int, default = 5)
parser.add_argument('--learning_rate', required=False, type=float, default =0.002)
#parser.add_argument('--n-threads', required=False, type=int, default=1, help="number of threads")
args = parser.parse_args()


tupleName = 'L2TauTrainTuple'
trainDir_path = '/TauMLTools/Training/python/output/'
datasetDir_path= "/L2SkimmedTuples/DataSetTraining/"
flatVars = ['nVertices','l1Tau_pt', 'l1Tau_eta', 'l1Tau_hwIso']
vecVars = {
        'caloRecHit_e':'energy',
        'caloRecHit_had':'energy',
        'patatrack':'pt'
}
dnn_path=""

if(args.machine =='local'):
    absolute_path = "/Users/valeriadamante/Desktop/Dottorato/"
    dnn_path = absolute_path+'/cmssimphase2/'+trainDir_path
elif (args.machine == 'gridui'):
    absolute_path = "/home/users/damante/"
    dnn_path = absolute_path+'/CMSSW_11_2_1_Patatrack/src/'+trainDir_path

dataset_path = absolute_path+datasetDir_path
print(dnn_path)

variables = {'evt','l1Tau_pt', 'l1Tau_hwIso'}

fileName_data = 'miniTuple_Data.root'
outfileName_data = 'miniEvtTuple_Data.root'
datasetFile_name_data = 'datasetData.npy'

fileName_test = 'DataSetTrainingWeight.root'
outfileName_test = 'miniEvtTuple.root'
datasetFile_name_test ='dataset.npy'

#check if npy file with dataset exists
# load dataset for data
if not os.path.exists(("{}/{}").format(dataset_path, datasetFile_name_data)):
    awkArray_data = GetAwkArrayData(dataset_path+fileName_data, tupleName)
    featMatrix_data = GetFeatureMatrix(awkArray_data, flatVars, vecVars, False)
    np.save(("{}/{}").format(dataset_path, datasetFile_name_data), featMatrix_data)
else:
    if(args.verbose>0):
        print(("file {}/{} exists and taking data from there").format(dataset_path, datasetFile_name_data))
    featMatrix_data = np.load(('{}/{}').format(dataset_path, datasetFile_name_data))

# load dataset for test
# dataset - x
if not os.path.exists(("{}/{}").format(dataset_path, datasetFile_name_test)):
    awkArray_test = GetAwkArrayData(dataset_path+fileName_test, tupleName)
    featMatrix_test = GetFeatureMatrix(awkArray_test, flatVars, vecVars, False)
    np.save(("{}/{}").format(dataset_path, datasetFile_name_test), featMatrix_test)
else:
    if(args.verbose>0):
        print(("file {}/{} exists and taking data from there").format(dataset_path, datasetFile_name_test))
    featMatrix_test = np.load(('{}/{}').format(dataset_path, datasetFile_name_test))

# dataset - MCtruth
if not os.path.exists(('{}/MCTruth.npy').format(dataset_path)):
    awkArray_test = GetAwkArray(dataset_path+fileName_test, tupleName)
    MCTruth = GetMCTruth(awkArray_test, 'genLepton_isTau')
else:
    if(args.verbose>0):
        print(("file {}/MCTruth.npy exists and taking data from there").format(dataset_path))
    MCTruth = np.load(('{}/MCTruth.npy').format(dataset_path))


params = {
    'activation_dense':'relu',
    'num_units_den_layers':int(2*len(featMatrix_test[0])/0.8),
    'l1Pt_position':1,
    'batch_size':200,
    'train_fraction': 0.6102414433536747,
    'validation_fraction': 0.19474661713982488,
    'opt_threshold' : 0.30631134733375802,
    'bigOrRate': 13603.37
}

params['dropout_rate_den_layers']= args.dropout
params['learning_rate'] = args.learning_rate
params['num_den_layers'] = args.num_den_layers
file_suffix = "model_{}Layers_{:.2f}Dropout_{:.3f}LearningRate".format(params['num_den_layers'],params['dropout_rate_den_layers'],params['learning_rate']).replace(".","p")
# dataset - weight
if not os.path.exists(("{}/weights_test.npy").format(dataset_path)):
    awkArray_test = GetAwkArrayData(dataset_path+fileName_test, tupleName)
    w_train, w_test, w_val = GetTrainTestWeight(awkArray_test, params['train_fraction'],params['validation_fraction'])
    np.save(("{}/weights_test.npy").format(dataset_path), w_test)
else:
    if(args.verbose>0):
        print(("file {}/weights_test.npy exists and taking data from there").format(dataset_path))
    w_test = np.load(("{}/weights_test.npy").format(dataset_path))

x_train, y_train, x_test, y_test, x_val, y_val = GetTrainTestFraction(MCTruth, featMatrix_test, params['train_fraction'],
                                                                      params['validation_fraction'], False)

model_path = dnn_path+ file_suffix +"/"+file_suffix

pt_bins = [20,30,35,40,50,60,70,80,90,350]
numerator = np.zeros(len(pt_bins))
print(len(pt_bins))
denominator = np.zeros(len(pt_bins))
algo_eff = np.zeros(len(pt_bins))
th = np.array(len(pt_bins))
eff = np.array(len(pt_bins))



for index in range(0, len(pt_bins)-1):
    print(("bin is [{},{}]").format(pt_bins[index], pt_bins[index+1]))
    # 2. PART-1 : get the algo efficiency per taus on testing dataset binned in Pt
    # 2.1 obtain dataset binned in pt bins : x_test_pt_20_30 ecc

    pt_low_threshold = pt_bins[index]
    pt_high_threshold = pt_bins[index+1]
    if(args.verbose>1):
        print("obtained all vars , now computing bin subsamples ")
    x_test_bin, y_test_bin, w_test_bin = GetTestSubSample(x_test,y_test, w_test, pt_low_threshold, pt_high_threshold, params['l1Pt_position'])
    #print(("{},{}").format(len(x_test_bin),len(x_test_bin[:])))

    # 2.2 make predictions for those datasets

    if(args.verbose>1):
        print(("loading model from {}").format(model_path))

    if not os.path.exists(model_path):
        print(("{} does not exist").format(model_path))

    model = keras.models.load_model(model_path)
    if(args.verbose>1):
        print(("model successfully loaded from {}").format(model_path))

    y_predict_test_bin = model.predict(x_test_bin)
    y_predict_test_bin.resize(len(y_test_bin))
    # 2.3 evaluate algo efficiency

    #numerator = 0
    #denominator = 0
    #th = 0.
    '''
    fpr, tpr, thresholds = roc_curve(y_test_bin, y_predict_test_bin,sample_weight=w_test_bin)
    for i in range(0, len(thresholds)):
        tpr_check = round(tpr[i], 2)
        fpr_check = round(fpr[i], 2)

        if(tpr_check>=0.95 ):
            print (fpr_check)
            if(fpr_check<0.5):
                print(("fpr = {}").format(fpr[i]))
                print(("tpr = {}, threshold = {}").format(tpr_check, thresholds[i]))
                th = (round(thresholds[i],2))
                break
    '''
    opt_threshold = 0.2141452357154776
    for i in range(0, len(y_predict_test_bin)):
        if y_test_bin[i] == 1:
            denominator[index] +=1
            if y_predict_test_bin[i]>opt_threshold:
                numerator[index] +=1
                #numerator[index] *=

    algo_eff[index] = numerator[index] / denominator[index]
    print(("algo efficiency is {}").format(algo_eff[index]))



# 3. PART-2 : get (algo) efficiency per event from the previously saved histogram


ROOT.gStyle.SetPaintTextFormat(".2f")
if(args.verbose>0):
    print("opening canvas, histogram")
#uniformTau = ROOT.TGraph(len(.Pt1Bins)-1,array('d',self.Pt1Bins),array('d',efficiencies))


pt_bins_array = np.array( pt_bins).astype(np.float)
full_name = "/Users/valeriadamante/Desktop/Dottorato/cmssimphase2/outputs/efficiency_algo_hltDoubleL2IsoTau26eta2p2.root"
readfile = ROOT.TFile( full_name, 'READ' )
print(type(readfile.Get("Efficiency")))
effHisto_count = ROOT.TEfficiency(readfile.Get("total_algo_clone"))
total_algohist = ROOT.TH1D(readfile.Get("total_algo"))
print(("efficiency histo count has {} bins").format(total_algohist.GetNbinsX()))
print(np.array( pt_bins))
hPassed = ROOT.TH1D("passed_algo","passed_algo",len(pt_bins)-1,  pt_bins_array)
hTotal = ROOT.TH1D("total_algo","total_algo",len(pt_bins)-1,pt_bins_array)
x = []
y = []
for index in range(0, len(pt_bins)-1):
    x.append(pt_bins[index])
    y.append(algo_eff[index])
    if(args.verbose>0):
         print(("x = [{},{}] \t y = {} ").format(pt_bins[index], pt_bins[index+1],y[index]))
    hPassed.SetBinContent((index+1), numerator[index])
    hTotal.SetBinContent((index+1), denominator[index])
n = len(x)

EffGraph = ROOT.TEfficiency(hPassed,hTotal)
for index in range(0, len(pt_bins)-1):
    print(("histcount efficiency = {}").format(effHisto_count.GetEfficiency(index)))
    print(("dNN efficiency = {}").format(EffGraph.GetEfficiency(index)))
canvas=ROOT.TCanvas()
canvas.cd()
EffGraph.SetMarkerColor( ROOT.kBlue )
EffGraph.SetMarkerStyle( 22 )
EffGraph.SetTitle("EfficiencyDNN;#tau P_{T_} (GeV);#epsilon")




#canvas2=ROOT.TCanvas()
#canvas2.cd()
effHisto_count.SetMarkerColor( ROOT.kRed )
effHisto_count.SetMarkerStyle( 21 )
effHisto_count.SetTitle("Efficiency;#tau P_{T_} (GeV);#epsilon")
effHisto_count.Draw("AP")
EffGraph.Draw("P SAME")
canvas.Update()
canvas.Print("plotEfficiencyDNN.png", "png")
#canvas3.SetLogx()
#canvas2.Update()
#canvas2.Print("plotEfficiencyCOUNT.png", "png")


readfile.Close()

full_name_DNN = "efficiency_DNN.root"
writefile = ROOT.TFile( full_name_DNN, 'RECREATE' )
EffGraph.Write()
writefile.Close()
#raw_input()


# 4. PART-3 : put both efficiencies on the same histogram
