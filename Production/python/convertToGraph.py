from tools import *
import os
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
absolute_path = '/afs/cern.ch/work/v/vdamante/public/CMSSW_11_2_1_Patatrack/src/'
path_of_model = absolute_path+'model/model_3D4CNN14CNN2_0p00Dropout_0p0010LearningRate/model_3D4CNN14CNN2_0p00Dropout_0p0010LearningRate/saved_model.pb'
final_path = absolute_path+'graph_model'
model = GetModelPath(args.machine, params)
save_graph()
