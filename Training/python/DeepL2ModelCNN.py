#!/usr/bin/env python
# coding: utf-8
import tensorflow as tf
import awkward as ak
import pickle
import numpy as np
import json
from produceGridDatasets import *
import argparse
import os
parser = argparse.ArgumentParser()
parser.add_argument('--machine', required=False, type=str, default="local", choices=["local", "cmssimphase2"]) #aggiungi pccms65
parser.add_argument('--effRate', required= False, type=str, default = 'test', choices=['test', 'eff','rate'])
parser.add_argument('--n_max_events', required=False, type=int, default=-1, help='max number of events to be processed')

parser.add_argument('--McTruth_file', required=False, type=str, default='MCTruth.npy', help='output file name')
parser.add_argument('--Weights_file', required=False, type=str, default='weights.npy', help='output file name')
parser.add_argument('--n_cellsX', required=False, type=int, default=5, help='number of cells along X dir')
parser.add_argument('--n_cellsY', required=False, type=int, default=5, help='number of cells along Y dir')
parser.add_argument('--verbose', required=False, type=int, default=0)
parser.add_argument('--gpu', required=False, type=bool, default=False)
args = parser.parse_args()

# ***** Get cell grid Matrix *****
kwArgs = {'n_max_events':args.n_max_events, 'n_cellsX':args.n_cellsX, 'n_cellsY':args.n_cellsY, 'timeInfo' : True, 'verbose' : args.verbose}
CellGridMatrix = GetCellGridNormMatrix(args.machine, args.effRate, **kwArgs)
# ***** Get MC and weights Matrices *****
MCTruth, weights = GetMCTruthWeights(args.machine, args.McTruth_file, args.Weights_file, **kwArgs)
# ****** Define CNN params ******
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
}
sepstr = '*'*80
# ****** Get train - test - validation samples ********
if(args.verbose)>0:
    print("preparing train/test/val samples")
number_of_batches = len(MCTruth)/params['batch_size']
x_train, y_train, w_train, x_test, y_test, w_test, x_val, y_val, w_val = GetTrainTestFraction(MCTruth, CellGridMatrix, weights, params['train_fraction'],  params['validation_fraction'], args.verbose)

# ******* parameter setting ******
# do =
# lr =
# nlayersdense =
# nlayerscnn2=
# nlayerscnn1=
# params['dropout_rate_den_layers']= do
# params['learning_rate'] = lr
# params['num_den_layers'] = ly

# ***** create model directory *****

model_path = GetModelPath(args.machine, params)
model_directory = GetModelDir(args.machine)

log_directory = model_directory+"/log"
if not os.path.exists(log_directory):
    os.mkdir(log_directory)

history_directory = model_directory+"/history"
if not os.path.exists(history_directory):
    os.mkdir(history_directory)


var_directory = model_directory+"/vars"
if not os.path.exists(var_directory):
    os.mkdir(var_directory)
# creating files to be written
file_suffix = GetFileSuffix(params)
logfile = open(log_directory+"/output_"+file_suffix+".log","w+")
history_file_path = history_directory+"/history_"+file_suffix


logfile.write(sepstr+"\n")
logfile.write("\t \t \t parameters \n")
logfile.write(sepstr+"\n")
for par in params:
    str1 = str(par)+" = "+str(params[par])
logfile.write(str1+"\n")
logfile.write(sepstr+"\n")
logfile.write("\t \t end of parameters \n")
logfile.write(sepstr+"\n")

if args.verbose>0:
    print("defining model ")

#  - loss - optimizer - compile model

# ******* define model *******
if os.path.exists(model_path):
    print(("{} already exists").format(model_path))
    #continue
model = CNNModel(params)
#predictions = model.call(x_train).numpy()
loss_fn = tf.keras.losses.BinaryCrossentropy(from_logits=False)
#loss_fn(y_train, predictions).numpy()
opt = tf.keras.optimizers.Adam(learning_rate=params['learning_rate'])
#model.build([ None ] + list(x_train.shape[1:]))
model(x_train[0:2,:,:,:])
model.compile(optimizer=opt,loss=loss_fn,metrics=['accuracy'])
#print(model(x_val[0:2, :, :, :]))
#model.summary()
logfile.write(sepstr+"\n")
logfile.write("\t \t \t model summary \n")
logfile.write(sepstr+"\n")
model.summary(print_fn=lambda x: logfile.write(x + '\n'))
logfile.write("\n"+sepstr+"\n")
logfile.write("\t \t end of model summary \n")
logfile.write(sepstr+"\n")

# ******* Train the model *******
if(args.verbose>0):
    print("preparing training")
csvfile=model_directory+"/"+file_suffix+".csv"
#model_path= model_directory+"/model_"+file_suffix+".h5"
model_path= model_directory+"/"+file_suffix
# prepare training
early_stopping = tf.keras.callbacks.EarlyStopping(monitor=params['monitor_es'], patience=params['patience_es'])
csv_logger = tf.keras.callbacks.CSVLogger(csvfile, separator=',', append=False )
save_model = SaveBestModel(output=model_path, verbose=2)
history = model.fit(x_train, y_train, batch_size = params['batch_size'], verbose=2, epochs = params['epochs'] , callbacks=[save_model, csv_logger,early_stopping], sample_weight=w_train, validation_data = (x_val, y_val, w_val))
# ***** Save everything ******
np.save(history_file_path+'.npy',history.history)
f2 = open(history_file_path+".pkl","wb")
pickle.dump(history.history,f2)
f2.close()
# make predictions for test sample
