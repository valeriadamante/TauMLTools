#!/usr/bin/env python
# coding: utf-8

import tensorflow as tf
import uproot
import awkward as ak
import pickle
import numpy as np
import json
from DataLoaderL2 import *
import matplotlib
matplotlib.use('agg')
import pylab as plt
from sklearn.metrics import classification_report, roc_auc_score, accuracy_score, f1_score, confusion_matrix
import argparse
import os
parser = argparse.ArgumentParser()
parser.add_argument('--machine', required=False, type=str, default="local", choices=["local", "cmssimphase2"]) #aggiungi pccms65
parser.add_argument('--n_max_events', required=False, type=int, default=-1, help='max number of events to be processed')
parser.add_argument('--input_file', required=False, type=str, default='DataSetTrainingWeight.root', help='input file name')
parser.add_argument('--input_tuple', required=False, type=str, default='L2TauTrainTuple', help='input tree name')
parser.add_argument('--featMatrix_file', required=False, type=str, default='FeatMatrix.npy', help='output file name')
parser.add_argument('--McTruth_file', required=False, type=str, default='MCTruth.npy', help='output file name')
parser.add_argument('--Weights_file', required=False, type=str, default='weights.npy', help='output file name')
parser.add_argument('--verbose', required=False, type=int, default=0)
args = parser.parse_args()
# ****** Define directories *******
dir_dict = {} # add also pccms65
dir_dict["cmssimphase2"]={}
dir_dict["cmssimphase2"]["data"]="/home/valeria/DataSetTraining/"
dir_dict["cmssimphase2"]["output"]="/home/valeria/output/"
dir_dict["local"]={}
dir_dict["local"]["data"]="/Users/valeriadamante/Desktop/Dottorato/L2SkimmedTuples/DataSetTraining/"
dir_dict["local"]["output"]="/Users/valeriadamante/Desktop/Dottorato/cmssimphase2/outputProva/"
# ******** Define file & tuple names ********
absolute_path =  dir_dict[args.machine]["data"]
treeName = args.input_tuple
inFile = absolute_path+args.input_file
featMatrixFile = absolute_path+args.featMatrix_file
MCTruthFile = absolute_path+args.McTruth_file
WeightsFile = absolute_path+args.Weights_file
outDir = dir_dict[args.machine]["output"]
if not os.path.exists(outDir):
    os.mkdir(outDir)
# ***** Load Feat Matrix ******
if not os.path.exists(featMatrixFile):
    flatVars = ['nVertices','l1Tau_pt', 'l1Tau_eta', 'l1Tau_hwIso']
    vecVars = {
            'caloRecHit_e':'energy',
            'caloRecHit_had':'energy',
            'patatrack':'pt'
    }
    awkArray = GetAwkArray(inFile, treeName)
    featMatrix = GetFeatureMatrix(awkArray, flatVars, vecVars, False)
else:
    print(("file {} exists and taking data from there").format(featMatrixFile))
    featMatrix = np.load(featMatrixFile)
# ***** Load MCTruth Matrix *******
if(args.verbose)>0:
    print("Loading MCTruthFile")
if not os.path.exists(MCTruthFile):
    awkArray = GetAwkArray(inFile, treeName)
    MCTruth = GetMCTruth(awkArray, 'genLepton_isTau')
    np.save(MCTruthFile, MCTruth)
else:
    if(args.verbose)>0:
        print(("file {} exists and taking data from there").format(MCTruthFile))
    MCTruth = np.load(MCTruthFile)
# ****** Load Weights matrix *******
if not os.path.exists(WeightsFile):
    if(args.verbose)>0:
        print("Loading WeightsFile")
    awkArray = GetAwkArray(inFile, treeName)
    weights = GetMCTruth(awkArray, 'weight')
    np.save(WeightsFile, weights)
else:
    if(args.verbose)>0:
        print(("file {} exists and taking data from there").format(WeightsFile))
    weights = np.load(WeightsFile)

# create all vars array, mctruth array and feature matrix

params = {
    'activation_dense':'relu',
    'num_units_den_layers':int(2*len(featMatrix[0])/0.8),
    #'dropout_rate_den_layers':0.2,
    'batch_size':200,
    'train_fraction': 0.6102414433536747,
    'validation_fraction': 0.19474661713982488,
    #'learning_rate':0.002,
    'monitor_es':'val_loss',
    #'patience_es':2,
    #'epochs':10,
    'patience_es':10,
    'epochs':100000,
}
sepstr = '*'*80
# ****** Get train - test - validation samples ********
if(args.verbose)>0:
    print("preparing train/test/val samples")
number_of_batches = len(MCTruth)/params['batch_size']
x_train, y_train, w_train, x_test, y_test, w_test, x_val, y_val, w_val = GetTrainTestFraction(MCTruth, featMatrix, weights, params['train_fraction'],  params['validation_fraction'], args.verbose) 
for ly in range(5,11):
    for lr in np.arange(0.001, 0.02, 0.002):
        for do in np.arange(0.2, 0.5, 0.1):
            #ly = 5
            #lr = 0.001
            #do = 0.2
            params['dropout_rate_den_layers']= do
            params['learning_rate'] = lr
            params['num_den_layers'] = ly
            # creating specific model files_directoriesmodel_directory="output/model_"+file_suffix


            file_suffix ="{}Layers_{:.2f}Dropout_{:.3f}LearningRate".format(params['num_den_layers'],params['dropout_rate_den_layers'],params['learning_rate']).replace(".","p")
            print(file_suffix)

            model_directory=outDir+"/model_"+file_suffix
            if not os.path.exists(model_directory):
                os.mkdir(model_directory)

            log_directory = model_directory+"/log"
            if not os.path.exists(log_directory):
                os.mkdir(log_directory)

            history_directory = model_directory+"/history"
            if not os.path.exists(history_directory):
                os.mkdir(history_directory)

            plot_directory = outDir+"/plots"
            if not os.path.exists(plot_directory):
                os.mkdir(plot_directory)

            var_directory = model_directory+"/vars"
            if not os.path.exists(var_directory):
                os.mkdir(var_directory)


            # creating files to be written

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

            print("defining model ")

            # define model - loss - optimizer - compile model

            mymodel_path = model_directory+"/model_"+file_suffix
            if os.path.exists(mymodel_path):
                print(("{} already exists").format(mymodel_path))
                continue
            model = BasicModel(params)
            predictions = model.call(x_train).numpy()

            loss_fn = tf.keras.losses.BinaryCrossentropy(from_logits=False)
            #loss_fn(y_train, predictions).numpy()

            opt = tf.keras.optimizers.Adam(
            learning_rate=params['learning_rate'],
            )

            model.compile(optimizer=opt,
                      loss=loss_fn,
                      metrics=['accuracy'])

            model.build(x_train.shape)
            model.summary()
            logfile.write(sepstr+"\n")
            logfile.write("\t \t \t model summary \n")
            logfile.write(sepstr+"\n")
            model.summary(print_fn=lambda x: logfile.write(x + '\n'))
            logfile.write("\n"+sepstr+"\n")
            logfile.write("\t \t end of model summary \n")
            logfile.write(sepstr+"\n")

            print("preparing training")
            csvfile=model_directory+"/model_"+file_suffix+".csv"
            model_path= model_directory+"/model_"+file_suffix+".h5"
            #mymodel_path= model_directory+"/model_"+file_suffix
            # prepare training
            early_stopping = tf.keras.callbacks.EarlyStopping(monitor=params['monitor_es'], patience=params['patience_es'])
            csv_logger = tf.keras.callbacks.CSVLogger(csvfile, separator=',', append=False )
            #model_checkpt =  tf.keras.callbacks.ModelCheckpoint(filepath=model_path, monitor=params['monitor_es'], verbose=1, save_best_only=True, save_weights_only=False, mode='auto')
            my_own_cb = SaveBestModel(output=mymodel_path, verbose=2)
            history = model.fit(x_train, y_train, batch_size = params['batch_size'], verbose=2 ,
                            epochs = params['epochs'] , callbacks=[my_own_cb, csv_logger,early_stopping,], sample_weight=w_train,
                validation_data = (x_val, y_val, w_val))

            np.save(history_file_path+'.npy',history.history)
            f2 = open(history_file_path+".pkl","wb")
            pickle.dump(history.history,f2)
            f2.close()
            # make predictions for test sample
            '''
            y_pred = model.predict(x_test)
            y_pred.resize(len(y_test))
            predictions = [round(value) for value in y_pred]
            # evaluate predictions
            accuracy = accuracy_score(y_test, predictions)
            print("Accuracy: %.2f%%" % (accuracy * 100.0))

            logfile.write("\n"+sepstr+"\n")
            logfile.write("\t \t end of accuracy \n")
            logfile.write(sepstr+"\n")
            logfile.write(str(accuracy))
            logfile.write("\n"+sepstr+"\n")
            logfile.write("\t \t end of accuracy \n")
            logfile.write(sepstr+"\n")

            logfile.close()
            #get_ipython().run_line_magic('matplotlib', 'inline')


            # plot loss vs epoch for validation and training
            fig = plt.figure(figsize=(20,20))
            #plt.figure()
            ax = fig.add_subplot(2, 2, 1)
            ax.plot(history.history['loss'], label='loss')
            ax.plot(history.history['val_loss'], label='val_loss')
            ax.legend(loc="upper right")
            ax.set_title('loss VS epoch')
            ax.set_xlabel('epoch')
            ax.set_ylabel('loss')

            # plot accuracy vs epoch for validation and training
            ax = fig.add_subplot(2, 2, 2)
            ax.plot(history.history['accuracy'], label='acc')
            ax.plot(history.history['val_accuracy'], label='val_acc')
            ax.legend(loc="lower left")
            ax.set_title('accuracy VS epoch')
            ax.set_xlabel('epoch')
            ax.set_ylabel('acc')

            # Plot ROC
            # generate a no skill prediction (majority class)
            ns_probs = [0 for _ in range(len(y_test))]
            y_predict = model.predict(x_test)
            print(len(y_predict))
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
            np.save(var_directory+'/fpr_'+file_suffix+'.npy', fpr) # save
            np.save(var_directory+'/tpr_'+file_suffix+'.npy', tpr) # save
            np.save(var_directory+'/thresholds_'+file_suffix+'.npy', thresholds) # save
            np.save(var_directory+'/y_predict_'+file_suffix+'.npy', y_predict) # save
            np.save(var_directory+'/y_test_'+file_suffix+'.npy', y_test) # save
            roc_auc = auc(fpr, tpr)
            ax = fig.add_subplot(2, 2, 3)
            ax.plot(fpr, tpr, lw=2, ms=2.1, color='cyan', label='auc = %.3f' % (roc_auc))
            ax.plot([0, 1], [0, 1], linestyle='--', lw=2, ms=2.1, color='k', label='random chance')
            ax.set_xlim([0, 1.0])
            ax.set_ylim([0, 1.0])
            ax.set_xlabel('false positive rate')
            ax.set_ylabel('true positive rate')
            ax.set_title('receiver operating curve')
            ax.legend(loc="lower right")

            # plot distribution signal and background for train and test
            decisions, low_high=compare_train_test(model, x_train, y_train, x_test, y_test)

            ax = fig.add_subplot(2,2,4)
            bins=30
            ax.hist(decisions[0], color='b', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True,
                label='S (train)')
            ax.hist(decisions[1], color='r', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True,
                label='B (train)')
            hist, bins = np.histogram(decisions[2], bins=bins, range=low_high, density=True)
            scale = len(decisions[2]) / sum(hist)
            err = np.sqrt(hist * scale) / scale
            width = (bins[1] - bins[0])
            center = (bins[:-1] + bins[1:]) / 2
            ax.errorbar(center, hist, yerr=err, fmt='o', c='b', label='S (test)')

            hist, bins = np.histogram(decisions[3], bins=bins, range=low_high, density=True)
            scale = len(decisions[3]) / sum(hist)
            err = np.sqrt(hist * scale) / scale

            ax.errorbar(center, hist, yerr=err, fmt='o', c='r', label='B (test)')
            ax.set_title('S,B distribution')
            ax.set_xlabel("Output")
            ax.set_ylabel("Arbitrary units")
            ax.legend(loc='best')
            #plt.show()
            for par in params:
                print(par , "= ", params[par])


            fig.savefig(plot_directory+'/plots_'+file_suffix+'.pdf')
            fig.savefig(plot_directory+'/plots_'+file_suffix+'.png')

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
            fig.savefig(plot_directory+'/ROC_val_'+file_suffix+'.pdf')
            fig.savefig(plot_directory+'/ROC_val_'+file_suffix+'.png')
            #model.save("output/model_"+file_suffix)\
            '''
