#count average Taus
import ROOT
import numpy as np
import pandas as pd
import root_pandas

#tt = root_pandas.read_root(inpath + tt_file, 'events_3j1t')#, where='Jet_score_best>0.7 || Jet_score_secondbest>0.7')
#wjets = root_pandas.read_root(inpath + wj_file, 'events_3j1t')#, where='Jet_score_best>0.7 || Jet_score_secondbest>0.7')
#stbb = root_pandas.read_root(inpath + stbb_file, 'events_3j1t' )#,where='Jet_score_best>0.7 || Jet_score_secondbest>0.7')
#stbq = root_pandas.read_root(inpath+stbq_file, 'events_3j1t', where='Jet_score_best>0.7 || Jet_score_secondbest>0.7')
#absolute_path = "/Users/valeriadamante/Desktop/Dottorato/gridui/L2SkimmedTuples/DataSetTraining/"
absolute_path = "/home/users/damante/L2SkimmedTuples/DataSetTraining/"
#print(ROOT.RDataFrame('L2TauTrainTuple',absolute_path+"/all_Data.root").GetColumnNames())

df_prova = root_pandas.read_root(absolute_path+"all_Data.root", 'L2TauTrainTuple', {'evt','l1Tau_pt', 'l1Tau_hwIso'})
#print(df_prova.head())
#print( type(df_prova))

#print(df_prova['l1Tau_pt'])

df = df_prova[(df_prova.l1Tau_pt>=32) & ((df_prova.l1Tau_hwIso>0) | (df_prova.l1Tau_pt>=70))].groupby('evt').count()
print("df shape")
print(df.shape[0])
print("df count>=2")
print(df[df.l1Tau_pt<2].shape[0])
#print(df['l1Tau_pt'])
print(pd.DataFrame.mean(df['l1Tau_pt']))
