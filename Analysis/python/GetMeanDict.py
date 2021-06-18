import ROOT

#def GetMeanValue(df):
absolute_path = "/Users/valeriadamante/Desktop/Dottorato/gridui/L2SkimmedTuples/DataSetTraining/"
df = ROOT.RDataFrame("L2TauTrainTuple", absolute_path+"miniTuple.root")
#print( "column \t\t mean \t\t std_dev")
print("{")
for i in df.GetColumnNames():
    if(df.StdDev(i).GetValue()!=0):
        print("\n\t\"%s\":{\"mean\":%.5f, \"std\":%.5f, \"lim_min\":-5, \"lim_max\":5}," % (i, df.Mean(i).GetValue(), df.StdDev(i).GetValue()))
print("\n}")
