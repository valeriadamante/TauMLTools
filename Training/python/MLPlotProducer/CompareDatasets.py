
logfile_norm = "/Users/valeriadamante/Desktop/Dottorato/gridui/CMSSW_11_2_1_Patatrack/src/cellGridStructure.log"
logfile_CMSSW = "/Users/valeriadamante/Desktop/Dottorato/public/CMSSW_11_2_1_Patatrack/src/cellGridStructureCMSSW.log"
with open(logfile_norm) as f1, open(logfile_CMSSW) as f2:
    k=0
    l=0
    for x, y in zip(f1, f2):
        k+=1
        x = float(x.replace(',', ''))
        y = float(y.replace(',', ''))
        if(round(y,2)!=round(x,2)):
            print("sono diversi {}\t{} alla linea {}".format(x, y,k))
            l+=1

print(("ci sono {} diversi").format(l))
