import os
command =""
for i in range(1, 91):
    command="root -b /Users/valeriadamante/Desktop/Dottorato/gridui/CMSSW_11_2_1_Patatrack/src/TauMLTools/Analysis/bin/GetHistogramsAndNumbers.cpp+O\\(1,"+str(i)+"\\) & \n wait"
    os.system(command)
