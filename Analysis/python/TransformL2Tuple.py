import os
import argparse


parser = argparse.ArgumentParser()
#parser.add_argument('machine', type=str, default="local", choices=["lxplus", "local"])
parser.add_argument('--reference', required=True, type=str, help= "input reference file")
parser.add_argument('--targets', required=True, type=str, help= "input target files to compare separated by comma")
parser.add_argument('--variables', required=False, type=str, default= "", help= "input variables separated by comma")
parser.add_argument('--txtname', required=False, type=str, default="all_paths", help= "txt file with all paths")
parser.add_argument('--dirsuffix', required=False, type=str, default="", help="additional suffix for base directory name")
parser.add_argument('--topdir', required=False, type=str, default="CMSSW_11_2_0_phase2", help="base directory name")
parser.add_argument('--norm', required=False, type=bool, default=False, help="Draw normalized plots")
args = parser.parse_args()

# ****** dictionary for in/out directories ******
dir_dict = {}
dir_dict["local"]={}
#dir_dict["local"]["input"]="/Users/valeriadamante/Desktop/Dottorato/public/CMSSW_11_2_0_pre9/src/RootFiles/"
dir_dict["local"]["input"] = "/Users/valeriadamante/Desktop/Dottorato/public/CMSSW_11_3_0_pre4/src/Validate_files/output/"
dir_dict["local"]["output"] = "/Users/valeriadamante/Desktop/Dottorato/public/CMSSW_11_3_0_pre4/src/ValidationPlots/"
#dir_dict["local"]["output"]="/Users/valeriadamante/Desktop/Dottorato/public/CMSSW_11_2_0_pre9/src/"

dir_dict["lxplus"]={}
dir_dict["lxplus"]["input"]="/afs/cern.ch/work/v/vdamante/public/CMSSW_11_3_0_pre4/src/Validate_files/output/"
dir_dict["lxplus"]["output"]="/afs/cern.ch/work/v/vdamante/public/CMSSW_11_3_0_pre4/src/ValidationPlots/"#"/eos/home-v/vdamante/www/phase2validation/"

input_dir = dir_dict[args.machine]["input"]
out_dir =  dir_dict[args.machine]['output']

# **** create txt file with all paths  ****
input_files_str = args.reference+","+args.targets
if not os.path.exists(args.txtname+".txt"):
    command="root -l -b make_histogram.cpp+O\(\\\""+input_files_str+"\\\",\\\""+args.txtname+"\\\"\)"
    print command
    os.system(command)

# **** create list of files ****
input_files = []
input_files.append(input_dir+args.reference)
for file in (input_dir+args.targets).split(","):
    input_files.append(file)

# **** initialize the class and its values ****
make_histogram = histomaker(input_files)
make_histogram.file_name=args.txtname+".txt"                            # file name with paths
make_histogram.top_level_dir = out_dir+args.topdir+args.dirsuffix       # base directory name
make_histogram.norm=args.norm
# **** eventually initialize variables ****
variables =[]
if args.variables != "":
    for var in (args.variables).split(","):
        variables.append(var)
    #print variables
    make_histogram.variables = variables

# **** draw histograms ****
