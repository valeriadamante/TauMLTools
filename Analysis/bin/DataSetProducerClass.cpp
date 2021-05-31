#include "/TauMLTools/Analysis/interface/DataSetProducer.h"

DataSetProducer::DataSetProducer(){
  i=0;
  //absolute_path = "/Users/valeriadamante/Desktop/Dottorato/gridui/L2SkimmedTuples";
  //DYFile.push_back(absolute_path+"all_DY.root");
  //VBFFile.push_back(absolute_path+"all_TT.root");
  //TTFile.push_back(absolute_path+"all_VBF.root");
  //WJFile.push_back(absolute_path+"all_WJetsToLNu.root");
  //QCDFile.push_back(absolute_path+"all_QCD.root");
  //RDF_Map[DYFile]=RDF("L2TauTrainTuple", DYFile);
}/*
DataSetProducer::DataSetProducer(const std::vector<std::string> signals, const std::vector<std::string> backgrounds){
    ROOT::RDataFrame dfSignal("L2TauTrainTuple", signals);
    ROOT::RDataFrame dfBackground("L2TauTrainTuple", backgrounds);

}
DataSetProducer::DataSetProducer(const std::vector<std::string> signals, const std::vector<std::string> backgrounds, const std::vector<std::string> backgroundsFiltered){
    ROOT::RDataFrame dfSignal("L2TauTrainTuple", signals);
    ROOT::RDataFrame dfBackground("L2TauTrainTuple", backgrounds);
    ROOT::RDataFrame dfBackgroundFiltered("L2TauTrainTuple", backgroundsFiltered);

}
DataSetProducer::DataSetProducer(const std::vector<std::string> DYFile,const std::vector<std::string> VBFFile,const std::vector<std::string> TTFile,const std::vector<std::string> WJFile,const std::vector<std::string> QCDFile){
        RDF df_DY("L2TauTrainTuple",DYFile);
        RDF df_VBF("L2TauTrainTuple",VBFFile);
        RDF df_TT("L2TauTrainTuple",TTFile);
        RDF df_WJetsToLNu("L2TauTrainTuple",WJFile);
        RDF dfQCD("L2TauTrainTuple",QCDFile);
        RDF df_DataSetTraining("L2TauTrainTuple","DataSetTraining.root");
}/*
void {

    RDF df_DY("L2TauTrainTuple",absolute_path+"all_DY.root");
    RDF df_TT("L2TauTrainTuple",absolute_path+"all_TT.root");
    RDF df_VBF("L2TauTrainTuple",absolute_path+"all_VBF.root");
    RDF df_WJetsToLNu("L2TauTrainTuple",absolute_path+"all_WJetsToLNu.root");
    RDF dfQCD("L2TauTrainTuple",absolute_path+"all_QCD.root");
    RDF dfQCDFiltered("L2TauTrainTuple",absolute_path+"QCDFiltered.root");
    RDF df_DataSetTraining("L2TauTrainTuple",absolute_path+"DataSetTraining.root");
}
*/
