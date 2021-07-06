#include "NewTuple.h"
#include "NewTuple.C"

void ValidateTuple(){
  NewTuple nt;
  nt.SetAbsolutePath("/Users/valeriadamante/Desktop/Dottorato/L2SkimmedTuples/DataSetTraining/");

  nt.SetPlotDir("/Users/valeriadamante/Desktop/Dottorato/plots/");

  //std::string plotDir = "/home/users/damante/CMSSW_11_2_1_Patatrack/src/plots/";
  //nt.SetAbsolutePath("/home/Valeria/Desktop/Dottorato/L2SkimmedTuples/DataSetTraining/");
  //nt.SetAbsolutePath("/home/users/damante/L2SkimmedTuples/DataSetTraining/");

  nt.SetOldTupleFile({nt.GetAbsolutePath()+"DataSetTraining.root"});
  nt.SetNewTupleFile({nt.GetAbsolutePath()+"miniTuple.root"});
  nt.GetWeight("Eph1_ForWeighting.root");
  nt.DrawOnlyHistos("Eph1_ForWeighting.root");
}
