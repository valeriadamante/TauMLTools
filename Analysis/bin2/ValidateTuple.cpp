#include "NewTuple.h"
#include "NewTuple.C"

void ValidateTuple(){
  NewTuple nt;
  nt.SetAbsolutePath("/Users/valeriadamante/Desktop/Dottorato/L2SkimmedTuples/DataSetTraining/");
  //nt.SetAbsolutePath("/home/Valeria/Desktop/Dottorato/L2SkimmedTuples/DataSetTraining/");
  nt.SetOldTupleFile({nt.GetAbsolutePath()+"DataSetTraining.root"});
  nt.SetNewTupleFile({nt.GetAbsolutePath()+"miniTuple.root"});
  //nt.GetWeight("all_Data.root");
  nt.DrawOnlyHistos("all_Data.root");
}
