#include <memory>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "ROOT/RDataFrame.hxx"
#include <boost/bimap.hpp>
#include "TauMLTools/Core/interface/RootExt.h"
#include "TauMLTools/Core/interface/program_main.h"

struct Arguments {
    run::Argument<std::string> input{"input", "input root file with tau tuple"};
    //run::Argument<std::string> output{"output", "output root file with training tuple"};
    //run::Argument<unsigned> n_threads{"n-threads", "number of threads", 1};
    //run::Argument<Long64_t> start_entry{"start-entry", "start entry", 0};
    //run::Argument<Long64_t> end_entry{"end-entry", "end entry", std::numeric_limits<Long64_t>::max()};
    //run::Argument<float> training_weight_factor{"training-weight-factor",
    //    "additional factor to the normalization of the training weights", 4.f};
    //run::Argument<int> parity{"parity", "take odd (parity=1), even (parity=0) or all (parity=-1) events", -1};
    //run::Argument<bool> isQCD{"isQCD", "is monte carlo sample", false};
};
// Semplice esempio di una classe C++
class HistoMaker
{
public:
    using DataSetMap = std::map<std::string, std::map<std::string, std::string>>;
    //using DatasetBiMap = boost::bimap<std::string, unsigned>;
    //using FileNames = std::set<std::string>;
    using RDF = ROOT::RDF::RNode;
    //int numerator = 0;
    //int denominator = 0;
    //TH1D *h_num;
    //TH1D *h_den;


    // costruttore
    HistoMaker(const Arguments& _args) ;
    // metodo
    void Run();

private:

  // metodo privato aprire files
  //  static std::shared_ptr<TFile> OpenFile(const std::string& file_name);
  // metodo privato aprire TTrees
  //  static std::shared_ptr<TTree> ReadTree(const std::shared_ptr<TFile>& file, const std::string& tree_name);


// quantita private
private:
    const Arguments args;
  //  std::shared_ptr<TFile> inputFile;
  //  std::shared_ptr<TTree> inputTree;
    //ROOT::RDataFrame dataFrame;
    //RDF df;
    //DatasetBiMap known_datasets;
    //FileNames input_files;
    std::string parent_directory;
    std::string parent_directory_2;
    DataSetMap datasets;
    //std::vector<float> pt_thresholds;
    //std::map<std::string, std::string> branch_types;
};
