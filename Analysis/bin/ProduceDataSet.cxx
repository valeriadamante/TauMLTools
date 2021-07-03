/*! Produce training tuple from tau tuple.
*/

#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/variadic.hpp>
#include <boost/math/constants/constants.hpp>

#include "TauMLTools/Core/interface/program_main.h"
#include "TauMLTools/Core/interface/AnalysisMath.h"
#include "TauMLTools/Core/interface/RootExt.h"
#include "TauMLTools/Analysis/interface/L2EventTuple.h"
#include "TauMLTools/Analysis/interface/L2TrainingTauTuple.h"
#include "TauMLTools/Core/interface/ProgressReporter.h"


struct Arguments {
    //run::Argument<std::string> input{"input", "input root file with tau tuple"};
    run::Argument<std::string> output{"output", "output root file with training tuple"};
    run::Argument<unsigned> n_threads{"n-threads", "number of threads", 1};
    run::Argument<Long64_t> first_iteration{"first-iteration", "first iterations", 0};
    run::Argument<Long64_t> max_iterations{"max-iterations", "max iterations", 3780};
    //run::Argument<Long64_t> end_entry{"end-entry", "end entry", std::numeric_limits<Long64_t>::max()};
    //run::Argument<float> training_weight_factor{"training-weight-factor",
    //    "additional factor to the normalization of the training weights", 4.f};
    //run::Argument<int> parity{"parity", "take odd (parity=1), even (parity=0) or all (parity=-1) events", -1};
    //run::Argument<bool> isQCD{"isQCD", "is monte carlo sample", false};
};


namespace analysis {

std::string absolute_path = "/home/users/damante/L2SkimmedTuples/DataSetTraining/";

class DataSetProducer {
public:
    using TrainingTau = tau_train_tuple::Tau;
    using TrainingTauTuple = tau_train_tuple::TauTrainTuple;
    //std::vector<std::string> inputFiles = {absolute_path+"all_WJetsToLNu.root", absolute_path+"all_TT.root" , absolute_path+"all_DY.root", absolute_path+"all_ZPrime.root", absolute_path+"QCDFiltered.root" };
    //std::vector<int> all_events = {26362, 139019, 182810, 405918};
    std::vector<int> all_events = {26363, 139059, 182811, 407783};
    std::vector<int> all_events_norm = {7, 37, 49, 107};
    const int filling_WJ = 7 ;
    const int filling_TT = 37 ;
    const int filling_DY = 49 ;
    //const int filling_ZP = 30 ;
    const int filling_QCD = 107 ;

    DataSetProducer(const Arguments& _args) :
        args(_args), WJetFile(root_ext::OpenRootFile(absolute_path+"all_WJetsToLNu.root")), TTFile(root_ext::OpenRootFile(absolute_path+"all_TT.root")), /*ZPrimeFile(root_ext::OpenRootFile(absolute_path+"all_ZPrime.root")),*/ DYFile(root_ext::OpenRootFile(absolute_path+"all_DY.root")), QCDFile(root_ext::OpenRootFile(absolute_path+"QCDFiltered.root")),outputFile(root_ext::CreateRootFile(absolute_path+args.output(), ROOT::kLZ4, 4)), WJTuple(WJetFile.get(), true), TTTuple( TTFile.get(), true), /*ZPrimeTuple( ZPrimeFile.get(), true),*/ DYTuple( DYFile.get(), true), QCDTuple( QCDFile.get(), true), outputTuple(outputFile.get(), false)
    {

        if(args.n_threads() > 1)
            ROOT::EnableImplicitMT(args.n_threads());

    }
    void Run()
    {
        const Long64_t first_iteration =0;
        const Long64_t n_iterations = 3780;
        int filling_counter = std::max(first_iteration, args.first_iteration());
        int WJ_counter =std::max(first_iteration, args.first_iteration()*filling_WJ);
        int TT_counter =std::max(first_iteration, args.first_iteration()*filling_TT);
        int DY_counter =std::max(first_iteration, args.first_iteration()*filling_DY);
        //int ZP_counter =std::max(first_iteration, args.first_iteration()*filling_ZP);
        int QCD_counter =std::max(first_iteration, args.first_iteration()*filling_QCD);
        int filling_counter_max = std::min(n_iterations, args.max_iterations());
        tools::ProgressReporter reporter(10, std::cout, "Creating training tuple...");
        reporter.SetTotalNumberOfEvents(filling_counter_max);
        while(filling_counter < filling_counter_max){
          /* WJets */
          if(WJ_counter>=WJTuple.GetEntries()){
            std::cout << "superate entries" << std::endl ;
            WJ_counter =0;
          }
          WriteTuple(WJTuple, WJ_counter, filling_WJ);
          WJ_counter += filling_WJ;
          //std::cout << "WJ counter " <<WJ_counter <<std::endl;
          //std::cout << "fill count " <<filling_counter <<std::endl;
          // if(filling_counter >=3000)
            // std::cout << "WJ filled " << std::endl;

          /* TT */
          if(TT_counter>=TTTuple.GetEntries()){
            std::cout << "superate entries" << std::endl ;
            TT_counter =0;
          }
          WriteTuple(TTTuple, TT_counter, filling_TT);
          TT_counter += filling_TT;
          //std::cout << "TT counter " <<TT_counter <<std::endl;
          //std::cout << "fill count " <<filling_counter <<std::endl;
          // if(filling_counter >=3000)
            // std::cout << "TT filled " << std::endl;

          /* DY */
          if(DY_counter>=DYTuple.GetEntries()){
            std::cout << "superate entries" << std::endl ;
            DY_counter =0;
          }
          WriteTuple(DYTuple, DY_counter, filling_DY);
          DY_counter += filling_DY;

          //std::cout << "DY counter " <<DY_counter <<std::endl;
          //std::cout << "fill count " <<filling_counter <<std::endl;
          // if(filling_counter >=3000)
            // std::cout << "DY filled " << std::endl;
          /* QCD */
          if(QCD_counter>=QCDTuple.GetEntries()){
            break;
          }
          WriteTuple(QCDTuple, QCD_counter, filling_QCD);
          QCD_counter += filling_QCD;
          //std::cout << "QCD counter " <<QCD_counter <<std::endl;
          //std::cout << "fill count " <<filling_counter <<std::endl;
          // if(filling_counter >=3000)
            // std::cout << "QCD filled " << std::endl;

          if(++filling_counter % 100 == 0)
              reporter.Report(filling_counter);
        }

        outputTuple.Write();
        std::cout << "Training tuples has been successfully stored in " << args.output() << "." << std::endl;
    }

private :

    void WriteTuple(TrainingTauTuple &Tuple, int Counter, int Fill){
      if(Counter+Fill>Tuple.GetEntries())
        return;
      for(Long64_t current_entry = Counter; current_entry < Counter+Fill; ++current_entry) {
          //std::cout << "entry = " << current_entry << std::endl;
          Tuple.GetEntry(current_entry);
          //std::cout << "TTree entries == " <<TTuple.GetEntries()<<std::endl;
          const auto& tau = Tuple.data();
          CopyBranches(tau);
      }
    }

    void CopyBranches(const TrainingTau& tau)
    {
        auto& out = outputTuple();
        out.run = tau.run;
        out.lumi = tau.lumi;
        out.evt = tau.evt;
        out.genEventWeight = tau.genEventWeight;
        out.sampleType= tau.sampleType;
        out.genLepton_index= tau.genLepton_index;
        out.genLepton_kind= tau.genLepton_kind;
        out.genLepton_charge= tau.genLepton_charge;
        out.l1Tau_index= tau.l1Tau_index;
        out.genLepton_vis_pt= tau.genLepton_vis_pt;
        out.genLepton_vis_eta= tau.genLepton_vis_eta;
        out.genLepton_vis_phi= tau.genLepton_vis_phi;
        out.genLepton_vis_mass= tau.genLepton_vis_mass;
        out.l1Tau_pt= tau.l1Tau_pt;
        out.l1Tau_eta= tau.l1Tau_eta;
        out.l1Tau_phi= tau.l1Tau_phi;
        out.l1Tau_mass= tau.l1Tau_mass;
        out.l1Tau_hwIso= tau.l1Tau_hwIso;
        out.l1Tau_hwQual= tau.l1Tau_hwQual;
        out.l1Tau_towerIEta= tau.l1Tau_towerIEta;
        out.l1Tau_towerIPhi= tau.l1Tau_towerIPhi;
        out.l1Tau_rawEt= tau.l1Tau_rawEt;
        out.l1Tau_isoEt= tau.l1Tau_isoEt;
        out.l1Tau_hasEM= tau.l1Tau_hasEM;
        out.l1Tau_isMerged= tau.l1Tau_isMerged;
        out.caloRecHit_eb_rho= tau.caloRecHit_eb_rho;
        out.caloRecHit_ee_rho= tau.caloRecHit_ee_rho;
        out.caloRecHit_eb_eta= tau.caloRecHit_eb_eta;
        out.caloRecHit_ee_eta= tau.caloRecHit_ee_eta;
        out.caloRecHit_eb_phi= tau.caloRecHit_eb_phi;
        out.caloRecHit_ee_phi= tau.caloRecHit_ee_phi;
        out.caloRecHit_eb_energy= tau.caloRecHit_eb_energy;
        out.caloRecHit_ee_energy= tau.caloRecHit_ee_energy;
        out.caloRecHit_eb_time= tau.caloRecHit_eb_time;
        out.caloRecHit_ee_time= tau.caloRecHit_ee_time;
        out.caloRecHit_eb_detId= tau.caloRecHit_eb_detId;
        out.caloRecHit_ee_detId= tau.caloRecHit_ee_detId;
        out.caloRecHit_eb_chi2= tau.caloRecHit_eb_chi2;
        out.caloRecHit_ee_chi2= tau.caloRecHit_ee_chi2;
        out.caloRecHit_eb_energyError= tau.caloRecHit_eb_energyError;
        out.caloRecHit_ee_energyError= tau.caloRecHit_ee_energyError;
        out.caloRecHit_eb_timeError= tau.caloRecHit_eb_timeError;
        out.caloRecHit_ee_timeError= tau.caloRecHit_ee_timeError;
        out.caloRecHit_eb_flagsBits= tau.caloRecHit_eb_flagsBits;
        out.caloRecHit_ee_flagsBits= tau.caloRecHit_ee_flagsBits;
        out.caloRecHit_eb_isRecovered= tau.caloRecHit_eb_isRecovered;
        out.caloRecHit_ee_isRecovered= tau.caloRecHit_ee_isRecovered;
        out.caloRecHit_eb_isTimeValid= tau.caloRecHit_eb_isTimeValid;
        out.caloRecHit_ee_isTimeValid= tau.caloRecHit_ee_isTimeValid;
        out.caloRecHit_eb_isTimeErrorValid= tau.caloRecHit_eb_isTimeErrorValid;
        out.caloRecHit_ee_isTimeErrorValid= tau.caloRecHit_ee_isTimeErrorValid;
        out.caloRecHit_hbhe_rho= tau.caloRecHit_hbhe_rho;
        out.caloRecHit_hbhe_eta= tau.caloRecHit_hbhe_eta;
        out.caloRecHit_hbhe_phi= tau.caloRecHit_hbhe_phi;
        out.caloRecHit_hbhe_energy= tau.caloRecHit_hbhe_energy;
        out.caloRecHit_hbhe_time= tau.caloRecHit_hbhe_time;
        out.caloRecHit_hbhe_detId= tau.caloRecHit_hbhe_detId;
        out.caloRecHit_hbhe_chi2= tau.caloRecHit_hbhe_chi2;
        out.caloRecHit_hbhe_flags= tau.caloRecHit_hbhe_flags;
        out.caloRecHit_ho_rho= tau.caloRecHit_ho_rho;
        out.caloRecHit_ho_eta= tau.caloRecHit_ho_eta;
        out.caloRecHit_ho_phi= tau.caloRecHit_ho_phi;
        out.caloRecHit_ho_energy= tau.caloRecHit_ho_energy;
        out.caloRecHit_ho_time= tau.caloRecHit_ho_time;
        out.caloRecHit_ho_detId= tau.caloRecHit_ho_detId;
        out.caloRecHit_ho_aux= tau.caloRecHit_ho_aux;
        out.caloRecHit_ho_flags= tau.caloRecHit_ho_flags;
        out.caloRecHit_hf_rho= tau.caloRecHit_hf_rho;
        out.caloRecHit_hf_eta= tau.caloRecHit_hf_eta;
        out.caloRecHit_hf_phi= tau.caloRecHit_hf_phi;
        out.caloRecHit_hf_energy= tau.caloRecHit_hf_energy;
        out.caloRecHit_hf_time= tau.caloRecHit_hf_time;
        out.caloRecHit_hf_detId= tau.caloRecHit_hf_detId;
        out.caloRecHit_hf_flags= tau.caloRecHit_hf_flags;
        out.caloRecHit_hf_timeFalling= tau.caloRecHit_hf_timeFalling;
        out.caloRecHit_hf_auxHF= tau.caloRecHit_hf_auxHF;
        out.caloRecHit_hf_aux= tau.caloRecHit_hf_aux;
        out.patatrack_pt= tau.patatrack_pt;
        out.patatrack_eta= tau.patatrack_eta;
        out.patatrack_phi= tau.patatrack_phi;
        out.patatrack_chi2= tau.patatrack_chi2;
        out.patatrack_ndof= tau.patatrack_ndof;
        out.patatrack_charge= tau.patatrack_charge;
        out.patatrack_quality= tau.patatrack_quality;
        out.patatrack_dxy= tau.patatrack_dxy;
        out.patatrack_dz= tau.patatrack_dz;
        out.patatrack_vertex_id= tau.patatrack_vertex_id;
        out.patavert_z= tau.patavert_z;
        out.patavert_weight= tau.patavert_weight;
        out.patavert_ptv2= tau.patavert_ptv2;
        out.patavert_chi2= tau.patavert_chi2;
        out.patavert_ndof= tau.patavert_ndof;
        outputTuple.Fill();
    }


private:
    const Arguments args;
    std::shared_ptr<TFile>  WJetFile, TTFile, /*ZPrimeFile,*/ DYFile, QCDFile, outputFile;
    TrainingTauTuple WJTuple, TTTuple, /*ZPrimeTuple,*/ DYTuple, QCDTuple, outputTuple;

};

} // namespace analysis

PROGRAM_MAIN(analysis::DataSetProducer, Arguments)
