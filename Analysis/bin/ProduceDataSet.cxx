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
    run::Argument<Long64_t> max_iterations{"max-iterations", "max iterations", 12152};
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
    //std::vector<int> all_events = {26362, 139019, 182810, 405918}; // 1st iteration
    // std::vector<int> all_events = {26363, 139059, 182811, 407783}; // 2nd iteration
    // TT - DY - WJ - QCD
    std::vector<int> all_events = {258559, 473078, 78699, 1620224}; // 3rd iteration
    //td::vector<int> all_events_norm = {7, 37, 49, 107}; // 1st and 2nd iterations
    std::vector<int> all_events_norm = {21, 39, 7, 133}; // 3rd iteration
    const int filling_WJ = 7; // 7 ;
    const int filling_TT = 21; // 37 ;
    const int filling_DY = 39; // 49 ;
    const int filling_QCD = 133; //107 ;

    DataSetProducer(const Arguments& _args) :
        args(_args), WJetFile(root_ext::OpenRootFile(absolute_path+"all_WJetsToLNu.root")), TTFile(root_ext::OpenRootFile(absolute_path+"all_TT.root")), /*ZPrimeFile(root_ext::OpenRootFile(absolute_path+"all_ZPrime.root")),*/ DYFile(root_ext::OpenRootFile(absolute_path+"all_DY.root")), QCDFile(root_ext::OpenRootFile(absolute_path+"QCDFiltered.root")),outputFile(root_ext::CreateRootFile(absolute_path+args.output(), ROOT::kLZ4, 4)), WJTuple(WJetFile.get(), true), TTTuple( TTFile.get(), true), /*ZPrimeTuple( ZPrimeFile.get(), true),*/ DYTuple( DYFile.get(), true), QCDTuple( QCDFile.get(), true), outputTuple(outputFile.get(), false)
    {

        if(args.n_threads() > 1)
            ROOT::EnableImplicitMT(args.n_threads());

    }
    void Run()
    {
        const Long64_t first_iteration =0;
        //const Long64_t n_iterations = 3770; // 1st iteration
        //const Long64_t n_iterations = 3780; // 2nd iteration
        const Long64_t n_iterations = 12152; // 3rd iteration
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
            std::cout << "superate entries per WJ" << std::endl ;
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
            std::cout << "superate entries per TT " << std::endl ;
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
            std::cout << "superate entries per DY" << std::endl ;
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
    #define CP_BR(name) outputTuple().name = tau.name;
    void CopyBranches(const TrainingTau& tau)
    {
        auto& out = outputTuple();
        CP_BR(evt);
        CP_BR(run);
        CP_BR(lumi);
        CP_BR(defaultDiTauPath_lastModuleIndex);
        CP_BR(defaultDiTauPath_result);
        CP_BR(genEventWeight);
        CP_BR(sampleType);
        CP_BR(l1Tau_index);
        CP_BR(l1Tau_pt);
        CP_BR(l1Tau_eta);
        CP_BR(l1Tau_phi);
        CP_BR(l1Tau_mass);
        CP_BR(l1Tau_hwIso);
        CP_BR(l1Tau_hwQual);
        CP_BR(l1Tau_towerIEta);
        CP_BR(l1Tau_towerIPhi);
        CP_BR(genLepton_index);
        CP_BR(genLepton_kind);
        CP_BR(genLepton_charge);
        CP_BR(genLepton_vis_eta);
        CP_BR(genLepton_vis_pt);
        CP_BR(genLepton_vis_mass);
        CP_BR(genLepton_vis_phi);
        CP_BR(l1Tau_rawEt);
        CP_BR(l1Tau_isoEt);
        CP_BR(l1Tau_hasEM);
        CP_BR(l1Tau_isMerged);

        for(long unsigned int i=0; i < tau.caloRecHit_eb_rho.size() ; i++){
          out.caloRecHit_eb_rho.push_back(static_cast<float>(tau.caloRecHit_eb_rho.at(i)));
          out.caloRecHit_eb_eta.push_back(static_cast<float>(tau.caloRecHit_eb_eta.at(i)));
          out.caloRecHit_eb_phi.push_back(static_cast<float>(tau.caloRecHit_eb_phi.at(i)));
          out.caloRecHit_eb_energy.push_back(static_cast<float>(tau.caloRecHit_eb_energy.at(i)));
          out.caloRecHit_eb_time.push_back(static_cast<float>(tau.caloRecHit_eb_time.at(i)));
          out.caloRecHit_eb_detId.push_back(static_cast<ulong>(tau.caloRecHit_eb_detId.at(i)));
          out.caloRecHit_eb_chi2.push_back(static_cast<float>(tau.caloRecHit_eb_chi2.at(i)));
          out.caloRecHit_eb_energyError.push_back(static_cast<float>(tau.caloRecHit_eb_energyError.at(i)));
          out.caloRecHit_eb_timeError.push_back(static_cast<float>(tau.caloRecHit_eb_timeError.at(i)));
          out.caloRecHit_eb_flagsBits.push_back(static_cast<uint32_t>(tau.caloRecHit_eb_flagsBits.at(i)));
          out.caloRecHit_eb_isRecovered.push_back(static_cast<bool>(tau.caloRecHit_eb_isRecovered.at(i)));
          out.caloRecHit_eb_isTimeValid.push_back(static_cast<bool>(tau.caloRecHit_eb_isTimeValid.at(i)));
          out.caloRecHit_eb_isTimeErrorValid.push_back(static_cast<bool>(tau.caloRecHit_eb_isTimeErrorValid.at(i)));
        }
        for(long unsigned int i=0; i < tau.caloRecHit_ee_rho.size() ; i++){
          out.caloRecHit_ee_rho.push_back(static_cast<float>(tau.caloRecHit_ee_rho.at(i)));
          out.caloRecHit_ee_eta.push_back(static_cast<float>(tau.caloRecHit_ee_eta.at(i)));
          out.caloRecHit_ee_phi.push_back(static_cast<float>(tau.caloRecHit_ee_phi.at(i)));
          out.caloRecHit_ee_energy.push_back(static_cast<float>(tau.caloRecHit_ee_energy.at(i)));
          out.caloRecHit_ee_time.push_back(static_cast<float>(tau.caloRecHit_ee_time.at(i)));
          out.caloRecHit_ee_detId.push_back(static_cast<ulong>(tau.caloRecHit_ee_detId.at(i)));
          out.caloRecHit_ee_chi2.push_back(static_cast<float>(tau.caloRecHit_ee_chi2.at(i)));
          out.caloRecHit_ee_energyError.push_back(static_cast<float>(tau.caloRecHit_ee_energyError.at(i)));
          out.caloRecHit_ee_timeError.push_back(static_cast<float>(tau.caloRecHit_ee_timeError.at(i)));
          out.caloRecHit_ee_flagsBits.push_back(static_cast<uint32_t>(tau.caloRecHit_ee_flagsBits.at(i)));
          out.caloRecHit_ee_isRecovered.push_back(static_cast<bool>(tau.caloRecHit_ee_isRecovered.at(i)));
          out.caloRecHit_ee_isTimeValid.push_back(static_cast<bool>(tau.caloRecHit_ee_isTimeValid.at(i)));
          out.caloRecHit_ee_isTimeErrorValid.push_back(static_cast<bool>(tau.caloRecHit_ee_isTimeErrorValid.at(i)));
        }
        for(long unsigned int i=0; i < tau.caloRecHit_hbhe_rho.size() ; i++){
          out.caloRecHit_hbhe_rho.push_back(static_cast<float>(tau.caloRecHit_hbhe_rho.at(i)));
          out.caloRecHit_hbhe_eta.push_back(static_cast<float>(tau.caloRecHit_hbhe_eta.at(i)));
          out.caloRecHit_hbhe_phi.push_back(static_cast<float>(tau.caloRecHit_hbhe_phi.at(i)));
          out.caloRecHit_hbhe_energy.push_back(static_cast<float>(tau.caloRecHit_hbhe_energy.at(i)));
          out.caloRecHit_hbhe_time.push_back(static_cast<float>(tau.caloRecHit_hbhe_time.at(i)));
          out.caloRecHit_hbhe_detId.push_back(static_cast<ulong>(tau.caloRecHit_hbhe_detId.at(i)));
          out.caloRecHit_hbhe_chi2.push_back(static_cast<float>(tau.caloRecHit_hbhe_chi2.at(i)));
          out.caloRecHit_hbhe_flags.push_back(static_cast<ulong>(tau.caloRecHit_hbhe_flags.at(i)));
        }
        for(long unsigned int i=0; i < tau.caloRecHit_ho_rho.size() ; i++){
        out.caloRecHit_ho_rho.push_back(static_cast<float>(tau.caloRecHit_ho_rho.at(i)));
        out.caloRecHit_ho_eta.push_back(static_cast<float>(tau.caloRecHit_ho_eta.at(i)));
        out.caloRecHit_ho_phi.push_back(static_cast<float>(tau.caloRecHit_ho_phi.at(i)));
        out.caloRecHit_ho_energy.push_back(static_cast<float>(tau.caloRecHit_ho_energy.at(i)));
        out.caloRecHit_ho_time.push_back(static_cast<float>(tau.caloRecHit_ho_time.at(i)));
        out.caloRecHit_ho_detId.push_back(static_cast<ulong>(tau.caloRecHit_ho_detId.at(i)));
        out.caloRecHit_ho_aux.push_back(static_cast<ulong>(tau.caloRecHit_ho_aux.at(i)));
        out.caloRecHit_ho_flags.push_back(static_cast<ulong>(tau.caloRecHit_ho_flags.at(i)));
      }
        for(long unsigned int i=0; i < tau.caloRecHit_hf_rho.size() ; i++){
        out.caloRecHit_hf_rho.push_back(static_cast<float>(tau.caloRecHit_hf_rho.at(i)));
        out.caloRecHit_hf_eta.push_back(static_cast<float>(tau.caloRecHit_hf_eta.at(i)));
        out.caloRecHit_hf_phi.push_back(static_cast<float>(tau.caloRecHit_hf_phi.at(i)));
        out.caloRecHit_hf_energy.push_back(static_cast<float>(tau.caloRecHit_hf_energy.at(i)));
        out.caloRecHit_hf_time.push_back(static_cast<float>(tau.caloRecHit_hf_time.at(i)));
        out.caloRecHit_hf_detId.push_back(static_cast<ulong>(tau.caloRecHit_hf_detId.at(i)));
        out.caloRecHit_hf_flags.push_back(static_cast<ulong>(tau.caloRecHit_hf_flags.at(i)));
        out.caloRecHit_hf_timeFalling.push_back(static_cast<float>(tau.caloRecHit_hf_timeFalling.at(i)));
        out.caloRecHit_hf_auxHF.push_back(static_cast<uint32_t>(tau.caloRecHit_hf_auxHF.at(i)));
        out.caloRecHit_hf_aux.push_back(static_cast<ulong>(tau.caloRecHit_hf_aux.at(i)));
      }
      for(long unsigned int i=0; i < tau.patatrack_pt.size() ; i++){
        out.patatrack_pt.push_back(static_cast<float>(tau.patatrack_pt.at(i)));
        out.patatrack_eta.push_back(static_cast<float>(tau.patatrack_eta.at(i)));
        out.patatrack_phi.push_back(static_cast<float>(tau.patatrack_phi.at(i)));
        out.patatrack_chi2.push_back(static_cast<float>(tau.patatrack_chi2.at(i)));
        out.patatrack_ndof.push_back(static_cast<float>(tau.patatrack_ndof.at(i)));
        out.patatrack_charge.push_back(static_cast<int>(tau.patatrack_charge.at(i)));
        out.patatrack_quality.push_back(static_cast<uint>(tau.patatrack_quality.at(i)));
        out.patatrack_dxy.push_back(static_cast<float>(tau.patatrack_dxy.at(i)));
        out.patatrack_dz.push_back(static_cast<float>(tau.patatrack_dz.at(i)));
        out.patatrack_vertex_id.push_back(static_cast<float>(tau.patatrack_vertex_id.at(i)));
      }
      for(long unsigned int i=0; i < tau.patavert_z.size() ; i++){
        out.patavert_z.push_back(static_cast<float>(tau.patavert_z.at(i)));
        out.patavert_weight.push_back(static_cast<float>(tau.patavert_weight.at(i)));
        out.patavert_ptv2.push_back(static_cast<float>(tau.patavert_ptv2.at(i)));
        out.patavert_chi2.push_back(static_cast<float>(tau.patavert_chi2.at(i)));
        out.patavert_ndof.push_back(static_cast<int>(tau.patavert_ndof.at(i)));
      }
        outputTuple.Fill();
    }
    #undef CP_BR


private:
    const Arguments args;
    std::shared_ptr<TFile>  WJetFile, TTFile, /*ZPrimeFile,*/ DYFile, QCDFile, outputFile;
    TrainingTauTuple WJTuple, TTTuple, /*ZPrimeTuple,*/ DYTuple, QCDTuple, outputTuple;

};

} // namespace analysis

PROGRAM_MAIN(analysis::DataSetProducer, Arguments)
