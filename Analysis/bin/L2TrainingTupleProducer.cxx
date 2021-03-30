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
#include "TauMLTools/Analysis/interface/GenLepton.h"

struct Arguments {
    run::Argument<std::string> input{"input", "input root file with tau tuple"};
    run::Argument<std::string> output{"output", "output root file with training tuple"};
    run::Argument<unsigned> n_threads{"n-threads", "number of threads", 1};
    run::Argument<Long64_t> start_entry{"start-entry", "start entry", 0};
    run::Argument<Long64_t> end_entry{"end-entry", "end entry", std::numeric_limits<Long64_t>::max()};
    run::Argument<float> training_weight_factor{"training-weight-factor",
        "additional factor to the normalization of the training weights", 4.f};
    run::Argument<int> parity{"parity", "take odd (parity=1), even (parity=0) or all (parity=-1) events", -1};
    run::Argument<bool> isQCD{"isQCD", "is monte carlo sample", false};
};

namespace analysis {


class L2TrainingTupleProducer {
public:
    using Tau = train_tuple::Tau;
    using L2TauTuple = train_tuple::TrainTuple;
    using TrainingTau = tau_train_tuple::Tau;
    using TrainingTauTuple = tau_train_tuple::TauTrainTuple;

    L2TrainingTupleProducer(const Arguments& _args) :
        args(_args), inputFile(root_ext::OpenRootFile(args.input())),
        outputFile(root_ext::CreateRootFile(args.output(), ROOT::kLZ4, 4)),
        l2tauTuple(inputFile.get(), true), trainingTauTuple(outputFile.get(), false)
    {
        if(args.n_threads() > 1)
            ROOT::EnableImplicitMT(args.n_threads());
    }

    void Run()
    {
        const Long64_t end_entry = std::min(l2tauTuple.GetEntries(), args.end_entry());
        size_t n_processed = 0, n_total = static_cast<size_t>(end_entry - args.start_entry());
        tools::ProgressReporter reporter(10, std::cout, "Creating training tuple...");
        reporter.SetTotalNumberOfEvents(n_total);
        for(Long64_t current_entry = args.start_entry(); current_entry < end_entry; ++current_entry) {
            // std::cout << "evento = " << current_entry << std::endl;
            l2tauTuple.GetEntry(current_entry);
            const auto& tau = l2tauTuple.data();
            FillTauBranches(tau, args.isQCD());
            trainingTauTuple.Fill();
            if(++n_processed % 1000 == 0)
                reporter.Report(n_processed);
        }
        reporter.Report(n_processed, true);

        trainingTauTuple.Write();
        std::cout << "Training tuples has been successfully stored in " << args.output() << "." << std::endl;
    }

private:
    static constexpr float pi = boost::math::constants::pi<float>();

    template<typename Scalar>
    static Scalar DeltaPhi(Scalar phi1, Scalar phi2)
    {
        static constexpr Scalar pi = boost::math::constants::pi<Scalar>();
        Scalar dphi = phi1 - phi2;
        if(dphi > pi)
            dphi -= 2*pi;
        else if(dphi <= -pi)
            dphi += 2*pi;
        return dphi;
    }

    float deltaR(Float_t phi1,Float_t eta1,Float_t phi2,Float_t eta2) {
      auto dphi = DeltaPhi(phi1, phi2);
      auto deta = eta1-eta2;
      //std::cout << "deltaR = " << (std::sqrt(deta * deta + dphi * dphi)) <<std::endl;
      return (std::sqrt(deta * deta + dphi * dphi));

     }


    std::vector<int> FindIndices(const std::string& prefix, float l1eta, float l1phi, float delta_R_threshold){
        std::vector<int> indices_vector;
        std::vector<float> obj_eta = l2tauTuple.get<std::vector<float>>(prefix+"eta");
        std::vector<float> obj_phi = l2tauTuple.get<std::vector<float>>(prefix+"phi");
        for(std::vector<float>::size_type i = 0; i < obj_eta.size(); i++){
            if(deltaR(obj_phi.at(i), obj_eta.at(i), l1phi, l1eta)< delta_R_threshold) indices_vector.push_back(i);
        }
        return indices_vector;
    }

    #define CP_BR(name) trainingTauTuple().name = tau.name;
    /*
    void FillL1Taus(sconst Tau& tau, TauTrain& out , std::vector<int> tauindices){

      for (auto& k : tauindices){
        CP_BR(evt);
        CP_BR(run);
        CP_BR(lumi);
        CP_BR(genEventWeight);
        CP_BR(sampleType);

        out.l1Tau_pt = tau.l1Tau_pt.at(k);
        out.l1Tau_eta = tau.l1Tau_eta.at(k);
        out.l1Tau_phi = tau.l1Tau_phi.at(k);
        out.l1Tau_mass = tau.l1Tau_mass.at(k);
        out.l1Tau_hwIso = tau.l1Tau_hwIso.at(k);
        out.l1Tau_hwQual = tau.l1Tau_hwQual.at(k);
        out.l1Tau_towerIEta = tau.l1Tau_towerIEta.at(k);
        out.l1Tau_towerIPhi = tau.l1Tau_towerIPhi.at(k);
        out.l1Tau_rawEt = tau.l1Tau_rawEt.at(k);
        out.l1Tau_isoEt = tau.l1Tau_isoEt.at(k);
        out.l1Tau_hasEM = tau.l1Tau_hasEM.at(k);
        out.l1Tau_isMerged = tau.l1Tau_isMerged.at(k);

        //Float_t phi1,Float_t eta1,Float_t phi2,Float_t eta2

        std::vector<int> caloRecHit_ee_indices = FindIndices("caloRecHit_ee_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
        for(auto j: caloRecHit_ee_indices ){
          out.caloRecHit_ee_rho.push_back(tau.caloRecHit_ee_rho.at(j));
          out.caloRecHit_ee_eta.push_back(tau.caloRecHit_ee_eta.at(j));
          out.caloRecHit_ee_phi.push_back(tau.caloRecHit_ee_phi.at(j));
          out.caloRecHit_ee_energy.push_back(tau.caloRecHit_ee_energy.at(j));
          out.caloRecHit_ee_time.push_back(tau.caloRecHit_ee_time.at(j));
          out.caloRecHit_ee_detId.push_back(tau.caloRecHit_ee_detId.at(j));
          out.caloRecHit_ee_chi2.push_back(tau.caloRecHit_ee_chi2.at(j));
          out.caloRecHit_ee_energyError.push_back(tau.caloRecHit_ee_energyError.at(j));
          out.caloRecHit_ee_timeError.push_back(tau.caloRecHit_ee_timeError.at(j));
          out.caloRecHit_ee_flagsBits.push_back(tau.caloRecHit_ee_flagsBits.at(j));
          out.caloRecHit_ee_isRecovered.push_back(tau.caloRecHit_ee_isRecovered.at(j));
          out.caloRecHit_ee_isTimeValid.push_back(tau.caloRecHit_ee_isTimeValid.at(j));
          out.caloRecHit_ee_isTimeErrorValid.push_back(tau.caloRecHit_ee_isTimeErrorValid.at(j));
        }

        std::vector<int> caloRecHit_eb_indices = FindIndices("caloRecHit_eb_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
        for(auto j: caloRecHit_eb_indices ){
          out.caloRecHit_eb_rho.push_back(tau.caloRecHit_eb_rho.at(j));
          out.caloRecHit_eb_eta.push_back(tau.caloRecHit_eb_eta.at(j));
          out.caloRecHit_eb_phi.push_back(tau.caloRecHit_eb_phi.at(j));
          out.caloRecHit_eb_energy.push_back(tau.caloRecHit_eb_energy.at(j));
          out.caloRecHit_eb_time.push_back(tau.caloRecHit_eb_time.at(j));
          out.caloRecHit_eb_detId.push_back(tau.caloRecHit_eb_detId.at(j));
          out.caloRecHit_eb_chi2.push_back(tau.caloRecHit_eb_chi2.at(j));
          out.caloRecHit_eb_energyError.push_back(tau.caloRecHit_eb_energyError.at(j));
          out.caloRecHit_eb_timeError.push_back(tau.caloRecHit_eb_timeError.at(j));
          out.caloRecHit_eb_flagsBits.push_back(tau.caloRecHit_eb_flagsBits.at(j));
          out.caloRecHit_eb_isRecovered.push_back(tau.caloRecHit_eb_isRecovered.at(j));
          out.caloRecHit_eb_isTimeValid.push_back(tau.caloRecHit_eb_isTimeValid.at(j));
          out.caloRecHit_eb_isTimeErrorValid.push_back(tau.caloRecHit_eb_isTimeErrorValid.at(j));
        }

        std::vector<int> caloRecHit_hbhe_indices = FindIndices("caloRecHit_hbhe_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
        for(auto j: caloRecHit_hbhe_indices ){
          out.caloRecHit_hbhe_rho.push_back(tau.caloRecHit_hbhe_rho.at(j));
          out.caloRecHit_hbhe_eta.push_back(tau.caloRecHit_hbhe_eta.at(j));
          out.caloRecHit_hbhe_phi.push_back(tau.caloRecHit_hbhe_phi.at(j));
          out.caloRecHit_hbhe_energy.push_back(tau.caloRecHit_hbhe_energy.at(j));
          out.caloRecHit_hbhe_time.push_back(tau.caloRecHit_hbhe_time.at(j));
          out.caloRecHit_hbhe_detId.push_back(tau.caloRecHit_hbhe_detId.at(j));
          out.caloRecHit_hbhe_chi2.push_back(tau.caloRecHit_hbhe_chi2.at(j));
          out.caloRecHit_hbhe_flags.push_back(tau.caloRecHit_hbhe_flags.at(j));
          out.caloRecHit_hbhe_eraw.push_back(tau.caloRecHit_hbhe_eraw.at(j));
          out.caloRecHit_hbhe_eaux.push_back(tau.caloRecHit_hbhe_eaux.at(j));
          out.caloRecHit_hbhe_timeFalling.push_back(tau.caloRecHit_hbhe_timeFalling.at(j));
          out.caloRecHit_hbhe_idFront.push_back(tau.caloRecHit_hbhe_idFront.at(j));
          out.caloRecHit_hbhe_rho_front.push_back(tau.caloRecHit_hbhe_rho_front.at(j));
          out.caloRecHit_hbhe_eta_front.push_back(tau.caloRecHit_hbhe_eta_front.at(j));
          out.caloRecHit_hbhe_phi_front.push_back(tau.caloRecHit_hbhe_phi_front.at(j));
          out.caloRecHit_hbhe_auxHBHE.push_back(tau.caloRecHit_hbhe_auxHBHE.at(j));
          out.caloRecHit_hbhe_auxPhase1.push_back(tau.caloRecHit_hbhe_auxPhase1.at(j));
          out.caloRecHit_hbhe_auxTDC.push_back(tau.caloRecHit_hbhe_auxTDC.at(j));
          out.caloRecHit_hbhe_isMerged.push_back(tau.caloRecHit_hbhe_isMerged.at(j));
        }

        std::vector<int> caloRecHit_ho_indices = FindIndices("caloRecHit_ho_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
        for(auto j: caloRecHit_ho_indices ){
          out.caloRecHit_ho_rho.push_back(tau.caloRecHit_ho_rho.at(j));
          out.caloRecHit_ho_eta.push_back(tau.caloRecHit_ho_eta.at(j));
          out.caloRecHit_ho_phi.push_back(tau.caloRecHit_ho_phi.at(j));
          out.caloRecHit_ho_energy.push_back(tau.caloRecHit_ho_energy.at(j));
          out.caloRecHit_ho_time.push_back(tau.caloRecHit_ho_time.at(j));
          out.caloRecHit_ho_detId.push_back(tau.caloRecHit_ho_detId.at(j));
          out.caloRecHit_ho_aux.push_back(tau.caloRecHit_ho_aux.at(j));
          out.caloRecHit_ho_flags.push_back(tau.caloRecHit_ho_flags.at(j));
        }
        std::vector<int> caloRecHit_hf_indices = FindIndices("caloRecHit_hf_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
        for(auto j: caloRecHit_hf_indices ){
          out.caloRecHit_hf_rho.push_back(tau.caloRecHit_hf_rho.at(j))
          out.caloRecHit_hf_eta.push_back(tau.caloRecHit_hf_eta.at(j))
          out.caloRecHit_hf_phi.push_back(tau.caloRecHit_hf_phi.at(j))
          out.caloRecHit_hf_energy.push_back(tau.caloRecHit_hf_energy.at(j))
          out.caloRecHit_hf_time.push_back(tau.caloRecHit_hf_time.at(j))
          out.caloRecHit_hf_detId.push_back(tau.caloRecHit_hf_detId.at(j))
          out.caloRecHit_hf_flags.push_back(tau.caloRecHit_hf_flags.at(j))
          out.caloRecHit_hf_timeFalling.push_back(tau.caloRecHit_hf_timeFalling.at(j))
          out.caloRecHit_hf_auxHF.push_back(tau.caloRecHit_hf_auxHF.at(j))
          out.caloRecHit_hf_aux.push_back(tau.caloRecHit_hf_aux.at(j))
        }

        std::vector<int> patatrack_indices = FindIndices("patatrack_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
        for(auto j: patatrack_indices ){
          out.patatrack_pt.push_back( tau.patatrack_pt.at(j));
          out.patatrack_eta.push_back( tau.patatrack_eta.at(j));
          out.patatrack_phi.push_back( tau.patatrack_phi.at(j));
          out.patatrack_chi2.push_back( tau.patatrack_chi2.at(j));
          out.patatrack_ndof.push_back( tau.patatrack_ndof.at(j));
          out.patatrack_charge.push_back( tau.patatrack_charge.at(j));
          out.patatrack_quality.push_back( tau.patatrack_quality.at(j));
          out.patatrack_dxy.push_back( tau.patatrack_dxy.at(j));
          out.patatrack_dz.push_back( tau.patatrack_dz.at(j));
          out.patatrack_vertex_id.push_back( tau.patatrack_vertex_id.at(j));
        }
        std::vector<int> patavert_indices = FindIndices("patavert_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
        for(auto j: patavert_indices ){
          out.patavert_z.push_back(tau.patavert_z.at(j));
          out.patavert_weight.push_back(tau.patavert_weight.at(j));
          out.patavert_ptv2.push_back(tau.patavert_ptv2.at(j));
          out.patavert_chi2.push_back(tau.patavert_chi2.at(j));
          out.patavert_ndof.push_back(tau.patavert_ndof.at(j));
        }

      } // end of loop over l1taus to fill
    }
    */
    void FillTauBranches(const Tau& tau, bool isQCD)
    {

          auto& out = trainingTauTuple();
          float l1pt_threshold= 20.0;
          float genlep_eta_threshold = 2.1 ;
          float genlep_VisPt_threshold = 15.0 ;
          float delta_R_threshold = 0.5;
          std::vector<int> taus_indices;
          std::cout << "è entrato " << std::endl;
          //const auto fill_branches_l1= [&](std::vector<int> indices){
          //  FillVar(tau, out, indices);
          //};


         if(isQCD){
            for(std::vector<int>::size_type k = 0; k != tau.l1Tau_pt.size(); k++){
                if (tau.l1Tau_pt.at(k) > l1pt_threshold && tau.l1Tau_hwIso.at(k) > 0 ) {
                  taus_indices.push_back(k);
                  out.genLepton_kind = 6;
                  out.genLepton_vis_pt = -100.;
                  out.genLepton_vis_eta = -100.;
                  out.genLepton_vis_phi = -100.;
                  out.genLepton_charge = -100;
                  out.genLepton_vis_mass = -100.;
                }
            }

            for (auto& k : taus_indices){
              CP_BR(evt);
              CP_BR(run);
              CP_BR(lumi);
              CP_BR(genEventWeight);
              CP_BR(sampleType);

              out.l1Tau_pt = tau.l1Tau_pt.at(k);
              out.l1Tau_eta = tau.l1Tau_eta.at(k);
              out.l1Tau_phi = tau.l1Tau_phi.at(k);
              out.l1Tau_mass = tau.l1Tau_mass.at(k);
              out.l1Tau_hwIso = tau.l1Tau_hwIso.at(k);
              out.l1Tau_hwQual = tau.l1Tau_hwQual.at(k);
              out.l1Tau_towerIEta = tau.l1Tau_towerIEta.at(k);
              out.l1Tau_towerIPhi = tau.l1Tau_towerIPhi.at(k);
              out.l1Tau_rawEt = tau.l1Tau_rawEt.at(k);
              out.l1Tau_isoEt = tau.l1Tau_isoEt.at(k);
              out.l1Tau_hasEM = tau.l1Tau_hasEM.at(k);
              out.l1Tau_isMerged = tau.l1Tau_isMerged.at(k);

              //Float_t phi1,Float_t eta1,Float_t phi2,Float_t eta2

              std::vector<int> caloRecHit_ee_indices = FindIndices("caloRecHit_ee_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
              for(auto j: caloRecHit_ee_indices ){
                out.caloRecHit_ee_rho.push_back(tau.caloRecHit_ee_rho.at(j));
                out.caloRecHit_ee_eta.push_back(tau.caloRecHit_ee_eta.at(j));
                out.caloRecHit_ee_phi.push_back(tau.caloRecHit_ee_phi.at(j));
                out.caloRecHit_ee_energy.push_back(tau.caloRecHit_ee_energy.at(j));
                out.caloRecHit_ee_time.push_back(tau.caloRecHit_ee_time.at(j));
                out.caloRecHit_ee_detId.push_back(tau.caloRecHit_ee_detId.at(j));
                out.caloRecHit_ee_chi2.push_back(tau.caloRecHit_ee_chi2.at(j));
                out.caloRecHit_ee_energyError.push_back(tau.caloRecHit_ee_energyError.at(j));
                out.caloRecHit_ee_timeError.push_back(tau.caloRecHit_ee_timeError.at(j));
                out.caloRecHit_ee_flagsBits.push_back(tau.caloRecHit_ee_flagsBits.at(j));
                out.caloRecHit_ee_isRecovered.push_back(tau.caloRecHit_ee_isRecovered.at(j));
                out.caloRecHit_ee_isTimeValid.push_back(tau.caloRecHit_ee_isTimeValid.at(j));
                out.caloRecHit_ee_isTimeErrorValid.push_back(tau.caloRecHit_ee_isTimeErrorValid.at(j));
              }

              std::vector<int> caloRecHit_eb_indices = FindIndices("caloRecHit_eb_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
              for(auto j: caloRecHit_eb_indices ){
                out.caloRecHit_eb_rho.push_back(tau.caloRecHit_eb_rho.at(j));
                out.caloRecHit_eb_eta.push_back(tau.caloRecHit_eb_eta.at(j));
                out.caloRecHit_eb_phi.push_back(tau.caloRecHit_eb_phi.at(j));
                out.caloRecHit_eb_energy.push_back(tau.caloRecHit_eb_energy.at(j));
                out.caloRecHit_eb_time.push_back(tau.caloRecHit_eb_time.at(j));
                out.caloRecHit_eb_detId.push_back(tau.caloRecHit_eb_detId.at(j));
                out.caloRecHit_eb_chi2.push_back(tau.caloRecHit_eb_chi2.at(j));
                out.caloRecHit_eb_energyError.push_back(tau.caloRecHit_eb_energyError.at(j));
                out.caloRecHit_eb_timeError.push_back(tau.caloRecHit_eb_timeError.at(j));
                out.caloRecHit_eb_flagsBits.push_back(tau.caloRecHit_eb_flagsBits.at(j));
                out.caloRecHit_eb_isRecovered.push_back(tau.caloRecHit_eb_isRecovered.at(j));
                out.caloRecHit_eb_isTimeValid.push_back(tau.caloRecHit_eb_isTimeValid.at(j));
                out.caloRecHit_eb_isTimeErrorValid.push_back(tau.caloRecHit_eb_isTimeErrorValid.at(j));
              }

              std::vector<int> caloRecHit_hbhe_indices = FindIndices("caloRecHit_hbhe_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
              for(auto j: caloRecHit_hbhe_indices ){
                out.caloRecHit_hbhe_rho.push_back(tau.caloRecHit_hbhe_rho.at(j));
                out.caloRecHit_hbhe_eta.push_back(tau.caloRecHit_hbhe_eta.at(j));
                out.caloRecHit_hbhe_phi.push_back(tau.caloRecHit_hbhe_phi.at(j));
                out.caloRecHit_hbhe_energy.push_back(tau.caloRecHit_hbhe_energy.at(j));
                out.caloRecHit_hbhe_time.push_back(tau.caloRecHit_hbhe_time.at(j));
                out.caloRecHit_hbhe_detId.push_back(tau.caloRecHit_hbhe_detId.at(j));
                out.caloRecHit_hbhe_chi2.push_back(tau.caloRecHit_hbhe_chi2.at(j));
                out.caloRecHit_hbhe_flags.push_back(tau.caloRecHit_hbhe_flags.at(j));
                out.caloRecHit_hbhe_eraw.push_back(tau.caloRecHit_hbhe_eraw.at(j));
                out.caloRecHit_hbhe_eaux.push_back(tau.caloRecHit_hbhe_eaux.at(j));
                out.caloRecHit_hbhe_timeFalling.push_back(tau.caloRecHit_hbhe_timeFalling.at(j));
                out.caloRecHit_hbhe_idFront.push_back(tau.caloRecHit_hbhe_idFront.at(j));
                out.caloRecHit_hbhe_rho_front.push_back(tau.caloRecHit_hbhe_rho_front.at(j));
                out.caloRecHit_hbhe_eta_front.push_back(tau.caloRecHit_hbhe_eta_front.at(j));
                out.caloRecHit_hbhe_phi_front.push_back(tau.caloRecHit_hbhe_phi_front.at(j));
                out.caloRecHit_hbhe_auxHBHE.push_back(tau.caloRecHit_hbhe_auxHBHE.at(j));
                out.caloRecHit_hbhe_auxPhase1.push_back(tau.caloRecHit_hbhe_auxPhase1.at(j));
                out.caloRecHit_hbhe_auxTDC.push_back(tau.caloRecHit_hbhe_auxTDC.at(j));
                out.caloRecHit_hbhe_isMerged.push_back(tau.caloRecHit_hbhe_isMerged.at(j));
              }

              std::vector<int> caloRecHit_ho_indices = FindIndices("caloRecHit_ho_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
              for(auto j: caloRecHit_ho_indices ){
                out.caloRecHit_ho_rho.push_back(tau.caloRecHit_ho_rho.at(j));
                out.caloRecHit_ho_eta.push_back(tau.caloRecHit_ho_eta.at(j));
                out.caloRecHit_ho_phi.push_back(tau.caloRecHit_ho_phi.at(j));
                out.caloRecHit_ho_energy.push_back(tau.caloRecHit_ho_energy.at(j));
                out.caloRecHit_ho_time.push_back(tau.caloRecHit_ho_time.at(j));
                out.caloRecHit_ho_detId.push_back(tau.caloRecHit_ho_detId.at(j));
                out.caloRecHit_ho_aux.push_back(tau.caloRecHit_ho_aux.at(j));
                out.caloRecHit_ho_flags.push_back(tau.caloRecHit_ho_flags.at(j));
              }
              std::vector<int> caloRecHit_hf_indices = FindIndices("caloRecHit_hf_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
              for(auto j: caloRecHit_hf_indices ){
                out.caloRecHit_hf_rho.push_back(tau.caloRecHit_hf_rho.at(j));
                out.caloRecHit_hf_eta.push_back(tau.caloRecHit_hf_eta.at(j));
                out.caloRecHit_hf_phi.push_back(tau.caloRecHit_hf_phi.at(j));
                out.caloRecHit_hf_energy.push_back(tau.caloRecHit_hf_energy.at(j));
                out.caloRecHit_hf_time.push_back(tau.caloRecHit_hf_time.at(j));
                out.caloRecHit_hf_detId.push_back(tau.caloRecHit_hf_detId.at(j));
                out.caloRecHit_hf_flags.push_back(tau.caloRecHit_hf_flags.at(j));
                out.caloRecHit_hf_timeFalling.push_back(tau.caloRecHit_hf_timeFalling.at(j));
                out.caloRecHit_hf_auxHF.push_back(tau.caloRecHit_hf_auxHF.at(j));
                out.caloRecHit_hf_aux.push_back(tau.caloRecHit_hf_aux.at(j));
              }

              std::vector<int> patatrack_indices = FindIndices("patatrack_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
              for(auto j: patatrack_indices ){
                out.patatrack_pt.push_back( tau.patatrack_pt.at(j));
                out.patatrack_eta.push_back( tau.patatrack_eta.at(j));
                out.patatrack_phi.push_back( tau.patatrack_phi.at(j));
                out.patatrack_chi2.push_back( tau.patatrack_chi2.at(j));
                out.patatrack_ndof.push_back( tau.patatrack_ndof.at(j));
                out.patatrack_charge.push_back( tau.patatrack_charge.at(j));
                out.patatrack_quality.push_back( tau.patatrack_quality.at(j));
                out.patatrack_dxy.push_back( tau.patatrack_dxy.at(j));
                out.patatrack_dz.push_back( tau.patatrack_dz.at(j));
                out.patatrack_vertex_id.push_back( tau.patatrack_vertex_id.at(j));
              }
              std::vector<int> patavert_indices = FindIndices("patavert_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
              for(auto j: patavert_indices ){
                out.patavert_z.push_back(tau.patavert_z.at(j));
                out.patavert_weight.push_back(tau.patavert_weight.at(j));
                out.patavert_ptv2.push_back(tau.patavert_ptv2.at(j));
                out.patavert_chi2.push_back(tau.patavert_chi2.at(j));
                out.patavert_ndof.push_back(tau.patavert_ndof.at(j));
              }

            } // end of loop over l1taus to fill
          }
          else{
            std::cout << "NON QCD  " << std::endl;
            for (std::vector<int>::size_type i = 0; i != tau.genLepton_kind.size(); i++){
                /* 0. define quantities we need */
                float previous_deltaR = 100.0;
                int l1tau_previous_index = 100;
                std::vector<int> current_tau_index;
                //std::cout << "il gentau... " << i << std::endl;
                //reco_tau::gen_truth::GenLepton::Kind::TauDecayedToHadrons = 5
                /* 1. ask genlepton to be an hadronic tau within tracker acceptance and with a visible pt threshold */
                if(tau.genLepton_kind.at(i) != 5 || tau.genLepton_vis_pt.at(i) < genlep_VisPt_threshold || std::abs(tau.genLepton_vis_eta.at(i))>genlep_eta_threshold ) continue;
                //std::cout << "analizzo il gentau... " << i << std::endl;

                /* 1-1. fill general quantities */
                out.genLepton_kind = tau.genLepton_kind.at(i);
                out.genLepton_vis_pt = tau.genLepton_vis_pt.at(i);
                out.genLepton_vis_eta = tau.genLepton_vis_eta.at(i);
                out.genLepton_vis_phi = tau.genLepton_vis_phi.at(i);
                out.genLepton_charge = tau.genLepton_charge.at(i);
                out.genLepton_vis_mass = tau.genLepton_vis_mass.at(i);

                /* 2. loop over taus to find the closest one */
                for (std::vector<int>::size_type j = 0; j != tau.l1Tau_pt.size(); j++){
                  /* 2-1 isolation requirement and pt threshold on l1tau */
                  //std::cout << "il tau... " << j << std::endl;
                  if (tau.l1Tau_pt.at(j) < l1pt_threshold || tau.l1Tau_hwIso.at(j) < 0 ) continue;

                  //std::cout << "analizzo il tau... " << j << std::endl;
                  /* 2-2  discharge previously assigned l1taus */
                  if (std::find(taus_indices.begin(),taus_indices.end(), j) != taus_indices.end()) continue;

                  /* let's calculate deltaR: if it's the minimum */
                  float current_deltaR = deltaR(tau.l1Tau_phi.at(j), tau.l1Tau_eta.at(j), tau.genLepton_vis_phi.at(i), tau.genLepton_vis_eta.at(i));
                  //std::cout << "tau.l1Tau_phi " << tau.l1Tau_phi.at(j)  << std::endl;
                  //std::cout << "tau.l1Tau_eta " << tau.l1Tau_eta.at(j)  << std::endl;
                  //std::cout << "tau.genLepton_vis_phi " << tau.genLepton_vis_phi.at(i)  << std::endl;
                  //std::cout << "tau.genLepton_vis_eta " << tau.genLepton_vis_eta.at(i) << std::endl;
                  if (current_deltaR <0.3 && current_deltaR < previous_deltaR ){
                    previous_deltaR = current_deltaR;
                    l1tau_previous_index = j;
                  }

                } // end of loop over l1taus to find the closest one wrt genlepton

                //std::cout << "best tau ? " << l1tau_previous_index << "\n" << std::endl;
                /* save index */
                if(l1tau_previous_index!=100) taus_indices.push_back(l1tau_previous_index);
                if(l1tau_previous_index!=100) current_tau_index.push_back(l1tau_previous_index);
                //std::cout << "size di current tau index == " << current_tau_index.size() << std::endl;

                for (auto& k : current_tau_index){

                  //std::cout << "è entrato su tauindex loop e sta sull'elemento: " << k << std::endl;
                  CP_BR(evt);
                  CP_BR(run);
                  CP_BR(lumi);
                  CP_BR(genEventWeight);
                  CP_BR(sampleType);

                  out.l1Tau_pt = tau.l1Tau_pt.at(k);
                  out.l1Tau_eta = tau.l1Tau_eta.at(k);
                  out.l1Tau_phi = tau.l1Tau_phi.at(k);
                  out.l1Tau_mass = tau.l1Tau_mass.at(k);
                  out.l1Tau_hwIso = tau.l1Tau_hwIso.at(k);
                  out.l1Tau_hwQual = tau.l1Tau_hwQual.at(k);
                  out.l1Tau_towerIEta = tau.l1Tau_towerIEta.at(k);
                  out.l1Tau_towerIPhi = tau.l1Tau_towerIPhi.at(k);
                  out.l1Tau_rawEt = tau.l1Tau_rawEt.at(k);
                  out.l1Tau_isoEt = tau.l1Tau_isoEt.at(k);
                  out.l1Tau_hasEM = tau.l1Tau_hasEM.at(k);
                  out.l1Tau_isMerged = tau.l1Tau_isMerged.at(k);

                  //Float_t phi1,Float_t eta1,Float_t phi2,Float_t eta2

                  std::vector<int> caloRecHit_ee_indices = FindIndices("caloRecHit_ee_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
                  if (caloRecHit_ee_indices.size()>0){
                    for(auto j: caloRecHit_ee_indices ){
                      //std::cout << "è entrato su caloRecHit_ee_index loop e sta sull'elemento: " << j << std::endl;
                      out.caloRecHit_ee_rho.push_back(tau.caloRecHit_ee_rho.at(j));
                      out.caloRecHit_ee_eta.push_back(tau.caloRecHit_ee_eta.at(j));
                      out.caloRecHit_ee_phi.push_back(tau.caloRecHit_ee_phi.at(j));
                      out.caloRecHit_ee_energy.push_back(tau.caloRecHit_ee_energy.at(j));
                      out.caloRecHit_ee_time.push_back(tau.caloRecHit_ee_time.at(j));
                      out.caloRecHit_ee_detId.push_back(tau.caloRecHit_ee_detId.at(j));
                      out.caloRecHit_ee_chi2.push_back(tau.caloRecHit_ee_chi2.at(j));
                      out.caloRecHit_ee_energyError.push_back(tau.caloRecHit_ee_energyError.at(j));
                      out.caloRecHit_ee_timeError.push_back(tau.caloRecHit_ee_timeError.at(j));
                      out.caloRecHit_ee_flagsBits.push_back(tau.caloRecHit_ee_flagsBits.at(j));
                      out.caloRecHit_ee_isRecovered.push_back(tau.caloRecHit_ee_isRecovered.at(j));
                      out.caloRecHit_ee_isTimeValid.push_back(tau.caloRecHit_ee_isTimeValid.at(j));
                      out.caloRecHit_ee_isTimeErrorValid.push_back(tau.caloRecHit_ee_isTimeErrorValid.at(j));
                    }
                  }


                  std::vector<int> caloRecHit_eb_indices = FindIndices("caloRecHit_eb_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
                  if (caloRecHit_eb_indices.size()>0){
                    for(auto j: caloRecHit_eb_indices ){
                      out.caloRecHit_eb_rho.push_back(tau.caloRecHit_eb_rho.at(j));
                      out.caloRecHit_eb_eta.push_back(tau.caloRecHit_eb_eta.at(j));
                      out.caloRecHit_eb_phi.push_back(tau.caloRecHit_eb_phi.at(j));
                      out.caloRecHit_eb_energy.push_back(tau.caloRecHit_eb_energy.at(j));
                      out.caloRecHit_eb_time.push_back(tau.caloRecHit_eb_time.at(j));
                      out.caloRecHit_eb_detId.push_back(tau.caloRecHit_eb_detId.at(j));
                      out.caloRecHit_eb_chi2.push_back(tau.caloRecHit_eb_chi2.at(j));
                      out.caloRecHit_eb_energyError.push_back(tau.caloRecHit_eb_energyError.at(j));
                      out.caloRecHit_eb_timeError.push_back(tau.caloRecHit_eb_timeError.at(j));
                      out.caloRecHit_eb_flagsBits.push_back(tau.caloRecHit_eb_flagsBits.at(j));
                      out.caloRecHit_eb_isRecovered.push_back(tau.caloRecHit_eb_isRecovered.at(j));
                      out.caloRecHit_eb_isTimeValid.push_back(tau.caloRecHit_eb_isTimeValid.at(j));
                      out.caloRecHit_eb_isTimeErrorValid.push_back(tau.caloRecHit_eb_isTimeErrorValid.at(j));
                    }
                  }


                  std::vector<int> caloRecHit_hbhe_indices = FindIndices("caloRecHit_hbhe_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
                  // std::cout << "calorechit_hbhe ha size " << tau.caloRecHit_hbhe_rho.size()<< std::endl;
                  // std::cout << "calorechit_hbhe index ha size " << caloRecHit_hbhe_indices.size()<< std::endl;
                  if (caloRecHit_hbhe_indices.size()>0){
                    for(auto j: caloRecHit_hbhe_indices ){
                      // std::cout << "calorechit_ index " << j<< std::endl;
                      out.caloRecHit_hbhe_rho.push_back(tau.caloRecHit_hbhe_rho.at(j));
                      out.caloRecHit_hbhe_eta.push_back(tau.caloRecHit_hbhe_eta.at(j));
                      out.caloRecHit_hbhe_phi.push_back(tau.caloRecHit_hbhe_phi.at(j));
                      out.caloRecHit_hbhe_energy.push_back(tau.caloRecHit_hbhe_energy.at(j));
                      out.caloRecHit_hbhe_time.push_back(tau.caloRecHit_hbhe_time.at(j));
                      out.caloRecHit_hbhe_detId.push_back(tau.caloRecHit_hbhe_detId.at(j));
                      out.caloRecHit_hbhe_chi2.push_back(tau.caloRecHit_hbhe_chi2.at(j));
                      out.caloRecHit_hbhe_flags.push_back(tau.caloRecHit_hbhe_flags.at(j));
                      // out.caloRecHit_hbhe_eraw.push_back(tau.caloRecHit_hbhe_eraw.at(j));
                      // out.caloRecHit_hbhe_eaux.push_back(tau.caloRecHit_hbhe_eaux.at(j));
                      // out.caloRecHit_hbhe_timeFalling.push_back(tau.caloRecHit_hbhe_timeFalling.at(j));
                      // out.caloRecHit_hbhe_idFront.push_back(tau.caloRecHit_hbhe_idFront.at(j));
                      // out.caloRecHit_hbhe_rho_front.push_back(tau.caloRecHit_hbhe_rho_front.at(j));
                      // out.caloRecHit_hbhe_eta_front.push_back(tau.caloRecHit_hbhe_eta_front.at(j));
                      // out.caloRecHit_hbhe_phi_front.push_back(tau.caloRecHit_hbhe_phi_front.at(j));
                      // out.caloRecHit_hbhe_auxHBHE.push_back(tau.caloRecHit_hbhe_auxHBHE.at(j));/*
                      // out.caloRecHit_hbhe_auxPhase1.push_back(tau.caloRecHit_hbhe_auxPhase1.at(j));
                      // out.caloRecHit_hbhe_auxTDC.push_back(tau.caloRecHit_hbhe_auxTDC.at(j));
                      // out.caloRecHit_hbhe_isMerged.push_back(tau.caloRecHit_hbhe_isMerged.at(j));*/
                    }
                  }
                  std::vector<int> caloRecHit_ho_indices = FindIndices("caloRecHit_ho_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
                  if (caloRecHit_ho_indices.size()>0){
                    for(auto j: caloRecHit_ho_indices ){
                      out.caloRecHit_ho_rho.push_back(tau.caloRecHit_ho_rho.at(j));
                      out.caloRecHit_ho_eta.push_back(tau.caloRecHit_ho_eta.at(j));
                      out.caloRecHit_ho_phi.push_back(tau.caloRecHit_ho_phi.at(j));
                      out.caloRecHit_ho_energy.push_back(tau.caloRecHit_ho_energy.at(j));
                      out.caloRecHit_ho_time.push_back(tau.caloRecHit_ho_time.at(j));
                      out.caloRecHit_ho_detId.push_back(tau.caloRecHit_ho_detId.at(j));
                      out.caloRecHit_ho_aux.push_back(tau.caloRecHit_ho_aux.at(j));
                      out.caloRecHit_ho_flags.push_back(tau.caloRecHit_ho_flags.at(j));
                    }
                  }

                  std::vector<int> caloRecHit_hf_indices = FindIndices("caloRecHit_hf_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
                  if (caloRecHit_hf_indices.size()>0){
                    for(auto j: caloRecHit_hf_indices ){
                      out.caloRecHit_hf_rho.push_back(tau.caloRecHit_hf_rho.at(j));
                      out.caloRecHit_hf_eta.push_back(tau.caloRecHit_hf_eta.at(j));
                      out.caloRecHit_hf_phi.push_back(tau.caloRecHit_hf_phi.at(j));
                      out.caloRecHit_hf_energy.push_back(tau.caloRecHit_hf_energy.at(j));
                      out.caloRecHit_hf_time.push_back(tau.caloRecHit_hf_time.at(j));
                      out.caloRecHit_hf_detId.push_back(tau.caloRecHit_hf_detId.at(j));
                      out.caloRecHit_hf_flags.push_back(tau.caloRecHit_hf_flags.at(j));
                      out.caloRecHit_hf_timeFalling.push_back(tau.caloRecHit_hf_timeFalling.at(j));
                      out.caloRecHit_hf_auxHF.push_back(tau.caloRecHit_hf_auxHF.at(j));
                      out.caloRecHit_hf_aux.push_back(tau.caloRecHit_hf_aux.at(j));
                    }
                  }

                  std::vector<int> patatrack_indices = FindIndices("patatrack_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
                  if (patatrack_indices.size()>0){
                    for(auto j: patatrack_indices ){
                      out.patatrack_pt.push_back( tau.patatrack_pt.at(j));
                      out.patatrack_eta.push_back( tau.patatrack_eta.at(j));
                      out.patatrack_phi.push_back( tau.patatrack_phi.at(j));
                      out.patatrack_chi2.push_back( tau.patatrack_chi2.at(j));
                      out.patatrack_ndof.push_back( tau.patatrack_ndof.at(j));
                      out.patatrack_charge.push_back( tau.patatrack_charge.at(j));
                      out.patatrack_quality.push_back( tau.patatrack_quality.at(j));
                      out.patatrack_dxy.push_back( tau.patatrack_dxy.at(j));
                      out.patatrack_dz.push_back( tau.patatrack_dz.at(j));
                      out.patatrack_vertex_id.push_back( tau.patatrack_vertex_id.at(j));
                    }
                  }
                  /*
                  std::vector<int> patavert_indices = FindIndices("patavert_",tau.l1Tau_eta.at(k), tau.l1Tau_phi.at(k) , delta_R_threshold);
                  if (patavert_indices.size()>0){
                    for(auto j: patavert_indices ){
                      out.patavert_z.push_back(tau.patavert_z.at(j));
                      out.patavert_weight.push_back(tau.patavert_weight.at(j));
                      out.patavert_ptv2.push_back(tau.patavert_ptv2.at(j));
                      out.patavert_chi2.push_back(tau.patavert_chi2.at(j));
                      out.patavert_ndof.push_back(tau.patavert_ndof.at(j));
                    }
                  }
                  */

                } // end of loop over l1taus to fill




            } // end of loop over genleptons
          } // end of else

    }
    #undef CP_BR

private:
    const Arguments args;
    std::shared_ptr<TFile> inputFile, outputFile;
    L2TauTuple l2tauTuple;
    TrainingTauTuple trainingTauTuple;

};

} // namespace analysis

PROGRAM_MAIN(analysis::L2TrainingTupleProducer, Arguments)
