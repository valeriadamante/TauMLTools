/*! Produce training tuple from tau tuple.
*/

#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/variadic.hpp>
#include <boost/math/constants/constants.hpp>
#include <ROOT/RDataFrame.hxx>
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
    //run::Argument<float> training_weight_factor{"training-weight-factor",
    //    "additional factor to the normalization of the training weights", 4.f};
    //run::Argument<int> parity{"parity", "take odd (parity=1), even (parity=0) or all (parity=-1) events", -1};
    run::Argument<int> isQCDDataVBF{"isQCDDataVBF", "0 = QCD; 1 = TT,DY,ZPrime; 2 = VBF; 3=Data", false};
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
            //std::cout << "entry = " << current_entry << std::endl;
            l2tauTuple.GetEntry(current_entry);
            const auto& tau = l2tauTuple.data();
            FillTauBranches(tau, args.isQCDDataVBF());
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

    void FillL1Taus(const Tau& tau, std::vector<int> tauindices, float delta_R_threshold){

      auto& out = trainingTauTuple();
      for (auto& k : tauindices){
        //std::cout << "copio il tau .. " << k << std::endl;

        // global quantities
        //std::cout << "evento .. " << tau.evt << std::endl;
        CP_BR(evt);
        CP_BR(run);
        CP_BR(lumi);
        CP_BR(genEventWeight);
        CP_BR(sampleType);
        CP_BR(defaultDiTauPath_lastModuleIndex)
        CP_BR(defaultDiTauPath_result)

        // L1 tau quantities

        out.l1Tau_index = k;
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

        // calo rec hit quantities in a specific Delta R cone

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
        if (caloRecHit_hbhe_indices.size()>0){
          for(auto j: caloRecHit_hbhe_indices ){
            out.caloRecHit_hbhe_rho.push_back(tau.caloRecHit_hbhe_rho.at(j));
            out.caloRecHit_hbhe_eta.push_back(tau.caloRecHit_hbhe_eta.at(j));
            out.caloRecHit_hbhe_phi.push_back(tau.caloRecHit_hbhe_phi.at(j));
            out.caloRecHit_hbhe_energy.push_back(tau.caloRecHit_hbhe_energy.at(j));
            out.caloRecHit_hbhe_time.push_back(tau.caloRecHit_hbhe_time.at(j));
            out.caloRecHit_hbhe_detId.push_back(tau.caloRecHit_hbhe_detId.at(j));
            out.caloRecHit_hbhe_chi2.push_back(tau.caloRecHit_hbhe_chi2.at(j));
            out.caloRecHit_hbhe_flags.push_back(tau.caloRecHit_hbhe_flags.at(j));
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

        // patatrack and patavertices

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

        for(std::vector<int>::size_type j = 0; j != tau.patavert_z.size(); j++){
          out.patavert_z.push_back(tau.patavert_z.at(j));
          out.patavert_weight.push_back(tau.patavert_weight.at(j));
          out.patavert_ptv2.push_back(tau.patavert_ptv2.at(j));
          out.patavert_chi2.push_back(tau.patavert_chi2.at(j));
          out.patavert_ndof.push_back(tau.patavert_ndof.at(j));
        }
        trainingTauTuple.Fill();
      } // end of loop over l1taus to fill
    }

    void FillTauBranches(const Tau& tau, int isQCDDataVBF)
    {

          auto& out = trainingTauTuple();
          float l1pt_threshold= 32.0;
          float genlep_eta_threshold = 2.1 ;
          float genlep_VisPt_threshold = 15.0 ;
          float delta_R_threshold = 0.5;
          std::vector<int> taus_indices;
          const auto Fill_L1_Taus= [&](std::vector<int> tauindices){
            FillL1Taus(tau, tauindices, delta_R_threshold);
          };

         // QCD
         if(isQCDDataVBF==0){
            for(std::vector<int>::size_type k = 0; k != tau.l1Tau_pt.size(); k++){
                if( tau.l1Tau_pt.at(k) >= l1pt_threshold  && (tau.l1Tau_hwIso.at(k) > 0 || tau.l1Tau_pt.at(k) >= 70)) {
                  taus_indices.push_back(k);
                  out.genLepton_kind = 6;
                  out.genLepton_vis_pt = -100.;
                  out.genLepton_vis_eta = -100.;
                  out.genLepton_vis_phi = -100.;
                  out.genLepton_charge = -100;
                  out.genLepton_vis_mass = -100.;
                }
            }
            Fill_L1_Taus(taus_indices);
          }
         // TT, DY and ZPrime
         else if(isQCDDataVBF==1){
            //std::cout<< "evento .. " << tau.evt << std::endl ;
            for (std::vector<int>::size_type i = 0; i != tau.genLepton_kind.size(); i++){
                /* 0. define useful quantities */
                float previous_deltaR = 100.0;
                int l1tau_previous_index = 100;
                std::vector<int> current_tau_index;

                /* 1. ask genlepton to be an hadronic tau within tracker acceptance and with a visible pt threshold */
                if(tau.genLepton_kind.at(i) != 5 || tau.genLepton_vis_pt.at(i) < genlep_VisPt_threshold || std::abs(tau.genLepton_vis_eta.at(i))>genlep_eta_threshold ) continue;
                //std::cout << "analizzo il gentau... " << i << std::endl;

                /* 1-1. fill general quantities */
                out.genLepton_index = i;
                out.genLepton_kind = tau.genLepton_kind.at(i);
                out.genLepton_vis_pt = tau.genLepton_vis_pt.at(i);
                out.genLepton_vis_eta = tau.genLepton_vis_eta.at(i);
                out.genLepton_vis_phi = tau.genLepton_vis_phi.at(i);
                out.genLepton_charge = tau.genLepton_charge.at(i);
                out.genLepton_vis_mass = tau.genLepton_vis_mass.at(i);

                /* 2. loop over taus to find the closest one */
                for (std::vector<int>::size_type j = 0; j != tau.l1Tau_pt.size(); j++){
                    /* 2-1.  discharge previously assigned l1taus */
                    if (std::find(taus_indices.begin(),taus_indices.end(), j) != taus_indices.end()) {
                        //std::cout << " do not consider .. " << j << std::endl;
                        continue;
                    }
                    /* 2-2. isolation requirement and pt threshold on l1tau */
                    if (tau.l1Tau_pt.at(j) >= l1pt_threshold  && (tau.l1Tau_hwIso.at(j) > 0 || tau.l1Tau_pt.at(j) >= 70)) {
                        //std::cout << "analizzo il tau ... " << j << std::endl;

                        /* 2-3. let's calculate deltaR: if it's the minimum */
                        float current_deltaR = deltaR(tau.l1Tau_phi.at(j), tau.l1Tau_eta.at(j), tau.genLepton_vis_phi.at(i), tau.genLepton_vis_eta.at(i));
                        if (current_deltaR <0.3 && current_deltaR < previous_deltaR ){
                          previous_deltaR = current_deltaR;
                          l1tau_previous_index = j;
                        }
                    }
                } // end of loop over l1taus to find the closest one wrt genlepton

                /* 3. save index */
                if(l1tau_previous_index!=100) taus_indices.push_back(l1tau_previous_index);
                if(l1tau_previous_index!=100) current_tau_index.push_back(l1tau_previous_index);
                //std::cout << "closest tau ... " << l1tau_previous_index << std::endl;
                //for (auto& p : taus_indices){
                  //std::cout << "next iteration do not consider tau ... " << p << std::endl;
                //}
                Fill_L1_Taus(current_tau_index);
            } // end of loop over genleptons
            //std::cout << std::endl;
         } // end of else
         // VBF
         else if(isQCDDataVBF==2) {
           // first let's check event passing double big or
              if(tau.defaultDiTauPath_lastModuleIndex>5){
                  for (std::vector<int>::size_type i = 0; i != tau.genLepton_kind.size(); i++){
                      /* 0. define useful quantities */
                      float previous_deltaR = 100.0;
                      int l1tau_previous_index = 100;
                      std::vector<int> current_tau_index;
                      /* 1. ask genlepton to be an hadronic tau within tracker acceptance and with a visible pt threshold */
                      if(tau.genLepton_kind.at(i) != 5 || tau.genLepton_vis_pt.at(i) < genlep_VisPt_threshold || std::abs(tau.genLepton_vis_eta.at(i))>genlep_eta_threshold ) continue;
                      /* 1-1. fill general quantities */
                      out.genLepton_index = i;
                      out.genLepton_kind = tau.genLepton_kind.at(i);
                      out.genLepton_vis_pt = tau.genLepton_vis_pt.at(i);
                      out.genLepton_vis_eta = tau.genLepton_vis_eta.at(i);
                      out.genLepton_vis_phi = tau.genLepton_vis_phi.at(i);
                      out.genLepton_charge = tau.genLepton_charge.at(i);
                      out.genLepton_vis_mass = tau.genLepton_vis_mass.at(i);
                      /* 2. loop over taus to find the closest one */
                      for (std::vector<int>::size_type j = 0; j != tau.l1Tau_pt.size(); j++){
                          /* 2-1.  discharge previously assigned l1taus */
                          if (std::find(taus_indices.begin(),taus_indices.end(), j) != taus_indices.end()) {
                              continue;
                          }

                          /* 2-2. isolation requirement and pt threshold on l1tau */
                          if (tau.l1Tau_pt.at(j) >= l1pt_threshold  && (tau.l1Tau_hwIso.at(j) > 0 || tau.l1Tau_pt.at(j) >= 70)) {
                              /* 2-3. let's calculate deltaR: if it's the minimum */
                              float current_deltaR = deltaR(tau.l1Tau_phi.at(j), tau.l1Tau_eta.at(j), tau.genLepton_vis_phi.at(i), tau.genLepton_vis_eta.at(i));
                              if (current_deltaR <0.3 && current_deltaR < previous_deltaR ){
                                previous_deltaR = current_deltaR;
                                l1tau_previous_index = j;
                              }
                          }
                      } // end of loop over l1taus to find the closest one wrt genlepton
                      /* 3. save index */
                      if(l1tau_previous_index!=100) taus_indices.push_back(l1tau_previous_index);
                      if(l1tau_previous_index!=100) current_tau_index.push_back(l1tau_previous_index);
                      Fill_L1_Taus(current_tau_index);
                  } // end of loop over genleptons

              } // end of condition of ditauPath to have last module index > hltL1sDoubleTauBigOR (== 5 -> check on gdoc or by running getPathNames.py !!)
          } // end of else
          // DATA
          else if(isQCDDataVBF==3) {
               if(tau.defaultDiTauPath_lastModuleIndex>5){
                 for(std::vector<int>::size_type k = 0; k != tau.l1Tau_pt.size(); k++){
                     if( tau.l1Tau_pt.at(k) >= l1pt_threshold  && (tau.l1Tau_hwIso.at(k) > 0 || tau.l1Tau_pt.at(k) >= 70)) {
                       taus_indices.push_back(k);
                       out.genLepton_kind = 6;
                       out.genLepton_vis_pt = -100.;
                       out.genLepton_vis_eta = -100.;
                       out.genLepton_vis_phi = -100.;
                       out.genLepton_charge = -100;
                       out.genLepton_vis_mass = -100.;
                     }
                 }
                 Fill_L1_Taus(taus_indices);
               }

          }
          // in other cases it should rise an exception!!
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
