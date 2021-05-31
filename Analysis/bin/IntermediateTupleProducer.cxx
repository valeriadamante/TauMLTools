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
#include "TauMLTools/Analysis/interface/L2IntermediateTuple.h"
#include "TauMLTools/Core/interface/ProgressReporter.h"
#include "TauMLTools/Analysis/interface/GenLepton.h"

struct Arguments {
    run::Argument<std::string> input{"input", "input root file with tau tuple"};
    run::Argument<std::string> output{"output", "output root file with training tuple"};
    run::Argument<unsigned> n_threads{"n-threads", "number of threads", 1};
    run::Argument<Long64_t> start_entry{"start-entry", "start entry", 0};
    run::Argument<Long64_t> end_entry{"end-entry", "end entry", std::numeric_limits<Long64_t>::max()};
};

namespace analysis {


class IntermediateTupleProducer {
public:
    using Tau = train_tuple::Tau;
    using L2TauTuple = train_tuple::TrainTuple;
    using TrainingTau = tau_train_tuple::Tau;
    using TrainingTauTuple = tau_train_tuple::TauTrainTuple;
    using IntermTau = interm_tuple::IntermTau;
    using IntermTuple = interm_tuple::TauTrainIntermediateTuple;

    IntermediateTupleProducer(const Arguments& _args) :
        args(_args), inputFile(root_ext::OpenRootFile(args.input())),
        outputFile(root_ext::CreateRootFile(args.output(), ROOT::kLZ4, 4)),
        trainingTauTuple(inputFile.get(), true), intermTuple(outputFile.get(), false)
    {
        if(args.n_threads() > 1)
            ROOT::EnableImplicitMT(args.n_threads());
    }

    void Run()
    { const Long64_t end_entry = std::min(trainingTauTuple.GetEntries(), args.end_entry());
    size_t n_processed = 0, n_total = static_cast<size_t>(end_entry - args.start_entry());
    tools::ProgressReporter reporter(10, std::cout, "Creating training tuple...");
    reporter.SetTotalNumberOfEvents(n_total);
    for(Long64_t current_entry = args.start_entry(); current_entry < end_entry; ++current_entry) {
        //std::cout << "entry = " << current_entry << std::endl;
        trainingTauTuple.GetEntry(current_entry);
        const auto& tau = trainingTauTuple.data();
        FillIntermediateTuple(tau);
        if(++n_processed % 1000 == 0)
            reporter.Report(n_processed);
    }
    reporter.Report(n_processed, true);

    intermTuple.Write();
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
      template<typename Scalar>
      static Scalar DeltaEta(Scalar eta1, Scalar eta2)
      {
        return (eta1-eta2);
      }

      float DeltaR(Float_t phi1,Float_t eta1,Float_t phi2,Float_t eta2) {
        auto dphi = DeltaPhi(phi1, phi2);
        auto deta = DeltaEta(eta1, eta2);
        //std::cout << "deltaR = " << (std::sqrt(deta * deta + dphi * dphi)) <<std::endl;
        return (std::sqrt(deta * deta + dphi * dphi));

       }

      std::vector<int> ReorderByEnergy(std::vector<float>& energy){
        std::vector<int> indices;
        while(indices.size()<energy.size()){
          float current_energy_value=0.;
          int current_index =0;
          for(std::vector<float>::size_type i=0; i < energy.size() ; i++){
            if(std::find(indices.begin(), indices.end(), i)!= indices.end()) continue;
            if(energy.at(i) >= current_energy_value ){
              current_energy_value = energy.at(i);
              current_index = i;
            }
          }
          indices.push_back(current_index);
        }
        return indices;
      }
      #define CP_BR(name) intermTuple().name = tau.name;

      void FillIntermediateTuple(const TrainingTau& tau){

        auto& out = intermTuple();
        // global quantities
        //std::cout << "evento .. " << tau.evt << std::endl;
        CP_BR(evt);
        CP_BR(run);
        CP_BR(lumi);
        CP_BR(l1Tau_index);
        CP_BR(l1Tau_pt);
        CP_BR(l1Tau_eta);
        CP_BR(l1Tau_phi);
        CP_BR(l1Tau_mass);
        CP_BR(l1Tau_hwIso);
        if(tau.genLepton_kind == 5){
          out.genLepton_isTau = true;
        }
        else{
          out.genLepton_isTau = false;
        }
        CP_BR(genLepton_charge);
        CP_BR(genLepton_vis_eta);
        CP_BR(genLepton_vis_pt);
        CP_BR(genLepton_vis_mass);
        CP_BR(genLepton_vis_phi);

        out.nVertices=tau.patavert_z.size();
        std::vector<bool> ecal_isEndCap;
        std::vector<float> ecal_rho;
        std::vector<float> ecal_DeltaEta;
        std::vector<float> ecal_DeltaPhi;
        std::vector<float> ecal_DeltaR;
        std::vector<float> ecal_energy;
        std::vector<float> ecal_chi2;
        std::vector<bool> ecal_isRecovered;

        for(std::vector<float>::size_type i=0; i < tau.caloRecHit_ee_energy.size() ; i++){
          if(tau.caloRecHit_ee_energy.at(i)==0) continue;
          ecal_isEndCap.push_back(true);
          ecal_rho.push_back(tau.caloRecHit_ee_rho.at(i));
          ecal_DeltaEta.push_back(DeltaEta(tau.caloRecHit_ee_eta.at(i), tau.l1Tau_eta));
          ecal_DeltaPhi.push_back(DeltaPhi(tau.caloRecHit_ee_phi.at(i), tau.l1Tau_phi));
          ecal_DeltaR.push_back(DeltaR(tau.caloRecHit_ee_phi.at(i), tau.caloRecHit_ee_eta.at(i), tau.l1Tau_phi, tau.l1Tau_eta));
          ecal_energy.push_back(tau.caloRecHit_ee_energy.at(i));
          ecal_chi2.push_back(tau.caloRecHit_ee_chi2.at(i));
          ecal_isRecovered.push_back(tau.caloRecHit_ee_isRecovered.at(i));
        }
        for(std::vector<float>::size_type i=0; i < tau.caloRecHit_eb_energy.size() ; i++){
          if(tau.caloRecHit_eb_energy.at(i)==0) continue;
          ecal_isEndCap.push_back(false);
          ecal_rho.push_back(tau.caloRecHit_eb_rho.at(i));
          ecal_DeltaEta.push_back(DeltaEta(tau.caloRecHit_eb_eta.at(i), tau.l1Tau_eta));
          ecal_DeltaPhi.push_back(DeltaPhi(tau.caloRecHit_eb_phi.at(i), tau.l1Tau_phi));
          ecal_DeltaR.push_back(DeltaR(tau.caloRecHit_eb_phi.at(i), tau.caloRecHit_eb_eta.at(i), tau.l1Tau_phi, tau.l1Tau_eta));
          ecal_energy.push_back(tau.caloRecHit_eb_energy.at(i));
          ecal_chi2.push_back(tau.caloRecHit_eb_chi2.at(i));
          ecal_isRecovered.push_back(tau.caloRecHit_eb_isRecovered.at(i));
        }
        std::vector<int> ECalIndices = ReorderByEnergy(ecal_energy);
        for(auto& index : ECalIndices ){
          // double check to be sure there is  no calorechit  energy == 0
          if(ecal_energy.at(index)==0) continue;
          out.caloRecHit_e_isEndCap.push_back(ecal_isEndCap.at(index));
          out.caloRecHit_e_rho.push_back(ecal_rho.at(index));
          out.caloRecHit_e_DeltaEta.push_back(ecal_DeltaEta.at(index));
          out.caloRecHit_e_DeltaPhi.push_back(ecal_DeltaPhi.at(index));
          out.caloRecHit_e_DeltaR.push_back(ecal_DeltaR.at(index));
          out.caloRecHit_e_energy.push_back(ecal_energy.at(index));
          out.caloRecHit_e_chi2.push_back(ecal_chi2.at(index));
          out.caloRecHit_e_isRecovered.push_back(ecal_isRecovered.at(index));
        } 
        std::vector<bool> hcal_HadronSubDet;
        std::vector<float> hcal_rho;
        std::vector<float> hcal_DeltaEta;
        std::vector<float> hcal_DeltaPhi;
        std::vector<float> hcal_DeltaR;
        std::vector<float> hcal_energy;
        std::vector<float> hcal_chi2;
        std::vector<bool> hcal_time;
        for(std::vector<float>::size_type i=0; i < tau.caloRecHit_hbhe_energy.size() ; i++){
          if(tau.caloRecHit_hbhe_energy.at(i)==0) continue;
          hcal_HadronSubDet.push_back(true);
          hcal_rho.push_back(tau.caloRecHit_hbhe_rho.at(i));
          hcal_DeltaEta.push_back(DeltaEta(tau.caloRecHit_hbhe_eta.at(i), tau.l1Tau_eta));
          hcal_DeltaPhi.push_back(DeltaPhi(tau.caloRecHit_hbhe_phi.at(i), tau.l1Tau_phi));
          hcal_DeltaR.push_back(DeltaR(tau.caloRecHit_hbhe_phi.at(i), tau.caloRecHit_hbhe_eta.at(i), tau.l1Tau_phi, tau.l1Tau_eta));
          hcal_energy.push_back(tau.caloRecHit_hbhe_energy.at(i));
          hcal_chi2.push_back(tau.caloRecHit_hbhe_chi2.at(i));
          hcal_time.push_back(tau.caloRecHit_hbhe_time.at(i));
        }
        for(std::vector<float>::size_type i=0; i < tau.caloRecHit_ho_energy.size() ; i++){
          if(tau.caloRecHit_ho_energy.at(i)==0) continue;
          hcal_HadronSubDet.push_back(false);
          hcal_rho.push_back(tau.caloRecHit_ho_rho.at(i));
          hcal_DeltaEta.push_back(DeltaEta(tau.caloRecHit_ho_eta.at(i), tau.l1Tau_eta));
          hcal_DeltaPhi.push_back(DeltaPhi(tau.caloRecHit_ho_phi.at(i), tau.l1Tau_phi));
          hcal_DeltaR.push_back(DeltaR(tau.caloRecHit_ho_phi.at(i), tau.caloRecHit_ho_eta.at(i), tau.l1Tau_phi, tau.l1Tau_eta));
          hcal_energy.push_back(tau.caloRecHit_ho_energy.at(i));
          hcal_chi2.push_back(-100);
          hcal_time.push_back(tau.caloRecHit_ho_time.at(i));
        }
        std::vector<int> HCalIndices = ReorderByEnergy(hcal_energy);
        for(auto& index : HCalIndices ){
          // double check to be sure there is  no calorechit  energy == 0
          if(hcal_energy.at(index)==0) continue;
          out.caloRecHit_had_HadronSubDet.push_back(hcal_HadronSubDet.at(index));
          out.caloRecHit_had_rho.push_back(hcal_rho.at(index));
          out.caloRecHit_had_DeltaEta.push_back(hcal_DeltaEta.at(index));
          out.caloRecHit_had_DeltaPhi.push_back(hcal_DeltaPhi.at(index));
          out.caloRecHit_had_DeltaR.push_back(hcal_DeltaR.at(index));
          out.caloRecHit_had_energy.push_back(hcal_energy.at(index));
          out.caloRecHit_had_chi2.push_back(hcal_chi2.at(index));
          out.caloRecHit_had_time.push_back(hcal_time.at(index));
        }
        std::vector<float> patapt;
        std::vector<float> pataeta;
        std::vector<float> pataphi;
        std::vector<float> pataDeltaEta;
        std::vector<float> pataDeltaPhi;
        std::vector<float> pataDeltaR;
        std::vector<float> patachi2;
        std::vector<int> patandof;
        std::vector<int> patacharge;
        std::vector<int> pataquality;
        std::vector<float> patadxy;
        std::vector<float> patadz;
        std::vector<bool> patahasVertex;
        std::vector<float> patavert_z;
        std::vector<float> patavert_weight;
        std::vector<float> patavert_ptv2;
        std::vector<float> patavert_chi2;
        std::vector<int> patavert_ndof;

        for(std::vector<float>::size_type i=0; i < tau.patatrack_pt.size() ; i++){

          patapt.push_back(tau.patatrack_pt.at(i));
          pataeta.push_back(tau.patatrack_eta.at(i));
          pataphi.push_back(tau.patatrack_phi.at(i));
          pataDeltaEta.push_back(DeltaEta(tau.patatrack_eta.at(i), tau.l1Tau_eta));
          pataDeltaPhi.push_back(DeltaPhi(tau.patatrack_phi.at(i), tau.l1Tau_phi));
          pataDeltaR.push_back(DeltaR(tau.patatrack_phi.at(i),tau.patatrack_eta.at(i), tau.l1Tau_phi, tau.l1Tau_eta));
          patachi2.push_back(tau.patatrack_chi2.at(i));
          patandof.push_back(tau.patatrack_ndof.at(i));
          patacharge.push_back(tau.patatrack_charge.at(i));
          pataquality.push_back(tau.patatrack_quality.at(i));
          patadxy.push_back(tau.patatrack_dxy.at(i));
          patadz.push_back(tau.patatrack_dz.at(i));
          int vertex_id = tau.patatrack_vertex_id.at(i) ;
          if(vertex_id != -1){
            patahasVertex.push_back(true);
            patavert_z.push_back(tau.patavert_z.at(vertex_id));
            patavert_weight.push_back(tau.patavert_weight.at(vertex_id));
            patavert_ptv2.push_back(tau.patavert_ptv2.at(vertex_id));
            patavert_chi2.push_back(tau.patavert_chi2.at(vertex_id));
            patavert_ndof.push_back(tau.patavert_ndof.at(vertex_id));
          }
          else{
            patahasVertex.push_back(false);
            patavert_z.push_back(-100);
            patavert_weight.push_back(-100);
            patavert_ptv2.push_back(-100);
            patavert_chi2.push_back(-100);
            patavert_ndof.push_back(-100);
          }
        }
        std::vector<int> PataIndices = ReorderByEnergy(patapt);
        for(auto& index : PataIndices ){
            out.patatrack_pt.push_back(patapt.at(index));
            out.patatrack_eta.push_back(pataeta.at(index));
            out.patatrack_phi.push_back(pataphi.at(index));
            out.patatrack_DeltaEta.push_back(pataDeltaEta.at(index));
            out.patatrack_DeltaPhi.push_back(pataDeltaPhi.at(index));
            out.patatrack_DeltaR.push_back(pataDeltaR.at(index));
            out.patatrack_chi2.push_back(patachi2.at(index));
            out.patatrack_ndof.push_back(patandof.at(index));
            out.patatrack_charge.push_back(patacharge.at(index));
            out.patatrack_quality.push_back(pataquality.at(index));
            out.patatrack_dxy.push_back(patadxy.at(index));
            out.patatrack_dz.push_back(patadz.at(index));
            out.patatrack_hasVertex.push_back(patahasVertex.at(index));
            out.patatrack_vert_z.push_back(patavert_z.at(index));
            out.patatrack_vert_weight.push_back(patavert_weight.at(index));
            out.patatrack_vert_ptv2.push_back(patavert_ptv2.at(index));
            out.patatrack_vert_chi2.push_back(patavert_chi2.at(index));
            out.patatrack_vert_ndof.push_back(patavert_ndof.at(index));
        }
        intermTuple.Fill();
      }
      #undef CP_BR

private:
    const Arguments args;
    std::shared_ptr<TFile> inputFile, outputFile;
    TrainingTauTuple trainingTauTuple;
    IntermTuple intermTuple;

};

} // namespace analysis

PROGRAM_MAIN(analysis::IntermediateTupleProducer, Arguments)
