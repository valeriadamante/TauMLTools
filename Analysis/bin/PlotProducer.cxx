/*! Produce training tuple from tau tuple.
*/

#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/variadic.hpp>
#include <boost/foreach.hpp>
#include <boost/math/constants/constants.hpp>
#include "TauMLTools/Core/interface/program_main.h"
#include "TauMLTools/Core/interface/AnalysisMath.h"
#include "TauMLTools/Core/interface/RootExt.h"
#include "TauMLTools/Analysis/interface/L2EventTuple.h"
#include "TauMLTools/Analysis/interface/L2TrainingTauTuple.h"
#include "TauMLTools/Core/interface/ProgressReporter.h"
#include "TauMLTools/Analysis/interface/GenLepton.h"
#include "TH1D.h"

struct Arguments {
    run::Argument<std::string> input{"input", "input root file with tau tuple"};
    run::Argument<std::string> output{"output", "output root file with training tuple"};
    run::Argument<unsigned> n_threads{"n-threads", "number of threads", 1};
    run::Argument<Long64_t> end_entry{"end-entry", "end entry", std::numeric_limits<Long64_t>::max()};
    run::Argument<Long64_t> start_entry{"start-entry", "start entry", 0};
    //run::Argument<float> training_weight_factor{"training-weight-factor",
    //    "additional factor to the normalization of the training weights", 4.f};
    //run::Argument<int> parity{"parity", "take odd (parity=1), even (parity=0) or all (parity=-1) events", -1};
    //run::Argument<bool> isQCD{"isQCD", "is monte carlo sample", false};
};


namespace analysis {

class PlotProducer {
public:
    using Tau = train_tuple::Tau;
    using L2TauTuple = train_tuple::TrainTuple;
    using GenLep = reco_tau::gen_truth::GenLepton;

    PlotProducer(const Arguments& _args) :
        args(_args), inputFile(root_ext::OpenRootFile(args.input())), outputFile(root_ext::CreateRootFile(args.output(), ROOT::kLZ4, 4)), inputTuple("taus",inputFile.get(), true)
    {

        if(args.n_threads() > 1)
            ROOT::EnableImplicitMT(args.n_threads());

    }
    void Run()
    {
        const Long64_t end_entry = std::min(inputTuple.GetEntries(), args.end_entry());
        TH1D *histo = new TH1D("Epi","Epi",100, -1.01, 1.01);
        size_t n_processed = 0, n_total = static_cast<size_t>(end_entry - args.start_entry());
        tools::ProgressReporter reporter(10, std::cout, "Creating histogram...");
        reporter.SetTotalNumberOfEvents(n_total);
        for(Long64_t current_entry = args.start_entry(); current_entry < end_entry; ++current_entry) {
            //std::cout << "entry = " << current_entry << std::endl;
            inputTuple.GetEntry(current_entry);
            if (n_processed == 14000 || n_processed == 14000-1 || n_processed == 14000+1 ) continue;
            const auto& tau = inputTuple.data();
            FillTauBranches(tau, *histo);
            if(++n_processed % 1000 == 0)
                reporter.Report(n_processed);
        }
        reporter.Report(n_processed, true);
        histo->Write();
        std::cout << "Histogram successfully stored in " << args.output() << "." << std::endl;
    }


private :

    void FillTauBranches(const Tau& tau, TH1D &histo)
    {
      //std::cout << "genLeptons = " << tau.genLepton_nParticles.size() << std::endl;
      //std::cout << "genParticles = " << tau.genParticle_pdgId.size() << std::endl;
      int particle_counter =0;
      for(std::vector<int>::size_type i =0; i<tau.genLepton_nParticles.size(); i++){
        int part_index = 0;
        int lastMotherIndex =0 ;
        if(tau.genLepton_vis_pt.at(i)<15){
          while (part_index < tau.genLepton_nParticles.at(i) ){
            particle_counter++;
            part_index++;
            }
        }
        else {
          //std::cout << "lepton " << i << std::endl;
          //std::cout << "lepton Pt " << tau.genLepton_vis_pt.at(i) << std::endl;
          std::vector<int> genParticle_pdgId;
          std::vector<Long64_t> genParticle_mother;
          std::vector<int> genParticle_charge;
          std::vector<int> genParticle_isFirstCopy;
          std::vector<int> genParticle_isLastCopy;
          std::vector<float> genParticle_pt;
          std::vector<float> genParticle_eta;
          std::vector<float> genParticle_phi;
          std::vector<float> genParticle_mass;
          std::vector<float> genParticle_vtx_x;
          std::vector<float> genParticle_vtx_y;
          std::vector<float> genParticle_vtx_z;
          float EPiPlus=0.;
          float EPiZero=0.;
          while (part_index < tau.genLepton_nParticles.at(i) ){
            //std::cout << "particle " << particle_counter  << std::endl;
            lastMotherIndex = tau.genLepton_lastMotherIndex.at(i);
            genParticle_pdgId.push_back(tau.genParticle_pdgId.at(particle_counter));
            genParticle_mother.push_back(tau.genParticle_mother.at(particle_counter));
            genParticle_charge.push_back(tau.genParticle_charge.at(particle_counter));
            genParticle_isFirstCopy.push_back(tau.genParticle_isFirstCopy.at(particle_counter));
            genParticle_isLastCopy.push_back(tau.genParticle_isLastCopy.at(particle_counter));
            genParticle_pt.push_back(tau.genParticle_pt.at(particle_counter));
            genParticle_eta.push_back(tau.genParticle_eta.at(particle_counter));
            genParticle_phi.push_back(tau.genParticle_phi.at(particle_counter));
            genParticle_mass.push_back(tau.genParticle_mass.at(particle_counter));
            genParticle_vtx_x.push_back(tau.genParticle_vtx_x.at(particle_counter));
            genParticle_vtx_y.push_back(tau.genParticle_vtx_y.at(particle_counter));
            genParticle_vtx_z.push_back(tau.genParticle_vtx_z.at(particle_counter));
            //std::cout << "pos " << particle_counter << std::endl;
            if(part_index == tau.genLepton_nParticles.at(i) -1){
              GenLep lepton = GenLep::fromRootTuple(lastMotherIndex, genParticle_pdgId, genParticle_mother, genParticle_charge, genParticle_isFirstCopy, genParticle_isLastCopy, genParticle_pt, genParticle_eta, genParticle_phi, genParticle_mass, genParticle_vtx_x, genParticle_vtx_y, genParticle_vtx_z);
              //std::cout << "A" << std::endl;
              if(lepton.nChargedHadrons() == 1 && lepton.nNeutralHadrons()==1){
                bool both_found = true;
                //std::cout << "interessante! "<< std::endl;
                auto hadrons = lepton.hadrons();
                if(hadrons.size()==2){
                  for ( auto& hadron: hadrons) {
                       if(hadron->charge>0){
                         EPiPlus = hadron->p4.energy();
                       }
                       else if(hadron->charge==0){
                         EPiZero = hadron->p4.energy();
                       }
                       else{
                         both_found = false;
                       };
                  }
                  if(EPiPlus+EPiZero != 0. && both_found){
                    float variable = (EPiPlus-EPiZero)/(EPiPlus+EPiZero);
                    histo.Fill(variable);
                  }
                }
              }
            }
            particle_counter++;
            part_index++;
            }
          }
        }


    }


private:
    const Arguments args;
    std::shared_ptr<TFile> inputFile, outputFile;
    L2TauTuple inputTuple;

};

} // namespace analysis

PROGRAM_MAIN(analysis::PlotProducer, Arguments)
