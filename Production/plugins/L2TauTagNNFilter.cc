/*
Filter for L2 hadronic tau selection
*/
// system include files
#include <memory>
#include <Math/VectorUtil.h>
#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/variadic.hpp>
#include <boost/math/constants/constants.hpp>
#include "Compression.h"
// user include files
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CommonTools/Utils/interface/StringToEnumValue.h"
// maybe I need them ? - utilities
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/HLTPathStatus.h"
#include "DataFormats/Common/interface/HLTGlobalStatus.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
// Geometry
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTowerTopology.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
// caloRecHit
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitDefs.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
//Tracks
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//vertices
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
// L1 Tau
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "TauMLTools/Core/interface/Tools.h"
#include "TauMLTools/Core/interface/TextIO.h"
#include "TauMLTools/Production/interface/L2TauTagNNFilter.h"
// other utilities

enum NNInputs {
  nVertices = 0,
  l1Tau_pt = 1,
  l1Tau_eta = 2,
  l1Tau_hwIso = 3,
  ECalEnergySum = 4,
  ECalSize = 5,
  ECalEnergyStdDev = 6,
  ECalDeltaEta = 7,
  ECalDeltaPhi = 8,
  ECalChi2 = 9,
  ECalEnergySumForPositiveChi2 = 10,
  ECalSizeForPositiveChi2 = 11,
  HCalEnergySum = 12,
  HCalSize = 13,
  HCalEnergyStdDev = 14,
  HCalDeltaEta = 15,
  HCalDeltaPhi = 16,
  HCalChi2 = 17,
  HCalEnergySumForPositiveChi2 = 18,
  HCalSizeForPositiveChi2 = 19,
  PatatrackPtSum = 20,
  PatatrackSize = 21,
  PatatrackSizeWithVertex = 22,
  PatatrackPtSumWithVertex = 23,
  PatatrackChargeSum = 24,
  PatatrackDeltaEta = 25,
  PatatrackDeltaPhi = 26,
  PatatrackChi2OverNdof = 27,
  PatatrackNdof = 28,
  PatatrackDxy = 29,
  PatatrackDz = 30
};

L2TauNNTag::L2TauNNTag(const edm::ParameterSet& cfg):
  processName(cfg.getParameter<std::string>("processName")),
  l1Taus_token(consumes<l1t::TauBxCollection>(cfg.getParameter<edm::InputTag>("l1taus"))),
  hbhe_token(consumes<HBHERecHitCollection>(cfg.getParameter<edm::InputTag>("hbheInput"))),
  ho_token(consumes<HORecHitCollection>(cfg.getParameter<edm::InputTag>("hoInput"))),
  hf_token( consumes<HFRecHitCollection>(cfg.getParameter<edm::InputTag>("hfInput"))),
  ecalLabels(cfg.getParameter<std::vector<edm::InputTag> >("ecalInputs")),
  Geometry_token(esConsumes<CaloGeometry,CaloGeometryRecord>()),
  pataVertices_token(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("pataVertices"))),
  pataTracks_token(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("pataTracks")))
  {
    const unsigned nLabels = ecalLabels.size();
    for (unsigned i = 0; i != nLabels; i++)
      ecal_tokens.push_back(consumes<EcalRecHitCollection>(ecalLabels[i]));

  }
L2TauNNTag::~L2TauNNTag() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}



float L2TauNNTag::DeltaPhi(Float_t phi1, Float_t phi2)
{
    static constexpr float pi = boost::math::constants::pi<float>();
    float dphi = phi1 - phi2;
    if(dphi > pi)
        dphi -= 2*pi;
    else if(dphi <= -pi)
        dphi += 2*pi;
    return dphi;
}

float L2TauNNTag::DeltaEta(Float_t eta1, Float_t eta2)
{
  return (eta1-eta2);
}

float L2TauNNTag::DeltaR(Float_t phi1,Float_t eta1,Float_t phi2,Float_t eta2) {
  float dphi = DeltaPhi(phi1, phi2);
  float deta = DeltaEta(eta1, eta2);
  //std::cout << "deltaR = " << (std::sqrt(deta * deta + dphi * dphi)) <<std::endl;
  return (std::sqrt(deta * deta + dphi * dphi));

}

std::vector<int> L2TauNNTag::ReorderByEnergy(std::vector<float>& energy){
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

int L2TauNNTag::FindVertexIndex(const reco::VertexCollection& vertices, const reco::Track& track){
  const reco::TrackBase *track_to_compare = &track;
  for(size_t n = 0; n< vertices.size() ; ++n){
    for(auto k = vertices.at(n).tracks_begin(); k != vertices.at(n).tracks_end(); ++k){
        if(&**k==track_to_compare) return n;
    }
  }
  return -1;
}

void L2TauNNTag::FindObjectsAroundL1Tau(const caloRecHitCollections& caloRecHits, const reco::TrackCollection& patatracks,const reco::VertexCollection& patavertices, const l1t::TauBxCollection& l1Taus){
  float dR_max = 0.5;
  int nCellEta = 5;
  int nCellPhi = 5;
  int nVars = 31;
  float dEta_width = 2*dR_max/static_cast<float>(nCellEta);
  float dPhi_width = 2*dR_max/static_cast<float>(nCellPhi);
  // 1. reorder taus by Pt
  std::cout<< "delta eta = "<< dEta_width;
  std::cout<< "\ndelta phi = "<< dPhi_width;
  std::cout<< "\nnVars = "<< nVars;
  std::vector<float> l1Taus_pt;
  std::vector<float> l1Taus_eta;
  std::vector<float> l1Taus_phi;
  std::vector<float> l1Taus_hwIso;
  for(auto iter = l1Taus.begin(0); iter != l1Taus.end(0); ++iter) {
    l1Taus_pt.push_back(static_cast<float>(iter->polarP4().pt()));
    l1Taus_eta.push_back(static_cast<float>(iter->polarP4().eta()));
    l1Taus_phi.push_back(static_cast<float>(iter->polarP4().phi()));
    l1Taus_phi.push_back(static_cast<float>(iter->hwIso()));
  }
  // 2. create matrix
  std::vector<int> pt_indices=ReorderByEnergy(l1Taus_pt);
  int nTaus = pt_indices.size();
  float cellGridMatrix[nTaus][nCellEta][nCellPhi][nVars];
  // 3. find objects  around l1 tau and fill matrix
  for (int tau_idx : pt_indices){
    // 3.1 fill matrix with global observables
    for (int eta_idx = 0; eta_idx<5 ; eta_idx++){
      for (int phi_idx = 0; phi_idx<5 ; phi_idx++){
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::nVertices]=static_cast<float>(patavertices.size());
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::l1Tau_pt] = l1Taus_pt.at(tau_idx);
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::l1Tau_eta] = l1Taus_eta.at(tau_idx);
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::l1Tau_hwIso] = l1Taus_hwIso.at(tau_idx);
      }
    }
  float tauEta = l1Taus_eta.at(tau_idx);
  float tauPhi = l1Taus_phi.at(tau_idx);
    for (auto & caloRecHit_ee : *caloRecHits.ee){
      if(caloRecHit_ee.energy()<=0) continue;
      const auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_ee.id())->getPosition();
      float eeCalEta = position.eta();
      float eeCalPhi = position.phi();
      float eeCalEn = caloRecHit_ee.energy();
      float eeCalChi2 = caloRecHit_ee.chi2();
      if(DeltaR(eeCalPhi,eeCalEta,tauPhi,tauEta)<dR_max){
        auto deta = DeltaEta(eeCalEta, tauEta);
        int eta_idx = static_cast<int>(floor((deta + dR_max) / dEta_width + dR_max));
        auto dphi = DeltaPhi(eeCalPhi, tauPhi);
        int phi_idx = static_cast<int>(floor((dphi + dR_max) / dPhi_width + dR_max));
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalEnergySum]+=eeCalEn; //
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalSize]+=1;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalDeltaEta]+=deta;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalDeltaPhi]+=dphi;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalChi2]+=eeCalChi2; //
        if(eeCalChi2>0){
          cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalEnergySumForPositiveChi2]+=eeCalEn; //
          cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalSizeForPositiveChi2]+=1;
        }
        //EcalEnergyStdDev = 6,
      }
    } // end of loop over calorechits ee

    for (auto & caloRecHit_eb : *caloRecHits.eb){
      if(caloRecHit_eb.energy()<=0) continue;
      const auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_eb.id())->getPosition();
      float ebCalEta = position.eta();
      float ebCalPhi = position.phi();
      float ebCalEn = caloRecHit_eb.energy();
      float ebCalChi2 = caloRecHit_eb.chi2();
      if(DeltaR(ebCalPhi,ebCalEta,tauPhi,tauEta)<dR_max){
        auto deta = DeltaEta(ebCalEta, tauEta);
        int eta_idx = static_cast<int>(floor((deta + dR_max) / dEta_width + dR_max));
        auto dphi = DeltaPhi(ebCalPhi, tauPhi);
        int phi_idx = static_cast<int>(floor((dphi + dR_max) / dPhi_width + dR_max));
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalEnergySum]+=ebCalEn; //
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalSize]+=1;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalDeltaEta]+=deta;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalDeltaPhi]+=dphi;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalChi2]+=ebCalChi2; //
        if(ebCalChi2>0){
          cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalEnergySumForPositiveChi2]+=ebCalEn; //
          cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::ECalSizeForPositiveChi2]+=1;
        }
        //EcalEnergyStdDev = 6,
      }
    } // end of loop over calorechits eb


    for (auto & caloRecHit_hbhe : *caloRecHits.hbhe){
      if(caloRecHit_hbhe.energy()<=0) continue;
      const auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_hbhe.id())->getPosition();
      float hbheCalEta = position.eta();
      float hbheCalPhi = position.phi();
      float hbheCalEn = caloRecHit_hbhe.energy();
      float hbheCalChi2 = caloRecHit_hbhe.chi2();
      if(DeltaR(hbheCalPhi,hbheCalEta,tauPhi,tauEta)<dR_max){
        auto deta = DeltaEta(hbheCalEta, tauEta);
        int eta_idx = static_cast<int>(floor((deta + dR_max) / dEta_width + dR_max));
        auto dphi = DeltaPhi(hbheCalPhi, tauPhi);
        int phi_idx = static_cast<int>(floor((dphi + dR_max) / dPhi_width + dR_max));
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalEnergySum]+=hbheCalEn; //
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalSize]+=1;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalDeltaEta]+=deta;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalDeltaPhi]+=dphi;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalChi2]+=hbheCalChi2; //
        if(hbheCalChi2>0){
          cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalEnergySumForPositiveChi2]+=hbheCalEn; //
          cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalSizeForPositiveChi2]+=1;
        }
        //EcalEnergyStdDev = 6,
      }
    } // end of loop over calorechits hbhe

    for (auto & caloRecHit_ho : *caloRecHits.ho){
      if(caloRecHit_ho.energy()<=0) continue;
      const auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_ho.id())->getPosition();
      float hoCalEta = position.eta();
      float hoCalPhi = position.phi();
      float hoCalEn = caloRecHit_ho.energy();
      //float hoCalChi2 = caloRecHit_ho.chi2();
      if(DeltaR(hoCalPhi,hoCalEta,tauPhi,tauEta)<dR_max){
        auto deta = DeltaEta(hoCalEta, tauEta);
        int eta_idx = static_cast<int>(floor((deta + dR_max) / dEta_width + dR_max));
        auto dphi = DeltaPhi(hoCalPhi, tauPhi);
        int phi_idx = static_cast<int>(floor((dphi + dR_max) / dPhi_width + dR_max));
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalEnergySum]+=hoCalEn; //
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalSize]+=1;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalDeltaEta]+=deta;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalDeltaPhi]+=dphi;
        //cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalChi2]+=hoCalChi2; //
        //if(hoCalChi2>0){
        //  cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalEnergySumForPositiveChi2]+=hoCalEn; //
        //  cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::HCalSizeForPositiveChi2]+=1;
        //}
        //EcalEnergyStdDev = 6,
      }
    } // end of loop over calorechits ho



    for (unsigned n = 0; n < patatracks.size(); ++n){
      if(patatracks.at(n).pt()<=0) continue;
      float patatrackEta = patatracks.at(n).eta();
      float patatrackPhi = patatracks.at(n).phi();
      float patatrackPt = patatracks.at(n).pt();
      float patatrackNdof =patatracks.at(n).ndof();
      float patatrackChi2OverNdof ;
      if(patatracks.at(n).ndof() == 0.){
        //PatatrackChi2OverNdof = patatracks.at(n).chi2()/patatracks.at(n).ndof();
        patatrackChi2OverNdof =0.;

      }
      else{
          patatrackChi2OverNdof = patatracks.at(n).chi2()/(patatracks.at(n).ndof());
      }
      float patatrackDxy = patatracks.at(n).dxy();
      float patatrackDz = patatracks.at(n).dz();

      if(DeltaR(patatrackPhi,patatrackEta,tauPhi,tauEta)<dR_max){
        auto deta = DeltaEta(patatrackEta, tauEta);
        int eta_idx = int(floor((deta + dR_max) / dEta_width + dR_max));
        auto dphi = DeltaPhi(patatrackPhi, tauPhi);
        int phi_idx = int(floor((dphi + dR_max) / dPhi_width + dR_max));
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::PatatrackPtSum]+=patatrackPt; //
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::PatatrackSize]+=1;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::PatatrackChargeSum]+=patatracks.at(n).charge();
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::PatatrackDeltaEta]+=deta;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::PatatrackDeltaPhi]+=dphi;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::PatatrackChi2OverNdof]+=patatrackChi2OverNdof;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::PatatrackNdof]+=patatrackNdof;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::PatatrackDxy]+=patatrackDxy;
        cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::PatatrackDz]+=patatrackDz;
        if(FindVertexIndex(patavertices, patatracks.at(n) ) != -1){
          cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::PatatrackPtSumWithVertex]+=patatrackPt; //
          cellGridMatrix[tau_idx][eta_idx][phi_idx][NNInputs::PatatrackSizeWithVertex]+=1;
        }
      }
      // missing std dev
    } // end of loop over patatracks

  } // end of loop over taus
  std::cout << "cosa " << cellGridMatrix[0][2][1][2]<<std::endl;
}

bool L2TauNNTag::filter(edm::Event& event, const edm::EventSetup& eventsetup) {

  bool result = true;
  edm::Handle<l1t::TauBxCollection> l1Taus;
  event.getByToken(l1Taus_token, l1Taus);
  edm::Handle<EcalRecHitCollection> ebHandle;
  edm::Handle<EcalRecHitCollection> eeHandle;
  for (std::vector<edm::EDGetTokenT<EcalRecHitCollection> >::const_iterator i = ecal_tokens.begin(); i != ecal_tokens.end(); i++) {
    edm::Handle<EcalRecHitCollection> ec_tmp;
    event.getByToken(*i, ec_tmp);
    if (ec_tmp->empty())
      continue;
    if ((ec_tmp->begin()->detid()).subdetId() == EcalBarrel) {
      ebHandle = ec_tmp;
    }
    else if ((ec_tmp->begin()->detid()).subdetId() == EcalEndcap) {
      eeHandle = ec_tmp;
    }
  }
  std::vector<edm::EDGetTokenT<EcalRecHitCollection> >::const_iterator i;
  for (i = ecal_tokens.begin(); i != ecal_tokens.end(); i++) {
    edm::Handle<EcalRecHitCollection> ec;
    event.getByToken(*i, ec);
  }
  edm::Handle<HBHERecHitCollection> hbhe;
  event.getByToken(hbhe_token, hbhe);
  edm::Handle<HORecHitCollection> ho;
  event.getByToken(ho_token, ho);
  edm::Handle<HFRecHitCollection> hf;
  event.getByToken(hf_token, hf);
  edm::Handle<reco::TrackCollection> pataTracks;
  event.getByToken(pataTracks_token, pataTracks);
  edm::Handle<reco::VertexCollection> pataVertices;
  event.getByToken(pataVertices_token, pataVertices);
  edm::ESHandle<CaloGeometry> Geometry = eventsetup.getHandle(Geometry_token);
  caloRecHitCollections caloRecHits;
  caloRecHits.hbhe= &*hbhe;
  caloRecHits.ho= &*ho;
  caloRecHits.hf= &*hf;
  caloRecHits.eb= &*ebHandle;
  caloRecHits.ee= &*eeHandle;
  caloRecHits.Geometry = &*Geometry;
  FindObjectsAroundL1Tau(caloRecHits, *pataTracks, *pataVertices,*l1Taus);

  return result;
}
// ------------ method called once each job just before starting event loop  ------------
void L2TauNNTag::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void L2TauNNTag::endJob() {}

//define this as a plug-in
DEFINE_FWK_MODULE(L2TauNNTag);
