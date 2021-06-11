/*
Filter for L2 hadronic tau selection
*/

#include "Compression.h"
// maybe I need them ? - utilities
//#include "DataFormats/Candidate/interface/Candidate.h"
//#include "DataFormats/Candidate/interface/CandidateFwd.h"
//#include "DataFormats/Common/interface/AssociationMap.h"
//#include "DataFormats/Candidate/interface/LeafCandidate.h"
//#include "DataFormats/Common/interface/Association.h"
//#include "DataFormats/Common/interface/AssociationVector.h"
//#include "DataFormats/Common/interface/HLTPathStatus.h"
//#include "DataFormats/Common/interface/HLTGlobalStatus.h"
//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
//#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/StringToEnumValue.h"

// other utilities
#include "TauMLTools/Analysis/interface/L2TauTagNNFilter.h"
//#include "TauMLTools/Analysis/interface/GenLepton.h"
//#include "TauMLTools/Core/interface/Tools.h"
//#include "TauMLTools/Core/interface/TextIO.h"
//#include "TauMLTools/Production/interface/TauAnalysis.h"
//#include "TauMLTools/Production/interface/TauJet.h"


L2TauNNTag::L2TauNNTag(const edm::ParameterSet& cfg):
    processName(cfg.getParameter<std::string>("processName")),
    //defaultDiTauPath(cfg.getParameter<std::string>("defaultDiTauPath")),
    //TR_token(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("TriggerResults"))),
    //puInfo_token(mayConsume<std::vector<PileupSummaryInfo>>(cfg.getParameter<edm::InputTag>("puInfo"))),
    l1Taus_token(consumes<l1t::TauBxCollection>(cfg.getParameter<edm::InputTag>("l1taus"))),
    hbhe_token(consumes<HBHERecHitCollection>(cfg.getParameter<edm::InputTag>("hbheInput"))),
    ho_token(consumes<HORecHitCollection>(cfg.getParameter<edm::InputTag>("hoInput"))),
    hf_token( consumes<HFRecHitCollection>(cfg.getParameter<edm::InputTag>("hfInput"))),
    ecalLabels(cfg.getParameter<std::vector<edm::InputTag> >("ecalInputs")),
    Geometry_token(esConsumes<CaloGeometry,CaloGeometryRecord>()),
    pataVertices_token(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("pataVertices"))),
    pataTracks_token(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("pataTracks"))),
    {
      const unsigned nLabels = ecalLabels.size();
      for (unsigned i = 0; i != nLabels; i++)
        ecal_tokens.push_back(consumes<EcalRecHitCollection>(ecalLabels[i]));
    }

bool L2TauNNTag::filter(const edm::Event& event, const edm::EventSetup& eventsetup) override
  {
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
      caloRecHitCollections AllCaloRecHits;
      AllCaloRecHits.hbhe= &*hbhe;
      AllCaloRecHits.ho= &*ho;
      AllCaloRecHits.hf= &*hf;
      AllCaloRecHits.eb= &*ebHandle;
      AllCaloRecHits.ee= &*eeHandle;
      AllCaloRecHits.Geometry = &*Geometry;

  }
virtual void endJob() override;

private:

template<typename Scalar>
static Scalar L2TauNNTag::DeltaPhi(Scalar phi1, Scalar phi2)
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
static Scalar L2TauNNTag::DeltaEta(Scalar eta1, Scalar eta2)
{
  return (eta1-eta2);
}

template<typename Scalar>
static Scalar L2TauNNTag::DeltaR(Scalar phi1,Scalar eta1,Scalar phi2,Scalar eta2) {
  auto float dphi = DeltaPhi(phi1, phi2);
  auto deta = DeltaEta(eta1, eta2);
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

void L2TauNNTag::FindObjectsAroundL1Tau(const AllCaloRecHits& caloRecHits, const reco::TrackCollection& patatracks,const reco::VertexCollection& patavertices, const l1t::TauBxCollection& l1Taus){

  float dR_max = 0.5;
  int nCellEta = 5;
  int nCellPhi = 5;
  int nVars = 31;
  float dEta_width = 2*dR_max/static_cast<float>(nCellEta);
  float dPhi_width = 2*dR_max/static_cast<float>(nCellPhi);
  // 1. reorder taus by Pt
  std::vector<float> l1Taus_pt;
  for(auto iter = l1Taus.begin(0); iter != l1Taus.end(0); ++iter) {
    l1Taus_pt.push_back(static_cast<float>(iter->polarP4().pt()));
  }
  // 2. create matrix
  std::vector<int> pt_indices=ReorderByEnergy(l1Taus_pt);
  int nTaus = pt_indices.size();
  float cellGridMatrix[nTaus][nCellEta][nCellPhi][nVars];
  // 2.1 fill matrix with global observables
  for (auto & pt_idx : pt_indices){
    for (int eta_idx = 0; eta_idx<5 ; eta_idx++){
      for (int phi_idx = 0; phi_idx<5 ; phi_idx++){
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::nVertices]=patavertices.size();
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::l1Tau_pt]=l1Taus.at(pt_idx)->polarP4().pt();
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::l1Tau_eta]=l1Taus.at(pt_idx)->polarP4().eta();
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::l1Tau_hwIso]=l1Taus.at(pt_idx).hwIso();
      }
    }
  }
  // 3. find objects  around l1 tau and fill matrix
  for (int tau_idx : pt_indices){
    float tauEta = l1Taus.at(pt_idx)->polarP4.eta();
    float tauPhi = l1Taus.at(pt_idx)->polarP4.phi();
    for (auto & caloRecHit_ee : *caloRecHits.ee){
      if(caloRecHit_ee.energy()<=0) continue;
      const auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_ee.id())->getPosition();
      float eeCalEta = position.eta();
      float eeCalPhi = position.phi();
      float eeCalEn = caloRecHit_ee.energy();
      float eeCalChi2 = caloRecHit_ee.chi2();
      if(DeltaR(eeCalPhi,eeCalEta,tauPhi,tauEta)<dR_max){
        float deta = DeltaEta(eeCalEta, tauEta);
        int eta_idx = int(math.floor((deta + dR_max) / dEta_width + dR_max));
        float dphi = DeltaPhi(eeCalPhi, tauPhi);
        int phi_idx = int(math.floor((dphi + dR_max) / dPhi_width + dR_max));
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalEnergySum]+=eeCalEn; //
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalSize]+=1;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalDeltaEta]+=deta;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalDeltaPhi]+=dphi;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalChi2]+=eeCalChi2; //
        if(eeCalChi2>0){
          cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalEnergySumForPositiveChi2]+=eeCalEn; //
          cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalSizeForPositiveChi2]+=1;
        }
        //EcalEnergyStdDev = 6,
      }
    }
    for (auto & caloRecHit_eb : *caloRecHits.eb){
      if(caloRecHit_eb.energy()<=0) continue;
      const auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_eb.id())->getPosition();
      float ebCalEta = position.eta();
      float ebCalPhi = position.phi();
      float ebCalEn = caloRecHit_eb.energy();
      float ebCalChi2 = caloRecHit_eb.chi2();
      if(DeltaR(ebCalPhi,ebCalEta,tauPhi,tauEta)<dR_max){
        float deta = DeltaEta(ebCalEta, tauEta);
        int eta_idx = int(math.floor((deta + dR_max) / dEta_width + dR_max));
        float dphi = DeltaPhi(ebCalPhi, tauPhi);
        int phi_idx = int(math.floor((dphi + dR_max) / dPhi_width + dR_max));
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalEnergySum]+=ebCalEn; //
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalSize]+=1;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalDeltaEta]+=deta;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalDeltaPhi]+=dphi;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalChi2]+=ebCalChi2; //
        if(ebCalChi2>0){
          cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalEnergySumForPositiveChi2]+=ebCalEn; //
          cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::EcalSizeForPositiveChi2]+=1;
        }
        //EcalEnergyStdDev = 6,
      }
    }

    for (auto & caloRecHit_hbhe : *caloRecHits.hbhe){
      if(caloRecHit_hbhe.energy()<=0) continue;
      const auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_hbhe.id())->getPosition();
      float hbheCalEta = position.eta();
      float hbheCalPhi = position.phi();
      float hbheCalEn = caloRecHit_hbhe.energy();
      float hbheCalChi2 = caloRecHit_hbhe.chi2();
      if(DeltaR(hbheCalPhi,hbheCalEta,tauPhi,tauEta)<dR_max){
        float deta = DeltaEta(hbheCalEta, tauEta);
        int eta_idx = int(math.floor((deta + dR_max) / dEta_width + dR_max));
        float dphi = DeltaPhi(hbheCalPhi, tauPhi);
        int phi_idx = int(math.floor((dphi + dR_max) / dPhi_width + dR_max));
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalEnergySum]+=hbheCalEn; //
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalSize]+=1;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalDeltaEta]+=deta;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalDeltaPhi]+=dphi;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalChi2]+=hbheCalChi2; //
        if(hbheCalChi2>0){
          cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalEnergySumForPositiveChi2]+=hbheCalEn; //
          cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalSizeForPositiveChi2]+=1;
        }
        //EcalEnergyStdDev = 6,
      }
    }

    for (auto & caloRecHit_ho : *caloRecHits.ho){
      if(caloRecHit_ho.energy()<=0) continue;
      const auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_ho.id())->getPosition();
      float hoCalEta = position.eta();
      float hoCalPhi = position.phi();
      float hoCalEn = caloRecHit_ho.energy();
      //float hoCalChi2 = caloRecHit_ho.chi2();
      if(DeltaR(hoCalPhi,hoCalEta,tauPhi,tauEta)<dR_max){
        float deta = DeltaEta(hoCalEta, tauEta);
        int eta_idx = int(math.floor((deta + dR_max) / dEta_width + dR_max));
        float dphi = DeltaPhi(hoCalPhi, tauPhi);
        int phi_idx = int(math.floor((dphi + dR_max) / dPhi_width + dR_max));
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalEnergySum]+=hoCalEn; //
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalSize]+=1;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalDeltaEta]+=deta;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalDeltaPhi]+=dphi;
        //cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalChi2]+=hoCalChi2; //
        //if(hoCalChi2>0){
        //  cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalEnergySumForPositiveChi2]+=eeCalEn; //
        //  cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::HCalSizeForPositiveChi2]+=1;
        //}
        //EcalEnergyStdDev = 6,
      }
    }

    for (unsigned n = 0; n < patatracks.size(); ++n){
      if(patatracks.at(n).pt()<=0) continue;
      float patatrackEta = patatracks.at(n).eta();
      float patatrackPhi = patatracks.at(n).phi();
      float patatrackPt = patatracks.at(n).pt();
      float patatrackNdof =patatracks.at(n).ndof();
      float patatrackChi2OverNdof =0;
      if(patatracks.at(n).ndof()!=0){PatatrackChi2OverNdof= patatracks.at(n).chi2()/ patatracks.at(n).ndof();}
      float patatrackDxy = patatracks.at(n).dxy();
      float patatrackDz = patatracks.at(n).dz();

      if(DeltaR(patatrackPhi,patatrackEta,tauPhi,tauEta)<dR_max){
        float deta = DeltaEta(hoCalEta, tauEta);
        int eta_idx = int(math.floor((deta + dR_max) / dEta_width + dR_max));
        float dphi = DeltaPhi(hoCalPhi, tauPhi);
        int phi_idx = int(math.floor((dphi + dR_max) / dPhi_width + dR_max));
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::PatatrackPtSum]+=patatrackPt; //
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::PatatrackSize]+=1;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::PatatrackChargeSum]+=patatracks.at(n).charge();
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::PatatrackDeltaEta]+=deta;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::PatatrackDeltaPhi]+=dphi;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::PatatrackChi2OverNdof]+=patatrackChi2OverNdof;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::PatatrackNdof]+=patatrackNdof;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::PatatrackDxy]+=patatrackDxy;
        cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::PatatrackDz]+=patatrackDz;
        if(FindVertexIndex(patavertices, patatracks.at(n) ) != -1){
          cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::PatatrackPtSumWithVertex]+=patatrackPt; //
          cellGridMatrix[pt_idx][eta_idx][phi_idx][NNInputs::PatatrackSizeWithVertex]+=1;
        }
      }
      // missing std dev
    }
  }

  // 2. find object around taus
  // 2.1
}



};


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L2TauNNTag);
