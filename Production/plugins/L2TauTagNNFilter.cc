/*
Filter for L2 hadronic tau selection
*/
// system include files
#include <memory>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <Math/VectorUtil.h>
#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/variadic.hpp>
#include <boost/math/constants/constants.hpp>
#include "Compression.h"
// user include files
#include "FWCore/Framework/interface/stream/EDFilter.h"
//#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
// utilities
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
// Geometry
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
// caloRecHit
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EcalDetIdCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitDefs.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
//Tracks
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//vertices
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
// L1 Tau
#include "DataFormats/L1Trigger/interface/Tau.h"
namespace pt = boost::property_tree;
enum NNInputs {
  nVertices = 0,
  l1Tau_pt = 1,
  l1Tau_eta = 2,
  l1Tau_hwIso = 3,
  EcalEnergySum = 4,
  EcalSize = 5,
  EcalEnergyStdDev = 6,
  EcalDeltaEta = 7,
  EcalDeltaPhi = 8,
  EcalChi2 = 9,
  EcalEnergySumForPositiveChi2 = 10,
  EcalSizeForPositiveChi2 = 11,
  HcalEnergySum = 12,
  HcalSize = 13,
  HcalEnergyStdDev = 14,
  HcalDeltaEta = 15,
  HcalDeltaPhi = 16,
  HcalChi2 = 17,
  HcalEnergySumForPositiveChi2 = 18,
  HcalSizeForPositiveChi2 = 19,
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
std::map<int, std::string> varNameMap = {
  {0,"nVertices"},
  {1,"l1Tau_pt"},
  {2,"l1Tau_eta"},
  {3,"l1Tau_hwIso"},
  {4,"EcalEnergySum"},
  {5,"EcalSize"},
  {6,"EcalEnergyStdDev"},
  {7,"EcalDeltaEta"},
  {8,"EcalDeltaPhi"},
  {9,"EcalChi2"},
  {10,"EcalEnergySumForPositiveChi2"},
  {11,"EcalSizeForPositiveChi2"},
  {12,"HcalEnergySum"},
  {13,"HcalSize"},
  {14,"HcalEnergyStdDev"},
  {15,"HcalDeltaEta"},
  {16,"HcalDeltaPhi"},
  {17,"HcalChi2"},
  {18,"HcalEnergySumForPositiveChi2"},
  {19,"HcalSizeForPositiveChi2"},
  {20,"PatatrackPtSum"},
  {21,"PatatrackSize"},
  {22,"PatatrackSizeWithVertex"},
  {23,"PatatrackPtSumWithVertex"},
  {24,"PatatrackChargeSum"},
  {25,"PatatrackDeltaEta"},
  {26,"PatatrackDeltaPhi"},
  {27,"PatatrackChi2OverNdof"},
  {28,"PatatrackNdof"},
  {29,"PatatrackDxy"},
  {30,"PatatrackDz"},
};

struct CacheData {
  CacheData() : graphDef(nullptr) {}
  std::atomic<tensorflow::GraphDef*> graphDef;
};

class L2TauNNTag : public edm::stream::EDFilter<edm::GlobalCache<CacheData>> {
public:
  struct caloRecHitCollections {
    const HBHERecHitCollection  *hbhe;
    const HORecHitCollection *ho;
    const EcalRecHitCollection *eb;
    const EcalRecHitCollection *ee;
    const CaloGeometry *Geometry;
  };
  explicit L2TauNNTag(const edm::ParameterSet&, const CacheData*);
  ~L2TauNNTag() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions&);

  static std::unique_ptr<CacheData> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const CacheData*);

private:
  void beginJob();
  std::vector<int> get_tensor_shape(tensorflow::Tensor& tensor);
  static constexpr float pi = boost::math::constants::pi<float>();
  float DeltaPhi(Float_t phi1, Float_t phi2);
  float DeltaEta(Float_t eta1, Float_t eta2) ;
  float DeltaR(Float_t phi1,Float_t eta1,Float_t phi2,Float_t eta2);
  std::vector<int> ReorderByEnergy(std::vector<float>& energy);
  int FindVertexIndex(const reco::VertexCollection& vertices, const reco::Track& track);
  void initializeTensor(tensorflow::Tensor& tensor);
  void checknan(tensorflow::Tensor& tensor, bool printoutTensor);
  void standardizeTensor(tensorflow::Tensor& tensor, std::string normDictPath);
  void FindObjectsAroundL1Tau(const caloRecHitCollections& caloRecHits, const reco::TrackCollection& patatracks,const reco::VertexCollection& patavertices, const l1t::TauBxCollection& l1Taus, int evt_id,bool &result);
  bool filter(edm::Event& event, const edm::EventSetup& eventsetup) ;
  void endJob() ;

private:
  std::string processName;
  edm::EDGetTokenT<l1t::TauBxCollection> l1Taus_token;
  edm::EDGetTokenT<HBHERecHitCollection> hbhe_token;
  edm::EDGetTokenT<HORecHitCollection> ho_token;
  std::vector<edm::InputTag> ecalLabels;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> Geometry_token;
  edm::EDGetTokenT<reco::VertexCollection> pataVertices_token;
  edm::EDGetTokenT<reco::TrackCollection> pataTracks_token;
  std::string normalizationDict_;
  float discr_threshold_;
  std::vector<edm::EDGetTokenT<EcalRecHitCollection>> ecal_tokens;
  std::string graphPath_;
  std::string inputTensorName_;
  std::string outputTensorName_;
  tensorflow::Session* session_;

};


std::unique_ptr<CacheData> L2TauNNTag::initializeGlobalCache(const edm::ParameterSet& cfg) {
  // this method is supposed to create, initialize and return a CacheData instance
  CacheData* cacheData = new CacheData();

  // load the graph def and save it
  std::string graphPath = cfg.getParameter<std::string>("graphPath");
  cacheData->graphDef = tensorflow::loadGraphDef(graphPath);

  // set tensorflow log leven to warning
  tensorflow::setLogging("2");

  return std::unique_ptr<CacheData>(cacheData);
}

void L2TauNNTag::globalEndJob(const CacheData* cacheData) {
  // reset the graphDef
  if (cacheData->graphDef != nullptr) {
    delete cacheData->graphDef;
  }
}

void L2TauNNTag::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // defining this function will lead to a *_cfi file being generated when compiling
  edm::ParameterSetDescription desc;
  desc.add<std::string>("processName");
  desc.add<edm::InputTag>("l1taus");
  desc.add<edm::InputTag>("hbheInput");
  desc.add<edm::InputTag>("hoInput");
  desc.add<std::vector<edm::InputTag> >("ecalInputs");
  desc.add<edm::InputTag>("pataVertices");
  desc.add<edm::InputTag>("pataTracks");
  desc.add<std::string>("normalizationDict");
  desc.add<double>("discr_threshold");
  desc.add<std::string>("graphPath");
  descriptions.addWithDefaultLabel(desc);
}


L2TauNNTag::L2TauNNTag(const edm::ParameterSet& cfg, const CacheData* cacheData):
      processName(cfg.getParameter<std::string>("processName")),
      l1Taus_token(consumes<l1t::TauBxCollection>(cfg.getParameter<edm::InputTag>("l1taus"))),
      hbhe_token(consumes<HBHERecHitCollection>(cfg.getParameter<edm::InputTag>("hbheInput"))),
      ho_token(consumes<HORecHitCollection>(cfg.getParameter<edm::InputTag>("hoInput"))),
      ecalLabels(cfg.getParameter<std::vector<edm::InputTag> >("ecalInputs")),
      Geometry_token(esConsumes<CaloGeometry,CaloGeometryRecord>()),
      pataVertices_token(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("pataVertices"))),
      pataTracks_token(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("pataTracks"))),
      normalizationDict_(cfg.getParameter<std::string>("normalizationDict")),
      discr_threshold_(cfg.getParameter<double>("discr_threshold")),
      graphPath_(cfg.getParameter<std::string>("graphPath")),
      inputTensorName_((cacheData->graphDef.load())->node(0).name()),
      outputTensorName_((cacheData->graphDef.load())->node(cacheData->graphDef.load()->node_size()-1).name()),
      session_(tensorflow::createSession(cacheData->graphDef))
      {
        const unsigned nLabels = ecalLabels.size();
        for (unsigned i = 0; i != nLabels; i++)
          ecal_tokens.push_back(consumes<EcalRecHitCollection>(ecalLabels[i]));
      }

void L2TauNNTag::beginJob() {}

void L2TauNNTag::endJob() {
  // close the session
  tensorflow::closeSession(session_);
}


std::vector<int> L2TauNNTag::get_tensor_shape(tensorflow::Tensor& tensor) {
    std::vector<int> shape;
    int num_dimensions = tensor.shape().dims();
    for(int ii_dim=0; ii_dim<num_dimensions; ii_dim++) {
        shape.push_back(tensor.shape().dim_size(ii_dim));
    }
    return shape;
}
float L2TauNNTag::DeltaPhi(Float_t phi1, Float_t phi2) {
    static constexpr float pi = boost::math::constants::pi<float>();
  float dphi = phi1 - phi2;
    if(dphi > pi)
        dphi -= 2*pi;
    else if(dphi <= -pi)
        dphi += 2*pi;
    return dphi;
}
float L2TauNNTag::DeltaEta(Float_t eta1, Float_t eta2)  {
  return (eta1-eta2);
}
float L2TauNNTag::DeltaR(Float_t phi1,Float_t eta1,Float_t phi2,Float_t eta2) {
  float dphi = DeltaPhi(phi1, phi2);
  float deta = DeltaEta(eta1, eta2);
  return (std::sqrt(deta * deta + dphi * dphi));

}
std::vector<int> L2TauNNTag::ReorderByEnergy(std::vector<float>& energy){
    std::vector<int> indices;
    while(indices.size()<energy.size()){
      float current_energy_value=0.;
      int current_index = -100;
      for(std::vector<float>::size_type i=0; i < energy.size() ; i++){
        if(std::find(indices.begin(), indices.end(), i)!= indices.end()) continue;
        if(energy.at(i) >= current_energy_value ){
          current_energy_value = energy.at(i);
          current_index = static_cast<int>(i);
        }
      }
      if(current_index>=0){
        indices.push_back(current_index);
      }
    }
    return indices;
 }
int L2TauNNTag::FindVertexIndex(const reco::VertexCollection& vertices, const reco::Track& track){
  const reco::TrackBase *track_to_compare = &track;
  for(size_t n = 0; n< vertices.size() ; ++n){
    for(auto k = vertices.at(n).tracks_begin(); k != vertices.at(n).tracks_end(); ++k){
      if(&**k==track_to_compare){
        return n;
      }
    }
  }
  return -1;
}
void L2TauNNTag::initializeTensor(tensorflow::Tensor& tensor){
  std::vector<int> tensor_shape = get_tensor_shape(tensor);
  for(int tau_idx =0; tau_idx < tensor_shape.at(0); tau_idx++){
    for(int eta_idx =0; eta_idx < tensor_shape.at(1); eta_idx++){
      for(int phi_idx =0; phi_idx < tensor_shape.at(2); phi_idx++){
        for(int var_idx =0; var_idx < tensor_shape.at(3); var_idx++){
          tensor.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, var_idx)=static_cast<float>(0);
        }
      }
    }
  }
}
void L2TauNNTag::checknan(tensorflow::Tensor& tensor, bool printoutTensor){
  std::vector<int> tensor_shape = get_tensor_shape(tensor);
  for(int tau_idx =0; tau_idx < tensor_shape.at(0); tau_idx++){
    for(int eta_idx =0; eta_idx < tensor_shape.at(1); eta_idx++){
      for(int phi_idx =0; phi_idx < tensor_shape.at(2); phi_idx++){
        for(int var_idx =0; var_idx < tensor_shape.at(3); var_idx++){
          auto nonstd_var = tensor.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, var_idx);
          if(std::isnan(nonstd_var)){
            std::cout << "var is nan \nvar name= " << varNameMap.at(var_idx) << "\t tau_idx = " << tau_idx << "\t eta_idx = " << eta_idx << "\t phi_idx = " << phi_idx << std::endl;
            std::cout << "other vars in same cell \n";
            if(var_idx+1 < tensor_shape.at(3)) std::cout << varNameMap.at(var_idx+1) << "\t = " <<tensor.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, var_idx+1)<< std::endl;
            if(var_idx+2 < tensor_shape.at(3)) std::cout << varNameMap.at(var_idx+2) << "\t = " <<tensor.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, var_idx+2)<< std::endl;
            if(var_idx+3 < tensor_shape.at(3)) std::cout << varNameMap.at(var_idx+3) << "\t = " <<tensor.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, var_idx+3)<< std::endl;
            if(var_idx+4 < tensor_shape.at(3)) std::cout << varNameMap.at(var_idx+4) << "\t = " <<tensor.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, var_idx+4)<< std::endl;
          }
        }
      }
    }
  }
  if(printoutTensor){ /* this is only for debugging */
    for(int tau_idx =0; tau_idx < tensor_shape.at(0); tau_idx++){
      for(int phi_idx =0; phi_idx < tensor_shape.at(1); phi_idx++){
        for(int eta_idx =0; eta_idx < tensor_shape.at(2); eta_idx++){
          for(int var_idx =0; var_idx < tensor_shape.at(3); var_idx++){
            if(tensor.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::l1Tau_pt)*256.0 == 113.5){
              std::cout << tensor.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, var_idx) <<",\n";
              //std::cout << "\nvar name= " << varNameMap.at(var_idx) << "\t tau_idx = " << tau_idx << "\t eta_idx = " << eta_idx << "\t phi_idx = " << phi_idx << " \tvalue =" << tensor.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, var_idx);
            }
          }
        }
      }
    }
  }
}
void L2TauNNTag::standardizeTensor(tensorflow::Tensor& tensor, std::string normDictPath){
  pt::ptree loadPtreeRoot;
  pt::read_json(normDictPath, loadPtreeRoot);
  std::vector<int> tensor_shape = get_tensor_shape(tensor);
  for(int tau_idx =0; tau_idx < tensor_shape.at(0); tau_idx++){
    for(int  phi_idx=0; phi_idx < tensor_shape.at(1); phi_idx++){
      for(int eta_idx =0; eta_idx < tensor_shape.at(2); eta_idx++){
        for(int var_idx =0; var_idx < tensor_shape.at(3); var_idx++){
          pt::ptree var = loadPtreeRoot.get_child(varNameMap.at(var_idx));
          float mean = var.get_child("mean").get_value<float>();
          float std =var.get_child("std").get_value<float>();
          float min =var.get_child("min").get_value<float>();
          float max =var.get_child("max").get_value<float>();
          float nonstd_var = tensor.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, var_idx);
          float std_var = static_cast<float>((nonstd_var-mean)/std);
          if(std_var > max ){
            std_var = static_cast<float>(max);
          }
          else if (std_var < min){
            std_var = static_cast<float>(min);
          }
          tensor.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, var_idx) = static_cast<float>(std_var);
        }
      }
    }
  }
}
void L2TauNNTag::FindObjectsAroundL1Tau(const caloRecHitCollections& caloRecHits, const reco::TrackCollection& patatracks,const reco::VertexCollection& patavertices, const l1t::TauBxCollection& l1Taus, int evt_id,bool &result){

  float dR_max = 0.5;
  int nCellEta = 5;
  int nCellPhi = 5;
  int nVars = 31;
  float dEta_width = 2*dR_max/static_cast<float>(nCellEta);
  float dPhi_width = 2*dR_max/static_cast<float>(nCellPhi);

  std::vector<float> l1Taus_pt;
  std::vector<float> l1Taus_eta;
  std::vector<float> l1Taus_phi;
  std::vector<float> l1Taus_hwIso;
  for(auto iter = l1Taus.begin(0); iter != l1Taus.end(0); ++iter) {
    l1Taus_pt.push_back(static_cast<float>(iter->polarP4().pt()));
    l1Taus_eta.push_back(static_cast<float>(iter->polarP4().eta()));
    l1Taus_phi.push_back(static_cast<float>(iter->polarP4().phi()));
    l1Taus_hwIso.push_back(static_cast<float>(iter->hwIso()));
  }

  std::vector<int> pt_indices=ReorderByEnergy(l1Taus_pt);
  int nTaus = pt_indices.size();

  // Create tensor
  tensorflow::Tensor cellGridMatrix(tensorflow::DT_FLOAT, { nTaus, nCellEta, nCellPhi, nVars });
  initializeTensor(cellGridMatrix);

  // find objects  around l1 tau and fill matrix
  for (auto& tau_idx : pt_indices){
    // fill tensor with global observables
    for (int eta_idx = 0; eta_idx<nCellEta ; eta_idx++){
      for (int phi_idx = 0; phi_idx<nCellPhi ; phi_idx++){
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::nVertices) = static_cast<float>(patavertices.size());
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::l1Tau_pt)  = static_cast<float>(l1Taus_pt.at(tau_idx));
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::l1Tau_eta)  = static_cast<float>(l1Taus_eta.at(tau_idx));
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::l1Tau_hwIso)  = static_cast<float>(l1Taus_hwIso.at(tau_idx));
      }
    }

    // define l1tau eta and phi to evaluate dR
    float tauEta = l1Taus_eta.at(tau_idx);
    float tauPhi = l1Taus_phi.at(tau_idx);

    // fill tensor with Ecal around delta R < 0.5 wrt l1taus
    for (auto & caloRecHit_ee : *caloRecHits.ee){
      if(caloRecHit_ee.energy()<=0) continue;
      auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_ee.id())->getPosition();
      float eeCalEta = position.eta();
      float eeCalPhi = position.phi();
      float eeCalEn = caloRecHit_ee.energy();
      float eeCalChi2 = caloRecHit_ee.chi2();

      // require them to fall into deltaR < deltaR_max = 0.5 wrt l1taus
      if(DeltaR(eeCalPhi,eeCalEta,tauPhi,tauEta)<dR_max){
        float deta = DeltaEta(eeCalEta, tauEta);
        int eta_idx = static_cast<int>(floor(((deta + dR_max) / dEta_width)));
        float dphi = DeltaPhi(eeCalPhi, tauPhi);
        int phi_idx = static_cast<int>(floor(((dphi + dR_max) / dPhi_width)));
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergySum)+=static_cast<float>(eeCalEn); //
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalSize)+=static_cast<float>(1);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergyStdDev)+= static_cast<float>(eeCalEn * eeCalEn) ;
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalDeltaEta)+=static_cast<float>(deta*eeCalEn);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalDeltaPhi)+=static_cast<float>(dphi*eeCalEn);
        if(eeCalChi2>=0){
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalChi2)+=static_cast<float>(eeCalChi2*eeCalEn); //
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergySumForPositiveChi2)+=static_cast<float>(eeCalEn); //
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalSizeForPositiveChi2)+=static_cast<float>(1);
        }
      }
    } // end of loop over calorechit_ee

    for (auto & caloRecHit_eb : *caloRecHits.eb){
      if(caloRecHit_eb.energy()<=0) continue;
      auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_eb.id())->getPosition();
      float ebCalEta = position.eta();
      float ebCalPhi = position.phi();
      float ebCalEn = caloRecHit_eb.energy();
      float ebCalChi2 = caloRecHit_eb.chi2();

      // require them to fall into deltaR < deltaR_max = 0.5 wrt l1taus
      if(DeltaR(ebCalPhi,ebCalEta,tauPhi,tauEta)<dR_max){
        float deta = DeltaEta(ebCalEta, tauEta);
        int eta_idx = static_cast<int>(floor(((deta + dR_max) / dEta_width)));
        float dphi = DeltaPhi(ebCalPhi, tauPhi);
        int phi_idx = static_cast<int>(floor(((dphi + dR_max) / dPhi_width)));
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergySum)+=static_cast<float>(ebCalEn); //
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalSize)+=static_cast<float>(1);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergyStdDev)+= static_cast<float>(ebCalEn * ebCalEn) ;
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalDeltaEta)+=static_cast<float>(deta * ebCalEn);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalDeltaPhi)+=static_cast<float>(dphi * ebCalEn);
        if(ebCalChi2>=0){
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalChi2)+=static_cast<float>(ebCalChi2 * ebCalEn); //
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergySumForPositiveChi2)+=static_cast<float>(ebCalEn); //
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalSizeForPositiveChi2)+=static_cast<float>(1);
        }
      }
    } // end of loop over calorechit_eb

    // divide by the sum of all and define std dev
    for (int eta_idx = 0; eta_idx<nCellEta ; eta_idx++){
      for (int phi_idx = 0; phi_idx<nCellPhi ; phi_idx++){
        if(cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergySum)>static_cast<float>(0)){
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalDeltaEta) =  cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalDeltaEta)/cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergySum);
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalDeltaPhi) =  cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalDeltaPhi)/cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergySum);
        }
        if(cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergySumForPositiveChi2)>static_cast<float>(0)){
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalChi2) =  cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalChi2)/cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergySumForPositiveChi2);
        }
        if(cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalSize)>static_cast<float>(1)){
          float a = cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergyStdDev);
          float b = cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergySum);
          float c = cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalSize);
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergyStdDev) = ( a - (b*b/c) ) / (c-1);
        }
        else{
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::EcalEnergyStdDev) =static_cast<float>(0);
        }
      }
    }


    // fill tensor with Hcal around delta R < 0.5 wrt l1taus
    for (auto & caloRecHit_hbhe : *caloRecHits.hbhe){
      if(caloRecHit_hbhe.energy()<=0) continue;
      auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_hbhe.id())->getPosition();
      float hbheCalEta = position.eta();
      float hbheCalPhi = position.phi();
      float hbheCalEn = caloRecHit_hbhe.energy();
      float hbheCalChi2 = caloRecHit_hbhe.chi2();

      // require them to fall into deltaR < deltaR_max = 0.5 wrt l1taus
      if(DeltaR(hbheCalPhi,hbheCalEta,tauPhi,tauEta)<dR_max){
        float deta = DeltaEta(hbheCalEta, tauEta);
        int eta_idx = static_cast<int>(floor(((deta + dR_max) / dEta_width)));
        float dphi = DeltaPhi(hbheCalPhi, tauPhi);
        int phi_idx = static_cast<int>(floor(((dphi + dR_max) / dPhi_width)));
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergySum)+=static_cast<float>(hbheCalEn); //
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergyStdDev)+=static_cast<float>(hbheCalEn*hbheCalEn); //
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalSize)+=static_cast<float>(1);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalDeltaEta)+=static_cast<float>(deta*hbheCalEn);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalDeltaPhi)+=static_cast<float>(dphi*hbheCalEn);
        if(hbheCalChi2>=0){
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalChi2)+=static_cast<float>(hbheCalChi2*hbheCalEn); //
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergySumForPositiveChi2)+=static_cast<float>(hbheCalEn); //
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalSizeForPositiveChi2)+=static_cast<float>(1);
        }
      }
    } // end of loop over calorechit_hbhe

    for (auto & caloRecHit_ho : *caloRecHits.ho){
      if(caloRecHit_ho.energy()<=0) continue;
      auto& position = caloRecHits.Geometry->getGeometry(caloRecHit_ho.id())->getPosition();
      float hoCalEta = position.eta();
      float hoCalPhi = position.phi();
      float hoCalEn = caloRecHit_ho.energy();

      // require them to fall into deltaR < deltaR_max = 0.5 wrt l1taus
      if(DeltaR(hoCalPhi,hoCalEta,tauPhi,tauEta)<dR_max){
        float deta = DeltaEta(hoCalEta, tauEta);
        int eta_idx = static_cast<int>(floor(((deta + dR_max) / dEta_width)));
        float dphi = DeltaPhi(hoCalPhi, tauPhi);
        int phi_idx = static_cast<int>(floor(((dphi + dR_max) / dPhi_width)));
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergySum)+=static_cast<float>(hoCalEn); //
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergyStdDev)+=static_cast<float>(hoCalEn*hoCalEn); //
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalSize)+=static_cast<float>(1);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalDeltaEta)+=static_cast<float>(deta*hoCalEn);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalDeltaPhi)+=static_cast<float>(dphi*hoCalEn);
      }
    } // end of loop over calorechit_ho

    for (int eta_idx = 0; eta_idx<nCellEta ; eta_idx++){
      for (int phi_idx = 0; phi_idx<nCellPhi ; phi_idx++){
        if(cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergySum)>static_cast<float>(0)){
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalDeltaEta) =  cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalDeltaEta)/cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergySum);
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalDeltaPhi) =  cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalDeltaPhi)/cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergySum);
        }
        if(cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergySumForPositiveChi2)>static_cast<float>(0)){
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalChi2) =  cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalChi2)/cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergySumForPositiveChi2);
        }
        if(cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalSize)>static_cast<float>(1)){
          float a = cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergyStdDev);
          float b = cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergySum);
          float c = cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalSize);
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergyStdDev) = ( a - (b*b/c) ) / (c-1);
        }
        else{
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::HcalEnergyStdDev) =static_cast<float>(0);
        }
      }
    }

    // fill tensor with Patatracks around delta R < 0.5 wrt l1taus
    for (unsigned n = 0; n < patatracks.size(); n++){
      if(patatracks.at(n).pt()<=0) continue;
      float patatrackEta = patatracks.at(n).eta();
      float patatrackPhi = patatracks.at(n).phi();
      float patatrackPt = patatracks.at(n).pt();
      float patatrackNdof =patatracks.at(n).ndof();
      float patatrackChi2OverNdof ;
      if(patatracks.at(n).ndof() == 0.){
        patatrackChi2OverNdof = 0.;
      }
      else{
          patatrackChi2OverNdof = patatracks.at(n).chi2()/(patatracks.at(n).ndof());
      }
      float patatrackDxy = patatracks.at(n).dxy();
      float patatrackDz = patatracks.at(n).dz();

      // require them to fall into deltaR < deltaR_max = 0.5 wrt l1taus
      if(DeltaR(patatrackPhi,patatrackEta,tauPhi,tauEta)<dR_max){
        float deta = DeltaEta(patatrackEta, tauEta);
        int eta_idx = static_cast<int>(floor(((deta + dR_max) / dEta_width)));
        float dphi = DeltaPhi(patatrackPhi, tauPhi);
        int phi_idx = static_cast<int>(floor(((dphi + dR_max) / dPhi_width)));
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackPtSum) +=static_cast<float>(patatrackPt);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackSize) +=static_cast<float>(1);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackChargeSum) +=static_cast<float>(patatracks.at(n).charge());
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDeltaEta) +=static_cast<float>(deta*patatrackPt);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDeltaPhi) +=static_cast<float>(dphi*patatrackPt);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackChi2OverNdof) +=static_cast<float>(patatrackChi2OverNdof*patatrackPt);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackNdof) +=static_cast<float>(patatrackNdof*patatrackPt);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDxy) +=static_cast<float>(patatrackDxy*patatrackPt);
        cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackDz) +=static_cast<float>(patatrackDz*patatrackPt);
        if(FindVertexIndex(patavertices, patatracks.at(n) ) != -1){
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackPtSumWithVertex) +=static_cast<float>(patatrackPt);
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx, NNInputs::PatatrackSizeWithVertex) +=static_cast<float>(1);
        }
      }

    } // end of loop over patatracks
    for (int eta_idx = 0; eta_idx<nCellEta ; eta_idx++){
      for (int phi_idx = 0; phi_idx<nCellPhi ; phi_idx++){
        if(cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackPtSum)>static_cast<float>(0)){
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackDeltaEta) =  cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackDeltaEta)/cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackPtSum);
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackDeltaPhi) =  cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackDeltaPhi)/cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackPtSum) ;
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackChi2OverNdof) =  cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackChi2OverNdof)/cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackPtSum) ;
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackNdof) =  cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackNdof)/cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackPtSum) ;
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackDxy) =  cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackDxy)/cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackPtSum) ;
          cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackDz) =  cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackDz)/cellGridMatrix.tensor<float, 4>()(tau_idx, phi_idx, eta_idx,NNInputs::PatatrackPtSum) ;
        }
      }
    }
  } // end of loop over taus

  std::vector<int> cellGridMatrixShape = get_tensor_shape(cellGridMatrix);
  bool printoutevt=false;
  checknan(cellGridMatrix,printoutevt);
  standardizeTensor(cellGridMatrix, normalizationDict_);

  /* for debugging uncomment  following lines */
  //if(evt_id==135012){
  //  printoutevt=true;
  //}
  //checknan(cellGridMatrix,printoutevt);

  std::vector<tensorflow::Tensor> pred_vector;
  tensorflow::run(session_, {{inputTensorName_, cellGridMatrix}}, {outputTensorName_}, &pred_vector);

  for (auto& tau_idx : pt_indices){
    std::cout << "evt == " << evt_id<<" outcome for tau with pt =  " << l1Taus_pt.at(tau_idx) << " is \t "<< pred_vector[0].matrix<float>()(tau_idx, 0)<< std::endl;

  }
  int n_TauPassed =0;
  /* for debugging uncomment  following lines */
  //for (auto& tau_idx : pt_indices){
  //  if(l1Taus_pt.at(tau_idx)>=32. && (l1Taus_pt.at(tau_idx)>=70. || l1Taus_hwIso.at(tau_idx)>0) && pred_vector[0].matrix<float>()(tau_idx, 0)>static_cast<float>(discr_threshold_))
  //  {
  //    n_TauPassed+=1;
  //  }
  //}
  if(n_TauPassed==2){
    result = true;
  }
}
bool L2TauNNTag::filter(edm::Event& event, const edm::EventSetup& eventsetup) {
  bool result = false;
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
  edm::Handle<reco::TrackCollection> pataTracks;
  event.getByToken(pataTracks_token, pataTracks);
  edm::Handle<reco::VertexCollection> pataVertices;
  event.getByToken(pataVertices_token, pataVertices);
  edm::ESHandle<CaloGeometry> Geometry = eventsetup.getHandle(Geometry_token);
  caloRecHitCollections caloRecHits;
  caloRecHits.hbhe= &*hbhe;
  caloRecHits.ho= &*ho;
  caloRecHits.eb= &*ebHandle;
  caloRecHits.ee= &*eeHandle;
  caloRecHits.Geometry = &*Geometry;
  int evt_id = event.id().event();

  FindObjectsAroundL1Tau(caloRecHits, *pataTracks, *pataVertices,*l1Taus,evt_id,result);
  /* uncomment following lines for debugging */
  if (result == true){
    std::cout << "event passed == " << evt_id << std::endl;
  }

  return result;
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L2TauNNTag);
