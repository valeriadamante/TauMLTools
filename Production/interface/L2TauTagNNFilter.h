#ifndef TauMLTools_Production_L2TauTagNNFilter_h
#define TauMLTools_Production_L2TauTagNNFilter_h
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
#include "FWCore/Framework/interface/EDFilter.h"
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


// forward declarations
struct caloRecHitCollections {
  const HBHERecHitCollection  *hbhe;
  const HORecHitCollection *ho;
  const EcalRecHitCollection *eb;
  const EcalRecHitCollection *ee;
  const CaloGeometry *Geometry;
};

class L2TauNNTag : public edm::EDFilter{
public:
  explicit L2TauNNTag(const edm::ParameterSet& cfg);
  //static void fillDescriptions(edm::ConfigurationDescriptions&);
  ~L2TauNNTag() override;

private:
  void beginJob() override;
  bool filter(edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  static constexpr float pi = boost::math::constants::pi<float>();
  float DeltaPhi(Float_t phi1, Float_t phi2);
  float DeltaEta(Float_t eta1, Float_t eta2);
  float DeltaR(Float_t phi1,Float_t eta1,Float_t phi2,Float_t eta2);
  std::vector<int> get_tensor_shape(tensorflow::Tensor& tensor);
  std::vector<int> ReorderByEnergy(std::vector<float>& energy);
  void initializeTensor(tensorflow::Tensor& tensor);
  void checknan(tensorflow::Tensor& tensor, bool printoutTensor);
  void standardizeTensor(tensorflow::Tensor& tensor, std::string normDictPath);
  int FindVertexIndex(const reco::VertexCollection& vertices, const reco::Track& track);
  void FindObjectsAroundL1Tau(const caloRecHitCollections& caloRecHits, const reco::TrackCollection& patatracks,const reco::VertexCollection& patavertices, const l1t::TauBxCollection& l1Taus, int evt_id);

private:
  std::string processName;
  edm::EDGetTokenT<l1t::TauBxCollection> l1Taus_token;
  edm::EDGetTokenT<HBHERecHitCollection> hbhe_token;
  edm::EDGetTokenT<HORecHitCollection> ho_token;
  std::vector<edm::InputTag> ecalLabels;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> Geometry_token;
  edm::EDGetTokenT<reco::VertexCollection> pataVertices_token;
  edm::EDGetTokenT<reco::TrackCollection> pataTracks_token;
  std::vector<edm::EDGetTokenT<EcalRecHitCollection>> ecal_tokens;
  std::string graphPath_;
  std::string inputTensorName_;
  std::string outputTensorName_;
  std::string normalizationDict_;
  std::unique_ptr<tensorflow::GraphDef> graphDef_;
  std::unique_ptr<tensorflow::Session> session_;
  //tensorflow::GraphDef* graphDef_;
  //tensorflow::Session* session_;
};




#endif
