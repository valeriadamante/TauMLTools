#ifndef TauMLTools_Production_L2TauTagNNFilter_h
#define TauMLTools_Production_L2TauTagNNFilter_h
// system include files
#include <memory>
#include <Math/VectorUtil.h>
#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/variadic.hpp>
#include <boost/math/constants/constants.hpp>
// user include files
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
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

#include "Compression.h"


#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/StringToEnumValue.h"
#include "RecoTauTag/RecoTau/interface/PFRecoTauClusterVariables.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "TauMLTools/Analysis/interface/GenLepton.h"
#include "TauMLTools/Analysis/interface/SummaryTuple.h"
#include "TauMLTools/Analysis/interface/TauIdResults.h"
#include "TauMLTools/Analysis/interface/L2EventTuple.h"
#include "TauMLTools/Core/interface/Tools.h"
#include "TauMLTools/Core/interface/TextIO.h"
#include "TauMLTools/Production/interface/GenTruthTools.h"
#include "TauMLTools/Production/interface/MuonHitMatch.h"
#include "TauMLTools/Production/interface/TauAnalysis.h"
#include "TauMLTools/Production/interface/TauJet.h"

#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTowerTopology.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "RecoLocalCalo/CaloTowersCreator/interface/EScales.h"
#include "RecoLocalCalo/CaloTowersCreator/src/CaloTowersCreator.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
// forward declarations
struct caloRecHitCollections {
  const HBHERecHitCollection  *hbhe;
  const HORecHitCollection *ho;
  const HFRecHitCollection *hf;
  const EcalRecHitCollection *eb;
  const EcalRecHitCollection *ee;
  const CaloGeometry *Geometry;
};

class L2TauNNTag : public edm::EDFilter{
public:
  explicit L2TauNNTag(const edm::ParameterSet& cfg);
  ~L2TauNNTag() override;

private:
  void beginJob() override;
  bool filter(edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  VOID fillDescriptions(edm:: ConfigurationDescriptions& descriptions);
  static constexpr float pi = boost::math::constants::pi<float>();
  float DeltaPhi(Float_t phi1, Float_t phi2);
  float DeltaEta(Float_t eta1, Float_t eta2);
  float DeltaR(Float_t phi1,Float_t eta1,Float_t phi2,Float_t eta2);
  std::vector<int> ReorderByEnergy(std::vector<float>& energy);
  int FindVertexIndex(const reco::VertexCollection& vertices, const reco::Track& track);
  void FindObjectsAroundL1Tau(const caloRecHitCollections& caloRecHits, const reco::TrackCollection& patatracks,const reco::VertexCollection& patavertices, const l1t::TauBxCollection& l1Taus);

private:
  std::string processName;
  edm::EDGetTokenT<l1t::TauBxCollection> l1Taus_token;
  edm::EDGetTokenT<HBHERecHitCollection> hbhe_token;
  edm::EDGetTokenT<HORecHitCollection> ho_token;
  edm::EDGetTokenT<HFRecHitCollection> hf_token;
  std::vector<edm::InputTag> ecalLabels;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> Geometry_token;
  edm::EDGetTokenT<reco::VertexCollection> pataVertices_token;
  edm::EDGetTokenT<reco::TrackCollection> pataTracks_token;
  std::vector<edm::EDGetTokenT<EcalRecHitCollection>> ecal_tokens;
  std::string graphPath_;
  std::string inputTensorName_;
  std::string outputTensorName_;
  tensorflow::GraphDef* graphDef_;
  tensorflow::Session* session_;
};




#endif
