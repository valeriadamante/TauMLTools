#ifndef TauMLTools_Production_L2TauTagNNFilter_h
#define TauMLTools_Production_L2TauTagNNFilter_h
// system include files
#include <memory>

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
// L1 Tau
#include "DataFormats/L1Trigger/interface/Tau.h"
// forward declarations
enum class NNInputs { nVertices = 0,
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

class L2TauNNTag : public edm::EDFilter {
public:
  explicit L2TauNNTag(const edm::ParameterSet& cfg);
  ~L2TauNNTag() override;

private:
  void beginJob() override;

  bool filter(const edm::Event& event, const edm::EventSetup& eventsetup) override;
  void endJob() override;

  bool isPassingL2NNTag(const unsigned int nTaus) const;

  static constexpr float pi = boost::math::constants::pi<float>();
  template<typename Scalar>
  static Scalar DeltaPhi(Scalar phi1, Scalar phi2);
  template<typename Scalar>
  static Scalar DeltaEta(Scalar eta1, Scalar eta2);
  template<typename Scalar>
  static Scalar DeltaR(Scalar phi1,Scalar eta1,Scalar phi2,Scalar eta2);
  std::vector<int> ReorderByEnergy(std::vector<float>& energy);
  int FindVertexIndex(const reco::VertexCollection& vertices, const reco::Track& track);
  void FindObjectsAroundL1Tau(const AllCaloRecHits& caloRecHits, const reco::TrackCollection& patatracks,const reco::VertexCollection& patavertices, const l1t::TauBxCollection& l1Taus);

  struct caloRecHitCollections{
    const HBHERecHitCollection  *hbhe;
    const HORecHitCollection *ho;
    const HFRecHitCollection *hf;
    const EcalRecHitCollection *eb;
    const EcalRecHitCollection *ee;
    const CaloGeometry *Geometry;
  };
  std::string processName;
  //std::string defaultDiTauPath;
  //edm::EDGetTokenT<edm::TriggerResults> TR_token;
  //edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfo_token;
  edm::EDGetTokenT<l1t::TauBxCollection> l1Taus_token;
  edm::EDGetTokenT<HBHERecHitCollection> hbhe_token;
  edm::EDGetTokenT<HORecHitCollection> ho_token;
  edm::EDGetTokenT<HFRecHitCollection> hf_token;
  std::vector<edm::InputTag> ecalLabels;
  std::vector<edm::EDGetTokenT<EcalRecHitCollection>> ecal_tokens;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> Geometry_token;
  edm::EDGetTokenT<reco::VertexCollection> pataVertices_token;
  edm::EDGetTokenT<reco::TrackCollection> pataTracks_token;

  };


#endif
