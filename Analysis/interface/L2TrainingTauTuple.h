/*! Definition of a tuple with all the l2 level information. */

#pragma once

#include "TauMLTools/Core/interface/SmartTree.h"
#include <Math/VectorUtil.h>

#ifndef HBHE_VAR
  #define HBHE_VAR(type, name) VAR(std::vector<type>, caloRecHit_hbhe_##name)
#endif
#ifndef HO_VAR
  #define HO_VAR(type, name) VAR(std::vector<type>, caloRecHit_ho_##name)
#endif
#ifndef HF_VAR
  #define HF_VAR(type, name) VAR(std::vector<type>, caloRecHit_hf_##name)
#endif
#ifndef ECAL_VAR
  #define ECAL_VAR(type, name) VAR(std::vector<type>, caloRecHit_ee_##name) VAR(std::vector<type, caloRecHit_eb_##name)
#endif
#define PATATRACK_VAR(type, name) VAR(std::vector<type>, patatrack_##name)
#define PATAVERT(type, name) VAR(std::vector<type>, patavert_##name)



#define VAR2(type, name1, name2) VAR(type, name1) VAR(type, name2)
#define VAR3(type, name1, name2, name3) VAR2(type, name1, name2) VAR(type, name3)
#define VAR4(type, name1, name2, name3, name4) VAR3(type, name1, name2, name3) VAR(type, name4)

#define L2TRAINTAU_DATA() \
    /* Event Variables */ \
    VAR(UInt_t, run) /* run number */ \
    VAR(UInt_t, lumi) /* lumi section */ \
    VAR(ULong64_t, evt) /* event number */ \
    VAR(Float_t, genEventWeight) /* gen event weight */ \
    VAR(Int_t, sampleType) /* type of the sample (MC, Embedded or Data) */ \
    /* Gen lepton with the full decay chain */ \
    VAR(Int_t, genLepton_kind) /* kind of the gen lepton:
                                              Electron = 1, Muon = 2, TauElectron = 3, TauMuon = 4, Tau = 5, Other = 6 */\
    VAR(Int_t, genLepton_charge) /* charge of the gen lepton */ \
    VAR4(Float_t, genLepton_vis_pt, genLepton_vis_eta, genLepton_vis_phi, genLepton_vis_mass) /* visible 4-momentum of
                                                                                                 the gen lepton */ \
    /* L1 objects */ \
    VAR(Float_t, l1Tau_pt) /* L1 pt candidate*/ \
    VAR(Float_t, l1Tau_eta) /* L1 eta candidate*/ \
    VAR(Float_t, l1Tau_phi) /* L1 phi candidate*/ \
    VAR(Float_t, l1Tau_mass) /* L1 mass candidate*/ \
    VAR(int, l1Tau_hwIso) /* L1 hwIso candidate*/ \
    VAR(int, l1Tau_hwQual) /* L1 quality candidate*/ \
    VAR(int, l1Tau_towerIEta) /* L1 towerIEta candidate*/ \
    VAR(int, l1Tau_towerIPhi) /* L1 towerIPhi candidate*/ \
    VAR(int, l1Tau_rawEt) /* L1 rawEt candidate*/ \
    VAR(int, l1Tau_isoEt) /* L1 isoEt candidate*/ \
    VAR(bool, l1Tau_hasEM) /* L1 hasEM candidate*/ \
    VAR(bool, l1Tau_isMerged) /* L1 isMerged candidate*/ \
    /* CaloRecHits candidates */ \
    ECAL_VAR(Float_t, rho) /* */ \
    ECAL_VAR(Float_t, eta) /* */ \
    ECAL_VAR(Float_t, phi) /* */ \
    ECAL_VAR(Float_t, energy) /* */ \
    ECAL_VAR(Float_t, time) /* */ \
    ECAL_VAR(ULong64_t, detId) /* */ \
    ECAL_VAR(Float_t, chi2) /* */ \
    ECAL_VAR(Float_t, energyError) /* */ \
    ECAL_VAR(Float_t, timeError) /* */ \
    ECAL_VAR(uint32_t, flagsBits) /* */ \
    ECAL_VAR(Bool_t, isRecovered) /* */ \
    ECAL_VAR(Bool_t, isTimeValid) /* */ \
    ECAL_VAR(Bool_t, isTimeErrorValid) /* */ \
    HBHE_VAR(Float_t, rho) /* */ \
    HBHE_VAR(Float_t, eta) /* */ \
    HBHE_VAR(Float_t, phi) /* */ \
    HBHE_VAR(Float_t, energy) /* */ \
    HBHE_VAR(Float_t, time) /* */ \
    HBHE_VAR(ULong64_t, detId) /* */ \
    HBHE_VAR(Float_t, chi2) /* */ \
    HBHE_VAR(ULong64_t, flags) /* */ \
    HBHE_VAR(Float_t, eraw) /* */ \
    HBHE_VAR(Float_t, eaux) /* */ \
    HBHE_VAR(Float_t, timeFalling) /* */ \
    HBHE_VAR(ULong64_t, idFront) /* */ \
    HBHE_VAR(Float_t, rho_front) /* */ \
    HBHE_VAR(Float_t, eta_front) /* */ \
    HBHE_VAR(Float_t, phi_front) /* */ \
    HBHE_VAR(UInt_t, auxHBHE) /* */ \
    HBHE_VAR(UInt_t, auxPhase1) /* */ \
    HBHE_VAR(UInt_t, auxTDC) /* */ \
    HBHE_VAR(Bool_t, isMerged) /* */ \
    HO_VAR(Float_t, rho) /* */ \
    HO_VAR(Float_t, eta) /* */ \
    HO_VAR(Float_t, phi) /* */ \
    HO_VAR(Float_t, energy) /* */ \
    HO_VAR(Float_t, time) /* */ \
    HO_VAR(ULong64_t, detId) /* */ \
    HO_VAR(ULong64_t, aux) /* */ \
    HO_VAR(ULong64_t, flags) /* */ \
    HF_VAR(Float_t, rho) /* */ \
    HF_VAR(Float_t, eta) /* */ \
    HF_VAR(Float_t, phi) /* */ \
    HF_VAR(Float_t, energy) /* */ \
    HF_VAR(Float_t, time) /* */ \
    HF_VAR(ULong64_t, detId) /* */ \
    HF_VAR(ULong64_t, flags) /* */ \
    HF_VAR(Float_t, timeFalling) /* */ \
    HF_VAR(uint32_t, auxHF) /* */ \
    HF_VAR(ULong64_t, aux) /* */ \
    /* Tracks candidates */ \
    PATATRACK_VAR(Float_t, pt) /* track pt candidate*/ \
    PATATRACK_VAR(Float_t, eta) /* track eta candidate*/ \
    PATATRACK_VAR(Float_t, phi) /* track phi candidate*/ \
    PATATRACK_VAR(Float_t, chi2) /* track chi2 candidate*/ \
    PATATRACK_VAR(Int_t, ndof) /* track ndof candidate*/ \
    PATATRACK_VAR(Int_t, charge) /* pixelTrack charge candidate*/ \
    PATATRACK_VAR(UInt_t, quality) /* pixelTrack qualityMask candidate*/ \
    PATATRACK_VAR(Float_t, dxy) /* track dxy candidate*/ \
    PATATRACK_VAR(Float_t, dz) /* track dz candidate*/ \
    PATATRACK_VAR(Int_t, vertex_id) /* track associated vertex id candidate*/ \
    /* PATAVERTICES */ \
    PATAVERT(Float_t, z) /* x positions of vertices */ \
    PATAVERT(Float_t, weight) /* output weight (1/error^2) on the above */ \
    PATAVERT(Float_t, ptv2) /* vertices pt^2 */ \
    PATAVERT(Float_t, chi2) /* chi^2 of the vertices (PV) */ \
    PATAVERT(Int_t, ndof) /* number of degrees of freedom of vertices (PV) */ \

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(tau_train_tuple, Tau, TauTrainTuple, L2TRAINTAU_DATA, "TauTrainTuple")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(tau_train_tuple, TauTrainTuple, L2TRAINTAU_DATA)
#undef VAR
#undef VAR2
#undef VAR3
#undef VAR4
#undef TAU_ID
#undef TAU_VAR
#undef CALO_TOWER_VAR
#undef CALO_TAU_VAR
#undef PATATRACK_VAR
#undef PATAVERT
