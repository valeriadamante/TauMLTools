/*! Definition of a tuple with all the l2 level information. */

#pragma once

#include "TauMLTools/Core/interface/SmartTree.h"
#include <Math/VectorUtil.h>

#define HVAR(type, name) VAR(std::vector<type>, caloRecHit_had_##name)
#define EVAR(type, name) VAR(std::vector<type>, caloRecHit_e_##name)
#define PVAR(type, name) VAR(std::vector<type>, patatrack_##name)
#define VAR2(type, name1, name2) VAR(type, name1) VAR(type, name2)
#define VAR3(type, name1, name2, name3) VAR2(type, name1, name2) VAR(type, name3)
#define VAR4(type, name1, name2, name3, name4) VAR3(type, name1, name2, name3) VAR(type, name4)

#define INTERMEDIATE_DATA() \
    /* Event Variables */ \
    VAR(UInt_t, run) /* run number */ \
    VAR(UInt_t, lumi) /* lumi section */ \
    VAR(ULong64_t, evt) /* event number */ \
    VAR(Int_t, nVertices)/*e*/\
    /* Gen lepton with the full decay chain */ \
    VAR(Bool_t, genLepton_isTau) /* genlep index candidate*/ \
    VAR(Int_t, genLepton_charge) /* charge of the gen lepton */ \
    VAR4(Float_t, genLepton_vis_pt, genLepton_vis_eta, genLepton_vis_phi, genLepton_vis_mass) /* visible 4-momentum of
                                                                                                 the gen lepton */ \
    /* L1 objects */ \
    VAR(Int_t, l1Tau_index) /* L1 index candidate*/ \
    VAR(Float_t, l1Tau_pt) /* L1 pt candidate ! */ \
    VAR(Float_t, l1Tau_eta) /* L1 eta candidate ! */ \
    VAR(Float_t, l1Tau_phi) /* L1 phi candidate ! */ \
    VAR(Float_t, l1Tau_mass) /* L1 mass candidate ?! */ \
    VAR(int, l1Tau_hwIso) /* L1 hwIso candidate ! */ \
    /* CaloRecHits candidates */ \
    EVAR(Bool_t, isEndCap) /**/ \
    EVAR(Float_t, rho) /* ! */ \
    EVAR(Float_t, DeltaEta) /* Delta ! */ \
    EVAR(Float_t, DeltaPhi) /* Delta !  */ \
    EVAR(Float_t, DeltaR) /* Delta !  */ \
    EVAR(Float_t, energy) /* ! */ \
    EVAR(Float_t, chi2) /* ! */ \
    EVAR(Bool_t, isRecovered) /* ved */ \
    HVAR(Bool_t, HadronSubDet) /*e*/ \
    HVAR(Float_t, rho) /* d*/ \
    HVAR(Float_t, DeltaEta) /* */ \
    HVAR(Float_t, DeltaPhi) /* */ \
    HVAR(Float_t, DeltaR) /* */ \
    HVAR(Float_t, energy) /* */ \
    HVAR(Float_t, chi2) /* */ \
    HVAR(Float_t, time) /* */ \
    /* Tracks candidates */ \
    PVAR(Float_t, pt) /* track pt candidate*/ \
    PVAR(Float_t, eta) /* track eta candidate*/ \
    PVAR(Float_t, phi) /* track phi candidate*/ \
    PVAR(Float_t, DeltaEta) /* Delta eta with l1tau*/ \
    PVAR(Float_t, DeltaPhi) /* Delta phi with l1tau*/ \
    PVAR(Float_t, DeltaR) /* Delta R with l1tau*/ \
    PVAR(Float_t, chi2) /* track chi2 candidate*/ \
    PVAR(Int_t, ndof) /* track ndof candidate*/ \
    PVAR(Int_t, charge) /* pixelTrack charge candidate*/ \
    PVAR(UInt_t, quality) /* pixelTrack qualityMask candidate*/ \
    PVAR(Float_t, dxy) /* track dxy candidate*/ \
    PVAR(Float_t, dz) /* track dz candidate*/ \
    PVAR(Bool_t, hasVertex) /*t*/\
    PVAR(Float_t, vert_z) /* x positions of vertices */ \
    PVAR(Float_t, vert_weight) /* output weight (1/error^2) on the above */ \
    PVAR(Float_t, vert_ptv2) /* vertices pt^2 */ \
    PVAR(Float_t, vert_chi2) /* chi^2 of the vertices (PV) */ \
    PVAR(Int_t, vert_ndof) /* number of degrees of freedom of vertices (PV) */ \

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(interm_tuple, IntermTau, TauTrainIntermediateTuple, INTERMEDIATE_DATA, "L2TauTrainTuple")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(interm_tuple, TauTrainIntermediateTuple, INTERMEDIATE_DATA)
#undef VAR
#undef VAR2
#undef VAR3
#undef VAR4
#undef PVAR
#undef EVAR
#undef HVAR
