import FWCore.ParameterSet.Config as cms
from HLTrigger.Configuration.customizeHLTforPatatrack import customizeHLTforPatatrackTriplets
from TauMLTools.Production.l2TauNNTag_cfi import *

def update(process):
    thWp = {
        'opt_threshold_3': 0.180858813224404,
        'opt_threshold_4': 0.12267940863785043,
        'opt_threshold_5': 0.08411243185219064,
    }
    graphPath = '/afs/cern.ch/work/v/vdamante/public/CMSSW_11_2_1_Patatrack/src/TauMLTools/Analysis/config/graph_model/graph_saved_model.pb'
    normalizationDict = '/afs/cern.ch/work/v/vdamante/public/CMSSW_11_2_1_Patatrack/src/TauMLTools/Analysis/config/NormalizationDict.json'
    rateValue = 3
    rateWP=("opt_threshold_{}").format(rateValue)


    process = customizeHLTforPatatrackTriplets(process)
    process.l2TauNNTagFilter = l2TauNNTag.clone(
        debugLevel = 0,
        processName = cms.string('MLProva'),
        L1TauTrigger=cms.InputTag('hltL1sDoubleTauBigOR'),
        #l1taus=cms.InputTag('hltGtStage2Digis','Tau'),
        ecalInputs =cms.VInputTag("hltEcalRecHit:EcalRecHitsEB", "hltEcalRecHit:EcalRecHitsEE"),
        hbheInput = cms.InputTag("hltHbhereco"),
        hoInput = cms.InputTag("hltHoreco"),
        pataVertices = cms.InputTag("hltTrimmedPixelVertices"),
        pataTracks = cms.InputTag("hltPixelTracks"),
        graphPath = cms.string(graphPath),
        normalizationDict = cms.string(normalizationDict),
        discr_threshold = cms.double(thWp[rateWP])
    )
    # L2 updated Sequence
    process.L2Sequence = cms.Sequence(process.HLTDoCaloSequence +  process.l2TauNNTagFilter)

    # Regional -> Global customization
    process.hltHpsPFTauTrackPt1DiscriminatorReg.PFTauProducer = cms.InputTag("hltHpsPFTauProducer")
    process.hltHpsDoublePFTau35Reg.inputTag = cms.InputTag( "hltHpsPFTauProducer")
    process.hltHpsSelectedPFTausTrackPt1Reg.src = cms.InputTag( "hltHpsPFTauProducer")
    process.hltHpsPFTauMediumAbsoluteChargedIsolationDiscriminatorReg.PFTauProducer = cms.InputTag( "hltHpsPFTauProducer" )
    process.hltHpsPFTauMediumAbsoluteChargedIsolationDiscriminatorReg.particleFlowSrc = cms.InputTag( "hltParticleFlow" )
    process.hltHpsPFTauMediumRelativeChargedIsolationDiscriminatorReg.PFTauProducer = cms.InputTag( "hltHpsPFTauProducer" )
    process.hltHpsPFTauMediumRelativeChargedIsolationDiscriminatorReg.particleFlowSrc = cms.InputTag( "hltParticleFlow" )
    process.hltHpsPFTauMediumAbsOrRelChargedIsolationDiscriminatorReg.PFTauProducer = cms.InputTag( "hltHpsPFTauProducer" )
    process.hltHpsSelectedPFTausTrackPt1MediumChargedIsolationReg.src = cms.InputTag( "hltHpsPFTauProducer" )

    # re-define path with l2 updated sequence
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4 = cms.Path(process.HLTBeginSequence + process.hltL1sDoubleTauBigOR +
    process.hltPreDoubleMediumChargedIsoPFTauHPS35Trk1eta2p1Reg +  process.L2Sequence + process.HLTGlobalPFTauHPSSequence +
    process.HLTHPSDoublePFTauPt35Eta2p1Trk1Reg + process.HLTHPSMediumChargedIsoPFTauSequenceReg +
    process.hltHpsSelectedPFTausTrackPt1MediumChargedIsolationReg + process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationReg +
    process.hltHpsL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationMatchReg +
    process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationL1HLTMatchedReg + process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg +
    process.HLTEndSequence, process.HLTDoLocalPixelTask, process.HLTRecoPixelTracksTask, process.HLTRecopixelvertexingTask)


    process.schedule = cms.Schedule(*[ process.HLTriggerFirstPath, process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4,
    process.HLTriggerFinalPath, process.endjob_step ], tasks=[process.patAlgosToolsTask])

    return process
