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
        processName = cms.string('MLProva'),
        l1taus=cms.InputTag('hltGtStage2Digis','Tau'),
        ecalInputs =cms.VInputTag("hltEcalRecHit:EcalRecHitsEB", "hltEcalRecHit:EcalRecHitsEE"),
        hbheInput = cms.InputTag("hltHbhereco"),
        hoInput = cms.InputTag("hltHoreco"),
        pataVertices = cms.InputTag("hltTrimmedPixelVertices"),
        pataTracks = cms.InputTag("hltPixelTracks"),
        graphPath = cms.string(graphPath),
        normalizationDict = cms.string(normalizationDict),
        discr_threshold = cms.double(thWp[rateWP])
    )


    #process.l2TauNNTag.discr_threshold = cms.double(thWp[rateWP])
    #process.MLPathTest = cms.Path(process.HLTBeginSequence + process.hltL1sDoubleTauBigOR , process.HLTDoLocalPixelTask, process.HLTRecoPixelTracksTask, process.HLTRecopixelvertexingTask)
    process.MLPathTest = cms.Path(process.HLTBeginSequence + process.hltL1sDoubleTauBigOR + process.HLTDoCaloSequence +process.l2TauNNTagFilter, process.HLTDoLocalPixelTask, process.HLTRecoPixelTracksTask, process.HLTRecopixelvertexingTask)

    process.schedule = cms.Schedule(*[ process.HLTriggerFirstPath, process.MLPathTest, process.HLTriggerFinalPath, process.endjob_step ], tasks=[process.patAlgosToolsTask])

    process.options.wantSummary = cms.untracked.bool(False)
    return process
