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

    process.L2Sequence = cms.Sequence(process.HLTDoCaloSequence +  process.l2TauNNTagFilter)

    process.hltHpsPFTauProducer = cms.EDProducer( "RecoTauPiZeroUnembedder",
        src = cms.InputTag( "hltHpsPFTauProducerSansRefs" )
    )
    process.hltHpsDoublePFTau35 = cms.EDFilter( "HLT1PFTau",
        saveTags = cms.bool( True ),
        MinPt = cms.double( 35.0 ),
        MinN = cms.int32( 2 ),
        MaxEta = cms.double( 2.1 ),
        MinEta = cms.double( -1.0 ),
        MinMass = cms.double( -1.0 ),
        inputTag = cms.InputTag( "hltHpsPFTauProducer" ),
        MinE = cms.double( -1.0 ),
        triggerType = cms.int32( 84 ),
        MaxMass = cms.double( -1.0 )
    )

    process.hltHpsPFTauTrackPt1Discriminator = cms.EDProducer( "PFRecoTauDiscriminationByLeadingObjectPtCut",
        MinPtLeadingObject = cms.double( 1.0 ),
        PFTauProducer = cms.InputTag( "hltHpsPFTauProducer" ),
        Prediscriminants = cms.PSet(  BooleanOperator = cms.string( "and" ) ),
        UseOnlyChargedHadrons = cms.bool( True )
    )
    process.hltHpsSelectedPFTausTrackPt1 = cms.EDFilter( "PFTauSelector",
        discriminators = cms.VPSet(
          cms.PSet(  discriminator = cms.InputTag( "hltHpsPFTauTrackPt1Discriminator" ),
            selectionCut = cms.double( 0.5 )
          )
        ),
        discriminatorContainers = cms.VPSet(
        ),
        cut = cms.string( "pt > 0" ),
        src = cms.InputTag( "hltHpsPFTauProducer" )
    )
    process.hltHpsDoublePFTau35TrackPt1 = cms.EDFilter( "HLT1PFTau",
        saveTags = cms.bool( True ),
        MinPt = cms.double( 35.0 ),
        MinN = cms.int32( 2 ),
        MaxEta = cms.double( 2.1 ),
        MinEta = cms.double( -1.0 ),
        MinMass = cms.double( -1.0 ),
        inputTag = cms.InputTag( "hltHpsSelectedPFTausTrackPt1" ),
        MinE = cms.double( -1.0 ),
        triggerType = cms.int32( 84 ),
        MaxMass = cms.double( -1.0 )
    )

    process.HLTHPSDoublePFTauPt35Eta2p1Trk1 = cms.Sequence( process.hltHpsDoublePFTau35 + process.hltHpsPFTauTrackPt1Discriminator + process.hltHpsSelectedPFTausTrackPt1 + process.hltHpsDoublePFTau35TrackPt1 )


    process.hltHpsPFTauMediumAbsoluteChargedIsolationDiscriminator = cms.EDProducer( "PFRecoTauDiscriminationByIsolation",
        applyRhoCorrection = cms.bool( False ),
        PFTauProducer = cms.InputTag( "hltHpsPFTauProducer" ),
        storeRawOccupancy = cms.bool( False ),
        maximumSumPtCut = cms.double( 3.7 ),
        qualityCuts = cms.PSet(
          vertexTrackFiltering = cms.bool( False ),
          isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double( 0.2 ),
            minTrackPt = cms.double( 0.5 ),
            minGammaEt = cms.double( 0.5 ),
            minTrackHits = cms.uint32( 3 ),
            minTrackPixelHits = cms.uint32( 0 ),
            maxTrackChi2 = cms.double( 100.0 ),
            maxTransverseImpactParameter = cms.double( 0.1 ),
            useTracksInsteadOfPFHadrons = cms.bool( False )
          ),
          primaryVertexSrc = cms.InputTag( "hltPixelVertices" ),
          recoverLeadingTrk = cms.bool( False ),
          signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double( 0.2 ),
            minTrackPt = cms.double( 0.0 ),
            minGammaEt = cms.double( 0.5 ),
            minTrackHits = cms.uint32( 3 ),
            minTrackPixelHits = cms.uint32( 0 ),
            maxTrackChi2 = cms.double( 1000.0 ),
            maxTransverseImpactParameter = cms.double( 0.2 ),
            useTracksInsteadOfPFHadrons = cms.bool( False ),
            minNeutralHadronEt = cms.double( 1.0 )
          ),
          vxAssocQualityCuts = cms.PSet(
            minTrackPt = cms.double( 0.0 ),
            minGammaEt = cms.double( 0.5 ),
            minTrackHits = cms.uint32( 3 ),
            minTrackPixelHits = cms.uint32( 0 ),
            maxTrackChi2 = cms.double( 1000.0 ),
            maxTransverseImpactParameter = cms.double( 0.2 ),
            useTracksInsteadOfPFHadrons = cms.bool( False )
          ),
          pvFindingAlgo = cms.string( "closestInDeltaZ" )
        ),
        minTauPtForNoIso = cms.double( -99.0 ),
        vertexSrc = cms.InputTag( "NotUsed" ),
        applySumPtCut = cms.bool( True ),
        rhoConeSize = cms.double( 0.357 ),
        ApplyDiscriminationByTrackerIsolation = cms.bool( True ),
        storeRawPhotonSumPt_outsideSignalCone = cms.bool( False ),
        rhoUEOffsetCorrection = cms.double( 0.0 ),
        rhoProducer = cms.InputTag( "NotUsed" ),
        maxRelPhotonSumPt_outsideSignalCone = cms.double( 0.1 ),
        deltaBetaFactor = cms.string( "0.38" ),
        applyFootprintCorrection = cms.bool( False ),
        UseAllPFCandsForWeights = cms.bool( False ),
        relativeSumPtCut = cms.double( 0.03 ),
        Prediscriminants = cms.PSet(  BooleanOperator = cms.string( "and" ) ),
        applyOccupancyCut = cms.bool( False ),
        applyDeltaBetaCorrection = cms.bool( False ),
        WeightECALIsolation = cms.double( 0.33333 ),
        applyRelativeSumPtCut = cms.bool( False ),
        storeRawPUsumPt = cms.bool( False ),
        applyPhotonPtSumOutsideSignalConeCut = cms.bool( False ),
        maximumOccupancy = cms.uint32( 0 ),
        deltaBetaPUTrackPtCutOverride = cms.bool( True ),
        ApplyDiscriminationByWeightedECALIsolation = cms.bool( False ),
        maxAbsPhotonSumPt_outsideSignalCone = cms.double( 1.0E9 ),
        footprintCorrections = cms.VPSet(
          cms.PSet(  offset = cms.string( "0.0" ),
            selection = cms.string( "decayMode() = 0" )
          ),
          cms.PSet(  offset = cms.string( "0.0" ),
            selection = cms.string( "decayMode() = 1 || decayMode() = 2" )
          ),
          cms.PSet(  offset = cms.string( "2.7" ),
            selection = cms.string( "decayMode() = 5" )
          ),
          cms.PSet(  offset = cms.string( "0.0" ),
            selection = cms.string( "decayMode() = 6" )
          ),
          cms.PSet(  offset = cms.string( "max(2.0, 0.22*pt() - 2.0)" ),
            selection = cms.string( "decayMode() = 10" )
          )
        ),
        deltaBetaPUTrackPtCutOverride_val = cms.double( 0.5 ),
        ApplyDiscriminationByECALIsolation = cms.bool( False ),
        isoConeSizeForDeltaBeta = cms.double( 0.3 ),
        storeRawSumPt = cms.bool( False ),
        verbosity = cms.int32( 0 ),
        storeRawFootprintCorrection = cms.bool( False ),
        relativeSumPtOffset = cms.double( 0.0 ),
        customOuterCone = cms.double( -1.0 ),
        particleFlowSrc = cms.InputTag( "hltParticleFlow" )
    )


    process.hltHpsPFTauMediumRelativeChargedIsolationDiscriminator = cms.EDProducer( "PFRecoTauDiscriminationByIsolation",
        applyRhoCorrection = cms.bool( False ),
        PFTauProducer = cms.InputTag( "hltHpsPFTauProducer" ),
        storeRawOccupancy = cms.bool( False ),
        maximumSumPtCut = cms.double( 2.0 ),
        qualityCuts = cms.PSet(
          vertexTrackFiltering = cms.bool( False ),
          isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double( 0.2 ),
            minTrackPt = cms.double( 0.5 ),
            minGammaEt = cms.double( 0.5 ),
            minTrackHits = cms.uint32( 3 ),
            minTrackPixelHits = cms.uint32( 0 ),
            maxTrackChi2 = cms.double( 100.0 ),
            maxTransverseImpactParameter = cms.double( 0.1 ),
            useTracksInsteadOfPFHadrons = cms.bool( False )
          ),
          primaryVertexSrc = cms.InputTag( "hltPixelVertices" ),
          recoverLeadingTrk = cms.bool( False ),
          signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double( 0.2 ),
            minTrackPt = cms.double( 0.0 ),
            minGammaEt = cms.double( 0.5 ),
            minTrackHits = cms.uint32( 3 ),
            minTrackPixelHits = cms.uint32( 0 ),
            maxTrackChi2 = cms.double( 1000.0 ),
            maxTransverseImpactParameter = cms.double( 0.2 ),
            useTracksInsteadOfPFHadrons = cms.bool( False ),
            minNeutralHadronEt = cms.double( 1.0 )
          ),
          vxAssocQualityCuts = cms.PSet(
            minTrackPt = cms.double( 0.0 ),
            minGammaEt = cms.double( 0.5 ),
            minTrackHits = cms.uint32( 3 ),
            minTrackPixelHits = cms.uint32( 0 ),
            maxTrackChi2 = cms.double( 1000.0 ),
            maxTransverseImpactParameter = cms.double( 0.2 ),
            useTracksInsteadOfPFHadrons = cms.bool( False )
          ),
          pvFindingAlgo = cms.string( "closestInDeltaZ" )
        ),
        minTauPtForNoIso = cms.double( -99.0 ),
        vertexSrc = cms.InputTag( "NotUsed" ),
        applySumPtCut = cms.bool( False ),
        rhoConeSize = cms.double( 0.5 ),
        ApplyDiscriminationByTrackerIsolation = cms.bool( True ),
        storeRawPhotonSumPt_outsideSignalCone = cms.bool( False ),
        rhoUEOffsetCorrection = cms.double( 1.0 ),
        rhoProducer = cms.InputTag( "hltFixedGridRhoFastjetAll" ),
        maxRelPhotonSumPt_outsideSignalCone = cms.double( 0.1 ),
        deltaBetaFactor = cms.string( "0.38" ),
        applyFootprintCorrection = cms.bool( False ),
        UseAllPFCandsForWeights = cms.bool( False ),
        relativeSumPtCut = cms.double( 0.05 ),
        Prediscriminants = cms.PSet(  BooleanOperator = cms.string( "and" ) ),
        applyOccupancyCut = cms.bool( False ),
        applyDeltaBetaCorrection = cms.bool( False ),
        WeightECALIsolation = cms.double( 1.0 ),
        applyRelativeSumPtCut = cms.bool( True ),
        storeRawPUsumPt = cms.bool( False ),
        applyPhotonPtSumOutsideSignalConeCut = cms.bool( False ),
        maximumOccupancy = cms.uint32( 0 ),
        deltaBetaPUTrackPtCutOverride = cms.bool( True ),
        ApplyDiscriminationByWeightedECALIsolation = cms.bool( False ),
        maxAbsPhotonSumPt_outsideSignalCone = cms.double( 1.0E9 ),
        footprintCorrections = cms.VPSet(
          cms.PSet(  offset = cms.string( "0.0" ),
            selection = cms.string( "decayMode() = 0" )
          ),
          cms.PSet(  offset = cms.string( "0.0" ),
            selection = cms.string( "decayMode() = 1 || decayMode() = 2" )
          ),
          cms.PSet(  offset = cms.string( "2.7" ),
            selection = cms.string( "decayMode() = 5" )
          ),
          cms.PSet(  offset = cms.string( "0.0" ),
            selection = cms.string( "decayMode() = 6" )
          ),
          cms.PSet(  offset = cms.string( "max(2.0, 0.22*pt() - 2.0)" ),
            selection = cms.string( "decayMode() = 10" )
          )
        ),
        deltaBetaPUTrackPtCutOverride_val = cms.double( 0.5 ),
        ApplyDiscriminationByECALIsolation = cms.bool( False ),
        isoConeSizeForDeltaBeta = cms.double( 0.3 ),
        storeRawSumPt = cms.bool( False ),
        verbosity = cms.int32( 0 ),
        storeRawFootprintCorrection = cms.bool( False ),
        relativeSumPtOffset = cms.double( 60.0 ),
        customOuterCone = cms.double( -1.0 ),
        particleFlowSrc = cms.InputTag( "hltParticleFlow" )
    )
    process.hltHpsPFTauMediumAbsOrRelChargedIsolationDiscriminator = cms.EDProducer( "PFTauDiscriminatorLogicalAndProducer",
        Prediscriminants = cms.PSet(
          BooleanOperator = cms.string( "or" ),
          discr1 = cms.PSet(
            cut = cms.double( 0.5 ),
            Producer = cms.InputTag( "hltHpsPFTauMediumAbsoluteChargedIsolationDiscriminator" )
          ),
          discr2 = cms.PSet(
            cut = cms.double( 0.5 ),
            Producer = cms.InputTag( "hltHpsPFTauMediumRelativeChargedIsolationDiscriminator" )
          )
        ),
        FailValue = cms.double( 0.0 ),
        PassValue = cms.double( 1.0 ),
        PFTauProducer = cms.InputTag( "hltHpsPFTauProducer" )
    )

    process.HLTHPSMediumChargedIsoPFTauSequence = cms.Sequence( process.hltHpsPFTauMediumAbsoluteChargedIsolationDiscriminator + process.hltHpsPFTauMediumRelativeChargedIsolationDiscriminator + process.hltHpsPFTauMediumAbsOrRelChargedIsolationDiscriminator )
    process.hltHpsSelectedPFTausTrackPt1MediumChargedIsolation = cms.EDFilter( "PFTauSelector",
        discriminators = cms.VPSet(
          cms.PSet(  discriminator = cms.InputTag( "hltHpsPFTauTrackPt1Discriminator" ),
            selectionCut = cms.double( 0.5 )
          ),
          cms.PSet(  discriminator = cms.InputTag( "hltHpsPFTauMediumAbsOrRelChargedIsolationDiscriminator" ),
            selectionCut = cms.double( 0.5 )
          )
        ),
        discriminatorContainers = cms.VPSet(
        ),
        cut = cms.string( "pt > 0" ),
        src = cms.InputTag( "hltHpsPFTauProducer" )
    )
    process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolation = cms.EDFilter( "HLT1PFTau",
        saveTags = cms.bool( True ),
        MinPt = cms.double( 35.0 ),
        MinN = cms.int32( 2 ),
        MaxEta = cms.double( 2.1 ),
        MinEta = cms.double( -1.0 ),
        MinMass = cms.double( -1.0 ),
        inputTag = cms.InputTag( "hltHpsSelectedPFTausTrackPt1MediumChargedIsolation" ),
        MinE = cms.double( -1.0 ),
        triggerType = cms.int32( 84 ),
        MaxMass = cms.double( -1.0 )
    )

    process.hltHpsL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationMatch = cms.EDProducer( "L1THLTTauMatching",
        JetSrc = cms.InputTag( "hltHpsSelectedPFTausTrackPt1MediumChargedIsolation" ),
        EtMin = cms.double( 0.0 ),
        L1TauTrigger = cms.InputTag( "hltL1sDoubleTauBigOR" )
    )
    process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationL1HLTMatched = cms.EDFilter( "HLT1PFTau",
        saveTags = cms.bool( True ),
        MinPt = cms.double( 35.0 ),
        MinN = cms.int32( 2 ),
        MaxEta = cms.double( 2.1 ),
        MinEta = cms.double( -1.0 ),
        MinMass = cms.double( -1.0 ),
        inputTag = cms.InputTag( "hltHpsL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationMatch" ),
        MinE = cms.double( -1.0 ),
        triggerType = cms.int32( 84 ),
        MaxMass = cms.double( -1.0 )
    )
    process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02 = cms.EDFilter( "HLTPFTauPairDzMatchFilter",
        saveTags = cms.bool( True ),
        TriggerType = cms.int32( 84 ),
        JetSrc = cms.InputTag( "hltHpsL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationMatch" ),
        JetMinPt = cms.double( 35.0 ),
        JetMaxDZ = cms.double( 0.2 ),
        JetMinDR = cms.double( 0.5 ),
        JetMaxEta = cms.double( 2.1 )
    )

    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4 = cms.Path(
    process.HLTBeginSequence +
    process.hltL1sDoubleTauBigOR +
    process.hltPreDoubleMediumChargedIsoPFTauHPS35Trk1eta2p1Reg +
    process.L2Sequence +
    process.HLTGlobalPFTauHPSSequence +
    process.HLTHPSDoublePFTauPt35Eta2p1Trk1 +
    process.HLTHPSMediumChargedIsoPFTauSequence +
    process.hltHpsSelectedPFTausTrackPt1MediumChargedIsolation +
    process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolation +
    process.hltHpsL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationMatch +
    process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationL1HLTMatched +
    process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02+
    process.HLTEndSequence,
    process.HLTDoLocalPixelTask,
    process.HLTRecoPixelTracksTask,
    process.HLTRecopixelvertexingTask
    )



    process.schedule = cms.Schedule(*[ process.HLTriggerFirstPath, process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4, process.HLTriggerFinalPath, process.endjob_step ], tasks=[process.patAlgosToolsTask])

    #process.options.wantSummary = cms.untracked.bool(True)
    process.options.wantSummary = cms.untracked.bool(False)
    return process
