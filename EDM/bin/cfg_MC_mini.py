import FWCore.ParameterSet.Config as cms

process = cms.Process("bphAnalysis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

process.MessageLogger.categories.extend(["GetManyWithoutRegistration","GetByLabelWithoutRegistration"])
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/RECOSIMstep_94X_mc2017_realistic_v10-v1/70000/FE222926-00EC-E711-A352-782BCB20FB9D.root')
                            #fileNames = cms.untracked.vstring('/store/mc/RunIIAutumn18MiniAOD/BsToJpsiPhi_JpsiPhiFilterDG0_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/90000/BF3580C1-D5E4-7645-AEAC-F6A60D3152A0.root',)
                            #fileNames = cms.untracked.vstring('/store/mc/RunIIFall17MiniAODv2/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/80000/F69E6795-7CCA-E811-B4C8-EC0D9A0B3320.root',)
                            fileNames = cms.untracked.vstring('/store/mc/RunIIFall17MiniAODv2/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/80000/F6035168-38CA-E811-A5EE-A0369F3102B6.root',
                                                              '/store/mc/RunIIFall17MiniAODv2/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/80000/DC177185-BEB7-E811-9616-0025905A60DE.root')
    #fileNames = cms.untracked.vstring('/store/mc/RunIIFall17MiniAODv2/BsToJpsiPhi_BMuonFilter_DGamma0_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/1400F57A-1C43-E811-A00A-7CD30ACE2445.root',)
)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

### vtxTagInfo

process.load('RecoBTag/SoftLepton/softLepton_cff')

#process.load('RecoBTag/SoftLepton/softPFMuonTagInfos_cfi')
#process.load('RecoBTag/SoftLepton/softPFElectronTagInfos_cfi')

process.load('RecoBTag/SecondaryVertex/pfInclusiveSecondaryVertexFinderTagInfos_cfi')
process.load('RecoBTag/ImpactParameter/pfImpactParameterTagInfos_cfi')

#process.load('RecoBTag/SecondaryVertex/secondaryVertexTagInfos_cfi')

process.softPFMuonsTagInfos.primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices")
process.softPFMuonsTagInfos.jets = cms.InputTag("slimmedJets")
process.softPFMuonsTagInfos.muons = cms.InputTag("slimmedMuons")

process.softPFElectronsTagInfos.primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices")
process.softPFElectronsTagInfos.jets = cms.InputTag("slimmedJets")
process.softPFElectronsTagInfos.electrons = cms.InputTag("slimmedElectrons")

process.pfImpactParameterTagInfos.primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pfImpactParameterTagInfos.jets = cms.InputTag("slimmedJets")
process.pfImpactParameterTagInfos.candidates = cms.InputTag("packedPFCandidates")

process.tagInfoProd = cms.Sequence(
    process.softPFMuonsTagInfos
    + process.softPFElectronsTagInfos
    + process.pfImpactParameterTagInfos
    * process.pfSecondaryVertexTagInfos
    )

### vtxTagInfo end


process.pdAnalyzer = cms.EDAnalyzer('PDNtuplizer',

    ## optional
    eventList = cms.string('evtlist'),
    verbose = cms.untracked.string('f'),
    evStamp = cms.untracked.string('f'),

    ## mandatory
    ## ntuple file name: empty string to drop ntuple filling
    ntuName = cms.untracked.string('ntu.root'),
    ## histogram file name
    histName = cms.untracked.string('his.root'),

    labelTrigResults  = cms.string('TriggerResults::HLT'), 
    labelTrigEvent    = cms.string(''),
    labelTrigObjects  = cms.string('slimmedPatTrigger'),
    labelFilterLabels = cms.string('slimmedPatTrigger:filterLabels:PAT'),
    labelBeamSpot     = cms.string('offlineBeamSpot'),
    labelMets         = cms.string('slimmedMETs'),
    labelMuons        = cms.string('slimmedMuons'),
    labelElectrons    = cms.string('slimmedElectrons'),
    labelConversions  = cms.string('reducedEgamma:reducedConversions'),
#    labelTaus         = cms.string('slimmedTaus'),
    labelTaus         = cms.string(''),

    labelJets         = cms.string('slimmedJets'),
    labelPCCandidates = cms.string('packedPFCandidates::PAT'),

    labelGeneralTracks = cms.string(''),
    labelPVertices    = cms.string('offlineSlimmedPrimaryVertices'),
    labelSVertices    = cms.string('pfSecondaryVertexTagInfos'),
    labelSVTagInfo    = cms.string(''),
    labelPUInfo       = cms.string(''),
    labelGen          = cms.string('prunedGenParticles'),
    labelGPJ          = cms.string('slimmedGenJets'),

    labelCSV          = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
    labelTCHE         = cms.string('trackCountingHighEffBJetTags'),
    labelTags         = cms.vstring('pfDeepCSVJetTags:probudsg',
                                    'pfDeepCSVJetTags:probc',
                                    'pfDeepCSVJetTags:probcc',
                                    'pfDeepCSVJetTags:probb',
                                    'pfDeepCSVJetTags:probbb'),

    vs = cms.VPSet(
      cms.PSet( type  = cms.string('svtK0short'),
                label = cms.string('slimmedKshortVertices::PAT') ),
      cms.PSet( type  = cms.string('svtLambda0'),
                label = cms.string('slimmedLambdaVertices::PAT') ),
    ),
    vertReco = cms.vstring('svtBuJPsiK','svtBdJPsiKx','svtBsJPsiPhi'),

    acceptNewTrigPaths = cms.string('f'),
    write_hltlist = cms.string('f'),

    selectAssociatedPF = cms.string('f'),
    selectAssociatedTk = cms.string('f'),
    recoverMuonTracks = cms.string('t'),
    writeAllPrimaryVertices = cms.string('t'),

    jetPtMin  = cms.double(  5.0 ),
    jetEtaMax = cms.double(  2.5 ),
    trkDzMax  = cms.double(  0.8 ),
    trkPtMin  = cms.double(  0.5 ),
    trkEtaMax = cms.double(  3.0 ),
    dRmatchHLT = cms.double( 0.5 ),
    dPmatchHLT = cms.double( 0.5 ),


    savedTriggerPaths = cms.vstring(
        '*'
    ),

    ## trigger objects to save on ntuple:
    savedTriggerObjects = cms.vstring(
        'hltJet',
        'hltMuon',
        'hltElectron',
        'hltTrack'
    ),

    ## jet user info to save on ntuple:
    savedJetInfo = cms.vstring(
#       'puBeta*'
    ),

    ## Electron user info to save on ntuple:
    savedEleInfo = cms.vstring(
       '*'
    )

)
    
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False, #as energy corrections are not yet availible for 2018
                       era='2018-Prompt')      

#process.evNumFilter = cms.EDFilter('EvNumFilter',
#    eventList = cms.string('evList')
#)

#############################################################
#### PATH definition
#############################################################
# Let it run
process.p = cms.Path(
#    process.evNumFilter *
    process.egammaPostRecoSeq *
    process.tagInfoProd *
    process.pdAnalyzer
    )

