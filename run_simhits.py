import FWCore.ParameterSet.Config as cms

process = cms.Process('TRACKANA')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hit.root')
)

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '90X_upgrade2017_realistic_v20', '')

# Input source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/w/wangx/workSpace/public/tracking2015/testSample/PbPbstep2_DIGI2017_9.root')
    #fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/w/wangx/workSpace/public/tracking2015/testSample/step3_8.root')
)



process.demo = cms.EDAnalyzer('simHitsAnalyzer',
	allowDifferentSimHitProcesses = cms.bool(True),
        tpHitsSrc = cms.InputTag('mix', 'MergedTrackTruth'),
        hitHighTofSrc = cms.InputTag('g4SimHits', 'TrackerHitsPixelBarrelHighTof'),
        hitLowTofSrc = cms.InputTag('g4SimHits', 'TrackerHitsPixelBarrelLowTof'),
        hitEndHighTofSrc = cms.InputTag('g4SimHits', 'TrackerHitsPixelEndcapHighTof'),
        hitEndLowTofSrc = cms.InputTag('g4SimHits',  'TrackerHitsPixelEndcapLowTof'),
        hitTECHighTofSrc = cms.InputTag('g4SimHits', 'TrackerHitsTECHighTof'),
        hitTECLowTofSrc = cms.InputTag('g4SimHits',  'TrackerHitsTECLowTof'),
        hitTIBHighTofSrc = cms.InputTag('g4SimHits', 'TrackerHitsTIBHighTof'),
        hitTIBLowTofSrc = cms.InputTag('g4SimHits',  'TrackerHitsTIBLowTof'),
        hitTIDHighTofSrc = cms.InputTag('g4SimHits', 'TrackerHitsTIDHighTof'),
        hitTIDLowTofSrc = cms.InputTag('g4SimHits',  'TrackerHitsTIDLowTof'),
        hitTOBHighTofSrc = cms.InputTag('g4SimHits', 'TrackerHitsTOBHighTof'),
        hitTOBLowTofSrc = cms.InputTag('g4SimHits',  'TrackerHitsTOBLowTof'),
        simTrackSrc = cms.InputTag('g4SimHits', ''),
)

process.p = cms.Path(
		      process.demo
)
