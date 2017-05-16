import FWCore.ParameterSet.Config as cms

process = cms.Process('TRACKANA')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
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


# Input source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/w/wangx/workSpace/public/tracking2015/testSample/PbPbstep2_DIGI2017_9.root')
    #fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/w/wangx/workSpace/public/tracking2015/testSample/step3_8.root')
)



process.demo = cms.EDAnalyzer('simHitsAnalyzer',
        tpHitsSrc = cms.InputTag('mix', 'MergedTrackTruth'),
        hitSrc = cms.InputTag('g4SimHits', 'TrackerHitsPixelBarrelHighTof')
)

process.p = cms.Path(
		      process.demo
)
