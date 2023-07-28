# cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile
#     https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookWriteFrameworkModule
# 
# Import CMS python class definitions such as Process, Source, and EDProducer
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# cf: https://twiki.cern.ch/twiki/bin/view/CMS/CommandLineOptions
#     FWCore/ParameterSet/python/VarParsing.py
options = VarParsing('analysis')
options.register('GTAG', '106X_upgrade2018_realistic_v11BasedCandidateTmp_2022_08_09_01_32_34',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Global Tag"
)
options.parseArguments()

# Set up a process, named MYANALYZER in this case
process = cms.Process("MYANALYZER")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.Services_cff')


# -1 means all events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Configure the object that reads the input file
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/user/d/dazhang/private/cms/cosmics/data/2ad63d9f-234b-4c68-8c18-49d485d42bc7.root'
                )
                            )

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GTAG, '')

# configure the analysis module named "MyAnalyzer"
process.muonAnalyzer = cms.EDAnalyzer("MyAnalyzer",
         muonCollection = cms.InputTag("muons"),
         cosmicMuonCollection = cms.InputTag("cosmicMuons")
                                        )

process.muon1LegAnalyzer = cms.EDAnalyzer("MyAnalyzer",
         muonCollection = cms.InputTag("muons1Leg"),
         cosmicMuonCollection = cms.InputTag("cosmicMuons1Leg")
                                        )

process.muonEndCapsOnlyAnalyzer = cms.EDAnalyzer("MyAnalyzer",
         muonCollection = cms.InputTag("muonsBeamHaloEndCapsOnly"),
         cosmicMuonCollection = cms.InputTag("cosmicMuonsEndCapsOnly")
                                        )

process.muonNoRPCAnalyzer = cms.EDAnalyzer("MyAnalyzer",
         muonCollection = cms.InputTag("muonsNoRPC"),
         cosmicMuonCollection = cms.InputTag("cosmicMuonsNoRPC")
                                        )

process.muonWitht0CorrectionAnalyzer = cms.EDAnalyzer("MyAnalyzer",
         muonCollection = cms.InputTag("muonsWitht0Correction"),
         cosmicMuonCollection = cms.InputTag("cosmicMuonsWitht0Correction")
                                        )

process.splitmuonAnalyzer = cms.EDAnalyzer("MyAnalyzer",
         muonCollection = cms.InputTag("splitMuons"),
         cosmicMuonCollection = cms.InputTag("cosmicMuons")
                                        )

process.lhcSTAmuonAnalyzer = cms.EDAnalyzer("MyAnalyzer",
         muonCollection = cms.InputTag("lhcSTAMuons"),
         cosmicMuonCollection = cms.InputTag("cosmicMuons")
                                        )


# configure the output file for the analysis results
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
        "/afs/cern.ch/user/d/dazhang/private/cms/CMSSW_13_0_10/src/CosmicsAnalyzer/MyAnalyzer/test/muon_ntuple.root"
        )
                                    )

# Configure a path to run the analyzer module
process.p1 = cms.Path(process.muonAnalyzer)
# process.p2 = cms.Path(process.muon1LegAnalyzer)
# process.p3 = cms.Path(process.muonEndCapsOnlyAnalyzer)
# process.p4 = cms.Path(process.muonNoRPCAnalyzer)
# process.p5 = cms.Path(process.muonWitht0CorrectionAnalyzer)
# process.p6 = cms.Path(process.splitmuonAnalyzer)
# process.p7 = cms.Path(process.lhcSTAmuonAnalyzer)