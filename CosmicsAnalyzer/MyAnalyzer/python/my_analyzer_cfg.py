# cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile
#     https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookWriteFrameworkModule
# 
# Import CMS python class definitions such as Process, Source, and EDProducer
import FWCore.ParameterSet.Config as cms

# Set up a process, named MYANALYZER in this case
process = cms.Process("MYANALYZER")

process.load("FWCore.MessageService.MessageLogger_cfi")

# -1 means all events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Configure the object that reads the input file
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/user/d/dazhang/private/cms/cosmics/data/2ad63d9f-234b-4c68-8c18-49d485d42bc7.root'
                )
                            )

# configure the analysis module named "MyAnalyzer"
process.muonPhiAnalyzer = cms.EDAnalyzer("MyAnalyzer",
         muonCollection = cms.InputTag("muons")
                                        )

# configure the output file for the analysis results
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
        "/afs/cern.ch/user/d/dazhang/private/cms/CMSSW_13_0_10/src/CosmicsAnalyzer/MyAnalyzer/test/muon_phi_ntuple.root"
        )
                                    )

# Configure a path to run the analyzer module
process.p = cms.Path(process.muonPhiAnalyzer)