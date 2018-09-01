import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    limit = cms.untracked.int32(-1)
#)
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

## ------------------------------------------------------
#  NOTE: you can use a bunch of core tools of PAT to
#  taylor your PAT configuration; for a few examples
#  uncomment the lines below
## ------------------------------------------------------
from PhysicsTools.PatAlgos.tools.coreTools import *
#removeMCMatching(process, ['All'])
#To be commented for MC

## remove certain objects from the default sequence
#removeAllPATObjectsBut(process, ['Muons','Electrons'])
#removeAllPATObjectsBut(process, ['Jets'])

# make sure to keep the created objects
#process.out.outputCommands += ['drop *']
#process.out.outputCommands += ['keep *_offlinePrimaryVertices_*_*']
process.out.outputCommands += ['keep *_pat*_*_*',]
process.out.outputCommands += ['keep *_selectedPatJets_*_*',]
process.out.outputCommands += ['keep *_genParticles_*_*',]
#process.out.outputCommands += ['keep *',]


from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process, 
                    cms.InputTag('ak5PFJets'),   
                    jetCorrLabel     = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute']),  
                    jetIdLabel   = "ak5",
                    btagdiscriminators=['jetBProbabilityBJetTags', 'jetProbabilityBJetTags', 'trackCountingHighPurBJetTags', 'trackCountingHighEffBJetTags', 'negativeOnlyJetBProbabilityJetTags', 'negativeOnlyJetProbabilityJetTags', 'negativeTrackCountingHighEffJetTags', 'negativeTrackCountingHighPurJetTags', 'positiveOnlyJetBProbabilityJetTags', 'positiveOnlyJetProbabilityJetTags', 'simpleSecondaryVertexHighEffBJetTags', 'simpleSecondaryVertexHighPurBJetTags', 'combinedSecondaryVertexBJetTags', 'combinedSecondaryVertexPositiveBJetTags', 'combinedSecondaryVertexMVABJetTags'] 
                    ) 

## ------------------------------------------------------
#  In addition you usually want to change the following
#  parameters:
## ------------------------------------------------------
#
#globaltag
#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')
#process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'
process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1.db')
process.GlobalTag.globaltag = 'START53_LV6A1::All'

#luminosity Only for data
#import FWCore.ParameterSet.Config as cms
#import FWCore.PythonUtilities.LumiList as LumiList
#myLumis = LumiList.LumiList(filename='Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt').getCMSSWString().split(',')
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
#process.source.lumisToProcess.extend(myLumis)

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        #        'file:myfile.root'
        #'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/ElectronHad/AOD/12Oct2013-v1/20001/001F9231-F141-E311-8F76-003048F00942.root'
        #'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/MuEG/AOD/12Oct2013-v1/20000/FAFFBF92-173E-E311-9B05-003048F179BA.root'
#        'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/BTag/AOD/12Oct2013-v1/20000/F8A6DAE8-013E-E311-87AC-002354EF3BDE.root'
#        'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/70000/FC6FD4D7-F0EC-E411-B4CA-0025905A6060.root'    
        'file:MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6-v1_00000_FE4D03BD-7DBD-E311-B0D1-0025905A6080.root'
        )
                            )

process.patElectrons.useParticleFlow = cms.bool(True)


process.selectedPatJets.cut = cms.string('pt > 20. && abs(eta) < 2.5 && numberOfDaughters > 1 && neutralHadronEnergyFraction < 0.99 && neutralEmEnergyFraction < 0.99 && (chargedEmEnergyFraction < 0.99 || abs(eta) > 2.4) && (chargedHadronEnergyFraction > 0. || abs(eta) >= 2.4) && (chargedMultiplicity > 0 || abs(eta) >= 2.4)')

process.selectedPatElectrons.cut = cms.string('pt > 10 && abs(eta) < 2.5 && (trackIso + caloIso)/pt < 0.2 && electronID("eidRobustTight") > 0.')

process.selectedPatMuons.cut = cms.string('isPFMuon && (isGlobalMuon || isTrackerMuon) && pt > 10 && abs(eta) < 2.5 && (trackIso + caloIso)/pt < 0.2')

process.countPatElectrons.maxNumber = cms.uint32(0)
#process.countPatMuons.maxNumber = cms.uint32(0)
#process.countPatJets.minNumber = cms.uint32(3)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000) )


#input file
#import FWCore.Utilities.FileUtils as FileUtils
#files2011data = FileUtils.loadListFromFile ('CMS_Run2011A_BTag_AOD_12Oct2013-v1_20000_file_index/CMS_Run2011A_BTag_AOD_12Oct2013-v1_20000_file_index_0019') 
#files2011data = FileUtils.loadListFromFile ('FilesName/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file_index/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file_index_0028')
#readFiles = cms.untracked.vstring( *files2011data )
#process.source.fileNames = readFiles

#output file
#process.out.fileName = 'file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file/DYJets_M_50_0028_2XXX.root' ##  (e.g. 'myTuple.root')
process.out.fileName = 'file://Jets_PAT_data_500files.root' ##  (e.g. 'myTuple.root')

process.demo = cms.EDAnalyzer('MyTreeProducer',
#                              histoname = cms.string("file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file/Histo_DYJets_0028_2XXXX.root")
                              histoname = cms.string("file://Histo_DY_0001.root"),
                              readFromEOS = cms.bool(False)
)


#from CommonTools.UtilAlgos.EventCountProducer import *
process.nEventsTotal = cms.EDProducer("EventCountProducer")
process.nEventsFiltered = cms.EDProducer("EventCountProducer")

process.p = cms.Path(process.nEventsTotal+process.patDefaultSequence+process.nEventsFiltered+process.demo )
#process.p = cms.Path(process.patDefaultSequence+process.demo)
