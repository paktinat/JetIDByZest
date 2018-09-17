import FWCore.ParameterSet.Config as cms
#N="119"
#N=$1
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
#import FWCore.PythonUtilities.LumiList as LumiList
#myLumis = LumiList.LumiList(filename='Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt').getCMSSWString().split(',')
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
#process.source.lumisToProcess.extend(myLumis)


process.countPatElectrons.maxNumber = cms.uint32(0)
#comment for DY, uncomment for QCD
#process.countPatMuons.maxNumber = cms.uint32(0)
#process.countPatJets.minNumber = cms.uint32(2)
#comment for QCD, uncomment for DY
process.countPatMuons.minNumber = cms.uint32(2)
process.countPatJets.minNumber = cms.uint32(1)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#input file
import FWCore.Utilities.FileUtils as FileUtils
#files2011data = FileUtils.loadListFromFile ('FilesName/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file_index/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file_index_00'+N)
files2011data = FileUtils.loadListFromFile ('FilesName/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file_index/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file_index_01'+N)#It goes up to 224


readFiles = cms.untracked.vstring( *files2011data )
process.source.fileNames = readFiles


process.demo = cms.EDAnalyzer('MyTreeCounter',
#                              histoname = cms.string("file:///cmsdata2/paktinat/JetIDByZest/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file/HistoTrksMuVertex0/HistoTrksMuVertex_00"+N+".root"),
                              histoname = cms.string("file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_QCD_Pt_120to170_TuneZ2_7TeV_pythia6_AODSIM_file/HistoTrksMuVertex0/HistoTrksMuVertex_01"+N+".root"),
                              readFromEOS = cms.bool(True)
)

 
process.nEventsTotal = cms.EDProducer("EventCountProducer")
process.nEventsFiltered = cms.EDProducer("EventCountProducer")

process.p = cms.Path(process.nEventsTotal+process.nEventsFiltered+process.demo )
