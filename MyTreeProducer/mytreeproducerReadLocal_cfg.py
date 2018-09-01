import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

#removeMCMatching(process, ['All'])
#To be commented for MC

## ------------------------------------------------------
#  In addition you usually want to change the following
#  parameters:
## ------------------------------------------------------
#
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#To be changed for MC
#globaltag
process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')
process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'
#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1.db')
#process.GlobalTag.globaltag = 'START53_LV6A1::All'

#luminosity Only for data not MC
#import FWCore.ParameterSet.Config as cms
#import FWCore.PythonUtilities.LumiList as LumiList
#myLumis = LumiList.LumiList(filename='Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt').getCMSSWString().split(',')
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
#process.source.lumisToProcess.extend(myLumis)

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
#        'file:Jets_PAT_data_500files2_PFCandSelectedTracks.root'
#        'file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_QCD_Pt_120to170_TuneZ2_7TeV_pythia6_AODSIM_file/QCD_Pt_120to170_0015.root'
        #'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/ElectronHad/AOD/12Oct2013-v1/20001/001F9231-F141-E311-8F76-003048F00942.root'
        #'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/MuEG/AOD/12Oct2013-v1/20000/FAFFBF92-173E-E311-9B05-003048F179BA.root'
#        'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/BTag/AOD/12Oct2013-v1/20000/F8A6DAE8-013E-E311-87AC-002354EF3BDE.root'
        #'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/70000/FC6FD4D7-F0EC-E411-B4CA-0025905A6060.root'    
#        'file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file/DYJets_M_50_0220.root',
#        'file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file/DYJets_M_50_0221.root',
#        'file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file/DYJets_M_50_0222.root',
#        'file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file/DYJets_M_50_0223.root',
#        'file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file/DYJets_M_50_0224.root'
#        'file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file/DYJets_M_50_0225.root',
#        'file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file/DYJets_M_50_0226.root',
#        'file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file/DYJets_M_50_0227.root',
#        'file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file/DYJets_M_50_0228.root',
#        'file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file/DYJets_M_50_0229.root'
        'file:///cmsdata2/paktinat/JetIDByZest/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file/DoubleMu_Zpeak_0000.root',
        'file:///cmsdata2/paktinat/JetIDByZest/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file/DoubleMu_Zpeak_0001.root',
        'file:///cmsdata2/paktinat/JetIDByZest/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file/DoubleMu_Zpeak_0002.root',
        'file:///cmsdata2/paktinat/JetIDByZest/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file/DoubleMu_Zpeak_0003.root',
        'file:///cmsdata2/paktinat/JetIDByZest/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file/DoubleMu_Zpeak_0004.root',
        'file:///cmsdata2/paktinat/JetIDByZest/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file/DoubleMu_Zpeak_0005.root',
        'file:///cmsdata2/paktinat/JetIDByZest/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file/DoubleMu_Zpeak_0006.root',
        'file:///cmsdata2/paktinat/JetIDByZest/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file/DoubleMu_Zpeak_0007.root',
        'file:///cmsdata2/paktinat/JetIDByZest/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file/DoubleMu_Zpeak_0008.root',
        'file:///cmsdata2/paktinat/JetIDByZest/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file/DoubleMu_Zpeak_0009.root'


        )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#input file
#import FWCore.Utilities.FileUtils as FileUtils
#files2011data = FileUtils.loadListFromFile ('CMS_Run2011A_BTag_AOD_12Oct2013-v1_20000_file_index/CMS_Run2011A_BTag_AOD_12Oct2013-v1_20000_file_index_0019') 
#files2011data = FileUtils.loadListFromFile ('CMS_MonteCarlo2011_Summer11LegDR_QCD_Pt_120to170_TuneZ2_7TeV_pythia6_AODSIM_file_index/CMS_MonteCarlo2011_Summer11LegDR_QCD_Pt_120to170_TuneZ2_7TeV_pythia6_AODSIM_file_index_0015')
#files2011data = FileUtils.loadListFromFile ('FilesName/CMS_MonteCarlo2011_Summer11LegDR_QCD_Pt_120to170_TuneZ2_7TeV_pythia6_AODSIM_cmsdata2_file_index.txt')
#readFiles = cms.untracked.vstring( *files2011data )
#process.source.fileNames = readFiles
#process.source = cms.Source('PoolSource', fileNames = readFiles)


#output file
#process.out.fileName = 'file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_QCD_Pt_120to170_TuneZ2_7TeV_pythia6_AODSIM_file/QCD_Pt_120to170_0015.root' ##  (e.g. 'myTuple.root')
#process.outputModules.fileName = 'file://Jets_PAT_data_500files.root' ##  (e.g. 'myTuple.root')

process.demo = cms.EDAnalyzer('MyTreeProducer',
#                              histoname = cms.string("file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_QCD_Pt_120to170_TuneZ2_7TeV_pythia6_AODSIM_file/Histo_QCD_Pt_120to170_0015.root")
#                              histoname = cms.string("file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_QCD_Pt_120to170_TuneZ2_7TeV_pythia6_AODSIM_file/Histo0.root"),
#                              histoname = cms.string("file:Histo22_PFCandSelectedTracks.root"),

#                              histoname = cms.string("file:///cmsdata2/paktinat/JetIDByZest/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_AODSIM_file/HistoTrks/Histo_DYJets_M_50_022.root"),
                              histoname = cms.string("file:///cmsdata2/paktinat/JetIDByZest/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file/HistoTrksMuVertex/HistoDoubleMu_Zpeak_000.root"),
                              readFromEOS = cms.bool(False)
)

#from CommonTools.UtilAlgos.EventCountProducer import *
#process.nEventsTotal = cms.EDProducer("EventCountProducer")
#process.nEventsFiltered = cms.EDProducer("EventCountProducer")

#process.p = cms.Path(process.nEventsTotal+process.patDefaultSequence+process.nEventsFiltered+process.demo )
#process.p = cms.Path(process.patDefaultSequence+process.demo)
process.p = cms.Path(process.demo)
