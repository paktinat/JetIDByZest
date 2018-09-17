// -*- C++ -*-
//
// Package:    MyTreeProducer
// Class:      MyTreeCounter
// 
/**\class MyTreeCounter MyTreeCounter.cc JetIDByZest/MyTreeCounter/src/MyTreeCounter.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Mon Jul 31 13:55:17 IRDT 2017
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Jet from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetAnalysis
//#include "DataFormats/JetReco/interface/PFJet.h"
//#include "DataFormats/JetReco/interface/PFJetCollection.h"

//#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
//Jet
//Histo
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include <DataFormats/Common/interface/MergeableCounter.h>
#include "DataFormats/TrackReco/interface/Track.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include <Math/VectorUtil.h>
//


//
// class declaration
//

class MyTreeCounter : public edm::EDAnalyzer {
public:
  explicit MyTreeCounter(const edm::ParameterSet&);
  ~MyTreeCounter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  std::string myHistoName;
  bool isReadFromEOS;
  TFile* myFile;
  TTree*  mytree;
  TH1F* myHisto;  
  //  TH1F* myHisto2;  
  Int_t nJets;
  Int_t nTrks;
  Int_t nMus;
  Int_t nGenPars;
  Int_t nGenTrks;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MyTreeCounter::MyTreeCounter(const edm::ParameterSet& iConfig)

{ //now do what ever initialization is needed
  // myFile = new TFile("Histo.root","RECREATE");
  myHistoName = iConfig.getParameter<std::string>("histoname");
  isReadFromEOS = iConfig.getParameter<bool>("readFromEOS");

  myFile = new TFile(myHistoName.c_str(),"RECREATE");
  myHisto = new TH1F("myHisto","myHisto", 2, 0, 2);

};

MyTreeCounter::~MyTreeCounter()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  myFile->cd();
  myHisto->Write();

  myFile->Write();
  myFile->Close();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyTreeCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

}

// ------------ method called once each job just before starting event loop  ------------
void 
MyTreeCounter::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyTreeCounter::endJob() 
{

}

// ------------ method called when starting to processes a run  ------------
void 
MyTreeCounter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MyTreeCounter::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MyTreeCounter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MyTreeCounter::endLuminosityBlock(const edm::LuminosityBlock & lumi, const edm::EventSetup & iSetup)
{
  int nEventsTotalCntr = 1;
  int nEventsFilteredCntr = 1;
  if(isReadFromEOS){
    // Total number of events is the sum of the events in each of these luminosity blocks
    edm::Handle<edm::MergeableCounter> nEventsTotalCounter;
    lumi.getByLabel("nEventsTotal", nEventsTotalCounter);
    nEventsTotalCntr = nEventsTotalCounter->value;
    //std::cout<<" nEventsTotal "<<nEventsTotalCntr<<std::endl;

    edm::Handle<edm::MergeableCounter> nEventsFilteredCounter;
    lumi.getByLabel("nEventsFiltered", nEventsFilteredCounter);
    nEventsFilteredCntr = nEventsFilteredCounter->value;
    //std::cout<<" nEventsFiltered "<<nEventsFilteredCntr<<std::endl;
  }

  myHisto->Fill(0.5, nEventsTotalCntr);
  myHisto->Fill(1.5, nEventsFilteredCntr);

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyTreeCounter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyTreeCounter);
