// -*- C++ -*-
//
// Package:    MyTreeProducer
// Class:      MyTreeProducer
// 
/**\class MyTreeProducer MyTreeProducer.cc JetIDByZest/MyTreeProducer/src/MyTreeProducer.cc

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
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include <Math/VectorUtil.h>
//


//
// class declaration
//

class MyTreeProducer : public edm::EDAnalyzer {
public:
  explicit MyTreeProducer(const edm::ParameterSet&);
  ~MyTreeProducer();

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
  TFile* myFile;
  TTree*  mytree;
  TH1F* myHisto;  
  Int_t nJets;
  Int_t nGenPars;
  static const Int_t maxObjects = 30;
  // ----------member data ---------------------------
  //Kinematic
  Float_t genParPt[maxObjects];
  Float_t genParEta[maxObjects];
  Float_t genParPhi[maxObjects];
  Float_t genParMass[maxObjects];
  //ID
  Int_t genParID[maxObjects];
  Int_t genParMotherID[maxObjects];
  Float_t genParMotherPz[maxObjects];


  //Kinematic
  Float_t jetPt[maxObjects];
  Float_t jetEta[maxObjects];
  Float_t jetPhi[maxObjects];
  Float_t jetMass[maxObjects];
  //BTag Info
  Float_t jetBProbability[maxObjects];
  Float_t jetProbability[maxObjects];
  Float_t trackCountingHighPur[maxObjects];
  Float_t trackCountingHighEff[maxObjects];
  Float_t negativeOnlyJetBProbability[maxObjects];
  Float_t negativeOnlyJetProbability[maxObjects];
  Float_t negativeTrackCountingHighEff[maxObjects];
  Float_t negativeTrackCountingHighPur[maxObjects];
  Float_t positiveOnlyJetBProbability[maxObjects];
  Float_t positiveOnlyJetProbability[maxObjects];
  Float_t simpleSecondaryVertexHighEff[maxObjects];
  Float_t simpleSecondaryVertexHighPur[maxObjects];
  Float_t combinedSecondaryVertex[maxObjects];
  Float_t combinedSecondaryVertexPositive[maxObjects];
  Float_t combinedSecondaryVertexMVA[maxObjects];
  //Zest
  Float_t jetZest[maxObjects];
  Float_t jetPPerp[maxObjects];
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
MyTreeProducer::MyTreeProducer(const edm::ParameterSet& iConfig)

{ //now do what ever initialization is needed
  // myFile = new TFile("Histo.root","RECREATE");
  myHistoName = iConfig.getParameter<std::string>("histoname");
  myFile = new TFile(myHistoName.c_str(),"RECREATE");
  myHisto = new TH1F("myHisto","myHisto", 2, 0, 2);
  mytree  = new TTree("JetInfo", "Jet Information");
  mytree->Branch("nJets",&nJets,"nJets/I");
  //Kinematic
  mytree->Branch("jetPt",&jetPt,"jetPt[nJets]/F");
  mytree->Branch("jetEta",&jetEta,"jetEta[nJets]/F");
  mytree->Branch("jetPhi",&jetPhi,"jetPhi[nJets]/F");
  mytree->Branch("jetMass",&jetMass,"jetMass[nJets]/F");
  //BTag Info
  mytree->Branch("jetBProbability", &jetBProbability,"jetBProbability[nJets]/F");
  mytree->Branch("jetProbability", &jetProbability,"jetProbability[nJets]/F");
  mytree->Branch("trackCountingHighPur", &trackCountingHighPur,"trackCountingHighPur[nJets]/F");
  mytree->Branch("trackCountingHighEff", &trackCountingHighEff,"trackCountingHighEff[nJets]/F");
  mytree->Branch("negativeOnlyJetBProbability", &negativeOnlyJetBProbability,"negativeOnlyJetBProbability[nJets]/F");
  mytree->Branch("negativeOnlyJetProbability", &negativeOnlyJetProbability,"negativeOnlyJetProbability[nJets]/F");
  mytree->Branch("negativeTrackCountingHighEff", &negativeTrackCountingHighEff,"negativeTrackCountingHighEff[nJets]/F");
  mytree->Branch("negativeTrackCountingHighPur", &negativeTrackCountingHighPur,"negativeTrackCountingHighPur[nJets]/F");
  mytree->Branch("positiveOnlyJetBProbability", &positiveOnlyJetBProbability,"positiveOnlyJetBProbability[nJets]/F");
  mytree->Branch("positiveOnlyJetProbability", &positiveOnlyJetProbability,"positiveOnlyJetProbability[nJets]/F");
  mytree->Branch("simpleSecondaryVertexHighEff", &simpleSecondaryVertexHighEff,"simpleSecondaryVertexHighEff[nJets]/F");
  mytree->Branch("simpleSecondaryVertexHighPur", &simpleSecondaryVertexHighPur,"simpleSecondaryVertexHighPur[nJets]/F");
  mytree->Branch("combinedSecondaryVertex", &combinedSecondaryVertex,"combinedSecondaryVertex[nJets]/F");
  mytree->Branch("combinedSecondaryVertexPositive", &combinedSecondaryVertexPositive,"combinedSecondaryVertexPositive[nJets]/F");
  mytree->Branch("combinedSecondaryVertexMVA", &combinedSecondaryVertexMVA,"combinedSecondaryVertexMVA[nJets]/F");
  //Zest vars
  mytree->Branch("jetZest",&jetZest,"jetZest[nJets]/F");
  mytree->Branch("jetPPerp",&jetPPerp,"jetPPerp[nJets]/F");

  mytree->Branch("nGenPars",&nGenPars,"nGenPars/I");
  //Kinematic gen Particles
  mytree->Branch("genParPt",&genParPt,"genParPt[nGenPars]/F");
  mytree->Branch("genParEta",&genParEta,"genParEta[nGenPars]/F");
  mytree->Branch("genParPhi",&genParPhi,"genParPhi[nGenPars]/F");
  mytree->Branch("genParMass",&genParMass,"genParMass[nGenPars]/F");
  //ID
  mytree->Branch("genParID",&genParID,"genParID[nGenPars]/I");
  mytree->Branch("genParMotherID",&genParMotherID,"genParMotherID[nGenPars]/I");
  mytree->Branch("genParMotherPz",&genParMotherPz,"genParMotherPz[nGenPars]/F");
}

MyTreeProducer::~MyTreeProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  myFile->cd();
  myHisto->Write();
  mytree->Write();
  myFile->Write();
  myFile->Close();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;



#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif

  bool HasGoodVtx = "False";
  edm::Handle<reco::VertexCollection> offlinePV;

  if(iEvent.eventAuxiliary().isRealData()){
    iEvent.getByLabel("offlinePrimaryVertices", offlinePV);
    for ( reco::VertexCollection::const_iterator pv = offlinePV->begin(); pv != offlinePV->end(); ++pv ) {
      if(!(pv->isFake()) && pv->ndof() > 4 && fabs(pv->z()) <= 24.0 && pv->position().Rho() <=2.0){
	HasGoodVtx = "True";
	break;
      }
    }
    // std::cout<<" pv.size() "<<offlinePV->size()<<endl;
  }
  else
    {
      HasGoodVtx = "True";
      nGenPars = 0;
      edm::Handle<reco::GenParticleCollection> genParticles;
      iEvent.getByLabel("genParticles", genParticles);
      for ( reco::GenParticleCollection::const_iterator gp = genParticles->begin(); gp != genParticles->end(); ++gp ){
	if(gp->status() == 3){
	  //cout<<nGenPars<<" pdgID "<<gp->pdgId()<<" status "<<gp->status()<<endl;//
	  cout<< gp->p4()<<endl;
	  genParPt[nGenPars] = gp->pt();
	  genParEta[nGenPars] = gp->eta();
	  genParPhi[nGenPars] = gp->phi();
	  genParMass[nGenPars] = gp->mass();
	  genParID[nGenPars] = gp->pdgId();
	  if(gp->mother() != 0)
	    {//cout<<" motherID "<<gp->mother()->pdgId()<<" motherPz "<<gp->mother()->p4().pz()<<endl;
	      genParMotherID[nGenPars] = gp->mother()->pdgId();
	      genParMotherPz[nGenPars] = gp->mother()->p4().pz();}
	  else{
	    genParMotherID[nGenPars] = 0;
	    genParMotherPz[nGenPars] = 0;
	  }

	  nGenPars++;
	}
      }
    }

  if(!HasGoodVtx){
    for ( reco::VertexCollection::const_iterator pv = offlinePV->begin(); pv != offlinePV->end(); ++pv ) 
      std::cout<<"  pv.SumPt "<<pv->p4().pt()<<"  pv.isFake "<<pv->isFake()<<"  pv.ndof "<<pv->ndof()<<"  pv.abs(z) "<<fabs(pv->z())<<"  pv.Rho "<<pv->position().Rho()<<std::endl;
  }

  int i = 0;
  edm::Handle<pat::ElectronCollection> eleHandle;
  iEvent.getByLabel("selectedPatElectrons", eleHandle);
  for ( pat::ElectronCollection::const_iterator ele = eleHandle->begin(); ele != eleHandle->end(); ++ele ) {
    double pt = ele->pt();
    //     if( ele->Lepton().puChargedHadronIso() != -1){
    //     if(ele->userIsolation(pat::PfPUChargedHadronIso) != -1){
    //std::cout<<" chargedHadronIso "<<ele->pfIsolationVariables().chargedHadronIso<<" neutralHadronIso "<<ele->pfIsolationVariables().neutralHadronIso<<" photonIso "<<ele->pfIsolationVariables().photonIso<<" puChargedHadronIso "<<ele->userIsolation(pat::User4Iso)<<" puChargedHadronIso "<<ele->userIsolation("User4Iso")<<" puChargedHadronIso "<<ele->userIsolation(pat::User5Iso)<<" puChargedHadronIso "<<ele->userIsolation("User5Iso")<<" puChargedHadronIso "<<ele->userIsolation(pat::PfAllParticleIso)<<" puChargedHadronIso "<<ele->userIsolation(pat::PfGammaIso)<<endl;
    std::cout<<" ele No. "<<i++<<" pt "<<pt<<" eta "<<ele->eta()<<" new relIso (trk + calo)/pt "<<(ele->trackIso() + ele->caloIso())/pt<<endl;
    //std::cout<<" ele No. "<<i++<<" pt "<<pt<<" eta "<<ele->eta()<<" isoRel "<<(ele->pfIsolationVariables().chargedHadronIso+ max(ele->pfIsolationVariables().neutralHadronIso + ele->pfIsolationVariables().photonIso - 0.5 * ele->userIsolation(pat::PfPUChargedHadronIso),0.0))/pt<<" id "<<ele->electronID("eidRobustLoose")<<std::endl;
     
     
  }//}



  i = 0;
  edm::Handle<pat::MuonCollection> muHandle;
  iEvent.getByLabel("selectedPatMuons", muHandle);
  for ( pat::MuonCollection::const_iterator mu = muHandle->begin(); mu != muHandle->end(); ++mu ) {
    double pt = mu->pt();
    //     if( ele->Lepton().puChargedHadronIso() != -1){                                                                                                                        
    //     if(ele->userIsolation(pat::PfPUChargedHadronIso) != -1){                                                                                                              
    //std::cout<<" mu No. "<<i<<" pt "<<pt<<" eta "<<mu->eta()<<" isPfIsoValid "<<mu->isPFIsolationValid()<<" isIsoValid "<<mu->isIsolationValid()<<" chargedHadronIso "<<mu->pfIsolationR04().sumChargedHadronPt<<" neutralHadronIso "<<mu->pfIsolationR04().sumNeutralHadronEt<<" photonIso "<<mu->pfIsolationR04().sumPhotonEt<<" puChargedHadronIso "<<mu->pfIsolationR04().sumPUPt<<endl;
    std::cout<<" mu No. "<<i++<<" pt "<<pt<<" eta "<<mu->eta()<<" new relIso (trk + calo)/pt "<<(mu->trackIso() + mu->caloIso())/pt<<endl;
    //std::cout<<" isoRel "<<(mu->pfIsolationR04().sumChargedHadronPt+ max(mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5 * mu->pfIsolationR04().sumPUPt,0.0))/pt<<std::endl;
  }//}                                                                                                                                                                           


  /* 
     int jetNo = 0;
    
     edm::Handle<reco::PFJetCollection> pfjetH;
     iEvent.getByLabel("ak5PFJets", pfjetH);

     for ( reco::PFJetCollection::const_iterator jet = pfjetH->begin(); jet != pfjetH->end(); ++jet ) {
     double pt = jet->pt();
     if(pt < 20)
     break;
     std::cout<<" Jet No. "<<jetNo++<<" pt "<<pt<<std::endl;
     jet->print();
     }

     // Get b tag information
     edm::Handle<reco::JetTagCollection> bTagHandle;
     iEvent.getByLabel("trackCountingHighPurBJetTags", bTagHandle);
     const reco::JetTagCollection & bTags = *(bTagHandle.product());

     // Loop over jets and study b tag info.
     for (unsigned int i = 0; i != bTags.size(); ++i) {
     if(bTags[i].first->pt() < 20)
     break;
     cout<<" Jet "<< i 
     <<" has b tag discriminator = "<<bTags[i].second
     << " and jet Pt = "<<bTags[i].first->pt()<<endl;
     }
  
     jetNo = 0;   */
  nJets = 0;
  edm::Handle<pat::JetCollection> jetHandle;
  iEvent.getByLabel("selectedPatJets", jetHandle);
  for ( pat::JetCollection::const_iterator jet = jetHandle->begin(); jet != jetHandle->end(); ++jet ) {
    double pt = jet->pt();
    if(pt < 20)
      break;
    // cout<<" Jet No. "<<jetNo++<<" pt "<<pt<<" bDisCombi "<<jet->bDiscriminator("combinedSecondaryVertexBJetTags")<<" isPF "<<jet->isPFJet()<<" Calo "<<jet->isCaloJet()<<endl;

    if(pt > 20){
      //Kinematic
      jetPt[nJets] = pt;
      jetEta[nJets] = jet->eta();
      jetPhi[nJets] = jet->phi();
      jetMass[nJets] = jet->mass();
      //BTag Info
      jetBProbability[nJets] = jet->bDiscriminator("jetBProbabilityBJetTags");
      jetProbability[nJets] = jet->bDiscriminator("jetProbabilityBJetTags");
      trackCountingHighPur[nJets] = jet->bDiscriminator("trackCountingHighPurBJetTags");
      trackCountingHighEff[nJets] = jet->bDiscriminator("trackCountingHighEffBJetTags");
      negativeOnlyJetBProbability[nJets] = jet->bDiscriminator("negativeOnlyJetBProbabilityJetTags");
      negativeOnlyJetProbability[nJets] = jet->bDiscriminator("negativeOnlyJetProbabilityJetTags");
      negativeTrackCountingHighEff[nJets] = jet->bDiscriminator("negativeTrackCountingHighEffJetTags");
      negativeTrackCountingHighPur[nJets] = jet->bDiscriminator("negativeTrackCountingHighPurJetTags");
      positiveOnlyJetBProbability[nJets] = jet->bDiscriminator("positiveOnlyJetBProbabilityJetTags");
      positiveOnlyJetProbability[nJets] = jet->bDiscriminator("positiveOnlyJetProbabilityJetTags");
      simpleSecondaryVertexHighEff[nJets] = jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
      simpleSecondaryVertexHighPur[nJets] = jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
      combinedSecondaryVertex[nJets] = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
      combinedSecondaryVertexPositive[nJets] = jet->bDiscriminator("combinedSecondaryVertexPositiveBJetTags");
      combinedSecondaryVertexMVA[nJets] = jet->bDiscriminator("combinedSecondaryVertexMVABJetTags");

      math::XYZVector jetVec = jet->momentum();
      float pPerp = 0.0;
      float zest = 0.0; //1706.03904
      //      float bib  = 0.0; //1706.03904 boost invariant boadening 
      std::vector<reco::PFCandidatePtr> pfCandidates = jet->getPFConstituents();
      for( std::vector<reco::PFCandidatePtr>::const_iterator pfCand = pfCandidates.begin(); pfCand != pfCandidates.end(); ++pfCand){
	
	//	math::XYZVector pfCandVec = (*pfCand)->momentum();
	float pPerpCand = ROOT::Math::VectorUtil::Perp((*pfCand)->momentum(),jetVec);
	pPerp += pPerpCand;
      }
      
      for( std::vector<reco::PFCandidatePtr>::const_iterator pfCand = pfCandidates.begin(); pfCand != pfCandidates.end(); ++pfCand){
	float pPerpCand = ROOT::Math::VectorUtil::Perp((*pfCand)->momentum(),jetVec);
	zest += exp(-pPerp/pPerpCand);
	//	bib = pPerp/jet->mass();
      }
      zest = -1.0/log(zest);
      //      myHisto->Fill(bib);
      
      jetZest[nJets] = zest;
      jetPPerp[nJets] = pPerp;
  
      nJets++;
    }
  }

  //jetNo = 0; 
  // edm::Handle<pat::JetCollection> myJetHandle;
  // iEvent.getByLabel("myPatJets", myJetHandle);
  // for ( pat::JetCollection::const_iterator jet = myJetHandle->begin(); jet != myJetHandle->end(); ++jet ) {
  //   double pt = jet->pt();
  //   if(pt < 20)
  //     break;
  //   cout<<" Jet No. "<<jetNo++<<" pt "<<pt<<" bDisCombi "<<jet->bDiscriminator("combinedSecondaryVertexBJetTags")<<" isPF "<<jet->isPFJet()<<" Calo "<<jet->isCaloJet()<<endl;
  //     //process.patJets.discriminatorSources = cms.VInputTag(cms.InputTag("combinedSecondaryVertexBJetTags"), cms.InputTag("combinedSecondaryVertexMVABJetTags"), cms.InputTag("jetBProbabilityBJetTags"), cms.InputTag("jetProbabilityBJetTags"), cms.InputTag("simpleSecondaryVertexHighEffBJetTags"), cms.InputTag("simpleSecondaryVertexHighPurBJetTags"), cms.InputTag("trackCountingHighEffBJetTags"), cms.InputTag("trackCountingHighPurBJetTags"), cms.InputTag("softPFElectronBJetTags"), cms.InputTag("ghostTrackBJetTags"), cms.InputTag("softPFMuonBJetTags"))
  // }

  mytree->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void 
MyTreeProducer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyTreeProducer::endJob() 
{

}

// ------------ method called when starting to processes a run  ------------
void 
MyTreeProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MyTreeProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MyTreeProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MyTreeProducer::endLuminosityBlock(const edm::LuminosityBlock & lumi, const edm::EventSetup & iSetup)
{// Total number of events is the sum of the events in each of these luminosity blocks
  edm::Handle<edm::MergeableCounter> nEventsTotalCounter;
  lumi.getByLabel("nEventsTotal", nEventsTotalCounter);
  int nEventsTotalCntr = nEventsTotalCounter->value;
  std::cout<<" nEventsTotal "<<nEventsTotalCntr<<std::endl;

  edm::Handle<edm::MergeableCounter> nEventsFilteredCounter;
  lumi.getByLabel("nEventsFiltered", nEventsFilteredCounter);
  int nEventsFilteredCntr = nEventsFilteredCounter->value;
  std::cout<<" nEventsFiltered "<<nEventsFilteredCntr<<std::endl;

  myHisto->Fill(0.5, nEventsTotalCntr);
  myHisto->Fill(1.5, nEventsFilteredCntr);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyTreeProducer);
