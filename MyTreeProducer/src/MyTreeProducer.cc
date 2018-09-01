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
#include "DataFormats/TrackReco/interface/Track.h"
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

  // static const Int_t maxObjects = 30;
  static const Int_t maxObjects = 1000;
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
  Float_t genTrkPt[maxObjects];
  Float_t genTrkEta[maxObjects];
  Float_t genTrkPhi[maxObjects];
  Float_t genTrkMass[maxObjects];
  //ID
  Int_t genTrkID[maxObjects];

  //Kinematic
  Float_t jetPt[maxObjects];
  Float_t jetEta[maxObjects];
  Float_t jetPhi[maxObjects];
  Float_t jetMass[maxObjects];
  //BTag Info
  /*  Float_t jetBProbability[maxObjects];
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
  Float_t combinedSecondaryVertexMVA[maxObjects];*/
  /*//Zest
  Float_t jetZest[maxObjects];
  Float_t jetPPerp[maxObjects];
  //JetPFConstituents
  Int_t nAllJetsContstituents;//number of all jetConstituents in the event
  Int_t nJetContstituents[maxObjects];//number of jetConstituents for a special jet
  Int_t indJetContstituent0[maxObjects];//number of total jetConstituents before this special jet
  Float_t jetContEta[maxObjects*40];
  Float_t jetContPhi[maxObjects*40];
  Float_t jetContPt[maxObjects*40];
  Float_t jetContMass[maxObjects*40];
  Int_t jetContChargex3[maxObjects*40];
*/
  //charged tracks and muon
  Float_t trkPt[maxObjects];
  Float_t trkEta[maxObjects];
  Float_t trkPhi[maxObjects];
  Float_t trkCharge[maxObjects];

  Float_t muPt[maxObjects];
  Float_t muEta[maxObjects];
  Float_t muPhi[maxObjects];
  Float_t muCharge[maxObjects];

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
  isReadFromEOS = iConfig.getParameter<bool>("readFromEOS");

  myFile = new TFile(myHistoName.c_str(),"RECREATE");
  myHisto = new TH1F("myHisto","myHisto", 2, 0, 2);
  //  myHisto2 = new TH1F("myHisto2","myHisto2", 100000, 0, 100);
  mytree  = new TTree("JetInfo", "Jet Information");
  mytree->Branch("nJets",&nJets,"nJets/I");
  //Kinematic
  mytree->Branch("jetPt",&jetPt,"jetPt[nJets]/F");
  mytree->Branch("jetEta",&jetEta,"jetEta[nJets]/F");
  mytree->Branch("jetPhi",&jetPhi,"jetPhi[nJets]/F");
  mytree->Branch("jetMass",&jetMass,"jetMass[nJets]/F");

  //JetPFConstituents
  /*  mytree->Branch("nAllJetsContstituents",&nAllJetsContstituents,"nAllJetsContstituents/I");//number of all jetConstituents in the event
  mytree->Branch("nJetContstituents",&nJetContstituents,"nJetContstituents[nJets]/I");//number of jetConstituents for a special jet
  mytree->Branch("indJetContstituent0",&indJetContstituent0,"indJetContstituent0[nJets]/I");//number of total jetConstituents before this special jet
  mytree->Branch("jetContEta",&jetContEta,"jetContEta[nAllJetsContstituents]/F");
  mytree->Branch("jetContPhi",&jetContPhi,"jetContPhi[nAllJetsContstituents]/F");
  mytree->Branch("jetContPt",&jetContPt,"jetContPt[nAllJetsContstituents]/F");
  mytree->Branch("jetContMass",&jetContMass,"jetContMass[nAllJetsContstituents]/F");
  mytree->Branch("jetContChargex3",&jetContChargex3,"jetContChargex3[nAllJetsContstituents]/I"); //charge times 3 known as qx3
 
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
*/
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

  mytree->Branch("nGenTrks",&nGenTrks,"nGenTrks/I");
  //Kinematic gen Particles
  mytree->Branch("genTrkPt",&genTrkPt,"genTrkPt[nGenTrks]/F");
  mytree->Branch("genTrkEta",&genTrkEta,"genTrkEta[nGenTrks]/F");
  mytree->Branch("genTrkPhi",&genTrkPhi,"genTrkPhi[nGenTrks]/F");
  mytree->Branch("genTrkMass",&genTrkMass,"genTrkMass[nGenTrks]/F");
  //ID
  mytree->Branch("genTrkID",&genTrkID,"genTrkID[nGenTrks]/I");

  mytree->Branch("nTrks",&nTrks,"nTrks/I");
  //Kinematic Tracks
  mytree->Branch("trkPt",&trkPt,"trkPt[nTrks]/F");
  mytree->Branch("trkEta",&trkEta,"trkEta[nTrks]/F");
  mytree->Branch("trkPhi",&trkPhi,"trkPhi[nTrks]/F");
  mytree->Branch("trkCharge",&trkCharge,"trkCharge[nTrks]/F");

  mytree->Branch("nMus",&nMus,"nMus/I");
  //Kinematic Muons
  mytree->Branch("muPt",&muPt,"muPt[nMus]/F");
  mytree->Branch("muEta",&muEta,"muEta[nMus]/F");
  mytree->Branch("muPhi",&muPhi,"muPhi[nMus]/F");
  mytree->Branch("muCharge",&muCharge,"muCharge[nMus]/F");


};

MyTreeProducer::~MyTreeProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  myFile->cd();
  myHisto->Write();
  //  myHisto2->Write();
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

  //if(iEvent.id().event()  == 283642){247684

  reco::GenParticleCollection::const_iterator gp;
  pat::MuonCollection::const_iterator mu0, mu1, mu;
  reco::VertexCollection::const_iterator pv;
  pat::JetCollection::const_iterator jet;
  reco::TrackCollection::const_iterator trk;

 nJets = 0;
 nTrks = 0;
 nMus = 0;
 nGenPars = 0;
 nGenTrks = 0;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif

  int i = 0;
  edm::Handle<pat::ElectronCollection> eleHandle;
  iEvent.getByLabel("cleanPatElectrons", eleHandle);
  for ( pat::ElectronCollection::const_iterator ele = eleHandle->begin(); ele != eleHandle->end(); ele++ ) {
    double pt = ele->pt();
    //     if( ele->Lepton().puChargedHadronIso() != -1){
    //     if(ele->userIsolation(pat::PfPUChargedHadronIso) != -1){
    //std::cout<<" chargedHadronIso "<<ele->pfIsolationVariables().chargedHadronIso<<" neutralHadronIso "<<ele->pfIsolationVariables().neutralHadronIso<<" photonIso "<<ele->pfIsolationVariables().photonIso<<" puChargedHadronIso "<<ele->userIsolation(pat::User4Iso)<<" puChargedHadronIso "<<ele->userIsolation("User4Iso")<<" puChargedHadronIso "<<ele->userIsolation(pat::User5Iso)<<" puChargedHadronIso "<<ele->userIsolation("User5Iso")<<" puChargedHadronIso "<<ele->userIsolation(pat::PfAllParticleIso)<<" puChargedHadronIso "<<ele->userIsolation(pat::PfGammaIso)<<endl;
    std::cout<<" ele No. "<<i++<<" pt "<<pt<<" eta "<<ele->eta()<<" new relIso (trk + calo)/pt "<<(ele->trackIso() + ele->caloIso())/pt<<endl;
    //std::cout<<" ele No. "<<i++<<" pt "<<pt<<" eta "<<ele->eta()<<" isoRel "<<(ele->pfIsolationVariables().chargedHadronIso+ max(ele->pfIsolationVariables().neutralHadronIso + ele->pfIsolationVariables().photonIso - 0.5 * ele->userIsolation(pat::PfPUChargedHadronIso),0.0))/pt<<" id "<<ele->electronID("eidRobustLoose")<<std::endl;
     
     
  }//}



  i = 0;
  //  cout<< "Here0 "<<endl;
  edm::Handle<pat::MuonCollection> muHandle;
  iEvent.getByLabel("cleanPatMuons", muHandle);
  // for ( pat::MuonCollection::const_iterator mu = muHandle->begin(); mu != muHandle->end(); ++mu ) {
  //   double pt = mu->pt();
  //   //     if( ele->Lepton().puChargedHadronIso() != -1){                                                                                                                        
  //   //     if(ele->userIsolation(pat::PfPUChargedHadronIso) != -1){                                                                                                              
  //   //std::cout<<" mu No. "<<i<<" pt "<<pt<<" eta "<<mu->eta()<<" isPfIsoValid "<<mu->isPFIsolationValid()<<" isIsoValid "<<mu->isIsolationValid()<<" chargedHadronIso "<<mu->pfIsolationR04().sumChargedHadronPt<<" neutralHadronIso "<<mu->pfIsolationR04().sumNeutralHadronEt<<" photonIso "<<mu->pfIsolationR04().sumPhotonEt<<" puChargedHadronIso "<<mu->pfIsolationR04().sumPUPt<<endl;
  //   std::cout<<" mu No. "<<i++<<" pt "<<pt<<" eta "<<mu->eta()<<" new relIso (trk + calo)/pt "<<(mu->trackIso() + mu->caloIso())/pt<<endl;
  //   //std::cout<<" isoRel "<<(mu->pfIsolationR04().sumChargedHadronPt+ max(mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5 * mu->pfIsolationR04().sumPUPt,0.0))/pt<<std::endl;
  // }//}     
  //cout<< "Here00 "<<endl;
  bool badEvent = true;
  if(muHandle->size() > 1){
     mu0 = muHandle->begin();
     mu1 = mu0+1;
     //     cout<<" mu0->charge() "<<endl;
     if(mu0->charge() != mu1->charge()){
       float mass = (mu0->p4() + mu1->p4()).mass();
       //      myHisto->Fill(mass);
       if(fabs(mass - 90) < 25)
   	badEvent = false;
     }
  }

  //badEvent = false;

  //cout<< "Here000 "<<endl;
  if(!badEvent){
    reco::TrackBase::Point muVertex = muHandle->begin()->vertex();
    
    for ( mu = muHandle->begin(); mu != muHandle->end(); mu++ ) {
      muPt[nMus] = mu->pt();
      muPhi[nMus] = mu->phi();
      muEta[nMus] = mu->eta();
      muCharge[nMus] = mu->charge();
      nMus++;
    }
    //cout<< "Here0001 "<<endl;
  edm::Handle<reco::VertexCollection> offlinePV;
  iEvent.getByLabel("offlinePrimaryVertices", offlinePV);
  /*  float vx = (offlinePV->begin())->position().x();
  float vy = (offlinePV->begin())->position().y();
  float vz = (offlinePV->begin())->position().z();
  */
  // myHisto2->Fill(fabs(vx - muVertex.x()));
  // myHisto2->Fill(fabs(vy - muVertex.y()) + 10);
  // myHisto2->Fill(fabs(vz - muVertex.z()) + 20);
  // cout<< "Here00011 "<<endl;

  // cout<<" vx "<<vx<<" vy "<<vy<<endl;
  if(iEvent.eventAuxiliary().isRealData()){
    
    for ( pv = offlinePV->begin(); pv != offlinePV->end(); pv++ ) {
      if(!(pv->isFake()) && pv->ndof() > 4 && fabs(pv->z()) <= 24.0 && pv->position().Rho() <=2.0){
	break;
      }
    }
    // std::cout<<" pv.size() "<<offlinePV->size()<<endl;
  }
  else
    {// cout<< "Here000111 "<<endl;
      int ii = 0;
      float genVtxX = 0.0;
      float genVtxY = 0.0;
      float genVtxZ = 0.0;

      edm::Handle<reco::GenParticleCollection> genParticles;
      iEvent.getByLabel("genParticles", genParticles);
      // for ( reco::GenParticleCollection::const_iterator gp = genParticles->begin(); gp != genParticles->end(); ++gp ){
      for ( gp = genParticles->begin(); gp != genParticles->end(); gp++ ){
	ii++;
	if(gp->status() == 3){
	  if(ii == 5){
	    genVtxX = gp->vertex().x();
	    genVtxY = gp->vertex().y();
	    genVtxZ = gp->vertex().z();    
	  }
	  // cout<<" genVtxX "<<genVtxX<<" genVtxY "<<genVtxY<<" genVtxZ "<<genVtxZ<<endl;
	  // cout<<ii<<" pdgID "<<gp->pdgId()<<" status "<<gp->status()<<" pz "<<gp->p4().pz()<<" gp->vertex () "<<gp->vertex();
	  // cout<< "Here00011123 "<<gp->mother()<<endl;
	  // if(gp->mother() != 0) cout<<" motherID "<<gp->mother()->pdgId()<<" motherPz "<<gp->mother()->p4().pz()<<endl;
	  // else 
	  //  cout<<endl;
	  //cout<< gp->p4()<<endl;
	  genParPt[nGenPars] = gp->pt();
	  genParEta[nGenPars] = gp->eta();
	  genParPhi[nGenPars] = gp->phi();
	  genParMass[nGenPars] = gp->mass();
	  genParID[nGenPars] = gp->pdgId();
	  //cout<< "Here00011123222 "<<endl;
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
	//	 cout<< "Here00014234 "<<endl;
	 if(gp->status() == 1 && gp->charge() != 0 &&  gp->pt() > 0.5 && fabs(gp->eta()) < 2.5){

	   // myHisto2->Fill(fabs(genVtxX - gp->vertex().x()) + 30);
	   // myHisto2->Fill(fabs(genVtxY - gp->vertex().y()) + 50);
	   // myHisto2->Fill(fabs(genVtxZ - gp->vertex().z()) + 70);

	   if(fabs(genVtxX - gp->vertex().x()) < 0.001 && fabs(genVtxY - gp->vertex().y()) < 0.001 && fabs(genVtxZ - gp->vertex().z()) < 0.001){
	     //   cout<<" genVtxX "<<genVtxX<<" genVtxY "<<genVtxY<<" genVtxZ "<<genVtxZ<<endl;
	     //   cout<<ii<<" pdgID "<<gp->pdgId()<<" status "<<gp->status()<<" pz "<<gp->p4().pz()<<" gp->vertex () "<<gp->vertex()<<endl;
	     // //if(genVertex != gp->vertex())
	     // //  cout<< gp->p4()<<endl;
	     // }
	     
	     genTrkPt[nGenTrks] = gp->pt();
	     genTrkEta[nGenTrks] = gp->eta();
	     genTrkPhi[nGenTrks] = gp->phi();
	     genTrkMass[nGenTrks] = gp->mass();
	     genTrkID[nGenTrks] = gp->pdgId();
	     nGenTrks++;
	   }// if(fabs(genVtxX - gp->vertex
	 }// if(gp->status() == 1 && 
      }// for ( reco::GenParticleCollectio
    }//else
  // for(int igen = 0; igen < nGenPars; igen++){
  //   cout<< igen <<" ID "<<genParID[igen]<<" Mother "<<genParMotherID[igen]<<" Pt "<<genParPt[igen]<<" Phi "<<genParPhi[igen]<<" Eta "<<genParEta[igen]<<endl;
  // }


  //int jetNo = 0;
     /* 
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
  //nAllJetsContstituents = 0;
  
  edm::Handle<pat::JetCollection> jetHandle;
  iEvent.getByLabel("selectedPatJets", jetHandle);
  for ( jet = jetHandle->begin(); jet != jetHandle->end(); jet++ ) {
    double pt = jet->pt();

    if(pt > 20){
      //int nJetConts = jet->getPFConstituents().size();
      //cout<<" Jet No. "<<jetNo++<<" pt "<<pt<<" bDisCombi "<<jet->bDiscriminator("combinedSecondaryVertexBJetTags")<<" isPF "<<jet->isPFJet()<<" Calo "<<jet->isCaloJet()<<" nConstituents "<<nJetConts<<endl;
      
      //Kinematic
      jetPt[nJets] = pt;
      jetEta[nJets] = jet->eta();
      jetPhi[nJets] = jet->phi();
      jetMass[nJets] = jet->mass();
      //cout<<nJets<<" jVx "<<jet->vx()<<endl;
      /* //BTag Info
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
      //JetPFConstituents
      indJetContstituent0[nJets] = nAllJetsContstituents;

     

      math::XYZVector jetVec = jet->momentum();
      float pPerp = 0.0;
      float zest = 0.0; //1706.03904
      //      float bib  = 0.0; //1706.03904 boost invariant boadening 
      std::vector<reco::PFCandidatePtr> pfCandidates = jet->getPFConstituents();
      for( std::vector<reco::PFCandidatePtr>::const_iterator pfCand = pfCandidates.begin(); pfCand != pfCandidates.end(); ++pfCand){

	//cout<<nJets<<" jCVx "<<(*pfCand)->vx()<<" jCVy "<<(*pfCand)->vy()<<" 3q "<<(*pfCand)->threeCharge()<<endl;
	
	//	math::XYZVector pfCandVec = (*pfCand)->momentum();
	
	jetContEta[nAllJetsContstituents] = (*pfCand)->eta();
	jetContPhi[nAllJetsContstituents] = (*pfCand)->phi();
	jetContPt[nAllJetsContstituents] = (*pfCand)->pt();
	jetContMass[nAllJetsContstituents] = (*pfCand)->mass();
	jetContChargex3[nAllJetsContstituents] = (*pfCand)->threeCharge();
	nAllJetsContstituents++;



	if((*pfCand)->threeCharge() != 0 && (fabs((*pfCand)->vx() - vx) > 0.1 || fabs((*pfCand)->vy() - vy) > 0.1 )) continue;
	if((*pfCand)->threeCharge() == 0 && (*pfCand)->pt() < 1) continue;
   
	float pPerpCand = ROOT::Math::VectorUtil::Perp((*pfCand)->momentum(),jetVec);
	pPerp += pPerpCand;

	//	cout<<" jetContEta[jetCon] "<<jetContEta[jetCon]<<endl;
      }
      nJetContstituents[nJets] = nJetConts;
      //cout<<" nAllJetsContstituents "<<nAllJetsContstituents<<" indJetContstituent0[nJets] "<<indJetContstituent0[nJets]<<" nJetContstituents[nJets] "<<nJetContstituents[nJets]<<endl;
      
      for( std::vector<reco::PFCandidatePtr>::const_iterator pfCand = pfCandidates.begin(); pfCand != pfCandidates.end(); ++pfCand){
	if((*pfCand)->threeCharge() != 0 && (fabs((*pfCand)->vx() - vx) > 0.1 || fabs((*pfCand)->vy() - vy) > 0.1 )) continue;
	if((*pfCand)->threeCharge() == 0 && (*pfCand)->pt() < 1) continue;
	
	float pPerpCand = ROOT::Math::VectorUtil::Perp((*pfCand)->momentum(),jetVec);
	zest += exp(-pPerp/pPerpCand);
	//	bib = pPerp/jet->mass();
      }
      zest = -1.0/log(zest);
      //      myHisto->Fill(bib);
      
      jetZest[nJets] = zest;
      jetPPerp[nJets] = pPerp;
      */
      nJets++;
    }    
  }//cout<< "Here00015675 "<<endl;
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

  reco::TrackBase::Point vert(muVertex.x(), muVertex.y(), muVertex.z());
  //reco::TrackBase::Point vert(vx, vy, vz);
  //cout<<" recVertex "<<muVertex<<" offPrimaryVtx "<<vert<<endl;
  //cout<< "Here "<<endl;
  edm::Handle<reco::TrackCollection> trkHandle;
  //cout<< "Here222222 "<<endl;
  iEvent.getByLabel("generalTracks", trkHandle);
  //cout<< "Here2 "<<endl;
  for (trk = trkHandle->begin(); trk != trkHandle->end(); trk++ ) {
    //cout<< "Here3 "<<endl;
    //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookTrackAnalysis
    if(trk->pt()<=0.5)continue;
    if(trk->charge()==0) continue;// NO neutral objects
    //    if(fabs(trk->pdgId())!=211) continue;//Due to the lack of the particle ID all the tracks for cms are pions(ID==211)
    //    if(!(trk->trackHighPurity())) continue; 
    if(fabs(trk->eta()) > 2.5)continue;
    
    //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTransientTracks#Examples_including_calculation_o
    double dxy = trk->dxy(vert);
    double dxy_error = trk->dxyError();
    double z0 = trk->dz(vert);
    double z0_error = trk->dzError();

    //Conditions from 1204.1411 CMS 7 TeV
    if(dxy_error == 0 || z0_error == 0)
      continue;
    
    if(fabs(dxy/dxy_error) > 3 || fabs(z0/z0_error) > 3) 
      continue;

    if(fabs(trk->ptError()/trk->pt()) > 0.05)
      continue;
    
    if (trk->quality(reco::TrackBase::highPurity)){
      trkPt[nTrks] = trk->pt();
      trkPhi[nTrks] = trk->phi();
      trkEta[nTrks] = trk->eta();
      trkCharge[nTrks] = trk->charge();
      
      nTrks++;
    }// cout<< "Here33 "<<endl;
  } //cout<< "Here222 "<<endl;
  }// if(!badEvent){
  
  mytree->Fill();
  //}//if(iEvent.id() == 327470){
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
MyTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyTreeProducer);
