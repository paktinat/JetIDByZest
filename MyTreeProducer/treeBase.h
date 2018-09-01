//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Sep 16 12:05:52 2017 by ROOT version 5.32/00
// from TTree JetInfo/Jet Information
// found on file: /cmsdata2/paktinat/JetIDByZest/Histo_0004_New.root
//////////////////////////////////////////////////////////

#ifndef treeBase_h
#define treeBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class treeBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nJets;
   Float_t         jetPt[7];   //[nJets]
   Float_t         jetEta[7];   //[nJets]
   Float_t         jetPhi[7];   //[nJets]
   Float_t         jetMass[7];   //[nJets]
   Float_t         jetBProbability[7];   //[nJets]
   Float_t         jetProbability[7];   //[nJets]
   Float_t         trackCountingHighPur[7];   //[nJets]
   Float_t         trackCountingHighEff[7];   //[nJets]
   Float_t         negativeOnlyJetBProbability[7];   //[nJets]
   Float_t         negativeOnlyJetProbability[7];   //[nJets]
   Float_t         negativeTrackCountingHighEff[7];   //[nJets]
   Float_t         negativeTrackCountingHighPur[7];   //[nJets]
   Float_t         positiveOnlyJetBProbability[7];   //[nJets]
   Float_t         positiveOnlyJetProbability[7];   //[nJets]
   Float_t         simpleSecondaryVertexHighEff[7];   //[nJets]
   Float_t         simpleSecondaryVertexHighPur[7];   //[nJets]
   Float_t         combinedSecondaryVertex[7];   //[nJets]
   Float_t         combinedSecondaryVertexPositive[7];   //[nJets]
   Float_t         combinedSecondaryVertexMVA[7];   //[nJets]
   Float_t         jetZest[7];   //[nJets]
   Float_t         jetPPerp[7];   //[nJets]

   // List of branches
   TBranch        *b_nJets;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetBProbability;   //!
   TBranch        *b_jetProbability;   //!
   TBranch        *b_trackCountingHighPur;   //!
   TBranch        *b_trackCountingHighEff;   //!
   TBranch        *b_negativeOnlyJetBProbability;   //!
   TBranch        *b_negativeOnlyJetProbability;   //!
   TBranch        *b_negativeTrackCountingHighEff;   //!
   TBranch        *b_negativeTrackCountingHighPur;   //!
   TBranch        *b_positiveOnlyJetBProbability;   //!
   TBranch        *b_positiveOnlyJetProbability;   //!
   TBranch        *b_simpleSecondaryVertexHighEff;   //!
   TBranch        *b_simpleSecondaryVertexHighPur;   //!
   TBranch        *b_combinedSecondaryVertex;   //!
   TBranch        *b_combinedSecondaryVertexPositive;   //!
   TBranch        *b_combinedSecondaryVertexMVA;   //!
   TBranch        *b_jetZest;   //!
   TBranch        *b_jetPPerp;   //!

   treeBase(TTree *tree=0);
   virtual ~treeBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef treeBase_cxx
treeBase::treeBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/cmsdata2/paktinat/JetIDByZest/Histo_0004_New.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/cmsdata2/paktinat/JetIDByZest/Histo_0004_New.root");
      }
      f->GetObject("JetInfo",tree);

   }
   Init(tree);
}

treeBase::~treeBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t treeBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t treeBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void treeBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetMass", jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetBProbability", jetBProbability, &b_jetBProbability);
   fChain->SetBranchAddress("jetProbability", jetProbability, &b_jetProbability);
   fChain->SetBranchAddress("trackCountingHighPur", trackCountingHighPur, &b_trackCountingHighPur);
   fChain->SetBranchAddress("trackCountingHighEff", trackCountingHighEff, &b_trackCountingHighEff);
   fChain->SetBranchAddress("negativeOnlyJetBProbability", negativeOnlyJetBProbability, &b_negativeOnlyJetBProbability);
   fChain->SetBranchAddress("negativeOnlyJetProbability", negativeOnlyJetProbability, &b_negativeOnlyJetProbability);
   fChain->SetBranchAddress("negativeTrackCountingHighEff", negativeTrackCountingHighEff, &b_negativeTrackCountingHighEff);
   fChain->SetBranchAddress("negativeTrackCountingHighPur", negativeTrackCountingHighPur, &b_negativeTrackCountingHighPur);
   fChain->SetBranchAddress("positiveOnlyJetBProbability", positiveOnlyJetBProbability, &b_positiveOnlyJetBProbability);
   fChain->SetBranchAddress("positiveOnlyJetProbability", positiveOnlyJetProbability, &b_positiveOnlyJetProbability);
   fChain->SetBranchAddress("simpleSecondaryVertexHighEff", simpleSecondaryVertexHighEff, &b_simpleSecondaryVertexHighEff);
   fChain->SetBranchAddress("simpleSecondaryVertexHighPur", simpleSecondaryVertexHighPur, &b_simpleSecondaryVertexHighPur);
   fChain->SetBranchAddress("combinedSecondaryVertex", combinedSecondaryVertex, &b_combinedSecondaryVertex);
   fChain->SetBranchAddress("combinedSecondaryVertexPositive", combinedSecondaryVertexPositive, &b_combinedSecondaryVertexPositive);
   fChain->SetBranchAddress("combinedSecondaryVertexMVA", combinedSecondaryVertexMVA, &b_combinedSecondaryVertexMVA);
   fChain->SetBranchAddress("jetZest", jetZest, &b_jetZest);
   fChain->SetBranchAddress("jetPPerp", jetPPerp, &b_jetPPerp);
   Notify();
}

Bool_t treeBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void treeBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t treeBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef treeBase_cxx
