//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Oct  8 15:26:47 2017 by ROOT version 5.32/00
// from TTree JetInfo/Jet Information
// found on file: Histo_DYJets_0028_2.root
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
   Float_t         jetPt[5];   //[nJets]
   Float_t         jetEta[5];   //[nJets]
   Float_t         jetPhi[5];   //[nJets]
   Float_t         jetMass[5];   //[nJets]
   Float_t         jetBProbability[5];   //[nJets]
   Float_t         jetProbability[5];   //[nJets]
   Float_t         trackCountingHighPur[5];   //[nJets]
   Float_t         trackCountingHighEff[5];   //[nJets]
   Float_t         negativeOnlyJetBProbability[5];   //[nJets]
   Float_t         negativeOnlyJetProbability[5];   //[nJets]
   Float_t         negativeTrackCountingHighEff[5];   //[nJets]
   Float_t         negativeTrackCountingHighPur[5];   //[nJets]
   Float_t         positiveOnlyJetBProbability[5];   //[nJets]
   Float_t         positiveOnlyJetProbability[5];   //[nJets]
   Float_t         simpleSecondaryVertexHighEff[5];   //[nJets]
   Float_t         simpleSecondaryVertexHighPur[5];   //[nJets]
   Float_t         combinedSecondaryVertex[5];   //[nJets]
   Float_t         combinedSecondaryVertexPositive[5];   //[nJets]
   Float_t         combinedSecondaryVertexMVA[5];   //[nJets]
   Float_t         jetZest[5];   //[nJets]
   Float_t         jetPPerp[5];   //[nJets]
   Int_t           nGenPars;
   Float_t         genParPt[12];   //[nGenPars]
   Float_t         genParEta[12];   //[nGenPars]
   Float_t         genParPhi[12];   //[nGenPars]
   Float_t         genParMass[12];   //[nGenPars]
   Int_t           genParID[12];   //[nGenPars]
   Int_t           genParMotherID[12];   //[nGenPars]
   Float_t         genParMotherPz[12];   //[nGenPars]

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
   TBranch        *b_nGenPars;   //!
   TBranch        *b_genParPt;   //!
   TBranch        *b_genParEta;   //!
   TBranch        *b_genParPhi;   //!
   TBranch        *b_genParMass;   //!
   TBranch        *b_genParID;   //!
   TBranch        *b_genParMotherID;   //!
   TBranch        *b_genParMotherPz;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Histo_DYJets_0028_2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Histo_DYJets_0028_2.root");
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
   fChain->SetBranchAddress("nGenPars", &nGenPars, &b_nGenPars);
   fChain->SetBranchAddress("genParPt", genParPt, &b_genParPt);
   fChain->SetBranchAddress("genParEta", genParEta, &b_genParEta);
   fChain->SetBranchAddress("genParPhi", genParPhi, &b_genParPhi);
   fChain->SetBranchAddress("genParMass", genParMass, &b_genParMass);
   fChain->SetBranchAddress("genParID", genParID, &b_genParID);
   fChain->SetBranchAddress("genParMotherID", genParMotherID, &b_genParMotherID);
   fChain->SetBranchAddress("genParMotherPz", genParMotherPz, &b_genParMotherPz);
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
