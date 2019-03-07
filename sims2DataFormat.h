//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar  6 11:01:20 2019 by ROOT version 6.12/06
// from TTree tree/lab setup entangled sim
// found on file: 2billion.root
//////////////////////////////////////////////////////////

#ifndef sims2DataFormat_h
#define sims2DataFormat_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

// Header file for the classes stored in the TTree if any.

class sims2DataFormat {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         theta1;
   Float_t         phi1;
   Float_t         E1_a;
   Float_t         X1_a;
   Float_t         Y1_a;
   Float_t         Z1_a;
   Float_t         E1_b;
   Float_t         X1_b;
   Float_t         Y1_b;
   Float_t         Z1_b;
   Float_t         theta2;
   Float_t         phi2;
   Float_t         E2_a;
   Float_t         X2_a;
   Float_t         Y2_a;
   Float_t         Z2_a;
   Float_t         E2_b;
   Float_t         X2_b;
   Float_t         Y2_b;
   Float_t         Z2_b;

   // List of branches
   TBranch        *b_theta1;   //!
   TBranch        *b_phi1;   //!
   TBranch        *b_E1_a;   //!
   TBranch        *b_X1_a;   //!
   TBranch        *b_Y1_a;   //!
   TBranch        *b_Z1_a;   //!
   TBranch        *b_E1_b;   //!
   TBranch        *b_X1_b;   //!
   TBranch        *b_Y1_b;   //!
   TBranch        *b_Z1_b;   //!
   TBranch        *b_theta2;   //!
   TBranch        *b_phi2;   //!
   TBranch        *b_E2_a;   //!
   TBranch        *b_X2_a;   //!
   TBranch        *b_Y2_a;   //!
   TBranch        *b_Z2_a;   //!
   TBranch        *b_E2_b;   //!
   TBranch        *b_X2_b;   //!
   TBranch        *b_Y2_b;   //!
   TBranch        *b_Z2_b;   //!

   sims2DataFormat(TTree *tree=0);
   virtual ~sims2DataFormat();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
};

#endif

#ifdef sims2DataFormat_cxx
sims2DataFormat::sims2DataFormat(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("2billion_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("2billion_1.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

sims2DataFormat::~sims2DataFormat()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t sims2DataFormat::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t sims2DataFormat::LoadTree(Long64_t entry)
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

void sims2DataFormat::Init(TTree *tree)
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

   fChain->SetBranchAddress("theta1", &theta1, &b_theta1);
   fChain->SetBranchAddress("phi1", &phi1, &b_phi1);
   fChain->SetBranchAddress("E1_a", &E1_a, &b_E1_a);
   fChain->SetBranchAddress("X1_a", &X1_a, &b_X1_a);
   fChain->SetBranchAddress("Y1_a", &Y1_a, &b_Y1_a);
   fChain->SetBranchAddress("Z1_a", &Z1_a, &b_Z1_a);
   fChain->SetBranchAddress("E1_b", &E1_b, &b_E1_b);
   fChain->SetBranchAddress("X1_b", &X1_b, &b_X1_b);
   fChain->SetBranchAddress("Y1_b", &Y1_b, &b_Y1_b);
   fChain->SetBranchAddress("Z1_b", &Z1_b, &b_Z1_b);
   fChain->SetBranchAddress("theta2", &theta2, &b_theta2);
   fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
   fChain->SetBranchAddress("E2_a", &E2_a, &b_E2_a);
   fChain->SetBranchAddress("X2_a", &X2_a, &b_X2_a);
   fChain->SetBranchAddress("Y2_a", &Y2_a, &b_Y2_a);
   fChain->SetBranchAddress("Z2_a", &Z2_a, &b_Z2_a);
   fChain->SetBranchAddress("E2_b", &E2_b, &b_E2_b);
   fChain->SetBranchAddress("X2_b", &X2_b, &b_X2_b);
   fChain->SetBranchAddress("Y2_b", &Y2_b, &b_Y2_b);
   fChain->SetBranchAddress("Z2_b", &Z2_b, &b_Z2_b);
   Notify();
}

Bool_t sims2DataFormat::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void sims2DataFormat::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t sims2DataFormat::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Int_t findAM(Float_t z);

Int_t findPixelID(Int_t am, Float_t x, Float_t y);

#endif // #ifdef sims2DataFormat_cxx
