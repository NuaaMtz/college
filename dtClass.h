//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov  3 15:07:21 2023 by ROOT version 6.28/06
// from TTree tree/phase space and bec events
// found on file: events.root
//////////////////////////////////////////////////////////

#ifndef dtClass_h
#define dtClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class dtClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           BEC;
   Double_t        eb;
   Double_t        ran;
   Double_t        P1[4];
   Double_t        P2[4];
   Double_t        P3[4];
   Double_t        theta1;
   Double_t        theta2;
   Double_t        theta3;
   Double_t        phi1;
   Double_t        phi2;
   Double_t        phi3;
   Double_t        M12;
   Double_t        M23;
   Double_t        M13;
   Double_t        Q12;
   Double_t        MM;
   Double_t        Pm1[4];
   Double_t        Pm2[4];
   Double_t        Pm3[4];
   Double_t        Pm4[4];
   Double_t        Pm5[4];

   // List of branches
   TBranch        *b_BEC;   //!
   TBranch        *b_eb;   //!
   TBranch        *b_ran;   //!
   TBranch        *b_P1;   //!
   TBranch        *b_P2;   //!
   TBranch        *b_P3;   //!
   TBranch        *b_theta1;   //!
   TBranch        *b_theta2;   //!
   TBranch        *b_theta3;   //!
   TBranch        *b_phi1;   //!
   TBranch        *b_phi2;   //!
   TBranch        *b_phi3;   //!
   TBranch        *b_M12;   //!
   TBranch        *b_M23;   //!
   TBranch        *b_M13;   //!
   TBranch        *b_Q12;   //!
   TBranch        *b_MM;   //!
   TBranch        *b_Pm1;   //!
   TBranch        *b_Pm2;   //!
   TBranch        *b_Pm3;   //!
   TBranch        *b_Pm4;   //!
   TBranch        *b_Pm5;   //!

   dtClass(TTree *tree=0);
   virtual ~dtClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef dtClass_cxx
dtClass::dtClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("events.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("events.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

dtClass::~dtClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t dtClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t dtClass::LoadTree(Long64_t entry)
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

void dtClass::Init(TTree *tree)
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

   fChain->SetBranchAddress("BEC", &BEC, &b_BEC);
   fChain->SetBranchAddress("eb", &eb, &b_eb);
   fChain->SetBranchAddress("ran", &ran, &b_ran);
   fChain->SetBranchAddress("P1", P1, &b_P1);
   fChain->SetBranchAddress("P2", P2, &b_P2);
   fChain->SetBranchAddress("P3", P3, &b_P3);
   fChain->SetBranchAddress("theta1", &theta1, &b_theta1);
   fChain->SetBranchAddress("theta2", &theta2, &b_theta2);
   fChain->SetBranchAddress("theta3", &theta3, &b_theta3);
   fChain->SetBranchAddress("phi1", &phi1, &b_phi1);
   fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
   fChain->SetBranchAddress("phi3", &phi3, &b_phi3);
   fChain->SetBranchAddress("M12", &M12, &b_M12);
   fChain->SetBranchAddress("M23", &M23, &b_M23);
   fChain->SetBranchAddress("M13", &M13, &b_M13);
   fChain->SetBranchAddress("Q12", &Q12, &b_Q12);
   fChain->SetBranchAddress("MM", &MM, &b_MM);
   fChain->SetBranchAddress("Pm1", Pm1, &b_Pm1);
   fChain->SetBranchAddress("Pm2", Pm2, &b_Pm2);
   fChain->SetBranchAddress("Pm3", Pm3, &b_Pm3);
   fChain->SetBranchAddress("Pm4", Pm4, &b_Pm4);
   fChain->SetBranchAddress("Pm5", Pm5, &b_Pm5);
   Notify();
}

Bool_t dtClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void dtClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t dtClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef dtClass_cxx
