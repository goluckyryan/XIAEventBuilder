//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec 13 18:37:46 2021 by ROOT version 6.24/06
// from TTree tree/tree
// found on file: efEu152b-000.root
//////////////////////////////////////////////////////////

#ifndef Analyzer_h
#define Analyzer_h

#include <fstream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include "mapping.h"
#include "armory/AnalysisLibrary.h"

// Header file for the classes stored in the TTree if any.

class Analyzer : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       evID;
   Double_t        e[NCRYSTAL];
   ULong64_t       t[NCRYSTAL];
   UShort_t        p[NCRYSTAL];
   Double_t        bgo[NBGO];
   Double_t        other[NOTHER];
   Int_t           multi;

   // List of branches
   TBranch        *b_event_ID;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_time;   //!
   TBranch        *b_pileup;   //!
   TBranch        *b_bgo;   //!
   TBranch        *b_other;   //!
   TBranch        *b_multiplicity;   //!

   Analyzer(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Analyzer() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(Analyzer,0);
};

#endif

#ifdef Analyzer_cxx
void Analyzer::Init(TTree *tree)
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
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evID",   &evID, &b_event_ID);
   fChain->SetBranchAddress("e",          e, &b_energy);
   fChain->SetBranchAddress("t",          t, &b_time);
   fChain->SetBranchAddress("p",          p, &b_pileup);
   fChain->SetBranchAddress("bgo",      bgo, &b_bgo);
   fChain->SetBranchAddress("other",  other, &b_other);
   fChain->SetBranchAddress("multi", &multi, &b_multiplicity);
}

Bool_t Analyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}



void Analyzer::SlaveBegin(TTree * /*tree*/){

   TString option = GetOption();

}


void Analyzer::SlaveTerminate(){

}

#endif // #ifdef Analyzer_cxx
