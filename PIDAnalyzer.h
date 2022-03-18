//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 18 18:03:07 2022 by ROOT version 6.26/00
// from TTree tree/tree
// found on file: haha.root
//////////////////////////////////////////////////////////

#ifndef PIDAnalyzer_h
#define PIDAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TStopwatch.h>

#include "mapping.h"

// Header file for the classes stored in the TTree if any.

#define MAX_MULTI 40

class PIDAnalyzer : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   ULong64_t       eventID;
   Int_t           runID;
   Int_t           multi;
   Int_t           multiGagg;
   Short_t         gammaID[MAX_MULTI];   //[multi]
   Double_t        gamma[MAX_MULTI];   //[multi]
   ULong64_t       gamma_t[MAX_MULTI];   //[multi]
   Int_t           gaggID[MAX_MULTI];   //[multiGagg]
   Double_t        gaggP[MAX_MULTI];   //[multiGagg]
   Double_t        gaggT[MAX_MULTI];   //[multiGagg]
   ULong64_t       gagg_t[MAX_MULTI];   //[multiGagg]

   // List of branches
   TBranch        *b_eventID;   //!
   TBranch        *b_runID;   //!
   TBranch        *b_multi;   //!
   TBranch        *b_multiGagg;   //!
   TBranch        *b_gammaID;   //!
   TBranch        *b_gamma;   //!
   TBranch        *b_gamma_t;   //!
   TBranch        *b_gaggID;   //!
   TBranch        *b_gaggP;   //!
   TBranch        *b_gaggT;   //!
   TBranch        *b_gagg_t;   //!

   PIDAnalyzer(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~PIDAnalyzer() { }
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

   ClassDef(PIDAnalyzer,0);
   
   ULong64_t totnumEntry;

   ULong64_t ProcessedEntries;
   Float_t Frac; ///Progress bar
   TStopwatch StpWatch;
};

#endif

#ifdef PIDAnalyzer_cxx
void PIDAnalyzer::Init(TTree *tree)
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

   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("runID", &runID, &b_runID);
   fChain->SetBranchAddress("multi", &multi, &b_multi);
   fChain->SetBranchAddress("multiGagg", &multiGagg, &b_multiGagg);
   fChain->SetBranchAddress("gammaID", gammaID, &b_gammaID);
   fChain->SetBranchAddress("gamma", gamma, &b_gamma);
   fChain->SetBranchAddress("gamma_t", gamma_t, &b_gamma_t);
   fChain->SetBranchAddress("gaggID", gaggID, &b_gaggID);
   fChain->SetBranchAddress("gaggP", gaggP, &b_gaggP);
   fChain->SetBranchAddress("gaggT", gaggT, &b_gaggT);
   fChain->SetBranchAddress("gagg_t", gagg_t, &b_gagg_t);
   
   printf("======================== Start processing....\n");
   StpWatch.Start();
   
   Frac = 0.1; 
   ProcessedEntries = 0; 

}

Bool_t PIDAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PIDAnalyzer::SlaveBegin(TTree * /*tree*/){

   TString option = GetOption();

}


void PIDAnalyzer::SlaveTerminate(){

}

#endif // #ifdef PIDAnalyzer_cxx
