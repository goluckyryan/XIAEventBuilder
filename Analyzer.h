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
#include <TStopwatch.h>

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
   ULong64_t       e_t[NCRYSTAL];
   UShort_t        p[NCRYSTAL];
   UShort_t        hit[NCRYSTAL];
   Double_t        bgo[NBGO];
   ULong64_t       bgo_t[NBGO];
   Double_t        other[NOTHER];
   Int_t           multi;

   // List of branches
   TBranch        *b_event_ID;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_time;   //!
   TBranch        *b_pileup;   //!
   TBranch        *b_hit;   //!
   TBranch        *b_bgo;   //!
   TBranch        *b_bgoTime;   //!
   TBranch        *b_other;   //!
   TBranch        *b_multi;   //!

   Analyzer(TTree * /*tree*/ =0) : fChain(0) { totnumEntry = 0; 
                                               Frac = 0.1; 
                                               ProcessedEntries = 0; 
                                               saveEV2 = true; 
                                               outEV2Name = "test.ev2";
                                             }
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

   ULong64_t totnumEntry;

   vector<vector<double>> eCorr;

   ULong64_t ProcessedEntries;
   Float_t Frac; ///Progress bar
   TStopwatch StpWatch;

   bool saveEV2;
   FILE * outEV2;
   TString outEV2Name;



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
   fChain->SetBranchAddress("e_t",      e_t, &b_time);
   fChain->SetBranchAddress("p",          p, &b_pileup);
   fChain->SetBranchAddress("hit",      hit, &b_hit);
   fChain->SetBranchAddress("bgo",      bgo, &b_bgo);
   fChain->SetBranchAddress("bgo_t",  bgo_t, &b_bgoTime);
   fChain->SetBranchAddress("other",  other, &b_other);
   fChain->SetBranchAddress("multi", &multi, &b_multi);

   if( saveEV2 ){
      printf("======================== Create output ev2 : %s\n", outEV2Name.Data());
      outEV2 = fopen(outEV2Name.Data(), "w+");
   }  


   printf("======================== Start processing....\n");
   StpWatch.Start();

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
