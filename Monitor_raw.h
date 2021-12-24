#ifndef Monitor_raw_h
#define Monitor_raw_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include "mapping.h"
#include "armory/AnalysisLibrary.h"

// Header file for the classes stored in the TTree if any.

class Monitor_raw : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Long64_t        evID;
   UShort_t        ID;
   UShort_t        e;
   ULong64_t       t;

   // List of branches
   TBranch        *b_data_ID;   //!
   TBranch        *b_ID;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_timestamp;   //!

   Monitor_raw(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Monitor_raw() { }
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

   ClassDef(Monitor_raw,0);
};

#endif

#ifdef Monitor_raw_cxx
void Monitor_raw::Init(TTree *tree)
{
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evID", &evID, &b_data_ID);
   fChain->SetBranchAddress("id",     &ID, &b_ID);
   fChain->SetBranchAddress("e",       &e, &b_energy);
   fChain->SetBranchAddress("t",       &t, &b_timestamp);
}

Bool_t Monitor_raw::Notify()
{
   return kTRUE;
}


void Monitor_raw::SlaveBegin(TTree * /*tree*/){
   TString option = GetOption();
}


void Monitor_raw::SlaveTerminate(){}


#endif // #ifdef Monitor_raw_cxx
