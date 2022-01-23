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

// Header file for the classes stored in the TTree if any.'

#define MAX_MULTI 200

class Analyzer : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       evID;
   Int_t           runID;
   Int_t           detID[MAX_MULTI];
   Double_t        e[MAX_MULTI];
   ULong64_t       e_t[MAX_MULTI];
   Int_t           qdc[MAX_MULTI][8];
   Int_t           multi;
   Int_t           multiCry;
   Int_t           multiGagg;

   // List of branches
   TBranch        *b_event_ID;   //!
   TBranch        *b_runID;   //!
   TBranch        *b_detID;     //!
   TBranch        *b_energy;   //!
   TBranch        *b_time;   //!
   TBranch        *b_qdc;   //!
   TBranch        *b_multi;   //!
   TBranch        *b_multiCry;   //!
   TBranch        *b_multiGagg;   //!

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

   void Save2ev2();
   bool saveEV2;
   FILE * outEV2;
   TString outEV2Name;
   
   
   double eCal[NCRYSTAL];
   double gamma[NCLOVER]; // added-back energy
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

   fChain->SetBranchAddress("evID",           &evID, &b_event_ID);
   fChain->SetBranchAddress("runID",         &runID, &b_runID);
   fChain->SetBranchAddress("detID",          detID, &b_detID);
   fChain->SetBranchAddress("e",                  e, &b_energy);
   fChain->SetBranchAddress("e_t",              e_t, &b_time);
   fChain->SetBranchAddress("qdc",              qdc, &b_qdc);
   fChain->SetBranchAddress("multi",         &multi, &b_multi);
   fChain->SetBranchAddress("multiCry",   &multiCry, &b_multiCry);
   fChain->SetBranchAddress("multiGagg", &multiGagg, &b_multiGagg);

   TString option = GetOption();
   if ( option != "" ) outEV2Name = option;

   if( saveEV2 ){
      printf("======================== Create output ev2 : \033[1;31m%s\033[0m\n", outEV2Name.Data());
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

void Analyzer::Save2ev2(){
   short *  out0 = new short[1];
   short *  outa = new short[1];
   short * outb = new short[1];

   int count = 0;
   for (int i = 0; i < NCRYSTAL ; i++){
      if( !TMath::IsNaN(eCal[i]) ) count++;
   }

   out0[0] = count;
   if( count == 0 ) return;
   
   fwrite(out0, 1, 1, outEV2);
   for( int i = 0; i < count; i++){
      if( TMath::IsNaN(eCal[i]) ) continue;
      outa[0] = i;
      fwrite(outa, 1, 1, outEV2); 
      outb[0] = TMath::Nint(eCal[i]);
      fwrite(outb, 2, 1, outEV2); 
   }

   fwrite(out0, 1, 1, outEV2); 
   
      
   /**
   int len = (int) gatedID.size();
   char out[2*len+2];
   out[0] = numGatedData;
   for( int i = 0; i < (int) gatedID.size(); i++){
      int id = gatedID[i];
      out[2*i+1] = id;
      out[2*i+2] = eCal[id];
   }
   out[2*len+1]=numGatedData;
   fwrite(out,3*len+2,sizeof(out),outEV2);
   */ 
}

#endif // #ifdef Analyzer_cxx
