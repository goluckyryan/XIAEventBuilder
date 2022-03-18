//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec 13 18:37:46 2021 by ROOT version 6.24/06
// from TTree tree/tree
// found on file: efEu152b-000.root
//////////////////////////////////////////////////////////

#ifndef PreAnalyzer_h
#define PreAnalyzer_h

#include <fstream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TStopwatch.h>
#include <TMacro.h>
#include <vector>

#include "mapping.h"
#include "armory/AnalysisLibrary.h"

// Header file for the classes stored in the TTree if any.'

#define MAX_MULTI 200

class PreAnalyzer : public TSelector {
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
   Bool_t          pileup[MAX_MULTI];

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
   TBranch        *b_pileup;    //!

   PreAnalyzer(TTree * /*tree*/ =0) : fChain(0) { totnumEntry = 0; 
                                               Frac = 0.1; 
                                               ProcessedEntries = 0; 
                                             }
   virtual ~PreAnalyzer() { }
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

   ClassDef(PreAnalyzer,0);

   ULong64_t totnumEntry;

   vector<vector<double>> eCorr;

   ULong64_t ProcessedEntries;
   Float_t Frac; ///Progress bar
   TStopwatch StpWatch;
   
   double eCal[NCRYSTAL];
   double    gamma[NCLOVER]; // added-back energy
   
   ///========================= for new root file
   TFile * saveFile;
   TTree * newTree;
   TString saveFileName;
   
   ///tree
   ULong_t eventID;
   Int_t runID_N;
   int   multi_N;
   int   multiGagg_N;
   int       gammaID[NCLOVER];
   double    gamma_N[NCLOVER]; 
   ULong64_t gamma_t[NCLOVER];
   int       gaggID[NGAGG];
   double    gagg_peak[NGAGG];
   double    gagg_tail[NGAGG];
   ULong64_t gagg_t[NGAGG];

};

#endif

#ifdef PreAnalyzer_cxx
void PreAnalyzer::Init(TTree *tree)
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
   fChain->SetBranchAddress("pileup",        pileup, &b_pileup);

   TString option = GetOption();
      
   printf("======================== Load parameters.\n");
   eCorr = LoadCorrectionParameters("correction_e.dat"); 
   
   ///======================== open a new file
   saveFileName = "haha.root";
   
   saveFile = new TFile( saveFileName,"recreate");
   
   TMacro e_corr("correction_e.dat");
   e_corr.Write("correction_e");
   
   newTree = new TTree("tree", "tree");
   
   eventID = -1;
   runID = 0;
   
   multi_N = 0;
   multiGagg_N = 0;
   
   newTree->Branch("eventID",       &eventID, "eventID/l");
   newTree->Branch("runID",         &runID_N, "runID/I");
   newTree->Branch("multi",         &multi_N, "multi/I"); 
   newTree->Branch("multiGagg", &multiGagg_N, "multiGagg/I"); 
   newTree->Branch("gammaID",        gammaID, "gammaID[multi]/S");
   newTree->Branch("gamma",          gamma_N, "gamma[multi]/D");
   newTree->Branch("gamma_t",        gamma_t, "gamma_t[multi]/l");
   newTree->Branch("gaggID",          gaggID, "gaggID[multiGagg]/I");
   newTree->Branch("gaggP",        gagg_peak, "gaggP[multiGagg]/D");
   newTree->Branch("gaggT",        gagg_tail, "gaggT[multiGagg]/D");
   newTree->Branch("gagg_t",          gagg_t, "gagg_t[multiGagg]/l");
   
   printf("======================== Start processing....\n");
   StpWatch.Start();

}

Bool_t PreAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}



void PreAnalyzer::SlaveBegin(TTree * /*tree*/){

   TString option = GetOption();

}


void PreAnalyzer::SlaveTerminate(){

}

#endif // #ifdef PreAnalyzer_cxx
