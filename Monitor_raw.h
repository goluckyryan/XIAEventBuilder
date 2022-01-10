#ifndef Monitor_raw_h
#define Monitor_raw_h

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TMacro.h>
#include <TSelector.h>
#include <TStopwatch.h>
#include <TTreeIndex.h>

#include "mapping.h"
#include "armory/AnalysisLibrary.h"

// Header file for the classes stored in the TTree if any.

class Monitor_raw : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Long64_t        evID;
   UShort_t        detID;
   UShort_t        energy;
   ULong64_t       energy_t;

   // List of branches
   TBranch        *b_data_ID;   //!
   TBranch        *b_ID;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_energy_timestamp;   //!

   Monitor_raw(TTree * /*tree*/ =0) : fChain(0) { isBuildEvent = false; 
                                                  timeWindow = 100;
                                                  saveFileName = "test.root";
                                                }
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
   
   void SetBuildEvent(bool buildOption = false){ isBuildEvent = buildOption;}

   ClassDef(Monitor_raw,0);
   
   //=================== progress
   ULong64_t ProcessedEntries = 0;
   Float_t Frac = 0.1; ///Progress bar
   TStopwatch StpWatch;
   
   //=================== correction parameters
   vector<vector<double>> eCorr;
   
   //=================== output tree;   
   bool isBuildEvent;
   int timeWindow;
   ULong64_t time0;  //time-0 for each event
   int       timeDiff; 
   Long64_t * index;  //!

   TFile * saveFile; //!
   TTree * newtree; //!
   TString saveFileName;
   Long64_t totnumEntry; /// of original root
   
   //tree  
   void ClearTreeData();
   void BuildEvent();
   
   Int_t eventID;
   double       e[NCRYSTAL];
   ULong64_t   e_t[NCRYSTAL];
   double    bgo[NBGO];
   ULong64_t bgo_t[NBGO];
   Short_t   other[NOTHER];
   Short_t  multi;
   
   
};

#endif

#ifdef Monitor_raw_cxx
void Monitor_raw::Init(TTree *tree)
{
   // Set branch addresses and branch pointers
   if (!tree) return;
   
   fChain = (TChain *) tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evID",    &evID, &b_data_ID);
   fChain->SetBranchAddress("id",     &detID, &b_ID);
   fChain->SetBranchAddress("e",     &energy, &b_energy);
   fChain->SetBranchAddress("e_t", &energy_t, &b_energy_timestamp);
   
   if( isBuildEventRoot || isBuildEvent){
   
      printf("======================== Buidling Index using the timestamp\n");
      tree->BuildIndex("e_t");
      TTreeIndex *in = (TTreeIndex*) tree->GetTreeIndex(); 
      index = in->GetIndex();     
      
      //for(int i = 0; i < 100; i++){
      //   printf(" %3d | %lld \n", i, index[i]);
      //} 
      
      
      TString option = GetOption();
      
      //printf("======================== Formation of output file name\n");
      //int numFile = fChain->GetListOfFiles()->GetLast() + 1;   //need input of TChain
      //
      //printf(".......... number of files : %d \n", numFile);
      //
      //if( numFile > 0 ) {
      //   int oldRunNum = -100;
      //   bool contFlag = false; // is runNumber continue;
      //   for( int i = 0; i < numFile ; i++){
      //      TString name = ((TChain *) fChain)->GetListOfFiles()->At(i)->GetTitle();
      //      int found = name.Last('/');
      //      name.Remove(0, found + 1 ); // this should give "XXX_run0XX.root"
      //      TString prefix = name;
      //      found = name.Last('.');
      //      name.Remove(found); // this should give "XXX_run0XX"
      //      found = name.Last('_');
      //      int runNum = name.Remove(0, found+4).Atoi(); // this should give the 3 digit run number 
      //
      //      if( i == 0 ) {
      //         found = prefix.First("_");
      //         prefix.Remove(found);
      //         saveFileName = prefix +  "_run";
      //      }
      //
      //      if( runNum == oldRunNum + 1 ){
      //         int kk = saveFileName.Sizeof();
      //         if( contFlag == false ){
      //            saveFileName.Remove(kk-2); //remove the "-"
      //            saveFileName += "-";
      //         }else{
      //            saveFileName.Remove(kk-5); //remove the runNum and "-"
      //         }
      //         contFlag = true;
      //      }
      //      if( runNum > oldRunNum + 1) contFlag = false;
      //      
      //      saveFileName += Form("%03d_", runNum);
      //      oldRunNum = runNum;
      //   }
      //   int kk = saveFileName.Sizeof();
      //   saveFileName.Remove(kk-2); // remove the last "-"
      //   saveFileName += ".root";
      //}else{
      //   gROOT->ProcessLine(".q");
      //}
      
      if( option != "" ) saveFileName = option;
      
      printf("save file name : %s \n", saveFileName.Data());
      printf("---------------------------------------------\n");
            
      //TODO;
      saveFile = new TFile(saveFileName, "recreate");
      saveFile->cd();
      newtree = new TTree("tree", "tree");

      printf("======================== Create output tree\n");
      
      newtree->Branch("evID", &eventID, "event_ID/l"); 
      newtree->Branch("e",         e, Form("e[%d]/D", NCRYSTAL));
      newtree->Branch("e_t",     e_t, Form("e_timestamp[%d]/l", NCRYSTAL));
      //newtree->Branch("p",    pileup, Form("pile_up_flag[%d]/s", NCRYSTAL));
      //newtree->Branch("hit",     hit, Form("hit[%d]/s", NCRYSTAL));

      newtree->Branch("bgo",     bgo, Form("BGO_e[%d]/D", NBGO));
      newtree->Branch("bgo_t", bgo_t, Form("BGO_timestamp[%d]/l", NBGO));

      newtree->Branch("other", other, Form("other_e[%d]/D", NOTHER));

      newtree->Branch("multi", &multi, "multiplicity_crystal/I");
      
      
      if( eCorr.size() > 0 ) {
         TMacro energyCorr(e_corr);
         energyCorr.Write("energyCorr");
      }
      
   }
   
   
   printf("======================== Start processing....\n");
   StpWatch.Start();
   eventID = 0;
      
}


void Monitor_raw::ClearTreeData(){
   
   for( int i = 0; i < NCRYSTAL; i++){
      e[i]      = TMath::QuietNaN();
      e_t[i]    = 0;
      //pileup[i] = 0;
      //hit[i]    = 0;
   }
   for( int i = 0; i < NBGO; i++) {
      bgo[i]   = TMath::QuietNaN();
      bgo_t[i] = 0 ;
   }
   for( int i = 0; i < NOTHER; i++) {
      other[i] = TMath::QuietNaN(); 
   }
   multi = 0;
}

void Monitor_raw::BuildEvent(){

   if( time0 == 0) time0 = energy_t;
   
   timeDiff = (int) (energy_t - time0);
   
   if( timeDiff < timeWindow ) {
      
      if ( detID < NCRYSTAL ){
         e[detID] = energy;
         e_t[detID] = energy_t;
         multi++;
      }
      if ( 100 <= detID && detID < 100 + NBGO ){
         bgo[detID-100]   = energy;
      bgo_t[detID-100] = energy_t;
      }
      if ( 200 <= detID && detID < 200 + NOTHER){
         other[detID-200] = energy;
      }
      
      //printf("%d | %3d   %6d  %10llu, %3d\n", multi, detID, energy, energy_t, timeDiff);
      
   }else{
      //---- end of event
      eventID ++;

      saveFile->cd();
      newtree->Fill();
      
      ClearTreeData();
      
      /// fill 1st data of an event                  
      time0 = energy_t;
      timeDiff = 0;
      
      if ( detID < NCRYSTAL ){
         e[detID] = energy;
         e_t[detID] = energy_t;
         multi = 1;
      }
      if ( 100 <= detID && detID < 100 + NBGO ){
         bgo[detID-100]   = energy;
      bgo_t[detID-100] = energy_t;
      }
      if ( 200 <= detID && detID < 200 + NOTHER){
         other[detID-200] = energy;
      }
      
   }
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
