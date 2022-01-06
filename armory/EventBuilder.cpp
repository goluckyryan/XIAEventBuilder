#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TBenchmark.h"
#include "TStopwatch.h"
#include "TTreeIndex.h"

#include "../mapping.h"

Int_t eventID = 0 ;
double       e[NCRYSTAL];
ULong64_t   e_t[NCRYSTAL];
double    bgo[NBGO];
ULong64_t bgo_t[NBGO];
Short_t   other[NOTHER];
Short_t  multi;

void ClearTreeData(){
   
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

int main(int argn, char **argv){
  printf("=====================================\n");
  printf("===          Event Builder        ===\n");
  printf("=====================================\n");  
    
  if (argn != 2 && argn != 3  && argn != 4 )    {
    printf("Usage :\n");
    printf("%s [_raw.root File] <timeWindows> <SaveFileName>\n", argv[0]);
    printf("         timeWindows : default = 100 \n");
    printf("        SaveFileName : default is *.root \n");
    return 1;
  }
    
  TString inFileName = argv[1]; // need to check name 
  int timeWindow = 100;
  if( argn >= 3 ) timeWindow = atoi(argv[2]);
  
  printf("======================== Opening input _raw.root \n");
  TFile * inFile = new TFile(inFileName, "READ");
  if( inFile->IsOpen() == false ) {
    printf("!!!! cannot open file %s \n", inFileName.Data());
    return 0;
  }
  
  TTree * tree = (TTree *) inFile->Get("tree");
  
  Long64_t        evID;
  UShort_t        detID;
  UShort_t        energy;
  ULong64_t       energy_t;

  TBranch        *b_data_ID;   //!
  TBranch        *b_ID;   //!
  TBranch        *b_energy;   //!
  TBranch        *b_energy_timestamp;   //!

  tree->SetBranchAddress("evID",    &evID, &b_data_ID);
  tree->SetBranchAddress("id",     &detID, &b_ID);
  tree->SetBranchAddress("e",     &energy, &b_energy);
  tree->SetBranchAddress("e_t", &energy_t, &b_energy_timestamp);
  
  Long64_t totnumEntry = tree->GetEntries();
  
  printf( "total Entry : %lld \n", totnumEntry);
  
  printf("======================== Buidling Index using the timestamp\n");
  tree->BuildIndex("e_t");
  TTreeIndex *in = (TTreeIndex*) tree->GetTreeIndex(); 
  Long64_t * index = in->GetIndex();
  
  ULong64_t time0;  //time-0 for each event
  int       timeDiff; 

  printf("======================== Create output tree\n");
  TString outFileName = inFileName;
  outFileName.Remove(inFileName.First('.'));
  outFileName.Append(".root");
  if( argn >=4 ) outFileName = argv[3];
  
  TFile * saveFile = new TFile(outFileName, "recreate");
  saveFile->cd();
  TTree * newtree = new TTree("tree", "tree");
  
  newtree->Branch("evID", &eventID, "event_ID/l"); 
  newtree->Branch("e",         e, Form("e[%d]/D", NCRYSTAL));
  newtree->Branch("e_t",     e_t, Form("e_timestamp[%d]/l", NCRYSTAL));
  //newtree->Branch("p",    pileup, Form("pile_up_flag[%d]/s", NCRYSTAL));
  //newtree->Branch("hit",     hit, Form("hit[%d]/s", NCRYSTAL));

  newtree->Branch("bgo",     bgo, Form("BGO_e[%d]/D", NBGO));
  newtree->Branch("bgo_t", bgo_t, Form("BGO_timestamp[%d]/l", NBGO));

  newtree->Branch("other", other, Form("other_e[%d]/D", NOTHER));

  newtree->Branch("multi", &multi, "multiplicity_crystal/I");
  
  ClearTreeData();
  
  printf("======================== Start processing....\n");
  Float_t Frac = 0.1; ///Progress bar
  TStopwatch StpWatch;
  StpWatch.Start();
  eventID = 0;
  
  for( Long64_t entry = 0; entry < totnumEntry; entry++){
    

    /*********** Progress Bar ******************************************/ 
    if (entry>totnumEntry*Frac-1) {
      TString msg; msg.Form("%llu", totnumEntry/1000);
      int len = msg.Sizeof();
      printf(" %3.0f%% (%*llu/%llu k) processed in %6.1f sec | expect %6.1f sec\n",
             Frac*100, len, entry/1000,totnumEntry/1000,StpWatch.RealTime(), StpWatch.RealTime()/Frac);
      StpWatch.Start(kFALSE);
      Frac+=0.1;
    }
    
    entry = index[entry];
    
    b_ID->GetEntry(entry);
    b_energy->GetEntry(entry);
    b_energy_timestamp->GetEntry(entry);
    
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

  printf("============================== finished.\n");

  saveFile->cd();
  newtree->Write();
  saveFile->Close();
  
  printf(" total number of event Built : %d \n", eventID);
  
}
