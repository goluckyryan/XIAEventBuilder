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

#define MAXMULTI 200

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
  
  printf(">>> Opening input %s \n", inFileName.Data());
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
  
  printf(" total Entry : %lld \n", totnumEntry);

  printf(" event Build window: %d tick = %d nsec \n", timeWindow, timeWindow * 10);
  
  printf(">>> Buidling Index using the timestamp\n");
  tree->BuildIndex("e_t");
  TTreeIndex *in = (TTreeIndex*) tree->GetTreeIndex(); 
  Long64_t * index = in->GetIndex();
  
  ULong64_t time0;  //time-0 for each event
  int       timeDiff; 
  
  TString outFileName = inFileName;
  outFileName.Remove(inFileName.First("_raw"));
  outFileName.Append(".root");
  if( argn >=4 ) outFileName = argv[3];
  
  printf(">>> out File name : \033[1,31m%s\033[m\n", outFileName.Data());
  
  printf(">>> Create output tree\n");
  TFile * saveFile = new TFile(outFileName, "recreate");
  saveFile->cd();
  TTree * newtree = new TTree("tree", "tree");
  
  Int_t multi = 0; /// this is total multipicilty for all detectors
  newtree->Branch("multi", &multi, "multi/I");  

  Int_t eventID = 0 ;
  newtree->Branch("evID", &eventID, "event_ID/l"); 
  
  Int_t multiCry = 0 ; /// thi is total multiplicity for all crystal
  newtree->Branch("multiCry", &multiCry, "multiplicity_crystal/I");  
  
  int id[MAXMULTI];
  double e[MAXMULTI];    
  ULong64_t e_t[MAXMULTI];
  newtree->Branch("id",       id, "id[multi]/I" );
  newtree->Branch("e",         e, "e[multi]/D" );
  newtree->Branch("e_t",     e_t, "e_timestamp[multi]/l");
  
  printf("================== Start processing....\n");
  Float_t Frac = 0.1; ///Progress bar
  TStopwatch StpWatch;
  StpWatch.Start();
  
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
    
    b_ID->GetEntry(entry, 0);
    b_energy->GetEntry(entry, 0);
    b_energy_timestamp->GetEntry(entry, 0);
    
    if( time0 == 0) {
      time0 = energy_t;
      multi = 0;
    }
    timeDiff = (int) (energy_t - time0);
  
    if( timeDiff < timeWindow ) {
      
      id[multi] = detID;
      e[multi] = energy;
      e_t[multi] = energy_t;
      multi ++;
      if( detID < NCRYSTAL ) multiCry++;
      
    }else{
      ///---- end of event
      saveFile->cd();
      newtree->Fill();
      eventID ++;
      
      ///---- clear data
      for( int i = 0; i < multi ; i ++){
        id[i] = 0;
        e[i] = TMath::QuietNaN();
        e_t[i] = 0;
      }
      multi = 0;
      multiCry = 0;
      
      
      /// fill 1st data of an event                  
      time0 = energy_t;
      timeDiff = 0;

      id[multi] = detID;
      e[multi] = energy;
      e_t[multi] = energy_t;
      multi++;
      if( detID < NCRYSTAL ) multiCry++;

    }
    
  }

  printf("============================== finished.\n");

  saveFile->cd();
  newtree->Write();
  saveFile->Close();
  
  printf(" total number of event Built : %d \n", eventID);
  
}
