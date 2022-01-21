#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TBenchmark.h"
#include <vector>

#define MAX_CRATES 2
#define MAX_BOARDS_PER_CRATE 13
#define MAX_CHANNELS_PER_BOARD 16
#define BOARD_START 2

#include "../mapping.h"
#include "../armory/DataBlock.h"
#include "../armory/evtReader.h"

//#############################################
//           main 
//#############################################
int main(int argn, char **argv) {
  
  printf("=====================================\n");
  printf("===      evt --> _raw.root        ===\n");
  printf("=====================================\n");  
    
  if (argn < 3  )    {
    printf("Usage :\n");
    printf("%s [outFile] [evt1] [evt2] [evt3] ..... \n", argv[0]);
    printf("e.g.: \n");
    printf("%s hahaha_raw.root  haha-000.evt  haha-001.evt  haha-002.evt\n", argv[0]);
    printf("%s hahaha_raw.root  `ls haha-*.evt`\n", argv[0]);
    return 1;
  }
  
  TString outFileName = argv[1];
  int nFiles = argn-2;
  TString inFileName[nFiles];
  for( int i = 0; i < nFiles ; i++){
    inFileName[i] = argv[i+2];
    printf(" in file - %2d: %s\n", i,  inFileName[i].Data());
  }
  
  printf("     out file: %s\n", outFileName.Data());


  evtReader * evt = new evtReader();
  DataBlock * data = evt->data;
  
  printf("====================================\n");

  //====== ROOT file
  TFile * outFile = new TFile(outFileName, "recreate");
  TTree * tree = new TTree("tree", "tree");

  tree->Branch("evID",    &data->eventID, "data_ID/L"); 
  tree->Branch("id",           &data->id, "ID/s");
  tree->Branch("e",        &data->energy, "crystal_energy/s");
  tree->Branch("e_t",        &data->time, "crystal_timestamp/l");

  
  TBenchmark gClock;
  gClock.Reset();
  gClock.Start("timer");


  //=========================================
  //=========================================
  //=========================================
  //=========================================
  for( int i = 0; i < nFiles; i++){
    
    evt->OpenFile(inFileName[i]);
    if( evt->IsOpen() == false ) continue;
    
    Long64_t measureCount = 0;
    printf("\033[1;31mProcessing file: %s\033[0m\n", inFileName[i].Data());
    TBenchmark clock2;
    clock2.Reset();
    clock2.Start("timer");

    //=============== Read File
    while( evt->IsEndOfFile() == false ){
      
      evt->ReadBlock();
      evt->PrintStatus(10000);
      
      //cern fill tree
      outFile->cd();
      tree->Fill();

    }
    
    clock2.Stop("timer");
    double time = clock2.GetRealTime("timer");
    float tempf = (float)evt->GetFilePos()/(1024.*1024.*1024.);
    printf("           measurements: \x1B[32m%lld \x1B[0m |  %.3f GB\n", evt->GetBlockID(), tempf);
    printf("           Time used:%3.0f min %5.2f sec\n",  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);
    printf("           Root file size so far: %.4f GB\n",  outFile->GetSize()/1024./1024./1024.);
    
  }
  
  gClock.Stop("timer");
  double time = gClock.GetRealTime("timer");
  gClock.Start("timer");
  float tempf = (float)evt->GetFilePos()/(1024.*1024.*1024.);
  printf("Total measurements: \x1B[32m%lld \x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\nTime used:%3.0f min %5.2f sec\033[A\r", 
                   evt->GetBlockID()+1, (100*evt->GetFilePos()/evt->GetFileSize()), tempf,  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);


  //cern save root
  outFile->cd();
  double totRootSize = outFile->GetSize()/1024./1024./1024.;
  tree->Write();
  outFile->Close();

  gClock.Stop("timer");
  time = gClock.GetRealTime("timer");
  printf("\n==================== finished.\r\n");
  printf("Total time spend : %3.0f min %5.2f sec\n", TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);
  printf(" File size of %s : %.3f GB \n", outFileName.Data(), totRootSize);

}
