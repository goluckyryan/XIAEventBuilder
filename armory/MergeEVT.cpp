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

class measurment{

public:
  UShort_t               ch;
  UShort_t             slot;
  UShort_t            crate;
  UShort_t     headerLength;  /// headerLength > 4, more data except tarce.
  UShort_t      eventLength;  /// eventLength = headerLength + trace 
  Bool_t             pileup;
  ULong64_t            time;
  UShort_t              cfd;
  UShort_t           energy;
  UShort_t     trace_length;
  Bool_t trace_out_of_range;
  
  Long64_t   timeDiff;  
  
  UShort_t id;
  
  measurment(){};
  
  void Clear(){
    ch = 0;
    slot = 0;
    crate = 0;
    eventLength = 0;
    pileup = false;
    time = 0;
    cfd = 0;
    energy = 0;
    trace_length = 0;
    trace_out_of_range = 0;
    timeDiff = 0;
    id = 0;
  }
  
  void Print(){
    printf("Crate: %d, Slot: %d, Ch: %d | id: %d \n", crate, slot, ch, id);
    printf("HeaderLength: %d, Event Length: %d, energy: %d, timeStamp: %llu\n", headerLength, eventLength, energy, time);
    printf("trace_length: %d, pile-up:%d\n", trace_length, pileup); 
  }
};
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

  Long64_t measureID = -1;  
  measurment data;
  
  printf("====================================\n");

  //====== ROOT file
  TFile * outFile = new TFile(outFileName, "recreate");
  TTree * tree = new TTree("tree", "tree");

  tree->Branch("evID",    &measureID, "data_ID/L"); 
  tree->Branch("id",        &data.id, "ID/s");
  tree->Branch("e",     &data.energy, "crystal_energy/s");
  tree->Branch("e_t",       &data.time, "crystal_timestamp/l");


  long int inFileSize = 0;
  long int inFilePos = 0;
  
  TBenchmark gClock;
  gClock.Reset();
  gClock.Start("timer");

  //=========================================
  //=========================================
  //=========================================
  //=========================================
  for( int i = 0; i < nFiles; i++){
    FILE * inFile = fopen(inFileName[i], "r");
    if( inFile == NULL ){
      printf("Cannot read file : %s \n", inFileName[i].Data());
      continue;
    }
    
    Long64_t measureCount = 0;
    printf("\033[1;31mProcessing file: %s\033[0m\n", inFileName[i].Data());
    TBenchmark clock2;
    clock2.Reset();
    clock2.Start("timer");

    //get file size
    fseek(inFile, 0L, SEEK_END);
    inFileSize = ftell(inFile);
    rewind(inFile); ///back to the File begining
    inFilePos = 0;
    
    unsigned int header[4]; //read 4 header, unsigned int = 4 byte = 32 bits.  
    unsigned long long nWords = 0;
    
    //ULong64_t timeLast = 0;

    //=============== Read File
    ///  while ( ! feof(inFile) ){
    while( inFilePos < inFileSize || feof(inFile) ){

      fread(header, sizeof(header), 1, inFile);
      inFilePos = ftell(inFile);
      measureID ++; 
      measureCount++;
      
      /// see the Pixie-16 user manual, Table4-2
      data.ch           =  header[0] & 0xF ;
      data.slot         = (header[0] >> 4) & 0xF;
      data.crate        = (header[0] >> 8) & 0xF;
      data.headerLength = (header[0] >> 12) & 0x1F;
      data.eventLength  = (header[0] >> 17) & 0x3FFF;
      data.pileup       =  header[0] >> 31 ;
      data.time         = ((ULong64_t)(header[2] & 0xFFFF) << 32) + header[1];
      data.cfd          =  header[2] >> 16 ; 
      data.energy       =  header[3] & 0xFFFF;
      data.trace_length = (header[3] >> 16) & 0x7FFF;
      data.trace_out_of_range =  header[3] >> 31;
      
      data.id = data.crate*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD + (data.slot-BOARD_START)*MAX_CHANNELS_PER_BOARD + data.ch;

      data.id = mapping[data.id];
      data.energy = data.energy / 2. ; // factor 2 is the rawe_rebin_factor;

      nWords += data.eventLength;
      
      //=== jump to next measurement
      if( data.eventLength > 4 ){
        fseek(inFile, sizeof(int) * (data.eventLength-4),  SEEK_CUR);
        inFilePos = ftell(inFile);
      }
    
      //event stats, print status every 10000 events
      if ( measureID % 10000 == 0 ) {
        float tempf = (float)inFileSize/(1024.*1024.*1024.);
        gClock.Stop("timer");
        double time = gClock.GetRealTime("timer");
        gClock.Start("timer");
        printf("Total measurements: \x1B[32m%lld \x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\nTime used:%3.0f min %5.2f sec\033[A\033[A\r", 
                     measureID +1 , (100*inFilePos/inFileSize), tempf,  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);
      }   
      
      //cern fill tree
      outFile->cd();
      tree->Fill();

    }
    inFilePos = ftell(inFile);
    fclose(inFile);
    
    clock2.Stop("timer");
    double time = clock2.GetRealTime("timer");
    float tempf = (float)inFileSize/(1024.*1024.*1024.);
    printf("           measurements: \x1B[32m%lld \x1B[0m |  %.3f GB\n", measureCount, tempf);
    printf("           Time used:%3.0f min %5.2f sec\n",  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);
    printf("           Root file size so far: %.4f GB\n",  outFile->GetSize()/1024./1024./1024.);
    
  }
  
  gClock.Stop("timer");
  double time = gClock.GetRealTime("timer");
  gClock.Start("timer");
  float tempf = (float)inFileSize/(1024.*1024.*1024.);
  printf("Total measurements: \x1B[32m%lld \x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\nTime used:%3.0f min %5.2f sec\033[A\r", 
                   measureID+1, (100*inFilePos/inFileSize), tempf,  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);


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
