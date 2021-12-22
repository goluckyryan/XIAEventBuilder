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
    
  if (argn != 2 && argn != 3 )    {
    printf("Usage :\n");
    printf("%s [evt File] [timeWindow] \n", argv[0]);
    printf("         timeWindow : number of tick, 1 tick = 10 ns. default = 100 \n");       
    return 1;
  }
  
  TString inFileName = argv[1];
  TString outFileName =  inFileName;
  outFileName.Remove(inFileName.First('.'));
  outFileName.Append("_raw.root");
  
  long int inFilePos = 0;
  TBenchmark gClock;
  gClock.Reset();
  gClock.Start("timer");
  
  Long64_t measureID = -1;
  
  measurment data;
  
  printf("====================================\n");

  FILE * inFile = fopen(inFileName, "r");
  if( inFile == NULL ){
    printf("Cannot read file : %s \n", inFileName.Data());
    return -404;
  }

  //get file size
  fseek(inFile, 0L, SEEK_END);
  long int inFileSize = ftell(inFile);
  rewind(inFile); ///back to the File begining

  printf(" in file: %s\n", inFileName.Data());
  printf("out file: %s\n", outFileName.Data());
  printf("--------------------------------\n");

  //====== ROOT file
  TFile * outFile = new TFile(outFileName, "recreate");
  TTree * tree = new TTree("tree", "tree");
  
  tree->Branch("evID",    &measureID, "data_ID/L"); 
  tree->Branch("detID",     &data.id, "det_ID/s");
  tree->Branch("e",     &data.energy, "energy/s");
  tree->Branch("t",       &data.time, "time_stamp/l");
  //tree->Branch("tdiff", &data.timeDiff, "time_Diff/L");
  
//=======TODO online event building

  unsigned int header[4]; //read 4 header, unsigned int = 4 byte = 32 bits.  
  unsigned long long nWords = 0;
  
  ULong64_t timeLast = 0;

  //=============== Read File
  ///  while ( ! feof(inFile) ){
  while ( inFilePos <= inFileSize ){

    fread(header, sizeof(header), 1, inFile);
    inFilePos = ftell(inFile);
    measureID ++; 
    
    /// see the Pixie-16 user manual, Table4-2
    data.ch =  header[0] & 0xF ;
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
    
    nWords += data.eventLength;
    
    //if( measureID == 0 ) {
    //  data.timeDiff = 0;
    //}else{
    //  data.timeDiff = (Long64_t) data.time - timeLast;
    //}
    //timeLast = data.time;
    
    //if( data.timeDiff == false ){
    //  printf("----------------------nWords: %llu, inFilePos: %lu\n", nWords, inFilePos);
    //  for(int i = 0; i < 4; i++){
    //    printf("  %x\n", header[i]);
    //  }
    //  data.Print();
    //}
    
    //=== jump to next measurement
    if( data.eventLength > 4 ){
      fseek(inFile, sizeof(int) * (data.eventLength-4),  SEEK_CUR);
      inFilePos = ftell(inFile);
    }
  
    //event stats, print status every 10000 events
    if ( measureID % 10000 == 0 ) {
      inFilePos = ftell(inFile);
      float tempf = (float)inFileSize/(1024.*1024.*1024.);
      gClock.Stop("timer");
      double time = gClock.GetRealTime("timer");
      gClock.Start("timer");
      printf("Total measurements: \x1B[32m%lld \x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\nTime used:%3.0f min %5.2f sec\033[A\033[A\r", 
                   measureID, (100*inFilePos/inFileSize), tempf,  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);
    }   
    
    //cern fill tree
    outFile->cd();
    tree->Fill();

  }
  
  inFilePos = ftell(inFile);
  gClock.Stop("timer");
  double time = gClock.GetRealTime("timer");
  gClock.Start("timer");
  float tempf = (float)inFileSize/(1024.*1024.*1024.);
  printf("Total measurements: \x1B[32m%lld \x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\nTime used:%3.0f min %5.2f sec\033[A\r", 
                   measureID, (100*inFilePos/inFileSize), tempf,  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);

  fclose(inFile);

  //cern save root
  outFile->cd();
  tree->Write();
  outFile->Close();

  gClock.Stop("timer");
  time = gClock.GetRealTime("timer");
  printf("\n==================== finished.\r\n");
  printf("Total time spend : %3.0f min %5.2f sec\n", TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);

}
