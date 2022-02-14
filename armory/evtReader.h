#ifndef EVTREADER_H
#define EVTREADER_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>

#include "TSystem.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TBenchmark.h"

#include "../mapping.h"
#include "../armory/DataBlock.h"

#define MAX_CRATES 2
#define MAX_BOARDS_PER_CRATE 13
#define MAX_CHANNELS_PER_BOARD 16
#define BOARD_START 2

class evtReader{

  public:
    DataBlock * data;

  private:
    FILE * inFile;
    
    long int  inFileSize;
    long int  inFilePos;
    bool      endOfFile;  
    bool      isOpened;
    Long64_t  blockID;
    long int  nBlock;
    
    unsigned int extraHeader[14];
    unsigned int traceBlock[4000];
    
    TBenchmark gClock;

  ///============================================ Methods
  public:

    evtReader();
    evtReader(TString inFileName);
    ~evtReader();
    
    void OpenFile(TString inFileName);
    
    void UpdateFileSize();
    bool IsEndOfFile();
    
    bool IsOpen()               {return isOpened;}         
    long int GetFilePos()       {return inFilePos;}
    long int GetFileSize()      {return inFileSize;}
    Long64_t GetBlockID()       {return blockID;}
    Long64_t GetNumberOfBlock() {return nBlock;}
    
    int ReadBlock(int opt = 0);  /// 0 = default, fill data
                                 /// 1 = no fill data

    void ScanNumberOfBlock();
    void PrintStatus(int mod);
    
};


//========================== implementation

evtReader::evtReader(){
  inFile     = 0;
  data       = new DataBlock();
  
  inFileSize = 0;
  inFilePos  = 0;
  
  nBlock     = 0;    
  blockID    = -1;
  endOfFile  = false;
  isOpened   = false;
}

  
evtReader::~evtReader(){
  fclose(inFile);
  delete inFile;
  delete data;
}
  

evtReader::evtReader(TString inFileName){ 
  inFile     = 0;
  data       = new DataBlock();
  
  inFileSize = 0;
  inFilePos  = 0;

  nBlock     = 0;    
  blockID    = -1;
  endOfFile  = false;
  isOpened   = false;
  
  OpenFile(inFileName);
}

void evtReader::OpenFile(TString inFileName){
  inFile = fopen(inFileName, "r");
  if( inFile == NULL ){
    printf("Cannot read file : %s \n", inFileName.Data());
  }else{
    fseek(inFile, 0L, SEEK_END);
    inFileSize = ftell(inFile);
    rewind(inFile); ///back to the File begining

    data->Clear();
    
    gClock.Reset();
    gClock.Start("timer");
  }
};

void evtReader::UpdateFileSize(){
  if( inFile == NULL ) return;
  fseek(inFile, 0L, SEEK_END);
  inFileSize = ftell(inFile);
  fseek(inFile, inFilePos, SEEK_SET);
}

bool evtReader::IsEndOfFile() {
  int haha = feof(inFile);
  return haha > 0 ? true: false;
}


int evtReader::ReadBlock(int opt = 0){

  if( feof(inFile) ) return -1;
  if( endOfFile ) return -1;

  unsigned int header[4]; ///read 4 header, unsigned int = 4 byte = 32 bits. 

  if ( fread(header, sizeof(header), 1, inFile) != 1 ) {
      endOfFile = true;
      return -1;
  }
  blockID ++;


  if( opt == 0 ){
    /// see the Pixie-16 user manual, Table4-2
    data->eventID      = blockID;
    data->ch           =  header[0] & 0xF ;
    data->slot         = (header[0] >> 4) & 0xF;
    data->crate        = (header[0] >> 8) & 0xF;
    data->headerLength = (header[0] >> 12) & 0x1F;
    data->eventLength  = (header[0] >> 17) & 0x3FFF;
    data->pileup       =  header[0] >> 31 ;
    data->time         = ((ULong64_t)(header[2] & 0xFFFF) << 32) + header[1];
    data->cfd          =  header[2] >> 16 ; 
    data->energy       = (header[3] & 0xFFFF ); 
    data->trace_length = (header[3] >> 16) & 0x7FFF;
    data->trace_out_of_range =  header[3] >> 31;

    data->id = data->crate*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD + (data->slot-BOARD_START)*MAX_CHANNELS_PER_BOARD + data->ch;
    data->detID = mapping[data->id];

    data->ClearQDC();

    ///======== read QDCsum
    if( data->headerLength >= 4 ){
      fread(extraHeader, sizeof(unsigned int) * (data->headerLength-4), 1, inFile);
      if( data->headerLength == 8 || data->headerLength == 16){
        data->trailing = extraHeader[0];
        data->leading  = extraHeader[1];
        data->gap      = extraHeader[2];
        data->baseline = extraHeader[3];
      }
      if( data->headerLength == 12 || data->headerLength == 16){
        for( int i = 0; i < 8; i++){
          int startID = 0;
          if( data->headerLength > 12) startID = 4; ///the 1st 4 words
          data->QDCsum[i] = extraHeader[i+startID];
        }
      }
    }else{
      for( int i = 0 ; i < 8; i++){ data->QDCsum[i] = 0;}
      data->trailing = 0;
      data->leading  = 0;
      data->gap      = 0;
      data->baseline = 0;
    }
    ///====== read trace
    if( data->eventLength > data->headerLength ){      
      fread(traceBlock, sizeof(unsigned int) * ( data->trace_length / 2 ), 1, inFile);
      
      for( int i = 0; i < data->trace_length/2 ; i++){
        data->trace[2*i+0] = traceBlock[i] & 0xFFFF ;
        data->trace[2*i+1] = (traceBlock[i] >> 16 ) & 0xFFFF ;
      }
      
      ///make QDC by trace
      if( data->headerLength == 4 || data->headerLength == 8 ) {
        for( int i = 0; i < 8; i++){ data->QDCsum[i] = 0;}
        for( int i = 0; i < data->trace_length; i++){
          if(   0 <= i && i <  31 ) data->QDCsum[0] += data->trace[i];
          if(  31 <= i && i <  60 ) data->QDCsum[1] += data->trace[i];
          if(  60 <= i && i <  75 ) data->QDCsum[2] += data->trace[i];
          if(  75 <= i && i <  95 ) data->QDCsum[3] += data->trace[i];
          if(  95 <= i && i < 105 ) data->QDCsum[4] += data->trace[i];
          if( 105 <= i && i < 160 ) data->QDCsum[5] += data->trace[i];
          if( 160 <= i && i < 175 ) data->QDCsum[6] += data->trace[i];
          if( 175 <= i && i < 200 ) data->QDCsum[7] += data->trace[i];
        }
      }
    }
  }
  
  if( opt == 1 ){
    
    data->headerLength = (header[0] >> 12) & 0x1F;
    data->eventLength  = (header[0] >> 17) & 0x3FFF;
    data->trace_length = (header[3] >> 16) & 0x7FFF;
    
    if( data->headerLength >= 4 ){
      fread(extraHeader, sizeof(unsigned int) * (data->headerLength-4), 1, inFile);
    }
    if( data->eventLength > data->headerLength ){     
      fread(traceBlock, sizeof(unsigned int) * ( data->trace_length / 2 ), 1, inFile);
    }
  }
  
  inFilePos = ftell(inFile);

  return 1; 
}


void evtReader::ScanNumberOfBlock(){
  
  nBlock = 0;
  while( ReadBlock(1) != -1 ){
    nBlock ++;
    PrintStatus(10000);
  }
  
  printf("\n\n\n");
  printf("scan complete: number of data Block : %ld\n", nBlock);
  
  inFilePos = 0;
  blockID = -1;
  
  rewind(inFile); ///back to the File begining
  endOfFile = false;
  
}


void evtReader::PrintStatus(int mod){
  
  ///==== event stats, print status every 10000 events
  if ( blockID % mod == 0 ) {
    UpdateFileSize();
    gClock.Stop("timer");
    double time = gClock.GetRealTime("timer");
    gClock.Start("timer");
    printf("Total measurements: \x1B[32m%llu \x1B[0m\nReading Pos: \x1B[32m %.3f/%.3f GB\x1B[0m\nTime used:%3.0f min %5.2f sec\033[A\033[A\r", 
                 blockID, inFilePos/(1024.*1024.*1024.), inFileSize/1024./1024./1024,  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);
  }   
  
}

#endif
