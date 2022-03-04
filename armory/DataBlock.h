
#ifndef DATABLOCK_H
#define DATABLOCK_H

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

#define MAX_TRACE_LENGHT 16000

class DataBlock{

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
  
  Int_t            trailing;
  Int_t             leading;
  Int_t                 gap;
  Int_t            baseline;
  Int_t           QDCsum[8];
  
  ULong64_t eventID;
  
  UShort_t  trace[MAX_TRACE_LENGHT];
  
  DataBlock(){
    Clear();
  };
  
  ~DataBlock(){};
  
  void Clear(){
    ch = 0;
    slot = 0;
    crate = 0;
    headerLength = 0;
    eventLength = 0;
    pileup = false;
    time = 0;
    cfd = 0;
    energy = 0;
    trace_length = 0;
    trace_out_of_range = 0;
    eventID = 0;
    ClearQDC();
    ClearTrace();
  }
  
  void ClearQDC(){
    trailing = 0;
    leading = 0;
    gap = 0;
    baseline = 0;
    for( int i = 0; i < 8; i++) QDCsum[i] = -1;
  }
  
  void ClearTrace(){
    for( int i = 0 ; i < MAX_TRACE_LENGHT; i++) trace[i] = 0;
  }
  
  
  void Print(){
    printf("============== eventID : %llu\n", eventID);
    printf("Crate: %d, Slot: %d, Ch: %d \n", crate, slot, ch);
    printf("HeaderLength: %d, Event Length: %d, energy: %d, timeStamp: %llu\n", headerLength, eventLength, energy, time);
    printf("trace_length: %d, pile-up:%d\n", trace_length, pileup); 
    if( headerLength > 4 ){
      if( headerLength > 12 ){
        printf(" trailing : %d\n", trailing);
        printf(" leading  : %d\n", leading);
        printf(" gap      : %d\n", gap);
        printf(" baseLine : %d\n", baseline);
      }
      printf(" QDCsum : \n");
      for( int i = 0; i < 8; i++) printf("    %-10d\n", QDCsum[i]);
    }
    if( eventLength > headerLength ){
      printf(" trace:\n");
      for( int i = 0 ; i < trace_length ; i++)printf("%3d|     %-10d\n",i, trace[i]);
    }
  }
  
};

#endif
