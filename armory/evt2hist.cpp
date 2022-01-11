#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TBenchmark.h"
#include "TH1F.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "../mapping.h"
#include "../armory/AnalysisLibrary.h"

#define MAX_CRATES 2
#define MAX_BOARDS_PER_CRATE 13
#define MAX_CHANNELS_PER_BOARD 16
#define BOARD_START 2

//#############################TODO
//  1) Get ADC data
//  2) Change to GUI
//  3) calibration gamma

int rawEnergyThreshold = 100;

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
  
  Int_t            trailing;
  Int_t             leading;
  Int_t                 gap;
  Int_t            baseline;
  Int_t              ADC[8];
  
  UShort_t id;
  Int_t  detID;
  ULong64_t eventID;
  
  UShort_t    trace[1024];
  
  measurment(){
    Clear();
  };
  
  ~measurment(){};
  
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
    id = 0;
    detID = -1;
    eventID = 0;
    trailing = 0;
    leading = 0;
    gap = 0;
    baseline = 0;
    for( int i = 0; i < 8; i++) ADC[i] = -1;
    for( int i = 0; i < 1024; i++) trace[i] = 0;
  }
  
  void Print(){
    printf("============== eventID : %llu\n", eventID);
    printf("Crate: %d, Slot: %d, Ch: %d | id: %d = detID : %d \n", crate, slot, ch, id, detID);
    printf("HeaderLength: %d, Event Length: %d, energy: %d, timeStamp: %llu\n", headerLength, eventLength, energy, time);
    printf("trace_length: %d, pile-up:%d\n", trace_length, pileup); 
    if( headerLength > 4 ){
      if( headerLength > 12 ){
        printf(" trailing : %d\n", trailing);
        printf(" leading  : %d\n", leading);
        printf(" gap      : %d\n", gap);
        printf(" baseLine : %d\n", baseline);
      }
      printf(" ADC : \n");
      for( int i = 0; i < 8; i++) printf("    %-10d\n", ADC[i]);
    }
    if( eventLength > headerLength ){
      printf(" trace:\n");
      for( int i = 0 ; i < trace_length ; i++)printf("%3d|     %-10d\n",i, trace[i]);
    }
  }
  
};

//#############################################
//           main 
//###############################################
int main(int argn, char **argv) {
    
  if (argn != 2 && argn != 3 )    {
    printf("Usage :\n");
    printf("%s [evt File]  [E corr]\n", argv[0]);
    printf("         [E corr] : correction file for gamma energy \n");       
    return 1;
  }
  
  TString inFileName = argv[1];

  TString corrFile = "";
  std::vector<std::vector<double>> eCorr;
  if( argn >= 3 ){
    corrFile = argv[2];
    eCorr.clear();
    eCorr = LoadCorrectionParameters(corrFile);
  }
  
  long int inFilePos;
  TBenchmark gClock;
  gClock.Reset();
  gClock.Start("timer");
  
  ULong64_t measureID = -1;
  
  measurment data;
  
  printf("====================================\n");

  FILE * inFile = fopen(inFileName, "r");
  if( inFile == NULL ){
    printf("Cannot read file : %s \n", inFileName.Data());
    return -404;
  }

  printf(" in file: \003[1;31m%s\033[m\n", inFileName.Data());
  printf(" Gamma energy correction file : %s\n", corrFile == "" ? "Not provided." : corrFile.Data());
  printf("--------------------------------\n");
  

  //================ get file size
  fseek(inFile, 0L, SEEK_END);
  long int inFileSize = ftell(inFile);
  rewind(inFile); ///back to the File begining
  unsigned long long fpos = 0;

  //================ Historgrams
  TH1F * he[NCRYSTAL];
  for( int i = 0 ; i < NCRYSTAL; i++){
    he[i] = new TH1F(Form("he%02d", i), Form("e-%2d", i), 1000, 0, 2000);
    switch (i % 4){
      case 0 : he[i]->SetLineColor(2); break;
      case 1 : he[i]->SetLineColor(4); break;
      case 2 : he[i]->SetLineColor(1); break;
      case 3 : he[i]->SetLineColor(kGreen+3); break;
    }
  }

  //================ Set Canvas
  TApplication * app = new TApplication ("app", &argn, argv);
  
  
  TCanvas * canvas = new TCanvas("fCanvas", "Online Spectrum", 1800, 1400);
  canvas->Divide(1, 9, 0);
  canvas->SetCrosshair(1);
  for( int i = 0; i < 9 ; i++){
    canvas->cd(i+1)->SetBottomMargin(0.1);
    canvas->cd(i+1)->SetRightMargin(0.002);
  }
  
  

  //=============== Read File
  unsigned int header[4]; //read 4 header, unsigned int = 4 byte = 32 bits.  

  while ( ! feof(inFile) ){

    fread(header, sizeof(header), 1, inFile);
    measureID ++;
    
    /// see the Pixie-16 user manual, Table4-2
    data.eventID      = measureID;
    data.ch           =  header[0] & 0xF ;
    data.slot         = (header[0] >> 4) & 0xF;
    data.crate        = (header[0] >> 8) & 0xF;
    data.headerLength = (header[0] >> 12) & 0x1F;
    data.eventLength  = (header[0] >> 17) & 0x3FFF;
    data.pileup       =  header[0] >> 31 ;
    data.time         = ((ULong64_t)(header[2] & 0xFFFF) << 32) + header[1];
    data.cfd          =  header[2] >> 16 ; 
    data.energy       = (header[3] & 0xFFFF ) /2; // I don;t know why it has to "rebin"
    data.trace_length = (header[3] >> 16) & 0x7FFF;
    data.trace_out_of_range =  header[3] >> 31;
    
    data.id = data.crate*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD + (data.slot-BOARD_START)*MAX_CHANNELS_PER_BOARD + data.ch;
    data.detID = mapping[data.id];
    
    ///======== read ADC
    if( data.headerLength >= 4 ){
      unsigned int extraHeader[data.headerLength-4];
      fread(extraHeader, sizeof(extraHeader), 1, inFile);
      //if( measureID < 10 ) for(int i = 0; i < data.headerLength - 4; i++) printf("  %x\n", extraHeader[i]);      
      if( data.headerLength > 12){
        data.trailing = extraHeader[0];
        data.leading  = extraHeader[1];
        data.gap      = extraHeader[2];
        data.baseline = extraHeader[3];
      }
      
      for( int i = 0; i < 8; i++){
        int startID = 0;
        if( data.headerLength > 12) startID = 4; ///the 1st 4 words
        data.ADC[i] = extraHeader[i+startID];
      }
    }
    ///====== read trace
    if( data.eventLength > data.headerLength ){      
      unsigned int traceBlock[data.trace_length / 2];
      fread(traceBlock, sizeof(traceBlock), 1, inFile);
      
      for( int i = 0; i < data.trace_length/2 ; i++){
        data.trace[2*i+0] = traceBlock[i] & 0xFFFF ;
        data.trace[2*i+1] = (traceBlock[i] >> 16 ) & 0xFFFF ;
      }
    }
    
    if( measureID < 10 ) {
      printf("----------------------event Length: %u, fpos: %llu byte (%lld words)\n", data.eventLength, fpos, fpos/4);
      for(int i = 0; i < 4; i++) printf("  %x\n", header[i]);
      data.Print();
    }
    
    //=== jump to next measurement. obsolete, we read the whole block
    ///if( data.eventLength > 4 ) {
    ///  if( fseek(inFile, sizeof(int) * (data.eventLength-data.headerLength),  SEEK_CUR) != 0 ) break;;
    ///}
    fpos = ftell(inFile);
    
    
    //==== Fill Histogram
    if( 0 <= data.detID && data.detID < 100 && data.energy > rawEnergyThreshold ){
      if( corrFile != ""){
        ///========= apply correction
        int order = (int) eCorr[data.detID].size();
        double eCal = 0;
        for( int k = 0; k < order ; k++){
           eCal += eCorr[data.detID][k] * TMath::Power(data.energy, k);
        }
        he[data.detID]->Fill(eCal);
      }else{
        he[data.detID]->Fill(data.energy);
      }
    }
    
    
    
    //==== event stats, print status every 10000 events
    if ( measureID % 10000 == 0 ) {
      /// update file size, slow down?
      fseek(inFile, 0L, SEEK_END);
      inFileSize = ftell(inFile);
      fseek(inFile, fpos, SEEK_SET);
      
      inFilePos = ftell(inFile);
      gClock.Stop("timer");
      double time = gClock.GetRealTime("timer");
      gClock.Start("timer");
      printf("Total measurements: \x1B[32m%llu \x1B[0m\nReading Pos: \x1B[32m %.3f/%.3f GB\x1B[0m\nTime used:%3.0f min %5.2f sec\033[A\033[A\r", 
                   measureID, inFilePos/(1024.*1024.*1024.), inFileSize/1024./1024./1024,  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);
    }   
    
    //==== Plot Canvas
    gClock.Stop("timer");
    int time = TMath::Floor(gClock.GetRealTime("timer")*1000); // in millisec
    gClock.Start("timer");
    if( time  % 1000 == 0 ){
      for( int i = 0; i < NCRYSTAL; i++){
        canvas->cd(i/4 +1);
        //canvas->cd(i/4 +1)->SetLogy();
        if( i % 4 == 0  ) {
          he[i]->Draw();
        }else{
          he[i]->Draw("same");
        }
      }
      canvas->Modified();
      canvas->Update(); 
      gSystem->ProcessEvents();
    }  
  }//---- end of file loop
  
  
  for( int i = 0; i < NCRYSTAL; i++){
    canvas->cd(i/4 +1);
    //canvas->cd(i/4 +1)->SetLogy();
    if( i % 4 == 0  ) {
      he[i]->Draw();
    }else{
      he[i]->Draw("same");
    }
  }
  canvas->Modified();
  canvas->Update(); 
  gSystem->ProcessEvents();
  
  
  inFilePos = ftell(inFile);
  gClock.Stop("timer");
  double time = gClock.GetRealTime("timer");
  gClock.Start("timer");
  printf("Total measurements: \x1B[32m%llu \x1B[0m\nReading Pos: \x1B[32m %.3f/%.3f GB\x1B[0m\nTime used:%3.0f min %5.2f sec\033[A\033[A\r", 
         measureID, inFilePos/(1024.*1024.*1024.), inFileSize/1024./1024./1024,  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);

  fclose(inFile);
  
  printf("\n============= reasched end of file\n");
  printf("\nCrtl+C to end program.\n");
  app->Run();


}
