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
  std::vector<std::vector<double>> corr;
  if( argn >= 3 ){
    corrFile = argv[2];
    corr.clear();
    corr = LoadCorrectionParameters(corrFile);
  }
  
  long int inFilePos;
  TBenchmark gClock;
  gClock.Reset();
  gClock.Start("timer");
  
  ULong64_t measureID = 0;
  
  measurment data;
  
  printf("====================================\n");

  FILE * inFile = fopen(inFileName, "r");
  if( inFile == NULL ){
    printf("Cannot read file : %s \n", inFileName.Data());
    return -404;
  }

  //get file size
  //fseek(inFile, 0L, SEEK_END);
  //long int inFileSize = ftell(inFile);
  //rewind(inFile); ///back to the File begining

  printf(" in file: %s\n", inFileName.Data());
  printf(" Gamma energy correction file : %s\n", corrFile == "" ? "Not provided." : corrFile.Data());
  printf("--------------------------------\n");
  
  //================ Historgrams
  TH1F * he[NCRYSTAL];
  for( int i = 0 ; i < NCRYSTAL; i++){
    he[i] = new TH1F(Form("he%02d", i), Form("e-%2d", i), 1000, 0, 6000);
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
    canvas->cd(i+1)->SetBottomMargin(0.06);
    canvas->cd(i+1)->SetRightMargin(0.002);
  }

  unsigned int header[4]; //read 4 header, unsigned int = 4 byte = 32 bits.  
  unsigned long long nWords = 0;
  unsigned long long fpos = 0;

  //=============== Read File
  while ( ! feof(inFile) ){

    fread(header, sizeof(header), 1, inFile);
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
    
    ///printf("----------------------nWords: %llu, fpos: %llu\n", nWords, fpos);
    ///for(int i = 0; i < 4; i++){
    ///  printf("  %x\n", header[i]);
    ///}
    ///data.Print();
    
    int detID = map[data.id];
    if( 0 <= detID && detID < 100 ){
      if( corrFile != ""){
        double eCal = ApplyCorrection(corr, detID, data.energy);      
        he[detID]->Fill(eCal);
      }else{
        he[detID]->Fill(data.energy);
      }
    }
    
    //=== jump to next measurement
    if( data.eventLength > 4 ){
      fseek(inFile, sizeof(int) * (data.eventLength-4),  SEEK_CUR);
      fpos += ftell(inFile);
    }
    
    //==== event stats, print status every 10000 events
    if ( measureID % 10000 == 0 ) {
      //fseek(inFile, 0L, SEEK_END);
      //inFileSize = ftell(inFile);
      //fseek(inFile, fpos, SEEK_SET);
      
      inFilePos = ftell(inFile);
      float tempf = (float)inFilePos/(1024.*1024.*1024.);
      gClock.Stop("timer");
      double time = gClock.GetRealTime("timer");
      gClock.Start("timer");
      printf("Total measurements: \x1B[32m%llu \x1B[0m\nData Read: \x1B[32m %.3f GB\x1B[0m\nTime used:%3.0f min %5.2f sec\033[A\033[A\r", 
                   measureID, tempf,  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);
    }   
    
    
    gClock.Stop("timer");
    int time = TMath::Floor(gClock.GetRealTime("timer")*1000); // in millisec
    gClock.Start("timer");
    if( time  % 1000 == 0 ){
      for( int i = 0; i < NCRYSTAL; i++){
        canvas->cd(i/4 +1);
        canvas->cd(i/4 +1)->SetLogy();
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
  }
  
  inFilePos = ftell(inFile);
  gClock.Stop("timer");
  double time = gClock.GetRealTime("timer");
  gClock.Start("timer");
  float tempf = (float)inFilePos/(1024.*1024.*1024.);
  printf("Total measurements: \x1B[32m%llu \x1B[0m\nData Read: \x1B[32m %.3f GB\x1B[0m\nTime used:%3.0f min %5.2f sec\033[A\033[A\r", 
                 measureID, tempf,  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);

  fclose(inFile);
  
  printf("\n============= reasched end of file\n");
  printf("Crtl+C to end program.\n");
  app->Run();


}
