#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <thread>

#include "TSystem.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TBenchmark.h"
#include "TH1F.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TClonesArray.h"

#include "../mapping.h"
#include "../armory/AnalysisLibrary.h"

#include "../armory/DataBlock.h"

#define MAX_CRATES 2
#define MAX_BOARDS_PER_CRATE 13
#define MAX_CHANNELS_PER_BOARD 16
#define BOARD_START 2

//#############################TODO
//  1) multiple file
//  2) Change to GUI
//  4) eventBuilding

DataBlock dataList[1000];
int dataCollected = 0;

void BuildEvent(){
  
  ///==== sort timestamp
  
  ///==== build events
  
  ///==== output to g-g hist
}

//#############################################
//           main 
//###############################################
int main(int argn, char **argv) {
    
  if (argn < 2 || argn > 7 )    {
    printf("Usage :\n");
    printf("%s [evt File]  [E corr] [raw E threshold] [Hist File] [Root File] [display data]\n", argv[0]);
    printf("          [E corr] : correction file for gamma energy \n");       
    printf(" [raw E threshold] : min raw E, default 10 ch \n");       
    printf("       [Hist File] : if provided, saved histograms in root. if value = 1, auto fileName. \n");       
    printf("       [Root File] : if provided, saved energy, timestamp, QDC, and trace in to root. if value = 1, auto fileName. \n");       
    printf("    [display data] : number of event from the first one to display. default 0. \n");
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

  int rawEnergyThreshold = 10;
  if( argn >= 4 ) rawEnergyThreshold = atoi(argv[3]);

  TString histFileName = ""; ///save gamma hist for calibration
  if( argn >= 5 ) histFileName = argv[4];
  if( histFileName == "1" ){
    histFileName = inFileName;
    histFileName.Remove(0, inFileName.Last('/')+1);
    histFileName.Remove(histFileName.First('.'));
    histFileName.Append("_hist.root");
  }
  if( histFileName == "0" ) histFileName = "";

  TString rootFileName = ""; ///save data into root
  if( argn >= 6 ) rootFileName = argv[5];
  if( rootFileName == "1" ){
    rootFileName = inFileName;
    rootFileName.Remove(0, inFileName.Last('/')+1);
    rootFileName.Remove(rootFileName.First('.'));
    rootFileName.Append("_raw.root");
  }
  if( rootFileName == "0" ) rootFileName = "";
  
  int maxDataDisplay = 0;
  if( argn >= 7 ) maxDataDisplay = atoi(argv[6]);

  long int inFilePos;
  TBenchmark gClock;
  gClock.Reset();
  gClock.Start("timer");
  
  ULong64_t measureID = -1;
  
  DataBlock data;
  
  printf("====================================\n");

  FILE * inFile = fopen(inFileName, "r");
  if( inFile == NULL ){
    printf("Cannot read file : %s \n", inFileName.Data());
    return -404;
  }

  printf(" in file: \033[1;31m%s\033[m\n", inFileName.Data());
  printf(" Gamma energy correction file : %s\n", corrFile == "" ? "Not provided." : corrFile.Data());
  printf(" raw E threshold : %d ch\n", rawEnergyThreshold);
  if( histFileName != "" ) printf(" Save histograms to %s\n", histFileName.Data()); 
  if( rootFileName != "" ) printf(" Save root to %s\n", rootFileName.Data()); 
  printf("--------------------------------\n");
  
  TFile * fFile = NULL;
  TTree * tree = NULL; 
  if( rootFileName != "" ){
    fFile = new TFile(rootFileName, "RECREATE");
    tree = new TTree("tree", "tree");
    
    
    tree->Branch("headerLenght", &data.headerLength, "HeaderLength/s");
    tree->Branch("detID",               &data.detID, "detID/s");
    tree->Branch("id",                     &data.id, "id/s");
    tree->Branch("e",                  &data.energy, "energy/s");
    tree->Branch("e_t",                  &data.time, "timestamp/l");
    tree->Branch("p",                  &data.pileup, "pileup/O");
    tree->Branch("qdc",                 data.QDCsum, "QDC_sum[8]/I");
    tree->Branch("trace_length", &data.trace_length, "trace_length/s");
    tree->Branch("trace",                data.trace, "trace[trace_length]/s"); 
  }

  //================ get file size
  fseek(inFile, 0L, SEEK_END);
  long int inFileSize = ftell(inFile);
  rewind(inFile); ///back to the File begining
  unsigned long long fpos = 0;

  //================ Historgrams
  TH1F * he[NCRYSTAL];
  for( int i = 0 ; i < NCRYSTAL; i++){
    he[i] = new TH1F(Form("he%02d", i), Form("e-%2d", i), 2000, 0, 2000);
    switch (i % 4){
      case 0 : he[i]->SetLineColor(2); break;
      case 1 : he[i]->SetLineColor(4); break;
      case 2 : he[i]->SetLineColor(1); break;
      case 3 : he[i]->SetLineColor(kGreen+3); break;
    }
  }
  
  TGraph * gTrace = new TGraph();
  TLatex text;
  text.SetNDC();
  text.SetTextFont(82);
  text.SetTextSize(0.04);
  text.SetTextColor(2);

  //================ Set Canvas
  TApplication * app = new TApplication ("app", &argn, argv);
  
  TCanvas * canvas = new TCanvas("fCanvas", "Online Spectrum", 1800, 2000);
  
  canvas->Divide(1, 9, 0);
  canvas->SetCrosshair(1);
  for( int i = 0; i < 9 ; i++){
    canvas->cd(i+1)->SetBottomMargin(0.1);
    canvas->cd(i+1)->SetRightMargin(0.002);
  }
  
  ///TCanvas * cTrace = new TCanvas("cTrace", "Trace", 100, 100, 1000, 500);
  

  //=============== Read File
  unsigned int header[4]; //read 4 header, unsigned int = 4 byte = 32 bits.  

  while ( ! feof(inFile) ){

    if ( fread(header, sizeof(header), 1, inFile) !=1 ) break;
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
    
    ///======== read QDCsum
    if( data.headerLength >= 4 ){
      unsigned int extraHeader[data.headerLength-4];
      fread(extraHeader, sizeof(extraHeader), 1, inFile);
      if( data.headerLength == 8 || data.headerLength == 16){
        data.trailing = extraHeader[0];
        data.leading  = extraHeader[1];
        data.gap      = extraHeader[2];
        data.baseline = extraHeader[3];
      }
      if( data.headerLength == 12 || data.headerLength == 16){
          for( int i = 0; i < 8; i++){
            int startID = 0;
            if( data.headerLength > 12) startID = 4; ///the 1st 4 words
            data.QDCsum[i] = extraHeader[i+startID];
          }
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
    
    if( measureID < maxDataDisplay ) {
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
    
    
    //===== Trace
    if( data.trace_length > 0 ) {
      gTrace->Clear();
      gTrace->Set(data.trace_length);
      gTrace->SetTitle(Form("eventID : %llu, detID: %d", data.eventID, data.detID));
      
      if( data.headerLength == 4 ) {
        for( int i = 0; i < 8; i++ ) data.QDCsum[i] = 0;
      }
      
      for( int i = 0; i < data.trace_length; i++){
        gTrace->SetPoint(i, i, data.trace[i]);

        ///if the header don't have ADC, make one
        if( data.headerLength < 12 ) {
          if(   0 <= i && i <  31 ) data.QDCsum[0] += data.trace[i];
          if(  31 <= i && i <  60 ) data.QDCsum[1] += data.trace[i];
          if(  60 <= i && i <  75 ) data.QDCsum[2] += data.trace[i];
          if(  75 <= i && i <  95 ) data.QDCsum[3] += data.trace[i];
          if(  95 <= i && i < 105 ) data.QDCsum[4] += data.trace[i];
          if( 105 <= i && i < 160 ) data.QDCsum[5] += data.trace[i];
          if( 160 <= i && i < 175 ) data.QDCsum[6] += data.trace[i];
          if( 175 <= i && i < 200 ) data.QDCsum[7] += data.trace[i];
        }
      }
    }
    
    if( rootFileName != "" ){
      fFile->cd();
      tree->Fill();
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
    if( time  % 1000 == 0 || time < 10){
      
      //==== for clover
      for( int i = 0; i < NCLOVER; i++){
        double maxY = 0;
        double y = 0;
        for( int j = 0; j < 4; j++){
          int mBin = he[4*i+j]->GetMaximumBin();
          y = he[4*i+j]->GetBinContent(mBin);
          if( maxY < y ) maxY = y;
        }
        for( int j = 0; j < 4; j++){
          canvas->cd(i+1);
          he[4*i+j]->GetYaxis()->SetRangeUser(0, maxY*1.2);
          if ( j ==  0) {
            he[4*i]->Draw();
          }else{
            he[4*i+j]->Draw("same");
          }
        }
      }
      canvas->Modified();
      canvas->Update(); 
      
      //==== for trace
      ///if( data.trace_length > 0 ){
      ///  cTrace->cd();
      ///  gTrace->Draw("AL");
      ///  
      ///  for( int i = 0; i < 8; i++){
      ///    text.DrawLatex(0.2, 0.8-0.05*i,  Form("%d", data.QDCsum[i]));
      ///  }
      ///  cTrace->Modified();
      ///  cTrace->Update();
      ///}
      
      
      gSystem->ProcessEvents();
    }  
  }//---- end of file loop
  
  
  for( int i = 0; i < NCLOVER; i++){
    double maxY = 0;
    double y = 0;
    for( int j = 0; j < 4; j++){
      int mBin = he[4*i+j]->GetMaximumBin();
      y = he[4*i+j]->GetBinContent(mBin);
      if( maxY < y ) maxY = y;
    }
    for( int j = 0; j < 4; j++){
      canvas->cd(i+1);
      he[4*i+j]->GetYaxis()->SetRangeUser(0, maxY*1.2);
      if ( j ==  0) {
        he[4*i]->Draw();
      }else{
        he[4*i+j]->Draw("same");
      }
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
  
  printf("\n\n\n============= reached end of file\n");
  
  if( histFileName != "" ) {
    printf(" save gamma histograms : \033[1;3m%s\033[m\n", histFileName.Data());
    TFile * fHist = new TFile(histFileName, "RECREATE");
    for( int i = 0; i < NCRYSTAL; i++) he[i]->Write("", TObject::kOverwrite);
    fHist->Close();
  }
  
  if( rootFileName != "" ){
    printf(" save into Root : \033[1;3m%s\033[m\n", rootFileName.Data());
    fFile->cd();
    tree->Write();
    fFile->Close();
  }
  
  printf("Crtl+C to end program.\n");
  
  app->Run();


}
