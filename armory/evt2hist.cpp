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
#include "TH2F.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TClonesArray.h"

#include "../mapping.h"
#include "../armory/AnalysisLibrary.h"

#include "../armory/DataBlock.h"
#include "../armory/evtReader.h"

#define sleepTime 5000  ///sleep for 5 sec
//#############################TODO
//  1) multiple file
//  2) Change to GUI
//  4) eventBuilding
//  3) last 10% data


//################################## user setting
int eRange[2] = {50, 1000}; ///min, max

bool PIDFlag = false;
int GAGGID = 209;
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

  TBenchmark gClock;
  gClock.Reset();
  gClock.Start("timer");
  
  printf("====================================\n");

  evtReader * evt = new evtReader();
  evt->OpenFile(inFileName);
  if( evt->IsOpen() == false ) return -404;
  DataBlock * data = evt->data;

  printf(" in file: \033[1;31m%s\033[m\n", inFileName.Data());
  printf(" Gamma energy correction file : %s\n", corrFile == "" ? "Not provided." : corrFile.Data());
  printf(" raw E threshold : %d ch\n", rawEnergyThreshold);
  if( histFileName != "" ) printf(" Save histograms to %s\n", histFileName.Data()); 
  if( rootFileName != "" ) printf(" Save root to %s\n", rootFileName.Data()); 
  printf("--------------------------------\n");
  printf("Scanning the evt file... \n");
  evt->ScanNumberOfBlock();
  printf("go to 90%% of data \n");
  evt->JumptoPrecent(9);
  
  //================ ROOT tree
  TFile * fFile = NULL;
  TTree * tree = NULL; 
  
  short detID;
  
  if( rootFileName != "" ){
    fFile = new TFile(rootFileName, "RECREATE");
    tree = new TTree("tree", "tree");
    
    tree->Branch("headerLenght", &data->headerLength, "HeaderLength/s");
    tree->Branch("detID",                     &detID, "detID/s");
    tree->Branch("e",                  &data->energy, "energy/s");
    tree->Branch("e_t",                  &data->time, "timestamp/l");
    tree->Branch("p",                  &data->pileup, "pileup/O");
    tree->Branch("qdc",                 data->QDCsum, "QDC_sum[8]/I");
    tree->Branch("trace_length", &data->trace_length, "trace_length/s");
    tree->Branch("trace",                data->trace, "trace[trace_length]/s"); 
  }

  //================ Historgrams
  TH1F * he[NCRYSTAL];
  for( int i = 0 ; i < NCRYSTAL; i++){
    he[i] = new TH1F(Form("he%02d", i), Form("e-%2d", i), eRange[1]-eRange[0], eRange[0], eRange[1]);
    switch (i % 4){
      case 0 : he[i]->SetLineColor(2); break;
      case 1 : he[i]->SetLineColor(4); break;
      case 2 : he[i]->SetLineColor(1); break;
      case 3 : he[i]->SetLineColor(kGreen+3); break;
    }
  }
  

  TH2F * hPID ;
  if( PIDFlag ) hPID = new TH2F(Form("hPID%d", GAGGID), Form("GAGG - %d; tail; peak", GAGGID), 400, -10, 600, 400, -50, 1000); 
  
  TGraph * gTrace = new TGraph();
  TLatex text;
  text.SetNDC();
  text.SetTextFont(82);
  text.SetTextSize(0.04);
  text.SetTextColor(2);

  //================ Set Canvas
  TApplication * app = new TApplication ("app", &argn, argv);
  
  TCanvas * canvas = new TCanvas("fCanvas", "Online Spectrum", 1800, 2000);
  
  canvas->Divide(3, TMath::Ceil(NCLOVER/3.), 0);
  canvas->SetCrosshair(1);
  for( int i = 0; i < 9 ; i++){
    canvas->cd(i+1)->SetBottomMargin(0.1);
    canvas->cd(i+1)->SetRightMargin(0.002);
  }
  
  ///TCanvas * cTrace = new TCanvas("cTrace", "Trace", 100, 100, 1000, 500);
  TCanvas * cPID;
  if( PIDFlag ) cPID = new TCanvas("cPID", "PID", 100, 100, 500, 500);
  
  //=============== Read File
  int sleepCount = 0;
  
  while ( true ){

    if( evt->ReadBlock()== -1 ) {
      break;
      //printf("\n\n\nReached the end of file, wait %d sec to see any update.\n", sleepTime); 
      //sleep( sleepTime );
      //evt->UpdateFileSize(); 
      //sleepCount ++;
      //if( sleepCount > 1 ) {
      //  printf("waited for %d sec. exit.\n", 4* sleepTime);
      //  break;
      //}else{
      //  continue;
      //}
    }
    sleepCount = 0;
    
    if( evt->GetBlockID() < maxDataDisplay ) {
      printf("----------------------event Length: %u, fpos: %lu byte (%ld words)\n", data->eventLength, evt->GetFilePos(), evt->GetFilePos()/4);
      data->Print();
    }
    
    //==== Fill Histogram
    
    int haha  = data->crate*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD + (data->slot-BOARD_START)*MAX_CHANNELS_PER_BOARD + data->ch;
    detID  = mapping[haha];
    
    if( 0 <= detID && detID < 100 && data->energy > rawEnergyThreshold ){
      if( corrFile != ""){
        ///========= apply correction
        int order = (int) eCorr[detID].size();
        double eCal = 0;
        for( int k = 0; k < order ; k++){
           eCal += eCorr[detID][k] * TMath::Power(data->energy, k);
        }
        he[detID]->Fill(eCal);
      }else{
        he[detID]->Fill(data->energy);
      }
    }
    
    ///============ QDC 
    if( PIDFlag && detID == GAGGID && (data->headerLength < data->eventLength) ){
        double bg   = (data->QDCsum[0] + data->QDCsum[1])/60.;
        double peak = data->QDCsum[3]/20. - bg;
        double tail = data->QDCsum[5]/55. - bg;
        hPID->Fill( tail , peak);
    }
    
    //===== Trace
    ///if( data->trace_length > 0 ) {
    ///  gTrace->Clear();
    ///  gTrace->Set(data->trace_length);
    ///  gTrace->SetTitle(Form("eventID : %llu, detID: %d", data->eventID, data->detID));
    ///  
    ///  for( int i = 0; i < data->trace_length; i++) gTrace->SetPoint(i, i, data->trace[i]);
    ///}
    
    if( rootFileName != "" ){
      fFile->cd();
      tree->Fill();
    }
    
    //==== event stats, print status every 10000 events
    evt->PrintStatus(10000);   
    
    //==== Plot Canvas
    gClock.Stop("timer");
    int time = TMath::Floor(gClock.GetRealTime("timer")*1000); // in millisec
    gClock.Start("timer");
    if( time  % 1000 == 0 || time < 10){
      
      ///==== for clover
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
      
      ///==== for trace
      ///if( data->trace_length > 0 ){
      ///  cTrace->cd();
      ///  gTrace->Draw("AL");
      ///  
      ///  for( int i = 0; i < 8; i++){
      ///    text.DrawLatex(0.2, 0.8-0.05*i,  Form("%d", data->QDCsum[i]));
      ///  }
      ///  cTrace->Modified();
      ///  cTrace->Update();
      ///}
      
      ///=== for GAGG PID
      if( PIDFlag ) {
        cPID->cd();
        cPID->SetLogz();
        hPID->Draw("colz");
        cPID->Modified();
        cPID->Update();
      } 
      
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
  
  evt->PrintStatus(1);
  
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
