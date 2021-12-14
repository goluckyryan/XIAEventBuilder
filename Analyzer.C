#define Analyzer_cxx

#include "Analyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <vector>
#include <TStopwatch.h>

//############################################ User setting

int rawEnergyRange[2] = {100,  2000};

double BGO_threshold = 100;

//############################################ end of user setting

ULong64_t NumEntries = 0;
ULong64_t ProcessedEntries = 0;
Float_t Frac = 0.1; ///Progress bar
TStopwatch StpWatch;

//############################################ histogram declaration

TH2F * heVID;
TH1F * he[NCLOVER];

TH2F * hgg[NCLOVER][NCLOVER];

TH2F * hcoin;

///----- after calibration
TH2F * heCalVID;
TH1F * heCal[NCLOVER]; // BGO veto
TH2F * hcoinBGO;

void Analyzer::Begin(TTree * tree){

   TString option = GetOption();

   NumEntries = tree->GetEntries();

   printf("======================== histogram declaration\n");
   heVID    = new TH2F("heVID",                                              "e vs ID; det ID; e [ch]", NCLOVER, 0, NCLOVER, rawEnergyRange[1] - rawEnergyRange[0], rawEnergyRange[0], rawEnergyRange[1]);
   heCalVID = new TH2F("heCalVID", Form("eCal vs ID (BGO veto > %.1f); det ID; e [ch]", BGO_threshold), NCLOVER, 0, NCLOVER, rawEnergyRange[1] - rawEnergyRange[0], rawEnergyRange[0], rawEnergyRange[1]);
   for( int i = 0; i < NCLOVER; i ++){
      he[i]    = new TH1F(   Form("he%02d", i),                                  Form("e -%02d", i), rawEnergyRange[1] - rawEnergyRange[0], rawEnergyRange[0], rawEnergyRange[1]);
      heCal[i] = new TH1F(Form("heCal%02d", i), Form("e -%02d (BGO veto > %.1f)", i, BGO_threshold), rawEnergyRange[1] - rawEnergyRange[0], rawEnergyRange[0], rawEnergyRange[1]);
   }
   
   for( int i = 0; i < NCLOVER; i++){
      for( int j = i+1; j < NCLOVER; j++){
         hgg[i][j] = new TH2F(Form("hgg%02d%02d", i, j), Form("e%02d vs e%02d; e%02d; e%02d", i, j, i, j), 
                 (rawEnergyRange[1] - rawEnergyRange[0])/2, rawEnergyRange[0], rawEnergyRange[1], 
                 (rawEnergyRange[1] - rawEnergyRange[0])/2, rawEnergyRange[0], rawEnergyRange[1]);
      }
   }
   
   hcoin = new TH2F("hcoin", "detector coin.; det ID; det ID", NCLOVER, 0, NCLOVER, NCLOVER, 0 , NCLOVER); 
   hcoinBGO = new TH2F("hcoinBGO", Form("detector coin. (BGO veto > %.1f); det ID; det ID", BGO_threshold), NCLOVER, 0, NCLOVER, NCLOVER, 0 , NCLOVER); 
   
   printf("======================== End of histograms Declaration\n");
   StpWatch.Start();

}

Bool_t Analyzer::Process(Long64_t entry){

   ProcessedEntries++;
   
   /*********** Progress Bar ******************************************/ 
   if (ProcessedEntries>NumEntries*Frac-1) {
      TString msg; msg.Form("%llu", NumEntries/1000);
      int len = msg.Sizeof();
      printf(" %3.0f%% (%*llu/%llu k) processed in %6.1f sec | expect %6.1f sec\n",
               Frac*100, len, ProcessedEntries/1000,NumEntries/1000,StpWatch.RealTime(), StpWatch.RealTime()/Frac);
      StpWatch.Start(kFALSE);
      Frac+=0.1;
   }

   b_energy->GetEntry(entry);
   b_time->GetEntry(entry);
   b_pileup->GetEntry(entry);
   b_bgo->GetEntry(entry);
   b_other->GetEntry(entry);
   b_multiplicity->GetEntry(entry);
   
   if( multi == 0 ) return kTRUE;

   for( int detID = 0; detID < NCLOVER ; detID ++){
      
      //======== baics gate when no energy or pileup 
      if( TMath::IsNaN(e[detID])) continue;
      //if( pileup[detID] == 1 ) continue;
      
      //======== Fill raw data
      heVID->Fill( detID, e[detID]);
      he[detID]->Fill(e[detID]);
      
      for( int detJ = detID +1; detJ < NCLOVER; detJ++) {
         if( TMath::IsNaN(e[detJ])) continue;
         hgg[detID][detJ]->Fill(e[detID], e[detJ]); // x then y
         hcoin->Fill(detID, detJ); 
      }
      
      //======== BGO veto
      for( int kk = 0; kk < NBGO ; kk++){
         if( TMath::IsNaN(bgo[kk]) ) continue;
         if( bgo[kk] > BGO_threshold ) {
            return kTRUE;
         }
      }
      
      //========= apply correction
      double eCal = e[detID];
      
      heCalVID->Fill( detID, eCal);
      heCal[detID]->Fill(eCal);
      
      for( int detJ = detID +1; detJ < NCLOVER; detJ++) {
         if( TMath::IsNaN(e[detJ])) continue;
         hcoinBGO->Fill(detID, detJ); 
      }
      
   }
   
   return kTRUE;
}


void listDraws(void) {
  printf("------------------- List of Plots -------------------\n");
  printf("      rawID() - Raw e vs ID\n");
  printf("      drawE() - Raw e for all %d detectors\n", NCLOVER);
  printf("-----------------------------------------------------\n");
}

void rawID(){
  TCanvas * cRawID = (TCanvas *) gROOT->FindObjectAny("cRawID");
  if( cRawID == NULL ) cRawID = new TCanvas("cRawID", "raw ID", 1000, 800);
  cRawID->cd(1)->SetGrid();
  heVID->Draw("colz");
}

void drawE(bool isLogy = false, bool cali = false){

   int numCol = NCLOVER / 4;

   TCanvas *cRawE = (TCanvas *) gROOT->FindObjectAny("cRawE");
   if( cRawE == NULL ) cRawE = new TCanvas("cRawE", cali ? "Cal e" : "Raw e", 1200, 800);
   cRawE->Clear();cRawE->Divide(numCol, 4);
   
   //cRawE->SetRightMargin(0);
   //cRawE->SetLeftMargin(0);
   //cRawE->SetTopMargin(0);
   //cRawE->SetBottomMargin(0);
   //cRawE->SetTicks(1,1);
   //cRawE->SetBorderMode(1);
   
   for (Int_t i = 0; i < 4; i++) {
      for( Int_t j = 0; j < numCol; j++){
         int canvasID = numCol * i + j + 1;
         cRawE->cd(canvasID); 
         cRawE->cd(canvasID)->SetGrid();       
         cRawE->cd(canvasID)->SetRightMargin(0.1);
         //cRawE->cd(canvasID)->SetLeftMargin(0);
         cRawE->cd(canvasID)->SetTopMargin(0);
         //cRawE->cd(canvasID)->SetBottomMargin(0);  
         //cRawE->cd(canvasID)->SetBorderMode(1);  
         if( isLogy ) cRawE->cd(canvasID)->SetLogy();
         int hID = 4*j+ i;
         if( cali ) {
            heCal[hID]->Draw("");
         }else{
            he[hID]->Draw("");
         }
      }
   }

}

void Analyzer::Terminate(){

   printf("============================== finishing.\n");
   gROOT->cd();

   int canvasXY[2] = {1200 , 1200} ;// x, y
   int canvasDiv[2] = {2,2};
   TCanvas *cCanvas  = new TCanvas("cCanvas", "" ,canvasXY[0],canvasXY[1]);
   cCanvas->Modified(); cCanvas->Update();
   cCanvas->cd(); cCanvas->Divide(canvasDiv[0],canvasDiv[1]);

   gStyle->SetOptStat("neiou");

   cCanvas->cd(1);
   cCanvas->cd(1)->SetLogz(1);
   heVID->Draw("colz");

   cCanvas->cd(2);
   cCanvas->cd(2)->SetLogz(1);
   heCalVID->Draw("colz");
   
   cCanvas->cd(3);
   cCanvas->cd(3)->SetLogz(1);
   hcoin->Draw("colz");
   
   cCanvas->cd(4);
   cCanvas->cd(4)->SetLogz(1);
   hcoinBGO->Draw("colz");
   
   listDraws();
   
   gROOT->ProcessLine(".L AutoFit.C");
   printf("=============== loaded AutoFit.C\n");
   

}
