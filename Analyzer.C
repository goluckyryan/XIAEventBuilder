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

int rawEnergyRange[2] = {100,  6000}; // in ch
int energyRange[2] = {100, 2000}; // keV

double BGO_threshold = 100; // in ch

//############################################ end of user setting

ULong64_t NumEntries = 0;
ULong64_t ProcessedEntries = 0;
Float_t Frac = 0.1; ///Progress bar
TStopwatch StpWatch;

//############################################ histogram declaration

TH2F * heVID;
TH1F * he[NCLOVER];

TH1F * h1, * h2;

TH2F * hgg[NCLOVER][NCLOVER];

TH2F * hcoin;

///----- after calibration and BGO veto
TH2F * heCalVID;
TH1F * heCal[NCLOVER]; 
TH2F * hcoinBGO;

void Analyzer::Begin(TTree * tree){

   TString option = GetOption();

   NumEntries = tree->GetEntries();

   printf("======================== histogram declaration\n");

   h1 = new TH1F("h1", "h1", 1900, 100, 2000);
   h2 = new TH1F("h2", "h2 BGO gated", 1900, 100, 2000);
   
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
   //b_pileup->GetEntry(entry);
   b_bgo->GetEntry(entry);
   //b_other->GetEntry(entry);
   //b_multiplicity->GetEntry(entry);
   
   if( multi == 0 ) return kTRUE;

   ///=========== Looping Crystals
   for( int detID = 0; detID < NCLOVER ; detID ++){
      
      //======== baics gate when no energy or pileup 
      if( TMath::IsNaN(e[detID])) continue;
      //if( pileup[detID] == 1 ) continue;
      
      //======== Fill raw data
      heVID->Fill( detID, e[detID]);

      he[detID]->Fill(e[detID]);

      if( 8 == detID ) {
	h1->Fill( e[detID]*0.307484 - 0.505163 );
	heCal[detID]->Fill( e[detID]*0.307484 - 0.505163 );
      }
      if( 9 == detID ) {
	h1->Fill( e[detID]*0.308628 + 0.672629 );
	heCal[detID]->Fill( e[detID]*0.308628 + 0.672629 );
      }
      if( 10 == detID ) {
	h1->Fill( e[detID]*0.308445 + 0.238095 );
	heCal[detID]->Fill( e[detID]*0.308445 + 0.238095 );
      }
      if( 11 == detID ) {
	h1->Fill( e[detID]*0.312665 + 0.359117 );
	heCal[detID]->Fill( e[detID]*0.312665 + 0.359117 );
      }
      
      
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


      
      if(  8 == detID ) h2->Fill( e[detID]*0.307484 - 0.505163 );
      if(  9 == detID ) h2->Fill( e[detID]*0.308628 + 0.672629 );
      if( 10 == detID ) h2->Fill( e[detID]*0.308445 + 0.238095 );
      if( 11 == detID ) h2->Fill( e[detID]*0.312665 + 0.359117 );

      
      
      //========= apply correction
      double eCal = e[detID];
      
      heCalVID->Fill( detID, eCal);
      //heCal[detID]->Fill(eCal);
      
      for( int detJ = detID +1; detJ < NCLOVER; detJ++) {
         if( TMath::IsNaN(e[detJ])) continue;
         hcoinBGO->Fill(detID, detJ); 
      }
      
   }
   
   return kTRUE;
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
   //cCanvas->cd(3)->SetLogz(1);
   //hcoin->Draw("colz");
   h1->Draw("");
   heCal[8]->SetLineColor(2);heCal[8]->Draw("same");
   heCal[9]->SetLineColor(4);heCal[9]->Draw("same");
   heCal[10]->SetLineColor(6);heCal[10]->Draw("same");
   heCal[11]->SetLineColor(7);heCal[11]->Draw("same");
   
   cCanvas->cd(4);
   //cCanvas->cd(4)->SetLogz(1);
   //hcoinBGO->Draw("colz");
   h2->Draw("");
   
   printf("=============== Analyzer Utility\n");
   gROOT->ProcessLine(".L Analyzer_Utilt.c");
   //gROOT->ProcessLine("listDraws()");
   
   printf("=============== loaded AutoFit.C\n");
   gROOT->ProcessLine(".L AutoFit.C");
   

}
