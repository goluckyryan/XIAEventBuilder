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

int rawEnergyRange[2] = {500,  6000}; // in ch
int energyRange[3] = {1, 100, 2000}; // keV {resol, min, max}

double BGO_threshold = 100; // in ch

TString e_corr = "correction_e.dat";

//############################################ end of user setting

ULong64_t NumEntries = 0;
ULong64_t ProcessedEntries = 0;
Float_t Frac = 0.1; ///Progress bar
TStopwatch StpWatch;

vector<vector<double>> eCorr;

//############################################ histogram declaration

TH2F * heVID;
TH1F * he[NCRYSTAL];

TH2F * hgg[NCRYSTAL][NCRYSTAL];

TH2F * hcoin;

///----- after calibration and BGO veto
TH2F * heCalVID;
TH1F * heCal[NCRYSTAL]; 
TH2F * hcoinBGO;

//############################################ BEGIN
void Analyzer::Begin(TTree * tree){

   TString option = GetOption();

   NumEntries = tree->GetEntries();

   printf("======================== Histograms declaration\n");
   
   heVID    = new TH2F("heVID",                                              "e vs ID; det ID; e [ch]", NCRYSTAL, 0, NCRYSTAL, rawEnergyRange[1] - rawEnergyRange[0], rawEnergyRange[0], rawEnergyRange[1]);
   heCalVID = new TH2F("heCalVID", Form("eCal vs ID (BGO veto > %.1f); det ID; Energy [keV]", BGO_threshold), NCRYSTAL, 0, NCRYSTAL, (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
   for( int i = 0; i < NCRYSTAL; i ++){
      he[i]    = new TH1F(   Form("he%02d", i),                                  Form("e -%02d", i), rawEnergyRange[1] - rawEnergyRange[0], rawEnergyRange[0], rawEnergyRange[1]);
      heCal[i] = new TH1F(Form("heCal%02d", i), Form("eCal -%02d (BGO veto > %.1f); Energy [keV];  count / %d keV", i, BGO_threshold, energyRange[0]), (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
   }
   
   for( int i = 0; i < NCRYSTAL; i++){
      for( int j = i+1; j < NCRYSTAL; j++){
         //hgg[i][j] = new TH2F(Form("hgg%02d%02d", i, j), Form("e%02d vs e%02d; e%02d; e%02d", i, j, i, j), 
         //        (rawEnergyRange[1] - rawEnergyRange[0])/2, rawEnergyRange[0], rawEnergyRange[1], 
         //        (rawEnergyRange[1] - rawEnergyRange[0])/2, rawEnergyRange[0], rawEnergyRange[1]);
      }
   }
   
   hcoin = new TH2F("hcoin", "detector coin.; det ID; det ID", NCRYSTAL, 0, NCRYSTAL, NCRYSTAL, 0 , NCRYSTAL); 
   hcoinBGO = new TH2F("hcoinBGO", Form("detector coin. (BGO veto > %.1f); det ID; det ID", BGO_threshold), NCRYSTAL, 0, NCRYSTAL, NCRYSTAL, 0 , NCRYSTAL); 
   
   printf("======================== End of histograms declaration\n");
   
   printf("======================== Load parameters.\n");
   
   eCorr = LoadCorrectionParameters(e_corr); 
   
   StpWatch.Start();
   printf("======================== Start processing....\n");

}
//############################################ PROCESS
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
   for( int detID = 0; detID < NCRYSTAL ; detID ++){
      
      //======== baics gate when no energy or pileup 
      if( TMath::IsNaN(e[detID])) continue;
      //if( pileup[detID] == 1 ) continue;
      
      //======== Fill raw data
      heVID->Fill( detID, e[detID]);
      he[detID]->Fill(e[detID]);
      
      
      for( int detJ = detID +1; detJ < NCRYSTAL; detJ++) {
         if( TMath::IsNaN(e[detJ])) continue;
         //hgg[detID][detJ]->Fill(e[detID], e[detJ]); // x then y
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
      double eCal = 0 ;
      int order = (int) eCorr[detID].size();
      for( int i = 0; i < order ; i++){
         eCal += eCorr[detID][i] * TMath::Power(e[detID], i);
      }
      
      
      heCalVID->Fill( detID, eCal);
      heCal[detID]->Fill(eCal);
      
      for( int detJ = detID +1; detJ < NCRYSTAL; detJ++) {
         if( TMath::IsNaN(e[detJ])) continue;
         hcoinBGO->Fill(detID, detJ); 
      }
      
   }
   
   return kTRUE;
}
//############################################ TERMINATE
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

   printf("=============== loaded AutoFit.C, try showFitMethos()\n");
   gROOT->ProcessLine(".L armory/AutoFit.C");   
   printf("=============== Analyzer Utility\n");
   gROOT->ProcessLine(".L armory/Analyzer_Utili.c");
   gROOT->ProcessLine("listDraws()");
   

   

}
