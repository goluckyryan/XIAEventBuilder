#define Analyzer_cxx

#include "Analyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TH1.h>
#include <TCutG.h>
#include <TCanvas.h>
#include <TMath.h>
#include <vector>
#include <stdio.h>

//############################################ User setting

int rawEnergyRange[2] = {0,  12000}; // in ch
int energyRange[3] = {1, 30, 800}; // keV {resol, min, max}

double BGO_threshold = 0; // in ch

int pidMaxRange[3] = {500, 300, 1800}; //nBin, tail, and peak

TString e_corr = "correction_e.dat";

TString cutFileName1 = "protonCut.root";
//TString cutFileName1 = "alphaCut.root";
//TString cutFileName1 = "LiCut.root";

bool save_ev2 = false;

//############################################ end of user setting


//############################################ histogram declaration
TH2F * heVID;
TH1F * he[NCRYSTAL];

//TH2F * hgg[NCRYSTAL][NCRYSTAL];

TH2F * hcoin;
TH2F * hcrystalBGO;

TH2F * hgg;


TH1F * hTDiff;



TH2F * hPID[NGAGG];
TH2F * hPID_A[NGAGG];
TH2F * hPID_B[NGAGG];
TH2F * hGAGGVID;
///----- after gamma gate
TH2F * hPID_A_g[NGAGG];


///----- after calibration and BGO veto
TH2F * heCalVID;
TH1F * heCal[NCRYSTAL]; 
TH2F * hcoinBGO;
TH2F * hcrystalBGO_G;
TH1F * hg[NCLOVER];

///----- after particle gate
TH1F * heCal_g[NCRYSTAL];
TH2F * heCalVID_g;

TH1F * hg_g[NCLOVER];



///============= cut
TCutG * cut1;

//############################################ BEGIN
void Analyzer::Begin(TTree * tree){

   TString option = GetOption();

   totnumEntry = tree->GetEntries();

   printf( "=========================================================================== \n");
   printf( "==========================  Analysis.C/h   ================================ \n");
   printf( "======  total Entry : %lld \n", totnumEntry);
   printf( "=========================================================================== \n");

   printf("======================== Histograms declaration\n");
   
   heVID      = new TH2F("heVID",                                              "e vs ID; det ID; e [ch]", NCRYSTAL, 0, NCRYSTAL, rawEnergyRange[1] - rawEnergyRange[0], rawEnergyRange[0], rawEnergyRange[1]);
   heCalVID   = new TH2F("heCalVID",   Form("eCal vs ID (BGO veto > %.1f); det ID; Energy [keV]", BGO_threshold), NCRYSTAL, 0, NCRYSTAL, (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
   heCalVID_g = new TH2F("heCalVID_g", Form("eCal vs ID (BGO veto > %.1f & particle); det ID; Energy [keV]", BGO_threshold), NCRYSTAL, 0, NCRYSTAL, (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
   
   hTDiff = new TH1F("hTDiff", "data time different within an event; tick [10 ns]", 110, 0, 110);
   
   heVID->SetNdivisions(-410, "X");
   heCalVID->SetNdivisions(-410, "X");
   
   for( int i = 0; i < NCRYSTAL; i ++){
      he[i]      = new TH1F(   Form("he%02d", i),                                  Form("e -%02d", i), rawEnergyRange[1] - rawEnergyRange[0], rawEnergyRange[0], rawEnergyRange[1]);
      heCal[i]   = new TH1F(Form("heCal%02d", i),   Form("eCal-%02d (BGO veto > %.1f); Energy [keV];  count / %d keV", i, BGO_threshold, energyRange[0]), (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
      heCal_g[i] = new TH1F(Form("heCal_g%02d", i), Form("eCal-%02d (BGO veto > %.1f & particle); Energy [keV];  count / %d keV", i, BGO_threshold, energyRange[0]), (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
   }
   
   for( int i = 0; i < NCLOVER; i++){
      hg[i]   = new TH1F(Form("hg%02d", i), Form("Clover-%02d (added-back)", i), (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
      hg_g[i] = new TH1F(Form("hg_g%02d", i), Form("Clover-%02d (added-back) particle", i), (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
   }
   
   hgg = new TH2F("hgg", "Gamma - Gamma", (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2], (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
   
   for( int i = 0; i < NGAGG; i++){
     hPID[i] = new TH2F(Form("hPID%02d", i), Form("PID-%2d; tail; peak ", i), pidMaxRange[0], -20, pidMaxRange[1], pidMaxRange[0], -50, pidMaxRange[2]);
     hPID_A[i] = new TH2F(Form("hPID_A%02d",i), Form("PID_A detID = %2d; tail; peak ", i) , pidMaxRange[0], -20, pidMaxRange[1], pidMaxRange[0], -50, pidMaxRange[2]);
     hPID_B[i] = new TH2F(Form("hPID_B%02d",i), Form("PID_B detID = %2d; tail; peak ", i) , pidMaxRange[0], -20, pidMaxRange[1], pidMaxRange[0], -50, pidMaxRange[2]);
     hPID_A_g[i] = new TH2F(Form("hPID_A_g%02d",i), Form("PID_A detID = %2d (gated); tail; peak ", i) , pidMaxRange[0], -20, pidMaxRange[1], pidMaxRange[0], -50, pidMaxRange[2]);
  }

  hGAGGVID = new TH2F("hGAGGVID", "GAGG V ID", 80, 0, 80, 400, -50, 2000);
   
   /**
   for( int i = 0; i < NCRYSTAL; i++){
      for( int j = i+1; j < NCRYSTAL; j++){
         //hgg[i][j] = new TH2F(Form("hgg%02d%02d", i, j), Form("e%02d vs e%02d; e%02d; e%02d", i, j, i, j), 
         //        (rawEnergyRange[1] - rawEnergyRange[0])/2, rawEnergyRange[0], rawEnergyRange[1], 
         //        (rawEnergyRange[1] - rawEnergyRange[0])/2, rawEnergyRange[0], rawEnergyRange[1]);
      }
   }*/
   
   hcoin = new TH2F("hcoin", "detector coin.; det ID; det ID", NCRYSTAL, 0, NCRYSTAL, NCRYSTAL, 0 , NCRYSTAL); 
   
   hcoinBGO = new TH2F("hcoinBGO", Form("detector coin. (BGO veto > %.1f); det ID; det ID", BGO_threshold), NCRYSTAL, 0, NCRYSTAL, NCRYSTAL, 0 , NCRYSTAL); 
   hcrystalBGO = new TH2F("hcrystalBGO", Form("crystal vs BGO ; det ID; BGO ID"), NCRYSTAL, 0, NCRYSTAL, NBGO, 0 , NBGO); 
   hcrystalBGO_G = new TH2F("hcrystalBGO_G", Form("crystal vs BGO (BGO veto); det ID; BGO ID"), NCRYSTAL, 0, NCRYSTAL, NBGO, 0 , NBGO); 
   
   printf("======================== Load parameters.\n");
   
   eCorr = LoadCorrectionParameters(e_corr); 
   
   
   if( cutFileName1 != ""){
      printf("======================== Load cuts.\n");
   
      TFile * cutFile1 = new TFile(cutFileName1);
      cut1 = (TCutG *) cutFile1->Get("CUTG");
      printf(" %s is loaded.\n", cut1->GetName());
   }
      
   saveEV2 = save_ev2;

}


void Analyzer::PID_calculation(int ID){

}


//############################################ PROCESS
Bool_t Analyzer::Process(Long64_t entry){

   ProcessedEntries++;
   
   /*********** Progress Bar ******************************************/ 
   if (ProcessedEntries>totnumEntry*Frac-1) {
      TString msg; msg.Form("%llu", totnumEntry/1000);
      int len = msg.Sizeof();
      printf(" %3.0f%% (%*llu/%llu k) processed in %6.1f sec | expect %6.1f sec\n",
               Frac*100, len, ProcessedEntries/1000,totnumEntry/1000,StpWatch.RealTime(), StpWatch.RealTime()/Frac);
      StpWatch.Start(kFALSE);
      Frac+=0.1;
   }

   b_energy->GetEntry(entry);
   b_time->GetEntry(entry);
   b_multi->GetEntry(entry);
   b_multiCry->GetEntry(entry);
   b_detID->GetEntry(entry);
   b_qdc->GetEntry(entry);
   b_pileup->GetEntry(entry);
   
   if( multi == 0 ) return kTRUE;
   
   for( int i = 0; i < NCRYSTAL; i++) eCal[i] = TMath::QuietNaN();
   for( int i = 0; i < NCLOVER; i++) gamma[i] = 0;

   ///printf("---------------------------- %d \n", multi);
   
   double bg[NGAGG][2]={0}, peak[NGAGG][2]={0}, tail[NGAGG][2] = {0};
   int count[NGAGG] = {0} ;

   ///=========== make the particle gate
   for( int i = 0; i < multi ; i ++){
      
      if( e_t[i] - e_t[0] > 20 ) continue;
      if( pileup[i] == 1 ) continue;
      if( e[i] < 100 ) continue;
      int id = detID[i];
      if( id < 200 || id >= 300 ) continue;
      id = id - 200;

      hGAGGVID->Fill(id, e[i]);
      
      // GAGG_A
      if( (0 <= id && id < 50) ) {

            bg[id][0] = (qdc[i][0] + qdc[i][1])/60.;
            peak[id][0] = qdc[i][3]/20. - bg[id][0];
            tail[id][0] = qdc[i][5]/55. - bg[id][0];
            
            hPID_A[id]->Fill( tail[id][0], peak[id][0]);
            count[id] ++;
      }

      // GAGG_B
      if( 50 <= id   ) {
            id = id - 50;
            
            bg[id][1] = (qdc[i][0] + qdc[i][1])/60.;
            peak[id][1] = qdc[i][3]/20. - bg[id][1];
            tail[id][1] = qdc[i][5]/55. - bg[id][1];
            
            hPID_B[id]->Fill( tail[id][1], peak[id][1]);
            count[id]++;
      }
      
   }


  ///#########################################################
  ///================ coincident gate between proton and gamma

  ///printf("======================\n");
  for ( int i = 0 ; i < NGAGG ; i++){   

    if( count[i] == 2 ){
      
      double tailAvg = (tail[i][0]+tail[i][1])/2.;
      double peakAvg = (peak[i][0]+peak[i][1])/2.;
       
      hPID[i]->Fill( tailAvg, peakAvg);
    }
  }
   
   
  ///=========== Looping data for the event   
  for( int i = 0; i < multi ; i ++){
      if( pileup[i] == 1 ) continue;
      int id = detID[i];
      
      ///printf("%d %f %llu\n", id, e[i], e_t[i]);

      //======== Fill raw data
      if( 0 <= id &&  id < NCRYSTAL ){  /// gamma data
         heVID->Fill( id, e[i]);
         he[id]->Fill(e[i]);
         
         for ( int j = i + 1; j < multi; j++){
            if( 100 <= detID[j] && detID[j] < 200 ) hcrystalBGO->Fill(id, detID[j]-100); /// crystal - BGO coincident 
            
            if( detID[j] < 100 ) hcoin->Fill(id, detID[j]); /// crystal-crystal coincident 
            
         }
      }
      
      if ( 100 < id && id < 200 ){ /// BGO data
         
      }
      
      if( 200 < id && id  < 300){ /// GAGG
         continue;
      }
      
      
      //======== TDiff veto
      //if( !(e_t[i] - e_t[0] < 20 || e_t[i] - e_t[0] > 35) ) continue;
      if( e_t[i] - e_t[0] > 20 ) continue;
      
      if ( i > 0 ) hTDiff->Fill( e_t[i] - e_t[0]);
      
     
      //======== BGO veto
      bool dropflag = false;
      if( id < NCRYSTAL && multi > 1) {
         for( int j =  i + 1; j < multi; j++){
            if( detID[j] >= 100 && (detID[j]-100)*4 <= id && id < (detID[j]-100 +1)*4) {
               dropflag = true;
               break;
            }
         }
      }
      if( dropflag ) continue;
      
      if( 0<= id && id < NCRYSTAL ) {
         if( e_corr == "" ){
            eCal[id] = e[i];
         }else{
            ///========= apply energy correction
            int order = (int) eCorr[id].size();
            eCal[id] = 0;
            for( int k = 0; k < order ; k++){
               eCal[id] += eCorr[id][k] * TMath::Power(e[i], k);
            }
         }
         
         if( id != 30 ) heCalVID->Fill( id,  eCal[id]);
         heCal[id]->Fill(eCal[id]);
      
         for ( int j = i + 1; j < multi; j++){
            if( 100 <= detID[j] && detID[j] < 200 ) hcrystalBGO_G->Fill(id, detID[j]-100); /// crystal - BGO coincident 
         }
         
         
         ///========== add back and remove cross talk
         int cloverID = id /4;                  
         
         if( eCal[id] > energyRange[1]/4. ) gamma[cloverID] += eCal[id];
         
      }
      
   }
   
   for( int i = 0; i < NCLOVER; i++){
    for( int j = 0; j < NCLOVER; j++){
      if( gamma[i] > 0 && gamma[j] > 0 && i != j ) hgg->Fill( gamma[i], gamma[j]);
    }
   }
   
   for( int i = 0 ; i < NCLOVER; i++){
      if( gamma[i] > 0  ) {
         
         hg[i]->Fill(gamma[i]);
         
         for( int gi = 0; gi < NGAGG ; gi ++){    
           if( cut1->IsInside(tail[gi][0], peak[gi][0]) )  {
            hg_g[i]->Fill(gamma[i]);  
           }    
         }
      }

      if( abs(gamma[i] - 1052 ) < 8 ){
        for( int i = 0; i < NGAGG ; i ++){        
          hPID_A_g[i]->Fill( tail[i][0], peak[i][0]);
        }
      }
   }
   
   
   if ( saveEV2) Save2ev2();
   
   
   return kTRUE;
}
//############################################ TERMINATE
void Analyzer::Terminate(){
   
  if(saveEV2) fclose(outEV2);

  printf("============================== finishing.\n");
  gROOT->cd();

  int canvasXY[2] = {1600 , 800} ;// x, y
  int canvasDiv[2] = {4,2};
  TCanvas *cCanvas  = new TCanvas("cCanvas", "" ,canvasXY[0],canvasXY[1]);
  cCanvas->Modified(); cCanvas->Update();
  cCanvas->cd(); cCanvas->Divide(canvasDiv[0],canvasDiv[1]);

  gStyle->SetOptStat("neiou");

  int padID = 0;

  //========================= canvas 1
  padID++;   
  cCanvas->cd(padID);
  cCanvas->cd(padID)->SetLogz(0);
  heCalVID->Draw("colz");

  //========================= canvas 2   
  padID++;
  cCanvas->cd(padID);
  cCanvas->cd(padID)->SetLogz(1);
  hGAGGVID->Draw("colz");

  //========================= canvas 3   
  padID++;
  cCanvas->cd(padID);
  cCanvas->cd(padID)->SetLogz(0);
  hTDiff->Draw();

  //========================= canvas 3   
  padID++;
  cCanvas->cd(padID);
  cCanvas->cd(padID)->SetLogz(0);
  hg[6]->Draw();

  //========================= canvas 5
  padID++;
  cCanvas->cd(padID);
  cCanvas->cd(padID)->SetLogz(1);
  //heCalVID->Draw("colz");

  hPID[9]->Draw("colz");

  //========================= canvas 6   
  padID++;
  cCanvas->cd(padID);
  cCanvas->cd(padID)->SetLogz(1);
  hPID_A[9]->Draw("colz");
  cut1->Draw("same");

  //========================= canvas 7   
  padID++;
  cCanvas->cd(padID);
  cCanvas->cd(padID)->SetLogz(1);
  hPID_B[9]->Draw("colz");

  //========================= canvas 8   
  padID++;
  cCanvas->cd(padID);
  cCanvas->cd(padID)->SetLogz(1);
  //hPID_A_g[9]->Draw("colz");
  hg_g[6]->Draw();

/*
   //========================= canvas 1   
   cCanvas->cd(4);
   cCanvas->cd(4)->SetLogy(1);
   //gROOT->ProcessLine(".x script.C");   
   //hcrystalBGO_G->Draw("colz");
   hg[0]->SetLineColor(2);
   hg[0]->Draw();
   hg_g[0]->Draw("same");
   //hcoinBGO->Draw("colz");
   
   TCanvas *cAux  = new TCanvas("cAux", "" ,1000, 0, 800,600);
   cAux->cd();
   
   TH1F * h0 = (TH1F *) heCal[0]->Clone("h0");
   h0->Add(heCal[1]);
   h0->Add(heCal[2]);
   h0->Add(heCal[3]);
   
   TH1F * h0_g = (TH1F *) heCal_g[0]->Clone("h0_g");
   h0_g->Add(heCal_g[1]);
   h0_g->Add(heCal_g[2]);
   h0_g->Add(heCal_g[3]);
   
   h0->SetLineColor(kGreen+3);
   h0->Draw();
   h0_g->Draw("same");
   hg[0]->Draw("same");
   
/**/
   printf("=============== loaded AutoFit.C, try showFitMethos()\n");
   gROOT->ProcessLine(".L armory/AutoFit.C");   
   printf("=============== Analyzer Utility\n");
   gROOT->ProcessLine(".L armory/Analyzer_Utili.c");
   gROOT->ProcessLine("listDraws()");
   

   

}
