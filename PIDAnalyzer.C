#define PIDAnalyzer_cxx


#include "PIDAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TH1.h>
#include <TCutG.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TObjArray.h>
#include <vector>
#include <stdio.h>

//############################################ User setting

int energyRange[3] = {1, 30, 800}; // keV {resol, min, max}

int pidMaxRange[3] = {500, 400, 1600}; //nBin, tail, and peak

TString cutFileName1 = "testCuts.root";
TString cutFileName2 = ""; //"alphaCut.root";
TString cutFileName3 = ""; //"tritonCut.root";

short timeGateFlag = 4; // 0 = off, 1 <, 2 >, 3 sandwish, 4 !sandwish
unsigned int timeGate[2] = {45, 65}; // if timeGateFlag < 3, only timeGate[0] use, else, {min, max}

//############################################ end of user setting

TH2F * hgg;

TH1F * hg[NCLOVER];
TH1F * hg_g1[NCLOVER];
TH1F * hg_g2[NCLOVER];
TH1F * hg_g3[NCLOVER];

TH2F * hPID[NGAGG];

///============= cut
TObjArray * cutList1;
TObjArray * cutList2;
TObjArray * cutList3;
TCutG * cut;

void PIDAnalyzer::Begin(TTree *tree){

   TString option = GetOption();
   
   totnumEntry = tree->GetEntries();

   printf( "=========================================================================== \n");
   printf( "==========================  Analysis.C/h   ================================ \n");
   printf( "======  total Entry : %lld \n", totnumEntry);
   printf( "=========================================================================== \n");

   printf("======================== Histograms declaration\n");
   
   
   for( int i = 0; i < NCLOVER; i++){
      hg[i]    = new TH1F(Form("hg%02d", i), Form("Clover-%02d (added-back)", i), (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
      hg_g1[i] = new TH1F(Form("hg_g1%02d", i), Form("Clover-%02d (added-back) particle", i), (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
      hg_g2[i] = new TH1F(Form("hg_g2%02d", i), Form("Clover-%02d (added-back) particle", i), (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
      hg_g3[i] = new TH1F(Form("hg_g3%02d", i), Form("Clover-%02d (added-back) particle", i), (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
   }
   
   hgg = new TH2F("hgg", "Gamma - Gamma", (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2], (energyRange[2] - energyRange[1])/energyRange[0], energyRange[1], energyRange[2]);
   
   
   for( int i = 0; i < NGAGG; i++){
     hPID[i] = new TH2F(Form("hPID%02d", i), Form("PID-%2d; tail; peak ", i), pidMaxRange[0], -20, pidMaxRange[1], pidMaxRange[0], -50, pidMaxRange[2]);
  }
  
  printf("======================== Load cuts.\n");   
   if( cutFileName1 != ""){
      TFile * cutFile1 = new TFile(cutFileName1);
         if( cutFile1->IsOpen()){
         cutList1 = (TObjArray *) cutFile1->FindObjectAny("cutList");
         int numCut1 = cutList1->GetEntries();
         printf("=========== found %d cutG in %s \n", numCut1, cutFile1->GetName());

         for(int i = 0; i < numCut1 ; i++){
            printf("cut name : %s , VarX: %s, VarY: %s, numPoints: %d \n",
               cutList1->At(i)->GetName(),
               ((TCutG*)cutList1->At(i))->GetVarX(),
               ((TCutG*)cutList1->At(i))->GetVarY(),
               ((TCutG*)cutList1->At(i))->GetN());
         }
      }else{
         cutList1 = NULL;
      }
   }

   if( cutFileName2 != ""){
      TFile * cutFile2 = new TFile(cutFileName2);
      if( cutFile2->IsOpen()){
         cutList2 = (TObjArray *) cutFile2->FindObjectAny("cutList");
         int numCut2 = cutList2->GetEntries();
         printf("=========== found %d cutG in %s \n", numCut2, cutFile2->GetName());

         for(int i = 0; 2 < numCut2 ; i++){
            printf("cut name : %s , VarX: %s, VarY: %s, numPoints: %d \n",
               cutList2->At(i)->GetName(),
               ((TCutG*)cutList2->At(i))->GetVarX(),
               ((TCutG*)cutList2->At(i))->GetVarY(),
               ((TCutG*)cutList2->At(i))->GetN());
         }
      }else{
         cutList2 = NULL;
      }
   }

   if( cutFileName3 != ""){
      TFile * cutFile3 = new TFile(cutFileName3);
      if( cutFile3->IsOpen()){
         cutList3 = (TObjArray *) cutFile3->FindObjectAny("cutList");
         int numCut3 = cutList3->GetEntries();
         printf("=========== found %d cutG in %s \n", numCut3, cutFile3->GetName());

         for(int i = 0; 2 < numCut3 ; i++){
            printf("cut name : %s , VarX: %s, VarY: %s, numPoints: %d \n",
               cutList3->At(i)->GetName(),
               ((TCutG*)cutList3->At(i))->GetVarX(),
               ((TCutG*)cutList3->At(i))->GetVarY(),
               ((TCutG*)cutList3->At(i))->GetN());
         }
      }else{
         cutList3 = NULL;
      }
   }

}

Bool_t PIDAnalyzer::Process(Long64_t entry){
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

   b_eventID->GetEntry(entry);
   b_runID->GetEntry(entry);
   b_multi->GetEntry(entry);
   b_multiGagg->GetEntry(entry);
   b_gammaID->GetEntry(entry);
   b_gamma->GetEntry(entry);
   b_gamma_t->GetEntry(entry);
   b_gaggID->GetEntry(entry);
   b_gaggP->GetEntry(entry);
   b_gaggT->GetEntry(entry);
   b_gagg_t->GetEntry(entry);

   ///################## Gagg
   bool fillFlag1 = false;
   bool fillFlag2 = false;
   bool fillFlag3 = false;

   for( int i = 0 ; i < multiGagg; i++){
      hPID[gaggID[i]]->Fill(gaggT[i],gaggP[i]);

      if( cutList1 != NULL && i < cutList1->GetEntries()) {
         cut = (TCutG *) cutList1->At(i);
         if( fillFlag1 == false && cut->IsInside(gaggT[i],gaggP[i]) )  {
            fillFlag1 = true; 
            break;
         }
      }
      if( cutList2 != NULL && i < cutList2->GetEntries()) {
         cut = (TCutG *) cutList2->At(i);
         if( fillFlag2 == false && cut->IsInside(gaggT[i],gaggP[i]) )  {
            fillFlag2 = true;  
            break;
         }
      }
      if( cutList3 != NULL && i < cutList3->GetEntries()) {
         cut = (TCutG *) cutList3->At(i);
         if( fillFlag2 == false && cut->IsInside(gaggT[i],gaggP[i]) )  {
            fillFlag3 = true; 
            break;
         }
      }
     
   } 
   
   
  ///################## Gamma data from Clover
  
  for( int i = 0; i < multi ; i ++){

      //======== TDiff veto
      //if( timeGateFlag == 1 && e_t[i] - e_t[0] > timeGate[0] ) continue;
      //if( timeGateFlag == 2 && e_t[i] - e_t[0] < timeGate[0] ) continue;
      //if( timeGateFlag == 3 && !(timeGate[0] < e_t[i] - e_t[0] && e_t[i] - e_t[0] < timeGate[1]) ) continue;
      //if( timeGateFlag == 4 &&   timeGate[0] < e_t[i] - e_t[0] && e_t[i] - e_t[0] < timeGate[1] ) continue;
      
      hg[gammaID[i]]->Fill(gamma[i]);

      if( fillFlag1 ) hg_g1[gammaID[i]]->Fill(gamma[i]);
      if( fillFlag2 ) hg_g2[gammaID[i]]->Fill(gamma[i]);
      if( fillFlag3 ) hg_g3[gammaID[i]]->Fill(gamma[i]);
      
      for( int j = 0 ; j < multi; j++){
         if( j != i ) hgg->Fill( gamma[i], gamma[j]);
      }
      
   }

   
   return kTRUE;
}



void PIDAnalyzer::Terminate(){

  printf("============================== finishing.\n");
  gROOT->cd();

  int canvasXY[2] = {1600 , 800} ;// x, y
  int canvasDiv[2] = {4,3};
  TCanvas *cCanvas  = new TCanvas("cCanvas", "" ,canvasXY[0],canvasXY[1]);
  if( !cCanvas->GetShowToolBar() ) cCanvas->ToggleToolBar();
  cCanvas->Modified(); cCanvas->Update();
  cCanvas->cd(); cCanvas->Divide(canvasDiv[0],canvasDiv[1]);

  gStyle->SetOptStat("neiou");

  int padID = 0;
  
   for( int i = 0; i < NCLOVER; i++){
      padID++;
      cCanvas->cd(padID);
      cCanvas->cd(padID)->SetLogz(0);
      hg[i]->Draw("");
      hg_g1[i]->SetLineColor(2);hg_g1[i]->Draw("same");
      hg_g2[i]->SetLineColor(4);hg_g2[i]->Draw("same");
      hg_g3[i]->SetLineColor(6);hg_g3[i]->Draw("same");
   }
  
   TCanvas *cPID  = new TCanvas("cPID", "" ,7*200,4*200);
   if( !cPID->GetShowToolBar() ) cPID->ToggleToolBar();
   cPID->Modified(); cPID->Update();
   cPID->cd(); cPID->Divide(7,4);
   
   for( int i = 0; i < NGAGG; i++){
      cPID->cd(i+1);
      hPID[i]->Draw("colz");
      if( cutList1 != NULL && i < cutList1->GetEntries()) ((TCutG*)cutList1->At(i))->Draw("same");
      if( cutList2 != NULL && i < cutList2->GetEntries()) ((TCutG*)cutList2->At(i))->Draw("same");
      if( cutList3 != NULL && i < cutList3->GetEntries()) ((TCutG*)cutList3->At(i))->Draw("same");
   }


   printf("=============== loaded AutoFit.C, try showFitMethos()\n");
   gROOT->ProcessLine(".L armory/AutoFit.C");   


}
