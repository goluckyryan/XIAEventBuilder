#define Monitor_raw_cxx


#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>

#include <TCanvas.h>
#include <TMath.h>
#include <vector>


//############################################ User setting

int rawEnergyRange[2] = {500,  6000}; // in ch, {min, max}
int    energyRange[3] = {1, 40,  2000}; // in keV, {resol, min, max}

TString e_corr = "correction_e.dat";

Long64_t maxEvent = -1;

bool isBuildEventRoot = true; // is true, no histogram


//############################################ end of user setting

#include "Monitor_raw.h"


//############################################ histogram declaration

TH1F * he[NCRYSTAL];
TH1F * heCal[NCRYSTAL];

TH2F * heCalvID;

//############################################ BEGIN
void Monitor_raw::Begin(TTree * tree){
   TString option = GetOption();

   totnumEntry = tree->GetEntries();

   printf( "=========================================================================== \n");
   printf( "========================== Monitor_raw.C/h ================================ \n");
   printf( "======  total Entry : %lld \n", totnumEntry);
   printf( "=========================================================================== \n");
   
   
   if( isBuildEvent == false && isBuildEventRoot == false) {
      printf("======================== Creating histograms\n");

      for( int i = 0 ; i < NCRYSTAL ; i++){
         he[i]    = new TH1F (Form("he%02d", i),            Form("e%02d", i),            rawEnergyRange[1]-rawEnergyRange[0], rawEnergyRange[0], rawEnergyRange[1]);  
         heCal[i] = new TH1F (Form("heCal%02d", i), Form("e%02d (Cali.)", i), (energyRange[2]-energyRange[1])/energyRange[0],    energyRange[1],    energyRange[2]);  
         
         switch (i%4){
            case 0: he[i]->SetLineColor(2);break;
            case 1: he[i]->SetLineColor(kYellow+3);break;
            case 2: he[i]->SetLineColor(kGreen+2);break;
            case 3: he[i]->SetLineColor(4);break;
         }
         switch (i%4){
            case 0: heCal[i]->SetLineColor(2);break;
            case 1: heCal[i]->SetLineColor(kYellow+3);break;
            case 2: heCal[i]->SetLineColor(kGreen+2);break;
            case 3: heCal[i]->SetLineColor(4);break;
         }
      }
      
      heCalvID = new TH2F("heCalvID", "ID vs Energy (Cali.); ID; Energy",   NCRYSTAL, 0, NCRYSTAL, (energyRange[2]-energyRange[1])/energyRange[0],    energyRange[1],    energyRange[2]);

      
      printf("------------------------ end of histograms creation.\n");
   }
   
   printf("======================== Load parameters.\n");
   eCorr = LoadCorrectionParameters(e_corr); 
   
   time0 = 0;
   timeDiff = 0;
   
   ClearTreeData();
   
   
}

//############################################ PROCESS
Bool_t Monitor_raw::Process(Long64_t entry){

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

   if ( maxEvent > 0 && entry > maxEvent ) {
      Abort(Form("user abort, > %lld \n", maxEvent));
      Terminate();
   }
   
   if( isBuildEventRoot || isBuildEvent ) entry = index[entry];
   
   b_ID->GetEntry(entry);
   b_energy->GetEntry(entry);
   b_energy_timestamp->GetEntry(entry);
   
   if ( isBuildEventRoot || isBuildEvent ) {
      BuildEvent();
   }else{
      if( detID < 100 ){
         he[detID]->Fill(energy);
         
         ///double eCal = ApplyCorrection(eCorr, detID, e); //slow, why?
         double eCal = eCorr[detID][0] + eCorr[detID][1]*energy;
         heCal[detID]->Fill(eCal);
         
         heCalvID->Fill(detID, eCal);
      }
   }

   return kTRUE;
}

//############################################ TERMINATE
void Monitor_raw::Terminate(){
   
   
   printf("============================== finished.\n");
   gROOT->cd();
   
   if ( isBuildEventRoot || isBuildEvent ) {
      saveFile->cd();
      newtree->Write();
      saveFile->Close();
   }else{
   
      int nCrystalPerClover = 4;
      int nClover = NCRYSTAL / nCrystalPerClover;

      TCanvas * cc = new TCanvas("cc", "cc", 2000, 2000);
      if( cc->GetShowEventStatus() == 0 ) cc->ToggleEventStatus();
      cc->Divide(1, 9, 0);
      
       for (Int_t i = 0; i < nClover; i++) {
         int canvasID =  i + 1;
         cc->cd(canvasID); 
         cc->cd(canvasID)->SetGrid();       
         cc->cd(canvasID)->SetTickx(2);   
         cc->cd(canvasID)->SetTicky(2);   
         cc->cd(canvasID)->SetBottomMargin(0.06);
         cc->cd(canvasID)->SetLogy();

         for( Int_t j = 0; j < nCrystalPerClover; j++){
            int hID = nCrystalPerClover*i+ j ;
            heCal[hID]->Draw("same");         
            //he[hID]->Draw("same");         
         }
      }
      cc->SetCrosshair(1);
      
      TCanvas * c1 = new TCanvas("c1", "c1", 1000, 1000);
      if( c1->GetShowEventStatus() == 0 ) c1->ToggleEventStatus();
      c1->SetLogz();
      c1->SetGridx();
      heCalvID->SetNdivisions(-409, "X");
      heCalvID->Draw("colz");
      
   }
   
}
