#define Monitor_raw_cxx

#include "Monitor_raw.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <vector>
#include <TStopwatch.h>

//############################################ User setting

int rawEnergyRange[2] = {500,  6000}; // in ch

TString e_corr = "correction_e.dat";

Long64_t maxEvent = -1;

//############################################ end of user setting

ULong64_t NumEntries = 0;
ULong64_t ProcessedEntries = 0;
Float_t Frac = 0.1; ///Progress bar
TStopwatch StpWatch;

vector<vector<double>> eCorr;

//############################################ histogram declaration

TH1F * hTDiff;
TH2F * hTDiffvEventID;

//############################################ BEGIN
void Monitor_raw::Begin(TTree * tree){
   TString option = GetOption();
   
   NumEntries = tree->GetEntries();
   
   printf("======================== Creating histograms\n");
   
   hTDiff = new TH1F("hTDiff", "time different between this and next event", 2000, -1000, 1000);
   hTDiffvEventID = new TH2F("hTDiffvEventID", "time different between this and next event; eventID; TDiff", 50, 0, maxEvent, 2000, -1000, 1000);
   
   printf("======================== end of histograms creation.\n");
   
   printf("======================== Load parameters.\n");
   eCorr = LoadCorrectionParameters(e_corr); 
   
   StpWatch.Start();
   printf("======================== Start processing....\n");
   
}

//############################################ PROCESS
Bool_t Monitor_raw::Process(Long64_t entry){

   ProcessedEntries++;
   
   /*********** Progress Bar ******************************************/ 
   if (ProcessedEntries>NumEntries*Frac-1) {
      TString msg; msg.Form("%llu", NumEntries/1000);
      int len = msg.Sizeof();
      //printf(" %3.0f%% (%*llu/%llu k) processed in %6.1f sec | expect %6.1f sec\033[A\n",
      printf(" %3.0f%% (%*llu/%llu k) processed in %6.1f sec | expect %6.1f sec\n",
               Frac*100, len, ProcessedEntries/1000,NumEntries/1000,StpWatch.RealTime(), StpWatch.RealTime()/Frac);
      StpWatch.Start(kFALSE);
      Frac+=0.1;
   }

   if ( maxEvent > 0 && entry > maxEvent ) {
      Abort(Form("user abort, > %lld \n", maxEvent));
      Terminate();
   }
   
   b_energy->GetEntry(entry);
   b_time_stamp->GetEntry(entry);
   
   ULong64_t t0 = t;
   
   if( entry < (Long64_t) NumEntries )  b_time_stamp->GetEntry(entry+1);;
   
   ULong64_t t1 = t;
  
   int tDiff = (int) t1 - t0;
  
   hTDiff->Fill( tDiff );
   hTDiffvEventID->Fill( entry, tDiff);

   return kTRUE;
}

//############################################ TERMINATE
void Monitor_raw::Terminate(){
   
   printf("============================== finishing.\n");
   gROOT->cd();
   
   TCanvas * cc = new TCanvas("cc", "cc", 2000, 1000);

   if( cc->GetShowEventStatus() == 0 ) cc->ToggleEventStatus();
 
   cc->Divide(2,1);
   
   cc->cd(1);
   hTDiff->Draw();
   
   cc->cd(2);
   cc->cd(2)->SetGrid();
   hTDiffvEventID->SetMarkerStyle(3);
   hTDiffvEventID->Draw("");
   
}
