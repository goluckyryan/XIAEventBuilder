#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TClonesArray.h>
#include <TBenchmark.h>
#include <TMath.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLine.h>
#include <TBox.h>
#include <iostream>

void readRawTrace(TString fileName, int minDetID = 0, int maxDetID = 1000){

/**///==============================================================   

   TFile * f1 = new TFile (fileName, "read");
   TTree * tree = (TTree *) f1->Get("tree");
   
   if( tree == NULL ) {
      printf("===== Are you using gen_runXXX.root ? please use runXXX.root\n");
      return;
   }
   
   int totnumEntry = tree->GetEntries();
   printf( "========== total Entry : %d \n", totnumEntry);
   
   TCanvas * cRead = new TCanvas("cRead", "Read Trace", 0, 1500, 800, 300);
   cRead->Divide(1,1);
   for( int i = 1; i <= 2 ; i++){
      cRead->cd(i)->SetGrid();
   }
   cRead->SetGrid();
   
   gStyle->SetOptFit(0);
   
/**///==============================================================   
   UShort_t   detID;
   UShort_t trace[16000];
   UShort_t  traceLength;
   Int_t      QDC[8];
   ULong64_t  time;
   tree->SetBranchAddress("detID", &detID);
   tree->SetBranchAddress("trace", trace);
   tree->SetBranchAddress("trace_length", &traceLength);
   tree->SetBranchAddress("qdc", QDC); 
   tree->SetBranchAddress("e_t", &time); 
   
   TLatex text ;
   text.SetNDC();
   text.SetTextFont(62);
   text.SetTextSize(0.05);
   text.SetTextColor(2);
   
   bool breakFlag = false;  
   bool lastEvFlag = false; 
   int oldEv  = 0;

   TGraph * g = new TGraph();

   int QDCsum[8] = {0}; ///directly sum from trace
   double x[9] = {0, 31, 60, 75, 95, 105, 160, 175, 200};
   Color_t color[8] = {4, kGreen+2, 4, kGreen+2, 4, kGreen+2, 4, kGreen+2};

   
   for( int ev = 0; ev < totnumEntry; ev++){
      
      if( lastEvFlag ) {         
         ev = oldEv;
         lastEvFlag = false;
      }
      tree->GetEntry(ev);
      
      if( !(minDetID <=  detID &&  detID <= maxDetID ) ) continue;

      printf("-------------------------------- ev : %d \n", ev);

      printf("id : %d, trace Length : %u, time: %16llu ( enter = next , q = stop, w = last)\n", detID, traceLength, time);
      
      g->Clear();
      g->Set(traceLength);

      int minY = 9999, maxY = 0;
      for( int i = 0 ; i < 8; i++) QDCsum[i] = 0;
      
      for( int k = 0; k < traceLength ; k++) {
         g->SetPoint(k, k, trace[k]);
         for( int i = 0; i < 8; i++) if( x[i] <= k && k < x[i+1] ) QDCsum[i] += trace[k];
         if( trace[k] > maxY ) maxY = trace[k];
         if( trace[k] < minY ) minY = trace[k];
      }
      
      g->SetTitle(Form("ev: %d, id : %d, trace Length : %u, time: %16llu\n", ev, detID, traceLength, time));
      
      cRead->cd(1);
      cRead->Clear();
      g->Draw("AL");

      /*
      text.SetTextColor(kRed);
      text.DrawLatex(0.12, 0.85,  Form("QDC, QDCtrace (diff)"));
      for( int i = 0; i < 8; i++) {
         text.SetTextColor(color[i]);
         text.DrawLatex(0.12, 0.8-0.05*i,  Form("%d, %d (%+-d)", QDC[i], QDCsum[i], QDC[i]-QDCsum[i]));

         TBox * box = new TBox(x[i], minY, x[i+1], maxY);
         box->SetFillColorAlpha(color[i], 0.1);
         box->Draw("same");
      } 
      */ 
      cRead->Update();         
      gSystem->ProcessEvents();
      
      
      char s[80];
      fgets(s, sizeof s, stdin); 

      if( s[0] == 113 ) { // 'q' = 113
         breakFlag = true;
         break;
      }else if( s[0] == 119 ) { // 'w' = 119
         
         if(  ev > 0  ) {
           lastEvFlag = true;
           oldEv = ev -1;
         }else{
           printf(" the first event!!! \n");
         }
      }
      
      if( breakFlag ) break;
   
   }
   
   //gROOT->ProcessLine(".q");

}
