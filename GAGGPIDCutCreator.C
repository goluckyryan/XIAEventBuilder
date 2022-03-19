#include <TH2F.h>
#include <TFile.h>
#include <TChain.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TString.h>
#include <TObjArray.h>

#include "mapping.h"

void GAGGPIDCutCreator(TString dataList, 
                       TString saveFileName = "GAGGPIDCuts.root", 
                       int detID = -100, // -100 = new, -X = append from X, X = only X
                       int peakRange=1200, 
                       int tailRange=400, 
                       bool isLogz = false, 
                       float frac = 0.3
                      ){
   
   printf("================ Graphic Cut Creator for GAGG ============== \n");
   
   TChain * chain = new TChain("tree");
   chain->Add(dataList);
   
   chain->GetListOfFiles()->Print();
   
   Long64_t numEntries = chain->GetEntries();
   
   printf(" total number of entries : %llu, using the first %.0f%% data.\n", numEntries, frac*100.);
   
   TString varX, varY, tag;

   gStyle->SetOptStat("neiou");

   TCanvas * cCutCreator = new TCanvas("cCutCreator", "GAGG PID Cut Creator", 100, 100, 800, 800);
   if( !cCutCreator->GetShowToolBar() ) cCutCreator->ToggleToolBar();

   cCutCreator->Update();
   if( isLogz ) cCutCreator->cd()->SetLogz();

   TString expression, gate;
   TH2F * h[NGAGG];

   TFile * cutFile;
   TCutG * cut = NULL;
   TObjArray * cutList;

   
   /// load the cutFile and load the cutList
   

   int startID = 0;
   int stopID = NGAGG;

   if( detID != -100 ) {
      cutFile = new TFile(saveFileName, "UPDATE");
      bool listExist = cutFile->GetListOfKeys()->Contains("cutList");
      if( !listExist ) {
         cutList = new TObjArray();
      }else{
         cutList = (TObjArray*) cutFile->FindObjectAny("cutList");
         int numCut = cutList->GetLast()+1;
         printf("----- found %d cuts \n", numCut);
         for( int k = 0; k < numCut; k++){
            if( cutList->At(k) != NULL ){
               printf("found a cut at %2d \n", k);
            }else{
               printf("     No cut at %2d \n", k);
            }
         }
      }
      
      if( detID < 0 ) startID = abs(detID);
      if( detID >= 0 ){
         startID = detID;
         stopID = detID+1;
      }
      
   }else{
      cutFile = new TFile(saveFileName, "recreate");
      cutList = new TObjArray();
   }
   printf("=================== plotting histogram, start at GAGG-%0d \n", startID);
   
   int count = 0;

   for (Int_t i = startID; i < stopID; i++) {

      varX.Form("gaggT"); varY.Form("gaggP");

      h[i] = new TH2F(Form("h%d", i), Form("GAGG-%d",i), 500, 0, tailRange, 500, 0, peakRange);

      expression.Form("%s:%s>>h%d", varY.Data(), varX.Data(),i);
      gate.Form("gaggID==%d", i);

      chain->Draw(expression, gate, "col", numEntries * frac);
      
      cut = (TCutG*) cutList->At(i);
      cut->Draw("same");
      
      if( h[i]->Integral() < 1000 ) {
         h[i]->SetMarkerStyle(20);
         h[i]->SetMarkerSize(0.4);
         h[i]->Draw("");
      }

      printf("======== make a graphic cut on the plot (double click to stop), %2d-th cut: ", i );

      if( detID != -100 ) 
      
      cCutCreator->Modified(); cCutCreator->Update();

      gPad->WaitPrimitive();

      cut = (TCutG*) gROOT->FindObject("CUTG");
      
      if( cut == NULL ){
         printf(" stopped by user.\n");
         break;
      }

      TString name; name.Form("cut%d", i);
      cut->SetName(name);
      cut->SetVarX(varX.Data());
      cut->SetVarY(varY.Data());
      cut->SetTitle(tag);
      cut->SetLineColor(i+1);
      if( detID != -100 ){
         cutList->AddAt(cut, abs(detID));
      }else{
         cutList->Add(cut);
      }
      printf(" cut-%d \n", i);
              
      count ++;
   }

   cutList->Write("cutList", TObject::kSingleKey);
   if( detID == -100) {
      printf("====> saved %d cuts into %s\n", count, saveFileName.Data());
   }else{
      if( detID < 0 ) {
         printf("====> ReDo %d cuts into %s, from GAGG-%2d\n", count, saveFileName.Data(), startID);
      }else{
         printf("====> ReDo %d cuts into %s, GAGG-%2d\n", count, saveFileName.Data(), startID);
      }
   }
   gROOT->ProcessLine(".q");

}
