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
                       int peakRange=1200, 
                       int tailRange=400, 
                       bool isLogz = false,
                       TString gate = "", 
                       TString treeName = "tree"
                      ){
   
   printf("================ Graphic Cut Creator for GAGG ============== \n");
   
   TChain * chain = new TChain(treeName);
   chain->Add(dataList);
   
   chain->GetListOfFiles()->Print();
   
   TString varX, varY, tag;

   gStyle->SetOptStat("neiou");

   TCanvas * cCutCreator = new TCanvas("cCutCreator", "GAGG PID Cut Creator", 100, 100, 800, 800);
   if( !cCutCreator->GetShowToolBar() ) cCutCreator->ToggleToolBar();

   cCutCreator->Update();
   if( isLogz ) cCutCreator->cd()->SetLogz();

   TCutG * cut = NULL;
   TObjArray * cutList = new TObjArray();

   TString expression;

   TH2F * h[NGAGG];

   int count = 0;
   for (Int_t i = 0; i < NGAGG; i++) {


      varX.Form("gaggT"); varY.Form("gaggP");

      h[i] = new TH2F(Form("h%d", i), Form("GAGG-%d",i), 500, 0, tailRange, 500, 0, peakRange);

      expression.Form("%s:%s>>h%d", varY.Data(), varX.Data(),i);
      gate.Form("gaggID==%d", i);

      chain->Draw(expression, gate, "col");
      
      if( h[i]->Integral() < 1000 ) {
         h[i]->SetMarkerStyle(20);
         h[i]->SetMarkerSize(0.4);
         h[i]->Draw("");
      }

      printf("======== make a graphic cut on the plot (double click to stop), %d-th cut: ", i );

      cCutCreator->Modified(); cCutCreator->Update();

      gPad->WaitPrimitive();

      cut = (TCutG*) gROOT->FindObject("CUTG");
      
      if( cut == NULL ){
         printf(" stopped by user. no file saved or changed. \n");
         break;
      }

      TString name; name.Form("cut%d", i);
      cut->SetName(name);
      cut->SetVarX(varX.Data());
      cut->SetVarY(varY.Data());
      cut->SetTitle(tag);
      cut->SetLineColor(i+1);
      cutList->Add(cut);

      printf(" cut-%d \n", i);
              
      count ++;
   }

   TFile * cutFile = new TFile(saveFileName, "recreate");
   cutList->Write("cutList", TObject::kSingleKey);
   
   printf("====> saved %d cuts into %s\n", count, saveFileName.Data());
   gROOT->ProcessLine(".q");

}
