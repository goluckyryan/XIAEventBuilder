void listDraws(void) {
  printf("------------------- List of Plots -------------------\n");
  printf("  newCanvas() - Create a new Canvas\n");
  printf("-----------------------------------------------------\n");
  printf("    rawEvID() - Raw e vs ID\n");
  printf("      drawE() - Raw e for all %d detectors\n", NCLOVER);
  printf("     drawGG() - Gamma - Gamma Coincident for all %d detectors\n", NCLOVER);
  printf("-----------------------------------------------------\n");
}

int nCanvas=0;
void newCanvas(int sizeX = 800, int sizeY = 600, int posX = 0, int posY = 0){
  TString name; name.Form("cNewCanvas%d | %s", nCanvas, canvasTitle.Data());
  TCanvas * cNewCanvas = new TCanvas(name, name, posX, posY, sizeX, sizeY);
  nCanvas++;
  cNewCanvas->cd();
}

void rawEvID(){
  TCanvas * cRawID = (TCanvas *) gROOT->FindObjectAny("cRawID");
  if( cRawID == NULL ) cRawID = new TCanvas("cRawID", "raw ID", 1000, 800);
  cRawID->cd(1)->SetGrid();
  heVID->Draw("colz");
}

void drawE(bool isLogy = false, bool cali = false){

   int nCrystal = 4;
   int numCol = NCLOVER / nCrystal;
   
   int size = 300;

   TCanvas *cRawE = (TCanvas *) gROOT->FindObjectAny("cRawE");
   if( cRawE == NULL ) cRawE = new TCanvas("cRawE", cali ? "Cal e" : "Raw e", size * numCol, size * nCrystal);
   cRawE->Clear();cRawE->Divide(numCol, 4);
   
   //cRawE->SetRightMargin(0);
   //cRawE->SetLeftMargin(0);
   //cRawE->SetTopMargin(0);
   //cRawE->SetBottomMargin(0);
   //cRawE->SetTicks(1,1);
   //cRawE->SetBorderMode(1);
   
   for (Int_t i = 0; i < nCrystal; i++) {
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
         int hID = nCrystal*j+ i;
         if( cali ) {
            heCal[hID]->Draw("");
         }else{
            he[hID]->Draw("");
         }
      }
   }

}

void drawGG(){

   int nCrystal = 4;
   int numCol = NCLOVER / nCrystal;
   
   int size = 300;

   TCanvas *cGG = (TCanvas *) gROOT->FindObjectAny("cGG");
   if( cGG == NULL ) cGG = new TCanvas("cGG", "Gamma - Gamma Coin.", size * NCLOVER, size * NCLOVER);
   cGG->Clear();cGG->Divide(NCLOVER, NCLOVER);
   
   for( int i = 0; i < NCLOVER; i ++){
      for( int j = i+; j < NCLOVER; j ++){
         cGG->cd( NCLOVER * i + j +1 );
         hgg[i][j]->Draw("colz");
      }
   }
   
}
