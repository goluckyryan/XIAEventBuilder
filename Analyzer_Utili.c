void listDraws(void) {
  printf("------------------- List of Plots -------------------\n");
  printf("      rawID() - Raw e vs ID\n");
  printf("      drawE() - Raw e for all %d detectors\n", NCLOVER);
  printf("-----------------------------------------------------\n");
}

void rawID(){
  TCanvas * cRawID = (TCanvas *) gROOT->FindObjectAny("cRawID");
  if( cRawID == NULL ) cRawID = new TCanvas("cRawID", "raw ID", 1000, 800);
  cRawID->cd(1)->SetGrid();
  heVID->Draw("colz");
}

void drawE(bool isLogy = false, bool cali = false){

   int numCol = NCLOVER / 4;

   TCanvas *cRawE = (TCanvas *) gROOT->FindObjectAny("cRawE");
   if( cRawE == NULL ) cRawE = new TCanvas("cRawE", cali ? "Cal e" : "Raw e", 1200, 800);
   cRawE->Clear();cRawE->Divide(numCol, 4);
   
   //cRawE->SetRightMargin(0);
   //cRawE->SetLeftMargin(0);
   //cRawE->SetTopMargin(0);
   //cRawE->SetBottomMargin(0);
   //cRawE->SetTicks(1,1);
   //cRawE->SetBorderMode(1);
   
   for (Int_t i = 0; i < 4; i++) {
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
         int hID = 4*j+ i;
         if( cali ) {
            heCal[hID]->Draw("");
         }else{
            he[hID]->Draw("");
         }
      }
   }

}
