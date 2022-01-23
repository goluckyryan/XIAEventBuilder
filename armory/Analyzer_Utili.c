void listDraws(void) {
  printf("------------------- List of Plots -------------------\n");
  printf("  newCanvas() - Create a new Canvas\n");
  printf("-----------------------------------------------------\n");
  printf("       eVID() - e vs ID\n");
  printf("      drawE() - e for all %d detectors\n", NCRYSTAL);
  //printf("     drawGG() - Gamma - Gamma Coincident for all %d detectors\n", NCRYSTAL);
  printf("-----------------------------------------------------\n");
  printf(" energyCalibration() - Calibrate energy \n");
  printf("-----------------------------------------------------\n");
}

int nCanvas=0;
void newCanvas(int sizeX = 800, int sizeY = 600, int posX = 0, int posY = 0){
  TString name; name.Form("cNewCanvas%d", nCanvas);
  TCanvas * cNewCanvas = new TCanvas(name, name, posX, posY, sizeX, sizeY);
  nCanvas++;
  cNewCanvas->cd();
}

void eVID(bool cal = false, bool logz = false){
  TCanvas * cRawID = (TCanvas *) gROOT->FindObjectAny("cRawID");
  if( cRawID == NULL ) cRawID = new TCanvas("cRawID", "raw ID", 1000, 800);
  
  if( cRawID->GetShowEventStatus() == 0 ) cRawID->ToggleEventStatus();
  
  cRawID->cd(1)->SetGrid();
  if( logz ) cRawID->cd(1)->SetLogz();
  cal ? heCalVID->Draw("colz") : heVID->Draw("colz");
}

void drawE(int CloverID = -1, bool cali = false, bool isLogy = false, double xMin = 0, double xMax = 0){

   int nCrystalPerClover = 4;
   int nClover = NCRYSTAL / nCrystalPerClover;
   
   if( CloverID >= nClover ) {
      printf("Clover-ID > nClover = %d. \n", nClover);
      return;
   }
   
   int size = 300;

   TCanvas *cRawE = (TCanvas *) gROOT->FindObjectAny("cRawE");
   if( cRawE == NULL ) cRawE = new TCanvas("cRawE", cali ? "Cal e" : "Raw e", size * nClover, size * nCrystalPerClover);
   
   if( cRawE->GetShowEventStatus() == 0 ) cRawE->ToggleEventStatus();
   
   cRawE->Clear();
   if( CloverID >= 0 ) {
      nClover = 1;
      cRawE->Divide(nClover, 1);
   }else{
      cRawE->Divide(nClover, nCrystalPerClover, 0);
   }
   
   ///find max y
   double maxY = 0;
   int nDet = nClover*nCrystalPerClover;
   for( int i = (CloverID < 0 ? 0 : nCrystalPerClover*CloverID) ; i < (CloverID < 0 ? nDet : nCrystalPerClover*CloverID  + nDet)  ; i++){
      int mBin = cali ? heCal[i]->GetMaximumBin() : he[i]->GetMaximumBin();
      double max = cali ? heCal[i]->GetBinContent(mBin) : he[i]->GetBinContent(mBin);
      if( max > maxY ) maxY = max; 
   }
   maxY = maxY * 1.1;
   ///printf("max Y : %f \n", maxY);

   for (Int_t i = 0; i < nClover; i++) {
      
      int hID = nCrystalPerClover * CloverID + i ;
      if( cali ) {
         heCal[hID]->SetMaximum(maxY);
         heCal[hID]->Draw("same");
      }else{
         he[hID]->SetMaximum(maxY);
         he[hID]->Draw("same");
      }
      
      
      
      ///for( Int_t j = 0; j < nCrystalPerClover; j++){
      ///   int canvasID = CloverID < 0 ? nClover*j+ i + 1 : j + 1;
      ///   cRawE->cd(canvasID); 
      ///   cRawE->cd(canvasID)->SetGrid();       
      ///   cRawE->cd(canvasID)->SetTickx(2);   
      ///   cRawE->cd(canvasID)->SetTicky(2);   
      ///   cRawE->cd(canvasID)->SetBottomMargin(0.06);
      ///   if ( i == nClover -1 )  cRawE->cd(canvasID)->SetRightMargin(0.002);
      ///   if( isLogy ) cRawE->cd(canvasID)->SetLogy();
      ///   int hID = CloverID < 0 ? nCrystalPerClover*i+ j : nCrystalPerClover * CloverID + j ;
      ///   if( cali ) {
      ///      if ( xMin != 0 || xMax != 0 ) heCal[hID]->GetXaxis()->SetRangeUser(xMin, xMax);
      ///      heCal[hID]->SetMaximum(maxY);
      ///      heCal[hID]->Draw("");
      ///   }else{
      ///      if ( xMin != 0 || xMax != 0 ) he[hID]->GetXaxis()->SetRangeUser(xMin, xMax);
      ///      he[hID]->SetMaximum(maxY);
      ///      he[hID]->Draw("");
      ///   }
      ///}
   }
   
   cRawE->SetCrosshair(1);
   
}

/**
void drawGG(){

   int nCrystal = 4;
   int numCol = NCRYSTAL / nCrystal;
   
   int size = 300;

   TCanvas *cGG = (TCanvas *) gROOT->FindObjectAny("cGG");
   if( cGG == NULL ) cGG = new TCanvas("cGG", "Gamma - Gamma Coin.", size * NCRYSTAL, size * NCRYSTAL);
   cGG->Clear();cGG->Divide(NCRYSTAL, NCRYSTAL);
   
   for( int i = 0; i < NCRYSTAL; i ++){
      for( int j = i+1; j < NCRYSTAL; j ++){
         cGG->cd( NCRYSTAL * i + j +1 );
         hgg[i][j]->Draw("colz");
      }
   }
   
}
*/

void energyCalibration(int detID = -1, int BG = 10, double threshold = 0.1, double sigmaMax = 5, int peakDensity = 10){
   
   TCanvas *cCal = (TCanvas *) gROOT->FindObjectAny("cCal");
   if( cCal == NULL ) cCal = new TCanvas("cCal", "Energy Calibration", 1000, 0, 1000, 600);
   cCal->Clear();

   cCal->Divide(2,1);
   cCal->SetGrid();
   
   vector<double> refEnergy = {121.738, 
                               244.699, 
                               344.281,
                               411.115,
                               443.965,
                               778.903,
                               867.390,
                               964.055,
                              1085.842,
                              ///1089.700,
                              1112.087,
                              1408.022}; 
   
   double a0[NCRYSTAL];
   double a1[NCRYSTAL];
   
   for( int i = 0 ; i < NCRYSTAL; i++){
      if( detID >= 0 &&  i != detID ) continue;
      
      cCal->cd(1);
      he[i]->Draw();
      vector<double> peaks = fitAuto(he[i], BG, threshold, sigmaMax, peakDensity);
      vector<vector<double>> output = FindMatchingPair(peaks, refEnergy);
      
      vector<double> haha1 = output[0];
      vector<double> haha2 = output[1];
      
      TGraph * graph = new TGraph(haha1.size(), &haha1[0], &haha2[0] );
      cCal->cd(2);
      graph->Draw("A*");
      
      TF1 * fit = new TF1("fit", "pol1" );
      graph->Fit("fit", "");
      
      a0[i] = fit->GetParameter(0);
      a1[i] = fit->GetParameter(1);
      
      if( detID < 0 ) {
         printf("%2d | a0: %14.10f, a1: %14.10f\n", i, a0[i], a1[i]);
      }else{
         printf("%2d | a0, a1 = %14.10f\t%14.10f\n", i, a0[i], a1[i]);
      }
   }
   
   if( detID < 0 ){
      FILE * paraOut;
      TString filename;
      filename.Form("correction_e_auto.dat");
      paraOut = fopen (filename.Data(), "w+");
      printf("=========== save e-correction parameters to %s \n", filename.Data());
      for( int i = 0; i < NCRYSTAL; i++){
         fprintf(paraOut, "%14.10f\t%14.10f\n", a0[i], a1[i]);
      }
      fflush(paraOut);
      fclose(paraOut);
   }
   
}
