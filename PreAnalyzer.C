#define PreAnalyzer_cxx

#include <TStyle.h>
#include <TMath.h>
#include <stdio.h>


#define TREESTRUCT  1 // if 0 = multi, 1 = array.

#include "PreAnalyzer.h"

//############################################ BEGIN
void PreAnalyzer::Begin(TTree * tree){

   TString option = GetOption();

   totnumEntry = tree->GetEntries();

   printf( "=========================================================================== \n");
   printf( "==========================  PreAnalysis.C/h   ============================= \n");
   printf( "======  total Entry : %lld \n", totnumEntry);
   printf( "=========================================================================== \n");


}



//############################################ PROCESS
Bool_t PreAnalyzer::Process(Long64_t entry){

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

   b_event_ID->GetEntry(entry);
   b_energy->GetEntry(entry);
   b_time->GetEntry(entry);
   b_multi->GetEntry(entry);
   b_multiCry->GetEntry(entry);
   b_detID->GetEntry(entry);
   b_qdc->GetEntry(entry);
   b_pileup->GetEntry(entry);
   b_runID->GetEntry(entry);
   
   multi_N = 0;
   multiGagg_N = 0;
   
   eventID = evID;
   runID_N = runID;
   
   for( int i = 0; i < NCLOVER; i++) {
      gammaID[i] = -1;
      gamma_N[i] = TMath::QuietNaN();
      gamma_t[i] = 0;
   }
   
   for( int i = 0; i < NGAGG; i++){
      gaggID[i] = -1;
      gagg_peak[i] = TMath::QuietNaN();
      gagg_tail[i] = TMath::QuietNaN();
      gagg_t[i] = 0;
   }
     
   double bg[NGAGG][2]={TMath::QuietNaN()}, peak[NGAGG][2]={TMath::QuietNaN()}, tail[NGAGG][2] = {TMath::QuietNaN()};
   ULong64_t gaggTime[NGAGG] = {0};
   ULong64_t gammaTime[NCLOVER] = {0};
   int count[NGAGG] = {0} ;   
   
  ///################## Gamma data from Clover
  for( int i = 0; i < NCRYSTAL; i++) eCal[i] = TMath::QuietNaN();
  for( int i = 0; i < NCLOVER; i++) gamma[i] = 0;

  for( int i = 0; i < multi ; i ++){
      if( pileup[i] == 1 ) continue;
      int id = detID[i];
         
      // GAGG_A
      if( (200 <= id && id < 250) ) {
            int id1 = id - 200;

            bg[id1][0] = (qdc[i][0] + qdc[i][1])/60.;
            peak[id1][0] = qdc[i][3]/20. - bg[id1][0];
            tail[id1][0] = qdc[i][5]/55. - bg[id1][0];
            
            if( gaggTime[id1] == 0 || e_t[i] < gaggTime[id1] ) gaggTime[id1] = e_t[i];

            count[id1] ++;
      }

      // GAGG_B
      if( 250 <= id && id  < 300  ) {
            int id2 = id - 250;
           
            bg[id2][1] = (qdc[i][0] + qdc[i][1])/60.;
            peak[id2][1] = qdc[i][3]/20. - bg[id2][1];
            tail[id2][1] = qdc[i][5]/55. - bg[id2][1];
            
            if( gaggTime[id2] == 0 || e_t[i] < gaggTime[id2] ) gaggTime[id2] = e_t[i];

            
            count[id2]++;
      }
      
      if( id >= 200 ) continue;
     
      //======== BGO veto
      bool dropflag = false;
      if( id < NCRYSTAL && multi > 1) {
         for( int j =  i + 1; j < multi; j++){
            if( 200 > detID[j] && detID[j] >= 100 && (detID[j]-100)*4 <= id && id < (detID[j]-100 +1)*4) {
               dropflag = true;
               break;
            }
         }
      }
      if( dropflag ) continue;
      
      if( 0<= id && id < NCRYSTAL ) {
         if( eCorr.size() == 0 ){
            eCal[id] = e[i];
         }else{
            ///========= apply energy correction
            int order = (int) eCorr[id].size();
            eCal[id] = 0;
            for( int k = 0; k < order ; k++){
               eCal[id] += eCorr[id][k] * TMath::Power(e[i], k);
            }
         }

         ///========== add back 
         int cloverID = id /4;                  
         if( eCal[id] > 10. ) {
            gamma[cloverID] += eCal[id];
            if( gammaTime[cloverID] == 0 || e_t[i] < gammaTime[cloverID] ) gammaTime[cloverID] = e_t[i];
            
         }
         ///========= remove cross talk
         
         ///========= doppler correction
         
      }

   }
   
   
   //################ Gamma-Paritcle 
   for( int i = 0 ; i < NCLOVER; i++){
      if( gamma[i] > 0  ) {
         
         if( TREESTRUCT == 0 ){
            gammaID[multi_N] = i;
            gamma_N[multi_N] = gamma[i];
            gamma_t[multi_N] = gammaTime[i];
         }else{
            gamma_N[i] = gamma[i];
            gamma_t[i] = gammaTime[i];
         }
         multi_N++; 
      }
   }
   
  for ( int i = 0 ; i < NGAGG ; i++){   
    if( count[i] == 2 ){
      
      if ( TREESTRUCT == 0 ){
         gaggID[multiGagg_N] = i;
         gagg_tail[multiGagg_N] = (tail[i][0]+tail[i][1])/2.;
         gagg_peak[multiGagg_N] = (peak[i][0]+peak[i][1])/2.;
         gagg_t[multiGagg_N]    = gaggTime[i];
      }else{
         gagg_tail[i] = (tail[i][0]+tail[i][1])/2.;
         gagg_peak[i] = (peak[i][0]+peak[i][1])/2.;
         gagg_t[i]    = gaggTime[i];
      }
      multiGagg_N++;
    }
  }
   
   if( multi_N == 0 && multiGagg_N == 0 ) return kTRUE;
   
   //############################ save
   saveFile->cd(); ///set focus on this file
   newTree->Fill();
   
   
   return kTRUE;
}
//############################################ TERMINATE
void PreAnalyzer::Terminate(){

   printf("============================== finishing.\n");

   saveFile->cd(); //set focus on this file
   newTree->Write(); 
   Long64_t nEntries = newTree->GetEntries();
   saveFile->Close();

   printf("-------------- done, saved in %s. number of entry : %lld\n", saveFileName.Data(), nEntries);

   gROOT->ProcessLine(".q");

}
