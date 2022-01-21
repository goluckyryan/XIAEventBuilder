/**********************************************************/
/*                                                        */
/*  Modified  by Ryan From                                */
/*                                                        */
/* PXI SCAN CODE -- J.M. Allmond (ORNL) -- July 2016      */
/*                                                        */
/**********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TBenchmark.h"

#define RAND ((float) rand() / ((unsigned int) RAND_MAX + 1))   // random number in interval (0,1)

#define MAX_CRATES 2
#define MAX_BOARDS_PER_CRATE 13
#define MAX_CHANNELS_PER_BOARD 16
#define BOARD_START 2

#define MAX_ID MAX_CRATES*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD

#define HEADER_LENGTH 4     //unit = words with 4 bytes per word
#define MAX_SUB_LENGTH 2016 //unit = words with 4 bytes per word ; 2004 --> 40 micro-second trace + 4 word header

#define RAWE_REBIN_FACTOR 2.0 // Rebin 32k pixie16 spectra to something smaller to fit better into 8k.

#include "../mapping.h"  

#include "../armory/evtReader.h"

unsigned long long int dataCount=0;
unsigned long long int pileUpCount=0;
unsigned long long int evtCount=0;

///////////////////////////////////
// START OF MAIN FUNCTION        //
///////////////////////////////////
int main(int argc, char **argv) {

  printf("=====================================\n");
  printf("===      evt.to --> root          ===\n");
  printf("=====================================\n");  

  // Check that the corrent number of arguments were provided.
  if (argc != 2 && argc != 3  && argc != 4 && argc != 5)    {
    printf("Incorrect number of arguments:\n");
    printf("%s [*.to File]  <timeWindow>  <fraction>  <saveFile>\n", argv[0]);
    printf("         timeWindow : number of tick, 1 tick = 10 ns. default = 100 \n");   
    printf("           fraction : 0 to 100, default 100, last fraction of evt.to to root  \n");   
    printf("                                             e.g. 10, last 10%% of the evt.to to root \n");   
    printf("           saveFile : default is replace evt.to with root  \n");   
    return 1;
  }
  
  //CERN ROOT things
  TString inFileName = argv[1];
  int timeWindow = 100; 
  if( argc >= 3 ) timeWindow = atoi(argv[2]);
  
  int frac = 100;
  if( argc >= 4 ) frac = abs(atoi(argv[3]));

  TString outFileName = inFileName;  
  outFileName.Remove(inFileName.First('.'));
  outFileName.Append(".root");
  int pos = inFileName.Last('/');
  outFileName.Remove(0, pos+1);
  if( argc >= 5) outFileName = argv[4];
  
  evtReader * evt = new evtReader();
  evt->OpenFile(inFileName);
  DataBlock * data = evt->data;
  
  printf("  in file : %s \n", inFileName.Data());
  printf(" our file : \033[1;31m%s\033[m\n", outFileName.Data());

  printf(" max number of detector channal: %d \n", MAX_ID);
  
  printf(" Skipping the frist \033[0;34m %d %% \033[m of data \n", 100 - frac);
  
  printf("------------------------ Event building time window : %d tics = %d nsec \n", timeWindow, timeWindow*10);

  
  TFile * outRootFile = new TFile(outFileName, "recreate");
  outRootFile->cd();
  TTree * tree = new TTree("tree", "tree");

  unsigned long long          evID = 0;
  int                        multi = 0;
  int                   id[MAX_ID] = {0};
  double                 e[MAX_ID] = {TMath::QuietNaN()};  
  unsigned long long   e_t[MAX_ID] = {0};
  int               qdc[MAX_ID][8] = {0};
  Int_t                   multiCry = 0 ; /// this is total multiplicity for all crystal

  //unsigned short  pileup[MAXMULTI];

  tree->Branch("evID",         &evID, "event_ID/l"); 
  tree->Branch("multi",       &multi, "multi/I"); 
  tree->Branch("detID",           id, "detID[multi]/I");
  tree->Branch("e",                e, "e[multi]/D");
  tree->Branch("e_t",            e_t, "e_timestamp[multi]/l");
  tree->Branch("qdc",            qdc, "qdc[multi][8]/I");
  tree->Branch("multiCry", &multiCry, "multiplicity_crystal/I");  
  
  
  int countGP = 0; //gamma-particle coincident

  /////////////////////
  // MAIN WHILE LOOP //
  /////////////////////
  while ( evt->IsEndOfFile() == false ) { //main while loop 
      
      long long int etime = -1; 
      long long int tdif = -1; 
      
      while (1) { //get subevents and event build for one "event" 
        
        if ( evt->ReadBlock() == -1) break;
        if ( evt->GetFilePos() < evt->GetFileSize() * ( 1. - frac/100.) ) continue; 
        
        //Set reference time for event building
        if (etime == -1) {
            etime = data->time;
            tdif  = 0;
            multi = 0;
        }else {
            tdif = data->time - etime;
        }    
      
        //Check for end of event, rewind, and break out of while loop
        if (tdif > timeWindow) {
            
            etime = data->time;
            tdif  = 0;
            multi = 0;
            multiCry = 0;
            
            id[multi]  = data->detID;
            e[multi]   = data->energy;
            e_t[multi] = data->time;
            for( int i = 0; i < 8; i++) qdc[multi][i] = data->QDCsum[i];
            multi++ ;
            if( data->detID < 100 ) multiCry ++;
            
            break;           
        }else{
          //if within time window, fill array;
          id[multi]  = data->detID;
          e[multi]   = data->energy;
          e_t[multi] = data->time;
          for( int i = 0; i < 8; i++) qdc[multi][i] = data->QDCsum[i];
          multi++ ;
          if( data->detID < 100 ) multiCry ++;
        }    
               
        // total pileups
        if (data->pileup == 1) {
            pileUpCount++;
        }
        
     
      } //end while loop for unpacking sub events and event building for one "event"
      if (multi==0) break; //end main WHILE LOOP when out of events 
      dataCount = evt->GetBlockID(); 
      evID ++;
     
      evt->PrintStatus(10000);
   
      // when no gagg, don't save
      if( multiCry == 0 ||  multi == multiCry )  continue;
      
      countGP ++;
   
      outRootFile->cd();
      tree->Fill();
          
  } // end main while loop 
  /////////////////////////
  // END MAIN WHILE LOOP //
  /////////////////////////
  evt->PrintStatus(1);

  outRootFile->cd();
  tree->Write();
  outRootFile->Close();

  printf("\n\n\n==================== finished.\r\n");
  
  printf(" number of Gamma - GAGG coincdient : %d\n", countGP);

  return 0;
}


