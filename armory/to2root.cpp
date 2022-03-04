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
  if (argc <= 3)    {
    printf("Incorrect number of arguments:\n");
    printf("%s [outFile] [timeWindow] [to File1]  [to File2] .... \n", argv[0]);
    printf("            outFile : output root file name\n");
    printf("         timeWindow : number of tick, 1 tick = 10 ns. default = 100 \n");   
    return 1;
  }
  
  TString outFileName = argv[1];
  int timeWindow = atoi(argv[2]);
  int nFile = argc - 3;
  TString inFileName[nFile];
  for( int i = 0 ; i < nFile ; i++){
    inFileName[i] = argv[i+3];
  }
  
  printf("====================================\n");
  
  evtReader * evt = new evtReader();
  DataBlock * data = evt->data;
  
  printf(" Number of input file : %d \n", nFile);
  printf(" out file : \033[1;31m%s\033[m\n", outFileName.Data());
  printf(" Event building time window : %d tics = %d nsec \n", timeWindow, timeWindow*10);


  TFile * outRootFile = new TFile(outFileName, "recreate");
  outRootFile->cd();
  TTree * tree = new TTree("tree", "tree");

  unsigned long long          evID = 0;
  int                        multi = 0;
  int                   id[MAX_ID] = {0};
  double                 e[MAX_ID] = {TMath::QuietNaN()};  
  unsigned long long   e_t[MAX_ID] = {0};
  int               qdc[MAX_ID][8] = {0};
  int                     multiCry = 0 ; /// this is total multiplicity for all crystal
  int                        runID = 0;  // date-run-fileID, Dec15-02-001 = 1502001
  int                    multiGagg = 0;
  //unsigned short  pileup[MAXMULTI];

  tree->Branch("evID",           &evID, "event_ID/l"); 
  tree->Branch("multi",         &multi, "multi/I"); 
  tree->Branch("detID",             id, "detID[multi]/I");
  tree->Branch("e",                  e, "e[multi]/D");
  tree->Branch("e_t",              e_t, "e_timestamp[multi]/l");
  tree->Branch("qdc",              qdc, "qdc[multi][8]/I");
  tree->Branch("multiCry",   &multiCry, "multiplicity_crystal/I");  
  tree->Branch("multiGagg", &multiGagg, "multiplicity_GAGG/I");  
  tree->Branch("runID",         &runID, "runID/I");  

  int countGP = 0; //gamma-particle coincident
  double totalDataSize = 0; 
  int outFileCount = 0;

  for( int i = 0; i < nFile; i++){
  
    evt->OpenFile(inFileName[i]);
    if( evt->IsOpen() == false ) continue;
    
    printf("==================================================== %d / %d\n", i+1, nFile);
    printf("\033[1;31m%s \033[m\n", inFileName[i].Data());
  
    int pos = inFileName[i].Last('/');
    TString temp = inFileName[i];
    temp.Remove(0, pos+1);
    temp.Remove(0, 3);
    temp.Remove(2, 1);
    temp.Remove(4, 1);
    temp.Remove(7);
    runID = atoi(temp.Data());
  
    long long int etime = -1; 
    long long int tdif = -1;


    while ( evt->IsEndOfFile() == false ) { //main while loop 
                  
        if ( evt->ReadBlock() == -1) break;
        
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
            
            //Gate
            //if( multiCry > 2 && multiGagg ==0  ) {
            //
            outRootFile->cd();
            tree->Fill();
            //
            //  countGP++;
            //}
            
            evID ++;
            
            etime = data->time;
            tdif  = 0;
            multi = 0;
            multiCry = 0;
            multiGagg = 0;

            int haha  = data->crate*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD + (data->slot-BOARD_START)*MAX_CHANNELS_PER_BOARD + data->ch;
            id[multi]  = mapping[haha];
            e[multi]   = data->energy;
            e_t[multi] = data->time;
            for( int i = 0; i < 8; i++) qdc[multi][i] = data->QDCsum[i];
            multi++ ;
            if( id[multi] < 100 ) multiCry ++;    
            if( id[multi] >= 200 ) multiGagg ++;            

        }else{
          //if within time window, fill array;
          int haha  = data->crate*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD + (data->slot-BOARD_START)*MAX_CHANNELS_PER_BOARD + data->ch;
          id[multi]  = mapping[haha];
          e[multi]   = data->energy;
          e_t[multi] = data->time;
          for( int i = 0; i < 8; i++) qdc[multi-1][i] = data->QDCsum[i];
          multi++ ;
          if( id[multi] < 100 ) multiCry ++;
          if( id[multi] >= 200 ) multiGagg ++;            

        }    
               
        // total pileups
        if (data->pileup == 1) {
            pileUpCount++;
        }
     
        evt->PrintStatus(10000);
     
    } // end main while loop 

    evt->PrintStatus(1);
    printf("\n\n\n");
    printf("             total number of event built : %llu\n", evID);    
    printf(" total number of Gamma - GAGG coincdient : %d (%.3f %%)\n", countGP, countGP*1.0/evID);    
    
    outRootFile->cd();
    tree->Write();
    
    totalDataSize += evt->GetFileSize()/1024./1024./1024.;
    
    double rootFileSize = outRootFile->GetSize()/1024./1024./1024. ; // in GB
    printf(" ----------- root file size : %.3f GB\n", rootFileSize);
    printf(" ---------- total read size : %.3f GB\n", totalDataSize);
    printf(" ----------- reduction rate : %.3f %%\n", rootFileSize*100./totalDataSize);
    
    if( rootFileSize > 3.0 ) {
      break;
    }
    
    ///try to open a new root file when file size > 2 GB
    /*if( rootFileSize > 2.0 ) {
      
      outRootFile->Close();
      delete outRootFile;
      delete tree;
      outFileCount += 1;
      
      if( outFileCount > 5 ) break;
      
      TString outFileName2 = outFileName;
      outFileName2.Insert(outFileName.Sizeof() - 6, Form("_%03d",outFileCount)); 
      outRootFile = new TFile( outFileName2, "recreate");
    
      tree = new TTree("tree", "tree");
      tree->Branch("evID",           &evID, "event_ID/l"); 
      tree->Branch("multi",         &multi, "multi/I"); 
      tree->Branch("detID",             id, "detID[multi]/I");
      tree->Branch("e",                  e, "e[multi]/D");
      tree->Branch("e_t",              e_t, "e_timestamp[multi]/l");
      tree->Branch("qdc",              qdc, "qdc[multi][8]/I");
      tree->Branch("multiCry",   &multiCry, "multiplicity_crystal/I");  
      tree->Branch("multiGagg", &multiGagg, "multiplicity_GAGG/I");  
      tree->Branch("runID",         &runID, "runID/I");  
    
    }*/
    
  }
  

  outRootFile->Close();

  printf("\n\n\n==================== finished.\r\n");
  
  printf(" number of Gamma - GAGG coincdient : %d\n", countGP);

  return 0;
}


