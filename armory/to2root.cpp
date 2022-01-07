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

#define MAXMULTI 100 

#include "../mapping.h"  

/////////////////////
// RAW EVENT TYPES //
/////////////////////
struct measurement
{
    int chn; 
    int sln;
    int crn;
    int id;
    int hlen;
    int elen;
    int trlen;          //number of samples
    int trwlen;         //number of words (two samples per word)
    int fcode;          //pileup flag
    long long int time;
    int ctime;
    int ctimef;
    int e;
    int extra;
    short tr[4096];
    int esum[4];
    int qsum[8];        
}; 
struct measurement data = {0};

unsigned long long int dataCount=0;
unsigned long long int pileUpCount=0;
unsigned long long int evtCount=0;

///////////////////////////////////
// START OF MAIN FUNCTION        //
///////////////////////////////////
int main(int argc, char **argv) {
  
  float tempf=0;

  //temp buffer for each sub event
  unsigned int sub[MAX_SUB_LENGTH];
  memset(sub, 0, sizeof(sub));

  printf("=====================================\n");
  printf("===      evt.to --> root          ===\n");
  printf("=====================================\n");  

  // Check that the corrent number of arguments were provided.
  if (argc != 2 && argc != 3 )    {
    printf("Incorrect number of arguments:\n");
    printf("%s [*.to File] [timeWindow] \n", argv[0]);
    printf("         timeWindow : number of tick, 1 tick = 10 ns. default = 100 \n");   
    return 1;
  }
  
  //CERN ROOT things
  TString inFileName = argv[1];
  TString outFileName = inFileName;
  
  int timeWindow = 100; 
  if( argc >= 3 ) timeWindow = atoi(argv[2]);
  
  outFileName.Remove(inFileName.First('.'));
  outFileName.Append(".root");
  
  printf("  in file : %s \n", inFileName.Data());
  printf(" our file : %s \n", outFileName.Data());

  printf(" max number of detector channal: %d \n", MAX_ID);
  printf("------------------------ Event building time window : %d tics = %d nsec \n", timeWindow, timeWindow*10);

  TFile * outRootFile = new TFile(outFileName, "recreate");
  outRootFile->cd();
  TTree * tree = new TTree("tree", "tree");

  unsigned long long          evID = 0;
  int                        multi = 0;
  int                 id[MAXMULTI] = {0};
  double               e[MAXMULTI] = {TMath::QuietNaN()};  
  unsigned long long e_t[MAXMULTI] = {0};
  Int_t                   multiCry = 0 ; /// this is total multiplicity for all crystal

  //unsigned short  pileup[MAXMULTI];
  //unsigned short     hit[MAXMULTI]; // number of hit in an event

  tree->Branch("evID", &evID, "event_ID/l"); 
  tree->Branch("multi", &multi, "multi/I"); 
  tree->Branch("id",    id, "id[multi]/I");
  tree->Branch("e",      e, "e[multi]/D");
  tree->Branch("e_t",  e_t, "e_timestamp[multi]/l");
  tree->Branch("multiCry", &multiCry, "multiplicity_crystal/I");  
  
  //open list-mode data file from PXI digitizer  
  FILE *fpr = fopen(argv[1], "r");
  long int fprsize,fprpos;
  if ( fpr == NULL) {
    fprintf(stderr, "Error, cannot open input file %s\n", argv[2]);
    return 1;
  }
 
  //get file size
  fseek(fpr, 0L, SEEK_END);
  fprsize = ftell(fpr);
  rewind(fpr);
  
  TBenchmark gClock;
  gClock.Reset();
  gClock.Start("timer");
  
  /////////////////////
  // MAIN WHILE LOOP //
  /////////////////////
  while (1) { //main while loop 

      /////////////////////////////////
      // UNPACK DATA AND EVENT BUILD //
      /////////////////////////////////
      
      long long int etime = -1; 
      long long int tdif = -1; 
      int sevtmult=0;  
      
      while (1) { //get subevents and event build for one "event" 
        
        if (fread(sub, sizeof(int)*HEADER_LENGTH, 1, fpr) != 1) break;
        
        data.chn    =  sub[0] & 0xF;     /// channel in digitizer
        data.sln    = (sub[0] & 0xF0) >> 4; /// digitizer ID
        data.crn    = (sub[0] & 0xF00) >> 8;  /// crate
        data.id     = data.crn*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD + (data.sln-BOARD_START)*MAX_CHANNELS_PER_BOARD + data.chn;   
        data.hlen   = (sub[0] & 0x1F000) >> 12;
        data.elen   = (sub[0] & 0x7FFE0000) >> 17;
        data.fcode  = (sub[0] & 0x80000000) >> 31;
        data.time   = ( (long long int)(sub[2] & 0xFFFF) << 32) + sub[1];
        data.ctime  = (sub[2] & 0x7FFF0000) >> 16;
        data.ctimef = (sub[2] & 0x80000000) >> 31;
        data.e      = (sub[3] & 0xFFFF);
        data.trlen  = (sub[3] & 0x7FFF0000) >> 16;
        data.trwlen =  data.trlen / 2;
        data.extra  = (sub[3] & 0x80000000) >> 31;   
        
        tempf = (float)data.e/RAWE_REBIN_FACTOR;// + RAND;
        data.e = (int)tempf;  
       
        //Set reference time for event building
        if (etime == -1) {
            etime = data.time;
            tdif = 0;
            multi = 0;
        }else {
            tdif = data.time - etime;
        }    
      
        //Check for end of event, rewind, and break out of while loop
        if (tdif > timeWindow) {
            fseek(fpr, -sizeof(int)*HEADER_LENGTH, SEEK_CUR); //fwrite/fread is buffered by system ; storing this in local buffer is no faster!
            break;           
        }else{
          //if within time window, fill array;
          int detID  = mapping[data.id];
          id[multi]  = detID;
          e[multi]   = data.e;
          e_t[multi] = data.time;
          multi++ ;
          if( detID < 100 ) multiCry ++;
        }    
               
        // total pileups
        if (data.fcode==1) {
            pileUpCount++;
        }
        
        
        //more data than just the header; read entire sub event, first rewind, then read data.elen
        fseek(fpr, -sizeof(int)*HEADER_LENGTH, SEEK_CUR);
        //if (fread(sub, sizeof(int)*dataBlock[sevtmult].elen, 1, fpr) != 1) break;
        if (fread(sub, sizeof(int)*data.elen, 1, fpr) != 1) break;
                              
        /**                      
        //trace
        k=0;
        for (i = dataBlock[sevtmult].hlen; i < dataBlock[sevtmult].elen; i++) {      
            dataBlock[sevtmult].tr[i - dataBlock[sevtmult].hlen + k] = sub[i] & 0x3FFF;
            dataBlock[sevtmult].tr[i - dataBlock[sevtmult].hlen + k + 1] = (sub[i]>>16) & 0x3FFF;
            k=k+1;
        } 
        
     // if (dataBlock[sevtmult].id == 4 && dataBlock[sevtmult].fcode == 1) DB(dataBlock[sevtmult].tr);
            
        //continue if no esum or qsum   
        if (dataBlock[sevtmult].hlen==HEADER_LENGTH) {
            sevtmult++;        
            continue;
        }
        
        //esum
        if (dataBlock[sevtmult].hlen==8 || dataBlock[sevtmult].hlen==16) { 
            for (i=4; i < 8; i++) {
                dataBlock[sevtmult].esum[i-4] = sub[i];
            }
        }
    
        //qsum
        if (dataBlock[sevtmult].hlen==12) { 
            for (i=4; i < 12; i++) {
                dataBlock[sevtmult].qsum[i-4] = sub[i];
            }
        }
    
        //qsum
        if (dataBlock[sevtmult].hlen==16) { 
            for (i=8; i < 16; i++) {
                dataBlock[sevtmult].qsum[i-8] = sub[i];
            }
        }    
        */
     
      } //end while loop for unpacking sub events and event building for one "event"
      if (multi==0) break; //end main WHILE LOOP when out of events 
      dataCount += multi;
      evID ++;
      
      /////////////////////////////////////
      // END UNPACK DATA AND EVENT BUILD //
      /////////////////////////////////////

      //event stats, print status every 10000 events
      if ( evID % 10000 == 0 ) {
        fprpos = ftell(fpr);
        tempf = (float)fprsize/(1024.*1024.*1024.);
        gClock.Stop("timer");
        double time = gClock.GetRealTime("timer");
        gClock.Start("timer");
        printf("Total dataBlock: \x1B[32m%llu \x1B[31m(%d%% pileup)\x1B[0m\nTotal Events: \x1B[32m%llu (%.1f <mult>)\x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\nTime used:%3.0f min %5.2f sec\033[3A\r", 
                   dataCount, (int)((100*pileUpCount)/dataCount), evID+1, (float)dataCount/((float)evID+1), (100*fprpos/fprsize), tempf,  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);
      }      

      
      outRootFile->cd();
      tree->Fill();

          
  } // end main while loop 
  /////////////////////////
  // END MAIN WHILE LOOP //
  /////////////////////////
  fprpos = ftell(fpr);
  tempf = (float)fprsize/(1024.*1024.*1024.);
  printf("Total SubEvents: \x1B[32m%llu \x1B[31m(%d%% pileup)\x1B[0m\nTotal Events: \x1B[32m%llu (%.1f <mult>)\x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\n\033[3A\n", 
  dataCount, (int)((100*pileUpCount)/dataCount), evID+1, (float)dataCount/((float)evID+1), (100*fprpos/fprsize), tempf);
           

  outRootFile->cd();
  tree->Write();
  outRootFile->Close();

  fclose(fpr);
  
  gClock.Stop("timer");
  double time = gClock.GetRealTime("timer");
  printf("\n\n==================== finished.\r\n");
  printf("Total time spend : %3.0f min %5.2f sec\n", TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);

  return 0;
}


