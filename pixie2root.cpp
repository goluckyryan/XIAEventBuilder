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

#include "mapping.h"  

/////////////////////
// RAW EVENT TYPES //
/////////////////////
struct subevent
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
    int energy;
    int extra;
    short tr[4096];
    int esum[4];
    int qsum[8];        
}; 
struct subevent subevt[MAX_ID]={0};
int sevtmult=0;
unsigned long long int sevtcount=0;
unsigned long long int pileupcount=0;
unsigned long long int evtcount=0;

int mult[1][4096]={0};

int tdifid[MAX_ID][8192]={0};

/****
int overwrite = 1;    

///////////////////////
// Write 2-byte data //
///////////////////////
void write_data2(char *filename, short *data, int xdim, int ydim, int overwrite) { //2byte per channel Write / Add to previous
    
    FILE *FP;
    int i;
    short *previous;
    
    if(!overwrite) {
        //allocate memory for 1d-array for reading in rows of 2d Radware matrix
        if ( ( previous = (short *)malloc(xdim * ydim * sizeof(short)) ) == NULL ) {
            printf("\nError, memory not allocated.\n");
            exit(1);
        }
      
        //open previous spectra file  
        if( (FP=fopen(filename, "r")) != NULL ){
            fread(previous, sizeof(short)*xdim*ydim, 1, FP);        
            fclose(FP);
            //update spectra
            for (i=0; i<xdim*ydim; i++) {
                if(previous[i] < (powf(2,sizeof(short)*8.0)-2)) 
                    data[i] = data[i] + previous[i];
                }   
            }
        else{
            printf("%s did not previously exist, creating ...\n", filename);       
        }     
   
        //Deallocate previous data
        free(previous);
    }
  
    FP=fopen(filename, "w");
    fwrite(data, sizeof(short)*xdim, ydim, FP);
    fclose(FP);            
}    


///////////////////////
// Write 4-byte data //
///////////////////////
void write_data4(char *filename, int *data, int xdim, int ydim, int overwrite) { //4byte per channel Write / Add to previous
    
    FILE *FP;
    int i;
    int *previous;
   
    if(!overwrite) {
        //allocate memory for 1d-array for reading in rows of 2d Radware matrix
        if ( ( previous = (int *)malloc(xdim * ydim * sizeof(int)) ) == NULL ) {
            printf("\nError, memory not allocated.\n");
            exit(1);
        }
      
        //open previous spectra file  
        if( (FP=fopen(filename, "r")) != NULL ){
            fread(previous, sizeof(int)*xdim*ydim, 1, FP);        
            fclose(FP);
            //update spectra
            for (i=0; i<xdim*ydim; i++) {
                if(previous[i] < (powf(2,sizeof(int)*8.0)-2)) 
                    data[i] = data[i] + previous[i];
            }   
        } 
        else{
            printf("%s did not previously exist, creating ...\n", filename);       
        }   
       
        //Deallocate previous data
        free(previous);
    }
  
    FP=fopen(filename, "w");
    fwrite(data, sizeof(int)*xdim, ydim, FP);
    fclose(FP);            
}    
******/

///////////////////////////////////
// START OF MAIN FUNCTION        //
///////////////////////////////////
int main(int argc, char **argv) {
  
  int i=0, j=0, k=0;
  float tempf=0;

  div_t e_div;
  lldiv_t lle_div;

  //temp buffer for each sub event
  unsigned int sub[MAX_SUB_LENGTH];
  memset(sub, 0, sizeof(sub));
  
  //Reference time and difference for event building
  long long int etime, tdif, idtime[MAX_ID]={0}, temptime;
  
  // Check that the corrent number of arguments were provided.
  if (argc != 2 && argc != 3 )    {
    printf("Incorrect number of arguments:\n");
    printf("%s [*.to File] [timeWindow] \n", argv[0]);
    printf("         timeWindow : number of tick, 1 tick = 10 ns. default = 100 \n");       
    return 1;
  }

  printf("=====================================\n");
  printf("===       evt --> root            ===\n");
  printf("=====================================\n");
  
  //CERN ROOT things
  TString inFileName = argv[1];
  TString outFileName = inFileName;
  
  int EVENT_BUILD_TIME = 100; 
  
  if( argc >= 3 ){
    EVENT_BUILD_TIME = atoi(argv[2]);
  }
  
  outFileName.Remove(inFileName.First('.'));
  outFileName.Append(".root");
  
  printf("  in file : %s \n", inFileName.Data());
  printf(" our file : %s \n", outFileName.Data());

  printf(" number of detector channal: %d \n", MAX_ID);
  printf("------------------------ Event building time window : %d tics = %d nsec \n", EVENT_BUILD_TIME, EVENT_BUILD_TIME*10);

  TFile * outRootFile = new TFile(outFileName, "recreate");
  outRootFile->cd();
  TTree * tree = new TTree("tree", "tree");

  unsigned long long evID = -1;
  double  energy[NCRYSTAL];  
  unsigned long long etimestamp[NCRYSTAL];
  double  bgo[NBGO]; 
  double  other[NOTHER];
  unsigned short pileup[NCRYSTAL];
  
  //const int maxMulti = 40;
  //double energy[maxMulti];  
  //unsigned timestamp[maxMulti];  
  //short detID[maxMulti];
  int multi;

  tree->Branch("evID", &evID, "event_ID/l"); 
  ///tree->Branch("detID", detID, Form("det ID[%d]/B", NCRYSTAL));
  tree->Branch("e", energy, Form("energy[%d]/D", NCRYSTAL));
  tree->Branch("t", etimestamp, Form("energy_time_stamp[%d]/l", NCRYSTAL));
  tree->Branch("p", pileup, Form("pile_up_flag[%d]/s", NCRYSTAL));
  
  tree->Branch("bgo", bgo, Form("BGO_energy[%d]/D", NBGO));
  tree->Branch("other", other, Form("other_energy[%d]/D", NOTHER));
  
  tree->Branch("multi", &multi, "multiplicity/I");
  
  //open list-mode data file from PXI digitizer  
  FILE *fpr;
  long int fprsize,fprpos;
  if ((fpr = fopen(argv[1], "r")) == NULL) {
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
      
      //CERN data clear
      for( int haha = 0; haha < NCRYSTAL; haha++){
         energy[haha] = TMath::QuietNaN();
         etimestamp[haha] = 0;
         pileup[haha] = 0;
      }
      for( int haha = 0; haha < NBGO; haha++) bgo[haha] = TMath::QuietNaN();
      for( int haha = 0; haha < NOTHER; haha++) other[haha] = TMath::QuietNaN();
      multi = 0;
      evID++;	
      
      etime=-1; tdif=-1; sevtmult=0;  
      //memset(&subevt, 0, sizeof(subevt)); //not needed since everything is redefined (except maybe trace on pileup evts)
      while (1) { //get subevents and event build for one "event" 
        
       // memset(&subevt[sevtmult], 0, sizeof(subevt[sevtmult])); //not needed since everything is redefined (except maybe trace on pileup evts)
        
        //read 4-byte header
        if (fread(sub, sizeof(int)*HEADER_LENGTH, 1, fpr) != 1) break;
        subevt[sevtmult].chn = sub[0] & 0xF;     /// channel in digitizer
        subevt[sevtmult].sln = (sub[0] & 0xF0) >> 4; /// digitizer ID
        subevt[sevtmult].crn = (sub[0] & 0xF00) >> 8;  /// crate
        subevt[sevtmult].id = subevt[sevtmult].crn*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD + (subevt[sevtmult].sln-BOARD_START)*MAX_CHANNELS_PER_BOARD + subevt[sevtmult].chn;   
        subevt[sevtmult].hlen = (sub[0] & 0x1F000) >> 12;
        subevt[sevtmult].elen = (sub[0] & 0x7FFE0000) >> 17;
        subevt[sevtmult].fcode = (sub[0] & 0x80000000) >> 31;
        subevt[sevtmult].time = ( (long long int)(sub[2] & 0xFFFF) << 32) + sub[1];
        subevt[sevtmult].ctime = (sub[2] & 0x7FFF0000) >> 16;
        subevt[sevtmult].ctimef = (sub[2] & 0x80000000) >> 31;
        subevt[sevtmult].energy = (sub[3] & 0xFFFF);
        subevt[sevtmult].trlen = (sub[3] & 0x7FFF0000) >> 16;
        subevt[sevtmult].trwlen = subevt[sevtmult].trlen / 2;
        subevt[sevtmult].extra = (sub[3] & 0x80000000) >> 31;       
 
        //rebin raw trap energy from 32k to ....            
        tempf = (float)subevt[sevtmult].energy/RAWE_REBIN_FACTOR;// + RAND;
        subevt[sevtmult].energy = (int)tempf;                        
 
        //check lengths (sometimes all of the bits for trace length are turned on ...)
       /* if (subevt[sevtmult].elen - subevt[sevtmult].hlen != subevt[sevtmult].trwlen) {
            printf("SEVERE ERROR: event, header, and trace length inconsistencies found\n");
            printf("event length = %d\n", subevt[sevtmult].elen);
            printf("header length = %d\n", subevt[sevtmult].hlen);
            printf("trace length = %d\n", subevt[sevtmult].trwlen);  
            printf("Extra = %d\n", subevt[sevtmult].extra); 
            printf("fcode = %d\n", subevt[sevtmult].fcode);              
            //sleep(1);          
            //return 0;
        } */ 
        
        
        //CERN fill tree
        
        ///========== need a mapping, can reduce the array size, speed up.
        
        int ch = map[subevt[sevtmult].id];
        if ( 0 <= ch && ch < NCRYSTAL ){
            energy[ch] = subevt[sevtmult].energy;
            etimestamp[ch] = subevt[sevtmult].time;
            pileup[ch] = subevt[sevtmult].fcode;
            multi++;
        }
        if ( 100 <= ch && ch < 100 + NBGO ){
            bgo[ch-100] = subevt[sevtmult].energy;
        }
        if ( 200 <= ch && ch < 200 + NOTHER){
            other[ch-200] = subevt[sevtmult].energy;
        }
       
        //Set reference time for event building
        if (etime == -1) {
            etime = subevt[sevtmult].time;
            tdif = 0;
        }
        else {
            tdif = subevt[sevtmult].time - etime;
            //if (tdif < 0) {
            //    printf("SEVERE ERROR: tdiff < 0, file must be time sorted\n");
            //    printf("etime = %lld, time = %lld, and tdif = %lld\n", etime, subevt[sevtmult].time, tdif);                
            //    return 0;   
            //}    
        }    
      
        //Check for end of event, rewind, and break out of while loop
        if (tdif > EVENT_BUILD_TIME) {
            fseek(fpr, -sizeof(int)*HEADER_LENGTH, SEEK_CUR); //fwrite/fread is buffered by system ; storing this in local buffer is no faster!
            break;           
        }    
               
        
        //time between sequential events for a single channel ; useful for determining optimal event building time
        temptime = (subevt[sevtmult].time - idtime[subevt[sevtmult].id])/100; //rebin to 1 micro-second
        if ( temptime >= 0 && temptime < 8192) {
            tdifid[subevt[sevtmult].id][temptime]++;
        }    
        idtime[subevt[sevtmult].id]=subevt[sevtmult].time; //store time for next subevent of channel    
    
        // total pileups
        if (subevt[sevtmult].fcode==1) {
            pileupcount++;
        }
        
        
        //more data than just the header; read entire sub event
        fseek(fpr, -sizeof(int)*HEADER_LENGTH, SEEK_CUR);
        if (fread(sub, sizeof(int)*subevt[sevtmult].elen, 1, fpr) != 1) break;
                              
        /*                      
        //trace
        k=0;
        for (i = subevt[sevtmult].hlen; i < subevt[sevtmult].elen; i++) {      
            subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k] = sub[i] & 0x3FFF;
            subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k + 1] = (sub[i]>>16) & 0x3FFF;
            k=k+1;
        } 
        
     // if (subevt[sevtmult].id == 4 && subevt[sevtmult].fcode == 1) DB(subevt[sevtmult].tr);
            
        //continue if no esum or qsum   
        if (subevt[sevtmult].hlen==HEADER_LENGTH) {
            sevtmult++;        
            continue;
        }
        
        //esum
        if (subevt[sevtmult].hlen==8 || subevt[sevtmult].hlen==16) { 
            for (i=4; i < 8; i++) {
                subevt[sevtmult].esum[i-4] = sub[i];
            }
        }
    
        //qsum
        if (subevt[sevtmult].hlen==12) { 
            for (i=4; i < 12; i++) {
                subevt[sevtmult].qsum[i-4] = sub[i];
            }
        }
    
        //qsum
        if (subevt[sevtmult].hlen==16) { 
            for (i=8; i < 16; i++) {
                subevt[sevtmult].qsum[i-8] = sub[i];
            }
        }    
        */
        sevtmult++;
     
      } //end while loop for unpacking sub events and event building for one "event"
      if (sevtmult==0) break; //end main WHILE LOOP when out of events 
      mult[0][sevtmult]++; //Histogram raw sub event multiplicity 
      sevtcount += sevtmult;
      evtcount++; //event-built number
      /////////////////////////////////////
      // END UNPACK DATA AND EVENT BUILD //
      /////////////////////////////////////

      //event stats, print status every 10000 events
      lle_div=lldiv(evtcount,10000);
      if ( lle_div.rem == 0 ) {
        fprpos = ftell(fpr);
        tempf = (float)fprsize/(1024.*1024.*1024.);
        gClock.Stop("timer");
        double time = gClock.GetRealTime("timer");
        gClock.Start("timer");
        printf("Total SubEvents: \x1B[32m%llu \x1B[31m(%d%% pileup)\x1B[0m\nTotal Events: \x1B[32m%llu (%.1f <mult>)\x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\nTime used:%3.0f min %5.2f sec\033[3A\r", 
                   sevtcount, (int)((100*pileupcount)/sevtcount), evtcount, (float)sevtcount/(float)evtcount, (100*fprpos/fprsize), tempf,  TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);
      }      


      //cern fill tree
      outRootFile->cd();
      tree->Fill();

          
  } // end main while loop 
  /////////////////////////
  // END MAIN WHILE LOOP //
  /////////////////////////
  fprpos = ftell(fpr);
  tempf = (float)fprsize/(1024.*1024.*1024.);
  printf("Total SubEvents: \x1B[32m%llu \x1B[31m(%d%% pileup)\x1B[0m\nTotal Events: \x1B[32m%llu (%.1f <mult>)\x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\n\033[3A\n", 
  sevtcount, (int)((100*pileupcount)/sevtcount), evtcount, (float)sevtcount/(float)evtcount, (100*fprpos/fprsize), tempf);
           
  //cern save root
  outRootFile->cd();
  tree->Write();
  outRootFile->Close();
  
  fclose(fpr);
  
  gClock.Stop("timer");
  double time = gClock.GetRealTime("timer");
  printf("\n==================== finished.\r\n");
  printf("Total time spend : %3.0f min %5.2f sec\n", TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);
  
  
  return 0;
}


