/**********************************************************/
/* PXI SCAN CODE -- J.M. Allmond (ORNL) -- July 2016      */
/*                                                        */
/* !unpak data from Pixie-16 digitizers, event build,     */
/* !and create detctors and user defined spectra          */
/*                                                        */
/* gcc -o pxi-scan pxi-scan.c                             */
/* ./pxi-scan -op datafile calibrationfile mapfile        */
/*                                                        */
/* ..... calibration file optional                        */
/* ..... map file optional                                */
/* ..... u for update spectra                             */
/* ..... o for overwrite spectra                          */
/* ..... p for print realtime stats                       */
/**********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define PRINT_CAL 1
#define PRINT_MAP 1

#define DB(x) fwrite(x, sizeof(x), 1, debugfile);
#define DEBUGFN "debug.mat"

#define RAND ((float) rand() / ((unsigned int) RAND_MAX + 1))   // random number in interval (0,1)
#define TRUE  1
#define FALSE 0

#define LINE_LENGTH 120

#define MAX_CRATES 2
#define MAX_BOARDS_PER_CRATE 13
#define MAX_CHANNELS_PER_BOARD 16
#define BOARD_START 2

#define MAX_ID MAX_CRATES*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD

#define HEADER_LENGTH 4     //unit = words with 4 bytes per word
#define MAX_SUB_LENGTH 2016 //unit = words with 4 bytes per word ; 2004 --> 40 micro-second trace + 4 word header

#define EVENT_BUILD_TIME 109 // 100 = 1 micro-second ; should be < L + G ~ 5.04 us (note 0.08 us scale factor in set file)

#define RAWE_REBIN_FACTOR 2.0 // Rebin 32k pixie16 spectra to something smaller to fit better into 8k.

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



//////////////////////////////////////////
// INPUT CALIBRATION AND MAP PARAMETERS //
//////////////////////////////////////////

float ecal[MAX_ID][2]={0};
float tcal[MAX_ID][2]={0};
float new_gain=1.0;

char map2type[MAX_ID]={0};
int map2det[MAX_ID]={0};
int map2deti[MAX_ID]={0};
float mapangles[MAX_ID][2]={0}; //theta and phi
float mapanglesi[MAX_ID][2]={0}; //theta and phi
float mapangles1[MAX_ID][2]={0}; //theta and phi
float mapangles2[MAX_ID][2]={0};//theta and phi

////////////////////
// DETECTOR TYPES //
////////////////////

//G = Ge
#define MAX_GE 16 // max number of Ge detectors
#define MAX_GE_XTL 4 // max number of crystals per Ge detector
#define MAX_GE_SEG 3 // max number of segments per Ge detector
#define MAX_GE_BGO 1 // max number of BGO PMTs per Ge detector 
#define GE_BGO_SUPPRESSION TRUE
struct gdetector
{    
    int xmult;
    int xid[MAX_GE_XTL]; 
    int xe[MAX_GE_XTL];
    long long int xt[MAX_GE_XTL]; 
    int xct[MAX_GE_XTL];   
    float xtheta[MAX_GE_XTL][4]; 
    float xphi[MAX_GE_XTL][4];
    bool xpileup[MAX_GE_XTL];
    bool xsuppress[MAX_GE_XTL];
    int xsubevtid[MAX_GE_XTL];

    int smult;
    int sid[MAX_GE_SEG]; 
    int se[MAX_GE_SEG];
    long long int st[MAX_GE_SEG];
    int sct[MAX_GE_SEG];
    float stheta[MAX_GE_SEG][4]; 
    float sphi[MAX_GE_SEG][4];
    bool spileup[MAX_GE_SEG];
    bool ssuppress[MAX_GE_SEG];
    int ssubevtid[MAX_GE_SEG];

    int bgomult;
    int bgoid[MAX_GE_BGO];
    int bgoe[MAX_GE_BGO];
    long long int bgot[MAX_GE_BGO];
    int bgoct[MAX_GE_BGO];
    float bgotheta[MAX_GE_XTL][4]; 
    float bgophi[MAX_GE_XTL][4];
    bool bgopileup[MAX_GE_BGO];
    int bgosubevtid[MAX_GE_BGO];

    int id;
    int energy;
    int edop;
    long long int time;   
    int ctime;  
    float theta[3]; //det, xtl, or seg angle
    float phi[3]; //det, xtl, or seg angle  
    bool suppress; //at least one xtl was suppressed by bgo
    bool pileup; //two or more unspressed xtls but at least one had pileup
    bool nonprompt; //two or more unspressed xtls but at least one was non-prompt with first xtl
    bool clean;
}; 
struct gdetector ge[MAX_GE]={0};
int gmult=0;
unsigned long long int gcount=0;

//S = Si
// .......

///////////////////////////////////////
// SPECTRA and FILE NAME DEFINITIONS // 
///////////////////////////////////////
//All spectra are considered two-dimensional arrays
//Must add "write spectra" at end of file 

//[Y-dim][X-dim]

////////////////
//Event Spectra
////////////////

#define HIT "hit.spn"
int hit[2][4096]={0}; //first for all hits, second for pilup hits

#define MULT "mult.spn" //total detector multiplicity for one event
int mult[1][4096]={0};

#define TDIFID "tdif.sec" //time diference between sequential events of a single channel ; 1 micro-second bins
int tdifid[MAX_ID][8192]={0};

#define E_RAW "e_raw.sec"
int e_raw[MAX_ID][8192]={0};

#define E_CAL "e_cal.sec"
int e_cal[MAX_ID][8192]={0};

#define TEVT_RAW "tevt_raw.sec" // 10 second bins
int tevt_raw[MAX_ID][8192]={0};

#define TEVT_CAL "tevt_cal.sec"
int tevt_cal[MAX_ID][8192]={0}; // 10 second bins

#define TCFD_RAW "tcfd_raw.sec"
int tcfd_raw[MAX_ID][8192]={0};

#define TDIF_RAW "tdif_raw.spn" // 10 nano-second bins
int tdif_raw[MAX_ID][4096]={0};

#define TDIF_CAL "tdif_cal.spn" // 10 nano-second bins
int tdif_cal[MAX_ID][4096]={0};

#define TDIF_CAL0_ETHRESH "tdif_cal0_ethresh.spn" //time difference relative to channel 0 ; 10 nano-second bins
int tdif_cal0_ethresh[MAX_ID][4096]={0};

#define IDID "idid.spn" //id vs id correlation hit pattern 
int idid[MAX_ID][MAX_ID]={0};

////////////////////////////
//Detector Processed Spectra
////////////////////////////

//Ge
#define GE_BGO_TDIF "ge_bgo_tdif.spn"
int ge_bgo_tdif[MAX_GE][4096]={0};

#define GE_XTL_TDIF "ge_xtl_tdif.spn"
int ge_xtl_tdif[MAX_GE][4096]={0};

#define GE_XTL_TDIF_ETHRESH "ge_xtl_tdif_ethresh.spn"
int ge_xtl_tdif_ethresh[MAX_GE][4096]={0};

#define GE_SPE_XTL "ge_spe_xtl.sec"
int ge_spe_xtl[MAX_GE*MAX_GE_XTL][8192]={0};

#define GE_SPE "ge_spe.sec"
int ge_spe[MAX_GE][8192]={0};

#define GE_SPE_CLEAN "ge_spe_clean.sec"
int ge_spe_clean[MAX_GE][8192]={0};

//trinity
#define PID "pid.spn"
int pid[4096][4096]={{0}};
#define PID_EVSP "pid_evsp.spn"
int pid_evsp[4096][4096]={{0}};
#define PID_EVST "pid_evst.spn"
int pid_evst[4096][4096]={{0}};
#define PID_EVSTAU "pid_evstau.spn"
int pid_evstau[4096][4096]={{0}};
#define PID_EVSR "pid_evsr.spn"
int pid_evsr[4096][4096]={{0}};
#define PID_TAUVSR "pid_tauvsr.spn"
int pid_tauvsr[4096][4096]={{0}};
#define GAGGDT "gaggdt.spn"
int gaggdt[4096][4096]={{0}};
#define SPTRACE "sptrace.spn"
int sptrace[1][4096]={{0}};

//////////////////////
//Final User Spectra
//////////////////////

//gamma-gamma
#define GG_TDIF "gg_tdif.spn"
int gg_tdif[4096][4096]={0};

#define GG_PROMPT "gg_prompt.spn"
int gg_prompt[4096][4096]={0};

#define GG_NPROMPT "gg_nprompt.spn"
int gg_nprompt[4096][4096]={0};




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




///////////////////////////////////
// START OF MAIN FUNCTION        //
///////////////////////////////////
int main(int argc, char **argv) {
  
  int i=0, j=0, k=0;
  float tempf=0;
  int max1=0, min1=0;
  int max2=0, min2=0;
  int maxid1=-1, minid1=-1;
  int maxid2=-1, minid2=-1;
  div_t e_div;
  lldiv_t lle_div;

  int overwrite = 1;    

  double etrace0,etrace1,btrace0,btrace1;
  double ptrace0,ptrace1,ttrace0,ttrace1,tautrace0,tautrace1;
  int dbcount = 0;
  long long int strace[500];
  memset(strace, 0, sizeof(strace));

  //temp buffer for each sub event
  unsigned int sub[MAX_SUB_LENGTH];
  memset(sub, 0, sizeof(sub));
  
  //Reference time and difference for event building
  long long int etime, tdif, idtime[MAX_ID]={0}, temptime;
  
   
  // Check that the corrent number of arguments were provided.
  if (argc<3)    {
    printf("Incorrect number of arguments:\n");
    printf("%s -op datafile calibrationfile mapfile \n", argv[0]);
    printf("\n .... calibration file is optional\n");
    printf(" .... map file is optional\n");
    printf(" .... o for overwrite spectra\n");
    printf(" .... u for update spectra\n");
    printf(" .... p for print realtime stats\n");
    return 1;
  }
    
  if(strstr(argv[1], "u") != NULL) {
    overwrite = 0;
    printf("Updating Spectra\n");
  }
  else {
    printf("Overwriting Spectra\n");
  }  
  
  //open list-mode data file from PXI digitizer  
  FILE *fpr;
  long int fprsize,fprpos;
  if ((fpr = fopen(argv[2], "r")) == NULL) {
    fprintf(stderr, "Error, cannot open input file %s\n", argv[2]);
    return 1;
  }
 
  //get file size
  fseek(fpr, 0L, SEEK_END);
  fprsize = ftell(fpr);
  rewind(fpr);
  
   
  //open debug file for streaming an 1d array
  FILE *debugfile;
  if ((debugfile = fopen(DEBUGFN, "w")) == NULL) {
    fprintf(stderr, "Error, cannot open %s\n", DEBUGFN);
    return 1;
  }   
   

  //buffer for reading in text files for calibrations and maps below
  char line[LINE_LENGTH];


  //open energy and time calibration file (e.g., *.ca3 file)  
  int calid=0; 
  float caloffset=0.0, calgain=0.0;
  int firstcal=0, necal=0, ntcal=0;
  
  FILE *fprcal;
  
  for (i=0; i<MAX_ID; i++) {
    ecal[i][0] = caloffset;
    ecal[i][1] = calgain;
    tcal[i][0] = caloffset;
    tcal[i][1] = calgain;     
  }  
  
  if (argc >= 4) {

    if ((fprcal = fopen(argv[3], "r")) == NULL) {
        fprintf(stderr, "Error, cannot open input file %s\n", argv[3]);
        return 1;
    } 
    
    printf("%s loaded!\n", argv[3]);

	while(fgets(line, LINE_LENGTH, fprcal) != NULL){
	    calid=0; caloffset=0; calgain=0;
        for(i=0; i<LINE_LENGTH; i++){
            if(line[i]=='#'){
                if(PRINT_CAL)printf("%s", line);
                break;
            }
            else if(line[i]>=0){
                if (firstcal==0) {
                    sscanf(line,"%f\n", &new_gain);    
	                if(PRINT_CAL) printf("%f\n", new_gain);                                       
                    firstcal=1;
                    break;
                }    
	            sscanf(line,"%d\t%f\t%f\n", &calid, &caloffset, &calgain);
	            if(PRINT_CAL) printf("%d\t%.4f\t%.4f\n", calid, caloffset, calgain);
		        if(calid >=0 && calid < MAX_ID) {
                    ecal[calid][0] = caloffset;
                    ecal[calid][1] = calgain;
                    necal++;
          		    break;
			    }
		        if(calid >=1000 && calid < 1000+MAX_ID) {
                    tcal[calid-1000][0] = caloffset;
                    tcal[calid-1000][1] = calgain;
                    ntcal++;
          		    break;
			    }			    
			    else {
				    printf("Error in reading %s : bad id or format\n", argv[3]);
				    return -1;
			    }
            }
	        else if(line[i]=='\n'){
                if(PRINT_CAL) printf("\n");
                break;
            }
            else {
                continue;
            }

        }
        	memset(line, 0, LINE_LENGTH);
 	  }
 	  
 	fclose(fprcal);
    printf("read %d energy calibrations\n", necal);
    printf("read %d time calibrations\n", ntcal);
  }
  


  //open ID->Detector map file (e.g., *.map file)  
  int mapid=0, detid=0, detidi; 
  char dettype=0;
  float theta=0, phi=0, thetai=0, phii=0, theta1=0, phi1=0, theta2=0, phi2=0; 
  FILE *fprmap;
  int nmmap=0;
  if (argc >= 5) {

    if ((fprmap = fopen(argv[4], "r")) == NULL) {
        fprintf(stderr, "Error, cannot open input file %s\n", argv[4]);
        return 1;
    } 

    printf("%s loaded!\n", argv[4]);

	while(fgets(line, LINE_LENGTH, fprmap) != NULL){
	    mapid=0; dettype=0; thetai=0; phii=0; theta=0; phi=0; theta1=0; phi1=0; theta2=0; phi2=0; 
        for(i=0; i<LINE_LENGTH; i++){
            if(line[i]=='#'){
                if(PRINT_MAP)printf("%s", line);
                break;
            }
            else if(line[i]>=0){   
	            sscanf(line,"%d\t%c\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", &mapid, &dettype, &detid, &detidi, &theta, &phi, &thetai, &phii, &theta1, &phi1, &theta2, &phi2);
	            if(PRINT_MAP) printf("%d\t%c\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", mapid, dettype, detid, detidi, theta, phi, thetai, phii, theta1, phi1, theta2, phi2);
		        if(mapid >=0 && mapid < MAX_ID) {
                    map2type[mapid] = dettype;
                    map2det[mapid] = detid;
                    map2deti[mapid] = detidi;
                    mapangles[mapid][0] = theta;
                    mapangles[mapid][1] = phi;
                    mapanglesi[mapid][0]=thetai; 
                    mapanglesi[mapid][1]=phii; 
                    mapangles1[mapid][0]= theta1; 
                    mapangles1[mapid][1]= phi1; 
                    mapangles2[mapid][0]= theta2; 
                    mapangles2[mapid][1]= phi2; 
                    nmmap++;
          		    break;
			    }			    
			    else {
				    printf("Error in reading %s : bad id or format\n", argv[4]);
				    return -1;
			    }
            }
	        else if(line[i]=='\n'){
                if(PRINT_MAP) printf("\n");
                break;
            }
            else {
                continue;
            }

        }
        	memset(line, 0, LINE_LENGTH);
 	  }
 	  
 	fclose(fprmap);
    printf("read %d id maps\n", nmmap);
  }




  /////////////////////
  // MAIN WHILE LOOP //
  /////////////////////
  while (1) { //main while loop 


      /////////////////////////////////
      // UNPACK DATA AND EVENT BUILD //
      /////////////////////////////////
      etime=-1; tdif=-1; sevtmult=0;  
      //memset(&subevt, 0, sizeof(subevt)); //not needed since everything is redefined (except maybe trace on pileup evts)
      while (1) { //get subevents and event build for one "event" 
        
       // memset(&subevt[sevtmult], 0, sizeof(subevt[sevtmult])); //not needed since everything is redefined (except maybe trace on pileup evts)
        
        //read 4-byte header
        if (fread(sub, sizeof(int)*HEADER_LENGTH, 1, fpr) != 1) break;
        subevt[sevtmult].chn = sub[0] & 0xF;
        subevt[sevtmult].sln = (sub[0] & 0xF0) >> 4;
        subevt[sevtmult].crn = (sub[0] & 0xF00) >> 8;
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
        
       
        //Set reference time for event building
        if (etime == -1) {
            etime = subevt[sevtmult].time;
            tdif = 0;
        }
        else {
            tdif = subevt[sevtmult].time - etime;
            if (tdif < 0) {
                printf("SEVERE ERROR: tdiff < 0, file must be time sorted\n");
                printf("etime = %lld, time = %lld, and tdif = %lld\n", etime, subevt[sevtmult].time, tdif);                
                return 0;   
            }    
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
        
        //Histogram raw spectra
        hit[0][subevt[sevtmult].id]++;
        if (subevt[sevtmult].fcode==1) 
        hit[1][subevt[sevtmult].id]++;     
  
        if (subevt[sevtmult].energy >= 0 && subevt[sevtmult].energy < 8192)
        e_raw[subevt[sevtmult].id][subevt[sevtmult].energy]++;
 
        if (subevt[sevtmult].time/1000000000 >= 0 && subevt[sevtmult].time/1000000000 < 8192) // rebin to 10 seconds
        tevt_raw[subevt[sevtmult].id][subevt[sevtmult].time/1000000000]++; // rebin to 10 seconds

        if (subevt[sevtmult].ctime >= 0 && subevt[sevtmult].ctime < 8192)
        tcfd_raw[subevt[sevtmult].id][subevt[sevtmult].ctime]++;

        if (tdif >= 0 && tdif < 4096 && sevtmult!=0)
        tdif_raw[subevt[sevtmult].id][tdif]++;


        //if CFD is enabled, ctime will be non-zero
        //tempf = (float)subevt[sevtmult].ctime*10.0/32768.0;
        //subevt[sevtmult].time = subevt[sevtmult].time + (long long int)tempf;    

        //Calibrate energy and time
        tempf = ((float)subevt[sevtmult].energy*ecal[subevt[sevtmult].id][1] + ecal[subevt[sevtmult].id][0])/new_gain;// + RAND;       
        subevt[sevtmult].energy = (int)tempf; 
        //subevt[sevtmult].time += (long long int)tcal[subevt[sevtmult].id][0];       
	
        //Histogram calibrated spectra
        if (subevt[sevtmult].energy >= 0 && subevt[sevtmult].energy < 8192)
        e_cal[subevt[sevtmult].id][subevt[sevtmult].energy]++;        
        
        if (subevt[sevtmult].time/1000000000 >= 0 && subevt[sevtmult].time/1000000000 < 8192)        
        tevt_cal[subevt[sevtmult].id][subevt[sevtmult].time/1000000000]++;
         
        //continue on if no trace, esum, or qsum
        if (subevt[sevtmult].hlen==HEADER_LENGTH && subevt[sevtmult].trwlen==0 ) {
            sevtmult++;
            continue;
        }
        
        //more data than just the header; read entire sub event
        fseek(fpr, -sizeof(int)*HEADER_LENGTH, SEEK_CUR);
        if (fread(sub, sizeof(int)*subevt[sevtmult].elen, 1, fpr) != 1) break;
                              
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
        
        sevtmult++;
     
      } //end while loop for unpacking sub events and event building for one "event"
      if (sevtmult==0) break; //end main WHILE LOOP when out of events 
      mult[0][sevtmult]++; //Histogram raw sub event multiplicity 
      sevtcount += sevtmult;
      evtcount++; //event-built number
      /////////////////////////////////////
      // END UNPACK DATA AND EVENT BUILD //
      /////////////////////////////////////

  	//int GAGG1=58;
	//int GAGG2=59;

  	int GAGG1=60;
	int GAGG2=70;

      if (sevtmult == 2 && (subevt[0].id == GAGG1 || subevt[0].id == GAGG2) && (subevt[1].id == GAGG1 || subevt[1].id == GAGG2) && ((subevt[0].time - subevt[1].time) + 2000 > 1992) && ((subevt[0].time - subevt[1].time) + 2000 < 2005) ) {

	if (evtcount < 200) DB(subevt[0].tr);

        etrace0 = 0;
        btrace0 = 0;
        etrace1 = 0;
        btrace1 = 0;
	ptrace0 = 0;
	ptrace1 = 0;
	ttrace0 = 0;
	ttrace1 = 0;
	tautrace0 = 0;
	tautrace1 = 0;

	for (i=0; i<60; i++) {
            btrace0 = btrace0 + subevt[0].tr[i];
            btrace1 = btrace1 + subevt[1].tr[i];
	}
	btrace0 = btrace0/60.;
        btrace1 = btrace1/60.;

	for (i=60; i<173; i++) { //180-230
            etrace0 = etrace0 + subevt[0].tr[i] - btrace0;
            etrace1 = etrace1 + subevt[1].tr[i] - btrace1;
	    tautrace0 = tautrace0 + (subevt[0].tr[i] - btrace0)*i;
            tautrace1 = tautrace1 + (subevt[1].tr[i] - btrace1)*i;
       	}
	tautrace0 = tautrace0 / etrace0;
        tautrace1 = tautrace1 / etrace1;
        etrace0 = etrace0 / 2.;
        etrace1 = etrace1 / 2.;


        //peak sum
        for (i=77; i<93; i++) { //180-223 // 197-213
            ptrace0 = ptrace0 + subevt[0].tr[i] - btrace0;
            ptrace1 = ptrace1 + subevt[1].tr[i] - btrace1;
        }
//        for (i=193; i<199; i++) { //180-223 // 197-213
//            ptrace0 = ptrace0 - (subevt[0].tr[i] - btrace0);
//            ptrace1 = ptrace1 - (subevt[1].tr[i] - btrace1);
//        }
        ptrace0 = ptrace0 / 2.;
        ptrace1 = ptrace1 / 2.;
	ptrace0 = (ptrace0 + ptrace1)/10;

        //tail sum
        for (i=101; i<160; i++) { //223-293 // 221-280
            ttrace0 = ttrace0 + subevt[0].tr[i] - btrace0;
            ttrace1 = ttrace1 + subevt[1].tr[i] - btrace1;
	}
        ttrace0 = ttrace0 / 2.;
        ttrace1 = ttrace1 / 2.;
	ttrace0 = (ttrace0 + ttrace1)/10;

        if ((int)ptrace0 > 0 && (int)ptrace0 < 4096 && (int)ttrace0 > 0 && (int)ttrace0 < 4096) pid[(int)ptrace0][(int)ttrace0]++;
        if ((int)etrace0 > 0 && (int)etrace0 < 4096 && (int)ptrace0 > 0 && (int)ptrace0 < 4096) pid_evsp[(int)etrace0][(int)ptrace0]++;
        if ((int)etrace0 > 0 && (int)etrace0 < 4096 && (int)ttrace0 > 0 && (int)ttrace0 < 4096) pid_evst[(int)etrace0][(int)ttrace0]++;
        if ((int)etrace0 > 0 && (int)etrace0 < 4096 && (int)tautrace0 > 0 && (int)tautrace0 < 4096) pid_evstau[(int)etrace0][(int)tautrace0]++;
        if ((int)(etrace0) > 0 && (int)(etrace0) < 4096 && (int)((100.*ttrace0)/ptrace0) > 0 && (int)((100.*ttrace0)/ptrace0) < 4096) pid_evsr[(int)(etrace0)][(int)((100.*ttrace0)/ptrace0)]++;
        if ((int)(tautrace0) > 0 && (int)(tautrace0) < 4096 && (int)((100.*ttrace0)/ptrace0) > 0 && (int)((100.*ttrace0)/ptrace0) < 4096) pid_tauvsr[(int)(tautrace0)][(int)((100.*ttrace0)/ptrace0)]++;



/*
	for (i=1; i<6; i++) {
	    etrace0 = etrace0 + subevt[0].qsum[i];
            etrace1 = etrace1 + subevt[1].qsum[i];
	    //printf("subevt[0].qsum[%d]=%d\n", i, subevt[0].qsum[i]);
	}
	etrace0 = etrace0 - subevt[0].qsum[0]*0.5/1.8;
        etrace1 = etrace1 - subevt[1].qsum[0]*0.5/1.8;

	etrace0 = etrace0 / 1;
        etrace1 = etrace1 / 1;
*/

	//return 0;

        //printf("etrace = %f\n", etrace);
        if ((int)etrace0 >0 && (int)etrace0 < 8192) e_raw[200-GAGG1+subevt[0].id][(int)etrace0]++;
        if ((int)etrace1 >0 && (int)etrace1 < 8192) e_raw[200-GAGG1+subevt[1].id][(int)etrace1]++;
        if ((int)etrace0 >0 && (int)etrace0 < 8192 && (int)ttrace0>500/10 && (int)ptrace0 > 978/10 && (int)ptrace0 < 1270/10 ) {
		e_raw[203][(int)etrace0]++;
        	//if (evtcount < 200) DB(subevt[0].tr);
	}
      //	if ((int)ttrace0<40 && (int)ptrace0 > 131) {
                if ((int)etrace0>1350 && (int)etrace0<1400 && (int)((100.*ttrace0)/ptrace0) < 200 ) {
//			DB(subevt[0].tr);
			for (i=0; i<500; i++) {
				strace[i]=strace[i]+subevt[0].tr[i];
			}
			dbcount++;
		}
      //  }

        etrace1 = etrace0 + etrace1;
        etrace1 = etrace1/2.0;
        if ((int)etrace1 >0 && (int)etrace1 < 8192) {
          e_raw[202][(int)etrace1]++;
        }

	if ((subevt[0].time - subevt[1].time) + 2000 > 0 && (subevt[0].time - subevt[1].time) + 2000 < 4096 && (int)etrace1 >0 && (int)etrace1 < 4096)
	gaggdt[(subevt[0].time - subevt[1].time) + 2000][(int)etrace1]++;

      }

      //skip detector building below if no map file
      if (argc >= 5) {        


      ////////////////////////////////////////
      // MAP SUB EVENTS INTO DETECTOR TYPES //
      ////////////////////////////////////////
      memset(&ge, 0, sizeof(ge)); //This is needed but could be replaced by setting suppress, pileup, nonprompt, clean, x/s/bmult to zero at start of loop!
      
      for (i=0; i<sevtmult; i++) { 
        
        for (j=0; j<sevtmult; j++) {
            if (i!=j)
            idid[subevt[i].id][subevt[j].id]++;
        }
        
        //printf("i=%d, sevtmult=%d, subevt[i].id=%d, map2type[subevt[i].id]=%c, subevt[i].energy=%d, subevt[i].time=%lld\n", i, sevtmult, subevt[i].id, map2type[subevt[i].id], subevt[i].energy, subevt[i].time);
        //fflush(stdout);
          
        //if CFD is enabled, ctime will be non-zero
        //tempf = (float)subevt[i].ctime*10.0/32768.0;
        //tempf = (float)subevt[i].ctime*1.0/32768.0; //should make no difference ; need decimals or convert time to 1 ns bins
        //printf("%lld e=%d ", subevt[i].time, subevt[i].energy);
        //subevt[i].time = subevt[i].time + (long long int)tempf;  
        //printf("%lld and %f\n", subevt[i].time, tempf); 
  
               
        //Histogram calibrated tdif spectra (Do here to keep subevents within an event time ordered during event build above)
        if (i==0) etime = subevt[i].time;
        tdif = abs(subevt[i].time - etime); 
        if (tdif >= 0 && tdif < 4096 && i!=0) tdif_cal[subevt[i].id][tdif]++;
        //if (tdif > 20 && tdif < 50) DB(subevt[i].tr);
        
        //tdif with respect to channel id 0
        if (i!=0 && subevt[0].id==0 && subevt[0].energy >= 10 && subevt[0].energy <= 8000 && subevt[i].energy >= 10  && subevt[i].energy <= 8000) {
            if (tdif >= 0 && tdif < 4096 ) tdif_cal0_ethresh[subevt[i].id][tdif]++;   
        }            
        

        /////////////////////
        // G = Ge Detector //
        /////////////////////              
        if ( map2type[subevt[i].id] == 'G' ) { //Keep G and ge or switch to C and clover/cl?
        
            if (map2deti[subevt[i].id] > 0 && map2deti[subevt[i].id] <= MAX_GE_XTL) { //Ge crystal 
               if ( ge[map2det[subevt[i].id]].xmult >= MAX_GE_XTL) {               
                    printf("SEVERE ERROR: Same Ge(xtl) twice within event build; Make event build time smaller!\n");    
                    printf("ge[map2det[subevt[i].id]].xmult=%d, ge[map2det[subevt[i].id]].smult=%d, ge[map2det[subevt[i].id]].bgomult=%d\n", ge[map2det[subevt[i].id]].xmult, ge[map2det[subevt[i].id]].smult,ge[map2det[subevt[i].id]].bgomult);
                    continue;
                    //return -1;
                }   
                ge[map2det[subevt[i].id]].xid[ge[map2det[subevt[i].id]].xmult] = map2deti[subevt[i].id];               
                ge[map2det[subevt[i].id]].xe[ge[map2det[subevt[i].id]].xmult] = subevt[i].energy;               
                ge[map2det[subevt[i].id]].xt[ge[map2det[subevt[i].id]].xmult] = subevt[i].time;                           
                ge[map2det[subevt[i].id]].xct[ge[map2det[subevt[i].id]].xmult] = subevt[i].ctime;  
                ge[map2det[subevt[i].id]].xpileup[ge[map2det[subevt[i].id]].xmult] = subevt[i].fcode;  
                ge[map2det[subevt[i].id]].xsubevtid[ge[map2det[subevt[i].id]].xmult] = i;  
                        
                ge[map2det[subevt[i].id]].xtheta[ge[map2det[subevt[i].id]].xmult][0] = mapangles[subevt[i].id][0];   
                ge[map2det[subevt[i].id]].xtheta[ge[map2det[subevt[i].id]].xmult][1] = mapanglesi[subevt[i].id][0];               
                ge[map2det[subevt[i].id]].xtheta[ge[map2det[subevt[i].id]].xmult][2] = mapangles1[subevt[i].id][0];               
                ge[map2det[subevt[i].id]].xtheta[ge[map2det[subevt[i].id]].xmult][3] = mapangles2[subevt[i].id][0];               
                ge[map2det[subevt[i].id]].xphi[ge[map2det[subevt[i].id]].xmult][0] = mapangles[subevt[i].id][1];               
                ge[map2det[subevt[i].id]].xphi[ge[map2det[subevt[i].id]].xmult][1] = mapanglesi[subevt[i].id][1];               
                ge[map2det[subevt[i].id]].xphi[ge[map2det[subevt[i].id]].xmult][2] = mapangles1[subevt[i].id][1];               
                ge[map2det[subevt[i].id]].xphi[ge[map2det[subevt[i].id]].xmult][3] = mapangles2[subevt[i].id][1]; 

                ge[map2det[subevt[i].id]].xmult++;   
            }                                                                
            if (map2deti[subevt[i].id] > MAX_GE_XTL && map2deti[subevt[i].id] <= MAX_GE_XTL + MAX_GE_SEG) { //Ge segment
               if ( ge[map2det[subevt[i].id]].smult >= MAX_GE_SEG ) {               
                    printf("SEVERE ERROR: Same Ge(seg) twice within event build; Make event build time smaller!\n");    
                    printf("ge[map2det[subevt[i].id]].xmult=%d, ge[map2det[subevt[i].id]].smult=%d, ge[map2det[subevt[i].id]].bgomult=%d\n", ge[map2det[subevt[i].id]].xmult, ge[map2det[subevt[i].id]].smult,ge[map2det[subevt[i].id]].bgomult);
                    continue;
                    //return -1;
                }   
                ge[map2det[subevt[i].id]].sid[ge[map2det[subevt[i].id]].smult] = map2deti[subevt[i].id];               
                ge[map2det[subevt[i].id]].se[ge[map2det[subevt[i].id]].smult] = subevt[i].energy;               
                ge[map2det[subevt[i].id]].st[ge[map2det[subevt[i].id]].smult] = subevt[i].time;                           
                ge[map2det[subevt[i].id]].sct[ge[map2det[subevt[i].id]].smult] = subevt[i].ctime;  
                ge[map2det[subevt[i].id]].spileup[ge[map2det[subevt[i].id]].smult] = subevt[i].fcode;  
                ge[map2det[subevt[i].id]].ssubevtid[ge[map2det[subevt[i].id]].smult] = i;  
                            
                ge[map2det[subevt[i].id]].smult++;               
            }            
            if (map2deti[subevt[i].id] > MAX_GE_XTL + MAX_GE_SEG) { //BGO 
               if ( ge[map2det[subevt[i].id]].bgomult >= MAX_GE_BGO ) {               
                    printf("SEVERE ERROR: Same Ge(bgo) twice within event build; Make event build time smaller!\n");    
                    printf("ge[map2det[subevt[i].id]].xmult=%d, ge[map2det[subevt[i].id]].smult=%d, ge[map2det[subevt[i].id]].bgomult=%d\n", ge[map2det[subevt[i].id]].xmult, ge[map2det[subevt[i].id]].smult,ge[map2det[subevt[i].id]].bgomult);
                    continue;
                    //return -1;
                }   
                ge[map2det[subevt[i].id]].bgoid[ge[map2det[subevt[i].id]].bgomult] = map2deti[subevt[i].id];               
                ge[map2det[subevt[i].id]].bgoe[ge[map2det[subevt[i].id]].bgomult] = subevt[i].energy;               
                ge[map2det[subevt[i].id]].bgot[ge[map2det[subevt[i].id]].bgomult] = subevt[i].time;                           
                ge[map2det[subevt[i].id]].bgoct[ge[map2det[subevt[i].id]].bgomult] = subevt[i].ctime;              
                ge[map2det[subevt[i].id]].bgopileup[ge[map2det[subevt[i].id]].bgomult] = subevt[i].fcode;
                ge[map2det[subevt[i].id]].bgosubevtid[ge[map2det[subevt[i].id]].bgomult] = i;  

                ge[map2det[subevt[i].id]].bgomult++;               
            } 
                 
                
        } //end G

        ///////////////////
        // S=Si Detector //
        ///////////////////  
        
      } // end i loop over sevtmult     
      ////////////////////////////////////////////
      // END MAP SUB EVENTS INTO DETECTOR TYPES //
      ////////////////////////////////////////////


      ////////////////////////////
      // PROCESS DETECTOR TYPES //
      ////////////////////////////
   

      /////////////////////
      // G = Ge Detector //
      /////////////////////      
      gmult = 0;
      for (i=0; i<MAX_GE; i++) {
        
        max1 = -1; max2 = -1;
        maxid1 = -1; maxid2 = -1;
        
        ge[i].id = i;
        //Addback and Compton Suppression 
	for (j=0; j<ge[i].xmult; j++) {
            
            //compton suppression per crystal 
            if (GE_BGO_SUPPRESSION == TRUE) {
                for (k=0; k<ge[i].bgomult; k++) {
                    tdif = abs( ge[i].xt[j] - ge[i].bgot[k] );
                    if (tdif >= 0 && tdif <= 4096) ge_bgo_tdif[i][tdif]++;
                    if ( (tdif < 50 && ge[i].bgoe[k] > 10) || ge[i].bgopileup[k]==TRUE) { //need to fix bgo pileup with trace analysis 
                        ge[i].xsuppress[j] = TRUE; 
                        ge[i].suppress = TRUE;                        
                    }  
                }
            }      
            
            //addback     
            if (ge[i].xsuppress[j] == FALSE) {       
                if (ge[i].xpileup[j] == TRUE) {
                    ge[i].pileup = TRUE;     
                    continue;
                }              
                //xtl spectra
                if (ge[i].xe[j] > 0 && ge[i].xe[j] < 8192 && ge[i].id >= 1) ge_spe_xtl[(ge[i].id-1)*MAX_GE_XTL + ge[i].xid[j]][ge[i].xe[j]]++;                        
                tdif = abs( ge[i].xt[j] - ge[i].time );
                if (tdif >= 0 && tdif <= 4096 && ge[i].time != 0 ) ge_xtl_tdif[i][tdif]++;           
                                       
                if (ge[i].xe[j] > 50 && ge[i].xe[j] < 4000) {                  
                    if (tdif >= 0 && tdif <= 4096 && ge[i].time != 0 ) ge_xtl_tdif_ethresh[i][tdif]++;           
                    if (tdif < 20 || ge[i].time == 0) {  
                        ge[i].energy = ge[i].energy + ge[i].xe[j];
                        if (max1 < ge[i].xe[j]) {
                            max1 = ge[i].xe[j];
                            maxid1 = j;
                            ge[i].time = ge[i].xt[j];
                            ge[i].ctime = ge[i].xct[j];
                        } 
                    }
                    else {
                        ge[i].nonprompt = TRUE; // the first time will become the adopted value / event   
                    }        
                } 
            }  
           
        }
        
        if (max1 == -1) continue;
        
        //Segmentation Position and Compton Suppression
        for (j=0; j<ge[i].smult; j++) {

            //compton suppression per segment 
            if (GE_BGO_SUPPRESSION == TRUE) {            
                for (k=0; k<ge[i].bgomult; k++) {
                    tdif = abs( ge[i].st[j] - ge[i].bgot[k] );
                    if (tdif < 50 && ge[i].bgoe[k] > 10) {
                        ge[i].ssuppress[j] = TRUE;    
                    }
                }
            }      

            //segment
            if (ge[i].ssuppress[j] == FALSE && ge[i].se[j] > 0 && ge[i].se[j] < 10000) {                    
                if (max2 < ge[i].se[j]) {
                    max2 = ge[i].se[j];
                    maxid2 = j;
                }                                           
            }  

        }    
        
        //Angle assignments
        ge[i].theta[0] = ge[i].xtheta[maxid1][0];                       //detector center
        ge[i].phi[0] = ge[i].xphi[maxid1][0];
        ge[i].theta[1] = ge[i].xtheta[maxid1][1];                       //crystal center
        ge[i].phi[1] = ge[i].xphi[maxid1][1];       
        if (ge[i].sid[maxid2] == 6) {                                   // side channel C  
            ge[i].theta[2] = ge[i].xtheta[maxid1][2];
            ge[i].phi[2] = ge[i].xphi[maxid1][2];  
        }
        else if (ge[i].sid[maxid2] == 5 || ge[i].sid[maxid2] == 7) {    // side channel L/R
            ge[i].theta[2] = ge[i].xtheta[maxid1][3];
            ge[i].phi[2] = ge[i].xphi[maxid1][3];             
        }  
        else {      
            ge[i].theta[2] = ge[i].xtheta[maxid1][1];                   // side channel failure --> crystal center
            ge[i].phi[2] = ge[i].xphi[maxid1][1];  
        }        
        
        
        //clean addback
        if (ge[i].suppress == FALSE && ge[i].pileup == FALSE && ge[i].nonprompt == FALSE) {
            ge[i].clean=TRUE; //maybe clean should not include suppress == FALSE?
        }    
        
        //Ge spectra
        if (ge[i].energy > 0 && ge[i].energy < 8192) ge_spe[ge[i].id][ge[i].energy]++;         
        if (ge[i].energy > 0 && ge[i].energy < 8192 && ge[i].clean == TRUE) ge_spe_clean[ge[i].id][ge[i].energy]++;
        
        //copy data and increment counters
        ge[gmult]=ge[i];
        gmult++;        
        gcount++;
        
      } //end G

     

      ///////////////////////////
      // END PROCESS DETECTORS //
      ///////////////////////////



      ////////////////////////
      // FINAL USER SPECTRA //
      ////////////////////////
        
        
      //gamma-gamma time and energy 
      for (i=0; i<gmult; i++) {
        for (j=0; j<gmult; j++) {
            if (i!=j) {
                tdif = ge[i].time - ge[j].time + 2000;
                //time difference matrix    
                if (tdif >= 0 && tdif < 4096 && ge[i].energy >= 0 && ge[i].energy < 4096 && ge[j].energy >= 0 && ge[j].energy < 4096) {
                    if (ge[i].energy < ge[j].energy)
                    gg_tdif[ge[i].energy][tdif]++;
                    else 
                    gg_tdif[ge[j].energy][tdif]++;
                }   
                //prompt gamma-gamma
                if (tdif >= 1920 && tdif <= 2080) {
                    if (ge[i].energy >= 0 && ge[i].energy < 4096 && ge[j].energy >= 0 && ge[j].energy < 4096) 
                    gg_prompt[ge[i].energy][ge[j].energy]++;
                }                   
                //non-prompt gamma-gamma (mult by 0.2555)
                if ( (tdif >= 1532 && tdif <= 1846) || (tdif >= 2138 && tdif <= 2452 ) ) {
                    if (ge[i].energy >= 0 && ge[i].energy < 4096 && ge[j].energy >= 0 && ge[j].energy < 4096) 
                    gg_nprompt[ge[i].energy][ge[j].energy]++;                
                }    
            }
        }
      }  


      ////////////////////////////
      // END FINAL USER SPECTRA //
      ////////////////////////////

      } //end argc >= 5 condition 

      //event stats, print status every 10000 events
      lle_div=lldiv(evtcount,10000);
      if ( lle_div.rem == 0 && strstr(argv[1], "p") != NULL) {
        fprpos = ftell(fpr);
        tempf = (float)fprsize/(1024.*1024.*1024.);
        printf("Total SubEvents: \x1B[32m%llu \x1B[31m(%d%% pileup)\x1B[0m\nTotal Events: \x1B[32m%llu (%.1f <mult>)\x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\n\033[3A\r", sevtcount, (int)((100*pileupcount)/sevtcount), evtcount, (float)sevtcount/(float)evtcount, (100*fprpos/fprsize), tempf);
      }      

          
  } // end main while loop 
  /////////////////////////
  // END MAIN WHILE LOOP //
  /////////////////////////
  fprpos = ftell(fpr);
  tempf = (float)fprsize/(1024.*1024.*1024.);
  printf("Total SubEvents: \x1B[32m%llu \x1B[31m(%d%% pileup)\x1B[0m\nTotal Events: \x1B[32m%llu (%.1f <mult>)\x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\n\033[3A\r", sevtcount, (int)((100*pileupcount)/sevtcount), evtcount, (float)sevtcount/(float)evtcount, (100*fprpos/fprsize), tempf);
           

  
  
  
  
  
  ////////////////////
  // WRITE SPECTRA //
  ///////////////////
  printf("\n\n\n\nWriting Spectra to Disk ...\n");
     
  //Event Spectra   
  write_data4(HIT, *hit, 4096, 2, overwrite);
  write_data4(MULT, *mult, 4096, 1, overwrite);
  write_data4(TDIFID, *tdifid, MAX_ID, 8192, overwrite);
  write_data4(E_RAW, *e_raw, MAX_ID, 8192, overwrite);
  write_data4(E_CAL, *e_cal, MAX_ID, 8192, overwrite);
  write_data4(TEVT_RAW, *tevt_raw, MAX_ID, 8192, overwrite);
  write_data4(TEVT_CAL, *tevt_cal, MAX_ID, 8192, overwrite);
  write_data4(TCFD_RAW, *tcfd_raw, MAX_ID, 8192, overwrite);
  write_data4(TDIF_RAW, *tdif_raw, MAX_ID, 4096, overwrite);
  write_data4(TDIF_CAL, *tdif_cal, MAX_ID, 4096, overwrite);
 // write_data4(TDIF_CAL0_ETHRESH, *tdif_cal0_ethresh, MAX_ID, 4096, overwrite);
//  write_data4(IDID, *idid, MAX_ID, MAX_ID, overwrite);

  //Detector Processed Spectra
  //Ge
//  write_data4(GE_BGO_TDIF, *ge_bgo_tdif, MAX_GE, 4096, overwrite);
//  write_data4(GE_XTL_TDIF, *ge_xtl_tdif, MAX_GE, 4096, overwrite);
//  write_data4(GE_XTL_TDIF_ETHRESH, *ge_xtl_tdif_ethresh, MAX_GE, 4096, overwrite);
//  write_data4(GE_SPE_XTL, *ge_spe_xtl, MAX_GE*MAX_GE_XTL, 8192, overwrite);
//  write_data4(GE_SPE, *ge_spe, MAX_GE, 8192, overwrite);
//  write_data4(GE_SPE_CLEAN, *ge_spe_clean, MAX_GE, 8192, overwrite);
  //trinity
  write_data4(PID, *pid, 4096, 4096, overwrite);
  write_data4(PID_EVSP, *pid_evsp, 4096, 4096, overwrite);
  write_data4(PID_EVST, *pid_evst, 4096, 4096, overwrite);
  write_data4(PID_EVSTAU, *pid_evstau, 4096, 4096, overwrite);
  write_data4(PID_EVSR, *pid_evsr, 4096, 4096, overwrite);
  write_data4(PID_TAUVSR, *pid_tauvsr, 4096, 4096, overwrite);
  write_data4(GAGGDT, *gaggdt, 4096, 4096, overwrite);

  for (i=0; i<500; i++) {
  	if (strace[i]>0) sptrace[0][i] = strace[i]/dbcount;
  }
  write_data4(SPTRACE, *sptrace, 1, 4096, overwrite);

  //Final User Spectra
  //gamma-gamma
//  write_data4(GG_TDIF, *gg_tdif, 4096, 4096, overwrite);
 // write_data4(GG_PROMPT, *gg_prompt, 4096, 4096, overwrite);
 // write_data4(GG_NPROMPT, *gg_nprompt, 4096, 4096, overwrite);

  
  fclose(fpr);
  fclose(debugfile);
  
  return 0;
}


