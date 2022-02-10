/*********************************************************/
/* PXI Time Order -- J.M. Allmond (ORNL) -- v1 Jul 2016  */
/*                                       -- v2 Feb 2018  */
/*                                       -- v3 Jun 2018  */
/*                                       -- v4 May 2019  */
/* modified for FSU T.L. Ryan Tang (FSU) -- v5 Feb 2022  */ 
/*                                                       */
/* !Time Order Events from Pixie-16 digitizers           */
/* !Max of:                                              */
/* !IDs = static, Evts = dynamic, data = dynamic         */
/*                                                       */
/* gcc -o pxi-time-order pxi-time-order.c                */
/* ./pxi-time-order datafile                             */
/*********************************************************/

/////////////////////////////////////////////////////////
//Code assumes that sequential sub events for a        //
//specific channel are time ordered; therefore,        //
//unmerge data into circular buffers on a per          //
//channel id basis and then remerge channels in        //
//time order.                                          //
/////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#define MAX_CRATES 2
#define MAX_BOARDS_PER_CRATE 13
#define MAX_CHANNELS_PER_BOARD 16
#define BOARD_START 2

#define MAX_ID MAX_CRATES*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD

#define HEADER_LENGTH 4       //unit = words with 4 bytes per word
#define MAX_SUB_LENGTH 2016   //unit = words with 4 bytes per word ; 2016 --> 40 micro-second trace + 4 word header + 12 extra header


#define DEF_SUB_EVENTS 100    //number of events for each dynamic buffer level
#define M1_SUB_EVENTS 1000    //manual input for irregular / non-linear / non-geometric progression 
#define M2_SUB_EVENTS 5000
#define M3_SUB_EVENTS 20000
#define M4_SUB_EVENTS 50000

#define M5_SUB_EVENTS 100000
#define MAX_SUB_EVENTS 200000

#define MAX_MALLOC 4000*1024*1024L //2GB 

#define MAXLONGLONGINT 9223372036854775807 

//TODO change to load file
#include "../mapping.h"

struct subevent
{
  long long int timestamp;
  int length; //unit = words with 4 bytes per word
  int detID;
  unsigned int *data;
};

struct subevent *subevents[MAX_ID]; 
int nevts[MAX_ID], iptr[MAX_ID];    
int maxevts[MAX_ID];               


void swap64(long long int * a, long long int *b){
  long long int t = *a;
  *a = *b;
  *b = t;
}

void swapInt(int * a, int *b){
  int t = *a;
  *a = *b;
  *b = t;
}

int partition(long long int arr[], int index [], int low, int high){
  
  long long int pivot = arr[high];
  int i = (low -1);
  
  for(int j = low; j <= high -1 ; j++){
      if( arr[j] < pivot ){
        i++;
        swap64(&arr[i], &arr[j]);
        swapInt(&index[i], &index[j]);
      }
  }
  swap64(&arr[i+1], &arr[high]);
  swapInt(&index[i+1], &index[high]);
  return i+1;
}

void quickSort(long long int arr[], int index[], int low, int high){
    if (low < high){
        int pi = partition(arr, index, low, high);
        quickSort(arr, index, low, pi - 1);
        quickSort(arr, index, pi + 1, high);
    }
}


int main(int argc, char **argv) {
          
  FILE *fpr, *fpw;
  long int fprsize=0, fprsize_orig=0, fprsize_old=-1, fprpos=0;

  int online = 0;

  unsigned int subhead[HEADER_LENGTH];
  memset(subhead, 0, sizeof(subhead));
  int pause=0;

  long long int nwords=0, evts_tot_read=0, evts_tot_write=0;
  long long int evts_tot_drop = 0;

  long long int time=0,time_old=0;
  int length=0;
  int chn=0; 
  int sln=0;
  int crn=0;
  int id=0;
  
  int idmax=0;
  int totmem=0;
  int outoforder=0;
  int evts_old=0;
  int evts_new=0;

  long long int timemin=0, timemin_old=0;
  int min_id = -1;
 
  memset(nevts, 0, sizeof(nevts));
  memset(iptr, 0, sizeof(iptr));     /// index of time
  
  long long int timeIndex[MAX_ID];
  int index[MAX_ID];
  int fillSize = 0;

  int i=0, j=0;

  //open input event file
  if ((fpr = fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "Error, cannot open input file %s\n", argv[1]);
    return 1;
  }

  //write time order file to current location, not location of event file
  char filenameto[80];
  char *filename = strrchr(argv[1], '/');
  if (filename == NULL) strcpy(filenameto,argv[1]);
  else strcpy(filenameto,filename+1);
  strcat(filenameto,".to");
  if ((fpw = fopen(filenameto, "w")) == NULL) {
    fprintf(stderr, "Error, cannot open output file %s\n", filenameto);
    return 1;
  }
  
  printf("open : \033[1;31m%s\033[m\n", argv[1]);
  printf("save : \033[1;34m%s\033[m\n", filenameto);
  
  int eventWindow = 0;
  
  if( argc >= 3 ) eventWindow = atoi(argv[2]);
  printf(" event build time window : %d ticks = %d ns \n", eventWindow, 10* eventWindow);
  

  //check for lockfile, active PID, and event file for auto "online" mode detection
  FILE *FPLOCK;
  char lockfile[1024];    
  strcpy(lockfile, getenv("HOME"));    
  strcat(lockfile, "/.Pixie16Lock");
  int lockpid;

  FILE *FPPATH;
  char pathfile[1024];
  char line[1024];
  char onlinefile[1024];
  strcpy(pathfile, getenv("HOME"));    
  strcat(pathfile, "/.Pixie16Path");  

  FPLOCK = fopen(lockfile, "r");
  if (FPLOCK != NULL) {
    fscanf(FPLOCK, "%d", &lockpid);
    fclose(FPLOCK);  
    //PID from lockfile matches a running PID; run timesort in "online" mode for now
    if (getpgid(lockpid) >= 0) {
      FPPATH = fopen(pathfile, "r");
      if (FPPATH == NULL) {
        online = 0;
      }else {
        fgets(line, 1024, FPPATH); //skip first line
        fgets(line, 1024, FPPATH); //need second line 
        sscanf(line,"%s\n", onlinefile);
        fclose(FPPATH);
        if (filename == NULL) {	
          if (strcmp(onlinefile,argv[1]) == 0)  online = 1;
        }else {
          if (strcmp(onlinefile,filename+1) == 0) online = 1;
        }
      }
    }
  }
  if (online == 1) printf("Auto Mode: \x1B[32mOnline\x1B[0m\n");
  else printf("Auto Mode: \x1B[32mOffline\x1B[0m\n");

  //check file size for auto "online" mode 
  fprpos = ftell(fpr);
  fseek(fpr, 0L, SEEK_END);
  fprsize = fprsize_orig = ftell(fpr);
  fseek(fpr, fprpos, SEEK_SET);

  //Get memory for default number of subevents per channel id
  for (i=0; i < MAX_ID; i++){
    subevents[i] = (struct subevent *) malloc(sizeof(struct subevent)*DEF_SUB_EVENTS);
    if (subevents[i] == NULL) {
      printf("malloc failed\n");
      return -1;
    }
    totmem += sizeof(struct subevent)*DEF_SUB_EVENTS;
    maxevts[i] = DEF_SUB_EVENTS;
    for (j=0; j<DEF_SUB_EVENTS; j++) {
      subevents[i][j].data = NULL;
      subevents[i][j].length = 0;
      subevents[i][j].timestamp = 0;      
      subevents[i][j].detID = -1;      
    }
  }
  
  int count = 0;
  int debugCount = 0;

  printf("Static Memory = %ld KB (cf. MAX_ID=%d)\n", sizeof(subevents)/1024, MAX_ID);
  while (1) { //main while loop

    /////////
    while (1) { //fill buffers until (A) maxevents or (maxevents and 2GB) is reached for any ID 
                //                   (B) EOF
                //                   (C) auto online mode will wait for updates and break out of fill buffers for narrow time window
                //                   read 4-byte header
      if (pause == 1) {
        pause = 0;   
      }else {

        //////////////
        //auto online
        while ( (fprsize - nwords*sizeof(int) < MAX_SUB_LENGTH*sizeof(int)) && online == 1) {
                    
          online = 0;
          usleep(100000); //wait 0.1 seconds before checking (prevents excessive cpu usage)
                    
          //check new file size 
          fprpos = ftell(fpr);
          fseek(fpr, 0L, SEEK_END);
          fprsize = ftell(fpr);
          fseek(fpr, fprpos, SEEK_SET);

          //check for lock file and active PID
          FPLOCK = fopen(lockfile, "r");
          if (FPLOCK != NULL) {
            fscanf(FPLOCK, "%d", &lockpid);
            fclose(FPLOCK);  
            if (getpgid(lockpid) >= 0) {
              FPPATH = fopen(pathfile, "r");
              if (FPPATH != NULL) {
                fgets(line, 1024, FPPATH); //skip first line
                fgets(line, 1024, FPPATH); //need second line 
                sscanf(line,"%s\n", onlinefile);
                fclose(FPPATH);	
                if (filename == NULL) {	
                  if (strcmp(onlinefile,argv[1]) == 0) online = 1;
                }else {
                  if (strcmp(onlinefile,filename+1) == 0) online = 1;
                }
              }
            }
          }
        } //end auto online


        //read 4-byte header
        if (fread(subhead, sizeof(subhead), 1, fpr) != 1) break; 
        nwords = nwords + HEADER_LENGTH;
        chn = subhead[0] & 0xF;
        sln = (subhead[0] & 0xF0) >> 4;
        crn = (subhead[0] & 0xF00) >> 8;
        id = crn*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD + (sln-BOARD_START)*MAX_CHANNELS_PER_BOARD + chn;   
        length = (subhead[0] & 0x7FFE0000) >> 17; //unit = words with 4 bytes per word
        time = ( (long long int)(subhead[2] & 0xFFFF) << 32) + subhead[1]; 
        if (id > idmax) idmax = id;

      }

      //check memory
      if (totmem > MAX_MALLOC) {printf("Error: Exceeded MAX_MALLOC"); return -1;}
        
      //Expand memory for more events (careful when final is to left of initial in circular buffer)
      if ( maxevts[id] - nevts[id] == 1 && totmem < MAX_MALLOC) {
                  
        if (maxevts[id] == DEF_SUB_EVENTS) {evts_old = DEF_SUB_EVENTS; evts_new = M1_SUB_EVENTS;}
        if (maxevts[id] == M1_SUB_EVENTS)  {evts_old = M1_SUB_EVENTS;  evts_new = M2_SUB_EVENTS;}
        if (maxevts[id] == M2_SUB_EVENTS)  {evts_old = M2_SUB_EVENTS;  evts_new = M3_SUB_EVENTS;}
        if (maxevts[id] == M3_SUB_EVENTS)  {evts_old = M3_SUB_EVENTS;  evts_new = M4_SUB_EVENTS;}
        if (maxevts[id] == M4_SUB_EVENTS)  {evts_old = M4_SUB_EVENTS;  evts_new = M5_SUB_EVENTS;}
        if (maxevts[id] == M5_SUB_EVENTS)  {evts_old = M5_SUB_EVENTS;  evts_new = MAX_SUB_EVENTS;}

        if (maxevts[id]==evts_old && totmem + (evts_new-evts_old)*(sizeof(struct subevent) + sizeof(unsigned int)*length) < MAX_MALLOC) {
          subevents[id] = (struct subevent *) realloc(subevents[id], sizeof(struct subevent)*evts_new);  
          if (subevents[id] == NULL) {
            printf("realloc failed\n");
            return -1;
          }
          totmem = totmem - sizeof(struct subevent)*evts_old + sizeof(struct subevent)*evts_new;
          maxevts[id] = evts_new;
          for (j=evts_old; j<evts_new; j++) { 
            subevents[id][j].data = NULL;
            subevents[id][j].length = 0;
            subevents[id][j].timestamp = 0; 
            subevents[id][j].detID = -1; 
          }
                   
          // if circular buffer is wrapped around (i.e., final is to left of intial, move data to right of initial)
          if (iptr[id] + nevts[id] > evts_old) {
            for (j=0; j<iptr[id] + nevts[id] - evts_old; j++) {
              if (subevents[id][evts_old+j].data == NULL) {
                subevents[id][evts_old+j].data = (unsigned int *)  malloc(sizeof(unsigned int)*subevents[id][j].length);
                if (subevents[id][evts_old+j].data == NULL) {
                  printf("malloc failed\n");
                  return -1;
                }
                totmem += sizeof(unsigned int)*subevents[id][j].length;
              }    
              subevents[id][evts_old+j].length = subevents[id][j].length;
              subevents[id][evts_old+j].timestamp = subevents[id][j].timestamp;
              subevents[id][evts_old+j].detID = subevents[id][j].detID;
              for (i=0; i<subevents[id][evts_old+j].length; i++) {
                subevents[id][evts_old+j].data[i]=subevents[id][j].data[i];
              }
              //free data memory until it's needed again
              free(subevents[id][j].data);
              subevents[id][j].data = NULL;
              totmem -= sizeof(unsigned int)*subevents[id][j].length; 
              subevents[id][j].length = 0;
              subevents[id][j].timestamp = 0;
              subevents[id][j].detID = -1;
            }
          }

        }  

      }  


      // time control of buffer filling for auto online mode (reset if initial value or large gap > 3.5 sec)
      //                                                      large gap could be from low rate or un/replug 
      if ( time_old == 0 || (time - time_old)/10000000 > 35 ) time_old = time; 
        
      //fill buffers until full (online mode will stop filling buffers after 2.5 sec lag betweeen output/input)
      if ( nevts[id] < maxevts[id] && ( (time - time_old)/10000000 < 25 || online == 0 ) ) {        
        j = nevts[id] + iptr[id];
        if (j >= maxevts[id]) j -= maxevts[id];
                  
        subevents[id][j].timestamp = time;
        subevents[id][j].detID = mapping[id];

        if (subevents[id][j].data == NULL) {
          subevents[id][j].data = (unsigned int *)  malloc(sizeof(unsigned int)*length);
          if (subevents[id][j].data == NULL) {
            printf("malloc failed\n");
            return -1;
          }
          totmem += sizeof(unsigned int)*length;
        }else if (length != subevents[id][j].length) { //not needed anymore since always free data after use now. Keep for future ...
          subevents[id][j].data = (unsigned int *)  realloc(subevents[id][j].data, sizeof(unsigned int)*length);
          if (subevents[id][j].data == NULL) {
            printf("realloc failed\n");
            return -1;
          }
          totmem = totmem - sizeof(unsigned int)*subevents[id][j].length + sizeof(unsigned int)*length;
        }

        subevents[id][j].length = length;

        if (length>HEADER_LENGTH) {
          if (fread(subevents[id][j].data + HEADER_LENGTH, (length-HEADER_LENGTH)*sizeof(int), 1, fpr) != 1) break;
          nwords = nwords + (length-HEADER_LENGTH);
        }
                  
        for (i=0; i < HEADER_LENGTH; i++) {
          subevents[id][j].data[i] = subhead[i];  
        }             

        nevts[id]++;
        evts_tot_read++;

      }else {
        pause = 1;
        break;
      }

    } // end while for fill buffers  maxevts[id]  
    /////////            
    
    
    //######################## FSU 
    // find group of event within timewindow, check is contain gamma. if not, throw away all. 
    
    if( eventWindow > 0 ){
    
      //quick sort of subevents[i][iptr[i]].timestamp, i = 0 , idmax +1
      for( i = 0 ; i < idmax + 1; i++) {
        
        if( nevts[i] > 0 ){
          timeIndex[i] = subevents[i][iptr[i]].timestamp; 
        }else{
          timeIndex[i] = MAXLONGLONGINT; 
        }
        index[i] = i;
        //printf("%3d , %llu, %d \n", i, timeIndex[i], index[i] );
      }
      quickSort(timeIndex, index, 0, idmax);
      for( i = 0 ; i < idmax + 1; i++) {
        //printf("%3d , %llu , %d\n", i, timeIndex[i], index[i]);
      }
      
      //reduce the index size
      fillSize = 1;
      for( i = 1; i < idmax + 1; i++){
        if( timeIndex[i] - timeIndex[0] < eventWindow) fillSize ++;
      }
      
      
      // display
      if ( count < debugCount) {
      //if ( count < debugCount || timeIndex[10] == MAXLONGLONGINT) {
        printf("===============================================\n");
        
        for( int i = 0; i < idmax+1; i++) {
          printf("%3d, %llu, %d \n", i, timeIndex[i], index[i]);
          if( i == fillSize - 1 ) printf("------------------- %d \n", fillSize);
        }

      }
      
      if( timeIndex[0] == MAXLONGLONGINT ) break;
      
      //CHeck if fill evt.to for the data.
      bool fillFlag = false;
      
      for( i = 0 ; i < fillSize; i++ ){
        if( count < debugCount ) printf("********************* %llu , detID : %d\n", timeIndex[0], subevents[index[i]][iptr[index[i]]].detID);
        if( subevents[index[i]][iptr[index[i]]].detID < 100 ) fillFlag = true;
      }
      if( count < debugCount ) printf("=============== fillFlag : %d \n", fillFlag);

      for( i = 0; i < fillSize; i++){
        min_id = index[i];

        if( fillFlag ){        
          fwrite(subevents[min_id][iptr[min_id]].data, sizeof(unsigned int)*subevents[min_id][iptr[min_id]].length, 1, fpw);  
          
          if( count < debugCount ) printf("filling  subevents[%d][iprt[%d]] \n", min_id, min_id);
          evts_tot_write++;  
        }else{
          evts_tot_drop++;
        }
        
        
        //free data memory up until it's needed again
        if( count < debugCount ) printf("    Free  subevents[%d][iprt[%d]] \n", min_id, min_id);
        
        free(subevents[min_id][iptr[min_id]].data); 
        subevents[min_id][iptr[min_id]].data = NULL; 
        totmem -= sizeof(unsigned int)*subevents[min_id][iptr[min_id]].length; 
        subevents[min_id][iptr[min_id]].length = 0; 
        subevents[min_id][iptr[min_id]].timestamp = 0;
        subevents[min_id][iptr[min_id]].detID = -1;
          
        nevts[min_id]--;  
        if (++iptr[min_id] >= maxevts[min_id]) iptr[min_id] -= maxevts[min_id];
            
      }
    
      count ++;
    
    }else{
      
      /////////      
      // write event with minimum time to file
      timemin_old = timemin; 
      timemin = -1;          
      for (i=0; i < idmax + 1; i++) {   //could be MAX_ID but limit ourselves to current max, idmax
        if (nevts[i] > 0) {
          if (timemin == -1) {
            timemin = subevents[i][iptr[i]].timestamp; 
            time_old = timemin;
            min_id = i;   
          }else if (subevents[i][iptr[i]].timestamp < timemin) {
            timemin = subevents[i][iptr[i]].timestamp;  
            time_old = timemin;
            min_id = i;  
          }    
        }
      } 

      if (timemin > -1) {  
        if (timemin < timemin_old) {
          printf("\nWarning!!! timemin = %lld and timemin_old = %lld and min_id = %d\n", timemin, timemin_old, min_id); 
          outoforder++;
        }  
        if (subevents[min_id][iptr[min_id]].data == NULL) {printf("Error: data = NULL\n"); return -1;}
        
        fwrite(subevents[min_id][iptr[min_id]].data, sizeof(unsigned int)*subevents[min_id][iptr[min_id]].length, 1, fpw);  

        //free data memory up until it's needed again
        free(subevents[min_id][iptr[min_id]].data); 
        subevents[min_id][iptr[min_id]].data = NULL; 
        totmem -= sizeof(unsigned int)*subevents[min_id][iptr[min_id]].length; 
        subevents[min_id][iptr[min_id]].length = 0; 
        subevents[min_id][iptr[min_id]].timestamp = 0;
          
        nevts[min_id]--;  
        if (++iptr[min_id] >= maxevts[min_id]) iptr[min_id] -= maxevts[min_id];
        evts_tot_write++;  
        
      }else break;
      /////////
      
    }

    //print statistics
    if( evts_tot_read % 10000 == 0 )
      printf("Malloc (%d MB) : evts in (\x1B[34m%lld\x1B[0m) : evts out (\x1B[32m%lld\x1B[0m) : evts drop (\x1B[32m%lld\x1B[0m) : diff (\x1B[31m%lld\x1B[0m)\r", 
        (totmem)/1024/1024, evts_tot_read, evts_tot_write, evts_tot_drop, evts_tot_read-evts_tot_write - evts_tot_drop); 
        
  } //end main while 
  

  //cleanup 
  fclose(fpr);
  fclose(fpw);
  for (i=0; i<MAX_ID; i++){
    free(subevents[i]);
    totmem -= sizeof(struct subevent)*maxevts[i];
  }
  
  //print statistics last time
  printf("\33[2K");
  printf("Malloc (%d MB) : evts in (\x1B[34m%lld\x1B[0m) : evts out (\x1B[32m%lld\x1B[0m) : diff (\x1B[31m%lld\x1B[0m)\n", (totmem)/1024/1024, evts_tot_read, evts_tot_write, evts_tot_read-evts_tot_write); 
  if (outoforder > 0) printf("\x1B[31mWarning, there are %d events out of time order\x1B[0m\n", outoforder);
  if (totmem != 0) printf("\x1B[31mError: total memory not conserved\x1B[0m\n");


  return 0;
}


