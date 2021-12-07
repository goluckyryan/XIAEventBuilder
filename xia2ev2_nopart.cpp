#include <stdio.h>
#include <iostream>
#include <fstream>
//for 64bit compiler:
#include <stdlib.h>
#include <string.h>

char outdata[70],outloga[70]; 

//int *outenergy=NULL;
int outsize=1;
//char *outcount=NULL;
//char *outdetnum=NULL;
char *outcount=new char[1];
char *outdetnum=new char[1];
//outenergy=new int[1];
int *outenergy=new int[1];
//outcount[0]=1;
int correctflag = 0,gateflag = 0,CntAll = 0,CntGate = 0;
long maxtime=100; //max number of 10 ns tics for an event
long maxeventime=0,mineventime=400000000,CntEvts[64];
int multiplicity[10],startdet=15,PileUp[64];
int evN=0;
int evStart = 0;
int evEnd = 1;
unsigned int *bufsam=NULL;
long *bufsiz=NULL;
long negativedelevtimecount=0,timezero=-30000;
int size1,size2,bsize,bsam=2,bread=1,chan,slot,chann,outevtcount=0;
int pu,energy,evtimehi,pileupcount=0,zerocount=0;
unsigned short cfd,cfdold;  
int detectorcount=0,mul=0;
short detenergy[256],detid[256];
unsigned long maskpu = 2147483648,evtimelo;
long long deltime=0LL,timeinevent=0LL;
unsigned long long evtime,multiplier=4294967296LL,oldtime=0LL;
long long zero=0LL;
int channel[65][100],timelo[65],timehi[65];
float correctable[65][100];
float fracx,correct,top,bottom,tdc,percent,averagetime=0;
int ipnt,iadc,itdc;
int pntcnt[65];
  
void InitializeVariables(char **args);
void WriteOutEvent(char **args,FILE *outdat);
void BuildEvent();

int main(int argn,char **argv) {
  if ((argn<3) || (argn>8)) {
    std::cerr<<"Usage:\n "<<argv[0]<<" filein.evt fileout correctfilename startdetnumber maxtime start_buf# stop_buf#\n  Converts physics buffers (type 30) from  \'file.evt\' to \'fileout.ev2\' from \'start_buf#\' until\n \'stop_buf#'  or whole file if no event numbers are given or \'start buf#\' = 0\n  startdetnum is the number of the dE ADC on a scale of 1 to n.\n  maxtime defaults to 100 if not specified\n  **Use \'NULL\' for the output file to avoid writing an ev2 file. \nNote that the start and end buffer numbers are for the input where each adc \nread forms another buffer.\n This version sets the time in event of the start detector to 100 and all other times relative to that. If a correct file is specified, it dorrects the detectors in the file for walk.  It does not write out events which do not contain the start detector.\nAlso it only writes out events in which every ADC times in within the limits specified in the correct.tab file."<<std::endl;
    exit(0);
  }

  FILE *infile=fopen(argv[1],"r");
  if (infile==NULL) {
    std::cerr<<"Problem opening "<<argv[1]<<std::endl;
    exit(0);
  }

  strcpy (outdata,argv[2]);
  strcat (outdata,".ev2");
  FILE *outdat=fopen(outdata,"wb");
  std::cout<<"Writing output data to "<<outdata<<std::endl; 


/**
  if (argn==3) {  // Only in and out files given
    evStart=0;
    evEnd=1;
  }
 if (argn == 4) {   // In, out, correct files given
    evStart=0;
    evEnd=1;
    correctflag = 1;
    }

 if (argn==5) {   // In, out, correct files and startdet given
    evStart=0;
    evEnd=1;
    startdet=atoi(argv[4]);
    correctflag = 1;
    }

 if (argn==6) {   // Also maxtime given
    evStart=0;
    evEnd=1;
    maxtime=atoi(argv[5]);
    startdet=atoi(argv[4]);
    correctflag = 1;
    }

  if (argn>6) {
    evStart=atoi(argv[6]);
    maxtime=atoi(argv[4]);
    startdet=atoi(argv[5]);
    correctflag = 1;
    if (evStart == 0) evEnd = 1;
    else {
      std::cerr<<"You must specify an end event number if the start event number is not 0.  Your start event number is "<<evStart<<std::endl;
      exit(0);
    }
    
  } else if (argn==7) {
    
    evStart=atoi(argv[5]);
    evEnd=atoi(argv[6]);
    maxtime=atoi(argv[4]);
    startdet=atoi(argv[3]);
  }
  * 
*/
  int sizes=0;
  long typebuf[100]={0},type;

  if (evEnd<evStart) {
    std::cerr<<"start event number must be less than last)."<<std::endl;
    fclose(infile);
    exit(0);
  }
// change this for 64bit compiler long *bufsam=NULL;

  std::cout<<"The start detector number is "<<startdet<<std::endl;

  InitializeVariables(argv);

  while ((evN<evEnd) || (evStart == 0)) {
    //std::cout<<evStart<<"  end"<<evEnd<<"  num"<<evN<<std::endl;

    bufsiz=new long[bsam];
    //change long -> int for 64 bit
    fread(bufsiz,sizeof(int),bread,infile);
    //std::cout<<bread<<"  "<<bufsiz[0]<<std::endl;

    sizes=bufsiz[0];
    if (feof(infile)) break;
    evN++;

    bsize = (sizes-4);
    // change for 64bit  bufsam=new long[bsize/4+1];
    bufsam=new unsigned int[bsize/4+1];
    fread((char*)bufsam,1,bsize,infile);

    //std::cout<<"("<<evN<<")  "<<bufsam[0]<<"  "<<bufsam[1]<<"  "<<bufsam[2]<<"  "<<bufsam[5]<<std::endl;

    //std::cout<<"evN "<<evN<<"  evStart "<<evStart<<"  evEnd "<<evEnd<<std::endl;
    if (evN>=evStart) {

      if (bufsam[0] == 30) {
        chan = (bufsam[2]) & (15);
        slot = ((bufsam[2]) & (240))/16;
        chann = (slot - 2)*16 + chan + 1;
        pu = ((bufsam[2]) & (maskpu))/maskpu;
        energy = ((bufsam[5]) & 65535);
        evtimehi = ((bufsam[4]) & 65535);
        evtimelo = bufsam[3];
        evtime = evtimelo + multiplier*evtimehi;
        cfd = bufsam[4]/65536;
        deltime = evtime - oldtime;
        if(deltime < zero) negativedelevtimecount++;
        if (energy==0) zerocount++;
        if(pu>0)pileupcount++;
        if ((pu > 0)&&(chann>0)&&(chann<65)) PileUp[chann-1]++;


        //	ignore pileups
            timeinevent=timeinevent+deltime;
            oldtime = evtime;

        if(timeinevent > maxtime) {
          //This detector starts another event.  First write out the previous event. 

          WriteOutEvent(argv,outdat);

          // Now store this signal which starts a new event
          timeinevent = 1;
          timezero = -30000;
          detectorcount = 0;
          mul = 0;
          averagetime = 0;

          BuildEvent();
        }else {
          //This is just another detector in the current event
          //Now we have to put the event just read into the buffer

          BuildEvent();

        } //end if(timeinevent

      } //end if bufsam[0]=30

    } //end if evN >= evStart
   
  } //end while evN < evEnd

  delete [] bufsiz;
  delete [] bufsam;

  if(feof(infile) && evN<evEnd && evStart != 0) {
    std::cerr<<"The requested event is out of range (too high)."<<std::endl;
  }

    //Write out last event

    WriteOutEvent(argv,outdat);

 //Print out summaries	
	std::cout<<negativedelevtimecount<<" negative delta event times;    "<<pileupcount<<"  pileups;   "<<zerocount<<" zeros"<<std::endl; 
  std::cout<<" max time in event = "<<maxeventime<<";   min time betweeen events = "<<mineventime<<"\n If either of these numbers is close to "<<maxtime<<" you may want to change maxtime"<<std::endl;       
  if (strcmp(argv[2],"NULL") != 0) std::cout<<"wrote "<<outevtcount-1<<" events to file "<<argv[2]<<std::endl;

  for(int is=0;is<10;is=is+2) {
    std::cout<<multiplicity[is]<<" mul "<<is+1<<"    "<<multiplicity[is+1]<<" mul "<<is+2<<std::endl;
  }
  
  for(int is=0; is<64;is++) {
    if (PileUp[is] > 0) std::cout<<"PU("<<is+1<<") = "<<PileUp[is]<<std::endl;
  }

  std::cout<<"Total coincidence events ="<<CntAll<<"  Number in time gates = "<<CntGate<<std::endl;
  percent = 0.;
  if (CntAll > 0) percent = 100. * float(CntGate)/ float(CntAll);
  std::cout<<"  Percent in time gate = "<<percent<<std::endl;
  
  fclose(infile);

//  Write summaries to file:

  strcpy (outloga,argv[2]);
  strcat (outloga,".log");
  FILE *outlog=fopen(outloga,"wb");
  std::cout<<"Writing conversion log to "<<outloga<<std::endl; 

  if (correctflag >0) fprintf(outlog,"Corrected time walk using file %s \n\n",argv[3]);
  if (strcmp(argv[2],"NULL") != 0) fprintf(outlog,"wrote %i events to %s using xia2ev2_part with particle detector number %i\n\n",outevtcount,outdata,startdet);

  fprintf(outlog,"negativedelevtimecount = %i, pileupcount = %i  zerocount = %i\n",negativedelevtimecount,pileupcount,zerocount); 

  fprintf(outlog,"max time (10 ns) in event = %i  min time betweeen events = %i\n\n",mineventime,maxtime); 
      

  for(int is=0;is<10;is=is+2) fprintf(outlog,"multiplicity %i: counts = %i;    multiplicity %i: counts = %i\n",is+1,multiplicity[is],is+2,multiplicity[is+1]);
  fprintf(outlog,"\n");
  
  for(int is=0; is<64;is++) {
    if (CntEvts[is] > 0) fprintf(outlog,"ADC(%i): Counts = %i  Pileups = %i\n",is+1,CntEvts[is],PileUp[is]);
  }

  percent = 0.;
  if (CntAll > 0) percent = 100. * float(CntGate)/ float(CntAll);
  fprintf(outlog,"Total coincidence events = %i   Number in time gates = %i\n   Percent in time gate = %5.2f\n",CntAll,CntGate,percent);


  if (strcmp(argv[2],"NULL") != 0) fclose(outdat);
  return 0;
}



void InitializeVariables(char **args) {

 //initialize correctable
 
  for (ipnt = 0; ipnt<65; ipnt++) pntcnt[ipnt]=0;
  for (int adc=0; adc<65; adc++){
    timelo[adc] = 0;  timehi[adc]=32767;
    for (ipnt=0; ipnt<100; ipnt++) {
      correctable[adc][ipnt] = -1.;
      channel[adc][ipnt]=-1;
    }
  }

  for (int is=0; is<10; is++) multiplicity[is]=0;
  for (int is=0; is<64;is++) PileUp[is]=0;
  for (int is=0; is<64; is++) CntEvts[is]=0; 


  if (correctflag > 0) {
// Read in correction table 
    std::cerr<<"correction table file = "<<args[3]<<std::endl;
    std::ifstream *in = new std::ifstream( args[3], std::ios::in);
    ipnt=0;
    float offset;
    *in >> iadc;
    *in >> offset;
    *in >> timelo[iadc];
    *in >> timehi[iadc];
    while (iadc >= 0) {
      if (startdet == -2) std::cout<<"iadc ="<<iadc<<"   offset = "<<offset<<"  timelo = " <<timelo[iadc]<<"  timehi = "<<timehi[iadc]<<std::endl;
      if (iadc > 0 && iadc < 65) {
        *in >> channel[iadc][ipnt]; 
        while (channel[iadc][ipnt] >=0) {
          *in >> correctable[iadc][ipnt];
	  correctable[iadc][ipnt] = correctable[iadc][ipnt] +offset;
	  pntcnt[iadc]++;
	  if (startdet == -2) std::cout<<"["<<iadc<<"]  "<<channel[iadc][ipnt]<<": "
	     <<correctable[iadc][ipnt]<<"\n";
          ipnt++;
          *in >> channel[iadc][ipnt];
//	      std::cout<<channel[iadc][ipnt]<<"\n";
          }
//        std::cout<<"pntcnt["<<iadc<<"] = "<<pntcnt[iadc]<<"\n";
       *in >> iadc;
       *in >> offset;
       *in >> timelo[iadc];
       *in >> timehi[iadc];

     } //end if (iadc
     ipnt=0;
   } //end while (iadc


     if (startdet == -2) return;

  } //end if (correctflag > 0

} //end void InitializaVariables

void WriteOutEvent(char **args,FILE *outdat) {

  if(timeinevent < mineventime) mineventime = timeinevent;
  if (mul>8) mul=9;
  multiplicity[mul]++;
   
if (strcmp(args[2],"NULL") != 0) {
    detectorcount = detectorcount & 255;
//  timezero will be > 0 if startdet fired in this event
//    if(detectorcount > 0 && timezero > -30000) {
    if(detectorcount > 0 && timezero == -30000) {
      CntAll ++;
//      printf("\n[%i] ",detectorcount);
      for (int icnt=1; icnt<detectorcount; icnt = icnt + 2) {
        iadc = detid[icnt-1];
//        printf("%i: %i %i   ",iadc,detenergy[icnt-1],detenergy[icnt]);
        if (pntcnt[iadc] > 0) {  //Is there a correct table for this ADC?
	  ipnt=0;

//        Linear interpolation between points in correct.tab
	  while (ipnt < pntcnt[iadc]  && (detenergy[icnt-1] >= channel[iadc][ipnt])) ipnt++;
          top = correctable[iadc][ipnt] - correctable[iadc][ipnt-1];
	  bottom = channel[iadc][ipnt] - channel[iadc][ipnt-1]; 
          fracx = top/bottom;
	  correct = correctable[iadc][ipnt-1] + fracx*(detenergy[icnt] - channel[iadc][ipnt-1]); 
          tdc = detenergy[icnt];
	  if (detid[icnt-1] != startdet) tdc = 200 - correct + tdc;
	  if (detid[icnt-1] == startdet) tdc = correct + tdc;
          if (tdc <1) tdc=1;  if (tdc > 32767) tdc=32767;
          detenergy[icnt] = tdc + 0.5;
          averagetime += tdc;
 	  if (detid[icnt-1] == startdet) timezero = (long) (tdc+0.5);
        }  //end if (pntcnt

        else if (timezero == -30000) {
          detenergy[icnt] = detenergy[icnt] + 99;
          averagetime += detenergy[icnt];
        }


      } //end for (int icnt=

//  Adjust for timezero
      averagetime = 2*averagetime/detectorcount;
      for (int icnt=1; icnt<detectorcount; icnt = icnt + 2) {
        tdc = detenergy[icnt];
//        tdc = tdc + 99 - timezero;
        tdc = tdc + 99 - averagetime;
        detenergy[icnt] = tdc;
      }
        

//  Check time gates
      gateflag = 0;
      for (int icnt=1; icnt<detectorcount; icnt = icnt + 2) {
        itdc = detenergy[icnt]; iadc = detid[icnt-1];
        if ((itdc < timelo[iadc])  || (itdc > timehi[iadc])) gateflag = 1;
      }

// write out event if in gate

      if (gateflag == 0) {
        CntGate ++;
        outevtcount++;
        outcount[0] = detectorcount;
	fwrite(outcount,sizeof(char),outsize,outdat);
        for (int icnt=0; icnt<detectorcount; icnt++) {
	  outdetnum[0]=detid[icnt];
          outenergy[0]=detenergy[icnt];
	  fwrite(outdetnum,sizeof(char),outsize,outdat);
	  fwrite(outenergy,sizeof(short),outsize,outdat);
        } // end for (int icnt=0
	fwrite(outcount,sizeof(char),outsize,outdat);
      } // end if (gateflag == 0)
	      
    } //end if(detectorcount > 0 && timezero > -30000)
  } //endif strcmp
//  timeinevent = 0;
//  timezero=-30000;

} //end void WriteOutEvent(

void BuildEvent() {

  if(timeinevent > maxeventime) maxeventime = timeinevent;
  mul++;
//first put energy in list
//suppress zero energies and those over range to not confuse gnuscope
  if(energy>10 && energy < 32000) {
    detid[detectorcount]=chann;
    detenergy[detectorcount]=energy;
    detectorcount++;
    if (chann == startdet) timezero = timeinevent + 99; 
//now put timeinevent in list as chann + 64
    detid[detectorcount]=chann+64;
    detenergy[detectorcount]=timeinevent;
    detectorcount++;
  } // end if (energy>10 ..

}  // end BuildEvent
