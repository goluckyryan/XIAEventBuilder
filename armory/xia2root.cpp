#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include <vector>

#define NUMDET 64  /// number of detector
#define STARTDETID 15

std::vector<std::string> SplitStr(std::string tempLine, std::string splitter, int shift = 0){

  std::vector<std::string> output;

  size_t pos;
  do{
    pos = tempLine.find(splitter); /// fine splitter
    if( pos == 0 ){ ///check if it is splitter again
      tempLine = tempLine.substr(pos+1);
      continue;
    }

    std::string secStr;
    if( pos == std::string::npos ){
      secStr = tempLine;
    }else{
      secStr = tempLine.substr(0, pos+shift);
      tempLine = tempLine.substr(pos+shift);
    }

    ///check if secStr is begin with space
    while( secStr.substr(0, 1) == " "){
      secStr = secStr.substr(1);
    };
    
    ///check if secStr is end with space
    while( secStr.back() == ' '){
      secStr = secStr.substr(0, secStr.size()-1);
    }

    output.push_back(secStr);
    //printf(" |%s---\n", secStr.c_str());
    
  }while(pos != std::string::npos );

  return output;
}

std::vector<std::vector<double>> LoadCorrectionParameters(TString corrFile){

  printf("==================== load correction parameters : %s", corrFile.Data());
  std::ifstream file;
  file.open(corrFile.Data());
  
  std::vector<std::vector<double>> corr;
  corr.clear();
  
  std::vector<double> detCorr;
  detCorr.clear();
  
  if( file.is_open() ){  
    
    while( file.good() ){
    
      std::string line;
      getline(file, line);
    
      if( line.substr(0,1) == "#" ) continue;
      if( line.substr(0,2) == "//" ) continue;
      if( line.size() == 0 ) continue;
      
      std::vector<std::string> temp = SplitStr(line, " ");
      
      detCorr.clear();
      for( int i = 0; i < (int) temp.size() ; i++){
        detCorr.push_back(std::stod(temp[i]));
      }
      corr.push_back(detCorr);
    }
    
    file.close();
    
    printf(".... done\n");
    printf("===== correction parameters \n");
    for( int i = 0; i < (int) corr.size(); i++){
      printf("det : %2d | ", i );
      int len = (int) corr[i].size();
      for( int j = 0; j < len - 1 ; j++){
        printf("%6.2f, ", corr[i][j]);
      }
      printf("%6.2f\n", corr[i][len-1]);
    }
  
  }else{
    printf(".... fail\n");
  }
  
  return corr;
}


//###################################################################################
//################                                        ###########################
//################        main                            ###########################
//################                                        ###########################
//###################################################################################
int main(int argn,char **argv) {
  if ( argn == 1 ) {
    printf("Usage: \n");
    printf("%s file_in.evt  raw_Opt  timeWidow   correctionFile\n", argv[0]);
    printf("                            |        |        |\n");
    printf("                            |        |        + correction file, row for det, col for order of correction\n");
    printf("                            |        |\n");
    printf("                            |        + when build event, event build window, 1 = 10 ns, default 100\n");
    printf("                            + default 0 = raw, 1 = event build \n");
    ///   std::cerr<<"Usage:\n "<<argv[0]<<" filein.evt fileout correctfilename STARTDETIDnumber maxtime start_buf# stop_buf#\n  
    ///    Converts physics buffers (type 30) from  \'file.evt\' to \'fileout.ev2\' from \'start_buf#\' until\n 
    ///   'stop_buf#'  or whole file if no event numbers are given or \'start buf#\' = 0\n  
    ///   STARTDETIDnum is the number of the dE ADC on a scale of 1 to n.\n  
    ///    maxtime defaults to 100 if not specified\n  
    ///   **Use \'NULL\' for the output file to avoid writing an ev2 file. \n
    ///   Note that the start and end buffer numbers are for the input where each adc \n
    ///   read forms another buffer.\n This version sets the time in event of the start detector to 100 and all other times relative to that. 
    ///   If a correct file is specified, it dorrects the detectors in the file for walk.  
    ///   It does not write out events which do not contain the start detector.\n
    ///   Also it only writes out events in which every ADC times in within the limits specified in the correct.tab file."<<std::endl;
    exit(0);
  }

  printf("=======================================================\n");
  printf("===             XIA evt file to CERN ROOT           ===\n");
  printf("=======================================================\n");
  printf("The start detector number is %d\n", STARTDETID);

  unsigned long long blockNum=0;  /// this is total block
  unsigned long long block30Num=0;  /// this is number of block-30
  unsigned long long eventID=0; /// this is eventID
  
  ///for( int i = 0 ; i < argn; i++) printf("............. %s \n", argv[i]);
  
  int rawOpt = 0;
  if ( argn >= 3 ) rawOpt = atoi(argv[2]);
  
  int timeWindow = 100;
  if ( argn >= 4 ) timeWindow = atoi(argv[3]);
  
  TString corrFileName = "";
  bool hasCorr = false;
  if ( argn == 5 ) {
    corrFileName = argv[4];
    hasCorr = true;
  }
  
  FILE *infile=fopen(argv[1],"r");
  if (infile==NULL) {
    printf("cannot open file : %s \n", argv[1]);
    ///std::cerr<<"Problem opening "<<argv[1]<<std::endl;
    exit(0);
  }

  TString inFileName = argv[1];
  TString outFileName = inFileName;
  outFileName.Remove(inFileName.First('.'));
  if( rawOpt == 0 ) outFileName.Append("_raw");
  outFileName.Append(".root");
  
  printf(" In file : %s \n",  inFileName.Data());
  printf("Out file : %s \n", outFileName.Data());

  std::vector<std::vector<double>> corr;
  if( hasCorr ) corr = LoadCorrectionParameters(corrFileName);

  int chan,slot,chann;
  int pu; /// pile up
  int energy;
  double cEnergy;
  unsigned long long evtime;
  unsigned short cfd;  

  int pileupcount = 0;
  int zerocount = 0;
  int PileUp[64];
  
  const unsigned long maskpu = 2147483648;
  const unsigned long multiplier = 4294967296LL;
  
  double energyA[NUMDET];
  double cEnergyA[NUMDET];
  unsigned long long  timeA[NUMDET];
  int puA[NUMDET];
  long long diffTimeA[NUMDET];
  unsigned short cfdA[NUMDET];
  int multi = 0; /// multipicilty in an event
  int detMulti[NUMDET]; /// multiplicity in a detector in an event

  TFile * outFile = new TFile(outFileName, "RECREATE");
  outFile->cd();
  TTree * tree = new TTree("tree", "tree");
    
  tree->Branch("eventID", &eventID, "event_number/l");

  if ( rawOpt == 0 ){ /// when save raw data
    tree->Branch("chan",     &chan, "chan/I");
    tree->Branch("slot",     &slot, "slot/I");
    tree->Branch("chann",   &chann, "channel number/I");
    tree->Branch("pu",         &pu, "pile-up/I");
    tree->Branch("energy", &energy, "energy/I");
    if( hasCorr) tree->Branch("cEnergy", &cEnergy, "corrected_energy/D");
    tree->Branch("time",   &evtime, "timestamp/l");
    tree->Branch("cfd",       &cfd, "cfd/s");
  }else{ /// when build event by time-window
    tree->Branch("energy",    energyA, Form("energy[%d]/D", NUMDET));
    if( hasCorr) tree->Branch("cEnergy",    cEnergyA, Form("corrected_energy[%d]/D", NUMDET));
    tree->Branch("time",        timeA, Form("timestamp[%d]/l", NUMDET));
    tree->Branch("dtime",   diffTimeA, Form("diff_time[%d]/L", NUMDET));
    tree->Branch("pu",            puA, Form("pile_up[%d]/I", NUMDET));
    tree->Branch("cfd",          cfdA, Form("cfd[%d]/I", NUMDET));
    tree->Branch("multi",      &multi, "multiplicity/I");
    tree->Branch("detMulti", detMulti, Form("det_multiplicity[%d]/I", NUMDET));
  }

  ///change this for 64bit compiler long *bufsam=NULL;
  
  //clear energy and time array
  for( int i = 0; i < NUMDET; i++){
    energyA[i] = TMath::QuietNaN();
    cEnergyA[i] = TMath::QuietNaN();
    timeA[i] = 0;
    diffTimeA[i] = -999;
    cfdA[i] = 0;
    puA[i] = -1;
    detMulti[i] = 0;
  }
  multi = 0;
  
  unsigned long long startTime = 0;
  long long diffTime = 0;

  int bread = 1; 
  int bsam = 2;
  long * bufsiz=new long[bsam];
  unsigned int *bufsam = NULL;

  printf("============ Start looping events | build event ? %s", rawOpt == 0 ? "No" : "Yes");
  if( rawOpt == 1 ) {
    printf(" time window : +- %d click\n", timeWindow);
  }else{
    printf("\n");
  }
  while ( !feof(infile) ) {
  
    // get buffer size
    ///change long -> int for 64 bit
    fread(bufsiz,sizeof(int),bread,infile); /// read 1 int (4 byte) from infile and save to bufsize
    int bsize = bufsiz[0] -4 ;
    if (feof(infile)) break;
    blockNum ++;

    ///change for 64bit  bufsam=new long[bsize/4+1];
    bufsam = new unsigned int[bsize/4+1];
    fread((char*)bufsam, 1, bsize, infile); /// read bsize of 1 byte from infile and save to char
    ///printf("============ bsize : %d \n", bsize);
    ///for( int i = 0; i < bsize; i++) printf("%d, ", bufsam[i]);
    ///printf("\n");

    if (bufsam[0] == 30) {
      block30Num ++;
      chan = (bufsam[2]) & (15);
      slot = ((bufsam[2]) & (240))/16;
      chann = (slot - 2)*16 + chan + 1;
      pu = ((bufsam[2]) & (maskpu))/maskpu;
      energy = ((bufsam[5]) & 65535);
      unsigned long long evtimehi = ((bufsam[4]) & 65535);
      unsigned long long evtimelo = bufsam[3];
      evtime = evtimelo + multiplier*evtimehi;
      cfd = bufsam[4]/65536;

      if ( energy == 0 ) zerocount++;
      if ( pu > 0 ) pileupcount++;
      if ((pu > 0 ) && ( chann > 0 ) && ( chann < 65 )) PileUp[chann-1]++;

      if( blockNum % 100000 == 0 ) printf(".");
      ///if( blockNum % 100000 == 0 ) printf("%llu \n", blockNum);
      ///if( block30Num < 50) printf("b30: %10llu, chan: %d, slot: %d, chann: %2d, pu: %2d, energy: %5d, evtime: %13llu, cfd: %d\n", block30Num, chan, slot, chann, pu, energy, evtime, cfd);

      /// energy correction
      if ( hasCorr ){
        cEnergy = 0;
        int order = (int) corr[chann-1].size();
        for( int i = 0; i < order ; i++){
          cEnergy += corr[chann-1][i] * TMath::Power((double)energy, i);
        }
      } 

      if ( rawOpt == 0 ) {
        eventID++;
        outFile->cd();
        tree->Fill();
      }else{ /// build event

        if ( startTime == 0 ) startTime = evtime;
      
        diffTime = evtime - startTime;
        
        if( -timeWindow < diffTime && diffTime < timeWindow ){

          if( !TMath::IsNaN(energyA[chann-1]) ) detMulti[chann-1] ++;
          energyA[chann-1] = energy;
          cEnergyA[chann-1] = cEnergy;
          timeA[chann-1] = evtime;
          diffTimeA[chann-1] = diffTime;
          puA[chann-1] = pu;
          detMulti[chann-1]++;
          multi++;
          
          
        }else{
          /// fill tree 
          eventID++;
          outFile->cd();
          tree->Fill();
          ///  clear energy and time array
          multi = 0;
          for( int i = 0; i < NUMDET; i++){
            energyA[i] = TMath::QuietNaN();
            cEnergyA[i] = TMath::QuietNaN();
            timeA[i] = 0;
            diffTimeA[i] = -999;
            puA[i] = -1;
            detMulti[i] = 0;
            cfdA[i] = 0;
          }
          
          /// fill the 1st data of a new event
          startTime = evtime;
          energyA[chann-1] = energy;
          cEnergyA[chann-1] = cEnergy;
          timeA[chann-1] = evtime;
          diffTimeA[chann-1] = 0;
          puA[chann-1] = pu;
          detMulti[chann-1]++;
          multi++;
        }
      }
      
    } ///end if bufsam[0]=30
    
  } 

  delete [] bufsiz;
  delete [] bufsam;
  fclose(infile);

  printf("\n============ end of event loop, totoal block read: %llu \n", blockNum);

  eventID++;
  outFile->cd();
  tree->Write();
  outFile->Close();
  
  
  //========================= Print summary
  printf("============================================\n");
  ///printf("         number of block: %llu\n", blockNum); 
  printf(" number of type 30 block: %llu\n", block30Num); 
  printf("             event built: %llu\n", eventID); 
  printf("============================================\n");
    
  return 0;
}

