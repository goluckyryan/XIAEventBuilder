#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>

using namespace std;

int main(int argc, char **argv){
    
  printf("=====================================\n");
  printf("===      ev2 --> txt              ===\n");
  printf("=====================================\n");
  
  if ( argc != 2 && argc != 3 ){
    printf("Usage:  \n");
    printf("%10s  [ev2 file]  [max number of event]\n", argv[0]);
  }
  
  string inFileName = argv[1];
  string outFileName = argv[1];
  int found = outFileName.find('.');
  int len =  outFileName.length();
  outFileName.erase(found, len-found);
  outFileName.append(".txt");
  
  long long maxNumEvent = 0; // if maxNumEvent <= 0, all event 
  if( argc >= 3 ) maxNumEvent = atoi(argv[2]);


  //open list-mode data file from PXI digitizer  
  FILE *inFile = fopen(argv[1], "r");
  long int inFileSize,inFilepos;
  if ( inFile == NULL) {
    printf("Error, cannot open input file %s\n", argv[1]);
    return 1;
  }
  
  //get file size
  fseek(inFile, 0L, SEEK_END);
  inFileSize = ftell(inFile);
  rewind(inFile);
  long int inFilePos = 0;

  printf("  in file : %s | file size : %f MB\n", inFileName.c_str(), inFileSize/1024./1024.);
  //printf(" out file : %s \n", outFileName.c_str());
  if( maxNumEvent <= 0 ){
    printf(" max number of event : ALL \n");
  }else{
    printf(" max number of event : %lld \n", maxNumEvent);
  }
  printf("-------------------------------------\n");
  
  long long nEvent = 0;
  short * multi = new short[1];
  short * id = new short[1];
  short * energy  = new short[1] ;
  
  while ( inFilePos < inFileSize || feof(inFile) ) {
    
    printf("------------------- %lld\n", nEvent);
    
    fread(multi, 1, 1, inFile);
    printf(" number of det : %d \n", multi[0]);

    for( int i = 0; i < multi[0] ; i++){
      
      fread(id, 1, 1, inFile);
      fread(energy, 2, 1, inFile);
      
      printf("%2d| %3d, %8d \n", i, id[0], energy[0]);
      
    }
    fread(multi, 1, 1, inFile);
    
    nEvent++;
    
    if( maxNumEvent > 0 &&  nEvent > maxNumEvent ) break;
    
  }

}
