#include "TChain.h"

void ChainAnalysis(){

  TChain * chain = new TChain("tree");

  /******************/

  //chain->Add("Dec15-07-00[0-9].root");

  chain->Add("ti74pt7a-*.root");

  /******************/
  printf("\033[31m");
  chain->GetListOfFiles()->Print();
  printf("\033[m");
  chain->Process("Analyzer.C+");


}
