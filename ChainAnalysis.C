#include "TChain.h"

void ChainAnalysis(){

  TChain * chain = new TChain("tree");

  chain->Add("Dec15-07-00[0-9].root");

  chain->GetListOfFiles()->Print();
  chain->Process("Analyzer.C+");


}
