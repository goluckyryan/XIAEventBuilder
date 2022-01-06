

void process_run(TString rootFile){
    
    TFile * file = new TFile(rootFile, "read");
    
    TTree * tree = (TTree *) file->Get("tree");
    
    TString ev2FileName = rootFile;
    ev2FileName.Remove(rootFile.First("."));
    ev2FileName.Append(".ev2");

    tree->Process("Analyzer.C+", ev2FileName);
    
}
