#include <vector>
#include <cstdlib>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TRootCanvas.h>
#include <TLatex.h>
#include <TH2D.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TClass.h>
#include <TTreeIndex.h>
#include "./include/structs.h"
#include "./include/time_cluster.h"	

using namespace std;


int main(int argc, char ** argv){
    //define time threshold
    Double_t time_threshold = 10; //[ns]

    if(argc < 3){
        cout << "Too few arguments: Usage: infile.root outfile.root" << endl;
        return -1;
    }

    char * infile_name = argv[1];
    char * outfile_name = argv[2];
    
    cout << "input file name: " << infile_name << " output file name: " << outfile_name << endl;
    
    //read mu3e tree
    TFile * infile = new TFile(infile_name,"READ");
    TTree * intree = (TTree*)(infile->Get("mu3e"));
    //intree->Print();

    //Deactivate all branches
    intree->SetBranchStatus("*", 0);
 
    //Activate only necessary branches
    for (auto activeBranchName : {"tilehit_time", "tilehit_tile", "tilehit_edep"}){
        intree->SetBranchStatus(activeBranchName, 1);
    }
    intree->Print();

    //number of entries
    uint32_t nentries = intree->GetEntries();
    cout << "Number of entries: " << nentries << endl;

    // Open root file for writing
	TFile * outfile = new TFile(outfile_name, "RECREATE");
    	if(!outfile || outfile->IsZombie()){
        cout << "Could not open root output file" << endl;
        return -1;
    }
    //make a new tree with the same structure as the input tree 
    TTree *outtree = intree->CloneTree(0); 
    //outtree->Write();
    //outtree->Print();


    //sort tree by time
    /*
    outtree->BuildIndex("0","tilehit_time");
    TTreeIndex* I = (TTreeIndex*)outtree->GetTreeIndex(); // get the tree index
    Long64_t* sort_index = I->GetIndex(); //create an array of entries in sorted order
    outtree->SetNotify(I); //set the tree to use the index?
    */

    //create time cluster branch
    vector<Int_t> time_cluster;
    outtree->Branch("time_cluster", &time_cluster);

    //gROOT->ProcessLine("#include <vector>");

    //initialize counters
    Int_t nclusters = 0;
    Double_t tmp_time_reference = 0;

    //initialize vectors
    vector<Double_t> *tilehit_times = 0;
    vector<Int_t> *tilehit_ids = 0;
    vector<Double_t> *tilehit_edeps = 0;

    //get all tile hit times, ids, and edeps
    intree->SetBranchAddress("tilehit_time", &tilehit_times);
    intree->SetBranchAddress("tilehit_tile", &tilehit_ids);
    intree->SetBranchAddress("tilehit_edep", &tilehit_edeps);


    //clustering loop
    for (uint32_t i=0; i < nentries; i++){
        intree->GetEntry(i);
        if (!tilehit_times->size() == 0){
            tmp_time_reference = tilehit_times->at(0);
            //for (uint32_t j=0; j < tilehit_times->size(); j++){
            for (auto time : *tilehit_times){
                if (abs(time - tmp_time_reference) < time_threshold){
                    time_cluster.push_back(nclusters);
                }
                else{
                    tmp_time_reference = time;
                    nclusters++;
                    time_cluster.push_back(nclusters);
                }

            }
        }
        outtree->Fill();
        time_cluster.clear();
    }

    
outtree->Write();
outfile->Close();

cout << "ROOT file has been written" << endl;

return 0;
}