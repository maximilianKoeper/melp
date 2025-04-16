#ifndef STRUCTS_HH
#define STRUCTS_HH

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include "TGraph.h"
#include "TCanvas.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::ofstream;


//cluster hit struct
struct ClusterHit{
    Int_t tile_id = 0;
    Int_t frame_id = -1;
    Double_t time = -1;
    Double_t edep = -1;
};

//cluster struct
struct Cluster{
    Int_t id;
    Int_t frame_id;
    Int_t master_id;
    vector<ClusterHit> hits;

    //returns the number of hits in the cluster
    Int_t GetNHits(){
        return hits.size();
    }
    
    //returns all tile ids in cluster
    vector<Int_t> GetTileIDs(){
        vector<Int_t> tile_ids;
        for (uint32_t i=1; i <= hits.size(); i ++){
            tile_ids.push_back(hits[i].tile_id);
        }
        return tile_ids;
    }
    
    //returns all times in cluster
    vector<Double_t> GetTimes(){
        vector<Double_t> times;
        for (uint32_t i=1; i <= hits.size(); i ++){
            times.push_back(hits[i].time);
        }
        return times;
    }
    
    //returns all edeps in cluster
    vector<Double_t> GetEdeps(){
        vector<Double_t> edeps;
        for (uint32_t i=1; i <= hits.size(); i ++){
            edeps.push_back(hits[i].edep);
        }
        return edeps;
    }

    //return all frame ids in cluster
    vector<Int_t> GetFrameIDs(){
        vector<Int_t> frame_ids;
        for (uint32_t i=1; i <= hits.size(); i ++){
            frame_ids.push_back(hits[i].frame_id);
        }
        return frame_ids;
    }

};



#endif /*STRUCTS_HH*/