// "mu3e" or "alignment/mu3e"
#define TTREE_LOC "mu3e"
int nbins = 20000;   //nbins
int lo = -10;        //ns
int hi = 10;         //ns

uint32_t get_neighbour(uint32_t tile_id, uint32_t direction);
int32_t get_index(std::vector<unsigned int> v, uint32_t K);

// --------------------------------------------
//
int build_histograms() {
  std::vector<TString> files;
  std::vector<char const*> env_vars = {"CONDOR_ATH_INPUT", "CONDOR_INPUT"};


  for (auto env_var : env_vars) {
    char const* content = std::getenv(env_var);
    if (content != NULL){
      std::stringstream ss((content));   // Insert the string into a stream
      std::string buf;                   // Have a buffer string
      while (ss >> buf)
        files.push_back(buf);
      break;
    }

  }

  // -----------------------------------------
  TFile *savefile = new TFile("hist.root", "RECREATE");
  std::cout << "Generating histograms" << '\n';
  std::map<uint32_t, TH1D> m;

  // Station 1
  for (int i=200000; i<202912; i++){
    std::string s_z = std::to_string(i) + "_z";
    std::string s_phi = std::to_string(i) + "_phi";
    char const *hist_name_z = s_z.c_str();
    char const *hist_name_phi = s_phi.c_str();
    m[i] = TH1D(hist_name_z, hist_name_z, nbins, lo, hi);
    m[i+200000] = TH1D(hist_name_phi, hist_name_phi, nbins, lo, hi);
  }

  // Station 2
  for (int i=300000; i<302912; i++){
    std::string s_z = std::to_string(i) + "_z";
    std::string s_phi = std::to_string(i) + "_phi";
    char const *hist_name_z = s_z.c_str();
    char const *hist_name_phi = s_phi.c_str();
    m[i] = TH1D(hist_name_z, hist_name_z, nbins, lo, hi);
    m[i+200000] = TH1D(hist_name_phi, hist_name_phi, nbins, lo, hi);
  }
  // -----------------------------------------
  // std::cout << m[202911].GetName() << '\n';
  // std::cout << m[202911+200000].GetName() << '\n';
  // std::cout << m[302911].GetName() << '\n';
  // std::cout << m[302911+200000].GetName() << '\n';


  for (auto f : files){
    std::cout << f << " -> # Frames: ";
    TFile *sourceFile = new TFile(f, "READ");
    TTree *tree = (TTree *)(sourceFile->Get(TTREE_LOC));
    if(!tree) return 1;

    vector<unsigned int> *tilehit_tile = 0;
    vector<double> *tilehit_time = 0;

    tree->SetBranchAddress("tilehit_tile", &tilehit_tile);
    tree->SetBranchAddress("tilehit_time", &tilehit_time);
    tree->SetBranchStatus("*", 1); // Activates reading of all branches

    uint32_t nEntries = tree->GetEntries();
    std::cout << nEntries << std::endl;

    for (uint32_t frame = 0; frame < nEntries; frame++){
      tree->GetEntry(frame);

      //for (unsigned int index = 0; index < tilehit_tile->size(); ++index)
      //{
      //  std::cout << &tilehit_tile[index] << "  "<< get_neighbour(&tilehit_tile[index], 2) << "  "<< get_neighbour(&tilehit_tile[index], 1) << "\n";
      //}
      //for(auto const& value: *tilehit_tile) {
      //  std::cout << value-*tilehit_tile->begin() << '\n';
      //  std::cout << value << "  "<< get_neighbour(value, 2) << "  "<< get_neighbour(value, 1) << "\n";
      //}
      for (int index = 0, size = tilehit_tile->size(); index < size; ++index)
      {
        uint32_t current_tile_id = tilehit_tile->at(index);
        //std::cout << current_tile_id << "  "<< get_neighbour(current_tile_id, 2) << "  "<< get_neighbour(current_tile_id, 1) << "\n";
        auto tile_id_z_neighbour = get_neighbour(current_tile_id, 2);
        auto tile_id_phi_neighbour = get_neighbour(current_tile_id, 1);

        auto index_z = get_index(*tilehit_tile, tile_id_z_neighbour);
        auto index_phi = get_index(*tilehit_tile, tile_id_phi_neighbour);

        if (index_z >= 0){
          
        }
        if (index_phi >= 0){

        }
      }
    }
  }


  // -----------------------------------------
  // Save everything
  std::cout << "Save File" << '\n';
  for (auto const& [key, val] : m){
    val.Write();
  }
  savefile->Close();
  return 0;
}

// --------------------------------------------
//
// direction:
//    1: up
//    2: right
uint32_t get_neighbour(uint32_t tile_id, uint32_t direction){
  uint32_t neighbour_id = tile_id;
  uint32_t tmp = tile_id - 200000;

  switch (direction) {
    case 1:
      if (tmp >= 100000){
        tmp -= 100000;
      }
      if ((tmp + 1) % 56 == 0){
        neighbour_id -= 55;
      }
      else{
        neighbour_id += 1;
      }
      return neighbour_id;
      break;
    // ---------
    case 2:
      neighbour_id += 56;
      return neighbour_id;
      break;
  }

  return 0;
}

// --------------------------------------------
//
int32_t get_index(std::vector<unsigned int> v, uint32_t K){
    auto it = find(v.begin(), v.end(), K);

    // If element was found
    if (it != v.end())
    {

        // calculating the index
        // of K
        uint32_t index = it - v.begin();

        // tile might be hit more than one time
        it++;
        it = find(it, v.end(), K);
        // if tile id has more than one occurrence in v return -2
        if (it == v.end()){
          return index;
        }
        else{
          return -2;
        }
    }
    else {
        // If the element is not
        // present in the vector
        return -1;
    }
}
