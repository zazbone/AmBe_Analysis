#include "../include/load_dataset.hpp"


// Generic function to load GEANT4 sim dataset
// Return tree and file for lifetime purpose
// Or may be to use tree or extend it
DataSetHolder ambeDataSet (char const* path) {
    // Lifetime in c++ is hell
    std::shared_ptr<TFile> ambe3file{TFile::Open(path)};
    std::shared_ptr<TTree> tree{ambe3file->Get<TTree>("T5")};
    // Add other friend manualy to tree field of DataSetHolder and create a new df
    tree->AddFriend("T1", path);
    tree->AddFriend("T9", path);
    ROOT::RDF::RNode df = ROOT::RDataFrame {*tree};
    DataSetHolder out {df, {tree}, {ambe3file}};
    return out;
}

// Load T5 tree from Geant4 dataset with S2 reconstruction data
// From Safran2 we extract CCubeE from ES and nsESClusCoins
// S2 does not report all G4 row, missing rows will be filled by empty RVec
DataSetHolder withS2(char const* G4path, char const* S2path) {
    std::shared_ptr<TFile> G4file{TFile::Open(G4path)};
    std::shared_ptr<TTree> G4tree{G4file->Get<TTree>("T5")};
    std::shared_ptr<TFile> S2file{TFile::Open(S2path)};
    std::shared_ptr<TTree> S2tree{S2file->Get<TTree>("DefaultReconstruction/ES")};
    S2tree->BuildIndex("cycleNum");
    G4tree->SetAlias("cycleNum", "LocalEntry$");
    G4tree->AddFriend(S2tree.get());
    ROOT::RDF::RNode df = ROOT::RDataFrame {*G4tree};
    return DataSetHolder {df, {G4tree, S2tree}, {G4file, S2file}};
}