typedef struct {
    ROOT::RDataFrame df;
    std::shared_ptr<TTree> tree;
    std::shared_ptr<TFile> file;
} DataSetHolder;

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
    ROOT::RDataFrame df {
        *tree,
        {"T1.pdg", "T1.energy",
        "pdg", "trackid", "parentid", "volid",
        "T9.pdg", "T9.trackid", "T9.parentid", "T9.parentPDG"} // TODO: find why collumn selection does not work
        };
    DataSetHolder out {df, tree, ambe3file};
    return out;
}