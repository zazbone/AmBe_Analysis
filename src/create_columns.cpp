#include "../include/create_columns.hpp"


// Create a new root file with all neccessary practical column for the analysis
// [G4path (char*)]:: the path to a source G4 ambe callibration 
// [outputfile (char*)]:: the name and path of the output file
void EQuenched(char const* G4Path, char const* outputFile) {
    using VecD = ROOT::RVec<double> const&;
    using VecI = ROOT::RVec<int> const&;
    using VecUI = ROOT::RVec<unsigned int> const&;
    using VecUL = ROOT::RVec<unsigned long> const&;
    using VectID = std::tuple<ROOT::RVec<int>, ROOT::RVec<double>> const&;

    std::unique_ptr<TFile> G4File {TFile::Open(G4Path)};

    std::unique_ptr<TTree> mainTree{G4File->Get<TTree>("T5")};
    mainTree->AddFriend("T1", G4Path);
    mainTree->AddFriend("T9", G4Path);
    ROOT::RDF::RNode df = ROOT::RDataFrame{*mainTree};
    
    df = filterT1Events(df, 1e-5);
    df = df.Define(
        MASK_T9,
        [](VecI pid, VecI tip, VecI pdg, VecUI CPID, VecD Ekin){
            return createT9Mask(pid, tip, pdg, CPID, Ekin);
        },
        {T9_PARENT_ID, T9_TRACK_ID, T9_PDG, T9_CREATOR_PROCESS_ID, T9_INITIAL_EKIN}
    );
    df = createT5Mask(df);
    df = df.Redefine(
        MASK_T5,
        [](VecI cmask, VecI pdg){
            return ROOT::VecOps::Where(cmask==maskNOISE && (pdg == 1000020040 || pdg == 1000010030), ROOT::RVec(cmask.size(), -1), cmask);
        },
        {MASK_T5, T5_PDG}
    );


    df = df.Define(VOLID_ALL,
        [](VecI volid, VecD Epvt, VecI chainMask){
            // chain mask that sghould be discard any time are labeled with -1
            ROOT::RVec<bool> mask = chainMask != -1;
            return volidEpvtTot(volid[mask], Epvt[mask]);
        },
        {T5_VOLID, T5_E_QUENCHED, MASK_T5}
    );
    df = df.Define(E_QUENCHED_ALL, [](VectID volid){return std::get<1>(volid);}, {VOLID_ALL});
    df = df.Redefine(VOLID_ALL, [](VectID volid){return std::get<0>(volid);}, {VOLID_ALL});
    df = df.Define("_indexRankedCube", [](VecD Epvt){return ROOT::VecOps::Argsort(Epvt, [](double x, double y) {return x > y;});}, {E_QUENCHED_ALL});
    df = df.Redefine(E_QUENCHED_ALL, [](VecUL index, VecD Epvt){return ROOT::VecOps::Take(Epvt, index);}, {"_indexRankedCube", E_QUENCHED_ALL});
    df = df.Redefine(VOLID_ALL, [](VecUL index, VecI volid){return ROOT::VecOps::Take(volid, index);}, {"_indexRankedCube", VOLID_ALL});

    df = df.Define(VOLID_NOISE,
        [](VecI volid, VecD Epvt, VecI chainMask){
            ROOT::RVec<bool> mask = chainMask == maskNOISE;
            return volidEpvtTot(volid[mask], Epvt[mask]);
        },
        {T5_VOLID, T5_E_QUENCHED, MASK_T5}
    );
    df = df.Define(E_QUENCHED_NOISE, [](VectID volid){return std::get<1>(volid);}, {VOLID_NOISE});
    df = df.Redefine(VOLID_NOISE, [](VectID volid){return std::get<0>(volid);}, {VOLID_NOISE});
    df = df.Redefine("_indexRankedCube", [](VecD Epvt){return ROOT::VecOps::Argsort(Epvt, [](double x, double y) {return x > y;});}, {E_QUENCHED_NOISE});
    df = df.Redefine(E_QUENCHED_NOISE, [](VecUL index, VecD Epvt){return ROOT::VecOps::Take(Epvt, index);}, {"_indexRankedCube", E_QUENCHED_NOISE});
    df = df.Redefine(VOLID_NOISE, [](VecUL index, VecI volid){return ROOT::VecOps::Take(volid, index);}, {"_indexRankedCube", VOLID_NOISE});

    df = df.Define(VOLID_GAMMA,
        [](VecI volid, VecD Epvt, VecI chainMask){
            ROOT::RVec<bool> mask = chainMask == maskGAMMA || chainMask == maskG_CHILD;
            return volidEpvtTot(volid[mask], Epvt[mask]);
        },
        {T5_VOLID, T5_E_QUENCHED, MASK_T5}
    );
    df = df.Define(E_QUENCHED_GAMMA, [](VectID volid){return std::get<1>(volid);}, {VOLID_GAMMA});
    df = df.Redefine(VOLID_GAMMA, [](VectID volid){return std::get<0>(volid);}, {VOLID_GAMMA});
    df = df.Redefine("_indexRankedCube", [](VecD Epvt){return ROOT::VecOps::Argsort(Epvt, [](double x, double y) {return x > y;});}, {E_QUENCHED_GAMMA});
    df = df.Redefine(E_QUENCHED_GAMMA, [](VecUL index, VecD Epvt){return ROOT::VecOps::Take(Epvt, index);}, {"_indexRankedCube", E_QUENCHED_GAMMA});
    df = df.Redefine(VOLID_GAMMA, [](VecUL index, VecI volid){return ROOT::VecOps::Take(volid, index);}, {"_indexRankedCube", VOLID_GAMMA});
    
    df = df.Define(VOLID_CE,
        [](VecI volid, VecD Epvt, VecI chainMask){
            ROOT::RVec<bool> mask = chainMask == maskCE_ELECT || chainMask == maskCE_CHILD;
            return volidEpvtTot(volid[mask], Epvt[mask]);
        },
        {T5_VOLID, T5_E_QUENCHED, MASK_T5}
    );
    df = df.Define(E_QUENCHED_CE, [](VectID volid){return std::get<1>(volid);}, {VOLID_CE});
    df = df.Redefine(VOLID_CE, [](VectID volid){return std::get<0>(volid);}, {VOLID_CE});
    df = df.Redefine("_indexRankedCube", [](VecD Epvt){return ROOT::VecOps::Argsort(Epvt, [](double x, double y) {return x > y;});}, {E_QUENCHED_CE});
    df = df.Redefine(E_QUENCHED_CE, [](VecUL index, VecD Epvt){return ROOT::VecOps::Take(Epvt, index);}, {"_indexRankedCube", E_QUENCHED_CE});
    df = df.Redefine(VOLID_CE, [](VecUL index, VecI volid){return ROOT::VecOps::Take(volid, index);}, {"_indexRankedCube", VOLID_CE});
    
    df = df.Define(VOLID_NOISE_GAMMA,
        [](VecI volid, VecD Epvt, VecI chainMask){
            ROOT::RVec<bool> mask = chainMask == maskNOISE_GAMMA || chainMask == maskNOISE_GCHILD;
            return volidEpvtTot(volid[mask], Epvt[mask]);
        },
        {T5_VOLID, T5_E_QUENCHED, MASK_T5}
    );
    df = df.Define(E_QUENCHED_NOISE_GAMMA, [](VectID volid){return std::get<1>(volid);}, {VOLID_NOISE_GAMMA});
    df = df.Redefine(VOLID_NOISE_GAMMA, [](VectID volid){return std::get<0>(volid);}, {VOLID_NOISE_GAMMA});
    df = df.Redefine("_indexRankedCube", [](VecD Epvt){return ROOT::VecOps::Argsort(Epvt, [](double x, double y) {return x > y;});}, {E_QUENCHED_NOISE_GAMMA});
    df = df.Redefine(E_QUENCHED_NOISE_GAMMA, [](VecUL index, VecD Epvt){return ROOT::VecOps::Take(Epvt, index);}, {"_indexRankedCube", E_QUENCHED_NOISE_GAMMA});
    df = df.Redefine(VOLID_NOISE_GAMMA, [](VecUL index, VecI volid){return ROOT::VecOps::Take(volid, index);}, {"_indexRankedCube", VOLID_NOISE_GAMMA});

    df = df.Define(VOLID_NOISE_CE,
        [](VecI volid, VecD Epvt, VecI chainMask){
            ROOT::RVec<bool> mask = chainMask == maskNOISE_CE || chainMask == maskNOISE_CE_CHILD;
            return volidEpvtTot(volid[mask], Epvt[mask]);
        },
        {T5_VOLID, T5_E_QUENCHED, MASK_T5}
    );
    df = df.Define(E_QUENCHED_NOISE_CE, [](VectID volid){return std::get<1>(volid);}, {VOLID_NOISE_CE});
    df = df.Redefine(VOLID_NOISE_CE, [](VectID volid){return std::get<0>(volid);}, {VOLID_NOISE_CE});
    df = df.Redefine("_indexRankedCube", [](VecD Epvt){return ROOT::VecOps::Argsort(Epvt, [](double x, double y) {return x > y;});}, {E_QUENCHED_NOISE_CE});
    df = df.Redefine(E_QUENCHED_NOISE_CE, [](VecUL index, VecD Epvt){return ROOT::VecOps::Take(Epvt, index);}, {"_indexRankedCube", E_QUENCHED_NOISE_CE});
    df = df.Redefine(VOLID_NOISE_CE, [](VecUL index, VecI volid){return ROOT::VecOps::Take(volid, index);}, {"_indexRankedCube", VOLID_NOISE_CE});

    df = df.Define(TOTAL_ALL, [](VecD Epvt){return ROOT::VecOps::Sum(Epvt[Epvt > ETOL]);}, {E_QUENCHED_ALL});
    df = df.Define(TOTAL_NOISE, [](VecD Epvt){return ROOT::VecOps::Sum(Epvt[Epvt > ETOL]);}, {E_QUENCHED_NOISE});
    df = df.Define(TOTAL_GAMMA, [](VecD Epvt){return ROOT::VecOps::Sum(Epvt[Epvt > ETOL]);}, {E_QUENCHED_GAMMA});
    df = df.Define(TOTAL_CE, [](VecD Epvt){return ROOT::VecOps::Sum(Epvt[Epvt > ETOL]);}, {E_QUENCHED_CE});
    df = df.Define(TOTAL_NOISE_GAMMA, [](VecD Epvt){return ROOT::VecOps::Sum(Epvt[Epvt > ETOL]);}, {E_QUENCHED_NOISE_GAMMA});
    df = df.Define(TOTAL_NOISE_CE, [](VecD Epvt){return ROOT::VecOps::Sum(Epvt[Epvt > ETOL]);}, {E_QUENCHED_NOISE_CE});


    df.Snapshot("A", outputFile, {
        MASK_T9, MASK_T5,
        VOLID_ALL, VOLID_NOISE, VOLID_GAMMA, VOLID_CE, VOLID_NOISE_GAMMA, VOLID_NOISE_CE,
        E_QUENCHED_ALL, E_QUENCHED_NOISE, E_QUENCHED_GAMMA, E_QUENCHED_CE, E_QUENCHED_NOISE_GAMMA, E_QUENCHED_NOISE_CE,
        TOTAL_ALL, TOTAL_NOISE, TOTAL_GAMMA, TOTAL_CE, TOTAL_NOISE_GAMMA, TOTAL_NOISE_CE
    });
    // unique_ptr call destructor here => dont need to Close the TFile
}