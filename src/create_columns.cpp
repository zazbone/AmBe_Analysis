#include "../include/create_columns.hpp"


// TODO: Define column name as constants and document it
void EQuenched(const char* G4Path, const char* outputFile) {
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
        "chainMaskT9",
        [](VecI pid, VecI tip, VecI pdg, VecUI CPID, VecD Ekin){
            return decayChainMask(pid, tip, pdg, CPID, Ekin);
        },
        {"T9.parentid", "T9.trackid", "T9.pdg", "T9.CreatorProcessID", "T9.initialEkin"}
    );
    df = createT5Mask(df);
    df = df.Redefine(
        "chainMaskT5",
        [](VecI cmask, VecI pdg){
            return ROOT::VecOps::Where(cmask==maskNOISE && (pdg == 1000020040 || pdg == 1000010030), ROOT::RVec(cmask.size(), -1), cmask);
        },
        {"chainMaskT5", "pdg"}
    );


    df = df.Define("volidAll",
        [](VecI volid, VecD Epvt, VecI chainMask){
            // chain mask that sghould be discard any time are labeled with -1
            ROOT::RVec<bool> mask = chainMask != -1;
            return volidEpvtTot(volid[mask], Epvt[mask]);
        },
        {"volid", "E_quenched", "chainMaskT5"}
    );
    df = df.Define("EpvtAll", [](VectID volid){return std::get<1>(volid);}, {"volidAll"});
    df = df.Redefine("volidAll", [](VectID volid){return std::get<0>(volid);}, {"volidAll"});
    df = df.Define("_indexRankedCube", [](VecD Epvt){return ROOT::VecOps::Argsort(Epvt, [](double x, double y) {return x > y;});}, {"EpvtAll"});
    df = df.Define("EpvtRankedCube", [](VecUL index, VecD Epvt){return ROOT::VecOps::Take(Epvt, index);}, {"_indexRankedCube", "EpvtAll"});
    df = df.Define("volidRankedCube", [](VecUL index, VecI volid){return ROOT::VecOps::Take(volid, index);}, {"_indexRankedCube", "volidAll"});

    df = df.Define("volidCE",
        [](VecI volid, VecD Epvt, VecI chainMask){
            ROOT::RVec<bool> mask = chainMask == maskCE_ELECT || chainMask == maskCE_CHILD;
            return volidEpvtTot(volid[mask], Epvt[mask]);
        },
        {"volid", "E_quenched", "chainMaskT5"}
    );
    df = df.Define("EpvtCE", [](VectID volid){return std::get<1>(volid);}, {"volidCE"});
    df = df.Redefine("volidCE", [](VectID volid){return std::get<0>(volid);}, {"volidCE"});
    df = df.Define("totalEnergyCE", [](VecD Epvt){
        return ROOT::VecOps::Sum(Epvt[Epvt > 0.1]);
    }, {"EpvtCE"});
    df = df.Define("_indexMECE", [](VecD Epvt){return indexMax(Epvt);}, {"EpvtCE"});
    df = df.Define("EpvtMECE", [](int i, VecD Epvt){return i == -1? 0.0: Epvt[i];}, {"_indexMECE", "EpvtCE"});
    df = df.Define("volidMECE", [](int i, VecI volid){return i == -1? -1: volid.at(i);}, {"_indexMECE", "volidCE"});

    df = df.Define("volidGamma",
        [](VecI volid, VecD Epvt, VecI chainMask){
            ROOT::RVec<bool> mask = chainMask == maskGAMMA || chainMask == maskG_CHILD;
            return volidEpvtTot(volid[mask], Epvt[mask]);
        },
        {"volid", "E_quenched", "chainMaskT5"}
    );
    df = df.Define("EpvtGamma", [](VectID volid){return std::get<1>(volid);}, {"volidGamma"});
    df = df.Redefine("volidGamma", [](VectID volid){return std::get<0>(volid);}, {"volidGamma"});

    df = df.Define("volidNeutron",
        [](VecI volid, VecD Epvt, VecI chainMask){
            ROOT::RVec<bool> mask = chainMask == maskNOISE;
            return volidEpvtTot(volid[mask], Epvt[mask]);
        },
        {"volid", "E_quenched", "chainMaskT5"}
    );
    df = df.Define("EpvtNeutron", [](VectID volid){return std::get<1>(volid);}, {"volidNeutron"});
    df = df.Redefine("volidNeutron", [](VectID volid){return std::get<0>(volid);}, {"volidNeutron"});

    df = df.Define("CEReachPvt",
        [](bool T1, VecI mask){return T1 && ROOT::VecOps::Any(mask == maskGAMMA) && ROOT::VecOps::Any(mask == maskCE_ELECT);},
        {"is_event", "chainMaskT5"}
    );
    
    df.Snapshot("A", outputFile, {
        "is_event", "chainMaskT9", "chainMaskT5", "CEReachPvt",
        "EpvtAll", "volidAll", "EpvtCE", "volidCE", "EpvtGamma", "volidGamma", "EpvtNeutron", "volidNeutron",
        "totalEnergyCE", "EpvtMECE", "volidMECE",
        "volidRankedCube", "EpvtRankedCube"
    });
    // G4File->Close(); // why closing the file cause a segfault ???
}