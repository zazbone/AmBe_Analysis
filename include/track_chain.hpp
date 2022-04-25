#pragma once

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>


const int maskNOISE = 0;
const int maskGAMMA = 1;
const int maskCE_ELECT = 2;
const int maskG_CHILD = 3;
const int maskCE_CHILD = 4;


void maskWalker(
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> & mask,
    int curentParentId,
    const int N,
    bool isCEChain
);

ROOT::RVec<int> decayChainMask(
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> const& pdg,
    ROOT::RVec<int> const& processId,
    ROOT::RVec<double> const& Ekin
);

double walker(
    //T9
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> const& parentid,
    // T5
    ROOT::RVec<int> const& trackidT5,
    ROOT::RVec<double> const& edep_pvt,
    // Index of the actual particle considered inside the T9 tab
    int partIndex,
    int N
);

int findPrimaryCE(
    ROOT::RVec<int> const& pdg,
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<double> const& initialEkin,
    int gammaId,
    int N
);

double track_entry(
    //T9
    ROOT::RVec<int> const& pdg,
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<double> const& initialEkin,
    // T5
    ROOT::RVec<int> const& trackidT5,
    ROOT::RVec<double> const& edep_pvt,
    bool onlyFirstCompton=false
);

ROOT::RVec<int> cubeWalker(
    //T9
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> const& parentid,
    // T5
    ROOT::RVec<int> const& trackidT5,
    ROOT::RVec<int> const& volid,
    // Index of the actual particle considered inside the T9 tab
    int partIndex,
    int N
);

ROOT::RVec<int> cubeTrackEntry(
    //T9
    ROOT::RVec<int> const& pdg,
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<double> const& initialEkin,
    // T5
    ROOT::RVec<int> const& trackidT5,
    ROOT::RVec<int> const& volid,
    bool onlyFirstCompton=false
);