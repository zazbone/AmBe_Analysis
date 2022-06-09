#pragma once

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>


// Event that should always be discarded (ex: Tritium & Alpha with wrong deposited energy)
const int maskInvalid = -1; 
// Initial neutron and Gamma identifier
const int maskNOISE = 0;
const int maskGAMMA = 1;

// Gamma child
const int maskG_CHILD = 3;
// Gamma compton edge electron and child
const int maskCE_ELECT = 2;
const int maskCE_CHILD = 4;
// NB: `maskCE_ELECT` and `maskCE_CHILD` are also `maskG_CHILD` but specific ones

// Identifier for 4.44 MeV gamma event and CE subprocess 
// Create by ^{12}C^* relaxing, after a n + ^{12}C inelastic process
// Relaxing 4.44 gamma
const int maskNOISE_GAMMA = 7;
// Is childs
const int maskNOISE_GCHILD = 8;
// Gamma compton edge electron and child
const int maskNOISE_CE = 5;
const int maskNOISE_CE_CHILD = 6;
// NB: `maskNOISE_CE` and `maskNOISE_CE_CHILD` are also `maskNOISE_GCHILD` but specific ones




// Sourcefile contaning the algorithm that apply the filter mask to T9 and T5


ROOT::RVec<int> createT9Mask(
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> const& pdg,
    ROOT::RVec<int> const& processId,
    ROOT::RVec<double> const& Ekin
);


void applyToChain(
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> & mask,
    int curentParentId,
    int maskToApply
);


// The recursive function that will recursively traverse the chain of particles created
// and progressivelly fill the mask vector
// Should not be use direcly, this function is used by `decayChainMask` that you should call insted
void maskWalker(
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> & mask,
    int curentParentId,
    const int N,
    bool isCEChain
) ;


// Function to apply to T9 that return a mask array alligned with the other T9 particles.
// Require:
//      [parentid (RVec<int>)]:: T9 parentid vector, use to build the decay tree
//      [trackid (RVec<int>)]:: T9 trackid vector,  use to build the decay tree
//      [pdg (RVec<int>)]:: T9 pdg vector, use to separate electron - gamma - neutron
//      [processId (RVec<int>)]:: T9 processId vector, use to find compton electron process
//      [Ekin (RVec<int>)]:: T9 initialEkin vector, kinetic energy at creation (not deposited energy)
// Return::RVec<int>
//  The mask for T9 event, each interger value correspond to a decay chain.
//  See: `track_chain.hpp` (Neutron chain, 4.44 MeV gamma, Compton edge)
ROOT::RVec<int> decayChainMask(
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> const& pdg,
    ROOT::RVec<int> const& processId,
    ROOT::RVec<double> const& Ekin
);


////                                                                                 ////
// Following function are now unused, they correspond to the algorithm's first version //
////                                                                                 ////


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