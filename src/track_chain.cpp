#include "../include/constants.hpp"
#include "../include/track_chain.hpp"
#include "../include/util.hpp"


// Sourcefile contaning the algorithm that apply the filter mask to T9 and T5


ROOT::RVec<int> createT9Mask(
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> const& pdg,
    ROOT::RVec<int> const& processId,
    ROOT::RVec<double> const& Ekin
) {
    const int N = parentid.size();
    ROOT::RVec<int> mask (N, maskNOISE);
    int initGammaIndex = -1;
    int initNeutronIndex = -1;
    for (int i = 0; i < N; i++) {
        if (parentid[i] == 0) {
            if (pdg[i] == GAMMA_PDG) {
                initGammaIndex = i;
                mask[i] = maskGAMMA;
                applyToChain(parentid, trackid, mask, trackid[i], maskG_CHILD);
            } else {
                initNeutronIndex = i;
            }
            if (
                initGammaIndex != -1 &&
                initNeutronIndex != -1
            ) { break; }
        }
    }

    for (int i = 0; i < N ; i++) {
        if (pdg[i] == GAMMA_PDG && parentid[i] == trackid[initNeutronIndex] && (std::abs(Ekin[i] - 4.438) < 1e-2)) {
            mask[i] = maskNOISE_GAMMA;
            applyToChain(parentid, trackid, mask, trackid[i], maskNOISE_GCHILD);
            break;
        }
    }

    int gammaCEIndex = -1;
    int neutronCEIndex = -1;
    double gammaCEMEkin = 0;
    double neutronCEMEkin = 0;
    for (int i = 0; i < N; i++) {
        if (pdg[i] == ELECTRON_PDG && processId[i] == CE_PROCESS_ID) {
            if (mask[i] == maskG_CHILD && (gammaCEIndex == -1 || gammaCEMEkin < Ekin[i])) {
                gammaCEIndex = i;
                gammaCEMEkin = Ekin[i];
            } else if (mask[i] == maskNOISE_GCHILD && (neutronCEIndex == -1 || neutronCEMEkin < Ekin[i])) {
                neutronCEIndex = i;
                neutronCEMEkin = Ekin[i];
            }
        }
    }
    
    if (gammaCEIndex != -1) {
        mask[gammaCEIndex] = maskCE_ELECT;
        applyToChain(parentid, trackid, mask, trackid[gammaCEIndex], maskCE_CHILD);
    }
    if (neutronCEIndex != -1) {
        mask[neutronCEIndex] = maskNOISE_CE;
        applyToChain(parentid, trackid, mask, trackid[neutronCEIndex], maskNOISE_CE_CHILD);
    }

    // Add alpha and tritium particles to discarded events.
    mask[pdg == 1000020040 || pdg == 1000010030] == -1;
    return mask;
}


// The recursive function that will recursively traverse the chain of particles created
// and progressivelly fill the mask vector
// Should not be use direcly, this function is used by `createT9Mask` that you should call insted
void applyToChain(
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> & mask,
    int curentParentId,
    int maskToApply
) {
    int N = parentid.size();
    for (int i = 0; i < N; i++) {
        if (parentid[i] == curentParentId) {
            mask[i] = maskToApply;
            applyToChain(parentid, trackid, mask, trackid[i], maskToApply);
        }
    }
}



////                                                                                 ////
// Following function are now unused, they correspond to the algorithm's first version //
////                                                                                 ////




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
) {
    for (int i = 0; i < N; i++) {
        if (parentid[i] == curentParentId) {
            if (isCEChain) {
                mask[i]= maskCE_CHILD;
            } else {
                mask[i] = maskG_CHILD;
            }
            maskWalker(parentid, trackid, mask, trackid[i], N, isCEChain);
        }
    }
}



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
) {
    // Vector may have inconsistante size, only well formed rows are considered
    const int N = std::max(
        {parentid.size(), trackid.size(), pdg.size(), Ekin.size()}
    );
    ROOT::RVec<int> mask (N, maskNOISE);

    int gammaId = -1;
    for (int i=0; i < N; i++){
        if (pdg[i]==GAMMA_PDG && std::abs(Ekin[i] - GAMMA_ENERGY) < 10e-6) {
            mask[i] = maskGAMMA;
            gammaId = trackid[i];
            break;
        }
    }
    if (gammaId < 0) {
        // Case where gamma does not reach the detector
        return mask;
    }

    int CEIndex = -1;
    double CEEkin = 0.;
    ROOT::RVec<int> childIndex {};
    for (int i=0; i < N; i++) {
        if (pdg[i] == ELECTRON_PDG && parentid[i] == gammaId) {
            childIndex.push_back(i);
            mask[i]=maskG_CHILD;
            if (CEEkin < Ekin[i] && processId[i] == CE_PROCESS_ID) {  // Higher energy electron found
                CEEkin = Ekin[i];
                if (-1 != CEIndex) {mask[CEIndex] = maskG_CHILD;}
                CEIndex = i;
                mask[CEIndex] = maskCE_ELECT;
            }
        }
    }
    if (childIndex.size() == 0) {
        // case where Gamma 4.44 does not create elecron
        return mask;
    }

    // Recursion into each gamma child
    for (int eIndex: childIndex) {
        if (eIndex == CEIndex) {
            maskWalker(parentid, trackid, mask, trackid[eIndex], N, true);
        } else {
            maskWalker(parentid, trackid, mask, trackid[eIndex], N, false);
        }
    }
    return mask;
}

double track_entry(
    //T9
    ROOT::RVec<int> const& pdg,
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<double> const& initialEkin,
    // T5
    ROOT::RVec<int> const& trackidT5,
    ROOT::RVec<double> const& edep_pvt,
    bool onlyFirstCompton
) {
    int N = std::min({
        pdg.size(),
        parentid.size(),
        trackid.size(),
        initialEkin.size(),
        });
    const std::vector<int> validId {1, 2, 3};
    for (int i = 0; i < N; i++) {
        if (
            (std::abs(initialEkin[i] - GAMMA_ENERGY) < 0.001) &&  // Had a valid energy
            (pdg[i] == 22) &&  // Ensure it's a photon
            (parentid[i] == 0) &&  // Should not have parents
            (std::find(validId.begin(), validId.end(), trackid[i]) != validId.end())
        ) {
            if (onlyFirstCompton) {
                int indexCE = findPrimaryCE(pdg, parentid, initialEkin, trackid[i], N);
                if (indexCE >= 0) {
                    // Bactracking starting at the compton edge electron
                    return walker(trackid, parentid, trackidT5, edep_pvt, indexCE, N);
                }
                return 0.; // No child elecron, no CE energy
            } else {
                // Full gamma energy chain
                return walker(trackid, parentid, trackidT5, edep_pvt, i, N);
            }
        }
    }
    return 0.;
}


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
) {
    std::vector<int> childIndex (0);
    int partId = trackid[partIndex];

    // For all particle look if parentId match the curent particle trackId
    // All trackId are unique, so this is the only check needed
    for (int i = 0; i < N; i++) {
        if (parentid[i] == partId) {
            childIndex.push_back(i);
        }
    }
    int nbChild = childIndex.size();
    int sizeT5 = std::min(trackidT5.size(), edep_pvt.size());
    if (nbChild == 0) {
        // recursion get-away
        for (int i = 0; i < sizeT5; i++) {
            if (trackidT5[i] == partId) {
                // End of chain particle (No child and deposit is energy)
                return edep_pvt[i];
            }
        }
        // Losed particle (Mo child + no deposited energy in pvt)
        return 0.0;  
    } else {
        // If it's not the end of the chain:
        //  -> Take the actual particle energy
        //  -> Recursively continu the tracking for all child and sum the energies
        double etot {0.};
        for (int i = 0; i < sizeT5; i++) {
            if (trackidT5[i] == partId) {
                // The actual particle deposit energy in pvt
                etot += edep_pvt[i];
                break;  // No way to find multiple time the same trackid
            }
        }
        for (int i = 0; i < nbChild; i++) {
            // Child become parent
            etot += walker(trackid, parentid, trackidT5, edep_pvt, childIndex[i], N);
        }
        // Case where the particle had childs, the energy returned is the particle energy in pvt (if there any)
        // and energy of all child in tree
        return etot;
    }
}


// For given gamma trackid (gammaId), search for most energetic electron child
// return the CE electron index in T9 arrays
// If returned index is -1 => no electron child case
int findPrimaryCE(
    ROOT::RVec<int> const& pdg,
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<double> const& initialEkin,
    int gammaId,
    int N
) {
    int indexCE {-1};
    double energyCE {0.};
    for (int i = 0; i < N; i++) {
        // Looking for primary electron
        if (parentid[i] == gammaId && pdg[i] == 11) {
            // Maximum energy check could be done in the first if but keep it clear
            if (energyCE < initialEkin[i]) {
                energyCE = initialEkin[i];
                indexCE = i;
            }
        }
    }
    return indexCE;
}

ROOT::RVec<int> cubeTrackEntry(
    //T9
    ROOT::RVec<int> const& pdg,
    ROOT::RVec<int> const& trackid,
    ROOT::RVec<int> const& parentid,
    ROOT::RVec<double> const& initialEkin,
    // T5
    ROOT::RVec<int> const& trackidT5,
    ROOT::RVec<int> const& volid,
    bool onlyFirstCompton
) {
    int N = std::min({
        pdg.size(),
        parentid.size(),
        trackid.size(),
        initialEkin.size(),
        });
    const std::vector<int> validId {1, 2, 3};
    for (int i = 0; i < N; i++) {
        if (
            (std::abs(initialEkin[i] - GAMMA_ENERGY) < 0.001) &&  // Had a valid energy
            (pdg[i] == 22) &&  // Ensure it's a photon
            (parentid[i] == 0) &&  // Should not have parents
            (std::find(validId.begin(), validId.end(), trackid[i]) != validId.end())
        ) {
            if (onlyFirstCompton) {
                int indexCE = findPrimaryCE(pdg, parentid, initialEkin, trackid[i], N);
                if (indexCE >= 0) {
                    // Bactracking starting at the compton edge electron
                    return cubeWalker(trackid, parentid, trackidT5, volid, indexCE, N);
                }
                return ROOT::RVec<int> {};; // No child elecron, no CE energy
            } else {
                // Full gamma energy chain
                return cubeWalker(trackid, parentid, trackidT5, volid, i, N);
            }
        }
    }
    // No gamma case, T1 gamma does not reach the detector
    return ROOT::RVec<int> {};
}


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
) {
    std::vector<int> childIndex (0);
    int partId = trackid[partIndex];

    // For all particle look if parentId match the curent particle trackId
    // All trackId are unique, so this is the only check needed
    for (int i = 0; i < N; i++) {
        if (parentid[i] == partId) {
            childIndex.push_back(i);
        }
    }
    int nbChild = childIndex.size();
    int sizeT5 = std::min(trackidT5.size(), volid.size());
    if (nbChild == 0) {
        // recursion get-away
        for (int i = 0; i < sizeT5; i++) {
            if (trackidT5[i] == partId) {
                // End of chain particle (No child and deposit is energy)
                return ROOT::RVec<int> {volid[i]};
            }
        }
        // Losed particle (Mo child + no deposited energy in pvt)
        return ROOT::RVec<int> {}; 
    } else {
        // If it's not the end of the chain:
        //  -> Take the actual particle energy
        //  -> Recursively continu the tracking for all child and sum the energies
        ROOT::RVec<int> cubes {};
        for (int i = 0; i < sizeT5; i++) {
            if (trackidT5[i] == partId) {
                // The actual particle deposit energy in pvt
                cubes.push_back(volid[i]);
                break;  // No way to find multiple time the same trackid
            }
        }
        for (int i = 0; i < nbChild; i++) {
            // Child become parent
            auto childCubes = cubeWalker(trackid, parentid, trackidT5, volid, childIndex[i], N);
            cubes.insert(cubes.end(), childCubes.begin(), childCubes.end());
        }
        // Case where the particle had childs, the energy returned is the particle energy in pvt (if there any)
        // and energy of all child in tree
        return cubes;
    }
}
