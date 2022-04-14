#include "../include/constants.hpp"


// Sourcefile contaning the algorithm that calculate the total energy deposit in pvt cube by gamme decay chain
// TODO: functions documentation

double walker(
    //T9
    ROOT::RVec<int>* trackid,
    ROOT::RVec<int>* parentid,
    // T5
    ROOT::RVec<int>* trackidT5,
    ROOT::RVec<double>* edep_pvt,
    // Index of the actual particle considered inside the T9 tab
    int partIndex,
    int N
);

int findPrimaryCE(
    ROOT::RVec<int>* pdg,
    ROOT::RVec<int>* parentid,
    ROOT::RVec<double>* initialEkin,
    int gammaId,
    int N
    );

double track_entry(
    //T9
    ROOT::RVec<int>* pdg,
    ROOT::RVec<int>* trackid,
    ROOT::RVec<int>* parentid,
    ROOT::RVec<double>* initialEkin,
    // T5
    ROOT::RVec<int>* trackidT5,
    ROOT::RVec<double>* edep_pvt,
    bool onlyFirstCompton=false
) {
    int N = std::min({
        pdg->size(),
        parentid->size(),
        trackid->size(),
        initialEkin->size(),
        });
    const std::vector<int> validId {1, 2, 3};
    for (int i = 0; i < N; i++) {
        if (
            (std::abs(initialEkin->at(i) - GAMMA_ENERGY) < 0.001) &&  // Had a valid energy
            (pdg->at(i) == 22) &&  // Ensure it's a photon
            (parentid->at(i) == 0) &&  // Should not have parents
            (std::find(validId.begin(), validId.end(), trackid->at(i)) != validId.end())
        ) {
            if (onlyFirstCompton) {
                int indexCE = findPrimaryCE(pdg, parentid, initialEkin, trackid->at(i), N);
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
    ROOT::RVec<int>* trackid,
    ROOT::RVec<int>* parentid,
    // T5
    ROOT::RVec<int>* trackidT5,
    ROOT::RVec<double>* edep_pvt,
    // Index of the actual particle considered inside the T9 tab
    int partIndex,
    int N
) {
    std::vector<int> childIndex (0);
    int partId = trackid->at(partIndex);

    // For all particle look if parentId match the curent particle trackId
    // All trackId are unique, so this is the only check needed
    for (int i = 0; i < N; i++) {
        if (parentid->at(i) == partId) {
            childIndex.push_back(i);
        }
    }
    int nbChild = childIndex.size();
    int sizeT5 = std::min(trackidT5->size(), edep_pvt->size());
    if (nbChild == 0) {
        // recursion get-away
        for (int i = 0; i < sizeT5; i++) {
            if (trackidT5->at(i) == partId) {
                // End of chain particle (No child and deposit is energy)
                return edep_pvt->at(i);
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
            if (trackidT5->at(i) == partId) {
                // The actual particle deposit energy in pvt
                etot += edep_pvt->at(i);
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
    ROOT::RVec<int>* pdg,
    ROOT::RVec<int>* parentid,
    ROOT::RVec<double>* initialEkin,
    int gammaId,
    int N
) {
    int indexCE {-1};
    double energyCE {0.};
    for (int i = 0; i < N; i++) {
        // Looking for primary electron
        if (parentid->at(i) == gammaId && pdg->at(i) == 11) {
            // Maximum energy check could be done in the first if but keep it clear
            if (energyCE < initialEkin->at(i)) {
                energyCE = initialEkin->at(i);
                indexCE = i;
            }
        }
    }
    return indexCE;
}

ROOT::RVec<int> cubeWalker(
    //T9
    ROOT::RVec<int>* trackid,
    ROOT::RVec<int>* parentid,
    // T5
    ROOT::RVec<int>* trackidT5,
    ROOT::RVec<int>* volid,
    // Index of the actual particle considered inside the T9 tab
    int partIndex,
    int N
);

ROOT::RVec<int> cubeTrackEntry(
    //T9
    ROOT::RVec<int>* pdg,
    ROOT::RVec<int>* trackid,
    ROOT::RVec<int>* parentid,
    ROOT::RVec<double>* initialEkin,
    // T5
    ROOT::RVec<int>* trackidT5,
    ROOT::RVec<int>* volid,
    bool onlyFirstCompton=false
) {
    int N = std::min({
        pdg->size(),
        parentid->size(),
        trackid->size(),
        initialEkin->size(),
        });
    const std::vector<int> validId {1, 2, 3};
    for (int i = 0; i < N; i++) {
        if (
            (std::abs(initialEkin->at(i) - GAMMA_ENERGY) < 0.001) &&  // Had a valid energy
            (pdg->at(i) == 22) &&  // Ensure it's a photon
            (parentid->at(i) == 0) &&  // Should not have parents
            (std::find(validId.begin(), validId.end(), trackid->at(i)) != validId.end())
        ) {
            if (onlyFirstCompton) {
                int indexCE = findPrimaryCE(pdg, parentid, initialEkin, trackid->at(i), N);
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
    ROOT::RVec<int>* trackid,
    ROOT::RVec<int>* parentid,
    // T5
    ROOT::RVec<int>* trackidT5,
    ROOT::RVec<int>* volid,
    // Index of the actual particle considered inside the T9 tab
    int partIndex,
    int N
) {
    std::vector<int> childIndex (0);
    int partId = trackid->at(partIndex);

    // For all particle look if parentId match the curent particle trackId
    // All trackId are unique, so this is the only check needed
    for (int i = 0; i < N; i++) {
        if (parentid->at(i) == partId) {
            childIndex.push_back(i);
        }
    }
    int nbChild = childIndex.size();
    int sizeT5 = std::min(trackidT5->size(), volid->size());
    if (nbChild == 0) {
        // recursion get-away
        for (int i = 0; i < sizeT5; i++) {
            if (trackidT5->at(i) == partId) {
                // End of chain particle (No child and deposit is energy)
                return ROOT::RVec<int> {volid->at(i)};
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
            if (trackidT5->at(i) == partId) {
                // The actual particle deposit energy in pvt
                cubes.push_back(volid->at(i));
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
