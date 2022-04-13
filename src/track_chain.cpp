#include "../include/constants.hpp"


// Sourcefile contaning the algorithm that calculate the total energy deposit in pvt cube by gamme decay chain


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


double track_entry(
    //T9
    ROOT::RVec<int>* pdg,
    ROOT::RVec<int>* trackid,
    ROOT::RVec<int>* parentid,
    ROOT::RVec<double>* initialEkin,
    // T5
    ROOT::RVec<int>* trackidT5,
    ROOT::RVec<double>* edep_pvt
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
            return walker(trackid, parentid, trackidT5, edep_pvt, i, N);
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