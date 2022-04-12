// Sourcefile contaning the algorithm that calculate the total energy deposit in pvt cube by gamme decay chain

const double MIN_E_LEVEL = 0.00001; // MEV
const double MAX_E_LEVEL = 4.5; //MEV


double walker(
    ROOT::RVec<int>* pdg,
    ROOT::RVec<double>* epvt,
    ROOT::RVec<int>* parentid,
    ROOT::RVec<int>* trackid,
    ROOT::RVec<int>* parentpdg,
    int parentindex
);


double track_entry(
    ROOT::RVec<int>* pdg,
    ROOT::RVec<double>* epvt,  // All pvt event from the dataset
    ROOT::RVec<int>* parentid, // Were i check if parent correspond to the tracking id
    ROOT::RVec<int>* trackid,  // Were i get the id of all child
    ROOT::RVec<int>* parentpdg
) {
    int n = std::min(parentid->size(), trackid->size());
    const std::vector<int> validid {0, 1, 2, 3};
    for (int i = 0; i < n; i++) {
        if (
            (MIN_E_LEVEL < epvt->at(i) && epvt->at(i) < MAX_E_LEVEL) &&  // Had a valid energy
            (pdg->at(i) == 22) &&  // Ensure it's a photon
            (parentid->at(i) == 0) &&  // Should not have parents
            (std::find(validid.begin(), validid.end(), trackid->at(i)) != validid.end())
        ) {
            return walker(pdg, epvt, parentid, trackid, parentpdg, i);
        }
    }
    return 0.;
}


double walker(
    ROOT::RVec<int>* pdg,  // Multiple same trackid so parent pdg need to be tracked
    ROOT::RVec<double>* epvt,  // All pvt event from the dataset
    ROOT::RVec<int>* parentid, // Were we check if parent correspond to the tracking id
    ROOT::RVec<int>* trackid,  // Were we get the id of all child
    ROOT::RVec<int>* parentpdg,
    int parentindex
) {
    std::vector<int> all_childindex (0);
    int n = std::min({pdg->size(), epvt->size(), parentid->size(), trackid->size(), parentpdg->size()});
    for (int i = 0; i < n; i++) {
        if (
            parentid->at(i) == trackid->at(parentindex) &&
            epvt->at(i) > MIN_E_LEVEL &&
            parentpdg->at(i) == pdg->at(parentindex)
        ) {
            all_childindex.push_back(i);
        }
    }
    int m = all_childindex.size();
    if (m == 0) {
        // recursion get-away
        return epvt->at(parentindex);
    } else {
        // If it's not the end of the chain => recursively continu it
        double etot {epvt->at(parentindex)};
        for (int i = 0; i < m; i++) {
            int index = all_childindex[i];
            etot += walker(pdg, epvt, parentid, trackid, parentpdg, index);  // child become parent
        }
        return etot;
    }
}