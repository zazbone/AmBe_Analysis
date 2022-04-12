{ 
    ROOT::RVec<int> pdg = {2102, 22, 2102, 11, 12, 22};
    ROOT::RVec<double> epvt = {100., 1.2, 100., 1.0, 0.8, 1.4};
    ROOT::RVec<int> trackid = {1, 1, 1, 2, 3, 4};
    ROOT::RVec<int> parentid = {0, 0, 0, 1, 1, 2};
    ROOT::RVec<int> parentpdg = {0, 0, 0, 22, 22, 11};
    double res = track_entry(&pdg, &epvt, &parentid, &trackid, &parentpdg);
    std::cout << res << endl;
}