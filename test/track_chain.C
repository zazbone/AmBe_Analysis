{ 
    ROOT::RVec<int> pdg = {22, 11, 11, 11};
    ROOT::RVec<int> trackid = {2, 27, 28, 29};
    ROOT::RVec<int> parentid = {0, 2, 2, 2};
    ROOT::RVec<double> initialEkin = {4.4400, 0, 0, 0};
    ROOT::RVec<int> trackidT5 = {29, 28, 27, 2};
    ROOT::RVec<double> epvt = {0.0044, 0.044, 0.44, 4.4};
    double res = track_entry(&pdg, &trackid, &parentid, &initialEkin, &trackidT5, &epvt);
    std::cout << "Test results: " << res << endl;
}