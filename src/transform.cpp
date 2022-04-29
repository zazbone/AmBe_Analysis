#include "../include/constants.hpp"
#include "../include/track_chain.hpp"
#include "../include/util.hpp"


// Contains all main function relative to our research that apply change to dataset

// Require:
//  [ambedf (RDataFrame)]:: Initial dataframe with `T1.pdg` && `T1.energy` required
//  [absTolerence (double)]:: floating point tolerence for approx_eq test of event energy
// Return the modified RDataFrame with two new collumns:
//  [gamma_event (RVec<int>)]:: Add flag to each individual particles created (1) if pdg == 22 && energy ~ 4.44 (MEV)
//  [is_event (bool)]:: One per row, true .if. one gamma_event in row .else. false
ROOT::RDF::RNode filterT1Events (ROOT::RDF::RNode ambedf, double absTol) {
    // Event filter expression
    std::ostringstream expr2;
    expr2 << "(T1.pdg == " << GAMMA_PDG << ")  && (ROOT::VecOps::abs(T1.energy - " << GAMMA_ENERGY << ") <" << absTol << ")";
    std::string filter_formula = expr2.str();
    auto dfT5view = ambedf.Define("gamma_event", filter_formula);
    return dfT5view.Define("is_event", "ROOT::VecOps::Sum(gamma_event) > 0");
}


ROOT::RVec<int> _fromMaskT9(
    ROOT::RVec<int> const& trackidT9,
    ROOT::RVec<int> const& maskT9,
    ROOT::RVec<int> const& trackidT5
) {
    int N = trackidT5.size();
    ROOT::RVec<int> maskT5 (N, 0);
    for (int i = 1; i < 5; i++) {
        ROOT::RVec<int> tid = trackidT9[maskT9==i];
        ROOT::RVec<bool> mask = in(tid, trackidT5);
        maskT5 = ROOT::VecOps::Where(mask, ROOT::RVec<int>(N, i), maskT5);
    }
    return maskT5;
}


ROOT::RDF::RNode createT5Mask(ROOT::RDF::RNode df) {
    return df.Define("chainMaskT5", "_fromMaskT9(T9.trackid, chainMaskT9, trackid)");
}

// Given volid and edep_pvt should have been already masked
std::tuple<ROOT::RVec<int>, ROOT::RVec<double>> volidEpvtTot(
    ROOT::RVec<int> const& volid, ROOT::RVec<double> const& edep_pvt
) {
    ROOT::RVec<int> totalVolid {};
    ROOT::RVec<double> totalEpvt {};
    const int N = volid.size();
    for (int i = 0; i < N; i++) {
        unsigned long j = 0;
        for (; j < totalVolid.size(); j++) {
            if (totalVolid[j] == volid[i]) {
                totalEpvt[j] += edep_pvt[i];
                break;
            }
        }
        if (j >= totalVolid.size()) {
            totalVolid.push_back(volid[i]);
            totalEpvt.push_back(edep_pvt[i]);
        }
    }
    return std::make_tuple(totalVolid, totalEpvt);
}


ROOT::RVec<double> euclidDist(int origine, ROOT::RVec<int> const& points) {
    int x = getCubeX(origine);
    int y = getCubeY(origine);
    int z = getCubeZ(origine);
    ROOT::RVec<int> X = ROOT::VecOps::Map(points, getCubeX);
    ROOT::RVec<int> Y = ROOT::VecOps::Map(points, getCubeY);
    ROOT::RVec<int> Z = ROOT::VecOps::Map(points, getCubeZ);
    ROOT::RVec<int> dx = X - x;
    ROOT::RVec<int> dy = Y - y;
    ROOT::RVec<int> dz = Z - z;
    return ROOT::VecOps::sqrt(dx * dx + dy * dy + dz * dz);
}

ROOT::RVec<int> hammingDist(int origine, ROOT::RVec<int> const& points) {
    int x = getCubeX(origine);
    int y = getCubeY(origine);
    int z = getCubeZ(origine);
    ROOT::RVec<int> X = ROOT::VecOps::Map(points, getCubeX);
    ROOT::RVec<int> Y = ROOT::VecOps::Map(points, getCubeY);
    ROOT::RVec<int> Z = ROOT::VecOps::Map(points, getCubeZ);
    ROOT::RVec<int> dx = ROOT::VecOps::Where(X - x > 0, X - x, x - X);
    ROOT::RVec<int> dy = ROOT::VecOps::Where(Y - y > 0, Y - y, y - Y);
    ROOT::RVec<int> dz = ROOT::VecOps::Where(Z - z > 0, Z - z, z - Z);
    return dx + dy + dz;
}

ROOT::RVec<int> neightborDist(int origine, ROOT::RVec<int> const& points) {
    int x = getCubeX(origine);
    int y = getCubeY(origine);
    int z = getCubeZ(origine);
    ROOT::RVec<int> X = ROOT::VecOps::Map(points, getCubeX);
    ROOT::RVec<int> Y = ROOT::VecOps::Map(points, getCubeY);
    ROOT::RVec<int> Z = ROOT::VecOps::Map(points, getCubeZ);
    ROOT::RVec<int> dx = X - x;
    ROOT::RVec<int> dy = Y - y;
    ROOT::RVec<int> dz = Z - z;
    return ROOT::VecOps::Map(dx, dy, dz, neightborDistCalc);
}


ROOT::RVec<int> zDist(int origine, ROOT::RVec<int> const& points) {
    int z = getCubeZ(origine);
    ROOT::RVec<int> Z = ROOT::VecOps::Map(points, getCubeZ);
    ROOT::RVec<int> dz = Z - z;
    return dz;
}