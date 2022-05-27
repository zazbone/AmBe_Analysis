#include "../include/constants.hpp"
#include "../include/track_chain.hpp"
#include "../include/util.hpp"
#include "../include/column.hpp"


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
    auto dfT5view = ambedf.Define("initial_gamma", filter_formula);
    return dfT5view.Define("is_event", "ROOT::VecOps::Sum(initial_gamma) > 0");
}


ROOT::RDF::RNode createT5Mask(ROOT::RDF::RNode df) {
    return df.Define("chainMaskT5", "_fromMaskT9(T9.trackid, chainMaskT9, trackid)");
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


ROOT::RDF::RNode sumContribution(
    ROOT::RDF::RNode & df,
    const char* name,
    const char* selector,
    const char* energy_column,
    const char* volid_column
) {
    std::ostringstream volidColStream {};
    volidColStream << "volid" << name;
    std::string volidCol = volidColStream.str();
    std::ostringstream epvtColStream {};
    epvtColStream << "Epvt" << name;
    std::string epvtCol = epvtColStream.str();
    std::ostringstream exprStream {};
    exprStream<< "volidEpvtTot(" << volid_column << "[" << selector << "]" << "," << energy_column << "[" << selector << "]" << ")";
    std::string expr = exprStream.str();
    ROOT::RDF::RNode out = df.Define(volidCol, expr);
    using CubeTuple = std::tuple<ROOT::RVec<int>, ROOT::RVecD>;
    out = out.Define(epvtCol, [](CubeTuple const& volid){return std::get<1>(volid);}, {volidCol});
    out = out.Redefine(volidCol, [](CubeTuple const& volid){return std::get<0>(volid);}, {volidCol});
    return out;
}


// Require: maskT5, T1 event, T9. tree, volidCE, EpvtCE, volidALL, EpvtAll
ROOT::RDF::RNode selectedCEEvents(ROOT::RDF::RNode & df) {
    ROOT::RDF::RNode out = df.Filter("ROOT::VecOps::Any(T9.initialEkin[chainMaskT9 == maskCE_ELECT] > 4)", "CE produced with Ekin > 4");
    out = out.Filter("(EpvtCE[EpvtCE > 0.1]).size() == 1", "1 Compton edge cube");
    out = out.Define("CECubeVolid", "volidCE[EpvtCE > 0.1][0]");
    out = out.Define("CECubeCEEpvt", "ROOT::VecOps::Sum(EpvtCE[EpvtCE > 0.1])");
    out = out.Filter("auto mask = EpvtCE > 0.1; auto cubeVolid = volid == volidCE[mask].at(0) ; return EpvtCE[mask][0] / E_quenched[cubeVolid][0] > 0.95;", "Signal purity above 95%");
    return out;
}