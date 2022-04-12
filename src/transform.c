// Contains all main function relative to our research that apply change to dataset


const int GAMMA_PDG = 22;
const double VALIDE_ENERGY = 4.44;

// Require:
//  [ambedf (RDataFrame)]:: Initial dataframe with `T1.pdg` && `T1.energy` required
//  [absTolerence (double)]:: floating point tolerence for approx_eq test of event energy
// Return the modified RDataFrame with two new collumns:
//  [gamma_event (RVec<int>)]:: Add flag to each individual particles created (1) if pdg == 22 && energy ~ 4.44 (MEV)
//  [is_event (bool)]:: One per row, true .if. one gamma_event in row .else. false
ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void> filterT1Events (ROOT::RDataFrame ambedf, double absTol) {
    // Event filter expression
    std::ostringstream expr2;
    expr2 << "(T1.pdg == " << GAMMA_PDG << ")  && (ROOT::VecOps::abs(T1.energy - " << VALIDE_ENERGY << ") <" << absTol << ")";
    std::string filter_formula = expr2.str();
    auto dfT5view = ambedf.Define("gamma_event", filter_formula);
    return dfT5view.Define("is_event", "ROOT::VecOps::Sum(gamma_event) > 0");
}
