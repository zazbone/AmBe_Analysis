#pragma once


#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>


ROOT::RDF::RNode filterT1Events (ROOT::RDF::RNode ambedf, double absTol);

ROOT::RVec<int> _fromMaskT9(
    ROOT::RVec<int> const& trackidT9,
    ROOT::RVec<int> const& maskT9,
    ROOT::RVec<int> const& trackidT5
);

ROOT::RDF::RNode createT5Mask(ROOT::RDF::RNode df);

std::tuple<ROOT::RVec<int>, ROOT::RVec<double>> volidEpvtTot(ROOT::RVec<int> const& volid, ROOT::RVec<double> const& edep_pvt);

ROOT::RVec<double> euclidDist(int origine, ROOT::RVec<int> const& points);

ROOT::RVec<int> hammingDist(int origine, ROOT::RVec<int> const& points);

ROOT::RVec<int> neightborDist(int origine, ROOT::RVec<int> const& points);

ROOT::RVec<int> zDist(int origine, ROOT::RVec<int> const& points);

ROOT::RDF::RNode sumContribution(
    ROOT::RDF::RNode & df,
    const char* name,
    const char* selector,
    const char* energy_column,
    const char* volid_column
);

ROOT::RDF::RNode selectedCEEvents(ROOT::RDF::RNode & df);