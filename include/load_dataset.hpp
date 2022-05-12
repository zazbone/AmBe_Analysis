#pragma once

#include "TTree.h"
#include "TFile.h"
#include <ROOT/RDataFrame.hxx>
#include <vector>
#include <memory>
 

typedef struct {
    ROOT::RDF::RNode df;
    std::vector<std::shared_ptr<TTree>> tree;
    std::vector<std::shared_ptr<TFile>> file;
} DataSetHolder;

DataSetHolder ambeDataSet (char const* path);

DataSetHolder withS2(char const* G4path, char const* S2path);