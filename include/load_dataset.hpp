#pragma once

#include "TTree.h"
#include <ROOT/RDataFrame.hxx>


typedef struct {
    ROOT::RDataFrame df;
    std::shared_ptr<TTree> tree;
    std::shared_ptr<TFile> file;
} DataSetHolder;

DataSetHolder ambeDataSet (char const* path);