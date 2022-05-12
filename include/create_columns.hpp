#pragma once

#include "TTree.h"
#include <ROOT/RDataFrame.hxx>


#include <memory>
#include <tuple>


#include "../include/transform.hpp"
#include "../include/track_chain.hpp"
#include "../include/util.hpp"


void EQuenched(const char* G4Path, const char* outputFile);