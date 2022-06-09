#pragma once

#include "TTree.h"
#include <ROOT/RDataFrame.hxx>


#include <memory>
#include <tuple>
#include <string_view>

#include "../include/transform.hpp"
#include "../include/track_chain.hpp"
#include "../include/util.hpp"
#include "../include/constants.hpp"

// Mask column name
constexpr char const* MASK_T9 = "chainMaskT9";
constexpr char const* MASK_T5 = "chainMaskT5";
// Base T9 column
constexpr char const* T9_PARENT_ID = "T9.parentid";
constexpr char const* T9_TRACK_ID = "T9.trackid";
constexpr char const* T9_PDG = "T9.pdg";
constexpr char const* T9_CREATOR_PROCESS_ID = "T9.CreatorProcessID";
constexpr char const* T9_INITIAL_EKIN = "T9.initialEkin";
// Base T5 column
constexpr char const* T5_PDG = "pdg";
constexpr char const* T5_VOLID = "volid";
constexpr char const* T5_E_QUENCHED = "E_quenched";
constexpr char const* T5_TRACK_ID = "trackid";

// Column groupby volid and sorted by E_quenched
// For all
constexpr char const* VOLID_ALL = "volidAll";
constexpr char const* E_QUENCHED_ALL = "E_quenchedAll";
// For Neutron
constexpr char const* VOLID_NOISE = "volidNoise";
constexpr char const* E_QUENCHED_NOISE = "E_quenchedNoise";
// For Gamma
constexpr char const* VOLID_GAMMA = "volidGamma";
constexpr char const* E_QUENCHED_GAMMA = "E_quenchedGamma";
// For CE
constexpr char const* VOLID_CE = "volidCE";
constexpr char const* E_QUENCHED_CE = "E_quenchedCE";
// For Noise Gamma
constexpr char const* VOLID_NOISE_GAMMA = "volidNoiseGamma";
constexpr char const* E_QUENCHED_NOISE_GAMMA = "E_quenchedNoiseGamma";
// For Noise CE
constexpr char const* VOLID_NOISE_CE = "volidNoiseCE";
constexpr char const* E_QUENCHED_NOISE_CE = "E_quenchedNoiseCE";

// Total Energy of a process
constexpr char const* TOTAL_ALL = "totalEnergyAll";
constexpr char const* TOTAL_NOISE = "totalEnergyNoise";
constexpr char const* TOTAL_GAMMA = "totalEnergyGamma";
constexpr char const* TOTAL_CE = "totalEnergyCE";
constexpr char const* TOTAL_NOISE_GAMMA = "totalEnergyNoiseGamma";
constexpr char const* TOTAL_NOISE_CE = "totalEnergyNoiseCE";


void EQuenched(const char* G4Path, const char* outputFile);