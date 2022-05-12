#pragma once


#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

std::set<int> RVecToSet(ROOT::RVec<int> const& vec);

ROOT::RVec<bool> in(ROOT::RVec<int> const& tid, ROOT::RVec<int> const& trackidT5);

ROOT::RVec<int> getIndex(ROOT::RVec<int> const& vec1, ROOT::RVec<int> const& vec2);

int neightborDistCalc(int dx, int dy, int dz);

template<typename T>
bool contain(ROOT::RVec<T> const& sequence, T const& value);

template<typename T>
int indexMax(ROOT::RVec<T> const& vec);

int getCubeX(long m_volID);

int getCubeY(long m_volID);

int getCubeZ(long m_volID);

int centerOfMass(ROOT::RVec<int> const& volid, ROOT::RVec<double> const& weight);

int calcVolid(int x, int y, int z);

double weightedMean(ROOT::RVec<int> const& value, ROOT::RVec<double> const& weight);

ROOT::RVec<bool> approxEq(ROOT::RVec<double> v, double value, double eps);