#include "../include/util.hpp"


std::set<int> RVecToSet(ROOT::RVec<int> const& vec) {
    std::set<int> out;
    for (int el: vec) {
        out.insert(el);
    }
    return out;
}


ROOT::RVec<bool> in(ROOT::RVec<int> const& tid, ROOT::RVec<int> const& trackidT5) {
    ROOT::RVec<bool > out(trackidT5.size(), false);
    for (int i: tid) {
        out = out || (trackidT5 == i);
    }
    return out;
}

template<typename T>
bool contain(ROOT::RVec<T> const& sequence, T const& value) {
    for (auto i: sequence) {
        if (i == value) {return true;}
    }
    return false;
}



ROOT::RVec<int> getIndex(ROOT::RVec<int> const& vec1, ROOT::RVec<int> const& vec2) {
    ROOT::RVec<int> out {};
    std::set<int> memory {};
    for (int el: vec1) {
        for (unsigned long index = 0; index < vec2.size(); index++) {
            if (el == vec2.at(index) && (memory.find(index) == memory.end())) {
                out.push_back(index);
                memory.emplace(index);
                break;         
            }
        }
    }
    return out;
}

int neightborDistCalc(int dx, int dy, int dz) {
    return std::max({std::abs(dx), std::abs(dy), std::abs(dz)});
}


// Return the index of the maximal value in the given vec
// Empty vec are supported, the function will return -1 insted of a valid index
template<typename T>
int indexMax(ROOT::RVec<T> const& vec) {
    const int N = vec.size();
    if (!N) {return -1;}
    int maxValue = vec[0];
    int index = 0;
    for (int i = 1; i < N; i++) {
        if (maxValue < vec[i]) {
            maxValue = vec[i];
            index = i;
        }
    }
    return index;
}
template int indexMax<int>(ROOT::RVec<int> const& vec);
template int indexMax<double>(ROOT::RVec<double> const& vec);

int getCubeX(long m_volID)
{
    if (m_volID >= 1e9)m_volID -= 1e9;
    return(((m_volID/10)%10000)/100);
}

int getCubeY(long m_volID)
{
    if (m_volID >= 1e9)m_volID -= 1e9;
    return((m_volID/10)%100);
}

int getCubeZ(long m_volID)
{
    if (m_volID >= 1e9)m_volID -= 1e9;
    return(m_volID/100000);
}

int calcVolid(int x, int y, int z) {
    return 1e9 + z * 100000 + x * 1000 + y * 10;
}


int centerOfMass(ROOT::RVec<int> const& volid, ROOT::RVec<double> const& edep_pvt) {
    if (volid.size() == 0) {return -1;}
    ROOT::RVec<int> X = ROOT::VecOps::Map(volid, getCubeX);
    ROOT::RVec<int> Y = ROOT::VecOps::Map(volid, getCubeY);
    ROOT::RVec<int> Z = ROOT::VecOps::Map(volid, getCubeZ);
    int meanX = std::lround(weightedMean(X, edep_pvt));
    int meanY = std::lround(weightedMean(Y, edep_pvt));
    int meanZ = std::lround(weightedMean(Z, edep_pvt));
    return calcVolid(meanX, meanY, meanZ);
}


double weightedMean(ROOT::RVec<int> const& value, ROOT::RVec<double> const& weight) {
    double iW = 1 /  ROOT::VecOps::Sum(weight);
    return ROOT::VecOps::Sum(value * weight) * iW;
}