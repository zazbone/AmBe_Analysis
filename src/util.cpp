std::set<int> RVecToSet(ROOT::RVec<int> vec) {
    std::set<int> out;
    for (int el: vec) {
        out.insert(el);
    }
    return out;
}