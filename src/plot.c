void savePlot(TObject * plt, char const* fileName) {
    auto pathBuffer = std::string("plots") + std::string(fileName);
    std::string pathPgf("plots/pgf");
    std::string pathPng("plots/png");
    TCanvas c;
    plt->DrawClone();
    c.Update();
    c.Print((std::string("plots/pdf/") + std::string(fileName) + std::string(".pdf")).c_str());
    c.Print((std::string("plots/pgf/") + std::string(fileName) + std::string(".pgf")).c_str());
    c.Print((std::string("plots/png/") + std::string(fileName) + std::string(".png")).c_str());
}