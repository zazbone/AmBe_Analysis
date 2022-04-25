#include "TObject.h"
#include "TCanvas.h"
#include <string>


void savePlot(TObject const* plt, std::string& dir, std::string& fileName) {
    TCanvas c;
    plt->DrawClone();
    c.Update();
    c.Print((dir + std::string{"/pdf/"} + fileName + std::string{".pdf"}).c_str());
    c.Print((dir + std::string{"/pgf/"} + fileName + std::string{".pgf"}).c_str());
    c.Print((dir + std::string{"/png/"} + fileName + std::string{".png"}).c_str());
}