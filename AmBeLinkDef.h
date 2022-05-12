// See: https://root.cern.ch/selecting-dictionary-entries-linkdefh
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ defined_in "include/AmBe.hpp";
#pragma link C++ defined_in "include/constants.hpp";
#pragma link C++ defined_in "include/create_columns.hpp";
#pragma link C++ defined_in "include/load_dataset.hpp";
#pragma link C++ defined_in "include/track_chain.hpp";
#pragma link C++ defined_in "include/transform.hpp";
#pragma link C++ defined_in "include/util.hpp";
#endif