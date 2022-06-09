import ROOT
from _path import *


def load_AmBe():
    ROOT.gInterpreter.AddIncludePath(str(INCLUDE_PATH))
    ROOT.gInterpreter.ProcessLine('#include "AmBe.hpp"')  
    ROOT.gSystem.Load(str(DLL_PATH))
    ROOT.gROOT.ProcessLine(f".x {SCRIPT_PATH / 'lhcbstyle.C'}")


def to_cpp_str(py_obj):
    return ROOT.std.string(str(py_obj))