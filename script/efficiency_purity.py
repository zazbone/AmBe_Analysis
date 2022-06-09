import ROOT 
from sys import argv
from utils import load_AmBe, to_cpp_str
from _path import *

load_AmBe()

df = ROOT.RDataFrame("A", argv[1])
Ntot = df.Count().GetValue()
df = df.Filter("E_quenchedAll[0] > 3.7")
Nsel = df.Count().GetValue()

df = df.Define("most_energ_cube", "volidAll[0]")

df = df.Define("most_energ_CEcube", "volidCE.size() ? volidCE[0] : -1")
df = df.Define("most_energ_Noisecube", "volidNoiseCE.size() ? volidNoiseCE[0] : -1")
Nsel_true = df.Filter("most_energ_CEcube == most_energ_cube || most_energ_Noisecube == most_energ_cube").Count().GetValue()

print("All event:                                  ", Ntot)
print("Number of selected event:                   ", Nsel)
print("Number of selected cube with a Compton edge:", Nsel_true)
print("Efficiency:", Nsel / Ntot)
print("Purity:    ", Nsel_true / Nsel)

print("######## Energy ########")

df = df.Define("ENCE", "ROOT::VecOps::Sum(E_quenchedNoiseCE[most_energ_cube == volidNoiseCE])")
df = df.Define("ECE", "ROOT::VecOps::Sum(E_quenchedCE[most_energ_cube == volidCE])")
df = df.Define("Energy", "std::max(ENCE, ECE)")
df = df.Define("delta_E", "Energy > 0.05 ? (Energy - E_quenchedAll[0]) / Energy : -1")

canvas = ROOT.TCanvas()
ROOT.gPad.SetLogy(False)
model = ROOT.RDF.TH1DModel("dN/dE", "dN/dE", 100, -1, 0.2)
hist = df.Histo1D(model, "delta_E")
hist.SetYTitle("Number of events")
hist.SetXTitle("dN/dE")
hist.SetStats(True)
hist.Draw()
canvas.Draw()
input()
canvas.SaveAs(f"delta_E.png")
canvas.Clear()