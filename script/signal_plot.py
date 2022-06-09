import ROOT 
from sys import argv
from utils import load_AmBe
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
df = df.Define("ENCE", "ROOT::VecOps::Sum(E_quenchedNoiseCE[most_energ_cube == volidNoiseCE])")
df = df.Define("ECE", "ROOT::VecOps::Sum(E_quenchedCE[most_energ_cube == volidCE])")
df = df.Define("Energy", "std::max(ENCE, ECE)")
df = df.Define("Noise", "E_quenchedAll - Energy")
columns = ["Energy", "Noise"]
color = {col: c for col, c in zip(columns, [ROOT.kRed, ROOT.kBlue])}

ROOT.gROOT.SetStyle("ATLAS")
canvas = ROOT.TCanvas("mainC", "Energy per cube distribution")
frame = canvas.DrawFrame(0.5, 0, 5, 140000)
frame.SetYTitle("# cube");
frame.SetXTitle("E (MeV)");
frame.GetYaxis().SetTitleOffset(1.4);
frame.GetXaxis().SetTitleOffset(1.4);
stack = ROOT.THStack("hist", "Energy per cube distribution")
stack.SetTitle(f"Energy per cube for {df.Count().GetValue()} events")
stack.SetNameTitle("Energy per cube", "Energy per cube distribution")
histo = [df.Histo1D(ROOT.RDF.TH1DModel(col, col, 100, 0.5, 5), col) for col in columns]
ROOT.RDF.RunGraphs(histo)
for col, h in zip(columns, histo):
    h.SetFillColor(color[col])
    stack.Add(h.GetValue())
stack.DrawClone()
canvas.Draw()
input()