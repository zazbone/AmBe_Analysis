import ROOT 
from sys import argv
from utils import load_AmBe, to_cpp_str
from _path import *

load_AmBe()

df = ROOT.RDataFrame("A", argv[1])


dfg = df.Filter(f"{ROOT.VOLID_GAMMA}[{ROOT.E_QUENCHED_GAMMA} > 0.5].size() != 0")
dfg = dfg.Define("gamma_dist", f"euclidDist({ROOT.VOLID_GAMMA}[{ROOT.E_QUENCHED_GAMMA} > 0.5][0], {ROOT.VOLID_GAMMA}[{ROOT.E_QUENCHED_GAMMA} > 0.5])")

canvas = ROOT.TCanvas()
ROOT.gPad.SetLogy(True)
model = ROOT.RDF.TH1DModel("number of cube", "Distance (cube)", 100, 0, 20)
hist = dfg.Histo1D(model, "gamma_dist")
hist.SetYTitle("Cube number")
hist.SetXTitle("Distance (cube)")
hist.SetStats(True)
hist.Draw()
canvas.Draw()
input()
canvas.SaveAs(f"gamma_dist.png")
canvas.Clear()


dfce = df.Filter(f"{ROOT.VOLID_CE}[{ROOT.E_QUENCHED_CE} > 0.5].size() != 0")
dfce = dfce.Define("CE_dist", f"euclidDist({ROOT.VOLID_CE}[{ROOT.E_QUENCHED_CE} > 0.5][0], {ROOT.VOLID_CE}[{ROOT.E_QUENCHED_CE} > 0.5])")

model = ROOT.RDF.TH1DModel("number of cube", "Distance (cube)", 100, 0, 10)
hist = dfce.Histo1D(model, "CE_dist")
hist.SetYTitle("Cube number")
hist.SetXTitle("Distance (cube)")
hist.SetStats(True)
hist.Draw()
canvas.Draw()
input()
canvas.SaveAs(f"CE_dist.png")