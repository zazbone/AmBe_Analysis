from tkinter import Canvas
import ROOT 
from sys import argv
from utils import load_AmBe
from _path import *

load_AmBe()

df = ROOT.RDataFrame("A", argv[1])


columns = [ROOT.E_QUENCHED_ALL, ROOT.E_QUENCHED_CE, ROOT.E_QUENCHED_NOISE_CE, ROOT.E_QUENCHED_GAMMA, ROOT.E_QUENCHED_NOISE_GAMMA, ROOT.E_QUENCHED_NOISE]
title = {
    ROOT.E_QUENCHED_ALL: "Energy per cube",
    ROOT.E_QUENCHED_NOISE: "Noise energy per cube",
    ROOT.E_QUENCHED_CE: "CE energy per cube",
    ROOT.E_QUENCHED_NOISE_CE: "Noise CE energy per cube",
    ROOT.E_QUENCHED_GAMMA: "Gamma energy per cube",
    ROOT.E_QUENCHED_NOISE_GAMMA: "Noise Gamme energy per cube"
}
color = {
    ROOT.E_QUENCHED_ALL: ROOT.kBlack,
    ROOT.E_QUENCHED_NOISE: ROOT.kBlue,
    ROOT.E_QUENCHED_CE: ROOT.kRed,
    ROOT.E_QUENCHED_NOISE_CE: ROOT.kGreen,
    ROOT.E_QUENCHED_GAMMA: ROOT.kViolet,
    ROOT.E_QUENCHED_NOISE_GAMMA: 15
}

canvas = ROOT.TCanvas("canvas1", "Energy distribution", 1080, 720)
#ROOT.gPad.SetLogy(True)
for col in columns:
    h = df.Histo1D(ROOT.RDF.TH1DModel(col, title[col], 100, 0.5, 5), col)
    h.SetName(title[col])
    h.SetYTitle("Cube number")
    h.SetXTitle("E (MeV)")
    h.SetLineColor(color[col])
    h.SetStats(True)
    h.Draw()
    canvas.Draw()
    #input()
    canvas.SaveAs(f"{col}.png")
    canvas.Clear()


stack = ROOT.THStack("hist", "Energy per cube distribution")
stack.SetNameTitle("Energy per cube", "Energy per cube distribution")
histo = [df.Histo1D(ROOT.RDF.TH1DModel(col, title[col], 100, 0.5, 5), col) for col in columns[1::]]
ROOT.RDF.RunGraphs(histo)
for col, h in zip(columns[1::], histo):
    h.SetLineColor(color[col])
    h.SetFillColor(color[col])
    stack.Add(h.GetValue())
#stack.Draw("nostack")
stack.Draw()
stack.GetHistogram().SetName(f"Energy per cube for {df.Count().GetValue()} events")
stack.GetHistogram().GetXaxis().SetTitle("E (MeV)")
stack.GetHistogram().GetYaxis().SetTitle("Cube number")
stack.GetHistogram().SetStats(True)
input()
canvas.SaveAs("FullStackFilled.png")
canvas.Clear()

stack2 = ROOT.THStack("hist", "Energy per cube distribution")
stack2.SetNameTitle("Neutron energy per cube", "Energy per cube distribution")
histo = [df.Histo1D(ROOT.RDF.TH1DModel(col, title[col], 100, 0.5, 5), col) for col in [ROOT.E_QUENCHED_NOISE_CE, ROOT.E_QUENCHED_NOISE_GAMMA, ROOT.E_QUENCHED_NOISE]]
ROOT.RDF.RunGraphs(histo)
for col, h in zip([ROOT.E_QUENCHED_NOISE_CE, ROOT.E_QUENCHED_NOISE_GAMMA, ROOT.E_QUENCHED_NOISE], histo):
    if col == ROOT.E_QUENCHED_NOISE:
        h.SetLineColor(ROOT.kBlack)
    else:
        h.SetLineColor(ROOT.kWhite)
    h.SetFillColor(ROOT.kWhite)
    stack.Add(h.GetValue())
#stack.Draw("nostack")
stack.Draw()

stack.GetHistogram().SetName(f"Energy per cube for {df.Count().GetValue()} events")
stack.GetHistogram().GetXaxis().SetTitle("E (MeV)")
stack.GetHistogram().GetYaxis().SetTitle("Cube number")
stack.GetHistogram().SetStats(True)
input()
canvas.SaveAs("NetronStack.png")
canvas.Clear()