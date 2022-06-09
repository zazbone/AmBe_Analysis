import ROOT
from sys import argv
from _path import *
import os


if len(argv) < 3:
    # TODO: Find the write exception to raise
    raise ("Wrong number of argument. Please provide the 2 imput file path [and output path (optional)]")

input_file1 = Path(argv[1])
input_file2 = Path(argv[2])
if len(argv) == 4:
    output_file = Path(argv[3])
else:
    print("Output path not provided, tmp folder created")
    os.mkdir("./tmp")
    name = input_file1.with_suffix(".txt").name
    output_file = Path("./tmp") / f"{name}"
input_file1 = str(input_file1.resolve())
input_file2 = str(input_file2.resolve())
output_file = str(output_file.resolve())

f = ROOT.TFile.Open(input_file1)
main_tree = f.Get("T5")
main_tree.AddFriend("T1", input_file1)
main_tree.AddFriend("T9", input_file1)
main_tree.AddFriend("A", input_file2)
df = ROOT.RDataFrame(main_tree)
with open(output_file, "w+") as out:
    out.write(df.Display(
        (ROOT.MASK_T9, ROOT.T9_PDG, ROOT.T9_INITIAL_EKIN, ROOT.T9_TRACK_ID, ROOT.T5_TRACK_ID, ROOT.MASK_T5, ROOT.T5_E_QUENCHED,
        ROOT.E_QUENCHED_NOISE_CE, ROOT.E_QUENCHED_NOISE_GAMMA, ROOT.TOTAL_NOISE_CE, ROOT.TOTAL_NOISE_GAMMA), 100, 1000).AsString())
f.Close()